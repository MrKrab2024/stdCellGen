#!/usr/bin/env python3
# MIT License
# Source/Drain orientation preprocessor for SPICE std-cell decks.

import argparse
import json
import os
import re
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Set, Optional

# -----------------------------
# Data models
# -----------------------------

@dataclass
class Device:
    name: str
    d: str
    g: str
    s: str
    b: str
    model: str
    params: str        # trailing text after model, preserved
    typ: str           # "pmos" or "nmos" or "unknown"
    line_idx: int      # original file line index
    subckt: str        # parent cell name

@dataclass
class Group:
    id: int
    dev_ids: List[int] = field(default_factory=list)
    ext_nets: Tuple[str, str] = ("", "")       # two external DS nets (ends)
    _internal_nets: Set[str] = field(default_factory=set)
    order: List[int] = field(default_factory=list)  # device indices in chain order
    merged_from: List[int] = field(default_factory=list)  # original group ids if merged

@dataclass
class Cell:
    name: str
    pins: List[str]
    start_line: int
    end_line: int
    devices: List[Device] = field(default_factory=list)

@dataclass
class OrientationResult:
    # key: "subckt:devname" -> (newD, newS)
    dev_to_oriented_ds: Dict[str, Tuple[str, str]] = field(default_factory=dict)
    groups_oriented: int = 0
    groups_ambiguous: int = 0
    devs_oriented: int = 0
    devs_ambiguous: int = 0
    sequences: List[List[str]] = field(default_factory=list)
    shared_nets: Set[str] = field(default_factory=set)
    notes: List[str] = field(default_factory=list)
    ambiguous_devices: List[Dict[str, str]] = field(default_factory=list)
    skipped_tg_devices: List[Dict[str, str]] = field(default_factory=list)
    debug_groups_pre: List[dict] = field(default_factory=list)
    debug_groups_post: List[dict] = field(default_factory=list)
    oriented_group_dirs: List[dict] = field(default_factory=list)
    components: List[dict] = field(default_factory=list)

# -----------------------------
# SPICE parsing
# -----------------------------

SUBCKT_RE = re.compile(r"^\s*\.SUBCKT\s+(\S+)\s*(.*)$", re.IGNORECASE)
ENDS_RE   = re.compile(r"^\s*\.ENDS\b", re.IGNORECASE)
# Matches: M<name> D G S B model [params...]
M_RE      = re.compile(r"^\s*M(\S*)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(.*)$", re.IGNORECASE)

def parse_spice(sp_path: str) -> Tuple[List[str], Dict[str, Cell]]:
    with open(sp_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    cells: Dict[str, Cell] = {}
    cur_cell: Optional[Cell] = None

    for i, raw in enumerate(lines):
        line = raw.strip()
        if not line or line.startswith("*"):
            continue
        m = SUBCKT_RE.match(line)
        if m:
            name = m.group(1)
            pins = m.group(2).split()
            cur_cell = Cell(name=name, pins=pins, start_line=i, end_line=-1, devices=[])
            cells[name] = cur_cell
            continue
        if ENDS_RE.match(line):
            if cur_cell:
                cur_cell.end_line = i
                cur_cell = None
            continue
        if cur_cell:
            mm = M_RE.match(line)
            if mm:
                mname = "M" + mm.group(1)
                d, g, s, b = mm.group(2), mm.group(3), mm.group(4), mm.group(5)
                model = mm.group(6)
                params = mm.group(7).strip()
                mlow = model.lower()
                typ = "pmos" if "pmos" in mlow else ("nmos" if "nmos" in mlow else "unknown")
                dev = Device(
                    name=mname, d=d, g=g, s=s, b=b, model=model, params=params,
                    typ=typ, line_idx=i, subckt=cur_cell.name
                )
                cur_cell.devices.append(dev)
    return lines, cells

# -----------------------------
# Utilities
# -----------------------------

class UnionFind:
    def __init__(self, n: int):
        self.p = list(range(n))
        self.r = [0]*n
    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.r[ra] < self.r[rb]:
            self.p[ra] = rb
        elif self.r[ra] > self.r[rb]:
            self.p[rb] = ra
        else:
            self.p[rb] = ra
            self.r[ra] += 1
        return True

# -----------------------------
# Series grouping (skip power/shared nets as links; merge parallel by ordered ends)
# -----------------------------

def build_series_groups(devs: List[Device], power_nets: Set[str], shared_nets: Optional[Set[str]] = None) -> Tuple[List[Group], Dict[int, int], List[str]]:
    """
    Build pure-series groups:
      - Union only via DS nets with degree==2, excluding power nets and shared_nets as linking points.
      - Endpoints arise naturally from degrees.
      - Extract an ordered chain within each component for orientation.
      - Merge parallel groups only when ordered external ends are identical (left,right must match).
    Returns:
      groups list, device_id -> group_id mapping, merge_notes
    """
    notes: List[str] = []
    if shared_nets is None:
        shared_nets = set()
    n = len(devs)
    if n == 0:
        return [], {}, notes

    # DS net -> device indices (exclude power for linking map)
    net_to_devs_nonpower: Dict[str, List[int]] = defaultdict(list)
    for i, d in enumerate(devs):
        if d.d not in power_nets:
            net_to_devs_nonpower[d.d].append(i)
        if d.s not in power_nets:
            net_to_devs_nonpower[d.s].append(i)

    uf = UnionFind(n)
    # Link only on nets with exactly two devices, and not in shared_nets
    for net, ids in net_to_devs_nonpower.items():
        if net in shared_nets:
            continue
        if len(ids) == 2:
            uf.union(ids[0], ids[1])

    comp: Dict[int, List[int]] = defaultdict(list)
    for i in range(n):
        comp[uf.find(i)].append(i)

    groups_raw: List[Group] = []
    dev_to_group_raw: Dict[int, int] = {}
    gid = 0

    for _, ids in comp.items():
        # Build DS incidence within this component (include power nets to allow endpoints)
        ds_count: Dict[str, int] = defaultdict(int)
        net_to_devs_all: Dict[str, List[int]] = defaultdict(list)
        for i in ids:
            d = devs[i]
            for net in (d.d, d.s):
                ds_count[net] += 1
                net_to_devs_all[net].append(i)

        # Natural endpoints selection
        degree1_nonpower = [n for n, cnt in ds_count.items() if cnt == 1 and n not in power_nets]
        power_present = [n for n in ds_count.keys() if n in power_nets]

        if len(degree1_nonpower) >= 2:
            ends = [degree1_nonpower[0], degree1_nonpower[1]]
        elif len(degree1_nonpower) == 1 and power_present:
            ends = [degree1_nonpower[0], power_present[0]]
        elif len(power_present) >= 2:
            ends = [power_present[0], power_present[1]]
        else:
            nets_all = list(ds_count.keys())
            if len(nets_all) >= 2:
                a = nets_all[0]
                b = next((x for x in nets_all[1:] if x != a), a)
                ends = [a, b]
            elif len(nets_all) == 1:
                ends = [nets_all[0], nets_all[0]]
            else:
                ends = ["", ""]

        # Build an ordered chain (device order) starting from ends[0]
        order: List[int] = []
        visited: Set[int] = set()
        current_net = ends[0]
        steps = 0
        max_steps = len(ids) + 5

        while steps < max_steps and len(visited) < len(ids):
            cand = [i for i in net_to_devs_all.get(current_net, []) if i in ids and i not in visited]
            if not cand:
                remaining = [i for i in ids if i not in visited]
                if not remaining:
                    break
                i0 = remaining[0]
                d0 = devs[i0]
                current_net = d0.d
                continue
            i = cand[0]
            order.append(i)
            visited.add(i)
            d = devs[i]
            next_net = d.s if d.d == current_net else d.d
            current_net = next_net
            steps += 1

        # Internal nets: degree==2 and not power
        internal_nets: Set[str] = set(n for n, cnt in ds_count.items() if cnt == 2 and n not in power_nets)

        g = Group(
            id=gid,
            dev_ids=list(ids),
            ext_nets=(ends[0], ends[1]),
            _internal_nets=internal_nets,
            order=order,
            merged_from=[gid]
        )
        groups_raw.append(g)
        for i in ids:
            dev_to_group_raw[i] = gid
        gid += 1

    # Merge parallel groups only when ordered external ends match exactly
    key_to_groups: Dict[Tuple[str, str], List[Group]] = defaultdict(list)
    for g in groups_raw:
        a, b = g.ext_nets
        key = (a, b)
        key_to_groups[key].append(g)

    groups: List[Group] = []
    dev_to_group: Dict[int, int] = {}
    new_gid = 0
    for key, glist in key_to_groups.items():
        if len(glist) == 1:
            g = glist[0]
            g.id = new_gid
            groups.append(g)
            for i in g.dev_ids:
                dev_to_group[i] = new_gid
            new_gid += 1
        else:
            # merge
            all_devs: List[int] = []
            internal_union: Set[str] = set()
            merged_from = []
            for g in glist:
                all_devs.extend(g.dev_ids)
                internal_union |= g._internal_nets
                merged_from.extend(g.merged_from if g.merged_from else [g.id])
            a, b = key
            mg = Group(
                id=new_gid,
                dev_ids=all_devs,
                ext_nets=(a, b),
                _internal_nets=internal_union,
                order=[],
                merged_from=merged_from,
            )
            groups.append(mg)
            for i in all_devs:
                dev_to_group[i] = new_gid
            notes.append(f"merged {len(glist)} parallel groups into one (ordered): ends=({a}, {b}), devs={len(all_devs)}")
            new_gid += 1

    return groups, dev_to_group, notes

# -----------------------------
# Shared nets (between PMOS/NMOS DS sets)
# -----------------------------

def find_shared_ds_nets(p_devs: List[Device], n_devs: List[Device], exclude: Set[str]) -> Set[str]:
    def ds_nets(devs: List[Device]) -> Set[str]:
        s = set()
        for d in devs:
            if d.d not in exclude: s.add(d.d)
            if d.s not in exclude: s.add(d.s)
        return s
    return ds_nets(p_devs) & ds_nets(n_devs)

# -----------------------------
# Group graph and components
# -----------------------------

def build_group_graph(groups: List[Group]) -> Dict[str, List[Tuple[int, str]]]:
    """ net -> list of (group_id, other_end_net) """
    adj: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    for g in groups:
        a, b = g.ext_nets
        if a:
            adj[a].append((g.id, b))
        if b:
            adj[b].append((g.id, a))
    return adj

def split_components(groups: List[Group]) -> List[Dict]:
    """Split into connected components on the net-group graph.
    Returns a list of dicts: {"nets": set(str), "group_ids": set(int)}
    """
    adj = build_group_graph(groups)
    all_nets: Set[str] = set()
    for g in groups:
        a, b = g.ext_nets
        if a: all_nets.add(a)
        if b: all_nets.add(b)

    visited: Set[str] = set()
    comps: List[Dict] = []

    for n0 in list(all_nets):
        if n0 in visited:
            continue
        if not n0:
            continue
        q = deque([n0])
        visited.add(n0)
        comp_nets: Set[str] = {n0}
        comp_groups: Set[int] = set()
        while q:
            u = q.popleft()
            for gid, v in adj.get(u, []):
                comp_groups.add(gid)
                if v and v not in visited:
                    visited.add(v)
                    comp_nets.add(v)
                    q.append(v)
        comps.append({"nets": comp_nets, "group_ids": comp_groups})
    return comps

# -----------------------------
# Transmission gate filter (skip groups whose both ends are middle nets)
# -----------------------------


def filter_transmission_groups(groups: List[Group], devs: List[Device], power_nets: Set[str], shared_nets: Set[str], skip_pairs: Optional[Set[frozenset]] = None) -> Tuple[List[Group], List[str], List[Dict[str, str]]]:
    notes: List[str] = []
    filtered: List[Group] = []
    skipped: List[Dict[str, str]] = []
    if skip_pairs is None:
        skip_pairs = set()
    new_gid = 0
    for g in groups:
        a, b = g.ext_nets
        if (a in shared_nets) and (b in shared_nets) and (frozenset((a, b)) in skip_pairs):
            notes.append('skip_tg_group gid={} ends=({}, {}) devs={}'.format(g.id, a, b, len(g.dev_ids)))
            for i in g.dev_ids:
                d = devs[i]
                skipped.append({
                    'device': d.name,
                    'group': str(g.id),
                    'ends': '{}--{}'.format(a, b),
                    'reason': 'transmission_gate_group'
                })
            continue
        g.id = new_gid
        filtered.append(g)
        new_gid += 1
    return filtered, notes, skipped

# -----------------------------
# Orientation helpers
# -----------------------------


def orient_chain_devices_in_group(group: Group,
                                  devs: List[Device],
                                  src_net: str,
                                  dst_net: str) -> Dict[str, Tuple[str, str]]:
    """
    Orient devices inside a series/parallel group.
    - If group.order exists (pure series), walk order from src_net to dst_net.
    - If order is empty (parallel-merged chains), enumerate simple device-paths
      from src_net to dst_net within this group, and orient each path along src->dst.
    Returns: { devname: (newD, newS) }
    """
    oriented: Dict[str, Tuple[str, str]] = {}

    # Helper: set orientation for a single step u->v over device index i
    def set_step(i: int, u: str, v: str):
        d = devs[i]
        # newS should be u, newD should be v
        newS, newD = u, v
        oriented[d.name] = (newD, newS)

    # Pure series with known device order
    if group.order:
        order = list(group.order)
        # Ensure first touches src_net; otherwise reverse
        def touches(i: int, net: str) -> bool:
            dd = devs[i]
            return (dd.d == net) or (dd.s == net)
        if order and not touches(order[0], src_net):
            order = list(reversed(order))
        cur = src_net
        for i in order:
            d = devs[i]
            # Determine the other end
            if d.d == cur:
                nxt = d.s
            elif d.s == cur:
                nxt = d.d
            else:
                # If the device does not touch current net, try aligning to dst
                if d.d == dst_net:
                    nxt = d.s
                    cur = d.d
                else:
                    nxt = d.d
                    cur = d.s
            set_step(i, cur, nxt)
            cur = nxt
        return oriented

    # Parallel-merged chains (order empty): enumerate device-simple paths from src to dst
    devset: Set[int] = set(group.dev_ids)
    # Build adjacency: net -> list of (device index, other net)
    net_to_edges: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    for i in group.dev_ids:
        d = devs[i]
        net_to_edges[d.d].append((i, d.s))
        net_to_edges[d.s].append((i, d.d))

    used: Set[int] = set()  # devices already oriented

    def dfs(u: str, visited_devs: Set[int], path: List[Tuple[int, str, str]]):
        # u: current net; path: list of (dev_index, from_net, to_net)
        if u == dst_net:
            # Orient this path (skip devices already done)
            for i, a, b in path:
                if i in used:
                    continue
                set_step(i, a, b)
                used.add(i)
            return True
        progressed = False
        for (i, v) in net_to_edges.get(u, []):
            if i in visited_devs or i in used:
                continue
            visited_devs.add(i)
            path.append((i, u, v))
            if dfs(v, visited_devs, path):
                progressed = True
            path.pop()
            visited_devs.remove(i)
        return progressed

    # Try to cover all devices by repeatedly finding disjoint paths src->dst
    # Prefer starting from src_net
    while True:
        before = len(used)
        dfs(src_net, set(), [])
        if len(used) == before:
            break

    # Any leftover devices (should not happen in clean parallel chains), attempt local fix:
    # orient them towards whichever endpoint they touch (src preferred)
    if len(used) < len(devset):
        for i in group.dev_ids:
            if i in used:
                continue
            d = devs[i]
            if d.d == src_net or d.s == src_net:
                other = d.s if d.d == src_net else d.d
                set_step(i, src_net, other)
            elif d.d == dst_net or d.s == dst_net:
                other = d.s if d.d == dst_net else d.d
                set_step(i, other, dst_net)
            else:
                # Fallback: arbitrary but stable
                set_step(i, d.s, d.d)
            used.add(i)

    return oriented

    adj = build_group_graph(groups)

    # Unordered pair -> group ids, for propagating orientation to parallel ends regardless of order
    pair_to_gids: Dict[frozenset, List[int]] = defaultdict(list)
    for gid in comp_group_ids:
        a, b = groups[gid].ext_nets
        pair_to_gids[frozenset((a, b))].append(gid)

    def neighbors(u: str):
        for gid, v in adj.get(u, []):
            if gid in comp_group_ids and v in comp_nets:
                yield gid, v

    group_oriented: Dict[int, Tuple[str, str]] = {}
    dev_oriented: Dict[str, Tuple[str, str]] = {}
    sequences: List[List[str]] = []
    conflicts: List[str] = []

    def confirm_path(seq: List[str]):
        for i_idx in range(1, len(seq)):
            a = seq[i_idx - 1]
            b = seq[i_idx]
            desired = (a, b)
            gids = pair_to_gids.get(frozenset((a, b)), [])
            for g2 in gids:
                if g2 in group_oriented:
                    prev = group_oriented[g2]
                    if prev != desired:
                        dev_names = [devs[i].name for i in groups[g2].dev_ids]
                        conflicts.append(
                            f"cell={subckt} group={g2} devs={dev_names} conflict: prev {prev} vs new {desired} on path {seq}"
                        )
                        continue
                else:
                    group_oriented[g2] = desired
                    oriented = orient_chain_devices_in_group(groups[g2], devs, desired[0], desired[1])
                    for dn, (newD, newS) in oriented.items():
                        key = f"{subckt}:{dn}"
                        dev_oriented[key] = (newD, newS)
        sequences.append(seq)

    # Queue entries are (current_net, path_nets)
    q = deque()
    q.append((start_net, [start_net]))

    # Light dedupe to avoid explosion: remember (tail_net, prev_net)
    seen_tail: Set[Tuple[str, str]] = set()

    while q:
        u, path = q.popleft()
        if len(group_oriented) == len(comp_group_ids):
            break
        for gid, v in neighbors(u):
            if v in path:
                notes.append(f"skip path-cycle {'->'.join(path)}->{v}")
                continue
            new_path = path + [v]
            if v in shared_nets:
                confirm_path(new_path)
                continue
            tail = (v, u)
            if tail in seen_tail:
                continue
            seen_tail.add(tail)
            q.append((v, new_path))

    # Diagnostics if never reached any shared
    comp_shared = shared_nets & comp_nets
    reached_targets = set(s[-1] for s in sequences) if sequences else set()
    if comp_shared and not (reached_targets & comp_shared):
        notes.append(
            f"no_path_to_shared_with_path_check: start={start_net}, comp_shared={sorted(list(comp_shared))}"
        )

    return group_oriented, dev_oriented, sequences, conflicts

# -----------------------------
# Per-side processing with ambiguous reasons
# -----------------------------


def orient_component(groups: List[Group],
                     devs: List[Device],
                     comp_nets: Set[str],
                     comp_group_ids: Set[int],
                     start_net: str,
                     shared_nets: Set[str],
                     subckt: str,
                     notes: List[str]) -> Tuple[Dict[int, Tuple[str, str]], Dict[str, Tuple[str, str]], List[List[str]], List[str]]:
    adj = build_group_graph(groups)

    # Unordered pair -> group ids, for propagating orientation to parallel ends regardless of order
    pair_to_gids: Dict[frozenset, List[int]] = defaultdict(list)
    for gid in comp_group_ids:
        a, b = groups[gid].ext_nets
        pair_to_gids[frozenset((a, b))].append(gid)

    def neighbors(u: str):
        for gid, v in adj.get(u, []):
            if gid in comp_group_ids and v in comp_nets:
                yield gid, v

    group_oriented: Dict[int, Tuple[str, str]] = {}
    dev_oriented: Dict[str, Tuple[str, str]] = {}
    sequences: List[List[str]] = []
    conflicts: List[str] = []

    def confirm_path(seq: List[str]):
        # Confirm each edge along seq
        for i in range(1, len(seq)):
            a = seq[i-1]
            b = seq[i]
            desired = (a, b)
            gids = pair_to_gids.get(frozenset((a, b)), [])
            for g2 in gids:
                prev = group_oriented.get(g2)
                if prev and prev != desired:
                    dev_names = [devs[i].name for i in groups[g2].dev_ids]
                    conflicts.append(f"cell={subckt} group={g2} devs={dev_names} conflict: prev {prev} vs new {desired} on path {seq}")
                    continue
                if not prev:
                    group_oriented[g2] = desired
                    oriented = orient_chain_devices_in_group(groups[g2], devs, desired[0], desired[1])
                    for dn, (newD, newS) in oriented.items():
                        key = f"{subckt}:{dn}"
                        dev_oriented[key] = (newD, newS)

    # BFS of paths from start_net; only prune cycles on the current path
    from collections import deque
    q = deque()
    q.append([start_net])
    seen_paths = 0
    while q:
        path = q.popleft()
        u = path[-1]
        for gid, v in neighbors(u):
            if v in path:
                # cycle on current path; skip
                continue
            new_path = path + [v]
            if v in shared_nets:
                confirm_path(new_path)
                sequences.append(new_path)
            else:
                q.append(new_path)
            seen_paths += 1

    return group_oriented, dev_oriented, sequences, conflicts
def orient_by_components(groups: List[Group],
                         devs: List[Device],
                         start_power_net: str,
                         shared_nets: Set[str],
                         subckt: str) -> OrientationResult:
    result = OrientationResult()
    comps = split_components(groups)
    # record components summary for debug
    result.components = [{"nets": sorted(list(c["nets"])), "group_ids": sorted(list(c["group_ids"]))} for c in comps]

    global_dev_oriented: Dict[str, Tuple[str, str]] = {}
    global_group_oriented: Dict[int, Tuple[str, str]] = {}

    comp_data = []  # per-comp: {group_ids, comp_shared, reached}

    for ci, comp in enumerate(comps):
        comp_nets: Set[str] = comp["nets"]
        comp_group_ids: Set[int] = comp["group_ids"]
        comp_shared = shared_nets & comp_nets
        if not comp_group_ids:
            continue
        adj = build_group_graph(groups)
        deg: Dict[str, int] = defaultdict(int)
        for n in comp_nets:
            for gid, v in adj.get(n, []):
                if gid in comp_group_ids and v in comp_nets:
                    deg[n] += 1
        endpoints = [n for n in comp_nets if deg[n] == 1]

        if start_power_net in comp_nets:
            start = start_power_net
        elif endpoints:
            start = endpoints[0]
            result.notes.append(f"component {ci} no power net; start at endpoint {start}")
        else:
            start = next(iter(comp_nets))
            result.notes.append(f"component {ci} no power/endpoint; start at arbitrary {start}")

        if not comp_shared:
            result.notes.append(f"component {ci} has no shared nets; orientation may remain partial")

        group_oriented, dev_oriented, seqs, conflicts = orient_component(
            groups, devs, comp_nets, comp_group_ids, start, comp_shared, subckt, result.notes
        )
        reached_shared = any(s and s[-1] in comp_shared for s in seqs)
        comp_data.append({
            "group_ids": set(comp_group_ids),
            "comp_shared": set(comp_shared),
            "reached": reached_shared,
        })

        for gid, ddir in group_oriented.items():
            if gid in global_group_oriented and global_group_oriented[gid] != ddir:
                result.groups_ambiguous += 1
                result.notes.append(f"cross-comp conflict gid={gid} prev={global_group_oriented[gid]} new={ddir}")
            else:
                global_group_oriented[gid] = ddir
        global_dev_oriented.update(dev_oriented)
        result.sequences.extend(seqs)
        for c in conflicts:
            result.groups_ambiguous += 1
            result.notes.append(c)

    result.dev_to_oriented_ds = global_dev_oriented
    # oriented group directions for debug
    result.oriented_group_dirs = [{"gid": gid, "from": ddir[0], "to": ddir[1]} for gid, ddir in global_group_oriented.items()]
    result.devs_oriented = len(global_dev_oriented)
    result.groups_oriented = len(global_group_oriented)

    # Ambiguous device details
    oriented_gids = set(global_group_oriented.keys())
    ambiguous_details: List[Dict[str, str]] = []

    gid_to_comp = {}
    for idx, cd in enumerate(comp_data):
        for gid in cd["group_ids"]:
            gid_to_comp[gid] = idx

    for gid, g in enumerate(groups):
        if gid in oriented_gids:
            continue
        ci = gid_to_comp.get(gid, None)
        if ci is None:
            reason = "group_not_in_any_component"
        else:
            cd = comp_data[ci]
            if not cd["comp_shared"]:
                reason = "component_no_shared_net"
            elif not cd["reached"]:
                reason = "no_path_to_shared_in_component"
            else:
                reason = "not_on_any_confirmed_path"
        for i in g.dev_ids:
            d = devs[i]
            ambiguous_details.append({
                "device": d.name,
                "group": str(gid),
                "ends": f"{g.ext_nets[0]}--{g.ext_nets[1]}",
                "reason": reason,
            })

    result.ambiguous_devices = ambiguous_details

    all_dev_keys = {f"{subckt}:{d.name}" for d in devs}
    result.devs_ambiguous = len(all_dev_keys - set(global_dev_oriented.keys()))
    return result

# -----------------------------
# Per-cell processing
# -----------------------------

def dump_groups(groups: List[Group], devs: List[Device]) -> List[dict]:
    out = []
    for g in groups:
        a, b = g.ext_nets
        dev_names = [devs[i].name for i in g.dev_ids]
        order_names = [devs[i].name for i in g.order] if g.order else []
        out.append({
            'gid': g.id,
            'ends': [a, b],
            'devices': dev_names,
            'order': order_names,
        })
    return out

def process_cell(cell: Cell, vdd: str, vss: str) -> Dict[str, OrientationResult]:
    p_devs = [d for d in cell.devices if d.typ == "pmos"]
    n_devs = [d for d in cell.devices if d.typ == "nmos"]
    power = {vdd, vss}

    shared = find_shared_ds_nets(p_devs, n_devs, exclude=power)

    p_groups, _, p_merge_notes = build_series_groups(p_devs, power_nets=power, shared_nets=shared)
    n_groups, _, n_merge_notes = build_series_groups(n_devs, power_nets=power, shared_nets=shared)

    # debug: groups before TG filtering
    p_debug_pre = dump_groups(p_groups, p_devs)
    n_debug_pre = dump_groups(n_groups, n_devs)
    # Compute end-pair intersection (unordered) across PMOS/NMOS; only these are true TGs
    def __end_pairs(_groups):
        _s = set()
        for _g in _groups:
            a, b = _g.ext_nets
            if a and b:
                _s.add(frozenset((a, b)))
        return _s
    skip_pairs = __end_pairs(p_groups) & __end_pairs(n_groups)
    # Skip TG-like groups only for shared end-pairs present on both sides

    p_groups, p_tg_notes, p_tg_skipped = filter_transmission_groups(p_groups, p_devs, power, shared, skip_pairs)
    n_groups, n_tg_notes, n_tg_skipped = filter_transmission_groups(n_groups, n_devs, power, shared, skip_pairs)
    p_res = orient_by_components(p_groups, p_devs, start_power_net=vdd, shared_nets=shared, subckt=cell.name)
    n_res = orient_by_components(n_groups, n_devs, start_power_net=vss, shared_nets=shared, subckt=cell.name)
    
    # debug: groups after TG filtering
    p_debug_post = dump_groups(p_groups, p_devs)
    n_debug_post = dump_groups(n_groups, n_devs)
    
    # attach debug groups
    p_res.debug_groups_pre = p_debug_pre
    p_res.debug_groups_post = p_debug_post
    n_res.debug_groups_pre = n_debug_pre
    n_res.debug_groups_post = n_debug_post
    
    # Attach TG-skipped devices
    p_res.skipped_tg_devices = p_tg_skipped
    n_res.skipped_tg_devices = n_tg_skipped

    def __candidate_names(groups, devs):
        idxs = set()
        for g in groups:
            for i in g.dev_ids:
                idxs.add(i)
        return {devs[i].name for i in idxs}

    p_candidates = __candidate_names(p_groups, p_devs)
    n_candidates = __candidate_names(n_groups, n_devs)

    def __oriented_names(res, subckt):
        keys = set(res.dev_to_oriented_ds.keys())
        names = set()
        for k in keys:
            parts = k.split(":", 1)
            if len(parts) == 2 and parts[0] == subckt:
                names.add(parts[1])
        return names

    p_oriented_names = __oriented_names(p_res, cell.name)
    n_oriented_names = __oriented_names(n_res, cell.name)

    p_res.devs_ambiguous = len(p_candidates - p_oriented_names)
    n_res.devs_ambiguous = len(n_candidates - n_oriented_names)


    # Attach shared nets and merge notes for reporting
    p_res.shared_nets = set(shared)
    n_res.shared_nets = set(shared)
    p_res.notes.extend(p_merge_notes)
    n_res.notes.extend(n_merge_notes)
    p_res.notes.extend(p_tg_notes)
    n_res.notes.extend(n_tg_notes)

    return {"pmos": p_res, "nmos": n_res}

# -----------------------------
# Rewrite SPICE with oriented D/S
# -----------------------------

def rewrite_spice(lines: List[str],
                  cells: Dict[str, Cell],
                  results: Dict[str, Dict[str, OrientationResult]],
                  out_path: str) -> None:
    dev_orient: Dict[str, Tuple[str, str]] = {}
    for cname, res_pair in results.items():
        for typ in ("pmos", "nmos"):
            dev_orient.update(res_pair[typ].dev_to_oriented_ds)

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as f:
        cur_subckt: Optional[str] = None
        for raw in lines:
            line = raw.rstrip("\n")
            s = line.strip()

            m_sub = SUBCKT_RE.match(s)
            if m_sub:
                cur_subckt = m_sub.group(1)
                f.write(line + "\n")
                continue
            if ENDS_RE.match(s):
                cur_subckt = None
                f.write(line + "\n")
                continue

            mm = M_RE.match(s)
            if cur_subckt and mm:
                mname = "M" + mm.group(1)
                d, g, s_net, b, model, params = (mm.group(2), mm.group(3), mm.group(4),
                                                 mm.group(5), mm.group(6), mm.group(7))
                key = f"{cur_subckt}:{mname}"
                if key in dev_orient:
                    newD, newS = dev_orient[key]
                    out = f"{mname} {newD} {g} {newS} {b} {model}{params}"
                    f.write(out + "\n")
                else:
                    f.write(line + "\n")
            else:
                f.write(line + "\n")

# -----------------------------
# Report / Log
# -----------------------------

def write_report(results: Dict[str, Dict[str, OrientationResult]], report_path: str) -> None:
    payload = {}
    for cname, res_pair in results.items():
        entry = {}
        for typ in ("pmos", "nmos"):
            r = res_pair[typ]
            entry[typ] = {
                "groups_oriented": r.groups_oriented,
                "groups_ambiguous": r.groups_ambiguous,
                "devices_oriented": r.devs_oriented,
                "devices_ambiguous": r.devs_ambiguous,
                "shared_nets": sorted(list(r.shared_nets)),
                "sequences": r.sequences,
                "notes": r.notes,
                "ambiguous_devices": r.ambiguous_devices,
                "skipped_tg_devices": r.skipped_tg_devices,
                "series_groups_pre": r.debug_groups_pre,
                "series_groups_post": r.debug_groups_post,
                "oriented_group_dirs": r.oriented_group_dirs,
                "components": r.components,
            }
        payload[cname] = entry
    report_dir = os.path.dirname(report_path)
    if report_dir:
        os.makedirs(report_dir, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

def write_log(results: Dict[str, Dict[str, OrientationResult]], cells: Dict[str, Cell], log_path: str) -> None:
    processed = []
    skipped = []
    for cname, res_pair in results.items():
        oriented = res_pair["pmos"].devs_oriented + res_pair["nmos"].devs_oriented
        ambiguous = res_pair["pmos"].devs_ambiguous + res_pair["nmos"].devs_ambiguous
        if oriented > 0:
            processed.append((cname, oriented, ambiguous, (len(res_pair["pmos"].skipped_tg_devices) + len(res_pair["nmos"].skipped_tg_devices))))
        else:
            p_notes = "; ".join(res_pair["pmos"].notes)
            n_notes = "; ".join(res_pair["nmos"].notes)
            skipped.append((cname, p_notes or n_notes or "no orientation applied"))
    log_dir = os.path.dirname(log_path)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("[orient_sd] processed cells\n")
        for cname, o, a, tgs in processed:
            f.write(f"  - {cname}: oriented={o}, ambiguous={a}, tg_skipped={tgs}\n")
        f.write("\n[orient_sd] skipped cells\n")
        for cname, reason in skipped:
            f.write(f"  - {cname}: {reason}\n")

# -----------------------------
# CLI
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Preprocess SPICE: determine MOS Source/Drain by multi-length series groups and path-first orientation to shared nets."
    )
    ap.add_argument("--sp", required=True, help="Input SPICE (e.g., data/AsAp7.sp)")
    ap.add_argument("--vdd", default="VDD", help="VDD net name")
    ap.add_argument("--vss", default="VSS", help="VSS net name")
    ap.add_argument("--cells", default="", help="Optional file with cell names to process (one per line); default: all")
    ap.add_argument("--out_sp", default="out/oriented.sp", help="Output oriented SPICE path")
    ap.add_argument("--report_json", default="out/orient_report.json", help="Output JSON report path")
    ap.add_argument("--log", default="out/orient_log.txt", help="Output log path")
    args = ap.parse_args()

    lines, cells = parse_spice(args.sp)

    if args.cells:
        with open(args.cells, "r", encoding="utf-8") as f:
            targets = {ln.strip() for ln in f if ln.strip()}
    else:
        targets = set(cells.keys())

    results: Dict[str, Dict[str, OrientationResult]] = {}
    for cname in sorted(cells.keys()):
        if cname not in targets:
            continue
        res = process_cell(cells[cname], vdd=args.vdd, vss=args.vss)
        results[cname] = res

    rewrite_spice(lines, cells, results, args.out_sp)
    write_report(results, args.report_json)
    write_log(results, cells, args.log)

    total = oriented = ambig = 0
    for _, res_pair in results.items():
        pd = res_pair["pmos"].devs_oriented + res_pair["nmos"].devs_oriented
        pa = res_pair["pmos"].devs_ambiguous + res_pair["nmos"].devs_ambiguous
        total += pd + pa
        oriented += pd
        ambig += pa
    print(f"[orient_sd] devices total={total}, oriented={oriented}, ambiguous={ambig}")
    print(f"[orient_sd] oriented SPICE: {args.out_sp}")
    print(f"[orient_sd] report JSON:   {args.report_json}")
    print(f"[orient_sd] log:           {args.log}")

if __name__ == "__main__":
    main()


