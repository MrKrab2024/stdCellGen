#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import argparse
from collections import defaultdict, deque

class Trans:
    __slots__ = ("id","d","g","s","b","model","type")
    def __init__(self, id, d, g, s, b, model, type_):
        self.id=id; self.d=d; self.g=g; self.s=s; self.b=b; self.model=model; self.type=type_

def is_pmos_model(model: str) -> bool:
    m = model.lower()
    return ("pmos" in m) or ("pfet" in m) or (m.startswith("p"))

def parse_spice_subckts(path: str):
    """
    Return dict: subckt_name -> list[Trans]
    Minimal parser: handles .SUBCKT ... to .ENDS; transistor lines 'M...' with D G S B model ...
    Joins '+' continuation lines; ignores comments starting '*'.
    """
    subckts = {}
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    # join continuations:
    joined = []
    buf = ""
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith('*'):
            if buf:
                joined.append(buf)
                buf = ""
            continue
        if line.startswith('+'):
            buf += " " + line[1:].strip()
        else:
            if buf:
                joined.append(buf)
            buf = line
    if buf:
        joined.append(buf)

    current = None
    in_subckt = False
    for line in joined:
        l = line.strip()
        if not l:
            continue
        low = l.lower()
        if low.startswith(".subckt"):
            parts = l.split()
            if len(parts) >= 2:
                current = parts[1]
                subckts[current] = []
                in_subckt = True
            continue
        if low.startswith(".ends"):
            in_subckt = False
            current = None
            continue
        if not in_subckt or current is None:
            continue
        # transistor line: Mxxx D G S B model ...
        if l[0] in ('M','m'):
            toks = l.split()
            if len(toks) < 6:
                continue
            mid = toks[0]
            d, g, s, b = toks[1], toks[2], toks[3], toks[4]
            model = toks[5]
            ttype = "PMOS" if is_pmos_model(model) else "NMOS"
            subckts[current].append(Trans(mid, d, g, s, b, model, ttype))
    return subckts

def build_ds_adjacency(trans_list):
    """
    Build net -> list of device indices for DS (source/drain) terminals.
    Return (net2devs, deg, devs_by_idx)
    """
    net2devs = defaultdict(list)
    for i, t in enumerate(trans_list):
        net2devs[t.s].append(i)
        net2devs[t.d].append(i)
    deg = {net: len(devs) for net, devs in net2devs.items()}
    return net2devs, deg

def extract_series_chains_same_polarity(trans_list):
    """
    For a given polarity list (only PMOS or only NMOS), compute all maximal DS series chains:
    A chain is a maximal sequence of devices connected by DS nets where all internal nets have degree == 2.
    Return list of dicts: { 'dev_idx': [i...], 'start_net': net, 'end_net': net, 'gate_seq': [g...], 'dev_ids':[id...] }
    Orientation rule: chain is oriented so that first device is traversed from its source to its drain,
    i.e., start_net == first_device.s and end_net == last_device.d after orientation.
    """
    net2devs, deg = build_ds_adjacency(trans_list)
    n = len(trans_list)
    assigned = [False]*n
    chains = []

    # Helper: extend from given end net outward while internal net degree == 2
    def extend_from(net_end, prev_dev_idx):
        devs = []
        curr_net = net_end
        prev_dev = prev_dev_idx
        seen_nets = set()
        while True:
            if curr_net in seen_nets:
                break
            seen_nets.add(curr_net)
            if deg.get(curr_net, 0) != 2:
                break
            # get the other device at curr_net
            candidates = net2devs.get(curr_net, [])
            next_dev = None
            for di in candidates:
                if di != prev_dev and not assigned[di]:
                    next_dev = di
                    break
            if next_dev is None:
                break
            devs.append(next_dev)
            # advance net: the other terminal of next_dev
            t = trans_list[next_dev]
            next_net = t.d if (t.s == curr_net) else (t.s if t.d == curr_net else None)
            if next_net is None:
                break
            prev_dev = next_dev
            curr_net = next_net
        return devs, curr_net

    for i in range(n):
        if assigned[i]:
            continue
        t0 = trans_list[i]
        # Expand both ends
        left_devs, left_net = extend_from(t0.s, i)
        right_devs, right_net = extend_from(t0.d, i)

        # Combine to full chain: left (reversed order) + [i] + right
        full_devs = list(reversed(left_devs)) + [i] + right_devs
        # Determine nets at ends
        # We can reconstruct by stepping from left_end towards right_end
        start_net = left_net  # this is the net just outside the leftmost device on the left side
        end_net = right_net

        # Now orient chain so that first device is traversed from source->drain
        first = trans_list[full_devs[0]]
        if start_net != first.s:
            # reverse chain and swap ends
            full_devs = list(reversed(full_devs))
            start_net, end_net = end_net, start_net
            first = trans_list[full_devs[0]]
            # after reversing, still might not orient if nets are odd (should be rare),
            # but in well-formed stacks start_net should now equal first.s

        # Build gate sequence in chain order
        gate_seq = [trans_list[di].g for di in full_devs]
        dev_ids = [trans_list[di].id for di in full_devs]

        # Mark assigned
        for di in full_devs:
            assigned[di] = True

        chains.append({
            'dev_idx': full_devs,
            'start_net': start_net,
            'end_net': end_net,
            'gate_seq': gate_seq,
            'dev_ids': dev_ids,
        })

    return chains

def find_split_pairs(chains):
    """
    Given list of chains (same polarity), find split pairs:
    Two chains with same len, same gate_seq, same (start_net,end_net), and disjoint device sets.
    Return list of tuples (idx1, idx2)
    """
    pairs = []
    # Bucket by (len, tuple(gate_seq), start, end)
    buckets = defaultdict(list)
    for idx, ch in enumerate(chains):
        key = (len(ch['dev_idx']), tuple(ch['gate_seq']), ch['start_net'], ch['end_net'])
        buckets[key].append(idx)
    for key, idxs in buckets.items():
        if len(idxs) < 2:
            continue
        # choose disjoint pairs
        for i in range(len(idxs)):
            for j in range(i+1, len(idxs)):
                a = chains[idxs[i]]; b = chains[idxs[j]]
                if set(a['dev_idx']).isdisjoint(set(b['dev_idx'])):
                    pairs.append((idxs[i], idxs[j]))
    return pairs

def analyze_sp_file(sp_path, show_details=False):
    subckts = parse_spice_subckts(sp_path)
    result = {}
    for name, devices in subckts.items():
        # Separate by polarity
        p_list = [t for t in devices if t.type == "PMOS"]
        n_list = [t for t in devices if t.type == "NMOS"]
        p_chains = extract_series_chains_same_polarity(p_list) if p_list else []
        n_chains = extract_series_chains_same_polarity(n_list) if n_list else []
        p_pairs = find_split_pairs(p_chains)
        n_pairs = find_split_pairs(n_chains)
        has_split = bool(p_pairs or n_pairs)
        if has_split:
            result[name] = {
                'pmos_pairs': p_pairs,
                'nmos_pairs': n_pairs,
                'pmos_chains': p_chains if show_details else None,
                'nmos_chains': n_chains if show_details else None,
            }
    return result

def main():
    ap = argparse.ArgumentParser(description="Detect split-structure cells in a SPICE library (.SUBCKT).")
    ap.add_argument("--sp", required=True, help="Path to ASAP7.sp (or other .sp)")
    ap.add_argument("--details", action="store_true", help="Show chains and pair details")
    args = ap.parse_args()

    res = analyze_sp_file(args.sp, show_details=args.details)
    if not res:
        print("No split-structure cells found.")
        return 0

    print(f"Found {len(res)} cell(s) with split structure:")
    for cell, info in res.items():
        print(f"- {cell}")
        if args.details:
            # Print PMOS
            if info['pmos_pairs']:
                print(f"  PMOS pairs: {len(info['pmos_pairs'])}")
                for (i,j) in info['pmos_pairs']:
                    ci = info['pmos_chains'][i]; cj = info['pmos_chains'][j]
                    print(f"    pair(P): len={len(ci['dev_idx'])} gates={ci['gate_seq']} endpoints=({ci['start_net']}->{ci['end_net']})")
            # Print NMOS
            if info['nmos_pairs']:
                print(f"  NMOS pairs: {len(info['nmos_pairs'])}")
                for (i,j) in info['nmos_pairs']:
                    ci = info['nmos_chains'][i]; cj = info['nmos_chains'][j]
                    print(f"    pair(N): len={len(ci['dev_idx'])} gates={ci['gate_seq']} endpoints=({ci['start_net']}->{ci['end_net']})")
    return 0

if __name__ == "__main__":
    sys.exit(main())
