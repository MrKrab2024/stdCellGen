#include "stdcell/structural_matcher.hpp"

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <queue>
#include <set>
#include <limits>

namespace stdcell {

struct Arena {
    std::vector<MergeGroup> nodes;
    int add(const MergeGroup& g){ nodes.push_back(g); return (int)nodes.size()-1; }
};

static bool is_rail(const std::string& n) {
    std::string u = n; for (auto& c : u) c = (char)toupper((unsigned char)c);
    return (u=="VDD"||u=="VSS"||u=="VCC"||u=="GND");
}

struct NetUsage { std::vector<int> as_src, as_dst; };

static void index_usage(const Arena& ar, const std::vector<int>& idxs, std::unordered_map<std::string, NetUsage>& usage) {
    usage.clear();
    for (int i : idxs) {
        const auto& g = ar.nodes[i];
        usage[g.src_net].as_src.push_back(i);
        usage[g.dst_net].as_dst.push_back(i);
    }
}

static int new_leaf(Arena& ar, TransType tp, const Transistor& t) {
    MergeGroup g; g.kind = MergeKind::Leaf; g.type = tp; g.level = 0; g.src_net = t.s; g.dst_net = t.d; g.gate_nets.insert(t.g); g.trans_ids.push_back(t.id);
    return ar.add(g);
}

static int merge_series_chain(Arena& ar, const std::vector<int>& chain, int level) {
    if (chain.size() == 1) return chain[0];
    const auto& A = ar.nodes[chain.front()]; const auto& B = ar.nodes[chain.back()];
    MergeGroup g; g.kind = MergeKind::Series; g.type = A.type; g.level = level; g.src_net = A.src_net; g.dst_net = B.dst_net;
    for (int i : chain) {
        const auto& n = ar.nodes[i];
        g.gate_nets.insert(n.gate_nets.begin(), n.gate_nets.end());
        g.children.push_back(i);
        g.trans_ids.insert(g.trans_ids.end(), n.trans_ids.begin(), n.trans_ids.end());
    }
    return ar.add(g);
}

static int merge_parallel_set(Arena& ar, const std::vector<int>& items, int level) {
    if (items.empty()) return -1;
    const auto& base = ar.nodes[items[0]];
    MergeGroup g; g.kind = MergeKind::Parallel; g.type = base.type; g.level = level; g.src_net = base.src_net; g.dst_net = base.dst_net;
    for (int i : items) {
        const auto& n = ar.nodes[i];
        g.gate_nets.insert(n.gate_nets.begin(), n.gate_nets.end());
        g.children.push_back(i);
        g.trans_ids.insert(g.trans_ids.end(), n.trans_ids.begin(), n.trans_ids.end());
    }
    return ar.add(g);
}

static int merge_parallel_all(Arena& ar, const std::vector<int>& items, int level) {
    if (items.empty()) return -1;
    const auto& base = ar.nodes[items[0]];
    MergeGroup g; g.kind = MergeKind::Parallel; g.type = base.type; g.level = level; g.src_net = base.src_net; g.dst_net = base.dst_net;
    for (int i : items) {
        const auto& n = ar.nodes[i];
        g.gate_nets.insert(n.gate_nets.begin(), n.gate_nets.end());
        g.children.push_back(i);
        g.trans_ids.insert(g.trans_ids.end(), n.trans_ids.begin(), n.trans_ids.end());
    }
    return ar.add(g);
}

struct BuildResult { int top = -1; bool ok = true; std::string reason; };

static BuildResult build_group_tree_subset(Arena& ar, const Netlist& nl, TransType tp, const std::unordered_set<int>& allow) {
    BuildResult br; br.top = -1; br.ok = true; br.reason.clear();
    std::vector<int> idxs;
    for (size_t i=0;i<nl.devices.size();++i) if (allow.count((int)i)) {
        const auto& t = nl.devices[i]; if (t.type == tp) idxs.push_back(new_leaf(ar, tp, t));
    }
    if (idxs.empty()) { br.top = -1; return br; }

    int level = 0; int guard = 0;
    while (idxs.size() > 1 && guard++ < 1000) {
        ++level;
        std::unordered_map<std::string, NetUsage> usage; index_usage(ar, idxs, usage);
        std::map<int,int> next_map; std::map<int,int> prev_map;
        for (const auto& kv : usage) {
            const auto& u = kv.second; const std::string& net = kv.first;
            if (u.as_dst.size() == 1 && u.as_src.size() == 1 && !is_rail(net)) {
                int a = u.as_dst[0]; int b = u.as_src[0];
                const auto& A = ar.nodes[a]; const auto& B = ar.nodes[b];
                if (A.type == B.type && A.dst_net == net && B.src_net == net) {
                    if (!next_map.count(a) && !prev_map.count(b)) { next_map[a] = b; prev_map[b] = a; }
                }
            }
        }
        std::vector<std::vector<int>> chains; std::set<int> in_chain;
        for (int i : idxs) {
            if (prev_map.count(i)) continue;
            std::vector<int> chain; int cur = i; bool has_edge = false; std::set<int> seen;
            while (next_map.count(cur)) {
                if (seen.count(cur)) { br.ok = false; br.reason = "Series cycle detected at level " + std::to_string(level); return br; }
                seen.insert(cur);
                int nxt = next_map[cur]; chain.push_back(cur); cur = nxt; has_edge = true;
            }
            if (has_edge) chain.push_back(cur);
            if (chain.size() >= 2) { chains.push_back(chain); for (int v : chain) in_chain.insert(v); }
        }
        std::map<std::pair<std::string,std::string>, std::vector<int>> buckets;
        for (int i : idxs) {
            const auto& n = ar.nodes[i]; buckets[std::make_pair(n.src_net, n.dst_net)].push_back(i);
        }
        std::vector<std::vector<int>> parallel_sets; std::set<int> in_parallel;
        for (auto& kv : buckets) {
            auto& items = kv.second; if ((int)items.size() >= 2) {
                bool clash = false; for (int i : items) if (in_chain.count(i)) { clash = true; break; }
                if (clash) { br.ok = false; br.reason = "Group appears in both series-chain and parallel at level " + std::to_string(level) + ", bucket (" + kv.first.first + "," + kv.first.second + ")"; return br; }
                bool overlap = false; for (int i : items) if (in_parallel.count(i)) { overlap = true; break; }
                if (!overlap) { parallel_sets.push_back(items); for (int i : items) in_parallel.insert(i); }
            }
        }
        if (chains.empty() && parallel_sets.empty()) {
            int m = merge_parallel_all(ar, idxs, level);
            idxs.clear(); idxs.push_back(m);
            break;
        }
        std::set<int> used = in_chain; used.insert(in_parallel.begin(), in_parallel.end());
        std::vector<int> survivors; for (int i : idxs) if (!used.count(i)) survivors.push_back(i);
        std::vector<int> created;
        for (auto& ch : chains) created.push_back(merge_series_chain(ar, ch, level));
        for (auto& setv : parallel_sets) created.push_back(merge_parallel_set(ar, setv, level));
        idxs.clear(); idxs.insert(idxs.end(), survivors.begin(), survivors.end()); idxs.insert(idxs.end(), created.begin(), created.end());
    }
    if (idxs.empty()) { br.ok = false; br.reason = "No groups built"; return br; }
    br.top = idxs[0]; return br;
}

static void collect_leaves_ids(const Arena& ar, int idx, std::vector<std::string>& out_ids) {
    const auto& g = ar.nodes[idx];
    if (g.kind == MergeKind::Leaf) { out_ids.insert(out_ids.end(), g.trans_ids.begin(), g.trans_ids.end()); return; }
    for (int c : g.children) collect_leaves_ids(ar, c, out_ids);
}

static int count_gate_matches(const std::vector<std::string>& a_ids, const std::vector<std::string>& b_ids,
                              const std::unordered_map<std::string, const Transistor*>& idx) {
    int n = std::min((int)a_ids.size(), (int)b_ids.size()); int score = 0;
    for (int i=0;i<n;++i) {
        auto ita = idx.find(a_ids[i]); auto itb = idx.find(b_ids[i]);
        if (ita!=idx.end() && itb!=idx.end() && ita->second->g == itb->second->g) ++score;
    }
    return score;
}

static void add_dummy(bool is_p, const std::vector<std::string>& ref_ids, const std::unordered_map<std::string, const Transistor*>& by_id,
                      PairMos& out) {
    out.id.clear(); out.is_dummy = true; out.type = (is_p?TransType::PMOS:TransType::NMOS);
    if (ref_ids.empty()) { out.g = out.s = out.d = out.b = "DUMMY"; return; }
    const auto* t = by_id.at(ref_ids[0]);
    out.g = t->g; out.s = t->s; out.d = t->d; out.b = t->b;
}

static PairGroup align_and_pair(int level,
                                const std::vector<std::string>& p_ids,
                                const std::vector<std::string>& n_ids,
                                const std::unordered_map<std::string, const Transistor*>& by_id,
                                bool p_fix_order,
                                bool n_fix_order) {
    PairGroup pg; pg.level = level;
    auto p_norm = p_ids; auto p_rev = p_ids; std::reverse(p_rev.begin(), p_rev.end());
    auto n_norm = n_ids; auto n_rev = n_ids; std::reverse(n_rev.begin(), n_rev.end());
    auto score_variant = [&](const std::vector<std::string>& a, const std::vector<std::string>& b){ return count_gate_matches(a,b,by_id); };

    int best_score = -1; bool use_p_rev = false; bool use_n_rev = false; bool pad_left = false;
    for (int pr=0; pr<=(p_fix_order?0:1); ++pr) {
        for (int nr=0; nr<=(n_fix_order?0:1); ++nr) {
            const auto& A = (pr? p_rev : p_norm);
            const auto& B = (nr? n_rev : n_norm);
            int m = (int)std::max(A.size(), B.size());
            // left pad
            {
                std::vector<std::string> ap=A, bp=B;
                if ((int)ap.size()<m) ap.insert(ap.begin(), m-(int)ap.size(), std::string());
                if ((int)bp.size()<m) bp.insert(bp.begin(), m-(int)bp.size(), std::string());
                int s = score_variant(ap,bp);
                if (s>best_score) { best_score=s; use_p_rev=pr; use_n_rev=nr; pad_left=true; }
            }
            // right pad
            {
                std::vector<std::string> ap=A, bp=B;
                if ((int)ap.size()<m) ap.insert(ap.end(), m-(int)ap.size(), std::string());
                if ((int)bp.size()<m) bp.insert(bp.end(), m-(int)bp.size(), std::string());
                int s = score_variant(ap,bp);
                if (s>best_score) { best_score=s; use_p_rev=pr; use_n_rev=nr; pad_left=false; }
            }
        }
    }
    const auto& A = (use_p_rev? p_rev : p_norm);
    const auto& B = (use_n_rev? n_rev : n_norm);
    int M = (int)std::max(A.size(), B.size());
    std::vector<std::string> pa=A, na=B;
    if (pad_left) {
        if ((int)pa.size()<M) pa.insert(pa.begin(), M-(int)pa.size(), std::string());
        if ((int)na.size()<M) na.insert(na.begin(), M-(int)na.size(), std::string());
    } else {
        if ((int)pa.size()<M) pa.insert(pa.end(), M-(int)pa.size(), std::string());
        if ((int)na.size()<M) na.insert(na.end(), M-(int)na.size(), std::string());
    }
    for (int i=0;i<M;++i) {
        Pair pair; pair.x=0; pair.y=0;
        if (!pa[i].empty()) { const auto* t = by_id.at(pa[i]); pair.pmos = {t->id, TransType::PMOS, t->g, t->s, t->d, t->b, false}; }
        else { std::vector<std::string> ref; if (pad_left){ for (int k=i+1;k<M;++k) if (!pa[k].empty()){ ref.push_back(pa[k]); break; }} else { for (int k=i-1;k>=0;--k) if (!pa[k].empty()){ ref.push_back(pa[k]); break; }} add_dummy(true, ref, by_id, pair.pmos); }
        if (!na[i].empty()) { const auto* t = by_id.at(na[i]); pair.nmos = {t->id, TransType::NMOS, t->g, t->s, t->d, t->b, false}; }
        else { std::vector<std::string> ref; if (pad_left){ for (int k=i+1;k<M;++k) if (!na[k].empty()){ ref.push_back(na[k]); break; } } else { for (int k=i-1;k>=0;--k) if (!na[k].empty()){ ref.push_back(na[k]); break; } } add_dummy(false, ref, by_id, pair.nmos); }
        pg.pairs.push_back(pair);
    }
    return pg;
}

static void collect_all_leaf_mos(const Arena& ar, int idx, TransType tp, const std::unordered_map<std::string, const Transistor*>& by_id, std::vector<PairMos>& out) {
    std::vector<std::string> ids; collect_leaves_ids(ar, idx, ids);
    for (auto& id : ids) { const auto* t = by_id.at(id); out.push_back({t->id, tp, t->g, t->s, t->d, t->b, false}); }
}

// Hungarian algorithm for maximum weight assignment. Pads to square.
static std::vector<int> hungarian_max(const std::vector<std::vector<double>>& w) {
    int n = (int)w.size(); if (n==0) return {};
    double M = 0.0; for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (w[i][j] > M) M = w[i][j];
    if (M <= 0) M = 1.0;
    // Build cost matrix for minimization: cost = M - weight
    std::vector<std::vector<double>> a(n+1, std::vector<double>(n+1, 0.0));
    for (int i=1;i<=n;++i) for (int j=1;j<=n;++j) a[i][j] = M - w[i-1][j-1];
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> u(n+1), v(n+1);
    std::vector<int> p(n+1), way(n+1);
    for (int i=1;i<=n;++i) {
        p[0] = i; int j0 = 0; std::vector<double> minv(n+1, INF); std::vector<char> used(n+1, false);
        do {
            used[j0] = true; int i0 = p[j0]; double delta = INF; int j1 = 0;
            for (int j=1;j<=n;++j) if (!used[j]) {
                double cur = a[i0][j] - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j=0;j<=n;++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else { minv[j] -= delta; }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0]; p[j0] = p[j1]; j0 = j1;
        } while (j0);
    }
    std::vector<int> ans(n, -1); // row i -> column ans[i]
    for (int j=1;j<=n;++j) if (p[j] > 0) ans[p[j]-1] = j-1;
    return ans;
}

static bool rec_match_groups(const Arena& P, int p_idx, const Arena& N, int n_idx,
                             const std::unordered_map<std::string, const Transistor*>& by_id,
                             std::vector<PairGroup>& out_groups,
                             std::vector<PairMos>& discrete) {
    const auto& PG = P.nodes[p_idx];
    const auto& NG = N.nodes[n_idx];
    if (PG.level != NG.level) return false;
    if (PG.level == 1 && NG.level == 1) {
        std::vector<std::string> p_ids; collect_leaves_ids(P, p_idx, p_ids);
        std::vector<std::string> n_ids; collect_leaves_ids(N, n_idx, n_ids);
        bool p_fix = (PG.kind == MergeKind::Series);
        bool n_fix = (NG.kind == MergeKind::Series);
        PairGroup pg = align_and_pair(PG.level, p_ids, n_ids, by_id, p_fix, n_fix);
        out_groups.push_back(std::move(pg));
        return true;
    }
    if (PG.level == 0 && NG.level == 0) {
        std::vector<std::string> p_ids; collect_leaves_ids(P, p_idx, p_ids);
        std::vector<std::string> n_ids; collect_leaves_ids(N, n_idx, n_ids);
        if (p_ids.empty() || n_ids.empty()) return false;
        PairGroup g = align_and_pair(0, {p_ids[0]}, {n_ids[0]}, by_id, true, true);
        out_groups.push_back(std::move(g));
        return true;
    }
    // Children bipartite matching via Hungarian (max weight); allow leftover leaves to discrete
    auto pch = PG.children; auto nch = NG.children;
    int L = (int)pch.size(), R = (int)nch.size();
    int n = std::max(L, R);
    std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));
    auto gate_overlap = [&](const MergeGroup& A, const MergeGroup& B){ int inter=0; for (auto& g : A.gate_nets) if (B.gate_nets.count(g)) ++inter; int denom=(int)std::max(A.gate_nets.size(), B.gate_nets.size()); return denom? (double)inter/denom : 0.0; };
    for (int i=0;i<L;++i) {
        for (int j=0;j<R;++j) {
            const auto& A = P.nodes[pch[i]]; const auto& B = N.nodes[nch[j]];
            if (A.level != B.level) { W[i][j] = 0.0; continue; }
            double s = 0.0; s += (A.kind==B.kind? 1.0 : -0.2); s += gate_overlap(A,B); s -= 0.02 * std::abs((int)A.trans_ids.size() - (int)B.trans_ids.size());
            if (s < 0) s = 0; W[i][j] = s;
        }
    }
    // pad rows/cols if needed (already zero)
    auto assign = hungarian_max(W);
    std::vector<bool> usedP(L,false), usedN(R,false);
    // Process matches
    for (int i=0;i<L;++i) {
        int j = (i < (int)assign.size()? assign[i] : -1);
        if (j<0 || j>=R) continue;
        if (W[i][j] <= 0) continue; // invalid/no-benefit
        int pi = pch[i]; int ni = nch[j];
        std::vector<PairGroup> tmp;
        if (!rec_match_groups(P, pi, N, ni, by_id, tmp, discrete)) return false;
        usedP[i]=true; usedN[j]=true;
        out_groups.insert(out_groups.end(), tmp.begin(), tmp.end());
    }
    // Leftovers: allow only leaves to go discrete
    for (int i=0;i<L;++i) if (!usedP[i]) {
        const auto& A = P.nodes[pch[i]];
        if (A.level == 0) {
            std::vector<std::string> ids; collect_leaves_ids(P, pch[i], ids);
            for (auto& id : ids) { const auto* t = by_id.at(id); discrete.push_back({t->id, TransType::PMOS, t->g, t->s, t->d, t->b, false}); }
        } else return false;
    }
    for (int j=0;j<R;++j) if (!usedN[j]) {
        const auto& B = N.nodes[nch[j]];
        if (B.level == 0) {
            std::vector<std::string> ids; collect_leaves_ids(N, nch[j], ids);
            for (auto& id : ids) { const auto* t = by_id.at(id); discrete.push_back({t->id, TransType::NMOS, t->g, t->s, t->d, t->b, false}); }
        } else return false;
    }
    return true;
}

static bool detect_inverter_pair(const Transistor& p, const Transistor& n, const std::string& vdd, const std::string& vss, std::string& gate_out) {
    if (p.g != n.g) return false;
    auto other = [&](const std::string& a, const std::string& b, const std::string& rail){ return (a==rail)? b : ((b==rail)? a : std::string()); };
    std::string y1 = other(p.s, p.d, vdd); if (y1.empty() || is_rail(y1)) return false;
    std::string y2 = other(n.s, n.d, vss); if (y2.empty() || is_rail(y2)) return false;
    if (y1 != y2) return false;
    gate_out = p.g; return true;
}

static bool is_same_undirected_netpair(const Transistor& a, const Transistor& b) {
    if ((a.s==b.s && a.d==b.d) || (a.s==b.d && a.d==b.s)) return true; return false;
}

StructuralMatchOutput match_structural_pairs(const Netlist& nl, const TechRules& tr) {
    StructuralMatchOutput out;
    std::unordered_map<std::string, const Transistor*> by_id; for (const auto& t : nl.devices) by_id[t.id]=&t;

    // 0) Preprocess: detect inverter and transmission-gate pairs; emit PairGroup and mark used
    std::vector<bool> used(nl.devices.size(), false);
    std::string vdd = nl.rails.size()>0? nl.rails[0] : std::string("VDD");
    std::string vss = nl.rails.size()>1? nl.rails[1] : std::string("VSS");
    std::vector<int> pidx, nidx; for (int i=0;i<(int)nl.devices.size();++i){ if (nl.devices[i].type==TransType::PMOS) pidx.push_back(i); else nidx.push_back(i);} 

    // Inverters
    for (int ip : pidx) if (!used[ip]) {
        for (int in : nidx) if (!used[in]) {
            std::string g;
            if (detect_inverter_pair(nl.devices[ip], nl.devices[in], vdd, vss, g)) {
                PairGroup gpg; gpg.level=0;
                const auto& tp = nl.devices[ip]; const auto& tn = nl.devices[in];
                Pair pr; pr.x=0; pr.y=0; pr.pmos={tp.id, TransType::PMOS, tp.g,tp.s,tp.d,tp.b,false}; pr.nmos={tn.id, TransType::NMOS, tn.g,tn.s,tn.d,tn.b,false};
                gpg.pairs.push_back(pr); out.groups.push_back(std::move(gpg));
                used[ip]=used[in]=true; break;
            }
        }
    }
    // Transmission-gate: PMOS+NMOS share same two non-rail nets (S/D swapped allowed); gate relation not required
    for (int ip : pidx) if (!used[ip]) {
        for (int in : nidx) if (!used[in]) {
            const auto& tp = nl.devices[ip]; const auto& tn = nl.devices[in];
            if (!is_same_undirected_netpair(tp, tn)) continue;
            if (is_rail(tp.s) || is_rail(tp.d)) continue;
            PairGroup gpg; gpg.level=0;
            Pair pr; pr.x=0; pr.y=0; pr.pmos={tp.id, TransType::PMOS, tp.g,tp.s,tp.d,tp.b,false}; pr.nmos={tn.id, TransType::NMOS, tn.g,tn.s,tn.d,tn.b,false};
            gpg.pairs.push_back(pr); out.groups.push_back(std::move(gpg));
            used[ip]=used[in]=true; break;
        }
    }

    // 1) Build D/S connectivity (ignore rails) on remaining devices; split into components
    std::unordered_map<std::string, std::vector<int>> net2devs;
    for (int i=0;i<(int)nl.devices.size();++i) if (!used[i]) {
        const auto& t = nl.devices[i]; if (!is_rail(t.s)) net2devs[t.s].push_back(i); if (!is_rail(t.d)) net2devs[t.d].push_back(i);
    }
    std::vector<int> remIdx; for (int i=0;i<(int)nl.devices.size();++i) if (!used[i]) remIdx.push_back(i);
    std::vector<std::vector<int>> comps;
    std::unordered_set<int> vis;
    for (int start : remIdx) if (!vis.count(start)) {
        std::vector<int> comp; std::queue<int>q; q.push(start); vis.insert(start);
        while(!q.empty()){
            int u=q.front(); q.pop(); comp.push_back(u);
            const auto& t = nl.devices[u];
            auto expand = [&](const std::string& net){ auto it=net2devs.find(net); if(it==net2devs.end()) return; for(int v:it->second) if(!vis.count(v)){ vis.insert(v); q.push(v);} };
            if (!is_rail(t.s)) expand(t.s); if (!is_rail(t.d)) expand(t.d);
        }
        comps.push_back(std::move(comp));
    }

    // 2) For each component, build P/N group trees and recursively match
    for (auto& comp : comps) {
        std::unordered_set<int> allow(comp.begin(), comp.end());
        Arena P, N; auto buildP = build_group_tree_subset(P, nl, TransType::PMOS, allow); auto buildN = build_group_tree_subset(N, nl, TransType::NMOS, allow);
        if (!buildP.ok || !buildN.ok || buildP.top<0 || buildN.top<0) {
            out.used_fallback = true; if (!buildP.ok) out.fail_reason=buildP.reason; else if (!buildN.ok) out.fail_reason=buildN.reason; else out.fail_reason="Empty P or N in component";
            if (buildP.top>=0) collect_all_leaf_mos(P, buildP.top, TransType::PMOS, by_id, out.discrete);
            if (buildN.top>=0) collect_all_leaf_mos(N, buildN.top, TransType::NMOS, by_id, out.discrete);
            continue;
        }
        std::vector<PairGroup> gvec;
        if (!rec_match_groups(P, buildP.top, N, buildN.top, by_id, gvec, out.discrete)) {
            out.used_fallback = true; out.fail_reason = "Iterative match failed in component";
            collect_all_leaf_mos(P, buildP.top, TransType::PMOS, by_id, out.discrete);
            collect_all_leaf_mos(N, buildN.top, TransType::NMOS, by_id, out.discrete);
            continue;
        }
        out.groups.insert(out.groups.end(), gvec.begin(), gvec.end());
    }

    return out;
}

MatchResult pairs_to_match_result(const StructuralMatchOutput& out) {
    MatchResult mr;
    auto emit_pair = [&](const Pair& pair){
        if (!pair.pmos.is_dummy && !pair.pmos.id.empty()) {
            FoldedDevice fd; fd.base_id = pair.pmos.id; fd.fold_index = 0; fd.type = TransType::PMOS; fd.w_um = 0; fd.g=pair.pmos.g; fd.s=pair.pmos.s; fd.d=pair.pmos.d; fd.b=pair.pmos.b; mr.p_row.push_back(fd);
        }
        if (!pair.nmos.is_dummy && !pair.nmos.id.empty()) {
            FoldedDevice fd; fd.base_id = pair.nmos.id; fd.fold_index = 0; fd.type = TransType::NMOS; fd.w_um = 0; fd.g=pair.nmos.g; fd.s=pair.nmos.s; fd.d=pair.nmos.d; fd.b=pair.nmos.b; mr.n_row.push_back(fd);
        }
    };
    for (const auto& g : out.groups) for (const auto& p : g.pairs) emit_pair(p);
    for (const auto& m : out.discrete) {
        if (m.id.empty()) continue;
        FoldedDevice fd; fd.base_id = m.id; fd.fold_index = 0; fd.type = m.type; fd.w_um = 0; fd.g=m.g; fd.s=m.s; fd.d=m.d; fd.b=m.b;
        if (m.type == TransType::PMOS) mr.p_row.push_back(fd); else mr.n_row.push_back(fd);
    }
    return mr;
}

MatchResult match_structural(const Netlist& nl, const TechRules& tr) {
    auto out = match_structural_pairs(nl, tr);
    return pairs_to_match_result(out);
}

} // namespace stdcell
