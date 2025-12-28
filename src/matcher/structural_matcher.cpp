#include "stdcell/structural_matcher.hpp"

#include <unordered_map>
#include <map>
#include <algorithm>
#include <queue>

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

static int merge_series(Arena& ar, int a, int b) {
    const auto& A = ar.nodes[a]; const auto& B = ar.nodes[b];
    MergeGroup g; g.kind = MergeKind::Series; g.type = A.type; g.level = std::max(A.level, B.level) + 1; g.src_net = A.src_net; g.dst_net = B.dst_net;
    g.gate_nets = A.gate_nets; g.gate_nets.insert(B.gate_nets.begin(), B.gate_nets.end());
    g.children = { a, b };
    g.trans_ids = A.trans_ids; g.trans_ids.insert(g.trans_ids.end(), B.trans_ids.begin(), B.trans_ids.end());
    return ar.add(g);
}

static int merge_parallel(Arena& ar, const std::vector<int>& items) {
    if (items.empty()) return -1;
    const auto& base = ar.nodes[items[0]];
    MergeGroup g; g.kind = MergeKind::Parallel; g.type = base.type; g.level = 0; g.src_net = base.src_net; g.dst_net = base.dst_net;
    for (int i : items) {
        const auto& n = ar.nodes[i];
        g.level = std::max(g.level, n.level + 1);
        g.gate_nets.insert(n.gate_nets.begin(), n.gate_nets.end());
        g.children.push_back(i);
        g.trans_ids.insert(g.trans_ids.end(), n.trans_ids.begin(), n.trans_ids.end());
    }
    return ar.add(g);
}

static std::vector<int> build_initial_groups(Arena& ar, const Netlist& nl, TransType tp) {
    std::vector<int> idxs;
    for (const auto& t : nl.devices) if (t.type == tp) idxs.push_back(new_leaf(ar, tp, t));
    return idxs;
}

static bool try_parallel_merges(Arena& ar, std::vector<int>& idxs) {
    std::map<std::pair<std::string,std::string>, std::vector<int>> buckets;
    for (int i : idxs) {
        const auto& n = ar.nodes[i];
        buckets[std::make_pair(n.src_net, n.dst_net)].push_back(i);
    }
    bool changed = false; std::vector<int> next;
    for (auto& kv : buckets) {
        auto& items = kv.second;
        if (items.size() >= 2) {
            int m = merge_parallel(ar, items);
            next.push_back(m);
            changed = true;
        } else {
            next.push_back(items[0]);
        }
    }
    if (changed) idxs.swap(next);
    return changed;
}

static bool try_series_merges(Arena& ar, std::vector<int>& idxs) {
    std::unordered_map<std::string, NetUsage> usage; index_usage(ar, idxs, usage);
    std::vector<std::pair<int,int>> pairs;
    for (const auto& kv : usage) {
        const auto& u = kv.second;
        if (u.as_dst.size() == 1 && u.as_src.size() == 1) {
            int a = u.as_dst[0]; int b = u.as_src[0];
            const auto& A = ar.nodes[a]; const auto& B = ar.nodes[b];
            if (A.type == B.type && A.dst_net == B.src_net && !is_rail(A.dst_net)) {
                pairs.emplace_back(a, b);
            }
        }
    }
    if (pairs.empty()) return false;
    std::set<int> used; std::vector<int> survivors; std::vector<int> created;
    for (auto& pr : pairs) {
        if (used.count(pr.first) || used.count(pr.second)) continue;
        used.insert(pr.first); used.insert(pr.second);
        int m = merge_series(ar, pr.first, pr.second);
        created.push_back(m);
    }
    for (int i : idxs) if (!used.count(i)) survivors.push_back(i);
    survivors.insert(survivors.end(), created.begin(), created.end());
    idxs.swap(survivors);
    return !created.empty();
}

static int build_group_tree(Arena& ar, const Netlist& nl, TransType tp) {
    auto idxs = build_initial_groups(ar, nl, tp);
    if (idxs.empty()) return -1;
    bool changed = true;
    for (int iter=0; iter<100 && changed; ++iter) {
        changed = false;
        changed = try_parallel_merges(ar, idxs) || changed;
        changed = try_series_merges(ar, idxs) || changed;
    }
    if (idxs.size() == 1) return idxs[0];
    int top = merge_parallel(ar, idxs);
    return top;
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
    const auto* t = by_id.at(ref_ids[0]); // copy nets from boundary
    out.g = t->g; out.s = t->s; out.d = t->d; out.b = t->b;
}

static PairGroup align_and_pair(MergeKind kind, int level,
                                const std::vector<std::string>& p_ids,
                                const std::vector<std::string>& n_ids,
                                const std::unordered_map<std::string, const Transistor*>& by_id) {
    PairGroup pg; pg.kind = kind; pg.level = level;
    // Try two orientations for N: normal vs reversed
    auto n_norm = n_ids; auto n_rev = n_ids; std::reverse(n_rev.begin(), n_rev.end());
    // Try padding on left or right; we choose based on maximizing gate matches
    auto score_variant = [&](const std::vector<std::string>& a, const std::vector<std::string>& b){ return count_gate_matches(a,b,by_id); };

    int best_score = -1; bool use_rev = false; bool pad_left = false;
    for (int rev=0; rev<2; ++rev) {
        const auto& b = (rev? n_rev : n_norm);
        int m = (int)std::max(p_ids.size(), b.size());
        // pad left
        {
            std::vector<std::string> ap = p_ids, bp = b;
            if (ap.size() < m) ap.insert(ap.begin(), m - ap.size(), std::string());
            if (bp.size() < m) bp.insert(bp.begin(), m - bp.size(), std::string());
            int s = score_variant(ap, bp);
            if (s > best_score) { best_score = s; use_rev = (rev!=0); pad_left = true; }
        }
        // pad right
        {
            std::vector<std::string> ap = p_ids, bp = b;
            if (ap.size() < m) ap.insert(ap.end(), m - ap.size(), std::string());
            if (bp.size() < m) bp.insert(bp.end(), m - bp.size(), std::string());
            int s = score_variant(ap, bp);
            if (s > best_score) { best_score = s; use_rev = (rev!=0); pad_left = false; }
        }
    }
    const auto& bn = (use_rev? n_rev : n_norm);
    int M = (int)std::max(p_ids.size(), bn.size());
    std::vector<std::string> pa = p_ids, na = bn;
    if (pad_left) {
        if (pa.size() < M) pa.insert(pa.begin(), M - pa.size(), std::string());
        if (na.size() < M) na.insert(na.begin(), M - na.size(), std::string());
    } else {
        if (pa.size() < M) pa.insert(pa.end(), M - pa.size(), std::string());
        if (na.size() < M) na.insert(na.end(), M - na.size(), std::string());
    }
    // Build pairs; missing ids are dummies placed only at edges as constructed
    for (int i=0;i<M;++i) {
        Pair pair; pair.x = 0; pair.y = 0;
        if (!pa[i].empty()) {
            const auto* t = by_id.at(pa[i]); pair.pmos = {t->id, TransType::PMOS, t->g, t->s, t->d, t->b, false};
        } else {
            // decide reference: if pad_left, take first real afterwards, else last real before
            std::vector<std::string> ref;
            if (pad_left) {
                for (int k=i+1;k<M;++k) if (!pa[k].empty()) { ref.push_back(pa[k]); break; }
            } else {
                for (int k=i-1;k>=0;--k) if (!pa[k].empty()) { ref.push_back(pa[k]); break; }
            }
            add_dummy(true, ref, by_id, pair.pmos);
        }
        if (!na[i].empty()) {
            const auto* t = by_id.at(na[i]); pair.nmos = {t->id, TransType::NMOS, t->g, t->s, t->d, t->b, false};
        } else {
            std::vector<std::string> ref;
            if (pad_left) {
                for (int k=i+1;k<M;++k) if (!na[k].empty()) { ref.push_back(na[k]); break; }
            } else {
                for (int k=i-1;k>=0;--k) if (!na[k].empty()) { ref.push_back(na[k]); break; }
            }
            add_dummy(false, ref, by_id, pair.nmos);
        }
        pg.pairs.push_back(pair);
    }
    return pg;
}

static void collect_all_leaf_mos(const Arena& ar, int idx, TransType tp, const std::unordered_map<std::string, const Transistor*>& by_id, std::vector<PairMos>& out) {
    std::vector<std::string> ids; collect_leaves_ids(ar, idx, ids);
    for (auto& id : ids) {
        const auto* t = by_id.at(id);
        out.push_back({t->id, tp, t->g, t->s, t->d, t->b, false});
    }
}

StructuralMatchOutput match_structural_pairs(const Netlist& nl, const TechRules& tr) {
    // Build arenas for P and N
    Arena P, N;
    int p_top = build_group_tree(P, nl, TransType::PMOS);
    int n_top = build_group_tree(N, nl, TransType::NMOS);

    StructuralMatchOutput out;
    // map id -> transistor
    std::unordered_map<std::string, const Transistor*> by_id; for (const auto& t : nl.devices) by_id[t.id]=&t;

    if (p_top < 0 || n_top < 0) {
        // fallback: nothing to match, mark all leaf mos as discrete
        out.used_fallback = true;
        if (p_top >= 0) collect_all_leaf_mos(P, p_top, TransType::PMOS, by_id, out.discrete);
        if (n_top >= 0) collect_all_leaf_mos(N, n_top, TransType::NMOS, by_id, out.discrete);
        return out;
    }

    const auto& PG = P.nodes[p_top];
    const auto& NG = N.nodes[n_top];
    if (PG.level != NG.level) {
        // non-leaf mismatch in level: fail per spec -> simple matching stub
        out.used_fallback = true;
        collect_all_leaf_mos(P, p_top, TransType::PMOS, by_id, out.discrete);
        collect_all_leaf_mos(N, n_top, TransType::NMOS, by_id, out.discrete);
        return out;
    }

    // If both level-0: directly pair into one group of size 1
    std::vector<std::string> p_ids; collect_leaves_ids(P, p_top, p_ids);
    std::vector<std::string> n_ids; collect_leaves_ids(N, n_top, n_ids);

    if (PG.level == 0 && NG.level == 0) {
        if (!p_ids.empty() && !n_ids.empty()) {
            PairGroup g = align_and_pair(MergeKind::Leaf, 0, {p_ids[0]}, {n_ids[0]}, by_id);
            out.groups.push_back(std::move(g));
        } else {
            // any lonely leaf goes to discrete
            if (!p_ids.empty()) { const auto* t = by_id.at(p_ids[0]); out.discrete.push_back({t->id, TransType::PMOS, t->g, t->s, t->d, t->b, false}); }
            if (!n_ids.empty()) { const auto* t = by_id.at(n_ids[0]); out.discrete.push_back({t->id, TransType::NMOS, t->g, t->s, t->d, t->b, false}); }
        }
        return out;
    }

    // Otherwise, create a PairGroup at this matched level.
    // For Series: order fixed except flip; For Parallel: free order ? we handle at leaf sequence alignment via reverse and edge padding.
    PairGroup pg = align_and_pair(PG.kind, PG.level, p_ids, n_ids, by_id);
    out.groups.push_back(std::move(pg));
    return out;
}

MatchResult pairs_to_match_result(const StructuralMatchOutput& out) {
    MatchResult mr;
    auto emit_pair = [&](const Pair& pair){
        if (!pair.pmos.is_dummy) {
            FoldedDevice fd; fd.base_id = pair.pmos.id; fd.fold_index = 0; fd.type = TransType::PMOS; fd.w_um = 0; fd.g=pair.pmos.g; fd.s=pair.pmos.s; fd.d=pair.pmos.d; fd.b=pair.pmos.b; mr.p_row.push_back(fd);
        }
        if (!pair.nmos.is_dummy) {
            FoldedDevice fd; fd.base_id = pair.nmos.id; fd.fold_index = 0; fd.type = TransType::NMOS; fd.w_um = 0; fd.g=pair.nmos.g; fd.s=pair.nmos.s; fd.d=pair.nmos.d; fd.b=pair.nmos.b; mr.n_row.push_back(fd);
        }
    };
    for (const auto& g : out.groups) {
        for (const auto& p : g.pairs) emit_pair(p);
    }
    // Append discrete devices at the end by polarity
    for (const auto& m : out.discrete) {
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
