#include "stdcell/matcher.hpp"

#include <cmath>
#include <algorithm>

namespace stdcell {

static void append_folds(const Transistor& t, const TechRules& tr,
                         std::vector<FoldedDevice>& out) {
    double maxw = std::max(1e-9, tr.max_diff_w_um);
    int folds = 1;
    double seg_w = t.w_um;
    if (tr.fold_enable && t.w_um > maxw) {
        folds = static_cast<int>(std::ceil(t.w_um / maxw));
        seg_w = t.w_um / folds;
    }
    for (int i = 0; i < folds; ++i) {
        FoldedDevice fd;
        fd.base_id = t.id;
        fd.fold_index = i;
        fd.type = t.type;
        fd.w_um = seg_w;
        fd.g = t.g; fd.s = t.s; fd.d = t.d; fd.b = t.b;
        out.push_back(fd);
    }
}

MatchResult match_and_fold(const Netlist& nl, const TechRules& tr) {
    MatchResult mr;
    // Separate and fold PMOS/NMOS; naive order preserves input order.
    for (const auto& t : nl.devices) {
        if (t.type == TransType::PMOS) append_folds(t, tr, mr.p_row);
    }
    for (const auto& t : nl.devices) {
        if (t.type == TransType::NMOS) append_folds(t, tr, mr.n_row);
    }

    // Simple heuristic: align by gate name order (stable sort by gate to encourage diffusion sharing)
    auto by_gate = [](const FoldedDevice& a, const FoldedDevice& b){ return a.g < b.g; };
    std::stable_sort(mr.p_row.begin(), mr.p_row.end(), by_gate);
    std::stable_sort(mr.n_row.begin(), mr.n_row.end(), by_gate);

    return mr;
}

} // namespace stdcell
