#include "stdcell/matchdata.hpp"
#include <sstream>
#include <unordered_map>
#include <iomanip>
#include <algorithm>

namespace stdcell {

static std::string str_type(TransType t) { return t == TransType::PMOS ? "P" : "N"; }

std::string MatchModelKey::key() const {
    std::ostringstream os; os << str_type(type) << "|L=" << std::fixed << std::setprecision(6) << l_um
        << "|WPF=" << w_per_fin_um << "|G=" << g << "|S=" << s << "|D=" << d << "|B=" << b;
    return os.str();
}

MatchConfig build_match_config(const TechRules& tr, const std::vector<std::string>& rails_hint) {
    MatchConfig mc; mc.tr = tr; mc.default_w_per_fin_um = tr.default_w_per_fin_um;
    if (!rails_hint.empty()) mc.vdd = rails_hint[0];
    if (rails_hint.size() > 1) mc.vss = rails_hint[1];
    return mc;
}

static int safe_round_to_int(double v) { if (v <= 0) return 0; return (int)llround(v); }

MatchData build_match_data(const Netlist& nl, const MatchConfig& mc) {
    MatchData md; md.cell_name = nl.cell_name;
    std::unordered_map<std::string, size_t> idx;
    for (const auto& t : nl.devices) {
        int nfin = t.nfin;
        double wpf = 0.0;
        if (nfin > 0) {
            wpf = (t.w_um > 0.0) ? (t.w_um / nfin) : mc.default_w_per_fin_um;
        } else if (mc.default_w_per_fin_um > 0.0 && t.w_um > 0.0) {
            nfin = std::max(1, safe_round_to_int(t.w_um / mc.default_w_per_fin_um));
            wpf = t.w_um / nfin;
        } else {
            nfin = (t.w_um > 0.0 ? 1 : 0);
            wpf = (t.w_um > 0.0 ? t.w_um : mc.tr.min_w_nmos_um);
        }
        int fingers = (t.fingers > 0 ? t.fingers : 1);

        MatchModelKey mk{t.type, t.l_um, wpf, t.g, t.s, t.d, t.b};
        std::string k = mk.key();
        auto it = idx.find(k);
        if (it == idx.end()) {
            MatchTrans mt; mt.name = t.id; mt.model = mk; mt.nfin = nfin; mt.fingers = fingers; mt.x = 0; mt.y = 0;
            idx.emplace(k, md.trans.size());
            md.trans.push_back(mt);
        } else {
            auto& mt = md.trans[it->second];
            mt.nfin += nfin;
            mt.fingers += fingers;
        }
    }
    return md;
}

MatchResult build_match_result(const MatchData& md, const MatchConfig& mc) {
    MatchResult mr;
    for (const auto& mt : md.trans) {
        double width = mt.model.w_per_fin_um * std::max(1, mt.nfin);
        FoldedDevice fd; fd.base_id = mt.name; fd.fold_index = 0; fd.type = mt.model.type; fd.w_um = width; fd.g = mt.model.g; fd.s = mt.model.s; fd.d = mt.model.d; fd.b = mt.model.b;
        if (fd.type == TransType::PMOS) mr.p_row.push_back(fd); else mr.n_row.push_back(fd);
    }
    auto by_gate = [](const FoldedDevice& a, const FoldedDevice& b){ return a.g < b.g; };
    std::stable_sort(mr.p_row.begin(), mr.p_row.end(), by_gate);
    std::stable_sort(mr.n_row.begin(), mr.n_row.end(), by_gate);
    return mr;
}

} // namespace stdcell
