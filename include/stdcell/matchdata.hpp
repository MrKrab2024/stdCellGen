#pragma once

#include "stdcell/core.hpp"
#include <string>
#include <vector>

namespace stdcell {

struct MatchConfig {
    TechRules tr;
    std::string vdd;
    std::string vss;
    double default_w_per_fin_um = 0.0; // 0 means derive if possible
};

struct MatchModelKey {
    TransType type;
    double l_um;
    double w_per_fin_um;
    std::string g, s, d, b;
    std::string key() const; // stable string key
};

struct MatchTrans {
    std::string name;   // first device id for this model group
    MatchModelKey model;
    int nfin = 0;       // accumulated fins
    int fingers = 0;    // accumulated fingers
    int x = 0;          // integer coords (to be filled by placer)
    int y = 0;
};

struct MatchData {
    std::string cell_name;
    std::vector<MatchTrans> trans;
};

MatchConfig build_match_config(const TechRules& tr, const std::vector<std::string>& rails_hint);
MatchData build_match_data(const Netlist& nl, const MatchConfig& mc);

// Convert MatchData to legacy MatchResult by expanding width = w_per_fin * nfin
MatchResult build_match_result(const MatchData& md, const MatchConfig& mc);

} // namespace stdcell
