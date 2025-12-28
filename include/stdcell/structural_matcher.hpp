#pragma once

#include "stdcell/core.hpp"
#include <string>
#include <vector>
#include <set>

namespace stdcell {

enum class MergeKind { Leaf, Series, Parallel };

struct MergeGroup {
    MergeKind kind = MergeKind::Leaf;
    TransType type = TransType::PMOS;
    int level = 0;
    std::string src_net;
    std::string dst_net;
    std::set<std::string> gate_nets;
    std::vector<std::string> trans_ids; // for leaves or flattened content
    std::vector<int> children;          // indices in arena
};

// Pairing result types
struct PairMos {
    std::string id;     // empty if dummy
    TransType type;     // PMOS/NMOS
    std::string g, s, d, b;
    bool is_dummy = false;
};

struct Pair {
    PairMos pmos;
    PairMos nmos;
    int x = 0; // integer placement coords (to be filled later)
    int y = 0; // y equal within a pair by construction
};

struct PairGroup {
    MergeKind kind = MergeKind::Parallel; // group structural kind
    int level = 0;                        // matched level
    std::vector<Pair> pairs;              // ordered pairs; for Series, order fixed except flip; for Parallel, order free
};

struct StructuralMatchOutput {
    std::vector<PairGroup> groups;   // matched groups
    std::vector<PairMos> discrete;   // unmatched level-0 devices
    bool used_fallback = false;      // true if structural matching failed at non-leaf level
};

// Build structural groups and produce pair groups and discrete devices per spec.
StructuralMatchOutput match_structural_pairs(const Netlist& nl, const TechRules& tr);

// Convenience: convert pair groups into legacy MatchResult ordering for current placer/router pipeline.
MatchResult pairs_to_match_result(const StructuralMatchOutput& out);

// Legacy API kept: returns MatchResult by delegating to pairs.
MatchResult match_structural(const Netlist& nl, const TechRules& tr);

} // namespace stdcell
