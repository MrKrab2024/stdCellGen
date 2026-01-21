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

// Unified device reference used by matcher outputs (no geometry, no net duplication)
struct DeviceRef {
    std::string id;       // refers to Netlist.devices[*].id
    bool is_dummy = false; // true when padding pair columns
};

struct Pair {
    DeviceRef p; // PMOS id (or dummy)
    DeviceRef n; // NMOS id (or dummy)
};

struct PairGroup {
    int level = 0;                 // matched level
    std::vector<Pair> pairs;       // ordered pairs
};

// Debug dump for failed structural matching
struct GroupDump {
    int level = 0;
    MergeKind kind = MergeKind::Leaf;
    std::string src_net;
    std::string dst_net;
    std::vector<std::string> trans_ids; // leaf devices inside this group
    std::vector<GroupDump> children;    // recursive children
};

struct FailureDump {
    std::string reason;  // reason string at failure point
    int at_level = -1;   // level where failure happened
    GroupDump left;      // PMOS-side snapshot
    GroupDump right;     // NMOS-side snapshot
};

struct StructuralMatchOutput {
    std::vector<PairGroup> groups;     // matched groups
    std::vector<DeviceRef> discrete;   // unmatched level-0 devices (ids only)
    bool used_fallback = false;        // true if structural matching failed at non-leaf level
    std::string fail_reason;           // description when used_fallback

    std::vector<FailureDump> failures; // detailed failure snapshots
};

// Build structural groups and produce pair groups and discrete devices per spec.
StructuralMatchOutput match_structural_pairs(const Netlist& nl, const TechRules& tr);

// Convert pair groups into legacy MatchResult ordering for current placer/router pipeline.
MatchResult pairs_to_match_result(const StructuralMatchOutput& out, const Netlist& nl);

// Legacy API kept: returns MatchResult by delegating to pairs.
MatchResult match_structural(const Netlist& nl, const TechRules& tr);

} // namespace stdcell