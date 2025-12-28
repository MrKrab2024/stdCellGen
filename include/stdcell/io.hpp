#pragma once

#include "stdcell/core.hpp"
#include <string>
#include <vector>

namespace stdcell {

// Parsing simple key=value config for tech rules
TechRules parse_tech_rules(const std::string& path);

// Parsing simple cell spec file (legacy)
Netlist parse_cell_spec(const std::string& path);

// Parse SPICE library and extract one subckt as Netlist
Netlist parse_spice_subckt(const std::string& lib_path, const std::string& subckt_name,
                           const std::vector<std::string>& rails_hint);

// Parse unified run config (-config)
RunConfig parse_run_config(const std::string& path);

// Apply overrides from RunConfig to TechRules
void apply_overrides(TechRules& tr, const RunConfig& rc);

// Serialize layout and routing to JSON
void write_layout_json(const std::string& path, const Routed& routed);

// Write brief text summary
void write_layout_txt(const std::string& path, const Routed& routed);

// MatchResult JSON IO for step gating
void write_match_json(const std::string& path, const MatchResult& mr);
bool read_match_json(const std::string& path, MatchResult& mr);

// Placement JSON IO (optional for route-only)
void write_placement_json(const std::string& path, const Placement& pl);
bool read_placement_json(const std::string& path, Placement& pl);

} // namespace stdcell
