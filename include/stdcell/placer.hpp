#pragma once

#include "stdcell/core.hpp"
#include <memory>
#include <string>

namespace stdcell {

struct IPlacer {
    virtual ~IPlacer() = default;
    virtual Placement place(const MatchResult& match, const TechRules& tr, const Netlist& nl) = 0;
};

std::unique_ptr<IPlacer> make_heuristic_placer();

// RL adapter: writes CSV, calls Python module, reads CSV back
std::unique_ptr<IPlacer> make_rl_adapter_placer(const std::string& python_module_dir);

} // namespace stdcell
