#pragma once

#include "stdcell/core.hpp"

namespace stdcell {

// Simple Manhattan router honoring M0 horizontal and M1 vertical preference
Routed route_simple(const MatchResult& match, Placement placement, const TechRules& tr, const Netlist& nl);

} // namespace stdcell
