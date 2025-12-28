#pragma once

#include "stdcell/core.hpp"

namespace stdcell {

MatchResult match_and_fold(const Netlist& nl, const TechRules& tr);

} // namespace stdcell
