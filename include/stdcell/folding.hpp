#pragma once

#include "stdcell/core.hpp"
#include <string>
#include <vector>

namespace stdcell {

// Orientation of a finger: left->right is Source->Drain or Drain->Source
enum class Orient { SD, DS };

// One finger placement result (for either PMOS or NMOS)
struct FingerPlacement {
    std::string base_id;   // original transistor id (e.g., M1)
    int finger_index = 0;  // index within the split sequence (0-based)
    int nfin = 0;          // fins in this finger
    int y = 0;           // physical track row (0,1,2,...)
    int x = 0;             // x of left S/D terminal in grid columns (0,2,4,...)
    Orient orient = Orient::SD;
    TransType type = TransType::PMOS; // keeps PMOS/NMOS
};

// A candidate variant for a pair (one PMOS + one NMOS)
struct PairVariant {
    std::vector<FingerPlacement> fingers; // combined PMOS and NMOS
};

// Compute effective fins from a device and tech rules; fallback to width-based estimate if needed.
int effective_nfin(const Transistor& t, const TechRules& tr);

// Enumerate all valid folding variants for a pair of transistors.
// - n: NMOS device (type must be NMOS)
// - p: PMOS device (type must be PMOS)
// - max_fins_per_row: each finger must have nfin <= this value
// - num_logical_rows: logical rows; each maps to two physical rows (y) per spec
std::vector<PairVariant> enumerate_pair_variants(
    const Transistor& n,
    const Transistor& p,
    int max_fins_per_row,
    int num_logical_rows,
    const TechRules& tr
);

} // namespace stdcell

