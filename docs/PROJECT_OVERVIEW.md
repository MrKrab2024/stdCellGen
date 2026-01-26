# Project Overview: StdCell Structural Matching & Folding

This document summarizes the goals, data model, pipeline, and I/O of the stdCell project so a new agent can continue development without re-deriving context.

## Scope and Goals
- Purpose: deterministic modeling of standard-cell transistor structure, pairing, and folding; no routing/DRC/optimization.
- Inputs: oriented SPICE netlist per library; technology options (multi-rail aliases, folding params).
- Outputs: reproducible logs describing matched pair groups, discrete MOS, and enumerated folding placements.
- Determinism: all passes are stable and order-independent given the same inputs (no random seeds).

## Repository & Key Paths
- Config: `data/asap7.config`
  - `library_spice`: points to `out/oriented_asap7.sp` (or another oriented netlist).
  - `vdd_net`, `vss_net`: comma-separated alias sets (multi-rail), e.g. `vdd,vddG,vddR` and `vss,vssG,vssD`.
- Matching core:
  - `include/stdcell/structural_matcher.hpp` (Pair, PairGroup, output types)
  - `src/matcher/structural_matcher.cpp` (iterative structural matching pipeline)
- Folding core:
  - `include/stdcell/folding.hpp` (Orient, finger placement, enumerator API)
  - `src/placer/folding.cpp` (enumeration implementation)
- Apps/tools:
  - `apps/stdcellgen.cpp` (runs matching; emits pairs bundle)
  - `apps/fold_dump.cpp` (reads bundle; emits folding logs)
- Main outputs:
  - Orientation: `out/orient_report.json`, `out/orient_log.txt`, oriented netlist `out/oriented_asap7.sp`
  - Matching debug: `out/analysis/debug/<Cell>.debug.json`
  - Matching summary: `out/analysis/asap7_match_summary.csv`
  - Per-cell pairs (human-readable): `out/analysis/pairs/<Cell>.pairs.json`
  - Per-library pairs bundle (machine): `out/oriented_asap7/pairs.bundle.json`
  - Folding logs: `out/analysis/folding/<Cell>.fold.json`

## Data Model (Essentials)
- Netlist devices: canonicalized MOS with fields `{id, type in {PMOS,NMOS}, g,s,d,b nets, nfin}`.
- DeviceRef: `{ id, is_dummy }` ? id-only reference used throughout pairing to keep structures small and stable.
- Pair: `{ p: DeviceRef, n: DeviceRef }` ? one PMOS and one NMOS, by id.
- PairGroup: ordered list of `Pair`; order is meaningful for later placement.
- StructuralMatchOutput: `{ groups: vector<PairGroup>, discrete: vector<DeviceRef> }`.
- TGGroup (extraction-time composite): recorded but its devices are excluded from normal matching and placed as an atomic pair.

## Matching Pipeline (Per Cell)
1) Build Inverter Complement Map (global per cell)
   - Detect simple inverter pairs; create a directed graph `inv_adj` mapping net -> inv(net).
   - Nets A and B are complementary if there exists an odd-length path via `inv_adj` (are_complementary).
2) Partition Into Components
   - Use connectivity (excluding rails) to find independent subgraphs to process.
3) Transmission Gate (TG) Pre-Pass (first step inside each component)
   - Definition: one PMOS and one NMOS whose drain nets are the same AND whose gate nets are complementary (per map above).
   - Action: mark both devices as consumed; create a `Pair` and place it in a dedicated `PairGroup` (TG as an atomic block). Exclude them from normal matching.
4) Iterative Structural Matching (unchanged logic)
   - Build series/parallel groups; merge parallel groups with identical ends; propagate sequences from rails to shared nets.
   - Do not conflate TG with ambiguous: TG is a recognized structure and recorded as a Pair; ambiguous devices are kept in `ambiguous_devices` and/or `discrete`.
5) Output
   - Per-component failures emit a snapshot with exact `groups` (by pairs) and `discrete` MOS.
   - Success emits per-cell pairs JSON and contributes to the library bundle.

Notes
- Multi-rail handling uses `vdd_net`/`vss_net` alias sets from config. Legacy `nl.rails[0/1]` is deprecated.
- Determinism is enforced by stable sorting when emitting pairs / groups.

## Debug & Logs (Shapes and Failures)
- `out/analysis/debug/<Cell>.debug.json`
  - If failure: includes `fail_reason`, `groups` (hierarchy levels with `pairs` by ids), and `discrete` (remaining MOS).
  - Use this to inspect why a component became asymmetric or empty on one side (e.g., "Empty P or N in component"). This often happens if pre-pass filters remove all devices from one polarity inside a component.
- `out/analysis/pairs/<Cell>.pairs.json`
  - Flattened list of `Pair` ids per cell for quick inspection.
- `out/oriented_asap7/pairs.bundle.json`
  - Machine-friendly aggregation of all cells' pairs; input to folding enumerator for placement studies.

## Folding Model (Deterministic Enumeration)
- Minimal-finger principle
  - Split MOS into fingers with `nfin <= max_fins_per_row`; use the unique minimal split (e.g., 7 with max 3 -> [3,3,1]).
  - Each finger named `{orig}_{idx}`; no alternative split vectors are enumerated.
- Orientations
  - Only strict alternation within each MOS sequence: either [SD,DS,SD,...] or [DS,SD,DS,...].
- Physical Rows
  - Use logical rows r in [0..num_logical_rows-1]. Map to physical y as:
    - r even: NMOS y = 2*r, PMOS y = 2*r + 1
    - r odd : PMOS y = 2*r, NMOS y = 2*r + 1
- Slot Model (pairing before placement)
  - Let p = #pMOS fingers, n = #nMOS fingers, M = max(p,n).
  - Create M slots. First min(p,n) slots are two-sided; the remaining are single-sided (one side missing). No separate "air" object.
- Folding Shapes
  - 1-row: single row of length M.
  - 2-row: parameters (l0,l1,dx) with l0>=1, l1>=1, l0+l1=M, and dx in [-(l1-1), (l0-1)]. Row1 is horizontally offset by `dx` slots.
  - Constraints: single-sided slots must live only at the endpoints and be contiguous (prefix or suffix in the linearized slot order).
  - Placement: fill row0 left->right, then row1 left->right; within a slot, p and n share the same column (same x), and must be in the same logical row (adjacent physical y). No SDB/DDB spacing inserted.
- Output
  - `out/analysis/folding/<Cell>.fold.json` contains, for each pair (by ids), `variants` count and a few `examples` of placements: each entry has `{id, idx, type, nfin, x, y, orient}` per finger.

## Typical Workflow
1) Set `data/asap7.config`:
   - `library_spice = out/oriented_asap7.sp`
   - `vdd_net = VDD,VDDG,VDDr` (example)
   - `vss_net = VSS,VSSG,VSSd` (example)
2) Build and run matching
   - Example: `cmake -B build -S . -DCMAKE_BUILD_TYPE=Release`
   - `cmake --build build --config Release`
   - `build/Release/stdcellgen.exe -c data/asap7.config`
   - Results under `out/analysis/pairs/`, `out/oriented_asap7/pairs.bundle.json`, and `out/analysis/debug/`.
3) Enumerate folding
   - `build/Release/fold_dump.exe out/oriented_asap7/pairs.bundle.json`
   - Inspect `out/analysis/folding/*.fold.json`.

## Troubleshooting
- Transmission gates not extracted in a component
  - Ensure complement map is built once per cell, and TG pre-pass runs as the first step inside each component.
  - A TG requires same drain net and complementary gates; if drains differ or complement cannot be proven, it remains in normal matching.
- "Empty P or N in component"
  - Pre-pass or filtering removed all devices of one polarity in a component; that component cannot form pairs. This is expected for some control-only islands.
- AO22x1 alignment questions
  - Series groups are merged only when their terminal ends match; shared internal nets (e.g., `net18`) do not force merging into a single series chain if not topologically justified.

## Design Principles (Keep in Mind)
- Keep `Pair` and `DeviceRef` id-only. When full device info is needed (for printing nets or types), resolve by id from the netlist.
- Do not special-case clock names; always infer complement via the inverter graph.
- Make generation, not filtering, enforce constraints (e.g., single-sided at endpoints, same logical row for n/p).
- Preserve existing matching semantics; TG pre-pass must not change non-TG behavior.

## Next Steps
- Maintain the pairs bundle as the stable interface between matching and folding.
- Keep regression: `out/analysis/asap7_match_summary.csv` should remain stable when refactoring.
- If changing `structural_matcher.cpp`, remove legacy fields (`PairMos`, `Pair.x/y`, `DeviceRef.{g,s,d,type}`) and look up nets/types by id only when emitting.

