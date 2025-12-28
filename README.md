# stdcellgen (demo)

C++ standard-cell generator scaffold inspired by AutoCellGen. Matching, placement, and routing are split into modules.

- Matching: folding and order alignment (src/matcher)
- Placement: heuristic or Python RL adapter (src/placer, python/rl_placer)
- Routing: minimal M0 horizontal rails and M1 vertical trunks (src/router)

Build

- Requirements: CMake (>=3.15) and a C++17 compiler
- Windows: run `scripts/build.ps1`

Run demo

- `scripts/run_demo.ps1` will generate `out/layout.json` and `out/layout.txt`
- Use `--use-rl` to invoke the Python stub RL placer

Tech constraints

- Tracks per row configurable (default 9)
- M0 horizontal only; M1 vertical only

License: MIT
