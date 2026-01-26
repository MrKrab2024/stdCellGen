#include "stdcell/placer.hpp"
#include "stdcell/matcher.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <chrono>

namespace stdcell {

namespace {

struct HeuristicPlacer : IPlacer {
    Placement place(const MatchResult& match, const TechRules& tr, const Netlist& nl) override {
        Placement pl;
        const double row_h = tr.tracks_per_row * tr.track_pitch_um;
        pl.cell_h_um = tr.rail_width_um * 2.0 + row_h * 2.0;

        // Decide cell width from larger of p/n row device counts
        const size_t np = match.p_row.size();
        const size_t nn = match.n_row.size();
        const size_t cols = std::max(np, nn);
        const double dev_pitch = tr.track_pitch_um; // 1 track per fold for demo
        const double spacing = 0.5 * tr.track_pitch_um;
        const double margin = tr.track_pitch_um; // space for pins near edges
        pl.cell_w_um = margin * 2.0 + cols * dev_pitch + (cols ? (cols - 1) * spacing : 0);

        // Row y positions (centered in their rows)
        const double n_row_y = tr.rail_width_um + row_h * 0.5; // center of NMOS row
        const double p_row_y = tr.rail_width_um + row_h + row_h * 0.5; // center of PMOS row

        // Place devices left-to-right
        auto place_row = [&](const std::vector<FoldedDevice>& row, double cy) {
            double x = margin;
            for (size_t i = 0; i < row.size(); ++i) {
                Rect r;
                r.name = row[i].base_id + std::string("_f") + std::to_string(row[i].fold_index);
                r.x = x; r.y = cy - (row_h * 0.5) * 0.6; // some height inside row
                r.w = dev_pitch; r.h = row_h * 0.6;
                pl.devices.push_back(r);
                x += dev_pitch + spacing;
            }
        };

        place_row(match.n_row, n_row_y);
        place_row(match.p_row, p_row_y);

        // Create pins: for demo, put all input pins on left, outputs on right
        const double pin_w = tr.track_pitch_um * 0.5;
        const double pin_h = tr.track_pitch_um * 1.5;
        for (const auto& pin : nl.pins) {
            Rect p;
            p.name = pin;
            bool is_output = (!pin.empty() && (pin[0] == 'Y' || pin[0] == 'Z' || pin[0] == 'Q'));
            p.x = is_output ? (pl.cell_w_um - margin * 0.5) : (margin * 0.5);
            p.y = tr.rail_width_um + row_h; // mid height between rows
            p.w = pin_w; p.h = pin_h;
            pl.pins.push_back(p);
        }
        // Also add rails as pins for reference
        // VSS rails at bottom
        for (const auto& rname : nl.vss_nets) {
            Rect p; p.name = rname; p.x = pl.cell_w_um * 0.5; p.w = pl.cell_w_um; p.h = tr.rail_width_um; p.y = 0; pl.pins.push_back(p);
        }
        for (const auto& rname : nl.vdd_nets) {
        }
        return pl;
    }
};

static std::string timestamp_tag() {
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    return std::to_string((long long)ms);
}

struct RLAdapterPlacer : IPlacer {
    std::string module_dir;
    explicit RLAdapterPlacer(std::string m) : module_dir(std::move(m)) {}

    Placement place(const MatchResult& match, const TechRules& tr, const Netlist& nl) override {
        // Fallback to heuristic if Python call fails
        HeuristicPlacer fallback;
        Placement base = fallback.place(match, tr, nl);

        // Prepare CSV in out/ folder
        std::string tag = timestamp_tag();
        std::string dev_csv = std::string("out/rl_devices_") + tag + ".csv";
        std::string out_csv = std::string("out/rl_place_") + tag + ".csv";
        {
            std::ofstream d(dev_csv);
            d << "id,fold,type,gate,src,drain,width_um\n";
            auto dump = [&](const std::vector<FoldedDevice>& row){
                for (const auto& fd : row) {
                    d << fd.base_id << "," << fd.fold_index << "," << (fd.type == TransType::PMOS ? "PMOS" : "NMOS")
                      << "," << fd.g << "," << fd.s << "," << fd.d << "," << fd.w_um << "\n";
                }
            };
            dump(match.n_row); dump(match.p_row);
        }

        // Build PYTHONPATH to include module_dir's parent (so rl_placer is importable)
        std::string cmd = "";
#ifdef _WIN32
        cmd = "set PYTHONPATH=%CD%/python&& python -m rl_placer.main --input " + dev_csv + " --output " + out_csv;
#else
        cmd = "PYTHONPATH=./python python3 -m rl_placer.main --input " + dev_csv + " --output " + out_csv;
#endif
        int rc = std::system(cmd.c_str());
        if (rc != 0) {
            return base; // fallback
        }

        // Read output CSV: id,fold,x
        std::ifstream in(out_csv);
        if (!in.good()) return base;
        std::string line; std::getline(in, line); // header
        std::map<std::string, double> xmap; // key: id_f
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            std::string id; std::string fold; std::string xs;
            std::getline(iss, id, ',');
            std::getline(iss, fold, ',');
            std::getline(iss, xs, ',');
            double x = std::atof(xs.c_str());
            xmap[id + std::string("_") + fold] = x;
        }
        // Apply x positions to devices; keep y,w,h from base
        for (auto& r : base.devices) {
            // r.name == baseId_fN
            auto pos = r.name.find("_f");
            if (pos == std::string::npos) continue;
            std::string baseId = r.name.substr(0, pos);
            std::string fold = r.name.substr(pos + 2);
            auto it = xmap.find(baseId + std::string("_") + fold);
            if (it != xmap.end()) r.x = it->second;
        }
        return base;
    }
};

} // namespace

std::unique_ptr<IPlacer> make_heuristic_placer() { return std::make_unique<HeuristicPlacer>(); }

std::unique_ptr<IPlacer> make_rl_adapter_placer(const std::string& python_module_dir) {
    return std::make_unique<RLAdapterPlacer>(python_module_dir);
}

} // namespace stdcell



