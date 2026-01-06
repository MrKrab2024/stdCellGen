#include "stdcell/io.hpp"
#include "stdcell/matcher.hpp"
#include "stdcell/matchdata.hpp"
#include "stdcell/structural_matcher.hpp"
#include "stdcell/placer.hpp"
#include "stdcell/router.hpp"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <numeric>
#include <iomanip>

using namespace stdcell;

static void usage() {
    std::cerr << "Usage: stdcellgen -config <config_file>\n";
}

static std::vector<std::string> split_csv_trim(const std::string& csv) {
    std::vector<std::string> out; std::stringstream ss(csv); std::string it;
    while (std::getline(ss, it, ',')) {
        size_t b = it.find_first_not_of(" \t"); if (b==std::string::npos) continue; size_t e = it.find_last_not_of(" \t");
        out.push_back(it.substr(b, e-b+1));
    }
    return out;
}

static void ensure_dirs(const RunConfig& rc) {
    auto default_dir = [&](const std::string& sub){ return (rc.out_dir.empty() ? std::string("out") : rc.out_dir) + "/" + sub; };
    if (rc.match_out_dir.empty()) (void)std::filesystem::create_directories(default_dir("match")); else (void)std::filesystem::create_directories(rc.match_out_dir);
    if (rc.place_out_dir.empty()) (void)std::filesystem::create_directories(default_dir("place")); else (void)std::filesystem::create_directories(rc.place_out_dir);
    if (rc.route_out_dir.empty()) (void)std::filesystem::create_directories(default_dir("route")); else (void)std::filesystem::create_directories(rc.route_out_dir);
}

struct MatchReportRow {
    std::string cell;
    bool used_fallback = false;
    std::string fail_reason;
    int group_count = 0;
    int pair_count = 0;
    int discrete_count = 0;
};

int main(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }

    std::string cfg_path;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "-config" || a == "--config") && i + 1 < argc) cfg_path = argv[++i];
    }
    if (cfg_path.empty()) { usage(); return 1; }

    RunConfig rc = parse_run_config(cfg_path);
    TechRules tr = parse_tech_rules(rc.tech_cfg);
    apply_overrides(tr, rc);

    // Rails
    std::vector<std::string> rails_hint;
    if (!rc.vdd_net.empty() && !rc.vss_net.empty()) { rails_hint = { rc.vdd_net, rc.vss_net }; }
    else { rails_hint = split_csv_trim(rc.rails_csv); }

    const bool do_match = rc.do_match;
    const bool do_place = rc.do_place;
    const bool do_route = rc.do_route;

    if (!do_match && !do_place && !do_route) {
        std::cerr << "Config error: no step enabled. Set match=true/place=true/route=true or use steps=...\n";
        return 2;
    }

    std::vector<std::string> cells;
    if (!rc.single_cell.empty()) { cells.push_back(rc.single_cell); }
    else if (!rc.cells_list.empty()) {
        std::ifstream f(rc.cells_list);
        std::string line; while (std::getline(f, line)) {
            std::string s = line; auto pos = s.find_first_of("#*;"); if (pos != std::string::npos) s = s.substr(0, pos);
            size_t b = s.find_first_not_of(" \t\r"); if (b==std::string::npos) continue; size_t e = s.find_last_not_of(" \t\r");
            cells.push_back(s.substr(b, e-b+1));
        }
    }

    auto first_enabled = [&](){ if (do_match) return std::string("match"); if (do_place) return std::string("place"); return std::string("route"); }();
    if (first_enabled == "match") {
        if (rc.library_spice.empty() || (cells.empty())) { std::cerr << "Config error: match requested but missing inputs.\n"; return 3; }
    } else if (first_enabled == "place") {
        if (rc.single_cell.size()) { if (rc.match_file.empty() && rc.match_in_dir.empty() && rc.match_out_dir.empty()) { std::cerr << "Config error: place-only needs match_file or match dir.\n"; return 4; } }
        else if (!cells.empty()) { if (rc.match_in_dir.empty() && rc.match_out_dir.empty()) { std::cerr << "Config error: place-only needs match_in_dir or match_out_dir.\n"; return 5; } }
        else { std::cerr << "Config error: place requested but no cells specified.\n"; return 6; }
    } else {
        if (rc.single_cell.size()) { if (rc.place_file.empty() && rc.place_in_dir.empty() && rc.place_out_dir.empty()) { std::cerr << "Config error: route-only needs place_file or place dir.\n"; return 7; } }
        else if (!cells.empty()) { if (rc.place_in_dir.empty() && rc.place_out_dir.empty()) { std::cerr << "Config error: route-only needs place_in_dir or place_out_dir.\n"; return 8; } }
        else { std::cerr << "Config error: route requested but no cells specified.\n"; return 9; }
    }

    ensure_dirs(rc);

    auto match_path_for = [&](const std::string& cell){ if (!rc.match_in_dir.empty()) return rc.match_in_dir + "/" + cell + ".match.json"; if (!rc.match_out_dir.empty()) return rc.match_out_dir + "/" + cell + ".match.json"; return std::string("out/match/") + cell + ".match.json"; };
    auto place_path_for = [&](const std::string& cell){ if (!rc.place_in_dir.empty()) return rc.place_in_dir + "/" + cell + ".place.json"; if (!rc.place_out_dir.empty()) return rc.place_out_dir + "/" + cell + ".place.json"; return std::string("out/place/") + cell + ".place.json"; };

    MatchConfig mc = build_match_config(tr, rails_hint);

    std::vector<MatchReportRow> report_rows;

    // analysis dir for CSV/debug
    std::string analysis_dir = (rc.out_dir.empty()? std::string("out") : rc.out_dir) + "/analysis";
    std::string debug_dir = analysis_dir + "/debug";
    std::error_code ec_mk;
    std::filesystem::create_directories(analysis_dir, ec_mk);
    std::filesystem::create_directories(debug_dir, ec_mk);

    auto json_escape = [](const std::string& s){ std::string o; o.reserve(s.size()+2); for(char c: s){ if(c=='\\'||c=='"') o.push_back('\\'), o.push_back(c); else if(c=='\n') { o += "\\n"; } else o.push_back(c);} return o; };

    auto run_cell = [&](const std::string& cell){
        Netlist nl;
        Placement pl;

        if (do_match) {
            nl = parse_spice_subckt(rc.library_spice, cell, rails_hint);
            if (nl.devices.empty()) {
                std::cerr << "Warning: subckt parse failed for " << cell << ". Skipping.\n"; return; }
            // structural matching with status
            StructuralMatchOutput mout = match_structural_pairs(nl, tr);
            // summarize
            MatchReportRow row; row.cell = cell; row.used_fallback = mout.used_fallback; row.fail_reason = mout.fail_reason; row.group_count = (int)mout.groups.size();
            int pairs_sum = 0; for (const auto& g : mout.groups) pairs_sum += (int)g.pairs.size(); row.pair_count = pairs_sum; row.discrete_count = (int)mout.discrete.size();
            report_rows.push_back(row);

            // write per-cell debug JSON only for fails
            if (mout.used_fallback) {
                std::string dbg = debug_dir + "/" + cell + ".debug.json";
                std::ofstream dj(dbg);
                dj << "{\n  \"cell\": \"" << json_escape(cell) << "\",\n  \"used_fallback\": " << (mout.used_fallback?"true":"false") << ",\n  \"fail_reason\": \"" << json_escape(mout.fail_reason) << "\",\n  \"groups\": [\n";
                for (size_t gi=0; gi<mout.groups.size(); ++gi) {
                    const auto& g = mout.groups[gi];
                    dj << "    {\"level\": " << g.level << ", \"pairs\": [\n";
                    for (size_t pi=0; pi<g.pairs.size(); ++pi) {
                        const auto& p = g.pairs[pi];
                        dj << "      {\"pmos\":{\"id\":\"" << json_escape(p.pmos.id) << "\",\"g\":\"" << json_escape(p.pmos.g) << "\",\"s\":\"" << json_escape(p.pmos.s) << "\",\"d\":\"" << json_escape(p.pmos.d) << "\"},"
                           << "\"nmos\":{\"id\":\"" << json_escape(p.nmos.id) << "\",\"g\":\"" << json_escape(p.nmos.g) << "\",\"s\":\"" << json_escape(p.nmos.s) << "\",\"d\":\"" << json_escape(p.nmos.d) << "\"}}";
                        if (pi+1 != g.pairs.size()) dj << ",";
                        dj << "\n";
                    }
                    dj << "    ]}";
                    if (gi+1 != mout.groups.size()) dj << ",";
                    dj << "\n";
                }
                dj << "  ],\n  \"discrete\": [\n";
                for (size_t di=0; di<mout.discrete.size(); ++di) {
                    const auto& m = mout.discrete[di];
                    dj << "    {\"type\":\"" << (m.type==TransType::PMOS?"PMOS":"NMOS") << "\",\"id\":\"" << json_escape(m.id) << "\",\"g\":\"" << json_escape(m.g) << "\",\"s\":\"" << json_escape(m.s) << "\",\"d\":\"" << json_escape(m.d) << "\"}";
                    if (di+1 != mout.discrete.size()) dj << ",";
                    dj << "\n";
                }
                dj << "  ]\n}\n";
            }

            // convert to legacy for downstream
            MatchResult mr = pairs_to_match_result(mout);
            std::string mpath = match_path_for(cell);
            write_match_json(mpath, mr);
        } else if (do_place || do_route) {
            // When not matching, we still need nets/rails at least for route
            nl.cell_name = cell; nl.rails = rails_hint;
        }

        // Place stage
        if (do_place) {
            MatchResult mr;
            if (!do_match) {
                std::string minpath = match_path_for(cell);
                if (!read_match_json(minpath, mr)) { std::cerr << "Error: cannot read match json: " << minpath << "\n"; return; }
            } else {
                std::string minpath = match_path_for(cell);
                read_match_json(minpath, mr);
            }
            std::unique_ptr<IPlacer> placer = rc.use_rl ? make_rl_adapter_placer("python/rl_placer") : make_heuristic_placer();
            pl = placer->place(mr, tr, nl);
            std::string ppath = place_path_for(cell);
            write_placement_json(ppath, pl);
        } else if (do_route) {
            std::string ppath = place_path_for(cell);
            if (!read_placement_json(ppath, pl)) { std::cerr << "Error: cannot read placement json: " << ppath << "\n"; return; }
        }

        // Route stage
        if (do_route) {
            if (nl.cell_name.empty()) { nl.cell_name = cell; nl.rails = rails_hint; }
            Routed routed = route_simple(MatchResult{}, pl, tr, nl);
            std::string outdir = (rc.route_out_dir.empty()? std::string("out/route") : rc.route_out_dir);
            std::string json_out = outdir + "/" + cell + ".layout.json";
            std::string txt_out = outdir + "/" + cell + ".layout.txt";
            write_layout_json(json_out, routed);
            write_layout_txt(txt_out, routed);
            std::cout << "Wrote " << json_out << " and " << txt_out << "\n";
        }
    };

    if (cells.empty()) { std::cerr << "No cells specified to run.\n"; return 10; }
    for (const auto& cell : cells) run_cell(cell);

    // Print match summary if we did match
    if (do_match && !report_rows.empty()) {
        int total = (int)report_rows.size();
        int fallback_cnt = 0; for (auto& r : report_rows) if (r.used_fallback) ++fallback_cnt;
        std::cout << "\nMatch Summary:" << std::endl;
        std::cout << "Total cells: " << total << ", Fallback: " << fallback_cnt
                  << " (" << std::fixed << std::setprecision(1) << (100.0 * fallback_cnt / std::max(1,total)) << "%)" << std::endl;
        std::cout << "Details:" << std::endl;
        for (const auto& r : report_rows) {
            std::cout << "- " << r.cell << ": " << (r.used_fallback ? "FAIL" : "OK")
                      << ", groups=" << r.group_count << ", pairs=" << r.pair_count << ", discrete=" << r.discrete_count;
            if (r.used_fallback && !r.fail_reason.empty()) std::cout << ", reason=" << r.fail_reason;
            std::cout << std::endl;
        }
        // Write CSV
        std::string csv_path = analysis_dir + "/asap7_match_summary.csv";
        std::ofstream csv(csv_path);
        csv << "cell,status,groups,pairs,discrete,reason\n";
        for (const auto& r : report_rows) {
            std::string status = r.used_fallback ? "FAIL" : "OK";
            // Quote reason and escape quotes
            std::string reason = r.fail_reason; for (char& c : reason) if (c=='"') c='\'';
            csv << r.cell << "," << status << "," << r.group_count << "," << r.pair_count << "," << r.discrete_count << ",\"" << reason << "\"\n";
        }
        std::cout << "CSV summary: " << csv_path << std::endl;
        std::cout << "Debug JSON (fails): " << debug_dir << std::endl;
    }

    return 0;
}
