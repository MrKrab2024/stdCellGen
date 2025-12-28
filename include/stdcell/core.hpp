#pragma once

#include <string>
#include <vector>
#include <map>
#include <optional>
#include <cstdint>

namespace stdcell {

enum class Dir { H, V };

inline const char* to_cstr(Dir d) { return d == Dir::H ? "H" : "V"; }

enum class TransType { PMOS, NMOS };

inline const char* to_cstr(TransType t) { return t == TransType::PMOS ? "PMOS" : "NMOS"; }

struct TechRules {
    int tracks_per_row = 9;
    double track_pitch_um = 0.19;
    double rail_width_um = 0.40;
    double min_w_pmos_um = 0.14;
    double min_w_nmos_um = 0.14;
    double max_diff_w_um = 0.50;
    bool fold_enable = true;
    int pin_pitch_tracks = 1;
    Dir preferred_m0_dir = Dir::H;
    Dir preferred_m1_dir = Dir::V;
    double default_w_per_fin_um = 0.0; // optional mapping from width to nfin
};

struct Transistor {
    std::string id;       // e.g. M1
    TransType type = TransType::PMOS;
    double w_um = 0.14;   // width
    double l_um = 0.05;   // length
    std::string g;        // gate net
    std::string s;        // source net
    std::string d;        // drain net
    std::string b;        // bulk net
    // SPICE extras
    std::string model_name; // model card/name
    int nfin = 0;           // number of fins (if provided or derived)
    int fingers = 1;        // number of gate fingers (nf)
};

struct Netlist {
    std::string cell_name;
    std::vector<std::string> rails;   // e.g. VDD,VSS
    std::vector<std::string> pins;    // e.g. A,Y
    std::vector<Transistor> devices;
};

struct FoldedDevice {
    std::string base_id;
    int fold_index = 0;
    TransType type = TransType::PMOS;
    double w_um = 0.0;
    std::string g, s, d, b;
};

struct MatchResult {
    std::vector<FoldedDevice> p_row;
    std::vector<FoldedDevice> n_row;
};

struct Rect {
    double x = 0, y = 0, w = 0, h = 0;
    std::string name;
};

struct WireSegment {
    std::string net;
    std::string layer;
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
};

struct Placement {
    double cell_w_um = 0.0;
    double cell_h_um = 0.0;
    std::vector<Rect> devices;
    std::vector<Rect> pins;
};

struct Routed {
    Placement placement;
    std::vector<WireSegment> wires;
};

struct RunConfig {
    std::string library_spice;
    std::string tech_cfg;
    std::string cells_list;
    std::string single_cell;
    std::string out_dir = "out";
    std::vector<std::string> steps;
    bool use_rl = false;
    std::string rails_csv = "VDD,VSS";
    std::string vdd_net;
    std::string vss_net;
    bool do_match = false;
    bool do_place = false;
    bool do_route = false;
    std::string match_file;
    std::string place_file;
    std::string match_in_dir;
    std::string match_out_dir;
    std::string place_in_dir;
    std::string place_out_dir;
    std::string route_out_dir;
    std::map<std::string, std::string> overrides;
};

} // namespace stdcell
