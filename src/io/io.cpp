#include "stdcell/io.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cctype>
#include <filesystem>

namespace stdcell {

static inline std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

static inline std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out; std::string cur; std::istringstream iss(s);
    while (std::getline(iss, cur, sep)) out.push_back(cur);
    return out;
}

static inline std::string tolower_str(std::string s) { for (auto& c : s) c = (char)std::tolower((unsigned char)c); return s; }
static inline std::string toupper_str(std::string s) { for (auto& c : s) c = (char)std::toupper((unsigned char)c); return s; }

TechRules parse_tech_rules(const std::string& path) {
    TechRules tr;
    std::ifstream f(path);
    if (!f) return tr; // defaults
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string k = trim(line.substr(0, pos));
        std::string v = trim(line.substr(pos + 1));
        auto to_bool = [](const std::string& s){ return s == "1" || s == "true" || s == "True" || s == "TRUE"; };
        if (k == "height_tracks" || k == "tracks_per_row") tr.tracks_per_row = std::stoi(v);
        else if (k == "track_pitch_um") tr.track_pitch_um = std::stod(v);
        else if (k == "rail_width_um") tr.rail_width_um = std::stod(v);
        else if (k == "min_w_um.pmos") tr.min_w_pmos_um = std::stod(v);
        else if (k == "min_w_um.nmos") tr.min_w_nmos_um = std::stod(v);
        else if (k == "max_diff_w_um") tr.max_diff_w_um = std::stod(v);
        else if (k == "fold_enable") tr.fold_enable = to_bool(v);
        else if (k == "pin_pitch_tracks") tr.pin_pitch_tracks = std::stoi(v);
        else if (k == "preferred_m0_dir") tr.preferred_m0_dir = (v == "H" ? Dir::H : Dir::V);
        else if (k == "preferred_m1_dir") tr.preferred_m1_dir = (v == "H" ? Dir::H : Dir::V);
        else if (k == "default_w_per_fin_um" || k == "w_per_fin_um") tr.default_w_per_fin_um = std::stod(v);
    }
    return tr;
}

Netlist parse_cell_spec(const std::string& path) {
    Netlist nl;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("cell=", 0) == 0) {
            nl.cell_name = trim(line.substr(5));
        } else if (line.rfind("rails=", 0) == 0) {
            auto vals = split(trim(line.substr(6)), ',');
            for (auto& s : vals) nl.rails.push_back(trim(s));
        } else if (line.rfind("pins=", 0) == 0) {
            auto vals = split(trim(line.substr(5)), ',');
            for (auto& s : vals) nl.pins.push_back(trim(s));
        } else if (line.rfind("trans=", 0) == 0) {
            Transistor t;
            std::string rest = trim(line.substr(6));
            auto toks = split(rest, ' ');
            for (auto& tok : toks) {
                if (tok.empty()) continue;
                auto eq = tok.find('=');
                if (eq == std::string::npos) { t.id = tok; continue; }
                std::string k = tok.substr(0, eq);
                std::string v = tok.substr(eq + 1);
                if (k == "type") t.type = (v == "PMOS" ? TransType::PMOS : TransType::NMOS);
                else if (k == "W" || k == "w") t.w_um = std::stod(v);
                else if (k == "L" || k == "l") t.l_um = std::stod(v);
                else if (k == "G") t.g = v;
                else if (k == "S") t.s = v;
                else if (k == "D") t.d = v;
                else if (k == "B") t.b = v;
            }
            if (t.id.empty()) t.id = "M" + std::to_string((int)nl.devices.size() + 1);
            nl.devices.push_back(t);
        }
    }
    return nl;
}

static double parse_metric_value_to_um(const std::string& v) {
    std::string s = trim(v);
    if (s.empty()) return 0.0;
    char last = (char)std::tolower((unsigned char)s.back());
    if (last == 'u') { return std::atof(s.c_str()); }
    return std::atof(s.c_str()) * 1e6;
}

Netlist parse_spice_subckt(const std::string& lib_path, const std::string& subckt_name,
                           const std::vector<std::string>& rails_hint) {
    Netlist nl; nl.cell_name = subckt_name; nl.rails = rails_hint;
    std::ifstream f(lib_path);
    if (!f.good()) return nl;
    std::string line; bool in = false; std::vector<std::string> pins;
    while (std::getline(f, line)) {
        std::string s = trim(line);
        if (s.empty() || s[0] == '*') continue;
        std::string ls = tolower_str(s);
        if (!in) {
            if (ls.rfind(".subckt", 0) == 0) {
                std::istringstream iss(s);
                std::string dot, name; iss >> dot >> name;
                if (name == subckt_name) {
                    in = true; pins.clear();
                    std::string pin;
                    while (iss >> pin) {
                        if (pin.find('=') != std::string::npos) break;
                        pins.push_back(pin);
                    }
                    nl.pins = pins;
                }
            }
        } else {
            if (ls.rfind(".ends", 0) == 0) break;
            if (s.empty() || s[0] == '+') continue;
            if (!s.empty() && (s[0] == 'M' || s[0] == 'm')) {
                std::istringstream iss(s);
                std::string id, d, g, sn, b, model;
                iss >> id >> d >> g >> sn >> b >> model;
                Transistor t; t.id = id; t.g = g; t.s = sn; t.d = d; t.b = b; t.model_name = model;
                std::string ml = tolower_str(model);
                if (ml.find("pmos") != std::string::npos || ml.find("pfet") != std::string::npos || (!ml.empty() && ml[0]=='p')) t.type = TransType::PMOS; else t.type = TransType::NMOS;
                std::string kv;
                while (iss >> kv) {
                    auto eq = kv.find('='); if (eq == std::string::npos) continue;
                    std::string k = tolower_str(kv.substr(0, eq));
                    std::string v = kv.substr(eq + 1);
                    if (k == "w") t.w_um = parse_metric_value_to_um(v);
                    else if (k == "l") t.l_um = parse_metric_value_to_um(v);
                    else if (k == "nfin" || k == "nfins") t.nfin = std::max(0, (int)std::round(std::atof(v.c_str())));
                    else if (k == "nf" || k == "f" || k == "fingers") t.fingers = std::max(1, (int)std::round(std::atof(v.c_str())));
                }
                if (t.fingers <= 0) t.fingers = 1;
                nl.devices.push_back(t);
            }
        }
    }
    return nl;
}

RunConfig parse_run_config(const std::string& path) {
    RunConfig rc;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string k = trim(line.substr(0, pos));
        std::string v = trim(line.substr(pos + 1));
        std::string kl = tolower_str(k);
        auto set_if = [&](const char* key, std::string& dst){ if (kl == key) dst = v; };
        set_if("library_spice", rc.library_spice);
        set_if("tech_cfg", rc.tech_cfg);
        set_if("cells_list", rc.cells_list);
        set_if("single_cell", rc.single_cell);
        set_if("out_dir", rc.out_dir);
        set_if("match_in_dir", rc.match_in_dir);
        set_if("match_out_dir", rc.match_out_dir);
        set_if("place_in_dir", rc.place_in_dir);
        set_if("place_out_dir", rc.place_out_dir);
        set_if("route_out_dir", rc.route_out_dir);
        set_if("rails", rc.rails_csv);
        set_if("vdd_net", rc.vdd_net);
        set_if("vss_net", rc.vss_net);
        set_if("match_file", rc.match_file);
        set_if("place_file", rc.place_file);
        if (kl == "steps") {
            rc.steps.clear(); auto vlist = split(v, ',');
            for (auto& s2 : vlist) rc.steps.push_back(tolower_str(trim(s2)));
        } else if (kl == "use_rl") {
            std::string vl = tolower_str(v); rc.use_rl = (vl == "1" || vl == "true" || vl == "yes");
        } else if (kl == "match" || kl == "do_match") {
            std::string vl = tolower_str(v); rc.do_match = (vl == "1" || vl == "true" || vl == "yes");
        } else if (kl == "place" || kl == "do_place") {
            std::string vl = tolower_str(v); rc.do_place = (vl == "1" || vl == "true" || vl == "yes");
        } else if (kl == "route" || kl == "do_route") {
            std::string vl = tolower_str(v); rc.do_route = (vl == "1" || vl == "true" || vl == "yes");
        } else if (kl != "library_spice" && kl != "tech_cfg" && kl != "cells_list" && kl != "single_cell" && kl != "out_dir" && kl != "match_in_dir" && kl != "match_out_dir" && kl != "place_in_dir" && kl != "place_out_dir" && kl != "route_out_dir" && kl != "rails" && kl != "vdd_net" && kl != "vss_net" && kl != "match_file" && kl != "place_file") {
            rc.overrides[k] = v;
        }
    }
    for (auto& s2 : rc.steps) { if (s2=="match") rc.do_match = true; else if (s2=="place") rc.do_place = true; else if (s2=="route") rc.do_route = true; }
    return rc;
}

void apply_overrides(TechRules& tr, const RunConfig& rc) {
    for (const auto& kv : rc.overrides) {
        std::string k = tolower_str(kv.first); std::string v = kv.second;
        if (k == "tracks_per_row" || k == "height_tracks") tr.tracks_per_row = std::stoi(v);
        else if (k == "track_pitch_um") tr.track_pitch_um = std::stod(v);
        else if (k == "rail_width_um") tr.rail_width_um = std::stod(v);
        else if (k == "min_w_um.pmos") tr.min_w_pmos_um = std::stod(v);
        else if (k == "min_w_um.nmos") tr.min_w_nmos_um = std::stod(v);
        else if (k == "max_diff_w_um") tr.max_diff_w_um = std::stod(v);
        else if (k == "fold_enable") tr.fold_enable = (tolower_str(v) == "true" || v == "1");
        else if (k == "pin_pitch_tracks") tr.pin_pitch_tracks = std::stoi(v);
        else if (k == "preferred_m0_dir") tr.preferred_m0_dir = (v == "H" ? Dir::H : Dir::V);
        else if (k == "preferred_m1_dir") tr.preferred_m1_dir = (v == "H" ? Dir::H : Dir::V);
        else if (k == "default_w_per_fin_um" || k == "w_per_fin_um") tr.default_w_per_fin_um = std::stod(v);
    }
}

static void ensure_parent_dir(const std::string& file) {
    std::error_code ec; std::filesystem::create_directories(std::filesystem::path(file).parent_path(), ec);
}

void write_layout_json(const std::string& path, const Routed& routed) {
    ensure_parent_dir(path);
    std::ofstream o(path);
    const auto& pl = routed.placement;
    o << "{\n";
    o << "  \"cell\": \"" << pl.devices.size() << "_devices\",\n";
    o << std::fixed << std::setprecision(6);
    o << "  \"width_um\": " << pl.cell_w_um << ",\n";
    o << "  \"height_um\": " << pl.cell_h_um << ",\n";
    o << "  \"devices\": [\n";
    for (size_t i = 0; i < pl.devices.size(); ++i) {
        const auto& r = pl.devices[i];
        o << "    {\"name\":\"" << r.name << "\",\"x\":" << r.x << ",\"y\":" << r.y
          << ",\"w\":" << r.w << ",\"h\":" << r.h << "}";
        if (i + 1 != pl.devices.size()) o << ",";
        o << "\n";
    }
    o << "  ],\n";
    o << "  \"pins\": [\n";
    for (size_t i = 0; i < pl.pins.size(); ++i) {
        const auto& r = pl.pins[i];
        o << "    {\"name\":\"" << r.name << "\",\"x\":" << r.x << ",\"y\":" << r.y
          << ",\"w\":" << r.w << ",\"h\":" << r.h << "}";
        if (i + 1 != pl.pins.size()) o << ",";
        o << "\n";
    }
    o << "  ],\n";
    o << "  \"wires\": [\n";
    for (size_t i = 0; i < routed.wires.size(); ++i) {
        const auto& w = routed.wires[i];
        o << "    {\"net\":\"" << w.net << "\",\"layer\":\"" << w.layer << "\",\"x1\":"
          << w.x1 << ",\"y1\":" << w.y1 << ",\"x2\":" << w.x2 << ",\"y2\":" << w.y2 << "}";
        if (i + 1 != routed.wires.size()) o << ",";
        o << "\n";
    }
    o << "  ]\n";
    o << "}\n";
}

void write_layout_txt(const std::string& path, const Routed& routed) {
    ensure_parent_dir(path);
    std::ofstream o(path);
    const auto& pl = routed.placement;
    o << "Cell size: " << pl.cell_w_um << " x " << pl.cell_h_um << " um\n";
    o << "Devices: " << pl.devices.size() << ", Pins: " << pl.pins.size() << ", Wires: " << routed.wires.size() << "\n";
}

void write_match_json(const std::string& path, const MatchResult& mr) {
    ensure_parent_dir(path);
    std::ofstream o(path);
    o << "{\n  \"p_row\": [\n";
    auto dump_row = [&](const std::vector<FoldedDevice>& row){
        for (size_t i = 0; i < row.size(); ++i) {
            const auto& fd = row[i];
            o << "    {\"base_id\":\"" << fd.base_id << "\",\"fold_index\":" << fd.fold_index
              << ",\"type\":\"" << (fd.type == TransType::PMOS ? "PMOS" : "NMOS") << "\",\"w_um\":" << fd.w_um
              << ",\"g\":\"" << fd.g << "\",\"s\":\"" << fd.s << "\",\"d\":\"" << fd.d << "\",\"b\":\"" << fd.b << "\"}";
            if (i + 1 != row.size()) o << ",";
            o << "\n";
        }
    };
    dump_row(mr.p_row);
    o << "  ],\n  \"n_row\": [\n";
    dump_row(mr.n_row);
    o << "  ]\n}\n";
}

static bool parse_match_obj(const std::string& line, FoldedDevice& out) {
    auto get_str = [&](const char* key)->std::optional<std::string>{
        std::string k = std::string("\"") + key + "\":";
        size_t p = line.find(k); if (p == std::string::npos) return std::nullopt; p += k.size();
        if (p < line.size() && line[p] == '"') { ++p; size_t q = line.find('"', p); if (q != std::string::npos) return line.substr(p, q - p); }
        return std::nullopt;
    };
    auto get_num = [&](const char* key)->std::optional<double>{
        std::string k = std::string("\"") + key + "\":";
        size_t p = line.find(k); if (p == std::string::npos) return std::nullopt; p += k.size();
        size_t q = p; while (q < line.size() && (std::isdigit((unsigned char)line[q]) || line[q]=='-' || line[q]=='.' || line[q]=='e' || line[q]=='E')) ++q;
        return std::atof(line.substr(p, q - p).c_str());
    };
    auto base = get_str("base_id"); auto type = get_str("type"); auto g = get_str("g"); auto s = get_str("s"); auto d = get_str("d"); auto b = get_str("b");
    auto fold = get_num("fold_index"); auto w = get_num("w_um");
    if (!base || !type || !g || !s || !d || !b || !fold || !w) return false;
    out.base_id = *base; out.fold_index = (int)(*fold); out.type = (*type == "PMOS" ? TransType::PMOS : TransType::NMOS);
    out.w_um = *w; out.g = *g; out.s = *s; out.d = *d; out.b = *b; return true;
}

bool read_match_json(const std::string& path, MatchResult& mr) {
    std::ifstream in(path); if (!in.good()) return false;
    std::string line; bool in_p = false, in_n = false;
    while (std::getline(in, line)) {
        std::string s = trim(line);
        if (s.find("\"p_row\"") != std::string::npos) { in_p = true; in_n = false; continue; }
        if (s.find("\"n_row\"") != std::string::npos) { in_n = true; in_p = false; continue; }
        if (s == "]" || s == "],") { in_p = in_n = false; continue; }
        if (s.size() && s[0] == '{') {
            FoldedDevice fd; if (parse_match_obj(s, fd)) { if (in_p) mr.p_row.push_back(fd); else if (in_n) mr.n_row.push_back(fd); }
        }
    }
    return (!mr.p_row.empty() || !mr.n_row.empty());
}

void write_placement_json(const std::string& path, const Placement& pl) {
    ensure_parent_dir(path);
    std::ofstream o(path);
    o << std::fixed << std::setprecision(6);
    o << "{\n  \"width_um\": " << pl.cell_w_um << ",\n  \"height_um\": " << pl.cell_h_um << ",\n  \"devices\": [\n";
    for (size_t i = 0; i < pl.devices.size(); ++i) {
        const auto& r = pl.devices[i];
        o << "    {\"name\":\"" << r.name << "\",\"x\":" << r.x << ",\"y\":" << r.y
          << ",\"w\":" << r.w << ",\"h\":" << r.h << "}";
        if (i + 1 != pl.devices.size()) o << ",";
        o << "\n";
    }
    o << "  ],\n  \"pins\": [\n";
    for (size_t i = 0; i < pl.pins.size(); ++i) {
        const auto& r = pl.pins[i];
        o << "    {\"name\":\"" << r.name << "\",\"x\":" << r.x << ",\"y\":" << r.y
          << ",\"w\":" << r.w << ",\"h\":" << r.h << "}";
        if (i + 1 != pl.pins.size()) o << ",";
        o << "\n";
    }
    o << "  ]\n}\n";
}

bool read_placement_json(const std::string& path, Placement& pl) {
    std::ifstream in(path); if (!in.good()) return false;
    std::string line; while (std::getline(in, line)) {
        auto s = trim(line);
        if (s.rfind("\"width_um\"", 0) == 0) { auto pos = s.find(':'); if (pos!=std::string::npos) pl.cell_w_um = std::atof(s.substr(pos+1).c_str()); }
        else if (s.rfind("\"height_um\"", 0) == 0) { auto pos = s.find(':'); if (pos!=std::string::npos) pl.cell_h_um = std::atof(s.substr(pos+1).c_str()); }
    }
    return true;
}

} // namespace stdcell
