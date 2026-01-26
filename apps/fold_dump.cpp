#include "stdcell/io.hpp"
#include "stdcell/folding.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>

using namespace stdcell;

static bool parse_pairs_bundle_cells(const std::string& path,
    std::string& netlist_path,
    std::unordered_map<std::string, std::vector<std::pair<std::string,std::string>>>& out)
{
    std::ifstream in(path); if (!in.good()) return false; std::string s((std::istreambuf_iterator<char>(in)), {});
    // find netlist
    {
        std::string key = "\"netlist\""; size_t p = s.find(key); if (p != std::string::npos) {
            size_t q1 = s.find('"', p + key.size()); if (q1 != std::string::npos) {
                size_t q2 = s.find('"', q1 + 1); if (q2 != std::string::npos) netlist_path = s.substr(q1 + 1, q2 - q1 - 1);
            }
        }
    }
    size_t cells_pos = s.find("\"cells\""); if (cells_pos == std::string::npos) return false;
    size_t brace = s.find('{', cells_pos); if (brace == std::string::npos) return false;
    size_t pos = brace + 1;
    while (true) {
        size_t q1 = s.find('"', pos); if (q1 == std::string::npos) break;
        size_t q2 = s.find('"', q1 + 1); if (q2 == std::string::npos) break;
        std::string cell = s.substr(q1 + 1, q2 - q1 - 1);
        size_t colon = s.find(':', q2 + 1); if (colon == std::string::npos) break;
        size_t cb = s.find('{', colon + 1); if (cb == std::string::npos) break;
        int depth = 1; size_t i = cb + 1; for (; i < s.size() && depth > 0; ++i) { if (s[i] == '{') ++depth; else if (s[i] == '}') --depth; }
        if (depth != 0) break; size_t ce = i - 1;
        std::string block = s.substr(cb, ce - cb + 1);
        std::vector<std::pair<std::string, std::string>> pairs;
        size_t bpos = 0;
        while (true) {
            size_t lb = block.find('[', bpos); if (lb == std::string::npos) break;
            size_t rb = block.find(']', lb + 1); if (rb == std::string::npos) break;
            std::string arr = block.substr(lb + 1, rb - lb - 1);
            size_t a1 = arr.find('"'); if (a1 == std::string::npos) { bpos = rb + 1; continue; }
            size_t a2 = arr.find('"', a1 + 1); if (a2 == std::string::npos) { bpos = rb + 1; continue; }
            size_t comma = arr.find(',', a2 + 1); if (comma == std::string::npos) { bpos = rb + 1; continue; }
            size_t b1 = arr.find('"', comma + 1); if (b1 == std::string::npos) { bpos = rb + 1; continue; }
            size_t b2 = arr.find('"', b1 + 1); if (b2 == std::string::npos) { bpos = rb + 1; continue; }
            std::string id1 = arr.substr(a1 + 1, a2 - a1 - 1);
            std::string id2 = arr.substr(b1 + 1, b2 - b1 - 1);
            // assume (PMOS, NMOS)
            if (!id1.empty() && !id2.empty()) pairs.emplace_back(id1, id2);
            bpos = rb + 1;
        }
        out[cell] = std::move(pairs);
        pos = ce + 1;
        size_t nextComma = s.find(',', pos);
        if (nextComma != std::string::npos && nextComma < s.find('"', pos)) pos = nextComma + 1;
    }
    return true;
}

int main(int argc, char** argv) {
    std::string bundle = (argc > 1 ? argv[1] : std::string("out/oriented_asap7/pairs.bundle.json"));
    std::string netlist_path;
    std::unordered_map<std::string, std::vector<std::pair<std::string,std::string>>> cells_pairs;
    if (!parse_pairs_bundle_cells(bundle, netlist_path, cells_pairs)) {
        std::cerr << "Failed to read pairs bundle: " << bundle << "\n"; return 1;
    }
    if (netlist_path.empty()) netlist_path = "out/oriented_asap7.sp";

    // Prepare out dir
    std::string analysis_dir = "out/analysis";
    std::string folding_dir = analysis_dir + "/folding";
    std::error_code ec; std::filesystem::create_directories(folding_dir, ec);

    TechRules tr = parse_tech_rules("data/tech.demo.cfg"); // default demo
    int max_fpr = 3; int num_logical_rows = 2;
    std::vector<std::string> vdd = { "VDD" }, vss = { "VSS" };

    int cells_ok = 0;
    for (const auto& kv : cells_pairs) {
        const std::string& cell = kv.first; const auto& plist = kv.second;
        // parse netlist for this subckt
        Netlist nl = parse_spice_subckt(netlist_path, cell, vdd, vss);
        if (nl.devices.empty()) { std::cerr << "Skip cell (no netlist): " << cell << "\n"; continue; }
        std::unordered_map<std::string, const Transistor*> by_id; for (const auto& t : nl.devices) by_id[t.id] = &t;
        std::string outjson = folding_dir + "/" + cell + ".fold.json";
        std::ofstream fj(outjson);
        fj << "{\n  \"cell\": \"" << cell << "\",\n  \"params\": {\"max_fins_per_row\": " << max_fpr << ", \"num_logical_rows\": " << num_logical_rows << "},\n  \"pairs\": [\n";
        for (size_t pi=0; pi<plist.size(); ++pi) {
            const auto& pr = plist[pi];
            const Transistor* tp = (by_id.count(pr.first)? by_id[pr.first] : nullptr);
            const Transistor* tn = (by_id.count(pr.second)? by_id[pr.second] : nullptr);
            if (!tp || !tn) continue;
            auto vars = enumerate_pair_variants(*tn, *tp, max_fpr, num_logical_rows, tr);
            // filter minimal-finger splits
            int nf_n = effective_nfin(*tn, tr);
            int nf_p = effective_nfin(*tp, tr);
            int minN = (nf_n + max_fpr - 1) / max_fpr;
            int minP = (nf_p + max_fpr - 1) / max_fpr;
            std::vector<PairVariant> filtered; filtered.reserve(vars.size());
            for (const auto& v : vars) {
                int cn=0, cp=0; for (const auto& f : v.fingers) { if (f.type==TransType::NMOS) ++cn; else ++cp; }
                if (cn==minN && cp==minP) filtered.push_back(v);
            }
            fj << "    {\"pmos_id\":\"" << pr.first << "\", \"nmos_id\":\"" << pr.second << "\", \"variants\": " << ((filtered.empty()? vars : filtered).size()) << ", \"examples\": [\n";
            auto& use = filtered.empty()? vars : filtered;
                // enforce adjacency: n/p must be on the same logical row (floor(y/2) equal)
                std::vector<PairVariant> aligned; aligned.reserve(use.size());
                for (const auto& v : use) {
                    std::set<int> rn, rp;
                    for (const auto& f : v.fingers) {
                        if (f.type==TransType::NMOS) rn.insert(f.y/2); else rp.insert(f.y/2);
                    }
                    if (rn.size()==1 && rp.size()==1 && *rn.begin()==*rp.begin()) aligned.push_back(v);
                }
                auto& use2 = aligned.empty()? use : aligned;
            size_t show = std::min<size_t>(use.size(), 3);
            for (size_t vi=0; vi<show; ++vi) {
                fj << "      [\n";
                const auto& fingers = use[vi].fingers;
                for (size_t fi=0; fi<fingers.size(); ++fi) {
                    const auto& f = fingers[fi];
                    fj << "        {\"id\":\"" << f.base_id << "\", \"idx\":" << f.finger_index
                       << ", \"type\":\"" << (f.type==TransType::PMOS?"PMOS":"NMOS")
                       << "\", \"nfin\":" << f.nfin << ", \"x\":" << f.x << ", \"y\":" << f.y
                       << ", \"orient\":\"" << (f.orient==Orient::SD?"SD":"DS") << "\"}";
                    if (fi+1 != fingers.size()) fj << ",";
                    fj << "\n";
                }
                fj << "      ]"; if (vi+1 != show) fj << ","; fj << "\n";
            }
            fj << "    ]}"; if (pi+1 != plist.size()) fj << ","; fj << "\n";
        }
        fj << "  ]\n}\n";
        ++cells_ok;
    }
    std::cout << "Folding dump written for " << cells_ok << " cells into " << folding_dir << "\n";
    return 0;
}
