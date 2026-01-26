#include "stdcell/router.hpp"
#include <algorithm>

namespace stdcell {

Routed route_simple(const MatchResult& match, Placement placement, const TechRules& tr, const Netlist& nl) {
    Routed r; r.placement = placement;

    // Create rails on M0 (horizontal)
    const double y_vss = tr.rail_width_um * 0.5;
    const double y_vdd = r.placement.cell_h_um - tr.rail_width_um * 0.5;
    // Bottom rails for all VSS-family nets
    for (const auto& rn : nl.vss_nets) {
        WireSegment v; v.net = rn; v.layer = "M0"; v.x1 = 0; v.x2 = r.placement.cell_w_um; v.y1 = v.y2 = y_vss; r.wires.push_back(v);
    }
    // Top rails for all VDD-family nets
    for (const auto& rn : nl.vdd_nets) {
        WireSegment v; v.net = rn; v.layer = "M0"; v.x1 = 0; v.x2 = r.placement.cell_w_um; v.y1 = v.y2 = y_vdd; r.wires.push_back(v);
    }

    // Vertical M1 stubs from each non-rail pin
    for (const auto& p : r.placement.pins) {
        bool is_rail = (std::find(nl.vdd_nets.begin(), nl.vdd_nets.end(), p.name) != nl.vdd_nets.end() ||
                        std::find(nl.vss_nets.begin(), nl.vss_nets.end(), p.name) != nl.vss_nets.end());
        if (is_rail) continue;
        WireSegment s; s.net = p.name; s.layer = "M1"; s.x1 = s.x2 = p.x; s.y1 = 0; s.y2 = r.placement.cell_h_um; r.wires.push_back(s);
    }

    return r;
}

} // namespace stdcell
