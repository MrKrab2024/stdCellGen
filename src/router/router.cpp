#include "stdcell/router.hpp"

#include <algorithm>

namespace stdcell {

Routed route_simple(const MatchResult& match, Placement placement, const TechRules& tr, const Netlist& nl) {
    Routed r; r.placement = placement;

    // Create VDD/VSS rails on M0 (horizontal only)
    const double y_vss = tr.rail_width_um * 0.5;
    const double y_vdd = r.placement.cell_h_um - tr.rail_width_um * 0.5;
    WireSegment vss; vss.net = "VSS"; vss.layer = "M0"; vss.x1 = 0; vss.x2 = r.placement.cell_w_um; vss.y1 = vss.y2 = y_vss; r.wires.push_back(vss);
    WireSegment vdd; vdd.net = "VDD"; vdd.layer = "M0"; vdd.x1 = 0; vdd.x2 = r.placement.cell_w_um; vdd.y1 = vdd.y2 = y_vdd; r.wires.push_back(vdd);

    // For demo, drop simple vertical M1 stubs from each pin to near cell middle
    for (const auto& p : r.placement.pins) {
        if (p.name == "VDD" || p.name == "VSS" || p.name == "GND") continue;
        WireSegment s; s.net = p.name; s.layer = "M1"; s.x1 = s.x2 = p.x; s.y1 = 0; s.y2 = r.placement.cell_h_um; r.wires.push_back(s);
    }

    return r;
}

} // namespace stdcell
