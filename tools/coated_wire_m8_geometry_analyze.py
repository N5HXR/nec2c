#!/usr/bin/env python3
"""Compare paired NEC/OpenEMS results for the M8 geometry battery."""

import argparse
import json
import math
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", nargs=5,
                        metavar=("NAME", "NEC_BARE", "NEC_COATED",
                                 "OPENEMS_BARE", "OPENEMS_COATED"),
                        required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path)
    return parser.parse_args()


def load_nec(path):
    return json.loads(Path(path).read_text())


def load_openems(path):
    path = Path(path)
    spec = json.loads((path/"spec.json").read_text())
    far = json.loads((path/"far-field.json").read_text())
    current = json.loads((path/"wire-current.json").read_text())
    return {"resonance": {"frequency_mhz": far["frequency_mhz"],
                            "resistance_ohm": far["resistance_ohm"]},
            "directivity_dbi": far["directivity_dbi"],
            "radiation_efficiency": far["radiation_efficiency"],
            "pattern": far["samples"], "current": current["samples"],
            "mesh": spec["mesh"], "mesh_pitch_mm": spec["mesh_pitch_mm"],
            "wall_seconds": spec.get("wall_seconds"),
            "geometry": spec["geometry"]}


def by_key(rows, fields):
    return {tuple(row[field] for field in fields): row for row in rows}


def pattern_error(left, right):
    lookup = by_key(left, ("theta_deg", "phi_deg"))
    values = []
    for row in right:
        key = (row["theta_deg"], row["phi_deg"])
        if key in lookup and row["relative_db"] >= -40:
            values.append(abs(lookup[key]["relative_db"]-row["relative_db"]))
    return max(values)


def current_error(left, right):
    lookup = {row["name"]: row for row in left}
    return max(abs(lookup[row["name"]]["relative_magnitude"]-
                   row["relative_magnitude"]) for row in right)


def direct(case):
    n = case["nec"]
    o = case["openems"]
    return {
        "resonance_error_percent": 100*(n["resonance"]["frequency_mhz"]/
                                         o["resonance"]["frequency_mhz"]-1),
        "resistance_error_percent": 100*(n["resonance"]["resistance_ohm"]/
                                          o["resonance"]["resistance_ohm"]-1),
        "resistance_error_ohm": (n["resonance"]["resistance_ohm"]-
                                  o["resonance"]["resistance_ohm"]),
        "directivity_error_db": n["directivity_dbi"]-o["directivity_dbi"],
        "efficiency_error_percentage_points": 100*(
            n["radiation_efficiency"]-o["radiation_efficiency"]),
        "pattern_maximum_db_above_minus_40": pattern_error(n["pattern"],
                                                              o["pattern"]),
        "current_maximum_normalized_delta": current_error(n["current"],
                                                             o["current"]),
    }


def effect(bare, coated):
    def shifts(source):
        b, c = bare[source], coated[source]
        return {
            "resonance_shift_percent": 100*(
                c["resonance"]["frequency_mhz"]/
                b["resonance"]["frequency_mhz"]-1),
            "resistance_change_ohm": (c["resonance"]["resistance_ohm"]-
                                       b["resonance"]["resistance_ohm"]),
            "directivity_change_db": c["directivity_dbi"]-b["directivity_dbi"],
            "efficiency_change_percentage_points": 100*(
                c["radiation_efficiency"]-b["radiation_efficiency"]),
            "pattern_change_db": pattern_error(c["pattern"], b["pattern"]),
            "current_change": current_error(c["current"], b["current"]),
        }
    nec, openems = shifts("nec"), shifts("openems")
    return {"nec": nec, "openems": openems,
            "resonance_shift_error_percentage_points":
                nec["resonance_shift_percent"]-
                openems["resonance_shift_percent"],
            "resistance_change_error_ohm":
                nec["resistance_change_ohm"]-
                openems["resistance_change_ohm"],
            "directivity_change_error_db":
                nec["directivity_change_db"]-
                openems["directivity_change_db"],
            "efficiency_change_error_percentage_points":
                nec["efficiency_change_percentage_points"]-
                openems["efficiency_change_percentage_points"],
            "pattern_change_error_db": abs(nec["pattern_change_db"]-
                                             openems["pattern_change_db"]),
            "current_change_error": abs(nec["current_change"]-
                                          openems["current_change"])}


def compact(case):
    """Retain numerical witnesses without duplicating full pattern grids."""
    result = {}
    for solver, source in case.items():
        result[solver] = {
            "resonance": source["resonance"],
            "directivity_dbi": source["directivity_dbi"],
            "radiation_efficiency": source["radiation_efficiency"],
            "current": source["current"],
        }
        for key in ("mesh", "mesh_pitch_mm", "wall_seconds", "geometry"):
            if key in source:
                result[solver][key] = source[key]
    return result


def plot(cases, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = {
        "dogleg-b6": "dogleg\nb6",
        "dogleg-b12": "dogleg\nb12",
        "staged90-b6": "staged\nfeed1.5",
        "staged90-b6-gap6-control": "staged\nfeed6",
        "groundplane45-b6": "ground\n45deg b6",
    }
    names = [labels.get(case["name"], case["name"]) for case in cases]
    x = list(range(len(names)))
    resonance = [case["coating_effect"][
        "resonance_shift_error_percentage_points"] for case in cases]
    resistance = [case["coating_effect"][
        "resistance_change_error_ohm"] for case in cases]
    current = [case["coating_effect"]["current_change_error"] for case in cases]
    figure, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].bar(x, resonance)
    axes[1].bar(x, resistance)
    axes[2].bar(x, current)
    axes[0].set_ylabel("Coating-shift error (percentage points)")
    axes[1].set_ylabel("Coating resistance-change error (ohm)")
    axes[2].set_ylabel("Coating current-change error")
    for axis in axes:
        axis.set_xticks(x, names, fontsize=8)
        axis.axhline(0, color="black", linewidth=.7)
        axis.grid(True, axis="y", alpha=.25)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160)
    plt.close(figure)


def main():
    a = arguments()
    cases = []
    for name, nb, nc, ob, oc in a.case:
        bare = {"nec": load_nec(nb), "openems": load_openems(ob)}
        coated = {"nec": load_nec(nc), "openems": load_openems(oc)}
        nec_coated = coated["nec"]
        q = (math.sqrt(nec_coated["epsilon_r"])*2*math.pi*
             coated["openems"]["resonance"]["frequency_mhz"]*1e6/
             299792458.0*nec_coated["coating_radius_mm"]/1000)
        cases.append({"name": name, "bare": compact(bare),
                      "coated": compact(coated),
                      "sqrt_er_k0_b": q,
                      "direct_errors": {"bare": direct(bare),
                                        "coated": direct(coated)},
                      "coating_effect": effect(bare, coated)})
    record = {"schema": "nec2c-coated-wire-m8-geometry-analysis",
              "method": ("Paired bare/coated comparisons separate ordinary "
                         "geometry and feed error from coating-model error."),
              "cases": cases}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(record, indent=2)+"\n")
    if a.figure:
        plot(cases, a.figure)
    for case in cases:
        print(case["name"], json.dumps(case["coating_effect"], indent=2))


if __name__ == "__main__":
    main()
