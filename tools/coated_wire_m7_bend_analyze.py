#!/usr/bin/env python3
"""Assemble the M7 bend-envelope evidence and comparison figure."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT/"tests/coated-wire/m7-bend-study.json"
MATLAB = ROOT/"tests/coated-wire/m7-bend-matlab-result.json"
OPENEMS = ROOT/"tests/coated-wire/m7-bend-openems.json"
OUTPUT = ROOT/"tests/coated-wire/m7-bend-analysis.json"
FIGURE = ROOT/"tests/coated-wire/m7-bend-envelope.png"


def main():
    study = json.loads(STUDY.read_text())
    cases = {row["name"]: row for row in study["cases"]}
    matlab = json.loads(MATLAB.read_text())
    matlab_rows = {row["name"]: row for row in matlab["cases"]
                   if row["status"] == "ok"}
    comparisons = []
    for name, row in matlab_rows.items():
        nec = cases[name]["nec"]["bare"]
        reference = row["resonance"]
        comparisons.append({
            "name": name,
            "nec_frequency_mhz": nec["frequency_mhz"],
            "matlab_frequency_mhz": reference["frequency_mhz"],
            "frequency_error_percent": 100*(
                nec["frequency_mhz"]/reference["frequency_mhz"]-1),
            "nec_resistance_ohm": nec["resistance_ohm"],
            "matlab_resistance_ohm": reference["resistance_ohm"],
            "resistance_error_ohm": (nec["resistance_ohm"]-
                                     reference["resistance_ohm"]),
        })
    eligible = []
    for row in study["cases"]:
        if (row["maximum_turn_deg"] <= 45 and
                row["minimum_chord_ratio"] >= 4):
            bare = row["nec"]["bare"]
            coated = row["nec"]["coated"]
            eligible.append({
                "name": row["name"],
                "maximum_turn_deg": row["maximum_turn_deg"],
                "minimum_chord_ratio": row["minimum_chord_ratio"],
                "nec_coating_shift_percent": 100*(
                    coated["frequency_mhz"]/bare["frequency_mhz"]-1),
            })
    openems = (json.loads(OPENEMS.read_text()) if OPENEMS.exists() else
               {"runs": []})
    openems_runs = {row["key"]: row for row in openems.get("runs", [])}
    rounded = rounded_validation(cases["arc-n1-r4b"], openems_runs)
    record = {
        "schema": "nec2c-coated-wire-m7-bend-analysis",
        "engineering_bands": {
            "high_confidence": {
                "resonance_percent": 2.0,
                "resistance_percent": 10.0,
                "coating_shift_percentage_points": 1.0,
            },
            "basically_right": {
                "resonance_percent": 5.0,
                "resistance_percent": 25.0,
                "coating_shift_percentage_points": 2.0,
            },
            "note": ("These are reporting bands, not promotion gates; "
                     "reference convergence is required first."),
        },
        "nec_matlab": {
            "comparisons": comparisons,
            "maximum_absolute_frequency_error_percent": max(
                abs(row["frequency_error_percent"]) for row in comparisons),
            "maximum_absolute_resistance_error_ohm": max(
                abs(row["resistance_error_ohm"]) for row in comparisons),
        },
        "x45_y4_eligible_cases": eligible,
        "openems": openems,
        "rounded_n1_r4_validation": rounded,
    }
    OUTPUT.write_text(json.dumps(record, indent=2)+"\n")
    plot(study["cases"], matlab_rows, openems.get("runs", []))
    print(OUTPUT)
    print(FIGURE)


def rounded_validation(nec_case, runs):
    pairs = []
    for label in ("2mm", "1mm"):
        bare = runs.get(f"bare-{label}-n1-r4")
        coated = runs.get(f"coated-{label}-n1-r4")
        if bare is None or coated is None:
            continue
        bare_resonance = bare["resonance"]
        coated_resonance = coated["resonance"]
        pairs.append({
            "mesh_mm": bare["mesh_mm"],
            "bare": bare_resonance,
            "coated": coated_resonance,
            "coating_frequency_shift_percent": 100*(
                coated_resonance["frequency_mhz"]/
                bare_resonance["frequency_mhz"]-1),
            "coating_resistance_change_percent": 100*(
                coated_resonance["resistance_ohm"]/
                bare_resonance["resistance_ohm"]-1),
        })
    nec_bare = nec_case["nec"]["bare"]
    nec_coated = nec_case["nec"]["coated"]
    nec_shift = 100*(nec_coated["frequency_mhz"]/
                     nec_bare["frequency_mhz"]-1)
    fine_bare = runs.get("bare-0.5mm-n1-r4")
    return {
        "geometry": "arc-n1-r4b",
        "nec_coating_frequency_shift_percent": nec_shift,
        "matched_mesh_pairs": pairs,
        "finest_matched_shift_error_percentage_points": (
            abs(pairs[-1]["coating_frequency_shift_percent"]-nec_shift)
            if pairs else None),
        "finest_bare_openems": (fine_bare["resonance"]
                                 if fine_bare else None),
        "note": ("Matched bare/coated grids isolate the coating effect; the "
                 "0.5 mm bare run separately establishes geometry convergence."),
    }


def plot(cases, matlab, openems):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    colors = plt.cm.viridis([0.05, 0.25, 0.45, 0.65, 0.85])
    by_chords = {}
    for row in cases:
        if row["chords_per_corner"]:
            by_chords.setdefault(row["chords_per_corner"], []).append(row)
    for color, (chords, rows) in zip(colors, sorted(by_chords.items())):
        rows.sort(key=lambda row: row["radius_ratio"])
        x = [row["radius_ratio"] for row in rows]
        label = f"NEC, {chords} chord{'s' if chords != 1 else ''}"
        axes[0].plot(x, [row["nec"]["bare"]["frequency_mhz"]
                         for row in rows], color=color, label=label)
        axes[1].plot(x, [row["nec"]["bare"]["resistance_ohm"]
                         for row in rows], color=color)
        axes[2].plot(x, [100*(row["nec"]["coated"]["frequency_mhz"]/
                              row["nec"]["bare"]["frequency_mhz"]-1)
                         for row in rows], color=color)
    for name, row in matlab.items():
        source = next(item for item in cases if item["name"] == name)
        x = source["radius_ratio"]
        axes[0].scatter(x, row["resonance"]["frequency_mhz"], marker="x",
                        color="black", zorder=5)
        axes[1].scatter(x, row["resonance"]["resistance_ohm"], marker="x",
                        color="black", zorder=5)
    # OpenEMS is deliberately drawn only after a point has completed; missing
    # points remain visually absent instead of being mistaken for agreement.
    for row in openems:
        if row.get("status") != "ok":
            continue
        x = row["radius_ratio"]
        axes[0].scatter(x, row["resonance"]["frequency_mhz"], facecolors="none",
                        edgecolors="red", marker="o", zorder=6)
        axes[1].scatter(x, row["resonance"]["resistance_ohm"], facecolors="none",
                        edgecolors="red", marker="o", zorder=6)
    by_key = {row["key"]: row for row in openems
              if row.get("status") == "ok"}
    for label in ("2mm", "1mm"):
        bare = by_key.get(f"bare-{label}-n1-r4")
        coated = by_key.get(f"coated-{label}-n1-r4")
        if bare is not None and coated is not None:
            shift = 100*(coated["resonance"]["frequency_mhz"]/
                         bare["resonance"]["frequency_mhz"]-1)
            axes[2].scatter(4, shift, facecolors="none", edgecolors="red",
                            marker="o", zorder=6)
    axes[0].set(ylabel="Bare resonance (MHz)")
    axes[1].set(ylabel="Bare resistance at resonance (ohm)")
    axes[2].set(ylabel="NEC coating-induced shift (%)")
    for axis in axes:
        axis.set_xlabel("Bend radius / coating radius")
        axis.grid(True, alpha=0.25)
    axes[0].scatter([], [], marker="x", color="black", label="MATLAB surface MoM")
    axes[0].scatter([], [], marker="o", facecolors="none", edgecolors="red",
                    label="OpenEMS")
    axes[0].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE, dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    main()
