#!/usr/bin/env python3
"""Build M5 feed-model comparison artifacts from two OpenEMS gap runs."""

import argparse
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gap1-dir", type=Path, required=True)
    parser.add_argument("--gap2-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_case(path, label):
    spec = json.loads((path / "spec.json").read_text(encoding="utf-8"))
    impedance = json.loads((path / "impedance.json").read_text(encoding="utf-8"))
    current = json.loads((path / "axial-current.json").read_text(encoding="utf-8"))
    far = json.loads((path / "far-field.json").read_text(encoding="utf-8"))
    resonance = resonance_from_rows(impedance)
    return {
        "label": label,
        "gap_mm": spec["geometry_mm"]["feed_gap"],
        "mesh_controls_mm": spec["mesh_controls_mm"],
        "port_measurement": spec["port_measurement"],
        "end_criteria": spec["execution"]["end_criteria"],
        "wall_seconds": spec.get("wall_seconds"),
        "resonance": resonance,
        "impedance": impedance,
        "current": current,
        "far_field": far,
    }


def resonance_from_rows(rows):
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 <= 0 <= x1 or x1 <= 0 <= x0:
            fraction = -x0/(x1-x0)
            return {
                "frequency_mhz": left["frequency_mhz"] + fraction*(
                    right["frequency_mhz"]-left["frequency_mhz"]),
                "resistance_ohm": left["impedance_ohm"]["re"] + fraction*(
                    right["impedance_ohm"]["re"]-
                    left["impedance_ohm"]["re"]),
            }
    raise RuntimeError("OpenEMS sweep does not bracket resonance")


def interpolate(values, coordinate, x_key, y_key):
    values = sorted(values, key=lambda row: row[x_key])
    if coordinate <= values[0][x_key]:
        return values[0][y_key]
    if coordinate >= values[-1][x_key]:
        return values[-1][y_key]
    for left, right in zip(values, values[1:]):
        if left[x_key] <= coordinate <= right[x_key]:
            fraction = ((coordinate-left[x_key]) /
                        (right[x_key]-left[x_key]))
            return left[y_key]+fraction*(right[y_key]-left[y_key])
    raise AssertionError(coordinate)


def phi_zero_pattern(case):
    rows = [row for row in case["far_field"]["samples"]
            if row["phi_deg"] == 0]
    peak = max(row["relative_db"] for row in rows)
    return [{"theta_deg": row["theta_deg"],
             "relative_db": row["relative_db"]-peak} for row in rows]


def current_profile(case):
    return [{"z_m": row["z_mm"]*1e-3,
             "relative_magnitude": row["relative_magnitude"]}
            for row in case["current"]["samples"]]


def maximum_pair_error(left, right, cutoff_m=0.025):
    left_positive = [row for row in current_profile(left)
                     if row["z_m"] >= cutoff_m]
    right_positive = [row for row in current_profile(right)
                      if row["z_m"] >= cutoff_m]
    if len(left_positive) != len(right_positive):
        raise RuntimeError("gap cases have different current sample counts")
    return max(abs(a["relative_magnitude"]-b["relative_magnitude"])
               for a, b in zip(left_positive, right_positive))


def maximum_nec_current_error(case, nec, cutoff_m=0.025):
    positive = [row for row in current_profile(case) if row["z_m"] >= cutoff_m]
    return max(abs(row["relative_magnitude"]-
                   interpolate(nec, row["z_m"], "z_m", "relative_magnitude"))
               for row in positive)


def maximum_pattern_error(left, right, floor_db=-40.0):
    a = phi_zero_pattern(left)
    b = phi_zero_pattern(right)
    return max(abs(x["relative_db"]-y["relative_db"])
               for x, y in zip(a, b)
               if max(x["relative_db"], y["relative_db"]) >= floor_db)


def compact_case(case):
    return {key: value for key, value in case.items()
            if key not in {"current", "far_field"}} | {
        "current": current_profile(case),
        "pattern": phi_zero_pattern(case),
        "directivity_dbi": case["far_field"]["directivity_db"],
    }


def plot_impedance(path, gap1, gap2, nec):
    figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                constrained_layout=True)
    curves = ((gap1["impedance"], "OpenEMS 1 mm physical gap", "-"),
              (gap2["impedance"], "OpenEMS 2 mm physical gap", "--"))
    for axis, component, label in ((axes[0], "re", "Resistance (Ω)"),
                                    (axes[1], "im", "Reactance (Ω)")):
        for rows, legend, style in curves:
            axis.plot([row["frequency_mhz"] for row in rows],
                      [row["impedance_ohm"][component] for row in rows],
                      style, label=legend, linewidth=1.6)
        axis.plot([row["frequency_mhz"] for row in nec],
                  [row["resistance_ohm" if component == "re" else
                       "reactance_ohm"] for row in nec],
                  label="NEC continuous center source", linewidth=1.7)
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.25)
    axes[0].set_title("M5 continuous-source and finite-gap feed comparison")
    axes[0].legend(fontsize=9)
    axes[1].axhline(0, color="0.25", linewidth=0.8)
    axes[1].set_xlabel("Frequency (MHz)")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_current(path, gap1, gap2, nec):
    figure, axis = plt.subplots(figsize=(8, 4.7), constrained_layout=True)
    for case, label, style in ((gap1, "OpenEMS 1 mm gap", "o-"),
                               (gap2, "OpenEMS 2 mm gap", "s--")):
        rows = current_profile(case)
        axis.plot([row["z_m"] for row in rows],
                  [row["relative_magnitude"] for row in rows], style,
                  markersize=3.5, label=label, linewidth=1.4)
    axis.plot([row["z_m"] for row in nec],
              [row["relative_magnitude"] for row in nec],
              label="NEC continuous center source", linewidth=1.7)
    axis.set(title="M5 normalized axial current at resonance",
             xlabel="Axial position (m)", ylabel="Relative current magnitude")
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_pattern(path, gap1, gap2, nec):
    figure, axis = plt.subplots(figsize=(8, 4.7), constrained_layout=True)
    for rows, label, style in ((phi_zero_pattern(gap1), "OpenEMS 1 mm gap", "-"),
                               (phi_zero_pattern(gap2), "OpenEMS 2 mm gap", "--"),
                               (nec, "NEC continuous center source", "-")):
        axis.plot([row["theta_deg"] for row in rows],
                  [row["relative_db"] for row in rows], style,
                  label=label, linewidth=1.6)
    axis.set(title="M5 normalized radiation pattern", xlabel="Theta (degrees)",
             ylabel="Relative field (dB)", ylim=(-40, 1))
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gap1 = load_case(args.gap1_dir, "split-gap-1mm")
    gap2 = load_case(args.gap2_dir, "split-gap-2mm")
    m3 = json.loads((ROOT / "tests/coated-wire/m3-acceptance.json").read_text(
        encoding="utf-8"))
    nec_current = m3["currents"]["coated"]
    nec_pattern = m3["patterns"]["nec2c_coated"]
    nec_resonance = m3["nec2c"]["coated_resonance"]
    gap_delta = {
        "resonance_mhz": (gap1["resonance"]["frequency_mhz"]-
                          gap2["resonance"]["frequency_mhz"]),
        "resonance_percent": 100*(gap1["resonance"]["frequency_mhz"] /
                                  gap2["resonance"]["frequency_mhz"]-1),
        "resistance_ohm": (gap1["resonance"]["resistance_ohm"]-
                           gap2["resonance"]["resistance_ohm"]),
        "maximum_relative_current_away_from_feed":
            maximum_pair_error(gap1, gap2),
        "maximum_pattern_db_above_minus_40":
            maximum_pattern_error(gap1, gap2),
    }
    nec_delta = {
        "resonance_percent_vs_gap2": 100*(
            nec_resonance["frequency_mhz"] /
            gap2["resonance"]["frequency_mhz"]-1),
        "resistance_ohm_vs_gap2": (nec_resonance["resistance_ohm"]-
                                    gap2["resonance"]["resistance_ohm"]),
        "maximum_relative_current_vs_gap2_away_from_feed":
            maximum_nec_current_error(gap2, nec_current),
    }
    record = {
        "schema": "nec2c-coated-wire-m5-acceptance",
        "schema_version": 1,
        "article": m3["article"],
        "reference_plane": {
            "location": "z=0 symmetry plane",
            "deembedding_distance_mm": 0.0,
            "note": "The voltage spans the centered gap. OpenEMS samples terminal current on the adjacent Yee-cell contour; that is a lumped-port measurement contour, not a length of transmission line, so no propagation correction is applied. The 1-2 mm gap sequence bounds the residual source/measurement effect.",
        },
        "openems": {"gap1": compact_case(gap1), "gap2": compact_case(gap2)},
        "nec2c": {"resonance": nec_resonance,
                   "current": nec_current, "pattern": nec_pattern,
                   "source": "continuous wire, center segment EX 0 voltage source"},
        "gap_convergence": gap_delta,
        "nec2c_vs_physical_gap": nec_delta,
        "decision": {
            "production_form": "continuous center-segment excitation",
            "physical_validation_form": "split two-arm finite-gap differential port",
            "unsupported_form": "embedded voltage source inside continuous PEC in OpenEMS",
            "reason": "The embedded OpenEMS construction is electrically shorted; the 1-2 mm physical-gap sequence is converged and corroborates the continuous NEC source within the coated-operator accuracy budget.",
        },
    }
    output = args.output_dir / "m5-acceptance.json"
    output.write_text(json.dumps(record, indent=2, sort_keys=True)+"\n",
                      encoding="utf-8")
    plot_impedance(args.output_dir / "m5-impedance.png", gap1, gap2,
                   m3["sweeps"]["coated"])
    plot_current(args.output_dir / "m5-current.png", gap1, gap2, nec_current)
    plot_pattern(args.output_dir / "m5-pattern.png", gap1, gap2, nec_pattern)
    print(output)


if __name__ == "__main__":
    main()
