#!/usr/bin/env python3
"""Summarize and plot completed coated-wire OpenEMS M1 cases."""

import argparse
import json
import math
from pathlib import Path


def parse_case(value):
    try:
        label, directory = value.split("=", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("case must be LABEL=DIRECTORY") from error
    return label, Path(directory)


def interpolate_zero(rows):
    if rows and max(math.hypot(row["impedance_ohm"]["re"],
                               row["impedance_ohm"]["im"])
                    for row in rows) < 1.0e-12:
        raise ValueError("degenerate zero terminal impedance")
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 == 0 or x0 * x1 < 0:
            fraction = 0 if x0 == 0 else -x0 / (x1 - x0)
            frequency = left["frequency_mhz"] + fraction * (
                right["frequency_mhz"] - left["frequency_mhz"])
            resistance = left["impedance_ohm"]["re"] + fraction * (
                right["impedance_ohm"]["re"] -
                left["impedance_ohm"]["re"])
            return frequency, resistance
    raise ValueError("frequency sweep has no reactance zero")


def load_case(label, directory):
    spec = json.loads((directory / "spec.json").read_text())
    impedance = json.loads((directory / "impedance.json").read_text())
    try:
        frequency, resistance = interpolate_zero(impedance)
        diagnostic = None
    except ValueError as error:
        frequency = resistance = None
        diagnostic = str(error)
    far_path = directory / "far-field.json"
    return {
        "label": label,
        "directory": str(directory),
        "specification": spec,
        "resonance_mhz": frequency,
        "resistance_at_resonance_ohm": resistance,
        "diagnostic": diagnostic,
        "impedance": impedance,
        "far_field": (json.loads(far_path.read_text())
                      if far_path.exists() else None),
    }


def maximum_pattern_delta(left, right, floor_db=-40.0):
    if not left or not right:
        return None
    left_samples = {(row["theta_deg"], row["phi_deg"]): row["relative_db"]
                    for row in left["samples"]}
    right_samples = {(row["theta_deg"], row["phi_deg"]): row["relative_db"]
                     for row in right["samples"]}
    common = [key for key in left_samples.keys() & right_samples.keys()
              if max(left_samples[key], right_samples[key]) >= floor_db]
    return (max(abs(left_samples[key] - right_samples[key]) for key in common)
            if common else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case,
                        required=True, metavar="LABEL=DIRECTORY")
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = [load_case(*item) for item in args.case]
    compact = []
    for case in cases:
        spec = case["specification"]
        controls = spec["mesh_controls_mm"]
        radial_cells = controls.get("coating_radial_cells")
        if (radial_cells is None and spec["material"] == "coated" and
                spec["coordinates"] == "cylindrical"):
            geometry = spec["geometry_mm"]
            thickness = (geometry["coating_outer_radius"] -
                         geometry["conductor_radius"])
            radial_cells = max(4, math.ceil(
                thickness / controls["local_transverse"]))
        compact.append({
            "label": case["label"],
            "coordinates": spec["coordinates"],
            "material": spec["material"],
            "feed": spec["feed"],
            "cylindrical_source": spec.get("cylindrical_source"),
            "coating_free_end": spec["geometry_mm"].get(
                "coating_free_end", "flush"),
            "feed_gap_mm": spec["geometry_mm"]["feed_gap"],
            "end_criteria": spec["execution"]["end_criteria"],
            "mesh_controls_mm": controls,
            "coating_radial_cells": radial_cells,
            "mesh_cells": spec["mesh"]["total_cells"],
            "resonance_mhz": case["resonance_mhz"],
            "resistance_at_resonance_ohm":
                case["resistance_at_resonance_ohm"],
            "valid_terminal_response": case["diagnostic"] is None,
            "diagnostic": case["diagnostic"],
            "far_field": (case["far_field"]
                          if case["diagnostic"] is None else None),
        })

    flush_by_mesh = {}
    for row in compact:
        controls = row["mesh_controls_mm"]
        key = (controls["antenna_axial"], controls["local_transverse"],
               row["coating_radial_cells"])
        if (row["valid_terminal_response"] and
                row["material"] == "coated" and
                row["coating_free_end"] == "flush"):
            flush_by_mesh[key] = row
    for row in compact:
        controls = row["mesh_controls_mm"]
        key = (controls["antenna_axial"], controls["local_transverse"],
               row["coating_radial_cells"])
        baseline = flush_by_mesh.get(key)
        if baseline is None or not row["valid_terminal_response"]:
            continue
        row["delta_from_same_mesh_flush"] = {
            "resonance_mhz": (row["resonance_mhz"] -
                              baseline["resonance_mhz"]),
            "resistance_at_resonance_ohm": (
                row["resistance_at_resonance_ohm"] -
                baseline["resistance_at_resonance_ohm"]),
            "directivity_db": (
                row["far_field"]["directivity_db"] -
                baseline["far_field"]["directivity_db"]
                if row["far_field"] and baseline["far_field"] else None),
            "maximum_sampled_pattern_db_above_minus_40": maximum_pattern_delta(
                row["far_field"], baseline["far_field"]),
        }

    output = args.output_prefix
    output.parent.mkdir(parents=True, exist_ok=True)
    output.with_suffix(".json").write_text(
        json.dumps({"schema": "nec2c-coated-wire-m1-summary",
                    "schema_version": 1, "cases": compact},
                   indent=2, sort_keys=True) + "\n")

    figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for case in cases:
        frequency = [row["frequency_mhz"] for row in case["impedance"]]
        resistance = [row["impedance_ohm"]["re"]
                      for row in case["impedance"]]
        reactance = [row["impedance_ohm"]["im"]
                     for row in case["impedance"]]
        axes[0].plot(frequency, resistance, label=case["label"])
        axes[1].plot(frequency, reactance, label=case["label"])
    axes[0].set_ylabel("Resistance (ohm)")
    axes[1].set_ylabel("Reactance (ohm)")
    axes[1].set_xlabel("Frequency (MHz)")
    axes[1].axhline(0, color="black", linewidth=0.6)
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(alpha=0.25)
    axes[1].grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output.with_name(output.name + "-impedance.png"), dpi=160)
    plt.close(figure)

    groups = {}
    for row in compact:
        if not row["valid_terminal_response"]:
            continue
        key = (row["coordinates"], row["material"], row["feed"],
               row["cylindrical_source"], row["coating_free_end"],
               row["feed_gap_mm"], row["end_criteria"])
        groups.setdefault(key, []).append(row)
    figure, axes = plt.subplots(1, 2, figsize=(9, 4))
    for key, rows in groups.items():
        if len(rows) < 2:
            continue
        rows.sort(key=lambda row: row["mesh_controls_mm"]["local_transverse"],
                  reverse=True)
        pitch = [row["mesh_controls_mm"]["local_transverse"] for row in rows]
        label = f"{key[4]}, gap {key[5]:g} mm"
        axes[0].plot(pitch, [row["resonance_mhz"] for row in rows],
                     marker="o", label=label)
        axes[1].plot(pitch,
                     [row["resistance_at_resonance_ohm"] for row in rows],
                     marker="o", label=label)
    axes[0].set_ylabel("Resonance (MHz)")
    axes[1].set_ylabel("Resistance at resonance (ohm)")
    for axis in axes:
        axis.set_xlabel("Local radial/transverse pitch (mm)")
        axis.invert_xaxis()
        axis.grid(alpha=0.25)
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(output.with_name(output.name + "-convergence.png"), dpi=160)
    plt.close(figure)

    cap_rows = [row for row in compact
                if row["material"] == "coated" and
                row["valid_terminal_response"]]
    cap_types = {row["coating_free_end"] for row in cap_rows}
    if len(cap_types) > 1:
        order = {name: index for index, name in enumerate(
            ("flush", "flat-overhang", "hemisphere", "recessed"))}
        cap_rows.sort(key=lambda row: order.get(row["coating_free_end"], 99))
        # Labels remain unique when selected endpoint cases are repeated on a
        # finer through-thickness mesh.
        labels = [row["label"] for row in cap_rows]
        figure, axes = plt.subplots(1, 2, figsize=(10, 4))
        deltas = [row["delta_from_same_mesh_flush"] for row in cap_rows]
        axes[0].bar(labels, [row["resonance_mhz"] for row in deltas])
        axes[1].bar(labels, [row["resistance_at_resonance_ohm"]
                            for row in deltas])
        axes[0].set_ylabel("Resonance shift from same-mesh flush (MHz)")
        axes[1].set_ylabel("Resistance shift from same-mesh flush (ohm)")
        axes[0].axhline(0, color="black", linewidth=0.6)
        axes[1].axhline(0, color="black", linewidth=0.6)
        for axis in axes:
            axis.tick_params(axis="x", rotation=25)
            axis.grid(axis="y", alpha=0.25)
        figure.tight_layout()
        figure.savefig(output.with_name(output.name + "-caps.png"), dpi=160)
        plt.close(figure)
    print(output.with_suffix(".json"))


if __name__ == "__main__":
    main()
