#!/usr/bin/env python3
"""Assemble compact M6 junction evidence and comparison figures."""

import argparse
import json
import math
from pathlib import Path


def parser():
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--nec", type=Path, required=True)
    value.add_argument("--coated-bare", type=Path, required=True)
    value.add_argument("--unequal", type=Path, required=True)
    value.add_argument("--output-dir", type=Path, required=True)
    return value


def load_case(directory, label):
    impedance = json.loads((directory / "impedance.json").read_text())
    current = json.loads((directory / "axial-current.json").read_text())
    far = json.loads((directory / "far-field.json").read_text())
    spec = json.loads((directory / "spec.json").read_text())
    frequency = far["frequency_mhz"]
    resistance = None
    for left, right in zip(impedance, impedance[1:]):
        if left["frequency_mhz"] <= frequency <= right["frequency_mhz"]:
            fraction = ((frequency-left["frequency_mhz"])/
                        (right["frequency_mhz"]-left["frequency_mhz"]))
            resistance = (left["impedance_ohm"]["re"] + fraction *
                          (right["impedance_ohm"]["re"]-
                           left["impedance_ohm"]["re"]))
            break
    if resistance is None:
        raise ValueError("resonance lies outside impedance sweep")
    samples = current["samples"]
    mirror_error = 0.0
    by_position = {round(row["z_mm"], 9): row for row in samples}
    for position, row in by_position.items():
        opposite = by_position.get(round(-position, 9))
        if opposite:
            mirror_error = max(mirror_error,
                               abs(row["relative_magnitude"]-
                                   opposite["relative_magnitude"]))
    return {
        "label": label,
        "resonance": {
            "frequency_mhz": frequency,
            "resistance_ohm": resistance,
        },
        "current_mirror_maximum_relative_delta": mirror_error,
        "directivity_dbi": far["directivity_db"],
        "end_criteria": spec["execution"]["end_criteria"],
        "wall_seconds": spec.get("wall_seconds"),
        "mesh": spec["mesh"],
        "positive_arm": {
            "material": spec["material"],
            **spec["material_parameters"],
        },
        "negative_arm": spec["negative_arm"],
        "impedance": impedance,
        "current": samples,
        "pattern": far["samples"],
    }


def plot(cases, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for case in cases:
        rows = case["impedance"]
        axes[0].plot([row["frequency_mhz"] for row in rows],
                     [row["impedance_ohm"]["im"] for row in rows],
                     label=case["label"])
        axes[1].plot([row["z_mm"] for row in case["current"]],
                     [row["relative_magnitude"] for row in case["current"]],
                     marker="o", label=case["label"])
        phi0 = [row for row in case["pattern"] if row["phi_deg"] == 0]
        axes[2].plot([row["theta_deg"] for row in phi0],
                     [max(-40, row["relative_db"]) for row in phi0],
                     label=case["label"])
    axes[0].axhline(0, color="black", linewidth=0.7)
    axes[0].set(xlabel="Frequency (MHz)", ylabel="Reactance (ohm)",
                xlim=(240, 330), ylim=(-150, 150))
    axes[1].set(xlabel="Axial position (mm)",
                ylabel="Normalized current magnitude", ylim=(0, 1.05))
    axes[2].set(xlabel="Theta (degrees)", ylabel="Relative field (dB)",
                ylim=(-40, 1))
    for axis in axes:
        axis.grid(True, alpha=0.25)
    axes[0].legend()
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def main():
    arguments = parser().parse_args()
    nec = json.loads(arguments.nec.read_text())
    cases = [load_case(arguments.coated_bare, "coated / bare"),
             load_case(arguments.unequal, "Er 3.2 / Er 2.2")]
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    plot(cases, arguments.output_dir / "m6-openems.png")
    compact_cases = []
    candidate_prefixes = ("coated_bare", "unequal_coating")
    for case, prefix in zip(cases, candidate_prefixes):
        compact = {key: value for key, value in case.items()
                   if key not in {"impedance", "current", "pattern"}}
        candidate_frequency = nec["experimental_candidates"][
            f"{prefix}_resonance_mhz"]
        candidate_resistance = nec["experimental_candidates"][
            f"{prefix}_resistance_ohm"]
        compact["experimental_nec_candidate_error"] = {
            "resonance_percent": 100*(candidate_frequency/
                                       case["resonance"]["frequency_mhz"]-1),
            "resistance_ohm": (candidate_resistance-
                                case["resonance"]["resistance_ohm"]),
        }
        compact_cases.append(compact)
    nec["openems_asymmetric_references"] = compact_cases
    nec["decision"]["openems_role"] = (
        "future transition-basis target; reciprocity rejects the present "
        "source-domain candidate before accuracy comparison")
    (arguments.output_dir / "m6-acceptance.json").write_text(
        json.dumps(nec, indent=2, sort_keys=True) + "\n")
    print(arguments.output_dir / "m6-acceptance.json")


if __name__ == "__main__":
    main()
