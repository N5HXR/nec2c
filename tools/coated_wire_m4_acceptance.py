#!/usr/bin/env python3
"""Generate M4 VSWR, loss, and power-balance acceptance artifacts."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--openems-dir", type=Path)
    parser.add_argument("--openems-end-energy-db", type=float)
    return parser.parse_args()


def run(binary, deck, directory):
    result = directory / f"{deck.stem}.json"
    report = directory / f"{deck.stem}.out"
    completed = subprocess.run([
        str(binary), "--backend", "builtin", "--jobs", "1", "--threads", "16",
        "--result-json", str(result), "-i", str(deck), "-o", str(report),
    ], cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(f"{deck} failed\n{completed.stdout}\n{completed.stderr}")
    return json.loads(result.read_text(encoding="utf-8"))


def rows(result):
    values = []
    for execution in result["jobs"][0]["executions"]:
        impedance = execution["feeds"][0]["impedance_ohm"]
        z = complex(impedance["re"], impedance["im"])
        gamma = abs((z-50)/(z+50))
        power = execution["power"]
        values.append({
            "frequency_mhz": execution["frequency_mhz"],
            "resistance_ohm": z.real,
            "reactance_ohm": z.imag,
            "vswr_50_ohm": (1+gamma)/(1-gamma),
            "accepted_w": power["accepted_w"],
            "radiated_w": power["air_radiated_w"],
            "dielectric_absorbed_w": power["dielectric_absorbed_w"],
            "radiation_efficiency": (power["air_radiated_w"] /
                                     power["accepted_w"]),
            "balance_defect_relative": power["balance_defect_relative"],
        })
    return values


def resonance(values):
    for left, right in zip(values, values[1:]):
        if left["reactance_ohm"] <= 0 <= right["reactance_ohm"]:
            fraction = -left["reactance_ohm"] / (
                right["reactance_ohm"]-left["reactance_ohm"])
            return {
                key: left[key]+fraction*(right[key]-left[key])
                for key in ("frequency_mhz", "resistance_ohm")
            }
    raise RuntimeError("sweep does not bracket resonance")


def interpolate(values, frequency, key):
    for left, right in zip(values, values[1:]):
        if left["frequency_mhz"] <= frequency <= right["frequency_mhz"]:
            fraction = ((frequency-left["frequency_mhz"]) /
                        (right["frequency_mhz"]-left["frequency_mhz"]))
            return left[key]+fraction*(right[key]-left[key])
    raise RuntimeError(f"frequency {frequency} is outside sweep")


def openems_rows(path):
    source = json.loads((path / "impedance.json").read_text(encoding="utf-8"))
    return [{"frequency_mhz": row["frequency_mhz"],
             "resistance_ohm": row["impedance_ohm"]["re"],
             "reactance_ohm": row["impedance_ohm"]["im"],
             "accepted_w": row["accepted_power_w"]} for row in source
            if 250.0 <= row["frequency_mhz"] <= 320.0]


def plot_vswr(path, low, lossy):
    figure, axis = plt.subplots(figsize=(8, 4.7), constrained_layout=True)
    for values, label in ((low, "tanδ = 0.002"), (lossy, "tanδ = 0.05")):
        axis.plot([row["frequency_mhz"] for row in values],
                  [row["vswr_50_ohm"] for row in values], label=label,
                  linewidth=1.8)
    axis.set(title="M4 coated-dipole 50 Ω VSWR sweep",
             xlabel="Frequency (MHz)", ylabel="VSWR", ylim=(1, 8))
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_power(path, low, lossy):
    figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                constrained_layout=True)
    for values, label in ((low, "tanδ = 0.002"), (lossy, "tanδ = 0.05")):
        frequencies = [row["frequency_mhz"] for row in values]
        axes[0].plot(frequencies,
                     [100*row["radiation_efficiency"] for row in values],
                     label=label, linewidth=1.8)
        axes[1].plot(frequencies,
                     [100*row["balance_defect_relative"] for row in values],
                     label=label, linewidth=1.8)
    axes[0].set(title="M4 lossy coating power accounting",
                ylabel="Radiation efficiency (%)")
    axes[1].set(xlabel="Frequency (MHz)",
                ylabel="Power-balance defect (%)")
    axes[1].axhline(0.3, color="0.3", linestyle="--", label="0.3% target")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_openems(path, nec, openems):
    figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                constrained_layout=True)
    for axis, key, label in ((axes[0], "resistance_ohm", "Resistance (Ω)"),
                             (axes[1], "reactance_ohm", "Reactance (Ω)")):
        axis.plot([row["frequency_mhz"] for row in nec],
                  [row[key] for row in nec], label="NEC M4", linewidth=1.8)
        axis.plot([row["frequency_mhz"] for row in openems],
                  [row[key] for row in openems], label="OpenEMS",
                  linewidth=1.6, linestyle="--")
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.25)
        axis.legend()
    axes[0].set_title("M4 constant-conductivity OpenEMS comparison")
    axes[1].set_xlabel("Frequency (MHz)")
    axes[1].axhline(0, color="0.2", linewidth=0.8)
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    decks = ROOT / "tests/coated-wire"
    with tempfile.TemporaryDirectory(prefix="nec2c-m4-report-") as temporary:
        work = Path(temporary)
        low = rows(run(args.binary, decks / "m4-low-loss-sweep.nec", work))
        lossy = rows(run(args.binary, decks / "m4-lossy-sweep.nec", work))
        equivalent = rows(run(args.binary,
                              decks / "m4-openems-equivalent-sweep.nec", work))
    record = {
        "schema": "nec2c-coated-wire-m4-acceptance",
        "schema_version": 1,
        "article": {"length_m": 0.475, "conductor_radius_m": 0.002,
                    "coating_outer_radius_m": 0.004, "er": 3.2},
        "sweeps": {"low_loss": low, "lossy": lossy,
                   "openems_equivalent": equivalent},
        "summary": {
            "low_loss_min_vswr": min(row["vswr_50_ohm"] for row in low),
            "lossy_min_vswr": min(row["vswr_50_ohm"] for row in lossy),
            "low_loss_min_efficiency": min(row["radiation_efficiency"] for row in low),
            "lossy_min_efficiency": min(row["radiation_efficiency"] for row in lossy),
            "maximum_balance_defect_relative": max(
                row["balance_defect_relative"] for row in low+lossy),
        },
    }
    if args.openems_dir:
        openems = openems_rows(args.openems_dir)
        far = json.loads((args.openems_dir / "far-field.json").read_text(
            encoding="utf-8"))
        nec_resonance = resonance(equivalent)
        openems_resonance = resonance(openems)
        accepted = interpolate(openems, far["frequency_mhz"], "accepted_w")
        nec_efficiency = (interpolate(equivalent,
                                     nec_resonance["frequency_mhz"],
                                     "radiated_w") /
                          interpolate(equivalent,
                                      nec_resonance["frequency_mhz"],
                                      "accepted_w"))
        record["openems"] = {
            "sweep": openems,
            "resonance": openems_resonance,
            "radiated_power_w": far["radiated_power_w"],
            "radiation_efficiency": far["radiated_power_w"]/accepted,
            "transient_end_energy_db": args.openems_end_energy_db,
        }
        record["nec2c_openems_comparison"] = {
            "nec2c_resonance": nec_resonance,
            "frequency_error_percent": 100*(
                nec_resonance["frequency_mhz"] /
                openems_resonance["frequency_mhz"]-1),
            "resistance_error_ohm": (nec_resonance["resistance_ohm"]-
                                     openems_resonance["resistance_ohm"]),
            "nec2c_radiation_efficiency": nec_efficiency,
            "efficiency_error_percentage_points": 100*(
                nec_efficiency-far["radiated_power_w"]/accepted),
        }
        plot_openems(args.output_dir / "m4-openems.png", equivalent, openems)
    output = args.output_dir / "m4-acceptance.json"
    output.write_text(json.dumps(record, indent=2, sort_keys=True)+"\n",
                      encoding="utf-8")
    plot_vswr(args.output_dir / "m4-vswr.png", low, lossy)
    plot_power(args.output_dir / "m4-power.png", low, lossy)
    print(output)


if __name__ == "__main__":
    main()
