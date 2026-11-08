#!/usr/bin/env python3
"""Generate the compact M3 coated-wire acceptance record and overlay plots."""

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
    return parser.parse_args()


def run(binary, deck, directory):
    result = directory / f"{deck.stem}.json"
    report = directory / f"{deck.stem}.out"
    completed = subprocess.run([
        str(binary), "--backend", "builtin", "--jobs", "1", "--threads", "1",
        "--result-json", str(result), "-i", str(deck), "-o", str(report),
    ], cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(f"{deck} failed\n{completed.stdout}\n{completed.stderr}")
    return json.loads(result.read_text())


def impedance_rows(result):
    rows = []
    for execution in result["jobs"][0]["executions"]:
        impedance = execution["feeds"][0]["impedance_ohm"]
        rows.append({"frequency_mhz": execution["frequency_mhz"],
                     "resistance_ohm": impedance["re"],
                     "reactance_ohm": impedance["im"]})
    return rows


def resonance(rows):
    for left, right in zip(rows, rows[1:]):
        if left["reactance_ohm"] <= 0 <= right["reactance_ohm"] or \
           right["reactance_ohm"] <= 0 <= left["reactance_ohm"]:
            fraction = -left["reactance_ohm"] / (
                right["reactance_ohm"]-left["reactance_ohm"])
            return {
                "frequency_mhz": left["frequency_mhz"] + fraction*(
                    right["frequency_mhz"]-left["frequency_mhz"]),
                "resistance_ohm": left["resistance_ohm"] + fraction*(
                    right["resistance_ohm"]-left["resistance_ohm"]),
            }
    raise RuntimeError("sweep does not bracket resonance")


def currents(result):
    values = []
    execution = result["jobs"][0]["executions"][0]
    for row in execution["wire_currents"]:
        current = row["current_a"]
        values.append({"z_m": row["center_m"][2],
                       "magnitude_a": math.hypot(current["re"], current["im"])})
    peak = max(row["magnitude_a"] for row in values)
    for row in values:
        row["relative_magnitude"] = row["magnitude_a"]/peak
    return values


def pattern(result):
    samples = result["jobs"][0]["executions"][0]["far_fields"][0]["samples"]
    values = []
    for row in samples:
        et, ep = row["r_e_theta"], row["r_e_phi"]
        magnitude = math.sqrt(et["re"]**2+et["im"]**2+
                              ep["re"]**2+ep["im"]**2)
        values.append({"theta_deg": row["theta_deg"], "magnitude": magnitude})
    peak = max(row["magnitude"] for row in values)
    for row in values:
        row["relative_db"] = (20*math.log10(row["magnitude"]/peak)
                              if row["magnitude"] else -300.0)
    return values


def openems_reference():
    cases = json.loads((ROOT / "tests/coated-wire/m1-openems-caps.json").read_text())["cases"]
    selected = next(case for case in cases
                    if case["coating_free_end"] == "flush" and
                    case["coating_radial_cells"] == 8)
    phi_zero = [sample for sample in selected["far_field"]["samples"]
                if sample["phi_deg"] == 0]
    peak = max(sample["relative_db"] for sample in phi_zero)
    return selected, [{"theta_deg": sample["theta_deg"],
                       "relative_db": sample["relative_db"]-peak}
                      for sample in phi_zero]


def plot_impedance(path, bare, coated, bare_res, coated_res, reference):
    figure, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                constrained_layout=True)
    for axis, key, label in ((axes[0], "resistance_ohm", "Resistance (Ω)"),
                             (axes[1], "reactance_ohm", "Reactance (Ω)")):
        axis.plot([row["frequency_mhz"] for row in bare], [row[key] for row in bare],
                  label="NEC bare", linewidth=1.8)
        axis.plot([row["frequency_mhz"] for row in coated], [row[key] for row in coated],
                  label="NEC coated", linewidth=1.8)
        axis.axvline(reference["bare_frequency_mhz"], color="0.35", linestyle=":",
                     label="OpenEMS bare resonance" if axis is axes[0] else None)
        axis.axvline(reference["coated_frequency_mhz"], color="0.55", linestyle="--",
                     label="OpenEMS coated resonance" if axis is axes[0] else None)
        axis.grid(True, alpha=0.25)
        axis.set_ylabel(label)
    axes[0].scatter([bare_res["frequency_mhz"], coated_res["frequency_mhz"]],
                    [bare_res["resistance_ohm"], coated_res["resistance_ohm"]], s=28)
    axes[0].legend(ncol=2, fontsize=9)
    axes[0].set_title("M3 dipole impedance and resonance shift")
    axes[1].axhline(0, color="0.2", linewidth=0.8)
    axes[1].set_xlabel("Frequency (MHz)")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_currents(path, bare, coated):
    figure, axis = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    axis.plot([row["z_m"] for row in bare], [row["relative_magnitude"] for row in bare],
              label="NEC bare at resonance", linewidth=1.8)
    axis.plot([row["z_m"] for row in coated], [row["relative_magnitude"] for row in coated],
              label="NEC coated at resonance", linewidth=1.8)
    axis.set(title="M3 normalized current distribution", xlabel="Axial position (m)",
             ylabel="Relative current magnitude")
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_pattern(path, nec, reference):
    figure, axis = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    axis.plot([row["theta_deg"] for row in nec], [row["relative_db"] for row in nec],
              label="NEC coated", linewidth=1.8)
    axis.plot([row["theta_deg"] for row in reference],
              [row["relative_db"] for row in reference], label="OpenEMS coated",
              linewidth=1.6, linestyle="--")
    axis.set(title="M3 coated-dipole normalized pattern", xlabel="Theta (degrees)",
             ylabel="Relative field (dB)", ylim=(-40, 1))
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    decks = ROOT / "tests/coated-wire"
    with tempfile.TemporaryDirectory(prefix="nec2c-m3-report-") as temporary:
        work = Path(temporary)
        bare_sweep_result = run(args.binary, decks / "m3-bare-sweep.nec", work)
        coated_sweep_result = run(args.binary, decks / "m3-coated-sweep.nec", work)
        bare_result = run(args.binary, decks / "m3-bare-dipole.nec", work)
        coated_result = run(args.binary, decks / "m3-coated-dipole.nec", work)

    bare_sweep = impedance_rows(bare_sweep_result)
    coated_sweep = impedance_rows(coated_sweep_result)
    bare_res = resonance(bare_sweep)
    coated_res = resonance(coated_sweep)
    bare_current = currents(bare_result)
    coated_current = currents(coated_result)
    coated_pattern = pattern(coated_result)
    openems, openems_pattern = openems_reference()
    reference = {
        "bare_frequency_mhz": 295.815987,
        "bare_resistance_ohm": 72.797845,
        "coated_frequency_mhz": openems["resonance_mhz"],
        "coated_resistance_ohm": openems["resistance_at_resonance_ohm"],
        "coated_directivity_dbi": openems["far_field"]["directivity_db"],
    }
    record = {
        "schema": "nec2c-coated-wire-m3-acceptance",
        "schema_version": 1,
        "article": {"length_m": 0.475, "conductor_radius_m": 0.002,
                    "coating_outer_radius_m": 0.004, "er": 3.2},
        "nec2c": {"bare_resonance": bare_res, "coated_resonance": coated_res,
                  "resonance_shift_percent":
                      100*(coated_res["frequency_mhz"]/bare_res["frequency_mhz"]-1)},
        "openems": reference,
        "errors": {
            "bare_frequency_percent": 100*(bare_res["frequency_mhz"]/
                                             reference["bare_frequency_mhz"]-1),
            "coated_frequency_percent": 100*(coated_res["frequency_mhz"]/
                                               reference["coated_frequency_mhz"]-1),
            "bare_resistance_ohm": bare_res["resistance_ohm"]-
                                   reference["bare_resistance_ohm"],
            "coated_resistance_ohm": coated_res["resistance_ohm"]-
                                      reference["coated_resistance_ohm"],
        },
        "sweeps": {"bare": bare_sweep, "coated": coated_sweep},
        "currents": {"bare": bare_current, "coated": coated_current},
        "patterns": {"nec2c_coated": coated_pattern,
                     "openems_coated": openems_pattern},
    }
    (args.output_dir / "m3-acceptance.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n")
    plot_impedance(args.output_dir / "m3-impedance.png", bare_sweep,
                   coated_sweep, bare_res, coated_res, reference)
    plot_currents(args.output_dir / "m3-current.png", bare_current, coated_current)
    plot_pattern(args.output_dir / "m3-pattern.png", coated_pattern, openems_pattern)
    print(args.output_dir / "m3-acceptance.json")


if __name__ == "__main__":
    main()
