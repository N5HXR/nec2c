#!/usr/bin/env python3
"""Collect or verify the M11A difficult-domain E-field measurement contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
DEFAULT_OUTPUT = ROOT / "tests/nec4-similarity/m11a-measurement.json"
DEFAULT_PLOT = ROOT / "tests/nec4-similarity/m11a-measurement.png"
COMPONENTS = ("rho_vertical", "z_vertical", "rho_horizontal",
              "phi_horizontal", "z_horizontal")
MATERIALS = ("average-soil", "low-loss-water", "seawater",
             "lossless-water")


def number(value: str) -> float:
    return float(value)


def optional_number(value: str) -> float | None:
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def collect(output: Path, plot: Path) -> None:
    executable = ROOT / "tools/nec4_m11_measure"
    if not executable.exists():
        raise SystemExit("build tools/nec4_m11_measure before collecting M11A")
    with tempfile.TemporaryDirectory(prefix="nec2c-m11a-") as name:
        csv_path = Path(name) / "measurement.csv"
        subprocess.run([str(executable), "--output", str(csv_path)],
                       cwd=ROOT, check=True)
        rows = []
        with csv_path.open(newline="", encoding="utf-8") as stream:
            for source in csv.DictReader(stream):
                row = {
                    "material": source["material"],
                    "family": source["family"],
                    "rho_lambda": number(source["rho_lambda"]),
                    "observer_z_lambda": number(source["z_lambda"]),
                    "source_depth_lg": number(source["depth_lg"]),
                    "production_method": source["method"],
                    "grazing_angle_degrees": number(
                        source["grazing_angle_deg"]),
                    "grazing_warning": bool(int(source["grazing_risk"])),
                    "medium_ok": bool(int(source["medium_ok"])),
                    "fine_ok": bool(int(source["fine_ok"])),
                    "alternate_ok": bool(int(source["alternate_ok"])),
                    "medium_to_fine_scaled_error": number(
                        source["medium_fine_error"]),
                    "alternate_to_fine_scaled_error": optional_number(
                        source["alternate_fine_error"]),
                    "production_to_fine_scaled_error": number(
                        source["production_fine_error"]),
                    "wall_seconds": {
                        "medium": number(source["medium_seconds"]),
                        "fine": number(source["fine_seconds"]),
                        "alternate": number(source["alternate_seconds"]),
                    },
                    "fine_direct_field": {
                        COMPONENTS[index]: {
                            "re": number(source[f"fine_c{index}_real"]),
                            "im": number(source[f"fine_c{index}_imag"]),
                        } for index in range(len(COMPONENTS))
                    },
                }
                rows.append(row)
    selected = [row["medium_to_fine_scaled_error"] for row in rows
                if math.isfinite(row["medium_to_fine_scaled_error"])]
    alternate = [row["alternate_to_fine_scaled_error"] for row in rows
                 if row["alternate_to_fine_scaled_error"] is not None]
    direct = [row["production_to_fine_scaled_error"] for row in rows
              if row["production_method"] != "asymptotic"]
    asymptotic = [row["production_to_fine_scaled_error"] for row in rows
                  if row["production_method"] == "asymptotic"]
    evidence = {
        "schema": "nec2c-nec4-m11a-difficult-domain",
        "schema_version": 1,
        "purpose": (
            "condition the public direct-contour oracle and map the accepted "
            "original-NEC4 production path before adaptive implementation"),
        "normalization": (
            "maximum component absolute difference divided by maximum "
            "fine-direct component magnitude"),
        "contour_sequences": {
            "medium": {"panel_subdivisions": 3,
                       "quadrature_order": 48, "cutoff_decay": 80.0},
            "fine": {"panel_subdivisions": 5,
                     "quadrature_order": 80, "cutoff_decay": 120.0},
            "alternate": (
                "fine settings with Bessel and Hankel contour ownership "
                "exchanged where both representations are available"),
        },
        "provisional_scaled_error_target": 1e-4,
        "summary": {
            "witness_count": len(rows),
            "selected_oracle_worst": max(selected),
            "alternate_oracle_worst": max(alternate),
            "accepted_direct_path_worst": max(direct),
            "current_asymptotic_worst": max(asymptotic),
            "current_asymptotic_median": sorted(asymptotic)[
                len(asymptotic) // 2],
            "current_asymptotic_above_target": sum(
                value > 1e-4 for value in asymptotic),
            "current_asymptotic_count": len(asymptotic),
        },
        "witnesses": rows,
    }
    output.write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n",
                      encoding="utf-8")
    plot_evidence(evidence, plot)


def plot_evidence(evidence: dict, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    switch, angle = axes
    colors = dict(zip(MATERIALS, ("tab:blue", "tab:orange", "tab:red",
                                  "tab:green")))
    for material in MATERIALS:
        points = [row for row in evidence["witnesses"]
                  if row["material"] == material and
                  row["family"] == "switch"]
        switch.semilogy([row["rho_lambda"] for row in points],
                        [row["production_to_fine_scaled_error"]
                         for row in points], marker="o", label=material,
                        color=colors[material])
        points = [row for row in evidence["witnesses"]
                  if row["material"] == material and
                  row["family"] == "grazing-angle"]
        angle.semilogy([row["grazing_angle_degrees"] for row in points],
                       [row["production_to_fine_scaled_error"]
                        for row in points], marker="o", label=material,
                       color=colors[material])
    for axis in axes:
        axis.axhline(evidence["provisional_scaled_error_target"],
                     color="black", linestyle="--", linewidth=1,
                     label="M11 target")
        axis.grid(True, which="both", alpha=.25)
        axis.set_ylabel("Scaled E-tensor error")
    switch.axvline(3., color="0.4", linestyle=":", linewidth=1)
    switch.set_xlabel("Horizontal separation (wavelengths)")
    switch.set_title("Current direct/asymptotic switch")
    angle.set_xlabel("Interface grazing angle (degrees)")
    angle.set_title("Current asymptotic path at 4 wavelengths")
    angle.legend(fontsize=8, loc="best")
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def verify(output: Path) -> None:
    evidence = json.loads(output.read_text(encoding="utf-8"))
    if evidence.get("schema") != "nec2c-nec4-m11a-difficult-domain":
        raise AssertionError("M11A evidence schema")
    rows = evidence["witnesses"]
    if len(rows) != 80 or {row["material"] for row in rows} != set(MATERIALS):
        raise AssertionError("M11A witness coverage")
    target = evidence["provisional_scaled_error_target"]
    summary = evidence["summary"]
    if summary["selected_oracle_worst"] > target:
        raise AssertionError("M11A selected-contour oracle conditioning")
    if summary["alternate_oracle_worst"] > target:
        raise AssertionError("M11A alternate-contour oracle conditioning")
    if summary["accepted_direct_path_worst"] > target:
        raise AssertionError("M11A accepted direct path")
    if summary["current_asymptotic_above_target"] < 20:
        raise AssertionError("M11A difficult-domain sensitivity coverage")
    if any(not row["fine_ok"] for row in rows):
        raise AssertionError("M11A fine direct-contour failure")
    print("M11A retained difficult-domain oracle and production map passed")
    print(f"  selected oracle worst: {summary['selected_oracle_worst']:.3e}")
    print(f"  alternate oracle worst: {summary['alternate_oracle_worst']:.3e}")
    print(f"  current asymptotic worst: {summary['current_asymptotic_worst']:.3e}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--plot-existing", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    arguments = parser.parse_args()
    if arguments.collect:
        collect(arguments.output, arguments.plot)
    verify(arguments.output)
    if arguments.plot_existing:
        plot_evidence(json.loads(arguments.output.read_text(encoding="utf-8")),
                      arguments.plot)


if __name__ == "__main__":
    main()
