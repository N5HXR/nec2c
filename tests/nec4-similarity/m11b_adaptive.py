#!/usr/bin/env python3
"""Collect or verify M11B adaptive total-E qualification evidence."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
DEFAULT_OUTPUT = ROOT / "tests/nec4-similarity/m11b-adaptive.json"
DEFAULT_PLOT = ROOT / "tests/nec4-similarity/m11b-adaptive.png"
MATERIALS = ("average-soil", "low-loss-water", "seawater",
             "lossless-water")
REGIONS = ("switch-volume", "interface-line", "grazing-sheet",
           "axis-volume")


def convert(source: dict[str, str]) -> dict:
    integer = ("cells", "leaves", "accepted_leaves", "unresolved_leaves",
               "deepest_level", "fitted_samples", "withheld_samples",
               "direct_failures", "adaptive_fallbacks")
    floating = ("worst_leaf_error", "worst_encountered_error",
                "independent_worst_error", "build_seconds",
                "evaluation_seconds")
    row = {"material": source["material"], "region": source["region"]}
    row.update({name: int(source[name]) for name in integer})
    row.update({name: float(source[name]) for name in floating})
    row["independent_witnesses"] = 128
    row["accelerated_witnesses"] = 128 - row["adaptive_fallbacks"]
    row["fallback_fraction"] = row["adaptive_fallbacks"] / 128
    return row


def collect(output: Path, plot: Path, workers: int) -> None:
    executable = ROOT / "tools/nec4_m11b_qualify"
    if not executable.exists():
        raise SystemExit("build tools/nec4_m11b_qualify before collecting M11B")
    with tempfile.TemporaryDirectory(prefix="nec2c-m11b-") as name:
        csv_path = Path(name) / "qualification.csv"
        subprocess.run([str(executable), "--workers", str(workers),
                        "--output", str(csv_path)], cwd=ROOT, check=True)
        with csv_path.open(newline="", encoding="utf-8") as stream:
            rows = [convert(row) for row in csv.DictReader(stream)]
    evidence = {
        "schema": "nec2c-nec4-m11b-adaptive-total-e",
        "schema_version": 1,
        "representation": (
            "tensor-cubic adaptive cells over total E after the accepted "
            "F2 phase and spreading transform"),
        "fit_nodes": [0.0, 1 / 3, 2 / 3, 1.0],
        "withheld_nodes": [1 / 6, 1 / 2, 5 / 6],
        "cell_tolerance": 5e-5,
        "external_accuracy_gate": 1e-4,
        "independent_points_per_region": 128,
        "cell_budget": {"ordinary": 512, "interface-line": 2048},
        "direct_settings": {"panel_subdivisions": 5,
                            "quadrature_order": 80,
                            "cutoff_decay": 120.0},
        "fallback_policy": (
            "unresolved or over-budget leaves return no interpolated value; "
            "the future runtime must evaluate the direct contour"),
        "results": rows,
    }
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    plot_evidence(evidence, plot)


def plot_evidence(evidence: dict, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    errors, fallback = axes
    x = list(range(len(REGIONS)))
    width = .19
    for offset, material in enumerate(MATERIALS):
        rows = [next(row for row in evidence["results"]
                     if row["material"] == material and row["region"] == region)
                for region in REGIONS]
        positions = [value + (offset - 1.5) * width for value in x]
        errors.bar(positions, [row["independent_worst_error"] for row in rows],
                   width=width, label=material)
        fallback.bar(positions, [100 * row["fallback_fraction"] for row in rows],
                     width=width, label=material)
    errors.axhline(evidence["external_accuracy_gate"], color="black",
                   linestyle="--", linewidth=1)
    errors.set_yscale("log")
    errors.set_ylabel("Worst independent scaled E error")
    errors.set_title("Accuracy where interpolation was accepted")
    fallback.set_ylabel("Direct fallback witnesses (%)")
    fallback.set_title("Representation rejection is explicit")
    for axis in axes:
        axis.set_xticks(x, REGIONS, rotation=20, ha="right")
        axis.grid(True, axis="y", alpha=.25)
    errors.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def verify(output: Path) -> None:
    evidence = json.loads(output.read_text(encoding="utf-8"))
    if evidence.get("schema") != "nec2c-nec4-m11b-adaptive-total-e":
        raise AssertionError("M11B evidence schema")
    rows = evidence["results"]
    if len(rows) != len(MATERIALS) * len(REGIONS):
        raise AssertionError("M11B coverage")
    if {(row["material"], row["region"]) for row in rows} != {
            (material, region) for material in MATERIALS for region in REGIONS}:
        raise AssertionError("M11B material/region matrix")
    gate = evidence["external_accuracy_gate"]
    for row in rows:
        if row["direct_failures"]:
            raise AssertionError(("M11B direct oracle failure", row))
        if row["accelerated_witnesses"] and row["independent_worst_error"] > gate:
            raise AssertionError(("M11B independent accuracy", row))
        if row["worst_leaf_error"] > evidence["cell_tolerance"] * 1.000001:
            raise AssertionError(("M11B accepted leaf accuracy", row))
        budget = evidence["cell_budget"][
            "interface-line" if row["region"] == "interface-line" else "ordinary"]
        if row["cells"] > budget:
            raise AssertionError(("M11B cell budget", row))
    for row in rows:
        if row["material"] != "lossless-water" and row["adaptive_fallbacks"]:
            raise AssertionError(("M11B unexpected lossy fallback", row))
    print("M11B retained adaptive total-E accuracy and fallback gates passed")
    print("  worst accepted independent error: "
          f"{max(row['independent_worst_error'] for row in rows):.3e}")
    lossless = [row for row in rows if row["material"] == "lossless-water"]
    print("  lossless fallback regions: " + ", ".join(
        f"{row['region']}={100 * row['fallback_fraction']:.1f}%"
        for row in lossless if row["adaptive_fallbacks"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--plot-existing", action="store_true")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    arguments = parser.parse_args()
    if arguments.collect:
        collect(arguments.output, arguments.plot, arguments.workers)
    verify(arguments.output)
    if arguments.plot_existing:
        plot_evidence(json.loads(arguments.output.read_text(encoding="utf-8")),
                      arguments.plot)


if __name__ == "__main__":
    main()
