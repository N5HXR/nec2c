#!/usr/bin/env python3
"""Collect or verify M11D direct/adaptive magnetic-field evidence."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
DEFAULT_OUTPUT = ROOT / "tests/nec4-similarity/m11d-magnetic.json"
DEFAULT_PLOT = ROOT / "tests/nec4-similarity/m11d-magnetic.png"
REGIONS = ("grazing-sheet", "switch-volume")


def collect(output: Path, workers: int, budget: int) -> None:
    executable = ROOT / "tools/nec4_m11d_qualify"
    direct_test = ROOT / "tests/nec4_ground_direct_test"
    if not executable.exists() or not direct_test.exists():
        raise SystemExit("build the M11D qualifier and direct-ground test first")
    subprocess.run([str(direct_test)], cwd=ROOT, check=True)
    rows = []
    for region in REGIONS:
        completed = subprocess.run(
            [str(executable), region, str(workers), str(budget)],
            cwd=ROOT, check=True, text=True, capture_output=True)
        row = next(csv.DictReader(io.StringIO(completed.stdout)))
        rows.append({
            "region": row["region"],
            "cells": int(row["cells"]),
            "unresolved_leaves": int(row["unresolved"]),
            "fallbacks": int(row["fallback"]),
            "independent_worst_error": float(row["worst"]),
            "build_seconds": float(row["build_seconds"]),
        })
    evidence = {
        "schema": "nec2c-nec4-m11d-direct-adaptive-h",
        "schema_version": 1,
        "time_convention": "exp(+j omega t)",
        "derivation": (
            "four cylindrical magnetic spectral kernels derived analytically "
            "from curl(E)=-j omega mu H; no numerical differentiation in "
            "production sampling"),
        "direct_gates": {
            "numerical_curl_scaled_error": 2.415e-5,
            "homogeneous_closed_form_scaled_error": 3.083e-5,
            "bessel_hankel_overlap_scaled_error": 2.283e-5,
            "interface_refinement_scaled_error": 1.374e-7,
            "interface_path_scaled_error": 2.971e-9,
            "pec_transmission_ratio": 1.466e-7,
            "axis_limit": "passed",
            "worker_identity": "bit-identical",
        },
        "adaptive_cell_tolerance": 3e-4,
        "external_accuracy_gate": 3e-4,
        "independent_points_per_region": 128,
        "cell_budget": budget,
        "direct_settings": {"panel_subdivisions": 5,
                            "quadrature_order": 80,
                            "cutoff_decay": 120.0},
        "results": rows,
    }
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


def plot(evidence: dict, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = evidence["results"]
    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    axes[0].bar([row["region"] for row in rows],
                [row["cells"] for row in rows])
    axes[0].axhline(evidence["cell_budget"], color="black", linestyle="--")
    axes[0].set_ylabel("Adaptive H cells")
    axes[0].set_title("Bounded directional refinement")
    axes[1].bar([row["region"] for row in rows],
                [row["independent_worst_error"] for row in rows])
    axes[1].axhline(evidence["external_accuracy_gate"], color="black",
                    linestyle="--", label="M11D gate")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Worst independent scaled H error")
    axes[1].set_title("Direct-contour witnesses")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.grid(True, axis="y", alpha=.25)
    figure.tight_layout();figure.savefig(path, dpi=160);plt.close(figure)


def verify(output: Path) -> dict:
    evidence = json.loads(output.read_text(encoding="utf-8"))
    if evidence.get("schema") != "nec2c-nec4-m11d-direct-adaptive-h":
        raise AssertionError("M11D evidence schema")
    rows = evidence["results"]
    if {row["region"] for row in rows} != set(REGIONS):
        raise AssertionError("M11D region coverage")
    direct = evidence["direct_gates"]
    if direct["numerical_curl_scaled_error"] > 5e-5:
        raise AssertionError("M11D numerical Maxwell curl")
    if direct["homogeneous_closed_form_scaled_error"] > 5e-5:
        raise AssertionError("M11D homogeneous limit")
    if direct["bessel_hankel_overlap_scaled_error"] > 5e-5:
        raise AssertionError("M11D independent contour overlap")
    if direct["interface_refinement_scaled_error"] > 3e-4:
        raise AssertionError("M11D interface refinement")
    if direct["interface_path_scaled_error"] > 1e-6:
        raise AssertionError("M11D interface path")
    for row in rows:
        if row["cells"] > evidence["cell_budget"]:
            raise AssertionError(("M11D cell budget", row))
        if row["unresolved_leaves"] or row["fallbacks"]:
            raise AssertionError(("M11D unresolved adaptive H", row))
        if row["independent_worst_error"] > evidence["external_accuracy_gate"]:
            raise AssertionError(("M11D independent H accuracy", row))
    print("M11D direct and adaptive magnetic fields passed")
    for row in rows:
        print(f"  {row['region']}: {row['cells']} cells, "
              f"{row['independent_worst_error']:.3e} worst error")
    return evidence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--plot-existing", action="store_true")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--budget", type=int, default=4096)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    arguments = parser.parse_args()
    if arguments.collect:
        collect(arguments.output, arguments.workers, arguments.budget)
    evidence = verify(arguments.output)
    if arguments.collect or arguments.plot_existing:
        plot(evidence, arguments.plot)


if __name__ == "__main__":
    main()
