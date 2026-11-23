#!/usr/bin/env python3
"""Verify or regenerate M11C interference-safe table evidence."""

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
DEFAULT_OUTPUT = ROOT / "tests/nec4-similarity/m11c-interference.json"
DEFAULT_PLOT = ROOT / "tests/nec4-similarity/m11c-interference.png"
REGIONS = ("grazing-sheet", "switch-volume")


def collect(output: Path, workers: int, budget: int) -> None:
    executable = ROOT / "tools/nec4_m11c_qualify"
    if not executable.exists():
        raise SystemExit("build tools/nec4_m11c_qualify before collecting M11C")
    retained = json.loads(output.read_text(encoding="utf-8"))
    results = []
    integer = ("baseline_cells", "baseline_unresolved", "baseline_fallback",
               "directional_budget", "directional_cells",
               "directional_unresolved", "directional_fallback",
               "both_fallback")
    floating = ("baseline_worst", "directional_worst",
                "maximum_cancellation", "baseline_seconds",
                "directional_seconds")
    for region in REGIONS:
        completed = subprocess.run(
            [str(executable), region, str(workers), str(budget)],
            cwd=ROOT, check=True, text=True, capture_output=True)
        row = next(csv.DictReader(io.StringIO(completed.stdout)))
        parsed = {name: int(row[name]) for name in integer}
        parsed.update({name: float(row[name]) for name in floating})
        results.append({
            "region": row["region"],
            "baseline_cells": parsed["baseline_cells"],
            "baseline_unresolved_leaves": parsed["baseline_unresolved"],
            "baseline_fallbacks": parsed["baseline_fallback"],
            "baseline_worst_error": parsed["baseline_worst"],
            "directional_cells": parsed["directional_cells"],
            "directional_unresolved_leaves": parsed["directional_unresolved"],
            "directional_fallbacks": parsed["directional_fallback"],
            "directional_worst_error": parsed["directional_worst"],
            "maximum_cancellation_ratio": parsed["maximum_cancellation"],
            "baseline_build_seconds": parsed["baseline_seconds"],
            "directional_build_seconds": parsed["directional_seconds"],
        })
    retained["directional_budget"] = budget
    retained["results"] = results
    output.write_text(json.dumps(retained, indent=2) + "\n", encoding="utf-8")


def plot(evidence: dict, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = evidence["results"]
    labels = [row["region"] for row in rows]
    x = list(range(len(rows)))
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    axes[0].bar([v - .18 for v in x],
                [100 * row["baseline_fallbacks"] / 256 for row in rows],
                width=.36, label="M11B tensor, 512 cells")
    axes[0].bar([v + .18 for v in x],
                [100 * row["directional_fallbacks"] / 256 for row in rows],
                width=.36, label="M11C directional, <=4096")
    axes[0].set_ylabel("Independent points requiring direct fallback (%)")
    axes[0].set_title("Difficult-domain coverage")
    axes[0].legend(fontsize=8)
    axes[1].bar(x, [row["directional_worst_error"] for row in rows],
                width=.55)
    axes[1].axhline(evidence["cell_tolerance"], color="black",
                    linestyle="--", linewidth=1, label="cell gate")
    axes[1].axhline(evidence["external_accuracy_gate"], color="gray",
                    linestyle=":", linewidth=1, label="external gate")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Worst independent scaled E error")
    axes[1].set_title("Accepted directional interpolation")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.set_xticks(x, labels)
        axis.grid(True, axis="y", alpha=.25)
    figure.tight_layout();figure.savefig(path, dpi=160);plt.close(figure)


def verify(output: Path) -> dict:
    evidence = json.loads(output.read_text(encoding="utf-8"))
    if evidence.get("schema") != "nec2c-nec4-m11c-interference-safe":
        raise AssertionError("M11C evidence schema")
    rows = evidence["results"]
    if {row["region"] for row in rows} != set(REGIONS):
        raise AssertionError("M11C region coverage")
    if evidence["spectral_partition"]["worst_recomposition_relative_error"] > 3e-13:
        raise AssertionError("M11C exact spectral recomposition")
    for row in rows:
        if row["directional_cells"] > evidence["directional_budget"]:
            raise AssertionError(("M11C cell budget", row))
        if row["directional_unresolved_leaves"] or row["directional_fallbacks"]:
            raise AssertionError(("M11C unresolved directional cell", row))
        if row["directional_worst_error"] > evidence["cell_tolerance"]:
            raise AssertionError(("M11C independent accuracy", row))
        if row["maximum_cancellation_ratio"] < 100:
            raise AssertionError(("M11C interference witness too weak", row))
    print("M11C bounded directional interference tables passed")
    for row in rows:
        print(f"  {row['region']}: {row['directional_cells']} cells, "
              f"{row['directional_worst_error']:.3e} worst error, "
              f"{row['maximum_cancellation_ratio']:.1f}x cancellation")
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
    if arguments.plot_existing or arguments.collect:
        plot(evidence, arguments.plot)


if __name__ == "__main__":
    main()
