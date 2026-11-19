#!/usr/bin/env python3
"""M9 production CLI, cache, direct-oracle, and warning acceptance gates."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
HERE = ROOT / "tests" / "nec4-similarity"


def run(source, directory, stem, model="grid", timing=False):
    report = directory / f"{stem}.out"
    result = directory / f"{stem}.json"
    command = [str(BINARY), "--wire-engine", "nec4", "--ground-engine", "nec4",
               "--nec4-ground-model", model, "--nec4-ground-cache",
               str(directory / "cache"), "--threads", "auto", "--result-json",
               str(result)]
    timing_path = directory / f"{stem}-timing.json"
    if timing:
        command += ["--timing-json", str(timing_path)]
    command += ["-i", str(source), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                               timeout=240)
    if completed.returncode:
        raise AssertionError((command, completed.returncode, completed.stderr,
                              report.read_text(errors="replace")[-4000:]
                              if report.exists() else ""))
    document = json.loads(result.read_text(encoding="utf-8"))
    timing_document = json.loads(timing_path.read_text()) if timing else None
    return document, report.read_text(errors="replace"), timing_document


def execution(document):
    return document["jobs"][0]["executions"][0]


def complex_pair(value):
    return complex(value["re"], value["im"])


def norm_delta(left, right):
    scale = math.sqrt(sum(abs(value) ** 2 for value in right))
    error = math.sqrt(sum(abs(a-b) ** 2 for a, b in zip(left, right)))
    return error / max(scale, 1e-300)


def phase(timing, name):
    rows = [row for row in timing["phases"] if row["name"] == name]
    if len(rows) != 1:
        raise AssertionError((name, timing["phases"]))
    return rows[0]["seconds"]


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m9-production-") as name:
        directory = Path(name)
        source = HERE / "m9-separated-media.nec"
        grid, cold_report, grid_timing = run(source, directory, "grid-cold",
                                             timing=True)
        repeat, warm_report, _ = run(source, directory, "grid-warm")
        direct, direct_report, direct_timing = run(source, directory, "direct",
                                                   model="direct", timing=True)
        if "0/4 cached tables" not in cold_report:
            raise AssertionError("cold grid build did not generate four tables")
        if "4/4 cached tables" not in warm_report:
            raise AssertionError("warm grid run did not load all four tables")
        if "3/3 cached tables" not in direct_report:
            raise AssertionError("direct mode did not reuse the three near tables")
        for key in ("feeds", "wire_currents", "patch_currents", "power",
                    "far_fields", "near_fields"):
            if execution(grid)[key] != execution(repeat)[key]:
                raise AssertionError(("warm-cache solve changed machine results",
                                      key))

        grid_execution = execution(grid)
        direct_execution = execution(direct)
        grid_z = complex_pair(grid_execution["feeds"][0]["impedance_ohm"])
        direct_z = complex_pair(direct_execution["feeds"][0]["impedance_ohm"])
        impedance_delta = abs(grid_z-direct_z)/max(abs(direct_z), 1e-300)
        grid_currents = [complex_pair(row["current_a"])
                         for row in grid_execution["wire_currents"]]
        direct_currents = [complex_pair(row["current_a"])
                           for row in direct_execution["wire_currents"]]
        current_delta = norm_delta(grid_currents, direct_currents)
        if impedance_delta > 2e-5 or current_delta > 2e-5:
            raise AssertionError(("grid/direct disagreement", impedance_delta,
                                  current_delta))

        grid_fill = phase(grid_timing, "matrix_fill")
        direct_fill = phase(direct_timing, "matrix_fill")
        if not (grid_fill > 0 and direct_fill > 0):
            raise AssertionError("missing positive matrix timings")

        grazing, _, _ = run(HERE / "m9-grazing-cli.nec", directory,
                            "grazing")
        warnings = [row for row in grazing["diagnostics"]
                    if row["code"] == "NEC4_M9_GRAZING_INTERFERENCE"]
        if len(warnings) != 1 or warnings[0]["severity"] != "warning":
            raise AssertionError(("grazing warning", grazing["diagnostics"]))

        evidence = {
            "schema": "nec2c-nec4-m9-production-acceptance",
            "schema_version": 1,
            "grid_direct": {
                "impedance_relative_delta": impedance_delta,
                "current_l2_relative_delta": current_delta,
                "grid_matrix_fill_seconds": grid_fill,
                "direct_matrix_fill_seconds": direct_fill,
                "matrix_fill_speedup": direct_fill / grid_fill,
            },
            "cache": {"cold_misses": 4, "warm_hits": 4,
                      "warm_results_exact": True},
            "grazing_warning": warnings[0],
        }
        (HERE / "m9-production-gate.json").write_text(
            json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
        print(f"M9 production grid/direct: impedance {impedance_delta:.3e}, "
              f"currents {current_delta:.3e} relative delta")
        print(f"M9 cached model-region matrix fill: {direct_fill/grid_fill:.2f}x "
              "versus direct contour on the witness")
        print("M9 cache identity and grazing-warning CLI gates passed")


if __name__ == "__main__":
    main()
