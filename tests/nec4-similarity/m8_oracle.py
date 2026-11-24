#!/usr/bin/env python3
"""M8 Gate D1: diagnostic point/line oracle interface."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()


def run(*arguments):
    with tempfile.NamedTemporaryFile(suffix=".json") as output:
        subprocess.run(
            [sys.executable, str(ROOT / "tools" / "nec4_m8_oracle.py"),
             *arguments, "--output", output.name],
            cwd=ROOT, check=True, text=True, stdout=subprocess.PIPE)
        return json.load(output)


def main():
    common = ("--method", "bessel", "--precision", "55", "--cutoff", "65",
              "--subdivisions", "2", "--order", "40")
    point = run("point", *common)
    if point["schema"] != "nec2c-nec4-m8-direct-oracle":
        raise AssertionError("point schema")
    if point["status"] != "reference-only; not connected to production matrix fill":
        raise AssertionError("production boundary")
    if len(point["points"]) != 1 or point["points"][0]["pairing"] != "12":
        raise AssertionError("point pairing")
    if float(point["points"][0]["worst_refinement_error"]) > 1e-9:
        raise AssertionError("point refinement")

    line = run("line", "--line-count", "5", *common)
    # z=0 is intentionally omitted rather than silently assigned to one side.
    if [entry["pairing"] for entry in line["points"]] != ["11", "11", "12", "12"]:
        raise AssertionError("line pairing sequence")
    if any(float(entry["worst_refinement_error"]) > 1e-9
           for entry in line["points"]):
        raise AssertionError("line refinement")
    if any(len(entry["cartesian_electric_dyadic"]) != 3
           for entry in line["points"]):
        raise AssertionError("dyadic shape")
    print("M8 Gate D1 diagnostic point/line oracle passed with explicit refinement errors")


if __name__ == "__main__":
    main()
