#!/usr/bin/env python3
"""M9 Gate E oracle-backed original-NEC4 grazing warning gate."""

import json
import os
from pathlib import Path


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()


def main():
    path = ROOT / "tests/nec4-similarity/m9-grazing-warning-gate.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    if report["schema"] != "nec2c-nec4-m9-grazing-warning-gate":
        raise AssertionError("grazing warning schema")
    if not report["warning_is_diagnostic_not_correction"]:
        raise AssertionError("warning must not change original NEC4 result")
    rows = report["rows"]
    if len(rows) != 12:
        raise AssertionError("qualification lattice")
    for row in rows:
        expected = row["direct_ray_angle_degrees"] <= 15.0
        if row["warning"] != expected:
            raise AssertionError(("warning ownership", row))
    witnessed = [row for row in rows if row["warning"]]
    if not all(row["method_errors"]["second_order"]
               ["maximum_component_relative_error"] > 0.5
               for row in witnessed):
        raise AssertionError("warning envelope no longer exposes severe error")
    unwarned = [row for row in rows if not row["warning"]]
    if not all(row["method_errors"]["second_order"]["vector_relative_error"] < .12
               for row in unwarned):
        raise AssertionError("warning boundary leaves a severe ordinary-saddle error")
    print("M9 grazing warning is oracle-backed and preserves legacy behavior")


if __name__ == "__main__":
    main()
