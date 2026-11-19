#!/usr/bin/env python3
"""M9 Gate F: verify the frozen direct-table layout evidence."""

import json
import os
from pathlib import Path


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
REPORT = ROOT / "tests/nec4-similarity/m9-table-layout.json"


def main():
    document = json.loads(REPORT.read_text(encoding="utf-8"))
    if document["schema"] != "nec2c-nec4-m9-table-layout-qualification":
        raise AssertionError("unexpected table-layout evidence schema")
    if document["schema_version"] != 1:
        raise AssertionError("unexpected table-layout evidence version")
    if document["reference_evaluator"] != "validated-bessel-representation":
        raise AssertionError("table evidence does not name its direct evaluator")
    contour_evidence = ROOT / document["contour_equivalence_evidence"]
    if not contour_evidence.is_file():
        raise AssertionError("missing independent Bessel/Hankel evidence")

    target = document["target_maximum_error"]
    if target != 0.005 or document["worst_error"] > target:
        raise AssertionError(("table interpolation error", document["worst_error"]))
    expected_steps = {
        "1": {"rho": 0.02, "observer_z": 0.02, "source_depth": 0.02},
        "2": {"rho": 0.03, "observer_z": 0.03, "source_depth": 0.015},
        "3": {"rho": 0.03, "observer_z": 0.03, "source_depth": 0.03},
    }
    if document["maximum_steps_over_Lg"] != expected_steps:
        raise AssertionError("qualified table steps changed without new evidence")

    materials = document["materials"]
    if set(materials) != {"lossless-er4", "average-soil", "seawater"}:
        raise AssertionError("material coverage changed")
    if len(document["cases"]) != 39:
        raise AssertionError("expected 39 off-grid qualification probes")
    coverage = {(case["material"], case["region"])
                for case in document["cases"]}
    if coverage != {(name, region) for name in materials for region in (1, 2, 3)}:
        raise AssertionError("material/region coverage is incomplete")
    for case in document["cases"]:
        if case["maximum_null_aware_relative_error"] > target:
            raise AssertionError((case["material"], case["region"],
                                  case["maximum_null_aware_relative_error"]))
        for coordinate, bounds in zip(case["probe_over_Lg"],
                                      case["cell_bounds_over_Lg"]):
            if not bounds[0] < coordinate < bounds[1]:
                raise AssertionError(("probe is not strictly inside its cell",
                                      case["material"], case["region"],
                                      coordinate, bounds))
    if max(entry["payload_bytes"] for entry in materials.values()) > 64 * 1024**2:
        raise AssertionError("qualified raw table payload exceeds 64 MiB")

    print("M9 direct-table layout passed: 39 off-grid probes, "
          f"worst error {document['worst_error']:.6e}")
    for name, entry in materials.items():
        print(f"{name}: {entry['payload_bytes'] / 1024**2:.2f} MiB")


if __name__ == "__main__":
    main()
