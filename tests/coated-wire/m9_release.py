#!/usr/bin/env python3
"""M9 release compatibility, order, evidence, and performance gates."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT/"nec2c")).resolve()
HERE = ROOT/"tests/coated-wire"


def run(source, directory, name, backend):
    result = directory/f"{name}.json"
    report = directory/f"{name}.out"
    completed = subprocess.run(
        (str(BINARY), "--backend", backend, "--jobs", "1", "--threads", "1",
         "--result-json", str(result), "-i", str(source), "-o", str(report)),
        cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stdout+completed.stderr)
    return json.loads(result.read_text()), report.read_text()


def execution_order(result):
    execution = result["jobs"][0]["executions"][0]
    return {
        "wire": len(execution["wire_currents"]),
        "patch": len(execution["patch_currents"]["rows"]),
        "volume": len(execution["volume_currents"] or []),
    }


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m9-release-") as name:
        directory = Path(name)
        bare, bare_report = run(HERE/"m9-bare-order.nec", directory,
                                "bare-builtin", "builtin")
        coated, coated_report = run(HERE/"m9-coated-order.nec", directory,
                                    "coated-builtin", "builtin")
        coated_lapack, _ = run(HERE/"m9-coated-order.nec", directory,
                               "coated-lapack", "lapack")

    if bare["schema"] != "nec2c-results" or coated["schema"] != bare["schema"]:
        raise AssertionError("coated job changed the result schema")
    if bare["schema_version"] != 1 or coated["schema_version"] != 1:
        raise AssertionError("coated job changed schema version 1")
    if "coated_wires" in bare["jobs"][0]:
        raise AssertionError("bare result acquired coated-only metadata")
    if bare["jobs"][0]["executions"][0]["material_evaluations"]:
        raise AssertionError("bare result acquired a material evaluation")
    metadata = coated["jobs"][0]["coated_wires"]
    if metadata["domains"][0]["segment_count"] != 401:
        raise AssertionError(metadata)
    evaluation = coated["jobs"][0]["executions"][0]["material_evaluations"][0]
    if evaluation["validity_status"] != "within-reference-envelope":
        raise AssertionError(evaluation)
    order = execution_order(bare)
    if order != {"wire": 401, "patch": 0, "volume": 0}:
        raise AssertionError(order)
    if execution_order(coated) != order or execution_order(coated_lapack) != order:
        raise AssertionError("coating or backend changed matrix order")
    for report in (bare_report, coated_report):
        if "TOTAL SEGMENTS USED: 401" not in report:
            raise AssertionError("traditional report lost segment-order witness")

    performance = json.loads((HERE/"m9-performance.json").read_text())
    expected_sizes = [95, 401, 801, 1601]
    if [case["segments"] for case in performance["cases"]] != expected_sizes:
        raise AssertionError("M9 performance size matrix changed")
    for case in performance["cases"]:
        if not case["matrix_order_equal"]:
            raise AssertionError(case)
        if case["bare"]["wire_unknowns"] != case["coated"]["wire_unknowns"]:
            raise AssertionError(case)
        if case["segments"] >= 401:
            overhead = case["overhead"]
            if not (3.0 <= overhead["matrix_fill_ratio"] <= 7.0):
                raise AssertionError(("unexpected coated fill cost", case))
            if overhead["wall_ratio"] > 4.0:
                raise AssertionError(("coated release overhead exceeded 4x", case))
            if (overhead["peak_rss_delta_kib"] > 2048 and
                    overhead["peak_rss_ratio"] > 1.08):
                raise AssertionError(("coating changed matrix memory materially", case))
            if abs(overhead["matrix_factor_ratio"]-1) > 0.10:
                raise AssertionError(("factor cost indicates changed order", case))

    acceptance = json.loads((HERE/"m9-acceptance.json").read_text())
    if acceptance["vswr"]["minimum"]["vswr_50_ohm"] > 1.5:
        raise AssertionError("manual VSWR article lost its match")
    if acceptance["vswr"]["validity_statuses"] != [
            "within-reference-envelope"]:
        raise AssertionError(acceptance["vswr"])
    if acceptance["pattern"]["symmetry_maximum_db_above_minus_40"] > 1e-6:
        raise AssertionError("manual pattern lost dipole symmetry")
    if acceptance["pattern"]["statuses"] != ["within-reference-envelope"]:
        raise AssertionError(acceptance["pattern"])

    evidence = (
        "m1-openems-caps.json", "m1-openems-audits.json",
        "m3-acceptance.json", "m4-acceptance.json", "m5-acceptance.json",
        "m6-openems.png", "m7-bend-openems.json",
        "m8-openems-thickness.json", "m8-geometry-analysis.json",
    )
    missing = [name for name in evidence if not (HERE/name).is_file()]
    if missing:
        raise AssertionError(("missing promoted evidence", missing))

    print("M9 bare schema/report compatibility and 401-unknown order passed")
    print("M9 builtin/LAPACK coated result shape passed")
    print("M9 retained 95/401/801/1601 performance and memory gates passed")
    print("M9 VSWR, pattern, and promoted-evidence bundle passed")


if __name__ == "__main__":
    main()
