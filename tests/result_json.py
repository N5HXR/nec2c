#!/usr/bin/env python3
"""Contract tests for the opt-in Antenna Turtle JSON interface."""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


BINARY = Path(os.environ.get("NEC2C", "./nec2c")).resolve()
ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()


def run(*args, expect=0):
    completed = subprocess.run([str(BINARY), *map(str, args)], cwd=ROOT,
                               text=True, capture_output=True)
    if completed.returncode != expect:
        raise AssertionError((completed.returncode, completed.stdout, completed.stderr))
    return completed


def load(path):
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


with tempfile.TemporaryDirectory() as directory:
    temporary = Path(directory)
    schema = load(ROOT / "docs/nec2c-results-v1.schema.json")
    assert schema["properties"]["schema"]["const"] == "nec2c-results"

    capabilities = json.loads(run("--capabilities-json").stdout)
    assert capabilities["schema"] == "nec2c-capabilities"
    assert capabilities["schema_version"] == 1
    assert capabilities["features"]["dielectric_resource_versions"] == [1, 2, 3]
    assert capabilities["features"]["result_json_versions"] == [1]
    assert capabilities["features"]["grounds"] == ["free-space", "pec", "sommerfeld-norton"]

    deck = ROOT / "tests/decks/patch-only.nec"
    plain = temporary / "plain.out"
    observed = temporary / "observed.out"
    classic_json = temporary / "classic.json"
    run("-i", deck, "-o", plain)
    run("--result-json", classic_json, "-i", deck, "-o", observed)
    normalize = lambda value: re.sub(r"TOTAL RUN TIME: \d+", "TOTAL RUN TIME: X", value)
    assert normalize(plain.read_text()) == normalize(observed.read_text())
    classic = load(classic_json)
    assert classic["completion"] == {"status": "ok"}
    job = classic["jobs"][0]
    assert job["dielectric"] is None and job["mesh_states"] == []
    execution = job["executions"][0]
    assert execution["mesh_state_id"] is None
    assert execution["solver"]["used"] == "legacy-nec-direct"
    for key in ("relative_tolerance", "maximum_iterations", "iterations",
                "reported_residual", "true_residual", "converged"):
        assert execution["solver"][key] is None
    assert execution["volume_currents"] is None
    assert len(execution["patch_currents"]["rows"]) == 4

    validation_json = temporary / "classic-validation.json"
    run("--validate-only", "--result-json", validation_json,
        "-i", deck, "-o", temporary / "classic-validation.out")
    validation = load(validation_json)
    assert validation["jobs"][0]["validation"] == "classic-input-only"
    assert validation["jobs"][0]["executions"] == []

    dielectric_json = temporary / "dielectric.json"
    run("--result-json", dielectric_json,
        "-i", ROOT / "tests/dielectric/m5-cube.nec",
        "-o", temporary / "dielectric.out")
    dielectric = load(dielectric_json)
    job = dielectric["jobs"][0]
    assert len(job["mesh_states"]) == 1 and len(job["executions"]) == 1
    execution = job["executions"][0]
    assert execution["status"] == "ok" and execution["mesh_state_id"] == 0
    assert execution["solver"]["used"] == "dense"
    assert execution["power"]["basis"] == "plane-wave-scattered-balance"
    assert len(execution["volume_currents"]["rows"]) == 1
    assert execution["volume_currents"]["rows"][0][0:4] == [1, 1, 1, 1]
    assert len(execution["far_fields"][0]["samples"]) == 6

    shutil.copy(ROOT / "tests/dielectric/m5-cube.diel", temporary / "m5-cube.diel")
    repeated = temporary / "repeated.nec"
    repeated.write_text("""CM RESULT STATE ID TEST
CE
GE
FR 0 1 0 0 100.0
EX 1 1 1 0 90.0 0.0 90.0 0.0 0.0 0.0
DF \"m5-cube.diel\"
RP 0 1 1 0000 0.0 0.0 0.0 0.0 0.0 0.0
FR 0 1 0 0 100.0
GN -1
RP 0 1 1 0000 0.0 0.0 0.0 0.0 0.0 0.0
EN
""", encoding="ascii")
    repeated_json = temporary / "repeated.json"
    run("--result-json", repeated_json, "-i", repeated,
        "-o", temporary / "repeated.out")
    repeated_result = load(repeated_json)["jobs"][0]
    assert [state["mesh_state_id"] for state in repeated_result["mesh_states"]] == [0, 1]
    assert [item["mesh_state_id"] for item in repeated_result["executions"]] == [0, 1]
    assert [item["trigger"]["request_index"] for item in repeated_result["executions"]] == [0, 1]
    assert [item["ground_state_id"] for item in repeated_result["executions"]] == [0, 1]

    invalid = temporary / "invalid.nec"
    invalid.write_text("CM INVALID RESOURCE\nCE\nGE\nDF \"missing.diel\"\nEN\n", encoding="ascii")
    invalid_json = temporary / "invalid.json"
    run("--result-json", invalid_json, "-i", invalid,
        "-o", temporary / "invalid.out", expect=255)
    failure = load(invalid_json)
    assert failure["completion"] == {"status": "validation-error"}
    assert failure["diagnostics"][0]["severity"] == "error"

print("result JSON contract tests passed")
