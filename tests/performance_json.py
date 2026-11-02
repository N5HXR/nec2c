#!/usr/bin/env python3
"""Contract tests for the opt-in performance timing artifact."""

import json
import os
from pathlib import Path
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
    legacy_schema = load(ROOT / "docs/nec2c-timing-v1.schema.json")
    assert legacy_schema["properties"]["schema_version"] == {"const": 1}
    assert "fft_backend" not in legacy_schema["properties"]
    schema = load(ROOT / "docs/nec2c-timing-v2.schema.json")
    assert schema["properties"]["schema"] == {"const": "nec2c-timing"}
    assert schema["properties"]["schema_version"] == {"const": 2}

    capabilities = json.loads(run("--capabilities-json").stdout)
    assert capabilities["features"]["timing_json_versions"] == [2]
    assert capabilities["features"]["default_dense_backend"] == "auto"
    assert capabilities["features"]["default_jobs"] == "auto"
    assert capabilities["features"]["default_threads"] == "auto"
    assert capabilities["features"]["dense_backends"][:2] == ["auto", "builtin"]
    assert capabilities["features"]["fft_backends"][:2] == ["auto", "builtin"]
    assert capabilities["features"]["default_fft_backend"] == "auto"
    assert capabilities["limits"]["lapack_auto_min_order"] >= 1

    deck = ROOT / "Input/EX1.nec"
    timing = temporary / "timing.json"
    run("--timing-json", timing, "-i", deck, "-o", temporary / "classic.out")
    artifact = load(timing)
    assert artifact["schema"] == "nec2c-timing"
    assert artifact["schema_version"] == 2
    assert artifact["completion"] == {"status": "ok"}
    assert Path(artifact["input"]["deck"]) == deck
    assert artifact["clock"] == {"wall": "monotonic", "cpu": "process"}
    assert artifact["dense_backend"] == {"requested": "auto", "used": "builtin"}
    assert artifact["fft_backend"] == {"requested": None, "used": None}
    assert artifact["process"]["wall_seconds"] >= 0.0
    assert artifact["process"]["cpu_seconds"] >= 0.0
    assert artifact["process"]["peak_rss_kib"] is None or artifact["process"]["peak_rss_kib"] >= 0
    phases = {phase["name"]: phase for phase in artifact["phases"]}
    for name in ("matrix_fill", "matrix_factor"):
        assert phases[name]["seconds"] >= 0.0
        assert phases[name]["calls"] >= 1
    assert not list(temporary.glob("timing.json.tmp.*"))

    builtin_timing = temporary / "builtin.json"
    run("--backend", "builtin", "--timing-json", builtin_timing,
        "-i", deck, "-o", temporary / "builtin.out")
    assert load(builtin_timing)["dense_backend"] == {
        "requested": "builtin", "used": "builtin"}
    assert "DENSE BACKEND: builtin (REQUESTED builtin)" in \
        (temporary / "builtin.out").read_text()

    if "lapack" in capabilities["features"]["dense_backends"]:
        lapack_timing = temporary / "lapack.json"
        run("--backend", "lapack", "--timing-json", lapack_timing,
            "-i", deck, "-o", temporary / "lapack.out")
        assert load(lapack_timing)["dense_backend"] == {
            "requested": "lapack", "used": "lapack"}

    run("--backend", "not-a-backend", "-i", deck,
        "-o", temporary / "invalid-backend.out", expect=1)

    run("--timing-json", "-", "-i", deck, "-o", temporary / "dash.out", expect=1)
    shared = temporary / "shared.json"
    run("--result-json", shared, "--timing-json", shared,
        "-i", deck, "-o", temporary / "shared.out", expect=1)
    run("--capabilities-json", "--timing-json", temporary / "unused.json", expect=1)
    assert not (temporary / "unused.json").exists()

    missing = temporary / "missing.nec"
    failed_timing = temporary / "failed.json"
    run("--timing-json", failed_timing, "-i", missing,
        "-o", temporary / "missing.out", expect=255)
    assert load(failed_timing)["completion"] == {"status": "error"}

print("performance timing contract tests passed")
