#!/usr/bin/env python3
"""Check deterministic dielectric matrix assembly across thread counts."""

import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[1])).resolve()
PROGRAM = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
CASES = (
    ("volume", ROOT / "tests/dielectric/m5-cube.nec", ("dielectric_volume_assembly",)),
    ("coupled", ROOT / "tests/dielectric/m7-mixed-cube.nec", (
        "dielectric_volume_assembly",
        "dielectric_volume_to_wire_assembly",
        "dielectric_volume_to_patch_assembly",
    )),
)


def normalize(document):
    for job in document["jobs"]:
        for execution in job["executions"]:
            execution["solver"]["elapsed_seconds"] = 0
    return document


def run(directory, label, deck, threads):
    result = directory / f"{label}-{threads}.json"
    timing = directory / f"{label}-{threads}-timing.json"
    report = directory / f"{label}-{threads}.out"
    subprocess.run(
        [str(PROGRAM), "-s", "dense", "--backend", "builtin",
         "--jobs", "1", "--threads", str(threads),
         "--result-json", str(result), "--timing-json", str(timing),
         "-i", str(deck), "-o", str(report)],
        cwd=ROOT, check=True, timeout=60,
    )
    return normalize(json.loads(result.read_text(encoding="utf-8"))), json.loads(
        timing.read_text(encoding="utf-8"))


with tempfile.TemporaryDirectory(prefix="nec2c-m5-assembly-") as temporary:
    directory = Path(temporary)
    for label, deck, required_phases in CASES:
        reference, serial_timing = run(directory, label, deck, 1)
        assert serial_timing["concurrency"]["threads_per_job"] == 1
        for threads in (2, 4):
            result, timing = run(directory, label, deck, threads)
            assert result == reference
            assert timing["concurrency"]["threads_per_job"] == threads
            phases = {phase["name"] for phase in timing["phases"]}
            assert set(required_phases) <= phases
        automatic, timing = run(directory, label, deck, "auto")
        assert automatic == reference
        concurrency = timing["concurrency"]
        assert concurrency["threads_requested"] == "auto"
        assert 1 <= concurrency["threads_per_job"] <= 16
        assert concurrency["threads_per_job"] <= concurrency["allowed_cpus"]

print("dielectric assembly is bit-identical at 1, 2, and 4 threads")
