#!/usr/bin/env python3
"""Verify the immutable pre-coated-wire numerical and text baseline."""

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
HERE = ROOT / "tests" / "coated-wire"
DECK = HERE / "m0-bare-dipole.nec"
GOLDEN = HERE / "m0-baseline.json"


def canonical_result(path):
    document = json.loads(path.read_text(encoding="utf-8"))
    # The deck path is invocation-dependent; its content hash remains covered.
    document["input"]["deck"] = "tests/coated-wire/m0-bare-dipole.nec"
    for job in document["jobs"]:
        for execution in job["executions"]:
            solver = execution.get("solver")
            if solver is not None:
                solver["elapsed_seconds"] = None
    return json.dumps(document, sort_keys=True, separators=(",", ":"))


def normalized_report(text):
    text = re.sub(r"TOTAL RUN TIME: \d+", "TOTAL RUN TIME: X", text)
    text = re.sub(r"FILL: \d+ msec  FACTOR: \d+ msec",
                  "FILL: X msec  FACTOR: X msec", text)
    return text


def digest(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


with tempfile.TemporaryDirectory(prefix="nec2c-coated-m0-") as directory:
    temporary = Path(directory)
    result = temporary / "result.json"
    report = temporary / "report.out"
    command = [str(BINARY), "--backend", "builtin", "--jobs", "1",
               "--threads", "1", "--result-json", str(result),
               "-i", str(DECK), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((completed.returncode, completed.stdout,
                              completed.stderr))

    observed = {
        "deck_sha256": hashlib.sha256(DECK.read_bytes()).hexdigest(),
        "canonical_result_sha256": digest(canonical_result(result)),
        "normalized_report_sha256": digest(
            normalized_report(report.read_text(encoding="utf-8"))),
    }
    expected = json.loads(GOLDEN.read_text(encoding="utf-8"))
    for key, value in observed.items():
        if expected[key] != value:
            raise AssertionError(f"{key}: expected {expected[key]}, got {value}")

print("coated-wire M0 bare-wire baseline is unchanged")
