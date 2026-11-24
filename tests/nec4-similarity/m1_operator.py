#!/usr/bin/env python3
"""M1 engine-boundary, oracle, and frozen-classic acceptance checks."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def run(command, expected=0):
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode != expected:
        raise AssertionError((command, completed.returncode,
                              completed.stdout, completed.stderr))
    return completed


capabilities = json.loads(run([str(BINARY), "--capabilities-json"]).stdout)
features = capabilities["features"]
assert features["wire_engines"][:2] == ["auto", "nec2"]
assert features["ground_engines"][:2] == ["auto", "nec2"]
assert "nec4ish" in features["ground_engines"]
assert features["planned_ground_engines"] == []
assert "nec4ish" in features["wire_engines"]

commands = [
        ([str(BINARY), "--ground-engine", "adaptive", "--capabilities-json"],
         "planned but not available"),
        ([str(BINARY), "--model-check", "maybe", "--capabilities-json"],
         "warn, strict, or off")]
for command, fragment in commands:
    completed = run(command, expected=1)
    assert fragment in completed.stderr

with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m1-") as directory:
    output = Path(directory)
    run([sys.executable, str(ROOT / "tools" / "nec4_m1_oracles.py"),
         "--self-test", "--output-dir", str(output)])
    assert (output / "m1-oracles.json").is_file()
    assert (output / "m1-kernel-error-map.png").is_file()
    assert (output / "m1-sommerfeld-error-map.png").is_file()

run([sys.executable,
     str(ROOT / "tests" / "nec4-similarity" / "m0_baseline.py"),
     "--binary", str(BINARY)])
print("M1 operator boundary, independent oracles, and classic identity passed")
