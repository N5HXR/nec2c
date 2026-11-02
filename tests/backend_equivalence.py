#!/usr/bin/env python3
"""Compare builtin and LAPACK result artifacts across representative models."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


BINARY = Path(os.environ.get("NEC2C", "./nec2c")).resolve()
ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
RELATIVE_TOLERANCE = 1.0e-8
ABSOLUTE_TOLERANCE = 1.0e-10


def invoke(*arguments, expect=0):
    completed = subprocess.run([str(BINARY), *map(str, arguments)], cwd=ROOT,
                               text=True, capture_output=True)
    if completed.returncode != expect:
        raise AssertionError((completed.returncode, completed.stdout, completed.stderr))
    return completed


def compare(left, right, path="$"):
    if isinstance(left, bool) or left is None or isinstance(left, str):
        assert left == right, (path, left, right)
    elif isinstance(left, (int, float)):
        assert isinstance(right, (int, float)), (path, left, right)
        assert math.isclose(left, right, rel_tol=RELATIVE_TOLERANCE,
                            abs_tol=ABSOLUTE_TOLERANCE), (path, left, right)
    elif isinstance(left, list):
        assert isinstance(right, list) and len(left) == len(right), path
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            compare(left_item, right_item, f"{path}[{index}]")
    else:
        assert isinstance(right, dict) and left.keys() == right.keys(), path
        for key in left:
            if path.endswith(".solver") and key in {"requested", "used", "elapsed_seconds"}:
                continue
            compare(left[key], right[key], f"{path}.{key}")


capabilities = json.loads(invoke("--capabilities-json").stdout)
if "lapack" not in capabilities["features"]["dense_backends"]:
    invoke("--backend", "lapack", "-i", ROOT / "Input/EX1.nec",
           "-o", os.devnull, expect=1)
    print("LAPACKE unavailable; explicit-backend rejection passed")
    raise SystemExit(0)

decks = [
    ROOT / "Input/EX1.nec",
    ROOT / "Input/EX4.nec",
    ROOT / "Input/EX5.nec",
    ROOT / "Input/EX7.nec",
    ROOT / "tests/decks/patch-only.nec",
    ROOT / "tests/decks/nt-network.nec",
    ROOT / "tests/performance/loop-yagi-438.nec",
    ROOT / "tests/performance/scaling-dipole-801.nec",
    ROOT / "tests/dielectric/m6-cube.nec",
]

with tempfile.TemporaryDirectory(prefix="nec2c-backend-equivalence-") as name:
    directory = Path(name)
    for index, deck in enumerate(decks):
        values = {}
        for backend in ("builtin", "lapack"):
            result = directory / f"{index}-{backend}.json"
            invoke("--backend", backend, "--result-json", result,
                   "-i", deck, "-o", directory / f"{index}-{backend}.out")
            values[backend] = json.loads(result.read_text(encoding="utf-8"))
        compare(values["builtin"], values["lapack"])

print(f"builtin/LAPACK equivalence passed for {len(decks)} models")
