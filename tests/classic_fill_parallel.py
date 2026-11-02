#!/usr/bin/env python3
"""Require exact classic matrix-fill results across process counts."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[1]))
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c"))
DECKS = (
    "tests/performance/scaling-dipole-801.nec",
    "tests/performance/m6-mixed-ground-651.nec",
    "tests/performance/m6-connected-patch-657.nec",
    "tests/performance/m6-loaded-network-402.nec",
    "tests/performance/m6-real-ground-201.nec",
    "tests/performance/m6-symmetry-604.nec",
    "tests/performance/m6-sweep-301.nec",
)


def without_timing(value):
    if isinstance(value, dict):
        return {
            key: without_timing(item)
            for key, item in value.items()
            if key != "elapsed_seconds"
        }
    if isinstance(value, list):
        return [without_timing(item) for item in value]
    return value


def run(deck, threads, directory):
    stem = f"{Path(deck).stem}-t{threads}"
    result = directory / f"{stem}.json"
    command = [
        str(BINARY), "--backend", "builtin", "--jobs", "1",
        "--threads", str(threads), "--result-json", str(result),
        "-i", str(ROOT / deck), "-o", str(directory / f"{stem}.out"),
    ]
    completed = subprocess.run(command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{' '.join(command)} failed ({completed.returncode})\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    normalized = json.dumps(without_timing(json.loads(result.read_text())),
        sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(normalized).hexdigest()


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-classic-fill-") as temporary:
        directory = Path(temporary)
        for deck in DECKS:
            serial = run(deck, 1, directory)
            parallel = run(deck, 4, directory)
            if serial != parallel:
                raise AssertionError(f"serial/parallel mismatch: {deck}")
    print("classic matrix fill is exact across interaction families and repeated builtin fills")


if __name__ == "__main__":
    main()
