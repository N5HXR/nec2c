#!/usr/bin/env python3
"""Exercise prepared network state across a multi-angle excitation sweep."""

import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[1]))
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c"))
DECK = ROOT / "tests/performance/m8-network-multi-angle.nec"


def scrub(value):
    if isinstance(value, dict):
        return {key: scrub(item) for key, item in value.items()
                if key != "elapsed_seconds"}
    if isinstance(value, list):
        return [scrub(item) for item in value]
    return value


def scrub_report(text):
    text = re.sub(r"FILL: \d+ msec  FACTOR: \d+ msec",
                  "FILL: <time>  FACTOR: <time>", text)
    return re.sub(r"TOTAL RUN TIME: \d+ msec", "TOTAL RUN TIME: <time>", text)


def run(threads, directory):
    stem = f"network-{threads}"
    result = directory / f"{stem}.json"
    timing = directory / f"{stem}-timing.json"
    report = directory / f"{stem}.out"
    command = [str(BINARY), "--backend", "builtin", "--jobs", "1",
        "--threads", str(threads), "--result-json", str(result),
        "--timing-json", str(timing), "-i", str(DECK), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed\n{completed.stderr}")
    return scrub(json.loads(result.read_text())), scrub_report(report.read_text()), \
        json.loads(timing.read_text())


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-m8-network-") as temporary:
        directory = Path(temporary)
        serial = run(1, directory)
        parallel = run(4, directory)
        if serial[:2] != parallel[:2]:
            raise AssertionError("prepared network sweep differs by thread count")
        report = serial[1]
        if report.count("STRUCTURE EXCITATION DATA AT NETWORK CONNECTION POINTS") != 6:
            raise AssertionError("multi-angle network sweep did not emit six solutions")
        phases = {item["name"]: item for item in serial[2]["phases"]}
        if phases.get("network_cache_miss", {}).get("calls") != 1:
            raise AssertionError("network preparation should miss exactly once")
        if phases.get("network_cache_hit", {}).get("calls") != 5:
            raise AssertionError("network preparation should hit five times")
    print("prepared network state is exact across six ordered plane-wave excitations")


if __name__ == "__main__":
    main()
