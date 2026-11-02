#!/usr/bin/env python3
"""Require deterministic ordered far fields from serial and parallel plans."""

import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[1]))
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c"))
CLASSIC = ROOT / "tests/performance/m8-field-101.nec"
DIELECTRIC = ROOT / "tests/performance/m8-volume-field-216.nec"
GROUND_CACHE = ROOT / "tests/performance/m8-ground-cache.nec"
GROUND_CACHE_MISS = ROOT / "tests/performance/m8-ground-cache-miss.nec"


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


def run(deck, threads, directory, solver=None, backend=None):
    stem = f"{deck.stem}-{threads}"
    result = directory / f"{stem}.json"
    timing = directory / f"{stem}-timing.json"
    report = directory / f"{stem}.out"
    command = [str(BINARY), "--jobs", "1", "--threads", str(threads),
        "--result-json", str(result), "--timing-json", str(timing)]
    if solver is not None:
        command += ["-s", solver]
    if backend is not None:
        command += ["--backend", backend]
    command += ["-i", str(deck), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed\n{completed.stderr}")
    document = scrub(json.loads(result.read_text()))
    timing_document = json.loads(timing.read_text())
    return document, scrub_report(report.read_text()), timing_document


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-m8-fields-") as temporary:
        directory = Path(temporary)
        classic_serial = run(CLASSIC, 1, directory)
        classic_auto = run(CLASSIC, "auto", directory)
        if classic_serial[:2] != classic_auto[:2]:
            raise AssertionError("classic automatic far fields differ from serial")
        if classic_auto[2]["concurrency"]["threads_per_job"] <= 1:
            raise AssertionError("classic automatic plan did not receive field workers")

        volume_serial = run(DIELECTRIC, 1, directory, "dense", "builtin")
        volume_parallel = run(DIELECTRIC, 4, directory, "dense", "builtin")
        if volume_serial[:2] != volume_parallel[:2]:
            raise AssertionError("volume-only parallel far fields differ from serial")
        if volume_parallel[2]["concurrency"]["threads_per_job"] != 4:
            raise AssertionError("volume-only plan did not receive four threads")

        cached = run(GROUND_CACHE, 1, directory)
        phases = {item["name"]: item for item in cached[2]["phases"]}
        if phases.get("ground_table_cache_miss", {}).get("calls") != 1:
            raise AssertionError("repeated real-ground setup did not record one cache miss")
        if phases.get("ground_table_cache_hit", {}).get("calls") != 1:
            raise AssertionError("repeated real-ground setup did not record one cache hit")

        invalidated = run(GROUND_CACHE_MISS, 1, directory)
        phases = {item["name"]: item for item in invalidated[2]["phases"]}
        if phases.get("ground_table_cache_miss", {}).get("calls") != 2:
            raise AssertionError("frequency change did not invalidate ground cache")
        if "ground_table_cache_hit" in phases:
            raise AssertionError("frequency change incorrectly reused ground cache")
    print("classic and volume-only far fields preserve exact ordered output under parallel evaluation; exact Sommerfeld setup reuse passed")


if __name__ == "__main__":
    main()
