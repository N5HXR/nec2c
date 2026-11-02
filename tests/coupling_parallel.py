#!/usr/bin/env python3
"""Require exact conductor-volume coupling across worker counts."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[1]))
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c"))
CASES = (
    (ROOT / "tests/performance/m5-coupled-216.nec", "fft"),
    (ROOT / "tests/performance/m7-connected-coupled-216.nec", "fft"),
    (ROOT / "tests/performance/m7-pec-coupled-27.nec", "fft"),
    (ROOT / "tests/performance/m7-real-coupled-27.nec", "dense"),
)


def scrub(value):
    if isinstance(value, dict):
        return {key: scrub(item) for key, item in value.items()
                if key != "elapsed_seconds"}
    if isinstance(value, list):
        return [scrub(item) for item in value]
    return value


def run(deck, solver, backend, threads, directory):
    stem = f"{deck.stem}-{threads}"
    result = directory / f"result-{stem}.json"
    timing = directory / f"timing-{stem}.json"
    command = [str(BINARY), "--jobs", "1", "--threads", str(threads),
        "-s", solver, "--fft-backend", backend,
        "--result-json", str(result), "--timing-json", str(timing),
        "-i", str(deck), "-o", str(directory / f"report-{stem}.out")]
    completed = subprocess.run(command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed\n{completed.stderr}")
    normalized = json.dumps(scrub(json.loads(result.read_text())),
        sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    phases = {item["name"]: item["seconds"]
              for item in json.loads(timing.read_text())["phases"]}
    if "dielectric_wire_to_volume_assembly" not in phases:
        raise AssertionError("missing conductor-volume coupling phase")
    report = (directory / f"report-{stem}.out").read_text()
    marker = "CONDUCTOR-VOLUME COUPLING WORKERS: "
    if marker not in report:
        raise AssertionError(f"missing coupling worker telemetry for {deck.name}")
    workers = int(report.split(marker, 1)[1].split()[0])
    return hashlib.sha256(normalized).hexdigest(), workers


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-coupling-") as temporary:
        directory = Path(temporary)
        for deck, solver in CASES:
            samples = {threads: run(deck, solver, "builtin", threads, directory)
                      for threads in (1, 2, 4)}
            hashes = {threads: sample[0] for threads, sample in samples.items()}
            if len(set(hashes.values())) != 1:
                raise AssertionError(f"{deck.name} results differ: {hashes}")
            if samples[4][1] <= 1:
                raise AssertionError(f"{deck.name} did not exercise parallel coupling: {samples}")
        auto = {threads: run(CASES[0][0], "fft", "auto", threads, directory)
                for threads in (1, 4)}
        if auto[1][0] != auto[4][0] or auto[4][1] <= 1:
            raise AssertionError(f"automatic FFT backend coupling differs: {auto}")
    print("free-space, connected patch, PEC, and real-ground coupling is "
          "bit-identical at 1, 2, and 4 workers; automatic FFT is exact")


if __name__ == "__main__":
    main()
