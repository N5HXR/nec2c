#!/usr/bin/env python3
"""Measure disabled and enabled timing overhead against an archived binary."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]


def normalized_hash(path):
    text = path.read_text(encoding="ascii", errors="replace")
    text = re.sub(r"FILL: \d+ msec  FACTOR: \d+ msec",
                  "FILL: TIME msec  FACTOR: TIME msec", text)
    text = re.sub(r"TOTAL RUN TIME: \d+ msec", "TOTAL RUN TIME: TIME msec", text)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def execute(binary, deck, directory, label, index, timing):
    output = directory / f"{label}-{index}.out"
    command = [str(binary), "-i", str(deck), "-o", str(output)]
    if timing:
        command[1:1] = ["--timing-json", str(directory / f"{label}-{index}.json")]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError((command, completed.returncode, completed.stdout, completed.stderr))
    return {"seconds": elapsed, "numerical_sha256": normalized_hash(output)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--deck", type=Path,
                        default=ROOT / "tests/performance/scaling-dipole-801.nec")
    parser.add_argument("--repeat", type=int, default=25)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.repeat < 3 or arguments.warmup < 0:
        parser.error("--repeat must be at least 3 and --warmup cannot be negative")
    binaries = {
        "baseline": arguments.baseline.resolve(),
        "disabled": arguments.candidate.resolve(),
        "enabled": arguments.candidate.resolve(),
    }
    samples = {label: [] for label in binaries}
    labels = list(binaries)
    with tempfile.TemporaryDirectory(prefix="nec2c-timing-overhead-") as temporary:
        directory = Path(temporary)
        for index in range(arguments.warmup):
            for label in labels:
                execute(binaries[label], arguments.deck.resolve(), directory,
                        f"warmup-{label}", index, label == "enabled")
        for index in range(arguments.repeat):
            order = labels[index % len(labels):] + labels[:index % len(labels)]
            for label in order:
                samples[label].append(execute(
                    binaries[label], arguments.deck.resolve(), directory,
                    label, index, label == "enabled"))
    medians = {label: statistics.median(item["seconds"] for item in values)
               for label, values in samples.items()}
    baseline = medians["baseline"]
    overhead = {label: 100.0 * (medians[label] / baseline - 1.0)
                for label in ("disabled", "enabled")}
    hashes = {item["numerical_sha256"] for values in samples.values() for item in values}
    report = {
        "schema": "nec2c-timing-overhead",
        "schema_version": 1,
        "deck": str(arguments.deck.resolve()),
        "repeat": arguments.repeat,
        "warmup": arguments.warmup,
        "median_seconds": medians,
        "overhead_percent": overhead,
        "limits_percent": {"disabled": 2.0, "enabled": 5.0},
        "stable_numerical_output": len(hashes) == 1,
        "passed": (overhead["disabled"] <= 2.0 and overhead["enabled"] <= 5.0 and
                   len(hashes) == 1),
        "samples": samples,
    }
    arguments.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                encoding="utf-8")
    print(f"disabled overhead: {overhead['disabled']:+.3f}%")
    print(f"enabled overhead: {overhead['enabled']:+.3f}%")
    print("numerical output: " + ("stable" if len(hashes) == 1 else "CHANGED"))
    print("overhead contract: " + ("PASS" if report["passed"] else "FAIL"))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
