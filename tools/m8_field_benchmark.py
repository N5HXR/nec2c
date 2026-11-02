#!/usr/bin/env python3
"""Benchmark M8 far-field, power-quadrature, and complete-solve scaling."""

import argparse
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "classic-field-801": {
        "deck": ROOT / "tests/performance/m8-field-801.nec",
        "arguments": [], "phase": "far_field_compute", "overall": "far_field"},
    "coupled-power": {
        "deck": ROOT / "tests/performance/m5-coupled-216.nec",
        "arguments": ["-s", "fft", "--fft-backend", "auto"],
        "phase": "dielectric_power_compute", "overall": "dielectric_solve"},
}


def parse_threads(value):
    try:
        result = [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("requires comma-separated integers") from error
    if not result or any(item < 1 for item in result) or len(set(result)) != len(result):
        raise argparse.ArgumentTypeError("requires unique positive integers")
    return result


def scrub(value):
    if isinstance(value, dict):
        return {key: scrub(item) for key, item in value.items()
                if key != "elapsed_seconds"}
    if isinstance(value, list):
        return [scrub(item) for item in value]
    return value


def relative_delta(reference, candidate, path="result"):
    if type(reference) is not type(candidate):
        raise AssertionError(f"type mismatch at {path}")
    if isinstance(reference, dict):
        if reference.keys() != candidate.keys():
            raise AssertionError(f"key mismatch at {path}")
        return max((relative_delta(reference[key], candidate[key],
                                   f"{path}/{key}") for key in reference), default=0.0)
    if isinstance(reference, list):
        if len(reference) != len(candidate):
            raise AssertionError(f"length mismatch at {path}")
        return max((relative_delta(a, b, f"{path}/{index}")
                    for index, (a, b) in enumerate(zip(reference, candidate))),
                   default=0.0)
    if isinstance(reference, float):
        if not math.isfinite(reference) or not math.isfinite(candidate):
            return 0.0 if reference == candidate else math.inf
        return abs(reference-candidate)/max(1.0, abs(reference), abs(candidate))
    if reference != candidate:
        raise AssertionError(f"value mismatch at {path}: {reference!r} != {candidate!r}")
    return 0.0


def execute(case, threads, index, warmup, directory, timeout):
    label = "warmup" if warmup else "run"
    stem = f"{case}-{threads}-{label}-{index}"
    result = directory / f"{stem}-result.json"
    timing = directory / f"{stem}-timing.json"
    report = directory / f"{stem}.out"
    definition = CASES[case]
    command = [str(ROOT / "nec2c"), "--backend", "auto", "--jobs", "1",
        "--threads", str(threads), *definition["arguments"],
        "--result-json", str(result), "--timing-json", str(timing),
        "-i", str(definition["deck"]), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed\n{completed.stderr}")
    timing_document = json.loads(timing.read_text())
    phases = {item["name"]: item["seconds"] for item in timing_document["phases"]}
    return {"index": index, "phase_seconds": phases[definition["phase"]],
        "overall_seconds": phases[definition["overall"]],
        "wall_seconds": timing_document["process"]["wall_seconds"],
        "peak_rss_kib": timing_document["process"]["peak_rss_kib"],
        "result": scrub(json.loads(result.read_text()))}


def median(runs, name):
    return statistics.median(run[name] for run in runs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=parse_threads,
        default=parse_threads("1,2,4,8,16,32"))
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--acceptance-threads", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if (1 not in arguments.threads or arguments.acceptance_threads not in arguments.threads or
            arguments.warmup < 0 or arguments.repeat < 1 or arguments.timeout <= 0):
        parser.error("thread list must include 1 and the acceptance count")

    output = {"schema": "nec2c-m8-field-benchmark", "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "settings": vars(arguments) | {"output": str(arguments.output)}, "cases": []}
    reference_results = {}
    with tempfile.TemporaryDirectory(prefix="nec2c-m8-fields-") as temporary:
        directory = Path(temporary)
        for case in CASES:
            item = {"id": case, "deck": str(CASES[case]["deck"].relative_to(ROOT)),
                    "phase": CASES[case]["phase"], "thread_counts": []}
            for threads in arguments.threads:
                for index in range(arguments.warmup):
                    execute(case, threads, index, True, directory, arguments.timeout)
                runs = [execute(case, threads, index, False, directory, arguments.timeout)
                        for index in range(arguments.repeat)]
                reference = reference_results.setdefault(case, runs[0]["result"])
                delta = max(relative_delta(reference, run["result"]) for run in runs)
                for run in runs:
                    del run["result"]
                item["thread_counts"].append({"threads": threads, "runs": runs,
                    "maximum_relative_numerical_delta": delta,
                    "phase_median_seconds": median(runs, "phase_seconds"),
                    "overall_median_seconds": median(runs, "overall_seconds"),
                    "wall_median_seconds": median(runs, "wall_seconds")})
            baseline = item["thread_counts"][0]
            for count in item["thread_counts"]:
                count["phase_speedup"] = baseline["phase_median_seconds"] / count["phase_median_seconds"]
                count["overall_speedup"] = baseline["overall_median_seconds"] / count["overall_median_seconds"]
            output["cases"].append(item)

    failures = []
    for case in output["cases"]:
        acceptance = next(item for item in case["thread_counts"]
                          if item["threads"] == arguments.acceptance_threads)
        if any(item["maximum_relative_numerical_delta"] > 1e-10
               for item in case["thread_counts"]):
            failures.append(f"{case['id']}: numerical delta")
        minimum_phase = 1.5 if case["id"] == "classic-field-801" else 1.2
        if acceptance["phase_speedup"] < minimum_phase:
            failures.append(f"{case['id']}: phase speedup")
        if case["id"] == "coupled-power" and acceptance["overall_speedup"] < 1.2:
            failures.append(f"{case['id']}: solve speedup")
    output["validation"] = {"failures": failures, "passed": not failures}
    target = arguments.output.resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
    for case in output["cases"]:
        print(case["id"] + ": " + ", ".join(
            f"{item['threads']}t={item['phase_median_seconds']:.6f}s/"
            f"{item['phase_speedup']:.2f}x phase, {item['overall_speedup']:.2f}x overall"
            for item in case["thread_counts"]))
    print(target)
    if failures:
        raise SystemExit("M8 field acceptance failed: " + ", ".join(failures))


if __name__ == "__main__":
    main()
