#!/usr/bin/env python3
"""Benchmark classic NEC matrix-fill scaling and deterministic results."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
CASES = (
    {
        "id": "wire-801",
        "deck": "tests/performance/scaling-dipole-801.nec",
        "description": "801-segment connected free-space wire fill",
        "provenance": "Established M1 deterministic scaling fixture",
    },
    {
        "id": "mixed-ground-651",
        "deck": "tests/performance/m6-mixed-ground-651.nec",
        "description": "201 extended-kernel wire segments and 225 patches over PEC ground",
        "provenance": "M6 deterministic mixed-interaction scaling fixture",
    },
)


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scrub_timing(value):
    if isinstance(value, dict):
        return {
            key: scrub_timing(item)
            for key, item in value.items()
            if key != "elapsed_seconds"
        }
    if isinstance(value, list):
        return [scrub_timing(item) for item in value]
    return value


def numerical_hash(path):
    value = scrub_timing(json.loads(path.read_text(encoding="utf-8")))
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def summary(values):
    if not values:
        return None
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "maximum": max(values),
    }


def parse_threads(value):
    try:
        result = [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "requires comma-separated positive integers"
        ) from error
    if not result or any(item < 1 for item in result) or len(set(result)) != len(result):
        raise argparse.ArgumentTypeError(
            "requires unique comma-separated positive integers"
        )
    return result


def execute(binary, case, threads, index, warmup, directory, environment, timeout):
    kind = "warmup" if warmup else "run"
    stem = f"{case['id']}-t{threads}-{kind}{index}"
    output = directory / f"{stem}.out"
    result = directory / f"{stem}-result.json"
    timing = directory / f"{stem}-timing.json"
    command = [
        str(binary),
        "--backend", "builtin",
        "--jobs", "1",
        "--threads", str(threads),
        "--result-json", str(result),
        "--timing-json", str(timing),
        "-i", str((ROOT / case["deck"]).resolve()),
        "-o", str(output),
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    external_wall = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(json.dumps({
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }, indent=2))
    timing_value = json.loads(timing.read_text(encoding="utf-8"))
    phases = {phase["name"]: phase["seconds"] for phase in timing_value["phases"]}
    return {
        "index": index,
        "command": command,
        "matrix_fill_seconds": phases.get("matrix_fill"),
        "external_wall_seconds": external_wall,
        "process_wall_seconds": timing_value["process"]["wall_seconds"],
        "process_cpu_seconds": timing_value["process"]["cpu_seconds"],
        "peak_rss_kib": timing_value["process"]["peak_rss_kib"],
        "numerical_sha256": numerical_hash(result),
    }


def summarize_runs(runs):
    return {
        "matrix_fill_seconds": summary([
            run["matrix_fill_seconds"] for run in runs
            if run["matrix_fill_seconds"] is not None
        ]),
        "external_wall_seconds": summary([run["external_wall_seconds"] for run in runs]),
        "process_wall_seconds": summary([run["process_wall_seconds"] for run in runs]),
        "process_cpu_seconds": summary([run["process_cpu_seconds"] for run in runs]),
        "peak_rss_kib": summary([
            run["peak_rss_kib"] for run in runs if run["peak_rss_kib"] is not None
        ]),
        "stable_within_thread_count": len({
            run["numerical_sha256"] for run in runs
        }) == 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--threads", type=parse_threads, default=parse_threads("1,2,4,8,16,32")
    )
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--minimum-speedup", type=float, default=2.0)
    parser.add_argument("--case", action="append", choices=[case["id"] for case in CASES])
    arguments = parser.parse_args()
    if (arguments.repeat < 1 or arguments.warmup < 0 or arguments.timeout <= 0 or
            not math.isfinite(arguments.minimum_speedup) or
            arguments.minimum_speedup < 0):
        parser.error("invalid repeat, warmup, timeout, or minimum speedup")
    if 1 not in arguments.threads:
        parser.error("--threads must include 1 as the speedup baseline")

    binary = arguments.binary.resolve()
    selected = [case for case in CASES if not arguments.case or case["id"] in arguments.case]
    environment = os.environ.copy()
    environment["OMP_DYNAMIC"] = "FALSE"
    report = {
        "schema": "nec2c-m6-fill-benchmark",
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": {
            "platform": platform.platform(),
            "logical_cpus": os.cpu_count(),
            "cpu_affinity": sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity") else None,
        },
        "binary": {"path": str(binary), "sha256": sha256(binary)},
        "settings": {
            "backend": "builtin",
            "jobs": 1,
            "threads": arguments.threads,
            "repeat": arguments.repeat,
            "warmup": arguments.warmup,
            "minimum_speedup": arguments.minimum_speedup,
            "omp_dynamic": environment["OMP_DYNAMIC"],
            "omp_proc_bind": environment.get("OMP_PROC_BIND"),
            "omp_places": environment.get("OMP_PLACES"),
        },
        "cases": [],
    }
    numerical_mismatches = []
    missing_fill = []
    with tempfile.TemporaryDirectory(prefix="nec2c-m6-fill-") as temporary:
        directory = Path(temporary)
        for case in selected:
            case_report = {
                **case,
                "deck_sha256": sha256(ROOT / case["deck"]),
                "thread_counts": [],
            }
            reference_hash = None
            for threads in arguments.threads:
                for index in range(arguments.warmup):
                    execute(
                        binary, case, threads, index, True, directory,
                        environment, arguments.timeout
                    )
                runs = [
                    execute(
                        binary, case, threads, index, False, directory,
                        environment, arguments.timeout
                    )
                    for index in range(arguments.repeat)
                ]
                hashes = {run["numerical_sha256"] for run in runs}
                if reference_hash is None:
                    reference_hash = next(iter(hashes))
                if hashes != {reference_hash}:
                    numerical_mismatches.append(f"{case['id']}@{threads}")
                if any(run["matrix_fill_seconds"] is None for run in runs):
                    missing_fill.append(f"{case['id']}@{threads}")
                case_report["thread_counts"].append({
                    "threads": threads,
                    "runs": runs,
                    "summary": summarize_runs(runs),
                })

            baseline = next(
                item["summary"]["matrix_fill_seconds"]
                for item in case_report["thread_counts"] if item["threads"] == 1
            )
            baseline_median = baseline["median"] if baseline else None
            for item in case_report["thread_counts"]:
                measured = item["summary"]["matrix_fill_seconds"]
                speedup = (
                    baseline_median / measured["median"]
                    if baseline_median is not None and measured and measured["median"] > 0
                    else None
                )
                item["matrix_fill_speedup"] = speedup
                item["matrix_fill_efficiency"] = (
                    speedup / item["threads"] if speedup is not None else None
                )
            speedups = [
                item["matrix_fill_speedup"] for item in case_report["thread_counts"]
                if item["matrix_fill_speedup"] is not None
            ]
            case_report["best_matrix_fill_speedup"] = max(speedups) if speedups else None
            case_report["meets_minimum_speedup"] = (
                case_report["best_matrix_fill_speedup"] is not None and
                case_report["best_matrix_fill_speedup"] >= arguments.minimum_speedup
            )
            case_report["numerically_identical_across_threads"] = not any(
                mismatch.startswith(case["id"] + "@")
                for mismatch in numerical_mismatches
            )
            report["cases"].append(case_report)

    speedup_failures = [
        case["id"] for case in report["cases"] if not case["meets_minimum_speedup"]
    ]
    report["validation"] = {
        "missing_matrix_fill": missing_fill,
        "numerical_mismatches": numerical_mismatches,
        "speedup_failures": speedup_failures,
        "passed": not missing_fill and not numerical_mismatches and not speedup_failures,
    }
    atomic_json(arguments.output.resolve(), report)
    for case in report["cases"]:
        values = []
        for item in case["thread_counts"]:
            timing = item["summary"]["matrix_fill_seconds"]
            values.append(
                f"{item['threads']}t=" +
                (f"{timing['median']:.6f}s/{item['matrix_fill_speedup']:.2f}x"
                 if timing else "missing")
            )
        print(f"{case['id']}: " + ", ".join(values))
    print(arguments.output.resolve())
    if not report["validation"]["passed"]:
        raise SystemExit("M6 matrix-fill acceptance failed")


if __name__ == "__main__":
    main()
