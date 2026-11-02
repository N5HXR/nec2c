#!/usr/bin/env python3
"""Benchmark and validate deterministic OpenMP dielectric assembly scaling."""

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
        "id": "volume-216",
        "deck": "tests/performance/m5-volume-216.nec",
        "resource": "tests/performance/m5-volume-216.diel",
        "arguments": ["-s", "dense", "--backend", "builtin"],
        "phases": ["dielectric_volume_assembly"],
        "description": "216-cell volume-only direct interaction assembly",
    },
    {
        "id": "coupled-rows-216",
        "deck": "tests/performance/m5-coupled-216.nec",
        "resource": "tests/performance/m5-coupled-216.diel",
        "arguments": ["-s", "fft", "--fft-backend", "builtin"],
        "phases": [
            "dielectric_volume_to_wire_assembly",
            "dielectric_volume_to_patch_assembly",
        ],
        "description": "216-cell mixed wire/patch destination-row assembly",
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


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summary(values):
    if not values:
        return None
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "maximum": max(values),
    }


def execute(binary, case, threads, index, warmup, directory, environment, timeout):
    label = "warmup" if warmup else "run"
    stem = f"{case['id']}-t{threads}-{label}{index}"
    report_path = directory / (stem + ".out")
    result_path = directory / (stem + "-result.json")
    timing_path = directory / (stem + "-timing.json")
    command = [
        str(binary),
        "--jobs", "1",
        "--threads", str(threads),
        *case["arguments"],
        "--result-json", str(result_path),
        "--timing-json", str(timing_path),
        "-i", str((ROOT / case["deck"]).resolve()),
        "-o", str(report_path),
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
        raise RuntimeError(
            json.dumps(
                {
                    "command": command,
                    "returncode": completed.returncode,
                    "stdout": completed.stdout,
                    "stderr": completed.stderr,
                },
                indent=2,
            )
        )
    timing = json.loads(timing_path.read_text(encoding="utf-8"))
    phases = {item["name"]: item["seconds"] for item in timing["phases"]}
    selected = {name: phases.get(name) for name in case["phases"]}
    available = [value for value in selected.values() if value is not None]
    return {
        "index": index,
        "command": command,
        "external_wall_seconds": external_wall,
        "process_wall_seconds": timing["process"]["wall_seconds"],
        "process_cpu_seconds": timing["process"]["cpu_seconds"],
        "peak_rss_kib": timing["process"]["peak_rss_kib"],
        "phases": selected,
        "assembly_seconds": sum(available) if len(available) == len(selected) else None,
        "numerical_sha256": numerical_hash(result_path),
    }


def summarize_runs(runs):
    phase_names = sorted({name for run in runs for name in run["phases"]})
    return {
        "assembly_seconds": summary(
            [run["assembly_seconds"] for run in runs if run["assembly_seconds"] is not None]
        ),
        "phases": {
            name: summary(
                [run["phases"][name] for run in runs if run["phases"][name] is not None]
            )
            for name in phase_names
        },
        "external_wall_seconds": summary([run["external_wall_seconds"] for run in runs]),
        "process_wall_seconds": summary([run["process_wall_seconds"] for run in runs]),
        "process_cpu_seconds": summary([run["process_cpu_seconds"] for run in runs]),
        "peak_rss_kib": summary(
            [run["peak_rss_kib"] for run in runs if run["peak_rss_kib"] is not None]
        ),
        "stable_within_thread_count": len({run["numerical_sha256"] for run in runs}) == 1,
    }


def parse_threads(value):
    try:
        result = [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("requires comma-separated positive integers") from error
    if not result or any(item < 1 for item in result) or len(set(result)) != len(result):
        raise argparse.ArgumentTypeError("requires unique comma-separated positive integers")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=parse_threads, default=parse_threads("1,2,4,8,16,32"))
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument(
        "--minimum-speedup",
        type=float,
        default=2.0,
        help="required best assembly speedup for each case; use 0 for a smoke run",
    )
    parser.add_argument("--case", action="append", choices=[case["id"] for case in CASES])
    parser.add_argument(
        "--allow-missing-phases",
        action="store_true",
        help="permit a pre-M5 binary to smoke-test the runner without assembly timers",
    )
    arguments = parser.parse_args()
    if (arguments.repeat < 1 or arguments.warmup < 0 or arguments.timeout <= 0 or
            not math.isfinite(arguments.minimum_speedup) or arguments.minimum_speedup < 0):
        parser.error("--repeat and --timeout must be positive; --warmup cannot be negative")
    if 1 not in arguments.threads:
        parser.error("--threads must include 1 as the speedup baseline")

    binary = arguments.binary.resolve()
    selected_cases = [case for case in CASES if not arguments.case or case["id"] in arguments.case]
    environment = os.environ.copy()
    environment["OMP_DYNAMIC"] = "FALSE"
    report = {
        "schema": "nec2c-m5-assembly-benchmark",
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": {
            "platform": platform.platform(),
            "logical_cpus": os.cpu_count(),
            "cpu_affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        },
        "binary": {"path": str(binary), "sha256": file_hash(binary)},
        "settings": {
            "threads": arguments.threads,
            "repeat": arguments.repeat,
            "warmup": arguments.warmup,
            "jobs": 1,
            "minimum_speedup": arguments.minimum_speedup,
            "omp_dynamic": environment["OMP_DYNAMIC"],
            "omp_proc_bind": environment.get("OMP_PROC_BIND"),
            "omp_places": environment.get("OMP_PLACES"),
        },
        "cases": [],
    }
    missing_phases = []
    numerical_mismatches = []
    with tempfile.TemporaryDirectory(prefix="nec2c-m5-assembly-") as temporary:
        directory = Path(temporary)
        for case in selected_cases:
            case_report = {
                "id": case["id"],
                "description": case["description"],
                "deck": case["deck"],
                "deck_sha256": file_hash(ROOT / case["deck"]),
                "resource": case["resource"],
                "resource_sha256": file_hash(ROOT / case["resource"]),
                "arguments": case["arguments"],
                "assembly_phases": case["phases"],
                "thread_counts": [],
            }
            reference_hash = None
            for threads in arguments.threads:
                for index in range(arguments.warmup):
                    execute(binary, case, threads, index, True, directory, environment, arguments.timeout)
                runs = [
                    execute(binary, case, threads, index, False, directory, environment, arguments.timeout)
                    for index in range(arguments.repeat)
                ]
                hashes = {run["numerical_sha256"] for run in runs}
                if reference_hash is None:
                    reference_hash = next(iter(hashes))
                if hashes != {reference_hash}:
                    numerical_mismatches.append(f"{case['id']}@{threads}")
                if any(run["assembly_seconds"] is None for run in runs):
                    missing_phases.append(f"{case['id']}@{threads}")
                case_report["thread_counts"].append({
                    "threads": threads,
                    "runs": runs,
                    "summary": summarize_runs(runs),
                })
            baseline = next(
                item["summary"]["assembly_seconds"]
                for item in case_report["thread_counts"]
                if item["threads"] == 1
            )
            baseline_median = baseline["median"] if baseline else None
            for item in case_report["thread_counts"]:
                measured = item["summary"]["assembly_seconds"]
                item["assembly_speedup"] = (
                    baseline_median / measured["median"]
                    if baseline_median is not None and measured and measured["median"] > 0
                    else None
                )
                item["assembly_efficiency"] = (
                    item["assembly_speedup"] / item["threads"]
                    if item["assembly_speedup"] is not None else None
                )
            measured_speedups = [
                item["assembly_speedup"] for item in case_report["thread_counts"]
                if item["assembly_speedup"] is not None
            ]
            case_report["best_assembly_speedup"] = (
                max(measured_speedups) if measured_speedups else None
            )
            case_report["meets_minimum_speedup"] = (
                case_report["best_assembly_speedup"] is not None and
                case_report["best_assembly_speedup"] >= arguments.minimum_speedup
            )
            case_report["numerically_identical_across_threads"] = not any(
                mismatch.startswith(case["id"] + "@") for mismatch in numerical_mismatches
            )
            report["cases"].append(case_report)

    speedup_failures = [
        case["id"] for case in report["cases"]
        if (case["best_assembly_speedup"] is not None and
            not case["meets_minimum_speedup"])
    ]
    report["validation"] = {
        "missing_phase_measurements": missing_phases,
        "numerical_mismatches": numerical_mismatches,
        "speedup_failures": speedup_failures,
        "passed": (not numerical_mismatches and not speedup_failures and
                   (not missing_phases or arguments.allow_missing_phases)),
        "missing_phases_allowed": arguments.allow_missing_phases,
    }
    atomic_json(arguments.output.resolve(), report)
    for case in report["cases"]:
        values = []
        for item in case["thread_counts"]:
            timing = item["summary"]["assembly_seconds"]
            values.append(
                f"{item['threads']}t=" +
                (f"{timing['median']:.6f}s/{item['assembly_speedup']:.2f}x" if timing else "missing")
            )
        print(f"{case['id']}: " + ", ".join(values))
    print(arguments.output.resolve())
    if not report["validation"]["passed"]:
        raise SystemExit("M5 assembly benchmark validation failed")


if __name__ == "__main__":
    main()
