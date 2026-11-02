#!/usr/bin/env python3
"""Benchmark deterministic conductor-to-volume coupling scaling."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
DECK = ROOT / "tests/performance/m5-coupled-216.nec"
PHASE = "dielectric_wire_to_volume_assembly"


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


def digest_json(path):
    encoded = json.dumps(scrub(json.loads(path.read_text())),
        sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def digest_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summary(values):
    return {"minimum": min(values), "median": statistics.median(values),
            "maximum": max(values)}


def execute(binary, backend, threads, index, warmup, directory, environment, timeout):
    label = "warmup" if warmup else "run"
    stem = f"coupled-216-t{threads}-{label}{index}"
    result = directory / f"{stem}-result.json"
    timing = directory / f"{stem}-timing.json"
    report_path = directory / f"{stem}.out"
    command = [str(binary), "--jobs", "1", "--threads", str(threads),
        "-s", "fft", "--fft-backend", backend,
        "--result-json", str(result), "--timing-json", str(timing),
        "-i", str(DECK), "-o", str(report_path)]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, env=environment, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout, check=False)
    external = time.perf_counter()-started
    if completed.returncode != 0:
        raise RuntimeError(json.dumps({"command": command,
            "returncode": completed.returncode, "stdout": completed.stdout,
            "stderr": completed.stderr}, indent=2))
    document = json.loads(timing.read_text())
    phases = {item["name"]: item["seconds"] for item in document["phases"]}
    match = re.search(r"CONDUCTOR-VOLUME COUPLING WORKERS: (\d+)([^\n]*)",
                      report_path.read_text())
    if not match:
        raise RuntimeError("missing conductor-volume worker telemetry")
    return {"index": index, "command": command,
        "coupling_seconds": phases.get(PHASE),
        "solve_seconds": phases.get("dielectric_solve"),
        "external_wall_seconds": external,
        "process_wall_seconds": document["process"]["wall_seconds"],
        "process_cpu_seconds": document["process"]["cpu_seconds"],
        "peak_rss_kib": document["process"]["peak_rss_kib"],
        "fft_backend_used": document["fft_backend"]["used"],
        "coupling_workers_used": int(match.group(1)),
        "coupling_serial_fallback": "FALLBACK" in match.group(2),
        "numerical_sha256": digest_json(result)}


def run_summary(runs):
    return {name: summary([run[name] for run in runs if run[name] is not None])
        for name in ("coupling_seconds", "solve_seconds", "external_wall_seconds",
                     "process_wall_seconds", "process_cpu_seconds", "peak_rss_kib")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=parse_threads,
        default=parse_threads("1,2,4,8,16,32"))
    parser.add_argument("--fft-backend", choices=("auto", "builtin", "fftw"),
        default="auto")
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--minimum-coupling-speedup", type=float, default=2.0)
    parser.add_argument("--minimum-solve-speedup", type=float, default=1.2)
    parser.add_argument("--acceptance-threads", type=int, default=16)
    arguments = parser.parse_args()
    if (1 not in arguments.threads or arguments.warmup < 0 or arguments.repeat < 1 or
            arguments.timeout <= 0 or
            not math.isfinite(arguments.minimum_coupling_speedup) or
            not math.isfinite(arguments.minimum_solve_speedup) or
            arguments.minimum_coupling_speedup < 0 or
            arguments.minimum_solve_speedup < 0 or
            arguments.acceptance_threads not in arguments.threads):
        parser.error("invalid thread, repetition, timeout, or speedup setting")

    binary = arguments.binary.resolve()
    environment = os.environ.copy()
    environment["OMP_DYNAMIC"] = "FALSE"
    report = {"schema": "nec2c-m7-coupling-benchmark", "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": {"platform": platform.platform(), "logical_cpus": os.cpu_count(),
            "cpu_affinity": sorted(os.sched_getaffinity(0))
                if hasattr(os, "sched_getaffinity") else None},
        "binary": {"path": str(binary), "sha256": digest_file(binary)},
        "fixture": {"deck": str(DECK.relative_to(ROOT)),
            "sha256": digest_file(DECK)},
        "settings": {"threads": arguments.threads, "jobs": 1,
            "fft_backend": arguments.fft_backend,
            "warmup": arguments.warmup, "repeat": arguments.repeat,
            "minimum_coupling_speedup": arguments.minimum_coupling_speedup,
            "minimum_solve_speedup": arguments.minimum_solve_speedup,
            "acceptance_threads": arguments.acceptance_threads,
            "omp_dynamic": environment["OMP_DYNAMIC"],
            "omp_proc_bind": environment.get("OMP_PROC_BIND"),
            "omp_places": environment.get("OMP_PLACES")}, "thread_counts": []}
    reference = None
    mismatches = []
    with tempfile.TemporaryDirectory(prefix="nec2c-m7-coupling-") as temporary:
        directory = Path(temporary)
        for threads in arguments.threads:
            for index in range(arguments.warmup):
                execute(binary, arguments.fft_backend, threads, index, True,
                        directory, environment, arguments.timeout)
            runs = [execute(binary, arguments.fft_backend, threads, index, False,
                            directory, environment, arguments.timeout)
                    for index in range(arguments.repeat)]
            hashes = {run["numerical_sha256"] for run in runs}
            if reference is None:
                reference = next(iter(hashes))
            if hashes != {reference}:
                mismatches.append(threads)
            report["thread_counts"].append({"threads": threads, "runs": runs,
                "coupling_workers_used": sorted({run["coupling_workers_used"]
                                                  for run in runs}),
                "fft_backends_used": sorted({run["fft_backend_used"]
                                              for run in runs}),
                "serial_fallback": any(run["coupling_serial_fallback"]
                                       for run in runs),
                "summary": run_summary(runs)})

    baseline = next(item for item in report["thread_counts"] if item["threads"] == 1)
    coupling_base = baseline["summary"]["coupling_seconds"]["median"]
    solve_base = baseline["summary"]["solve_seconds"]["median"]
    for item in report["thread_counts"]:
        coupling = item["summary"]["coupling_seconds"]["median"]
        solve = item["summary"]["solve_seconds"]["median"]
        item["coupling_speedup"] = coupling_base/coupling
        item["solve_speedup"] = solve_base/solve
    best_coupling = max(item["coupling_speedup"] for item in report["thread_counts"])
    best_solve = max(item["solve_speedup"] for item in report["thread_counts"])
    acceptance = next(item for item in report["thread_counts"]
                      if item["threads"] == arguments.acceptance_threads)
    telemetry_failures = [item["threads"] for item in report["thread_counts"]
        if item["threads"] > 1 and
           (item["coupling_workers_used"] == [1] or item["serial_fallback"])]
    report["validation"] = {"numerical_mismatches": mismatches,
        "best_coupling_speedup": best_coupling, "best_solve_speedup": best_solve,
        "acceptance_threads": arguments.acceptance_threads,
        "acceptance_coupling_speedup": acceptance["coupling_speedup"],
        "acceptance_solve_speedup": acceptance["solve_speedup"],
        "telemetry_failures": telemetry_failures,
        "passed": not mismatches and not telemetry_failures and
            acceptance["coupling_speedup"] >= arguments.minimum_coupling_speedup and
            acceptance["solve_speedup"] >= arguments.minimum_solve_speedup}
    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name+".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    temporary.replace(output)
    for item in report["thread_counts"]:
        coupling=item["summary"]["coupling_seconds"]["median"]
        solve=item["summary"]["solve_seconds"]["median"]
        print(f"{item['threads']:2d} workers: coupling {coupling:.6f}s/"
              f"{item['coupling_speedup']:.2f}x, solve {solve:.6f}s/"
              f"{item['solve_speedup']:.2f}x")
    print(output)
    if not report["validation"]["passed"]:
        raise SystemExit("M7 coupling acceptance failed")


if __name__ == "__main__":
    main()
