#!/usr/bin/env python3
"""Measure built-in FFT scaling with a fixed one-iteration workload."""

import argparse
import json
import os
import platform
import re
import statistics
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_counts(value):
    counts = []
    for item in value.split(","):
        count = int(item)
        if count < 1:
            raise argparse.ArgumentTypeError("thread counts must be positive")
        counts.append(count)
    return counts


def fields(output):
    return dict(re.findall(r"([a-z0-9_]+)=([^\s]+)", output))


def sample(binary, model, threads, maximum_iterations, backend):
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(threads)
    completed = subprocess.run(
        [str(binary), model, "1e-8", str(maximum_iterations), backend],
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=60,
        check=False,
    )
    expected_returncode = 1 if maximum_iterations == 1 else 0
    if completed.returncode != expected_returncode:
        raise RuntimeError(
            f"benchmark returned {completed.returncode}: {completed.stderr}"
        )
    parsed = fields(completed.stdout)
    required = (
        "build_seconds",
        "matvec_seconds",
        "solve_seconds",
        "peak_rss_kib",
        "matvec_checksum_real",
        "matvec_checksum_imag",
        "matvec_bitwise_hash",
        "residual",
        "threads",
        "backend",
    )
    missing = [name for name in required if name not in parsed]
    if missing:
        raise RuntimeError(f"missing {missing}: {completed.stdout}")
    captured_threads = int(parsed["threads"])
    if threads > 1 and captured_threads == 1:
        raise RuntimeError(
            "the FFT benchmark was built without OpenMP; run "
            "`autoreconf -fi && ./configure --with-lapack=auto && make -j8`, "
            "and verify configure reports `parallel built-in FFT: enabled`"
        )
    if captured_threads != threads:
        raise RuntimeError(
            f"requested {threads} threads, plan captured {parsed['threads']}"
        )
    if parsed["backend"] != backend and backend != "auto":
        raise RuntimeError(
            f"requested {backend} backend, plan captured {parsed['backend']}"
        )
    return {
        "build_seconds": float(parsed["build_seconds"]),
        "matvec_seconds": float(parsed["matvec_seconds"]),
        "solve_seconds": float(parsed["solve_seconds"]),
        "peak_rss_kib": int(parsed["peak_rss_kib"]),
        "signature": [
            parsed["matvec_checksum_real"],
            parsed["matvec_checksum_imag"],
            parsed["matvec_bitwise_hash"],
            parsed["residual"],
        ],
        "backend": parsed["backend"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "tools/m9_fft_benchmark")
    parser.add_argument("--model", choices=("tree", "skeleton"), default="tree")
    parser.add_argument("--threads", type=parse_counts, default=parse_counts("1,2,4,8,16,32"))
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--max-iterations", type=int, default=1)
    parser.add_argument("--backend", choices=("auto", "builtin", "fftw"), default="auto")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.warmup < 0 or arguments.repeat < 1 or arguments.max_iterations < 1:
        parser.error("warmup must be nonnegative; repeat and max-iterations must be positive")
    evidence = {
        "schema": "nec2c-fft-benchmark-v1",
        "model": arguments.model,
        "system": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "affinity_cpus": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        },
        "warmup": arguments.warmup,
        "repeat": arguments.repeat,
        "maximum_iterations": arguments.max_iterations,
        "backend_requested": arguments.backend,
        "results": [],
    }
    reference = None
    serial_median = None
    for count in arguments.threads:
        for _ in range(arguments.warmup):
            sample(arguments.binary, arguments.model, count, arguments.max_iterations, arguments.backend)
        samples = [
            sample(arguments.binary, arguments.model, count, arguments.max_iterations, arguments.backend)
            for _ in range(arguments.repeat)
        ]
        signatures = {tuple(item["signature"]) for item in samples}
        if len(signatures) != 1:
            raise RuntimeError(f"{count}-thread numerical result is not repeatable")
        signature = next(iter(signatures))
        if reference is None:
            reference = signature
        elif signature != reference:
            raise RuntimeError(f"{count}-thread numerical result differs from serial")
        medians = {
            name: statistics.median(item[name] for item in samples)
            for name in ("build_seconds", "matvec_seconds", "solve_seconds", "peak_rss_kib")
        }
        if serial_median is None:
            serial_median = medians["matvec_seconds"]
        result = {
            "threads": count,
            "median": medians,
            "matvec_speedup": serial_median / medians["matvec_seconds"],
            "samples": samples,
        }
        evidence["results"].append(result)
        print(
            f"{count:2d} threads: matvec {medians['matvec_seconds']:.6f}s, "
            f"speedup {result['matvec_speedup']:.2f}x"
        )
    if arguments.output:
        arguments.output.write_text(json.dumps(evidence, indent=2) + "\n")
        print(arguments.output)


if __name__ == "__main__":
    main()
