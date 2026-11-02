#!/usr/bin/env python3
"""Interleave built-in and FFTW runs to compare backend latency fairly."""

import argparse
import json
import math
import statistics
from pathlib import Path

from fft_benchmark import ROOT, parse_counts, sample


def medians(samples):
    result = {
        name: statistics.median(item[name] for item in samples)
        for name in ("build_seconds", "matvec_seconds", "solve_seconds",
                     "peak_rss_kib")
    }
    result["build_plus_solve_seconds"] = statistics.median(
        item["build_seconds"] + item["solve_seconds"] for item in samples
    )
    return result


def equivalent(left, right):
    for index in (0, 1):
        if not math.isclose(float.fromhex(left[index]), float.fromhex(right[index]),
                            rel_tol=2e-12, abs_tol=2e-12):
            return False
    return math.isclose(float(left[3]), float(right[3]),
                        rel_tol=2e-12, abs_tol=2e-12)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "tools/m9_fft_benchmark")
    parser.add_argument("--model", choices=("tree", "skeleton"), default="tree")
    parser.add_argument("--threads", type=parse_counts, default=parse_counts("1,16"))
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--max-iterations", type=int, default=300)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.warmup < 0 or arguments.repeat < 1 or arguments.max_iterations < 1:
        parser.error("warmup must be nonnegative; repeat and max-iterations must be positive")

    evidence = {
        "schema": "nec2c-fft-backend-comparison-v1",
        "model": arguments.model,
        "warmup_pairs": arguments.warmup,
        "measured_pairs": arguments.repeat,
        "maximum_iterations": arguments.max_iterations,
        "results": [],
    }
    for threads in arguments.threads:
        for iteration in range(arguments.warmup):
            order = ("builtin", "fftw") if iteration % 2 == 0 else ("fftw", "builtin")
            for backend in order:
                sample(arguments.binary, arguments.model, threads,
                       arguments.max_iterations, backend)
        samples = {"builtin": [], "fftw": []}
        for iteration in range(arguments.repeat):
            order = ("builtin", "fftw") if iteration % 2 == 0 else ("fftw", "builtin")
            for backend in order:
                samples[backend].append(sample(
                    arguments.binary, arguments.model, threads,
                    arguments.max_iterations, backend
                ))
        signatures = {}
        for backend in samples:
            observed = {tuple(item["signature"]) for item in samples[backend]}
            if len(observed) != 1:
                raise RuntimeError(f"{backend} result is not repeatable at {threads} threads")
            signatures[backend] = next(iter(observed))
        if not equivalent(signatures["builtin"], signatures["fftw"]):
            raise RuntimeError(f"backend numerical mismatch at {threads} threads")
        summary = {backend: medians(samples[backend]) for backend in samples}
        ratio = {
            name: summary["builtin"][name] / summary["fftw"][name]
            for name in ("build_seconds", "matvec_seconds", "solve_seconds",
                         "build_plus_solve_seconds")
        }
        evidence["results"].append({
            "threads": threads,
            "median": summary,
            "builtin_over_fftw": ratio,
            "samples": samples,
        })
        print(
            f"{threads:2d} threads: FFTW matvec {ratio['matvec_seconds']:.2f}x, "
            f"solve {ratio['solve_seconds']:.2f}x, "
            f"build+solve {ratio['build_plus_solve_seconds']:.2f}x"
        )
    if arguments.output:
        arguments.output.write_text(json.dumps(evidence, indent=2) + "\n")
        print(arguments.output)


if __name__ == "__main__":
    main()
