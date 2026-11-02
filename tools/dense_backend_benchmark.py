#!/usr/bin/env python3
"""Measure builtin/LAPACK crossover, accuracy, and LAPACK thread scaling."""

import argparse
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]


def parse_numbers(value):
    try:
        numbers = [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not numbers or any(number < 1 for number in numbers):
        raise argparse.ArgumentTypeError("values must be positive")
    return numbers


def write_deck(path, segments):
    if segments % 2 == 0:
        raise ValueError("segment counts must be odd for a center feed")
    path.write_text(
        "CM M2 SYNTHETIC DENSE BACKEND SCALING DIPOLE\n"
        "CE\n"
        f"GW 1 {segments} -5.05 0.0 0.0 5.05 0.0 0.0 0.0001\n"
        "GE 0\n"
        "FR 0 1 0 0 14.1\n"
        f"EX 0 1 {(segments + 1) // 2} 0 1.0 0.0\n"
        "XQ\n"
        "EN\n", encoding="ascii")


def phase(timing, name):
    return next(item["seconds"] for item in timing["phases"]
                if item["name"] == name)


def run_once(binary, deck, backend, threads, directory, label):
    result_path = directory / f"{label}-result.json"
    timing_path = directory / f"{label}-timing.json"
    output_path = directory / f"{label}.out"
    environment = os.environ.copy()
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS"):
        environment[name] = str(threads)
    command = [str(binary), "--backend", backend,
               "--result-json", str(result_path),
               "--timing-json", str(timing_path),
               "-i", str(deck), "-o", str(output_path)]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, env=environment,
                               text=True, capture_output=True)
    external = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(json.dumps({"command": command,
            "returncode": completed.returncode, "stdout": completed.stdout,
            "stderr": completed.stderr}, indent=2))
    timing = json.loads(timing_path.read_text(encoding="utf-8"))
    result = json.loads(result_path.read_text(encoding="utf-8"))
    feed = result["jobs"][0]["executions"][0]["feeds"][0]
    impedance = feed["impedance_ohm"]
    return {
        "factor_seconds": phase(timing, "matrix_factor"),
        "solve_seconds": phase(timing, "matrix_solve"),
        "wall_seconds": timing["process"]["wall_seconds"],
        "external_wall_seconds": external,
        "backend": timing["dense_backend"],
        "impedance_ohm": [impedance["re"], impedance["im"]],
    }


def summarize(runs, backend, threads):
    return {
        "backend_requested": backend,
        "backend_used": sorted({run["backend"]["used"] for run in runs}),
        "threads": threads,
        "repeat": len(runs),
        "factor_seconds": statistics.median(run["factor_seconds"] for run in runs),
        "solve_seconds": statistics.median(run["solve_seconds"] for run in runs),
        "wall_seconds": statistics.median(run["wall_seconds"] for run in runs),
        "impedance_ohm": runs[0]["impedance_ohm"],
    }


def measure(binary, deck, backend, threads, repeat, warmup, directory, prefix):
    for index in range(warmup):
        run_once(binary, deck, backend, threads, directory,
                 f"{prefix}-{backend}-warmup-{index}")
    runs = [run_once(binary, deck, backend, threads, directory,
                     f"{prefix}-{backend}-run-{index}")
            for index in range(repeat)]
    return summarize(runs, backend, threads)


def measure_pair(binary, deck, threads, repeat, warmup, directory, prefix):
    backends = ("builtin", "lapack")
    for index in range(warmup):
        for backend in backends[index % 2:] + backends[:index % 2]:
            run_once(binary, deck, backend, threads, directory,
                     f"{prefix}-{backend}-warmup-{index}")
    runs = {backend: [] for backend in backends}
    for index in range(repeat):
        order = backends[index % 2:] + backends[:index % 2]
        for backend in order:
            runs[backend].append(run_once(binary, deck, backend, threads, directory,
                                          f"{prefix}-{backend}-run-{index}"))
    return tuple(summarize(runs[backend], backend, threads) for backend in backends)


def relative_complex_error(left, right):
    a = complex(*left)
    b = complex(*right)
    return abs(a - b) / max(abs(a), 1.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", type=parse_numbers,
                        default=parse_numbers("31,63,95,127,129,191,255,401,801"))
    parser.add_argument("--threads", type=parse_numbers,
                        default=parse_numbers("1,2,4,8,16,32"))
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--warmup", type=int, default=2)
    arguments = parser.parse_args()
    if arguments.repeat < 1 or arguments.warmup < 0:
        parser.error("--repeat must be positive and --warmup cannot be negative")
    if any(size % 2 == 0 for size in arguments.sizes):
        parser.error("all --sizes values must be odd")

    binary = arguments.binary.resolve()
    report = {
        "schema": "nec2c-dense-backend-benchmark",
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": {"platform": platform.platform(), "logical_cpus": os.cpu_count()},
        "settings": {"sizes": arguments.sizes, "threads": arguments.threads,
                     "repeat": arguments.repeat, "warmup": arguments.warmup},
        "crossover": [], "thread_scaling": [],
    }
    with tempfile.TemporaryDirectory(prefix="nec2c-dense-backend-") as name:
        directory = Path(name)
        for size in arguments.sizes:
            deck = directory / f"dipole-{size}.nec"
            write_deck(deck, size)
            builtin, lapack = measure_pair(binary, deck, 1, arguments.repeat,
                                           arguments.warmup, directory, f"n{size}")
            automatic = measure(binary, deck, "auto", 1, 1, 0, directory,
                                f"n{size}-auto")
            report["crossover"].append({
                "order": size, "builtin": builtin, "lapack": lapack,
                "automatic_backend": automatic["backend_used"][0],
                "factor_speedup": builtin["factor_seconds"] / lapack["factor_seconds"],
                "wall_speedup": builtin["wall_seconds"] / lapack["wall_seconds"],
                "relative_impedance_error": relative_complex_error(
                    builtin["impedance_ohm"], lapack["impedance_ohm"]),
            })
        largest = max(arguments.sizes)
        deck = directory / f"dipole-{largest}.nec"
        write_deck(deck, largest)
        for threads in arguments.threads:
            report["thread_scaling"].append(measure(
                binary, deck, "lapack", threads, arguments.repeat,
                arguments.warmup, directory, f"scale-n{largest}-t{threads}"))

    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(output)
    for row in report["crossover"]:
        print(f"n={row['order']:4d} auto={row['automatic_backend']:7s} "
              f"factor={row['factor_speedup']:.3f}x "
              f"error={row['relative_impedance_error']:.3e}")
    print(output)


if __name__ == "__main__":
    main()
