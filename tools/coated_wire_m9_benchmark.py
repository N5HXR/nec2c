#!/usr/bin/env python3
"""Measure coated-wire runtime and memory overhead without changing order."""

import argparse
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


def integer_list(value):
    try:
        result = [int(item) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not result or any(item < 3 or item % 2 == 0 for item in result):
        raise argparse.ArgumentTypeError("sizes must be odd integers of at least 3")
    return result


def deck(segments, coated):
    lines = [
        "CM M9 COATED-WIRE PERFORMANCE PAIR", "CE",
        f"GW 1 {segments} 0 0 -0.2375 0 0 0.2375 0.0002", "GE 0",
    ]
    if coated:
        lines.append("LD 7 1 0 0 3.2 0.0004 0.015 0")
    lines.extend(("FR 0 1 0 0 276.65 0",
                  f"EX 0 1 {(segments+1)//2} 0 1 0", "XQ 0", "EN", ""))
    return "\n".join(lines)


def phase(timing, name):
    row = next((item for item in timing["phases"] if item["name"] == name), None)
    return row["seconds"] if row else 0.0


def run_once(binary, source, coated, threads, directory, label):
    result_path = directory/f"{label}-result.json"
    timing_path = directory/f"{label}-timing.json"
    report_path = directory/f"{label}.out"
    environment = os.environ.copy()
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS"):
        environment[name] = str(threads)
    command = [str(binary), "--backend", "auto", "--jobs", "1",
               "--threads", str(threads), "--result-json", str(result_path),
               "--timing-json", str(timing_path), "-i", str(source),
               "-o", str(report_path)]
    completed = subprocess.run(command, cwd=ROOT, env=environment,
                               text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(json.dumps({"command": command,
            "returncode": completed.returncode, "stdout": completed.stdout,
            "stderr": completed.stderr}, indent=2))
    timing = json.loads(timing_path.read_text())
    result = json.loads(result_path.read_text())
    execution = result["jobs"][0]["executions"][0]
    evaluations = execution["material_evaluations"]
    if coated != bool(evaluations):
        raise AssertionError((coated, evaluations))
    domain_count = len(result["jobs"][0].get("coated_wires", {}).get(
        "domains", []))
    return {
        "matrix_fill_seconds": phase(timing, "matrix_fill"),
        "matrix_factor_seconds": phase(timing, "matrix_factor"),
        "matrix_solve_seconds": phase(timing, "matrix_solve"),
        "wall_seconds": timing["process"]["wall_seconds"],
        "peak_rss_kib": timing["process"]["peak_rss_kib"],
        "backend_used": timing["dense_backend"]["used"],
        "wire_unknowns": len(execution["wire_currents"]),
        "patch_unknown_rows": len(execution["patch_currents"]["rows"]),
        "volume_unknowns": (len(execution["volume_currents"] or [])),
        "coating_domains": domain_count,
    }


def median(rows, key):
    return statistics.median(row[key] for row in rows)


def summarize(rows):
    result = {"repeat": len(rows)}
    for key in ("matrix_fill_seconds", "matrix_factor_seconds",
                "matrix_solve_seconds", "wall_seconds", "peak_rss_kib"):
        result[key] = median(rows, key)
    for key in ("backend_used", "wire_unknowns", "patch_unknown_rows",
                "volume_unknowns", "coating_domains"):
        values = sorted({row[key] for row in rows})
        if len(values) != 1:
            raise AssertionError((key, values))
        result[key] = values[0]
    return result


def ratio(coated, bare, key):
    return coated[key]/bare[key] if bare[key] else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT/"nec2c")
    parser.add_argument("--sizes", type=integer_list,
                        default=integer_list("95,401,801"))
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.threads < 1 or arguments.warmup < 0 or arguments.repeat < 1:
        parser.error("threads/repeat must be positive and warmup nonnegative")

    binary = arguments.binary.resolve()
    cases = []
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m9-benchmark-") as name:
        directory = Path(name)
        for segments in arguments.sizes:
            sources = {}
            for coated in (False, True):
                source = directory/f"n{segments}-{'coated' if coated else 'bare'}.nec"
                source.write_text(deck(segments, coated), encoding="ascii")
                sources[coated] = source
            for index in range(arguments.warmup):
                order = (False, True) if index % 2 == 0 else (True, False)
                for coated in order:
                    run_once(binary, sources[coated], coated, arguments.threads,
                             directory, f"n{segments}-warmup-{index}-{coated}")
            runs = {False: [], True: []}
            for index in range(arguments.repeat):
                order = (False, True) if index % 2 == 0 else (True, False)
                for coated in order:
                    runs[coated].append(run_once(
                        binary, sources[coated], coated, arguments.threads,
                        directory, f"n{segments}-run-{index}-{coated}"))
            bare, coated = summarize(runs[False]), summarize(runs[True])
            if (bare["wire_unknowns"] != segments or
                    coated["wire_unknowns"] != segments or
                    bare["patch_unknown_rows"] != coated["patch_unknown_rows"] or
                    bare["volume_unknowns"] != coated["volume_unknowns"]):
                raise AssertionError("coating changed the matrix order")
            overhead = {
                "matrix_fill_ratio": ratio(coated, bare, "matrix_fill_seconds"),
                "matrix_factor_ratio": ratio(coated, bare, "matrix_factor_seconds"),
                "wall_ratio": ratio(coated, bare, "wall_seconds"),
                "peak_rss_ratio": ratio(coated, bare, "peak_rss_kib"),
                "peak_rss_delta_kib": coated["peak_rss_kib"]-bare["peak_rss_kib"],
            }
            cases.append({"segments": segments, "bare": bare, "coated": coated,
                          "matrix_order_equal": True, "overhead": overhead})
            print(f"{segments}: fill {overhead['matrix_fill_ratio']:.2f}x, "
                  f"wall {overhead['wall_ratio']:.2f}x, "
                  f"RSS {overhead['peak_rss_ratio']:.2f}x")

    report = {
        "schema": "nec2c-coated-wire-m9-performance",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": {"platform": platform.platform(), "logical_cpus": os.cpu_count()},
        "settings": {"threads": arguments.threads, "warmup": arguments.warmup,
                     "repeat": arguments.repeat, "sizes": arguments.sizes,
                     "backend": "auto"},
        "cases": cases,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2)+"\n")
    print(arguments.output)


if __name__ == "__main__":
    main()
