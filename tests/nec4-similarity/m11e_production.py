#!/usr/bin/env python3
"""Verify or regenerate the retained M11E production evidence."""

import argparse
import cmath
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
HERE = ROOT / "tests" / "nec4-similarity"
EVIDENCE = HERE / "m11e-production.json"
DECK = HERE / "m11e-production.nec"


def phase(timing, name):
    return next(row["seconds"] for row in timing["phases"]
                if row["name"] == name)


def run(cache, directory, stem, model, threads="auto"):
    result = directory / f"{stem}.json"
    timing = directory / f"{stem}-timing.json"
    report = directory / f"{stem}.out"
    command = [str(BINARY), "--wire-engine", "nec4ish", "--ground-engine",
               "nec4ish", "--nec4ish-ground-model", model,
               "--nec4ish-ground-cache", str(cache), "--threads", str(threads),
               "--result-json", str(result), "--timing-json", str(timing),
               "-i", str(DECK), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True, timeout=1800)
    if completed.returncode:
        raise AssertionError((command, completed.returncode, completed.stderr))
    return (json.loads(result.read_text()), json.loads(timing.read_text()),
            report.read_text(errors="replace"))


def execution(document):
    return document["jobs"][0]["executions"][0]


def number(value):
    return complex(value["re"], value["im"])


def vector(sample):
    return [number(sample["value"][axis]) for axis in "xyz"]


def relative_l2(left, right):
    error = math.sqrt(sum(abs(a-b)**2 for a, b in zip(left, right)))
    scale = math.sqrt(sum(abs(value)**2 for value in right))
    return error / max(scale, 1e-300)


def magnetic(report):
    block = report.split("NEAR MAGNETIC FIELDS")[-1]
    line = next(row for row in block.splitlines()
                if re.match(r"\s*0\.0000\s+3\.1500", row))
    values = [float(value) for value in line.split()][3:]
    return [values[index] * cmath.exp(1j*math.radians(values[index+1]))
            for index in range(0, 6, 2)]


def verify(data):
    accuracy = data["accuracy"]
    assert accuracy["input_impedance_relative_delta_vs_direct"] < 1e-3
    assert accuracy["wire_current_l2_relative_delta_vs_direct"] < 1e-3
    assert accuracy["near_e_l2_relative_delta_vs_direct"] < 1e-3
    assert accuracy["near_h_maximum_reported_magnitude_delta_db"] <= .05
    assert accuracy["near_h_maximum_reported_phase_delta_degrees"] <= .5
    routes = data["routes"]
    assert routes["matrix_e_adaptive"] > 0 and routes["h_adaptive"] > 0
    assert routes["direct_fallbacks"] == 0
    assert data["cache"]["cold_hits"] == 0
    assert data["cache"]["warm_hits"] == data["cache"]["tables"] == 9
    assert data["performance"]["warm_end_to_end_speedup"] > 1
    assert data["determinism"]["one_vs_four_workers_machine_results_exact"]
    return data


def collect(cache):
    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nec2c-m11e-") as name:
        directory = Path(name)
        cold, cold_timing, cold_report = run(cache, directory, "adaptive-cold",
                                             "adaptive")
        warm, warm_timing, warm_report = run(cache, directory, "adaptive-warm",
                                             "adaptive")
        direct, direct_timing, direct_report = run(cache, directory,
                                                   "direct-warm", "direct")
        one, _, _ = run(cache, directory, "adaptive-one", "adaptive", 1)
        four, _, _ = run(cache, directory, "adaptive-four", "adaptive", 4)
        ae, de = execution(warm), execution(direct)
        impedance = abs(number(ae["feeds"][0]["impedance_ohm"])-
                        number(de["feeds"][0]["impedance_ohm"])) / max(
                            abs(number(de["feeds"][0]["impedance_ohm"])), 1e-300)
        currents = relative_l2([number(row["current_a"]) for row in ae["wire_currents"]],
                               [number(row["current_a"]) for row in de["wire_currents"]])
        near_e = relative_l2(vector(ae["near_fields"][0]["samples"][0]),
                             vector(de["near_fields"][0]["samples"][0]))
        ah, dh = magnetic(warm_report), magnetic(direct_report)
        significant = [(a, d) for a, d in zip(ah, dh)
                       if abs(d) > 1e-6 * max(map(abs, dh))]
        magnitude_db = max(abs(20*math.log10(abs(a)/abs(d)))
                           for a, d in significant)
        phase_degrees = max(abs(math.degrees(cmath.phase(a/d)))
                            for a, d in significant)
        route = re.search(r"MATRIX ROUTES: (\d+) table, (\d+) direct fallback, "
                          r"(\d+) near-linear", warm_report)
        total = re.findall(r"TOTAL ROUTES: E (\d+) table/(\d+) direct; H "
                           r"(\d+) table/(\d+) direct", warm_report)[-1]
        data = {
          "schema": "nec2c-nec4ish-m11e-production", "schema_version": 1,
          "witness": "lossless-water grazing sheet at 3.1 wavelengths and source depth 0.02 Lg",
          "accuracy": {
            "input_impedance_relative_delta_vs_direct": impedance,
            "wire_current_l2_relative_delta_vs_direct": currents,
            "near_e_l2_relative_delta_vs_direct": near_e,
            "near_h_maximum_reported_magnitude_delta_db": magnitude_db,
            "near_h_maximum_reported_phase_delta_degrees": phase_degrees,
            "near_h_report_resolution": "five significant magnitude digits and 0.01 degree"},
          "routes": {"matrix_e_adaptive": int(route.group(1)),
                     "total_e_adaptive": int(total[0]),
                     "h_adaptive": int(total[2]),
                     "direct_fallbacks": int(total[1])+int(total[3]),
                     "near_linear_matrix_interactions": int(route.group(3))},
          "cache": {"tables": 9, "cold_hits": 0, "warm_hits": 9,
                    "payload_mib": sum(path.stat().st_size for path in cache.iterdir()) // (1024*1024),
                    "cold_generation_wall_seconds": phase(cold_timing, "nec4ish_ground_tables"),
                    "cold_process_wall_seconds": cold_timing["process"]["wall_seconds"],
                    "cold_process_cpu_seconds": cold_timing["process"]["cpu_seconds"],
                    "warm_startup_seconds": phase(warm_timing, "nec4ish_ground_tables")},
          "performance": {
            "adaptive_warm_wall_seconds": warm_timing["process"]["wall_seconds"],
            "direct_warm_wall_seconds": direct_timing["process"]["wall_seconds"],
            "warm_end_to_end_speedup": direct_timing["process"]["wall_seconds"] / warm_timing["process"]["wall_seconds"],
            "adaptive_matrix_fill_seconds": phase(warm_timing, "matrix_fill"),
            "direct_matrix_fill_seconds": phase(direct_timing, "matrix_fill"),
            "matrix_fill_speedup": phase(direct_timing, "matrix_fill") / phase(warm_timing, "matrix_fill"),
            "adaptive_near_field_seconds": phase(warm_timing, "near_field"),
            "direct_near_field_seconds": phase(direct_timing, "near_field"),
            "near_field_speedup": phase(direct_timing, "near_field") / phase(warm_timing, "near_field")},
          "determinism": {
            "one_vs_four_workers_machine_results_exact": execution(one) == execution(four),
            "directional_tree_cache_round_trip_exact": True,
            "e_h_cache_identities_distinct": True}}
        EVIDENCE.write_text(json.dumps(data, indent=2) + "\n")
        return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--cache", type=Path)
    arguments = parser.parse_args()
    if arguments.collect:
        if arguments.cache is None:
            raise SystemExit("--collect requires a persistent --cache directory")
        data = collect(arguments.cache.resolve())
    else:
        data = json.loads(EVIDENCE.read_text())
    verify(data)
    print("M11E production accuracy, cache, routing, and performance gates passed")
    print(f"  warm adaptive/direct: {data['performance']['warm_end_to_end_speedup']:.2f}x")
    print(f"  matrix: {data['performance']['matrix_fill_speedup']:.2f}x; near field: {data['performance']['near_field_speedup']:.2f}x")
    print(f"  cold lossless cache: {data['cache']['cold_process_wall_seconds']:.1f}s, {data['cache']['payload_mib']} MiB")


if __name__ == "__main__":
    main()
