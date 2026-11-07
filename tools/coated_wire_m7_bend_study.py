#!/usr/bin/env python3
"""Run NEC and prepare MATLAB cases for the M7 bend-envelope study."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

from coated_wire_bend_geometry import (COATING_RADIUS_M, WIRE_RADIUS_M,
                                       insert_feed_point, study_cases)


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT/"nec2c"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT/"tests/coated-wire/m7-bend-study.json")
    parser.add_argument("--write-matlab-config", action="store_true")
    parser.add_argument("--run-nec", action="store_true")
    return parser.parse_args()


def segment_count(length):
    # About 5 mm per NEC segment, with an odd count where possible.
    count = max(1, round(length/0.005))
    return count


def nec_deck(case, coated):
    points = case["points_m"]
    geometry = []
    source = None
    tags = []
    for tag, (left, right) in enumerate(zip(points, points[1:]), 1):
        length = math.dist(left, right)
        count = segment_count(length)
        horizontal_feed = (abs(left[2]) < 1e-12 and
                           abs(right[2]) < 1e-12 and left[0] < 0 < right[0])
        if horizontal_feed and count % 2 == 0:
            count += 1
        geometry.append(
            "GW %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g" %
            (tag, count, *left, *right, WIRE_RADIUS_M))
        tags.append(tag)
        if horizontal_feed:
            source = (tag, (count+1)//2)
    if source is None:
        raise RuntimeError(f"{case['name']} has no center feed segment")
    lines = ["CM M7 BEND ENVELOPE", "CE", *geometry, "GE 0"]
    if coated:
        lines.extend(f"LD 7 {tag} 0 0 3.2 {COATING_RADIUS_M:.17g} 0 0"
                     for tag in tags)
    lines.extend(("FR 0 241 0 0 240 0.5",
                  f"EX 0 {source[0]} {source[1]} 0 1 0", "XQ", "EN", ""))
    return "\n".join(lines)


def resonance(result):
    rows = result["jobs"][0]["executions"]
    for left, right in zip(rows, rows[1:]):
        zl = complex(**{("real" if key == "re" else "imag"): value
                        for key, value in left["feeds"][0]["impedance_ohm"].items()})
        zr = complex(**{("real" if key == "re" else "imag"): value
                        for key, value in right["feeds"][0]["impedance_ohm"].items()})
        if zl.imag == 0 or zl.imag*zr.imag < 0:
            part = 0.0 if zl.imag == 0 else -zl.imag/(zr.imag-zl.imag)
            return {
                "frequency_mhz": left["frequency_mhz"] + part*(
                    right["frequency_mhz"]-left["frequency_mhz"]),
                "resistance_ohm": zl.real+part*(zr.real-zl.real),
            }
    raise RuntimeError("resonance not bracketed")


def run_nec(case, coated, directory):
    name = f"{case['name']}-{'coated' if coated else 'bare'}"
    source, report, result = (directory/f"{name}{suffix}"
                              for suffix in (".nec", ".out", ".json"))
    source.write_text(nec_deck(case, coated), encoding="ascii")
    command = [str(BINARY), "--backend", "auto", "--jobs", "1",
               "--threads", "1", "--result-json", str(result)]
    if coated:
        command.extend(("--coated-wire-junctions", "experimental"))
    command.extend(("-i", str(source), "-o", str(report)))
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
    return resonance(json.loads(result.read_text()))


def matlab_selection(cases):
    # Broad enough to identify angle/radius trends without spending hours on
    # thirty-one surface-MoM meshes. OpenEMS refinement is selected later at
    # the observed pass/fail boundary.
    wanted = {"sharp", "arc-n1-r2b", "arc-n1-r4b", "arc-n1-r8b",
              "arc-n2-r4b", "arc-n2-r8b", "arc-n3-r4b",
              "arc-n4-r4b", "arc-n6-r8b"}
    return [case for case in cases if case["name"] in wanted]


def main():
    args = arguments()
    cases = study_cases()
    record = {
        "schema": "nec2c-coated-wire-m7-bend-study",
        "wire_radius_m": WIRE_RADIUS_M,
        "coating_radius_m": COATING_RADIUS_M,
        "constant_centerline_length_m": 0.474,
        "cases": cases,
    }
    if args.run_nec:
        with tempfile.TemporaryDirectory(prefix="nec2c-m7-bends-") as path:
            directory = Path(path)
            for case in cases:
                case["nec"] = {
                    "bare": run_nec(case, False, directory),
                    "coated": run_nec(case, True, directory),
                }
                print(case["name"], case["nec"])
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    if args.write_matlab_config:
        selected = matlab_selection(cases)
        config = {
            "operation": "study",
            "output_filename": "m7-bend-matlab-result.json",
            "conductor_radius_m": WIRE_RADIUS_M,
            "maximum_edge_m": 0.002,
            "frequency_mhz": list(range(260, 341, 10)),
            "cases": [{"name": case["name"],
                       "points_m": insert_feed_point(case["points_m"])}
                      for case in selected],
        }
        path = ROOT/"tests/coated-wire/m7-bend-matlab-config.json"
        path.write_text(json.dumps(config, indent=2)+"\n")
        print(path)
    print(args.output)


if __name__ == "__main__":
    main()
