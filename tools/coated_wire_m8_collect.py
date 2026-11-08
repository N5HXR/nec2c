#!/usr/bin/env python3
"""Collect compact M8 NEC/OpenEMS thickness-map evidence."""

import argparse
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
C0 = 299792458.0


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", default=[],
                        metavar="ER,RADIUS_MM,DIRECTORY")
    parser.add_argument("--nec-map", action="append", default=[],
                        metavar="ER=PATH")
    parser.add_argument("--include-m1-b4", action="store_true")
    parser.add_argument("--output", type=Path, default=
                        ROOT/"tests/coated-wire/m8-openems-thickness.json")
    return parser.parse_args()


def resonance(rows):
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 == 0 or x0*x1 < 0:
            part = 0.0 if x0 == 0 else -x0/(x1-x0)
            return {
                "frequency_mhz": left["frequency_mhz"] + part*(
                    right["frequency_mhz"]-left["frequency_mhz"]),
                "resistance_ohm": left["impedance_ohm"]["re"] + part*(
                    right["impedance_ohm"]["re"]-
                    left["impedance_ohm"]["re"]),
            }
    raise RuntimeError("resonance is not bracketed")


def key(er, radius):
    return f"er{er:g}-b{radius:g}mm"


def pattern_error(nec_evaluation, openems_far):
    if not nec_evaluation or not openems_far:
        return None
    reference = {}
    for sample in openems_far["samples"]:
        if sample["phi_deg"] == 0:
            reference[sample["theta_deg"]] = sample["relative_db"]
    errors = []
    for sample in nec_evaluation["pattern"]:
        level = sample["relative_db"]
        expected = reference.get(sample["theta_deg"])
        if level is not None and expected is not None and expected >= -40:
            errors.append(abs(level-expected))
    return max(errors) if errors else None


def main():
    args = arguments()
    nec = {}
    for item in args.nec_map:
        er_text, path_text = item.split("=", 1)
        er = float(er_text)
        source = json.loads(Path(path_text).read_text())
        for row in source["cases"]:
            nec[(er, float(row["outer_radius_mm"]))] = row

    cases = []
    for item in args.case:
        er_text, radius_text, directory_text = item.split(",", 2)
        er, radius = float(er_text), float(radius_text)
        directory = Path(directory_text)
        spec = json.loads((directory/"spec.json").read_text())
        result = resonance(json.loads((directory/"impedance.json").read_text()))
        far_path = directory/"far-field.json"
        far = json.loads(far_path.read_text()) if far_path.exists() else None
        cases.append({
            "key": key(er, radius), "source": "openems-run",
            "epsilon_r": er, "outer_radius_mm": radius,
            "mesh": spec["mesh"], "wall_seconds": spec.get("wall_seconds"),
            "resonance": result,
            "directivity_db": far.get("directivity_db") if far else None,
            "reference_directory": str(directory),
        })

    if args.include_m1_b4:
        source = json.loads((ROOT/"tests/coated-wire/m1-openems-caps.json").read_text())
        row = next(item for item in source["cases"]
                   if item["label"] == "flush-8cell")
        cases.append({
            "key": key(3.2, 4), "source": "m1-retained-flush-8cell",
            "epsilon_r": 3.2, "outer_radius_mm": 4.0,
            "mesh": {"total_cells": row["mesh_cells"],
                     "coating_radial_cells": row["coating_radial_cells"]},
            "wall_seconds": None,
            "resonance": {"frequency_mhz": row["resonance_mhz"],
                          "resistance_ohm": row["resistance_at_resonance_ohm"]},
            "directivity_db": row["far_field"]["directivity_db"],
            "reference_directory": "tests/coated-wire/m1-openems-caps.json",
        })

    by_key = {row["key"]: row for row in cases}
    cases = [by_key[name] for name in sorted(by_key)]
    for case in cases:
        lookup = (case["epsilon_r"], case["outer_radius_mm"])
        if lookup not in nec:
            raise RuntimeError(f"missing NEC result for {lookup}")
        predicted = nec[lookup]
        nec_evaluation = predicted.get("resonance_evaluation")
        reference = case["resonance"]
        radius_m = case["outer_radius_mm"]/1000
        k0b = 2*math.pi*reference["frequency_mhz"]*1e6/C0*radius_m
        case.update({
            "thickness_mm": case["outer_radius_mm"]-2,
            "radius_ratio_b_over_a": case["outer_radius_mm"]/2,
            "k0_b_at_openems_resonance": k0b,
            "sqrt_er_k0_b": math.sqrt(case["epsilon_r"])*k0b,
            "nec": {"frequency_mhz": predicted["resonance_mhz"],
                    "resistance_ohm": predicted["resistance_ohm"],
                    "directivity_db": ((nec_evaluation or {}).get(
                        "directivity_db")),
                    "radiation_efficiency": ((nec_evaluation or {}).get(
                        "radiation_efficiency"))},
            "error": {
                "frequency_percent": 100*(predicted["resonance_mhz"]/
                                           reference["frequency_mhz"]-1),
                "resistance_ohm": (predicted["resistance_ohm"]-
                                   reference["resistance_ohm"]),
                "resistance_percent": 100*(predicted["resistance_ohm"]/
                                            reference["resistance_ohm"]-1),
                "directivity_db": (
                    None if not nec_evaluation or not far else
                    nec_evaluation["directivity_db"]-far["directivity_db"]),
                "maximum_pattern_db_above_minus_40_db": pattern_error(
                    nec_evaluation, far),
            },
        })

    record = {
        "schema": "nec2c-coated-wire-m8-thickness-evidence",
        "article": {"length_mm": 475, "conductor_radius_mm": 2,
                    "feed_gap_mm": 2, "coating_end": "flush",
                    "loss_tangent": 0},
        "interpretation": {
            "strict_reference_band_percent": 1.0,
            "candidate_engineering_band_percent": 5.0,
            "note": ("No cutoff is inferred until the reference grid spans "
                     "thickness, permittivity, loss, and end geometry."),
        },
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    for case in cases:
        print(f"{case['key']}: {case['error']['frequency_percent']:.3f}% "
              f"resonance, {case['error']['resistance_ohm']:.3f} ohm")
    print(args.output)


if __name__ == "__main__":
    main()
