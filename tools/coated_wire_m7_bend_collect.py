#!/usr/bin/env python3
"""Retain compact records from the expensive M7 bend OpenEMS runs."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/"tests/coated-wire/m7-bend-openems.json"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sharp-2mm", type=Path)
    parser.add_argument("--sharp-1mm", type=Path)
    parser.add_argument("--sharp-half-by-1mm", type=Path)
    parser.add_argument("--rounded-bare-2mm", type=Path)
    parser.add_argument("--rounded-bare-1mm", type=Path)
    parser.add_argument("--rounded-bare-half-mm", type=Path)
    parser.add_argument("--rounded-coated-2mm", type=Path)
    parser.add_argument("--rounded-coated-1mm", type=Path)
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
    raise RuntimeError("resonance not bracketed")


def collect(key, directory, stage="sharp-convergence"):
    spec = json.loads((directory/"spec.json").read_text())
    bend = spec.get("bend", {})
    return {
        "key": key, "status": "ok", "stage": stage,
        "material": spec["material"], "mesh_mm": spec["mesh_pitch_mm"],
        "z_mesh_mm": spec.get("z_mesh_pitch_mm", spec["mesh_pitch_mm"]),
        "chords_per_corner": bend.get("chords_per_corner", 0),
        "radius_ratio": bend.get("radius_ratio", 0),
        "maximum_turn_deg": bend.get("maximum_turn_deg", 90),
        "minimum_chord_ratio": bend.get("minimum_chord_ratio", 0),
        "mesh": spec["mesh"], "solver_wall_seconds": spec.get("wall_seconds"),
        "resonance": resonance(json.loads(
            (directory/"impedance.json").read_text())),
    }


def main():
    args = arguments()
    record = (json.loads(OUTPUT.read_text()) if OUTPUT.exists() else
              {"schema": "nec2c-coated-wire-m7-bend-openems", "runs": []})
    by_key = {row["key"]: row for row in record["runs"]}
    entries = (("bare-2mm-sharp", args.sharp_2mm),
               ("bare-1mm-sharp", args.sharp_1mm),
               ("bare-0.5x1mm-sharp", args.sharp_half_by_1mm))
    for key, path in entries:
        if path is not None:
            by_key[key] = collect(key, path)
    rounded = (
        ("bare-2mm-n1-r4", args.rounded_bare_2mm),
        ("bare-1mm-n1-r4", args.rounded_bare_1mm),
        ("bare-0.5mm-n1-r4", args.rounded_bare_half_mm),
        ("coated-2mm-n1-r4", args.rounded_coated_2mm),
        ("coated-1mm-n1-r4", args.rounded_coated_1mm),
    )
    for key, path in rounded:
        if path is not None:
            by_key[key] = collect(key, path, "rounded-n1-r4-validation")
    # Replace the two early pilot names with the canonical validation names.
    if args.rounded_bare_2mm is not None:
        by_key.pop("bare-2mm-n1-r4b", None)
    if args.rounded_bare_1mm is not None:
        by_key.pop("bare-1mm-n1-r4b", None)
    record["runs"] = sorted(by_key.values(), key=lambda row: row["key"])
    OUTPUT.write_text(json.dumps(record, indent=2)+"\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
