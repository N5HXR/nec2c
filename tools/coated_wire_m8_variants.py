#!/usr/bin/env python3
"""Collect M8 mesh, loss, and endpoint sensitivity evidence."""

import argparse
import bisect
import json
import math
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", default=[],
                        metavar="LABEL,OPENEMS_DIRECTORY,NEC_MAP")
    parser.add_argument("--binary", type=Path, default=ROOT/"nec2c")
    parser.add_argument("--output", type=Path, default=
                        ROOT/"tests/coated-wire/m8-variant-evidence.json")
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
                "accepted_power_w": left["accepted_power_w"] + part*(
                    right["accepted_power_w"]-left["accepted_power_w"]),
            }
    raise RuntimeError("resonance is not bracketed")


def nec_current_profile(binary, spec, nec_map, frequency_mhz):
    article = nec_map["article"]
    radius_m = spec["geometry_mm"]["coating_outer_radius"]/1000
    text = "\n".join((
        "CM M8 CURRENT PROFILE", "CE",
        "GW 1 95 0 0 -0.2375 0 0 0.2375 0.002", "GE 0",
        f"LD 7 1 0 0 {article['epsilon_r']:.17g} {radius_m:.17g} "
        f"{article.get('tan_delta', 0):.17g} "
        f"{article.get('conductivity_s_per_m', 0):.17g}",
        f"FR 0 1 0 0 {frequency_mhz:.17g} 0",
        "EX 0 1 48 0 1 0", "XQ", "EN", ""))
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m8-current-") as temp:
        directory = Path(temp)
        source, report, result = (directory/f"current{suffix}"
                                  for suffix in (".nec", ".out", ".json"))
        source.write_text(text, encoding="ascii")
        completed = subprocess.run(
            (str(binary), "--result-json", str(result), "-i", str(source),
             "-o", str(report)), cwd=ROOT, text=True, capture_output=True)
        if completed.returncode:
            raise RuntimeError(completed.stdout+completed.stderr)
        execution = json.loads(result.read_text())["jobs"][0]["executions"][0]
    rows = []
    for row in execution["wire_currents"]:
        value = row["current_a"]
        rows.append((1000*row["center_m"][2],
                     math.hypot(value["re"], value["im"])))
    peak = max(value for _, value in rows)
    return [(-237.5, 0.0), *[(z, value/peak) for z, value in rows],
            (237.5, 0.0)]


def interpolate(rows, position):
    positions = [row[0] for row in rows]
    right = bisect.bisect_left(positions, position)
    if right == 0:
        return rows[0][1]
    if right == len(rows):
        return rows[-1][1]
    left = right-1
    part = ((position-rows[left][0]) /
            (rows[right][0]-rows[left][0]))
    return rows[left][1]+part*(rows[right][1]-rows[left][1])


def current_comparison(binary, directory, spec, nec_map, frequency_mhz):
    path = directory/"axial-current.json"
    if not path.exists():
        return None
    reference = json.loads(path.read_text())["samples"]
    nec = nec_current_profile(binary, spec, nec_map, frequency_mhz)
    samples = []
    for row in reference:
        expected = row["relative_magnitude"]
        predicted = interpolate(nec, row["z_mm"])
        samples.append({"z_mm": row["z_mm"], "openems": expected,
                        "nec": predicted, "absolute_delta":
                        abs(predicted-expected)})
    half = spec["geometry_mm"]["length"]/2
    def region_max(predicate):
        values = [row["absolute_delta"] for row in samples
                  if predicate(abs(row["z_mm"])/half)]
        return max(values) if values else None
    return {
        "samples": samples,
        "maximum_normalized_delta": max(
            row["absolute_delta"] for row in samples),
        "near_feed_maximum_delta": region_max(lambda value: value <= 0.25),
        "interior_maximum_delta": region_max(
            lambda value: 0.25 < value < 0.75),
        "near_end_maximum_delta": region_max(lambda value: value >= 0.75),
    }


def main():
    args = arguments()
    cases = []
    for item in args.case:
        label, directory_text, nec_text = item.split(",", 2)
        directory = Path(directory_text)
        spec = json.loads((directory/"spec.json").read_text())
        impedance = json.loads((directory/"impedance.json").read_text())
        far = json.loads((directory/"far-field.json").read_text())
        reference = resonance(impedance)
        nec_map = json.loads(Path(nec_text).read_text())
        radius = spec["geometry_mm"]["coating_outer_radius"]
        nec = next(row for row in nec_map["cases"]
                   if row["outer_radius_mm"] == radius)
        nec_evaluation = nec.get("resonance_evaluation")
        openems_efficiency = (far["radiated_power_w"] /
                              reference["accepted_power_w"])
        currents = current_comparison(args.binary, directory, spec, nec_map,
                                      reference["frequency_mhz"])
        record = {
            "label": label,
            "openems_directory": str(directory),
            "nec_map": str(Path(nec_text)),
            "geometry_mm": spec["geometry_mm"],
            "material_parameters": spec["material_parameters"],
            "mesh_controls_mm": spec["mesh_controls_mm"],
            "mesh": spec["mesh"],
            "wall_seconds": spec.get("wall_seconds"),
            "openems": {
                "resonance": reference,
                "radiation_efficiency": openems_efficiency,
                "directivity_db": far["directivity_db"],
                "axisymmetric_postprocess": far.get(
                    "axisymmetric_postprocess"),
            },
            "nec": {
                "resonance_mhz": nec["resonance_mhz"],
                "resistance_ohm": nec["resistance_ohm"],
                "radiation_efficiency": (nec_evaluation or {}).get(
                    "radiation_efficiency"),
            },
            "error": {
                "frequency_percent": 100*(nec["resonance_mhz"] /
                                           reference["frequency_mhz"]-1),
                "resistance_percent": 100*(nec["resistance_ohm"] /
                                            reference["resistance_ohm"]-1),
                "efficiency_percentage_points": (
                    None if nec_evaluation is None else
                    100*(nec_evaluation["radiation_efficiency"]-
                         openems_efficiency)),
            },
        }
        if currents is not None:
            record["current_comparison"] = currents
        cases.append(record)

    output = {
        "schema": "nec2c-coated-wire-m8-variant-evidence",
        "interpretation": {
            "mesh_refinement": ("A small change confirms the OpenEMS "
                                "reference rather than expanding the NEC "
                                "validity envelope."),
            "loss": ("Constant NEC conductivity is used to match OpenEMS's "
                     "frequency-independent material kappa."),
            "endpoints": ("Non-flush cases measure physical endpoint-model "
                          "spread, not reduced-operator error."),
        },
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2)+"\n")
    for case in cases:
        error = case["error"]
        print(f"{case['label']}: f {error['frequency_percent']:.3f}%, "
              f"R {error['resistance_percent']:.3f}%", end="")
        if error["efficiency_percentage_points"] is not None:
            print(f", efficiency {error['efficiency_percentage_points']:.3f} pp")
        else:
            print()
    print(args.output)


if __name__ == "__main__":
    main()
