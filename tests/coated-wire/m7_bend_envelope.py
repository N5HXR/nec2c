#!/usr/bin/env python3
"""Retained-evidence checks for the M7 practical bend envelope."""

import json
import math
import os
from pathlib import Path


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()


def length(points):
    return sum(math.dist(left, right) for left, right in zip(points, points[1:]))


def main():
    study = json.loads((ROOT/"tests/coated-wire/m7-bend-study.json").read_text())
    analysis = json.loads(
        (ROOT/"tests/coated-wire/m7-bend-analysis.json").read_text())
    openems = json.loads(
        (ROOT/"tests/coated-wire/m7-bend-openems.json").read_text())
    cases = {row["name"]: row for row in study["cases"]}
    for row in cases.values():
        if abs(length(row["points_m"])-0.474) > 1e-12:
            raise AssertionError(f"centerline length changed: {row['name']}")
    comparison = analysis["nec_matlab"]
    if comparison["maximum_absolute_frequency_error_percent"] > 1.0:
        raise AssertionError("NEC/MATLAB bare resonance agreement regressed")
    if comparison["maximum_absolute_resistance_error_ohm"] > 1.0:
        raise AssertionError("NEC/MATLAB bare resistance agreement regressed")
    eligible = analysis["x45_y4_eligible_cases"]
    if not eligible:
        raise AssertionError("the proposed X=45, Y=4 envelope has no cases")
    if any(row["maximum_turn_deg"] > 45 or
           row["minimum_chord_ratio"] < 4 for row in eligible):
        raise AssertionError("invalid case entered X=45, Y=4 envelope")

    runs = {row["key"]: row for row in openems["runs"]}
    sharp = [runs[name] for name in (
        "bare-2mm-sharp", "bare-1mm-sharp", "bare-0.5x1mm-sharp")]
    frequencies = [row["resonance"]["frequency_mhz"] for row in sharp]
    if not frequencies[0] < frequencies[1] < frequencies[2]:
        raise AssertionError("sharp-bend OpenEMS refinement is not monotonic")
    nec = cases["sharp"]["nec"]["bare"]
    matlab = next(row for row in comparison["comparisons"]
                  if row["name"] == "sharp")
    fine = sharp[-1]["resonance"]
    if abs(100*(fine["frequency_mhz"]/nec["frequency_mhz"]-1)) > 2.0:
        raise AssertionError("fine OpenEMS/NEC resonance exceeds relaxed band")
    if abs(100*(fine["frequency_mhz"]/
               matlab["matlab_frequency_mhz"]-1)) > 2.0:
        raise AssertionError("fine OpenEMS/MATLAB resonance exceeds relaxed band")
    if abs(100*(fine["resistance_ohm"]/nec["resistance_ohm"]-1)) > 5.0:
        raise AssertionError("fine OpenEMS/NEC resistance exceeds relaxed band")

    rounded = analysis["rounded_n1_r4_validation"]
    pairs = sorted(rounded["matched_mesh_pairs"],
                   key=lambda row: row["mesh_mm"], reverse=True)
    if [row["mesh_mm"] for row in pairs] != [2.0, 1.0]:
        raise AssertionError("rounded matched-mesh OpenEMS evidence is incomplete")
    errors = [abs(row["coating_frequency_shift_percent"]-
                  rounded["nec_coating_frequency_shift_percent"])
              for row in pairs]
    if not errors[1] < errors[0]:
        raise AssertionError("rounded coating shift did not converge toward NEC")
    if errors[1] > 1.0:
        raise AssertionError("rounded coating shift exceeds high-confidence band")
    rounded_nec = cases["arc-n1-r4b"]["nec"]["bare"]
    rounded_fine = rounded["finest_bare_openems"]
    if abs(100*(rounded_fine["frequency_mhz"]/
               rounded_nec["frequency_mhz"]-1)) > 2.0:
        raise AssertionError("rounded fine OpenEMS/NEC resonance exceeds band")
    if abs(100*(rounded_fine["resistance_ohm"]/
               rounded_nec["resistance_ohm"]-1)) > 10.0:
        raise AssertionError("rounded fine OpenEMS/NEC resistance exceeds band")

    shifts = [row["nec_coating_shift_percent"] for row in eligible]
    print("M7 bend centerlines preserve 474 mm length across "
          f"{len(cases)} cases")
    print("M7 NEC/MATLAB maximum differences: "
          f"{comparison['maximum_absolute_frequency_error_percent']:.3f}% and "
          f"{comparison['maximum_absolute_resistance_error_ohm']:.3f} ohm")
    print("M7 sharp OpenEMS convergence: " + ", ".join(
        f"{row['mesh_mm']:g}x{row['z_mesh_mm']:g}mm="
        f"{row['resonance']['frequency_mhz']:.3f}MHz/"
        f"{row['resonance']['resistance_ohm']:.3f}ohm" for row in sharp))
    print(f"M7 X=45/Y=4 NEC coating shifts span "
          f"{min(shifts):.3f}% to {max(shifts):.3f}%")
    print("M7 rounded n1/r4 OpenEMS coating shifts: " + ", ".join(
        f"{row['mesh_mm']:g}mm="
        f"{row['coating_frequency_shift_percent']:.3f}%"
        for row in pairs))
    print("M7 rounded n1/r4 finest coating-shift disagreement: "
          f"{errors[1]:.3f} percentage points")


if __name__ == "__main__":
    main()
