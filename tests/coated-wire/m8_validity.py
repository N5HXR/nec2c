#!/usr/bin/env python3
"""M8 validity-envelope, retained-reference, and warning-policy checks."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT/"nec2c")).resolve()
CVEL = 299.8


def deck(epsilon_r, radius_m, frequency_mhz, count=1, step_mhz=0):
    return "\n".join((
        "CM M8 VALIDITY POLICY", "CE",
        "GW 1 95 0 0 -0.2375 0 0 0.2375 0.002", "GE 0",
        f"LD 7 1 0 0 {epsilon_r:.17g} {radius_m:.17g} 0 0",
        f"FR 0 {count} 0 0 {frequency_mhz:.17g} {step_mhz:.17g}",
        "EX 0 1 48 0 1 0", "XQ", "EN", ""))


def run(name, text):
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m8-test-") as temp:
        directory = Path(temp)
        source = directory/f"{name}.nec"
        report = directory/f"{name}.out"
        result = directory/f"{name}.json"
        source.write_text(text, encoding="ascii")
        completed = subprocess.run(
            (str(BINARY), "--result-json", str(result), "-i", str(source),
             "-o", str(report)), cwd=ROOT, text=True, capture_output=True)
        if completed.returncode:
            raise RuntimeError(completed.stdout+completed.stderr)
        return json.loads(result.read_text())


def evaluation(result, index=0):
    return result["jobs"][0]["executions"][index][
        "material_evaluations"][0]


def close(left, right, absolute):
    if abs(left-right) > absolute:
        raise AssertionError((left, right, absolute))


def main():
    below = run("below", deck(3.2, 0.006, 265))
    below_evaluation = evaluation(below)
    if below_evaluation["validity_status"] != "within-reference-envelope":
        raise AssertionError(below_evaluation)
    if any(row["code"] == "COATED_WIRE_VALIDITY_ENVELOPE_EXCEEDED"
           for row in below["diagnostics"]):
        raise AssertionError("in-envelope case emitted a validity warning")

    above = run("above", deck(10, 0.012, 210, count=3, step_mhz=15))
    warnings = [row for row in above["diagnostics"]
                if row["code"] == "COATED_WIRE_VALIDITY_ENVELOPE_EXCEEDED"]
    if len(warnings) != 1:
        raise AssertionError(("expected one warning for the sweep", warnings))
    if any(evaluation(above, index)["validity_status"] !=
           "warning-outside-reference-envelope" for index in range(3)):
        raise AssertionError("outside-envelope execution was not stamped")

    boundary_frequency = 0.12*CVEL/(math.sqrt(3.2)*2*math.pi*0.006)
    lower = evaluation(run("lower-boundary", deck(
        3.2, 0.006, boundary_frequency*(1-1e-6))))
    upper_result = run("upper-boundary", deck(
        3.2, 0.006, boundary_frequency*(1+1e-6)))
    upper = evaluation(upper_result)
    if lower["validity_status"] != "within-reference-envelope":
        raise AssertionError(lower)
    if upper["validity_status"] != "warning-outside-reference-envelope":
        raise AssertionError(upper)

    evidence = json.loads((ROOT/
        "tests/coated-wire/m8-openems-thickness.json").read_text())
    analysis = json.loads((ROOT/
        "tests/coated-wire/m8-thickness-analysis.json").read_text())
    variants = json.loads((ROOT/
        "tests/coated-wire/m8-variant-evidence.json").read_text())
    geometry = json.loads((ROOT/
        "tests/coated-wire/m8-geometry-analysis.json").read_text())
    ground_convergence = json.loads((ROOT/
        "tests/coated-wire/m8-groundplane45-convergence.json").read_text())
    boundary = analysis["strict_one_percent_boundary"]
    if not (0.10 < boundary["largest_passing_sqrt_er_k0_b"] < 0.12):
        raise AssertionError(boundary)
    if not (0.12 < boundary["first_failing_sqrt_er_k0_b"] < 0.14):
        raise AssertionError(boundary)
    fit = analysis["affine_fit_frequency_error_vs_sqrt_er_k0_b"]
    if fit["r_squared"] < 0.97:
        raise AssertionError(fit)
    radiation = analysis["radiation_metrics"]
    if radiation["maximum_absolute_directivity_error_db"] > 0.1:
        raise AssertionError(radiation)
    if radiation["maximum_pattern_error_db_above_minus_40_db"] > 1.0:
        raise AssertionError(radiation)

    by_key = {row["key"]: row for row in evidence["cases"]}
    by_label = {row["label"]: row for row in variants["cases"]}
    coarse = by_key["er10-b12mm"]["resonance"]
    refined = by_label["refined-er10-b12"]["openems"]["resonance"]
    if abs(100*(refined["frequency_mhz"]/coarse["frequency_mhz"]-1)) > 0.02:
        raise AssertionError("refined OpenEMS frequency did not converge")
    if abs(100*(refined["resistance_ohm"]/coarse["resistance_ohm"]-1)) > 0.01:
        raise AssertionError("refined OpenEMS resistance did not converge")

    loss = by_label["loss-er3p2-b6"]["error"]
    if abs(loss["frequency_percent"]) > 1 or abs(loss["resistance_percent"]) > 1:
        raise AssertionError(loss)
    if abs(loss["efficiency_percentage_points"]) > 0.1:
        raise AssertionError(loss)

    current = by_label["current-er3p2-b12"]["current_comparison"]
    if current["near_feed_maximum_delta"] > 0.01:
        raise AssertionError(current)
    if current["interior_maximum_delta"] > 0.02:
        raise AssertionError(current)
    if current["near_end_maximum_delta"] < 0.05:
        raise AssertionError(("expected endpoint-sensitive current sample",
                              current))

    flush = by_key["er3.2-b12mm"]["resonance"]["frequency_mhz"]
    flat = by_label["flat-overhang-er3p2-b12"]["openems"]["resonance"][
        "frequency_mhz"]
    cap_shift = 100*(flat/flush-1)
    if abs(cap_shift) < 2:
        raise AssertionError(("endpoint sensitivity unexpectedly small",
                              cap_shift))
    close(analysis["strict_one_percent_boundary"][
        "provisional_warning_threshold"], 0.12, 0.0)

    geometry_by_name = {row["name"]: row for row in geometry["cases"]}
    expected_geometry = {"dogleg-b6", "dogleg-b12", "staged90-b6",
                         "staged90-b6-gap6-control", "groundplane45-b6"}
    if set(geometry_by_name) != expected_geometry:
        raise AssertionError(set(geometry_by_name))
    for name, row in geometry_by_name.items():
        effect = row["coating_effect"]
        gates = {
            "resonance": (abs(effect[
                "resonance_shift_error_percentage_points"]), 0.25),
            "resistance": (abs(effect[
                "resistance_change_error_ohm"]), 1.0),
            "directivity": (abs(effect[
                "directivity_change_error_db"]), 0.015),
            "efficiency": (abs(effect[
                "efficiency_change_error_percentage_points"]), 0.5),
            "pattern": (abs(effect["pattern_change_error_db"]), 0.25),
            "current": (abs(effect["current_change_error"]), 0.005),
        }
        failed = {metric: values for metric, values in gates.items()
                  if values[0] > values[1]}
        if failed:
            raise AssertionError((name, failed))
        if row["sqrt_er_k0_b"] > 0.12:
            raise AssertionError((name, row["sqrt_er_k0_b"]))
    staged = geometry_by_name["staged90-b6"]["coating_effect"]
    control = geometry_by_name[
        "staged90-b6-gap6-control"]["coating_effect"]
    if abs(staged["resonance_shift_error_percentage_points"]-
           control["resonance_shift_error_percentage_points"]) > 0.15:
        raise AssertionError("staged-bend coating result was feed-gap sensitive")
    ground = geometry_by_name["groundplane45-b6"]
    if ground["bare"]["openems"]["geometry"]["radial_slope_deg"] != -45:
        raise AssertionError("ground-plane radials are not 45 degrees down")
    meshes = [row["mesh_pitch_mm"] for row in ground_convergence["openems"]]
    if meshes != [2.0, 1.5, 1.0]:
        raise AssertionError(meshes)
    convergence = ground_convergence["interpretation"]
    if abs(convergence["finest_mesh_shift_error_percentage_points"]) > 0.25:
        raise AssertionError(convergence)
    if convergence["maximum_adjacent_mesh_shift_change_percentage_points"] < 0.8:
        raise AssertionError("ground-plane mesh uncertainty was lost")
    close(convergence["reference_uncertainty_percentage_points"], 0.5, 0.0)

    print("M8 warning boundary is exact and emits once per frequency sweep")
    print("M8 lossless validity fit, radiation metrics, and refined reference passed")
    print(f"M8 lossy errors {loss['frequency_percent']:.3f}%/"
          f"{loss['resistance_percent']:.3f}%; efficiency "
          f"{loss['efficiency_percentage_points']:.3f} percentage point")
    print(f"M8 extended-flat endpoint resonance shift {cap_shift:.3f}%")
    print("M8 normalized-current deltas: "
          f"feed {current['near_feed_maximum_delta']:.3f}, "
          f"interior {current['interior_maximum_delta']:.3f}, "
          f"end {current['near_end_maximum_delta']:.3f}")
    print("M8 refined dogleg, staged-90, source-control, and 45-degree "
          "ground-plane coating effects passed")
    print("M8 ground-plane Cartesian resonance-shift uncertainty retained "
          "at 0.5 percentage point")


if __name__ == "__main__":
    main()
