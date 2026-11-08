#!/usr/bin/env python3
"""M6 topology checks for coated-wire junctions."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def deck(split=False, lossy=False, source_tag=1, source_segment=48,
         right="equal", sweep=True, half_length=0.24):
    geometry = ([f"GW 1 48 0 0 {-half_length:.17g} 0 0 0 0.002",
                 f"GW 2 48 0 0 0 0 0 {half_length:.17g} 0.002"] if split else
                [f"GW 1 96 0 0 {-half_length:.17g} 0 0 "
                 f"{half_length:.17g} 0.002"])
    loss = "0.015 0" if lossy else "0 0"
    coating = [f"LD 7 1 0 0 3.2 0.004 {loss}"]
    if split and right == "equal":
        coating.append(f"LD 7 2 0 0 3.2 0.004 {loss}")
    elif split and right == "unequal":
        coating.append(f"LD 7 2 0 0 2.2 0.003 {loss}")
    frequency = ("FR 0 101 0 0 250 0.5" if sweep else
                 "FR 0 1 0 0 280 0")
    return "\n".join([
        "CM M6 UNIFORM COATING JUNCTION", "CE", *geometry, "GE 0",
        *coating, frequency,
        f"EX 0 {source_tag} {source_segment} 0 1 0", "XQ", "EN", "",
    ])


def invoke(directory, name, text, experimental=False):
    source = directory / f"{name}.nec"
    report = directory / f"{name}.out"
    result = directory / f"{name}.json"
    source.write_text(text, encoding="ascii")
    command = [
        str(BINARY), "--backend", "builtin", "--jobs", "1", "--threads", "1",
        "--result-json", str(result), "-i", str(source), "-o", str(report)]
    if experimental:
        command[1:1] = ["--coated-wire-junctions", "experimental"]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    parsed = json.loads(result.read_text()) if result.exists() else None
    return completed, parsed


def run(directory, name, text, experimental=False):
    completed, parsed = invoke(directory, name, text, experimental)
    if completed.returncode:
        diagnostic = (parsed.get("diagnostics", [{}])[0].get("code")
                      if parsed else None)
        raise RuntimeError(f"{name}: {diagnostic}\n{completed.stdout}\n"
                           f"{completed.stderr}")
    return parsed


def complex_value(value):
    return complex(value["re"], value["im"])


def execution_rows(result):
    return result["jobs"][0]["executions"]


def resonance(result):
    rows = execution_rows(result)
    for left, right in zip(rows, rows[1:]):
        z0 = complex_value(left["feeds"][0]["impedance_ohm"])
        z1 = complex_value(right["feeds"][0]["impedance_ohm"])
        if z0.imag == 0 or z0.imag*z1.imag < 0:
            fraction = 0 if z0.imag == 0 else -z0.imag/(z1.imag-z0.imag)
            return (left["frequency_mhz"] + fraction*(right["frequency_mhz"]-
                                                       left["frequency_mhz"]),
                    z0.real + fraction*(z1.real-z0.real))
    raise AssertionError("resonance not bracketed")


def assert_same_solution(single, split, tolerance=2e-11):
    left_rows = execution_rows(single)
    right_rows = execution_rows(split)
    if len(left_rows) != len(right_rows):
        raise AssertionError("different sweep lengths")
    maximum = 0.0
    power_maximum = 0.0
    for left, right in zip(left_rows, right_rows):
        zl = complex_value(left["feeds"][0]["impedance_ohm"])
        zr = complex_value(right["feeds"][0]["impedance_ohm"])
        maximum = max(maximum, abs(zl-zr)/max(abs(zl), abs(zr), 1.0))
        for lc, rc in zip(left["wire_currents"], right["wire_currents"]):
            il = complex_value(lc["current_a"])
            ir = complex_value(rc["current_a"])
            maximum = max(maximum, abs(il-ir)/max(abs(il), abs(ir), 1.0e-12))
        for key in ("accepted_w", "air_radiated_w", "dielectric_absorbed_w",
                    "balance_defect_w", "balance_defect_relative"):
            lp = left["power"].get(key)
            rp = right["power"].get(key)
            if lp is not None or rp is not None:
                if lp is None or rp is None:
                    raise AssertionError(f"power field {key} changed presence")
                power_maximum = max(
                    power_maximum,
                    abs(lp-rp)/max(abs(lp), abs(rp), 1.0e-15))
    if maximum > tolerance:
        raise AssertionError(f"tag split changed solution by {maximum:.3e}")
    if power_maximum > 1.0e-5:
        raise AssertionError(f"tag split changed power by {power_maximum:.3e}")
    return maximum, power_maximum


def experimental_reciprocity(directory, kind):
    left = run(directory, f"{kind}-left",
               deck(split=True, right=kind, source_tag=1, source_segment=24,
                    sweep=False), experimental=True)
    right = run(directory, f"{kind}-right",
                deck(split=True, right=kind, source_tag=2, source_segment=25,
                     sweep=False), experimental=True)
    ilr = complex_value(left["jobs"][0]["executions"][0]
                        ["wire_currents"][72]["current_a"])
    irl = complex_value(right["jobs"][0]["executions"][0]
                        ["wire_currents"][23]["current_a"])
    defect = abs(ilr-irl)/max(abs(ilr), abs(irl), 1.0e-30)
    for result in (left, right):
        diagnostics = result.get("diagnostics", [])
        if not any(item["severity"] == "warning" and
                   item["code"] == "COATED_WIRE_EXPERIMENTAL_JUNCTION"
                   for item in diagnostics):
            raise AssertionError("experimental junction warning is missing")
        if (result["jobs"][0]["coated_wires"]["operator_status"] !=
                "experimental-junction"):
            raise AssertionError("experimental junction status is missing")
    return defect


def uniform_reciprocity(directory):
    left = run(directory, "equal-left",
               deck(split=True, source_tag=1, source_segment=24,
                    sweep=False))
    right = run(directory, "equal-right",
                deck(split=True, source_tag=2, source_segment=25,
                     sweep=False))
    ilr = complex_value(left["jobs"][0]["executions"][0]
                        ["wire_currents"][72]["current_a"])
    irl = complex_value(right["jobs"][0]["executions"][0]
                        ["wire_currents"][23]["current_a"])
    return abs(ilr-irl)/max(abs(ilr), abs(irl), 1.0e-30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        help="write the deterministic acceptance metrics")
    arguments = parser.parse_args()
    frozen_path = ROOT / "tests/coated-wire/m6-acceptance.json"
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m6-") as temporary:
        root = Path(temporary)
        single = run(root, "single", deck())
        split = run(root, "split", deck(split=True))
        single_lossy = run(root, "single-lossy", deck(lossy=True))
        split_lossy = run(root, "split-lossy", deck(split=True, lossy=True))
        lossless_error, lossless_power_error = assert_same_solution(
            single, split, 1.0e-8)
        lossy_error, lossy_power_error = assert_same_solution(
            single_lossy, split_lossy, 1.0e-8)
        strict_cases = {
            "bare": "COATED_WIRE_BARE_JUNCTION_UNSUPPORTED",
            "unequal": "COATED_WIRE_JUNCTION_MISMATCH",
        }
        for kind, expected in strict_cases.items():
            completed, parsed = invoke(root, f"strict-{kind}",
                                       deck(split=True, right=kind,
                                            sweep=False))
            observed = (parsed.get("diagnostics", [{}])[0].get("code")
                        if parsed else None)
            if completed.returncode == 0 or observed != expected:
                raise AssertionError((kind, expected, observed))
        equal_defect = uniform_reciprocity(root)
        if equal_defect > 1.0e-12:
            raise AssertionError(f"uniform junction reciprocity {equal_defect:.3e}")
        bare_defect = experimental_reciprocity(root, "bare")
        unequal_defect = experimental_reciprocity(root, "unequal")
        if bare_defect <= 1.0e-3 or unequal_defect <= 1.0e-3:
            raise AssertionError("an experimental transition crossed the M6 "
                                 "reciprocity rejection threshold; review policy")
        bare_candidate = resonance(run(
            root, "bare-candidate-sweep",
            deck(split=True, right="bare", half_length=0.2375),
            experimental=True))
        unequal_candidate = resonance(run(
            root, "unequal-candidate-sweep",
            deck(split=True, right="unequal", half_length=0.2375),
            experimental=True))

    if frozen_path.exists():
        frozen = json.loads(frozen_path.read_text(encoding="utf-8"))
        observed = frozen["uniform_junction"]
        for key, value in {
                "maximum_relative_delta_lossless": lossless_error,
                "maximum_relative_delta_lossy": lossy_error,
                "reciprocity_relative_defect": equal_defect}.items():
            if not math.isclose(observed[key], value, rel_tol=1e-9,
                                abs_tol=1e-15):
                raise AssertionError((key, observed[key], value))
        for reference in frozen["openems_asymmetric_references"]:
            if reference["end_criteria"] != 1.0e-4:
                raise AssertionError("M6 OpenEMS end criterion changed")
            if reference["mesh"]["minimum_steps"][2] < 0.5:
                raise AssertionError("M6 OpenEMS axial sliver cell returned")
            if abs(reference["experimental_nec_candidate_error"]
                   ["resistance_ohm"]) <= 2.0:
                raise AssertionError("experimental junction unexpectedly met "
                                     "the M3 resistance budget")

    status = split["jobs"][0]["coated_wires"]["operator_status"]
    lossy_status = split_lossy["jobs"][0]["coated_wires"]["operator_status"]
    if status != "lossless-uniform-junction":
        raise AssertionError(status)
    if lossy_status != "lossy-uniform-junction":
        raise AssertionError(lossy_status)
    frequency, resistance = resonance(split)
    if not (math.isfinite(frequency) and math.isfinite(resistance)):
        raise AssertionError("nonfinite split-tag resonance")
    print("M6 identical coated-tag junction preserves the single-tag solution")
    print(f"M6 maximum lossless relative delta {lossless_error:.3e}; "
          f"lossy {lossy_error:.3e}")
    print(f"M6 maximum power delta: lossless {lossless_power_error:.3e}, "
          f"lossy {lossy_power_error:.3e}")
    print(f"M6 split-tag resonance {frequency:.6f} MHz/{resistance:.6f} ohm")
    print(f"M6 identical-junction reciprocity defect {equal_defect:.3e}")
    print(f"M6 experimental reciprocity defects: coated/bare {bare_defect:.3e}, "
          f"unequal coating {unequal_defect:.3e}")
    print(f"M6 experimental candidate resonances: coated/bare "
          f"{bare_candidate[0]:.6f} MHz, unequal {unequal_candidate[0]:.6f} MHz")
    print("M6 strict mode retains stable coated/bare and unequal-coating errors")
    if arguments.output:
        record = {
            "schema": "nec2c-coated-wire-m6-junctions",
            "schema_version": 1,
            "uniform_junction": {
                "maximum_relative_delta_lossless": lossless_error,
                "maximum_relative_delta_lossy": lossy_error,
                "maximum_relative_power_delta_lossless": lossless_power_error,
                "maximum_relative_power_delta_lossy": lossy_power_error,
                "reciprocity_relative_defect": equal_defect,
                "resonance_mhz": frequency,
                "resistance_ohm": resistance,
                "operator_status_lossless": status,
                "operator_status_lossy": lossy_status,
            },
            "experimental_candidates": {
                "coated_bare_reciprocity_relative_defect": bare_defect,
                "unequal_coating_reciprocity_relative_defect": unequal_defect,
                "coated_bare_resonance_mhz": bare_candidate[0],
                "coated_bare_resistance_ohm": bare_candidate[1],
                "unequal_coating_resonance_mhz": unequal_candidate[0],
                "unequal_coating_resistance_ohm": unequal_candidate[1],
            },
            "decision": {
                "identical_coating_junction": "supported",
                "coated_bare_junction": "rejected-by-default",
                "unequal_coating_junction": "rejected-by-default",
                "experimental_cli": "--coated-wire-junctions experimental",
            },
        }
        arguments.output.write_text(json.dumps(record, indent=2,
                                                sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
