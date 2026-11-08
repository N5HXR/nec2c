#!/usr/bin/env python3
"""M4 acceptance tests for lossy coated-wire materials and power."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
EPSILON0 = 8.8541878128e-12


def deck(tan_delta=0.0, sigma=0.0, frequency=276.65,
         count=1, step=0.0, pattern=False):
    return "\n".join([
        "CM M4 LOSS TEST", "CE",
        "GW 1 95 0 0 -0.2375 0 0 0.2375 0.002", "GE 0",
        f"LD 7 1 0 0 3.2 0.004 {tan_delta:.17g} {sigma:.17g}",
        f"FR 0 {count} 0 0 {frequency:.17g} {step:.17g}",
        "EX 0 1 48 0 1 0",
        "RP 0 19 1 1000 0 0 5 0 0 0" if pattern else "XQ 0",
        "EN", "",
    ])


def run(directory, name, text, backend="builtin", threads=1, jobs=1):
    source = directory / f"{name}.nec"
    report = directory / f"{name}.out"
    result = directory / f"{name}.json"
    source.write_text(text, encoding="ascii")
    completed = subprocess.run(
        [str(BINARY), "--backend", backend, "--jobs", str(jobs),
         "--threads", str(threads), "--result-json", str(result),
         "-i", str(source), "-o", str(report)],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if completed.returncode:
        raise RuntimeError(f"{name} failed\n{completed.stdout}\n{completed.stderr}")
    return json.loads(result.read_text(encoding="utf-8"))


def close(left, right, relative=1e-12, absolute=1e-14):
    if not math.isclose(left, right, rel_tol=relative, abs_tol=absolute):
        raise AssertionError((left, right))


def impedance(execution):
    value = execution["feeds"][0]["impedance_ohm"]
    return complex(value["re"], value["im"])


def main():
    frequency = 276.65
    tan_delta = 0.05
    # Match NEC's historical TP constant so both cards evaluate to the same
    # binary complex permittivity, not merely the same rounded material value.
    equivalent_sigma = (6.283185308*frequency*1e6*EPSILON0*3.2*tan_delta)
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m4-") as temporary:
        root = Path(temporary)
        zero = run(root, "zero", deck())
        tangent = run(root, "tangent", deck(tan_delta=tan_delta, pattern=True))
        conductivity = run(root, "conductivity",
                           deck(sigma=equivalent_sigma, pattern=True))
        tangent_execution = tangent["jobs"][0]["executions"][0]
        conductivity_execution = conductivity["jobs"][0]["executions"][0]
        close(impedance(tangent_execution).real,
              impedance(conductivity_execution).real)
        close(impedance(tangent_execution).imag,
              impedance(conductivity_execution).imag)
        for key in ("accepted_w", "air_radiated_w", "dielectric_absorbed_w",
                    "balance_defect_w", "balance_defect_relative"):
            close(tangent_execution["power"][key],
                  conductivity_execution["power"][key], relative=2e-12)

        sweep_tangent = run(root, "sweep-tangent",
                            deck(tan_delta=0.02, frequency=200,
                                 count=3, step=100))
        sweep_sigma = run(root, "sweep-sigma",
                          deck(sigma=equivalent_sigma, frequency=200,
                               count=3, step=100))
        acceptance_sweep = run(root, "acceptance-sweep",
                               deck(tan_delta=tan_delta, frequency=250,
                                    count=15, step=5))
        parallel_sweep = run(root, "parallel-sweep",
                             deck(tan_delta=0.02, frequency=200,
                                  count=3, step=100), jobs=3)
        tangent_losses = []
        sigma_losses = []
        for execution in sweep_tangent["jobs"][0]["executions"]:
            epsilon = execution["material_evaluations"][0][
                "complex_relative_permittivity"]
            tangent_losses.append(-epsilon["im"])
        for execution in sweep_sigma["jobs"][0]["executions"]:
            epsilon = execution["material_evaluations"][0][
                "complex_relative_permittivity"]
            sigma_losses.append(-epsilon["im"])
        for loss in tangent_losses:
            close(loss, 3.2*0.02)
        close(sigma_losses[0]/sigma_losses[1], 1.5)
        close(sigma_losses[1]/sigma_losses[2], 4/3)
        for serial, parallel in zip(
                sweep_tangent["jobs"][0]["executions"],
                parallel_sweep["jobs"][0]["executions"]):
            close(impedance(serial).real, impedance(parallel).real,
                  relative=0.0, absolute=0.0)
            close(impedance(serial).imag, impedance(parallel).imag,
                  relative=0.0, absolute=0.0)
            if serial["material_evaluations"] != parallel["material_evaluations"]:
                raise AssertionError("parallel material evaluation changed")

        lapack = run(root, "lapack", deck(tan_delta=tan_delta),
                     backend="lapack", threads=4)
        close(impedance(tangent_execution).real,
              impedance(lapack["jobs"][0]["executions"][0]).real,
              relative=1e-8, absolute=1e-10)
        close(impedance(tangent_execution).imag,
              impedance(lapack["jobs"][0]["executions"][0]).imag,
              relative=1e-8, absolute=1e-10)

        for result in (tangent, conductivity, sweep_tangent, sweep_sigma,
                       acceptance_sweep):
            for execution in result["jobs"][0]["executions"]:
                power = execution["power"]
                if power["dielectric_absorbed_w"] < 0:
                    raise AssertionError("passive coating reported negative loss")
        for execution in acceptance_sweep["jobs"][0]["executions"]:
            power = execution["power"]
            if power["balance_defect_relative"] > 0.003:
                raise AssertionError(power)

        zero_execution = zero["jobs"][0]["executions"][0]
        if zero_execution["power"]["dielectric_absorbed_w"] is not None:
            raise AssertionError("zero loss changed the legacy power record")

        acceptance = json.loads((ROOT / "tests/coated-wire/m4-acceptance.json").read_text(
            encoding="utf-8"))
        comparison = acceptance["nec2c_openems_comparison"]
        if abs(comparison["frequency_error_percent"]) > 1.0:
            raise AssertionError(comparison)
        if abs(comparison["resistance_error_ohm"]) > 2.0:
            raise AssertionError(comparison)
        if abs(comparison["efficiency_error_percentage_points"]) > 1.0:
            raise AssertionError(comparison)

        power = tangent_execution["power"]
        print("M4 tangent/conductivity equivalence and per-frequency evaluation passed")
        print(f"M4 lossy impedance: {impedance(tangent_execution).real:.6f} "
              f"{impedance(tangent_execution).imag:+.6f}j ohm")
        print(f"M4 dielectric absorption: {power['dielectric_absorbed_w']:.9f} W; "
              f"power defect {100*power['balance_defect_relative']:.4f}%")
        print(f"M4 OpenEMS: resonance error {comparison['frequency_error_percent']:.4f}%, "
              f"resistance error {comparison['resistance_error_ohm']:.4f} ohm, "
              f"efficiency error {comparison['efficiency_error_percentage_points']:.4f} point")
        print("M4 passivity, backend/process equivalence, and zero-loss compatibility passed")


if __name__ == "__main__":
    main()
