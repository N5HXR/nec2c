#!/usr/bin/env python3
"""M7 validation for conducting-patch/dielectric coupling."""

import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


FLOAT = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")


def run(binary, root, name):
    deck = root / "tests" / "dielectric" / name
    with tempfile.NamedTemporaryFile(suffix=".out") as output:
        process = subprocess.run(
            [str(binary), "-i", str(deck), "-o", output.name],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        text = Path(output.name).read_text(encoding="utf-8")
    if process.returncode != 0:
        raise AssertionError(f"{name} exited {process.returncode}\n{text}\n{process.stderr}")
    return text


def coupled_patches(text):
    section = text.split("COUPLED SURFACE PATCH CURRENTS (A/M)", 1)[1]
    section = section.split("COUPLED VOLUME POLARIZATION", 1)[0]
    values = []
    for line in section.splitlines():
        numbers = [float(x) for x in FLOAT.findall(line)]
        if len(numbers) == 7:
            values.append(tuple(complex(numbers[i], numbers[i + 1]) for i in (1, 3, 5)))
    return values


def legacy_patches(text):
    section = text.split("SURFACE PATCH CURRENTS", 1)[1]
    values = []
    for line in section.splitlines():
        numbers = [float(x) for x in FLOAT.findall(line)]
        if len(numbers) == 14:
            values.append(tuple(complex(numbers[i], numbers[i + 1]) for i in (8, 10, 12)))
        if values and "DATA CARD" in line:
            break
    return values


def metric(text, label):
    match = re.search(re.escape(label) + r"\s*([-+0-9.Ee]+)", text)
    if not match:
        raise AssertionError(f"missing {label!r}")
    return float(match.group(1))


def feed_impedance(text, coupled):
    heading = "COUPLED ANTENNA INPUT PARAMETERS" if coupled else "ANTENNA INPUT PARAMETERS"
    section = text.split(heading, 1)[1]
    for line in section.splitlines():
        numbers = [float(x) for x in FLOAT.findall(line)]
        if len(numbers) >= (8 if coupled else 10):
            return complex(numbers[6], numbers[7])
    raise AssertionError(f"missing feed impedance after {heading}")


def norm(vector):
    return math.sqrt(sum(abs(component) ** 2 for row in vector for component in row))


def subtract(left, right):
    return [tuple(a - b for a, b in zip(lrow, rrow)) for lrow, rrow in zip(left, right)]


def main():
    binary = Path(os.environ.get("NEC2C", "./nec2c")).resolve()
    root = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
    legacy = run(binary, root, "m7-patch-legacy.nec")
    vacuum = run(binary, root, "m7-patch-vacuum.nec")
    coarse = run(binary, root, "m7-patch-cube-coarse.nec")
    er2 = run(binary, root, "m7-patch-cube-er2.nec")
    fine = run(binary, root, "m7-patch-cube.nec")
    far = run(binary, root, "m7-patch-cube-far.nec")
    mixed_legacy = run(binary, root, "m7-mixed-legacy.nec")
    mixed_vacuum = run(binary, root, "m7-mixed-vacuum.nec")
    mixed_loaded = run(binary, root, "m7-mixed-cube.nec")
    surface_vacuum = run(binary, root, "m7-patch4-vacuum.nec")
    surface_loaded = run(binary, root, "m7-patch4-cube.nec")

    legacy_j = legacy_patches(legacy)
    vacuum_j = coupled_patches(vacuum)
    if len(legacy_j) != 4 or len(vacuum_j) != 4:
        raise AssertionError("patch vacuum identity did not produce four currents")
    for expected, actual in zip(legacy_j, vacuum_j):
        for want, got in zip(expected, actual):
            if abs(got - want) > 5.0e-6:
                raise AssertionError(f"patch vacuum identity mismatch: {got} versus {want}")

    z_legacy = feed_impedance(mixed_legacy, False)
    z_vacuum = feed_impedance(mixed_vacuum, True)
    z_loaded = feed_impedance(mixed_loaded, True)
    if abs(z_legacy - z_vacuum) > 0.02:
        raise AssertionError(f"mixed vacuum feed mismatch: {z_vacuum} versus {z_legacy}")
    with (root / "tests" / "openems-m7-mixed-10mm.csv").open(newline="") as stream:
        openems = next(csv.DictReader(stream))
    nec_delta = z_loaded - z_vacuum
    openems_delta = complex(float(openems["delta_r_ohm"]), float(openems["delta_x_ohm"]))
    if abs(nec_delta.real - openems_delta.real) > 2.0 or abs(nec_delta.imag - openems_delta.imag) > 2.0:
        raise AssertionError(f"M7 openEMS perturbation mismatch: {nec_delta} versus {openems_delta}")

    coarse_delta = subtract(coupled_patches(coarse), vacuum_j)
    fine_delta = subtract(coupled_patches(fine), vacuum_j)
    far_delta = subtract(coupled_patches(far), vacuum_j)
    if norm(fine_delta) == 0.0:
        raise AssertionError("dielectric produced no patch-current perturbation")
    if norm(subtract(coarse_delta, fine_delta)) / norm(fine_delta) > 0.05:
        raise AssertionError("patch response did not converge under volume refinement")
    if not norm(far_delta) < norm(coarse_delta):
        raise AssertionError("patch dielectric effect did not decrease with separation")
    c2v_er2 = metric(er2, "PATCH-TO-VOLUME BLOCK NORM:")
    c2v_er3 = metric(coarse, "PATCH-TO-VOLUME BLOCK NORM:")
    v2c_er2 = metric(er2, "VOLUME-TO-PATCH BLOCK NORM:")
    v2c_er3 = metric(coarse, "VOLUME-TO-PATCH BLOCK NORM:")
    if not math.isclose(c2v_er3 / c2v_er2, 2.0, rel_tol=2.0e-4):
        raise AssertionError("patch-to-volume block did not scale with susceptibility")
    if not math.isclose(v2c_er3, v2c_er2, rel_tol=2.0e-4):
        raise AssertionError("volume-to-patch Green block changed with material")

    # Compare the dielectric perturbation of the area-integrated surface
    # current.  The 2x2 and 4x4 patches have areas 0.04 and 0.01 m^2.
    def integrated_delta(loaded, unloaded, area):
        delta = subtract(coupled_patches(loaded), coupled_patches(unloaded))
        return tuple(area * sum(row[axis] for row in delta) for axis in range(3))

    surface2 = integrated_delta(coarse, vacuum, 0.04)
    surface4 = integrated_delta(surface_loaded, surface_vacuum, 0.01)
    surface4_norm = math.sqrt(sum(abs(value) ** 2 for value in surface4))
    surface_difference = math.sqrt(sum(abs(a - b) ** 2 for a, b in zip(surface2, surface4)))
    if surface4_norm == 0.0 or surface_difference / surface4_norm > 0.15:
        raise AssertionError("patch response did not converge under surface refinement")

    for label in (
        "WIRE-TO-VOLUME BLOCK NORM:", "VOLUME-TO-WIRE BLOCK NORM:",
        "PATCH-TO-VOLUME BLOCK NORM:", "VOLUME-TO-PATCH BLOCK NORM:",
    ):
        if not metric(mixed_loaded, label) > 0.0:
            raise AssertionError(f"mixed model did not exercise {label}")
    for text in (vacuum, coarse, er2, fine, far, mixed_vacuum, mixed_loaded,
                 surface_vacuum, surface_loaded):
        if metric(text, "RELATIVE RESIDUAL:") > 1.0e-10:
            raise AssertionError("M7 coupled residual exceeded 1e-10")

    print("PASS: patch-only and mixed wire/patch vacuum identities")
    print("PASS: all four conductor/volume block directions are exercised")
    print("PASS: patch cross blocks have the required material scaling")
    print("PASS: patch response converges under surface and volume refinement")
    print("PASS: patch dielectric effect decreases with separation")
    print("PASS: M7 coupled residuals are at most 1e-10")
    print("PASS: mixed PEC/dielectric impedance perturbation agrees with openEMS within 2 ohms")
    return 0


if __name__ == "__main__":
    sys.exit(main())
