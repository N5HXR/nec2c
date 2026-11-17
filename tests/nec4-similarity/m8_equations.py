#!/usr/bin/env python3
"""M8 Gate A: convention and Sommerfeld-equation transcription tests."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import (  # noqa: E402
    CrossGeometry,
    assemble_cartesian,
    differentiated_field_densities,
    expanded_field_densities,
    material,
    outgoing_gamma,
    relative_error,
)


def complex_text(value):
    return {"real": mp.nstr(mp.re(value), 45),
            "imag": mp.nstr(mp.im(value), 45)}


def main():
    mp.mp.dps = 80
    frequency = mp.mpf("7.1e6")
    vacuum = material(frequency, 1)
    soil = material(frequency, 13, conductivity=mp.mpf("0.005"))

    # Convention and square-root sheet: at lambda=0, gamma=+j*k gives
    # exp(gamma*z' - gamma*z)=exp(-j*k*(z-z')) for z'<0<z.
    gamma0 = outgoing_gamma(mp.mpf(0), vacuum.wave_number)
    convention_error = relative_error(gamma0, 1j * vacuum.wave_number)
    geometry = CrossGeometry(mp.mpf("3.25"), mp.mpf("-0.7"), mp.mpf("0.4"))
    phase_error = relative_error(
        mp.exp(gamma0 * geometry.source_z - gamma0 * geometry.observer_z),
        mp.exp(-1j * vacuum.wave_number *
               (geometry.observer_z - geometry.source_z)))
    evanescent = outgoing_gamma(2 * vacuum.wave_number, vacuum.wave_number)
    if not (mp.re(evanescent) > 0 and abs(mp.im(evanescent)) < mp.mpf("1e-70")):
        raise AssertionError(("evanescent sheet", evanescent))

    # Bessel derivative identities are used only in the expanded path; the
    # differentiated path below uses mp.diff directly.
    bessel_errors = []
    for argument in (mp.mpf("0.125"), mp.mpc("2.75", "0.4"),
                     mp.mpc("11.0", "-0.75")):
        first = mp.diff(lambda x: mp.besselj(0, x), argument)
        second = mp.diff(lambda x: mp.besselj(0, x), argument, 2)
        bessel_errors += [relative_error(first, -mp.besselj(1, argument)),
                          relative_error(second, (mp.besselj(2, argument) -
                                                  mp.besselj(0, argument)) / 2)]

    # Compare the derivative form (5-2) to the independently expanded form
    # (5-7) at complex spectral points representative of the rising and
    # horizontal portions of the Figure 5-2 contour.
    witnesses = [
        mp.mpc("0.03125", "0.03125"),
        mp.mpc("0.2", "0.03125"),
        mp.mpc("0.75", "0.125"),
    ]
    identity_errors = {}
    witness_values = {}
    for index, spectral_point in enumerate(witnesses):
        expanded = expanded_field_densities(
            spectral_point, geometry, soil, vacuum, frequency)
        differentiated = differentiated_field_densities(
            spectral_point, geometry, soil, vacuum, frequency)
        witness_values[str(index)] = {
            "lambda": complex_text(spectral_point),
            "expanded": {name: complex_text(value)
                         for name, value in expanded.items()},
        }
        for name in expanded:
            identity_errors[f"{index}:{name}"] = relative_error(
                expanded[name], differentiated[name])

    # Equation (5-1): at phi=0 the radial/vertical plane is x-z and the
    # transverse horizontal response appears as -E_phi^H in yy.
    components = expanded_field_densities(
        witnesses[1], geometry, soil, vacuum, frequency)
    zero = assemble_cartesian(components, mp.mpf(0))
    assembly_errors = {
        "xx": relative_error(zero[0, 0], components["rho_h"]),
        "yy": relative_error(zero[1, 1], -components["phi_h"]),
        "xz": relative_error(zero[0, 2], components["rho_v"]),
        "zx": relative_error(zero[2, 0], components["z_h"]),
        "zz": relative_error(zero[2, 2], components["z_v"]),
        "xy": abs(zero[0, 1]),
        "yz": abs(zero[1, 2]),
    }

    worst = max([convention_error, phase_error] + bessel_errors +
                list(identity_errors.values()) + list(assembly_errors.values()))
    target = mp.mpf("1e-60")
    if worst > target:
        raise AssertionError(("M8 Gate A identity failure", mp.nstr(worst, 12)))

    evidence = {
        "schema": "nec2c-nec4-m8-equation-lock",
        "schema_version": 1,
        "convention": "exp(+j omega t), outgoing exp(-j k R), gamma(0)=+j k",
        "equations": ["5-1", "5-2a-e", "5-3", "5-4", "5-7a-e"],
        "precision_decimal_digits": mp.mp.dps,
        "target": mp.nstr(target, 8),
        "worst_relative_or_scaled_error": mp.nstr(worst, 16),
        "convention_error": mp.nstr(convention_error, 16),
        "phase_error": mp.nstr(phase_error, 16),
        "bessel_identity_errors": [mp.nstr(value, 16) for value in bessel_errors],
        "field_identity_errors": {name: mp.nstr(value, 16)
                                  for name, value in identity_errors.items()},
        "assembly_errors": {name: mp.nstr(value, 16)
                            for name, value in assembly_errors.items()},
        "witnesses": witness_values,
    }
    output = ROOT / "tests" / "nec4-similarity" / "m8-equation-lock.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M8 Gate A equation/convention lock passed; worst error "
          f"{mp.nstr(worst, 6)} at {mp.mp.dps} decimal digits")


if __name__ == "__main__":
    main()
