#!/usr/bin/env python3
"""M8 Gate B3: Appendix H cancellation and singularity subtraction."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import (  # noqa: E402
    CrossGeometry,
    MU0,
    appendix_h_difference_terms,
    appendix_h_remainder_field_densities,
    appendix_h_residual_singular_fields,
    appendix_h_v0,
    appendix_h_v0_fields,
    appendix_h_vprime_densities,
    expanded_field_densities,
    material,
    outgoing_gamma,
    potential_densities,
    relative_error,
)
from nec4_m8_direct import appendix_h_bessel_contour, bessel_contour  # noqa: E402


COMPONENTS = ("rho_v", "z_v", "rho_h", "phi_h", "z_h")


def raw_vprime(lambda_value, geometry, lower, upper):
    """Literal V12-V0 density from the first equation in Appendix H."""
    gamma1 = outgoing_gamma(lambda_value, lower.wave_number)
    gamma2 = outgoing_gamma(lambda_value, upper.wave_number)
    j0_measure = mp.besselj(0, lambda_value * geometry.rho) * lambda_value
    exponential = mp.exp(-gamma2 *
                         (geometry.observer_z - geometry.source_z))
    v12 = potential_densities(lambda_value, geometry, lower, upper)[0]
    v0 = (2 * exponential /
          (gamma2 * (lower.wave_number ** 2 + upper.wave_number ** 2)) *
          j0_measure)
    return v12 - v0


def differentiated_v0_fields(frequency, lower, upper, geometry):
    """Independent H-1 derivative construction from the scalar V0."""
    omega = 2 * mp.pi * frequency
    q = 1j * omega * MU0 / (4 * mp.pi)

    def value(rho, source_z, observer_z):
        return appendix_h_v0(
            frequency, lower, upper,
            CrossGeometry(rho, source_z, observer_z))

    rho = geometry.rho
    source_z = geometry.source_z
    observer_z = geometry.observer_z
    d_rho = mp.diff(lambda r: value(r, source_z, observer_z), rho)
    d_rho_rho = mp.diff(lambda r: value(r, source_z, observer_z), rho, 2)
    d_z_z = mp.diff(lambda z: value(rho, source_z, z), observer_z, 2)
    d_rho_z = mp.diff(
        lambda r: mp.diff(lambda z: value(r, source_z, z), observer_z), rho)
    d_rho_source = mp.diff(
        lambda r: mp.diff(lambda z: value(r, z, observer_z), source_z), rho)
    return {
        "rho_v": -q * d_rho_z,
        "z_v": -q * (d_z_z + upper.wave_number ** 2 *
                       value(rho, source_z, observer_z)),
        "rho_h": -q * d_rho_rho,
        "phi_h": q * d_rho / rho,
        "z_h": q * d_rho_source,
    }


def main():
    mp.mp.dps = 90
    frequency = mp.mpf("7.1e6")
    lower = material(frequency, 13, conductivity=mp.mpf("0.005"))
    upper = material(frequency, 1)
    geometry = CrossGeometry(mp.mpf("0.27"), mp.mpf("-0.018"),
                             mp.mpf("0.021"))

    # H-3 is exactly the literal V12-V0 integrand where direct subtraction is
    # still well conditioned.  H-4 is checked independently by differentiating
    # that literal expression with respect to source height.
    density_errors = {}
    for index, lambda_value in enumerate((mp.mpc("0.2", "0.03"),
                                          mp.mpc("1.7", "0.2"),
                                          mp.mpc("17", "0.5"))):
        stable, stable_derivative = appendix_h_vprime_densities(
            lambda_value, geometry, lower, upper)
        literal = raw_vprime(lambda_value, geometry, lower, upper)
        literal_derivative = mp.diff(
            lambda source_z: raw_vprime(
                lambda_value,
                CrossGeometry(geometry.rho, source_z, geometry.observer_z),
                lower, upper), geometry.source_z)
        density_errors[f"{index}:H-3"] = relative_error(stable, literal)
        density_errors[f"{index}:H-4"] = relative_error(
            stable_derivative, literal_derivative)

    # At very large lambda, verify the stable T terms against a much-higher
    # precision construction.  This exercises the cancellation regime that
    # motivated Appendix H, without treating a low-precision raw subtraction
    # as an oracle.
    lambda_large = mp.mpc("1e28", "2e27")
    with mp.workdps(180):
        gamma1_reference = outgoing_gamma(lambda_large, lower.wave_number)
        gamma2_reference = outgoing_gamma(lambda_large, upper.wave_number)
        reference_terms = appendix_h_difference_terms(
            gamma1_reference, gamma2_reference, geometry.source_z,
            lower.wave_number ** 2, upper.wave_number ** 2)
    with mp.workdps(65):
        gamma1 = outgoing_gamma(lambda_large, lower.wave_number)
        gamma2 = outgoing_gamma(lambda_large, upper.wave_number)
        stable_terms = appendix_h_difference_terms(
            gamma1, gamma2, geometry.source_z,
            lower.wave_number ** 2, upper.wave_number ** 2)
    cancellation_errors = [relative_error(value, reference)
                           for value, reference in zip(stable_terms,
                                                       reference_terms)]

    # H-2a-e are compared with a separate numerical differentiation of V0.
    analytic_fields = appendix_h_v0_fields(
        frequency, lower, upper, geometry)
    derivative_fields = differentiated_v0_fields(
        frequency, lower, upper, geometry)
    field_errors = {name: relative_error(analytic_fields[name],
                                         derivative_fields[name])
                    for name in COMPONENTS}

    # H-1 remainder densities plus the independently differentiated spectral
    # V0 contribution must reconstruct the original Section 5 density.
    remainder_density_errors = {}
    for index, lambda_value in enumerate((mp.mpc("0.4", "0.05"),
                                          mp.mpc("4", "0.4"))):
        full = expanded_field_densities(
            lambda_value, geometry, lower, upper, frequency)
        remainder = appendix_h_remainder_field_densities(
            lambda_value, geometry, lower, upper, frequency)

        def v0_density(rho, source_z, observer_z):
            gamma2 = outgoing_gamma(lambda_value, upper.wave_number)
            return (2 * mp.exp(-gamma2 * (observer_z - source_z)) /
                    (gamma2 * (lower.wave_number ** 2 +
                               upper.wave_number ** 2)) *
                    mp.besselj(0, lambda_value * rho) * lambda_value)

        q = 1j * 2 * mp.pi * frequency * MU0 / (4 * mp.pi)
        rho0, source0, observer0 = (geometry.rho, geometry.source_z,
                                    geometry.observer_z)
        v0_here = v0_density(rho0, source0, observer0)
        spectral_v0 = {
            "rho_v": -q * mp.diff(
                lambda r: mp.diff(lambda z: v0_density(r, source0, z),
                                  observer0), rho0),
            "z_v": -q * (mp.diff(
                lambda z: v0_density(rho0, source0, z), observer0, 2) +
                upper.wave_number ** 2 * v0_here),
            "rho_h": -q * mp.diff(
                lambda r: v0_density(r, source0, observer0), rho0, 2),
            "phi_h": q / rho0 * mp.diff(
                lambda r: v0_density(r, source0, observer0), rho0),
            "z_h": q * mp.diff(
                lambda r: mp.diff(lambda z: v0_density(r, z, observer0),
                                  source0), rho0),
        }
        for name in COMPONENTS:
            remainder_density_errors[f"{index}:{name}"] = relative_error(
                full[name], remainder[name] + spectral_v0[name])

    # Lock the second subtraction formulas to finite values over representative
    # approach directions.  Their cancellation against E' is tested by the
    # near-coincident direct-field qualification, not by self-comparison here.
    residual_witnesses = {}
    for label, point in {
            "radial": CrossGeometry(mp.mpf("2e-4"), mp.mpf("-1e-4"),
                                     mp.mpf("1e-4")),
            "oblique": CrossGeometry(mp.mpf("1.5e-4"), mp.mpf("-2e-4"),
                                      mp.mpf("0.5e-4")),
    }.items():
        values = appendix_h_residual_singular_fields(
            frequency, lower, upper, point)
        if not all(mp.isfinite(value) for value in values.values()):
            raise AssertionError(("non-finite Appendix H residual", label))
        residual_witnesses[label] = {
            name: {"real": mp.nstr(mp.re(value), 24),
                   "imag": mp.nstr(mp.im(value), 24)}
            for name, value in values.items()}

    # Exercise the subtraction in an actual quadrature, not only pointwise
    # algebra.  The full Section 5 integral and E'+S agree near coincidence.
    near = CrossGeometry(mp.mpf("1e-3"), mp.mpf("-5e-4"),
                         mp.mpf("5e-4"))
    subtracted = appendix_h_bessel_contour(
        frequency, lower, upper, near, precision_digits=65,
        cutoff_decay=80, panel_subdivisions=3, quadrature_order=52)
    unsubtracted = bessel_contour(
        frequency, lower, upper, near, precision_digits=65,
        cutoff_decay=80, panel_subdivisions=3, quadrature_order=52)
    integrated_errors = {
        name: relative_error(subtracted.components[name],
                             unsubtracted.components[name])
        for name in COMPONENTS}

    # The second subtraction must remove 1/R growth.  Along a fixed approach
    # angle every E'' component remains bounded as the geometry halves.
    finite_bounds = {}
    for scale_text in ("1e-3", "5e-4", "2.5e-4"):
        scale = mp.mpf(scale_text)
        result = appendix_h_bessel_contour(
            frequency, lower, upper,
            CrossGeometry(scale, -scale / 2, scale / 2),
            precision_digits=60, cutoff_decay=75,
            panel_subdivisions=2, quadrature_order=48)
        finite_bounds[scale_text] = max(
            abs(value) for value in result.finite_remainder.values())

    worst_density = max(density_errors.values())
    worst_cancellation = max(cancellation_errors)
    worst_field = max(field_errors.values())
    worst_remainder_density = max(remainder_density_errors.values())
    worst_integrated = max(integrated_errors.values())
    worst_finite = max(finite_bounds.values())
    identity_target = mp.mpf("1e-70")
    cancellation_target = mp.mpf("1e-55")
    if (worst_density > identity_target or worst_field > identity_target or
            worst_remainder_density > identity_target):
        raise AssertionError(("Appendix H transcription", worst_density,
                              worst_field, worst_remainder_density))
    if worst_cancellation > cancellation_target:
        raise AssertionError(("Appendix H cancellation", worst_cancellation))
    if worst_integrated > mp.mpf("1e-9"):
        raise AssertionError(("Appendix H integrated reconstruction",
                              worst_integrated))
    if worst_finite > mp.mpf("10"):
        raise AssertionError(("Appendix H E'' is not bounded", worst_finite))

    evidence = {
        "schema": "nec2c-nec4-m8-singularity-gate",
        "schema_version": 1,
        "equations": ["H-1a-e", "H-2a-e", "H-3", "H-4"],
        "identity_target": mp.nstr(identity_target, 8),
        "cancellation_target": mp.nstr(cancellation_target, 8),
        "worst_density_identity": mp.nstr(worst_density, 16),
        "worst_large_lambda_cancellation": mp.nstr(worst_cancellation, 16),
        "worst_v0_field_identity": mp.nstr(worst_field, 16),
        "worst_remainder_field_density_identity": mp.nstr(
            worst_remainder_density, 16),
        "worst_integrated_reconstruction": mp.nstr(worst_integrated, 16),
        "largest_finite_remainder": mp.nstr(worst_finite, 16),
        "density_errors": {name: mp.nstr(value, 16)
                           for name, value in density_errors.items()},
        "large_lambda_term_errors": [mp.nstr(value, 16)
                                     for value in cancellation_errors],
        "v0_field_errors": {name: mp.nstr(value, 16)
                            for name, value in field_errors.items()},
        "remainder_field_density_errors": {
            name: mp.nstr(value, 16)
            for name, value in remainder_density_errors.items()},
        "residual_singularity_witnesses": residual_witnesses,
        "integrated_reconstruction_errors": {
            name: mp.nstr(value, 16) for name, value in integrated_errors.items()},
        "finite_remainder_approach": {
            scale: mp.nstr(value, 16) for scale, value in finite_bounds.items()},
    }
    output = ROOT / "tests" / "nec4-similarity" / "m8-singularity-gate.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M8 Gate B3 Appendix H passed: density "
          f"{mp.nstr(worst_density, 5)}, large-lambda {mp.nstr(worst_cancellation, 5)}, "
          f"V0 field {mp.nstr(worst_field, 5)}, remainder "
          f"{mp.nstr(worst_remainder_density, 5)}, integrated "
          f"{mp.nstr(worst_integrated, 5)}")


if __name__ == "__main__":
    main()
