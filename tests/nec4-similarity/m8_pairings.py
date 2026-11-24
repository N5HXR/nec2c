#!/usr/bin/env python3
"""M8 Gate C1: Appendix G transformations in the identical-media limit."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, transformed_pairing  # noqa: E402
from nec4_m8_equations import (  # noqa: E402
    appendix_f_21_field_densities, homogeneous_components, material,
    relative_error,
)


def main():
    mp.mp.dps = 70
    frequency = mp.mpf("7.1e6")
    media = {
        "vacuum": material(frequency, 1),
        "lossless-er4": material(frequency, 4),
        "lossy-er13": material(frequency, 13, conductivity=mp.mpf("0.005")),
    }
    samples = {
        "12": (mp.mpf("-0.7"), mp.mpf("0.4")),
        "21": (mp.mpf("0.4"), mp.mpf("-0.7")),
        "11": (mp.mpf("-0.7"), mp.mpf("-0.4")),
        "22": (mp.mpf("0.7"), mp.mpf("0.4")),
    }
    rho = mp.mpf("1")
    errors = {}
    for medium_name, medium in media.items():
        for pairing, (source_z, observer_z) in samples.items():
            transformed = transformed_pairing(
                frequency, medium, medium, rho, source_z, observer_z,
                method="bessel", precision_digits=60, cutoff_decay=75,
                panel_subdivisions=2, quadrature_order=48)
            analytic = homogeneous_components(
                frequency, medium, rho, source_z, observer_z)
            for name in COMPONENTS:
                errors[f"{medium_name}:{pairing}:{name}"] = relative_error(
                    transformed[name], analytic[name])

    # A non-identical interface is evaluated directly with Appendix F's 21
    # integrands.  This path does not call Appendix G and therefore checks the
    # signs in reciprocity equations (G-6)-(G-10) independently.
    lower = media["lossy-er13"]
    upper = media["vacuum"]
    source_z = mp.mpf("0.7")
    observer_z = mp.mpf("-0.4")
    scale = min(1 / rho, 1 / (source_z - observer_z))
    p0 = mp.mpc(scale, scale)
    cutoff = mp.re(p0) + mp.mpf(80) / (source_z - observer_z)
    panel_width = mp.pi * scale / 2
    nodes, weights = mp.gauss_quadrature(56, "legendre")
    appendix_f = {name: mp.mpc(0) for name in COMPONENTS}

    def integrate_segment(start, end):
        midpoint = (start + end) / 2
        half_width = (end - start) / 2
        for node, weight in zip(nodes, weights):
            spectral = midpoint + half_width * node
            density = appendix_f_21_field_densities(
                spectral, rho, source_z, observer_z,
                lower, upper, frequency)
            for name in COMPONENTS:
                appendix_f[name] += half_width * weight * density[name]

    integrate_segment(mp.mpc(0), p0)
    left = mp.re(p0)
    while left < cutoff:
        right = min(left + panel_width, cutoff)
        integrate_segment(mp.mpc(left, mp.im(p0)),
                          mp.mpc(right, mp.im(p0)))
        left = right
    transformed_21 = transformed_pairing(
        frequency, lower, upper, rho, source_z, observer_z,
        method="bessel", precision_digits=65, cutoff_decay=80,
        panel_subdivisions=2, quadrature_order=56)
    appendix_f_errors = {
        name: relative_error(appendix_f[name], transformed_21[name])
        for name in COMPONENTS}

    # Source-free interface conditions are evaluated as one-sided limits.
    # Linear Richardson extrapolation removes the O(delta) displacement error
    # that would otherwise dominate a comparison of two finite-height points.
    delta = mp.mpf("5e-5")
    limit_options = dict(method="bessel", precision_digits=60,
                         cutoff_decay=80, panel_subdivisions=2,
                         quadrature_order=52)
    below_far = transformed_pairing(
        frequency, lower, upper, rho, mp.mpf("-0.7"), -delta,
        **limit_options)
    below_near = transformed_pairing(
        frequency, lower, upper, rho, mp.mpf("-0.7"), -delta / 2,
        **limit_options)
    above_far = transformed_pairing(
        frequency, lower, upper, rho, mp.mpf("-0.7"), delta,
        **limit_options)
    above_near = transformed_pairing(
        frequency, lower, upper, rho, mp.mpf("-0.7"), delta / 2,
        **limit_options)
    below_limit = {name: 2 * below_near[name] - below_far[name]
                   for name in COMPONENTS}
    above_limit = {name: 2 * above_near[name] - above_far[name]
                   for name in COMPONENTS}

    def true_relative(left, right):
        return abs(left - right) / max(abs(left), abs(right), mp.mpf("1e-80"))

    boundary_errors = {
        # Tangential E for vertical, radial-horizontal, and phi-horizontal
        # source orientations.
        "tangential:rho_v": true_relative(
            below_limit["rho_v"], above_limit["rho_v"]),
        "tangential:rho_h": true_relative(
            below_limit["rho_h"], above_limit["rho_h"]),
        "tangential:phi_h": true_relative(
            below_limit["phi_h"], above_limit["phi_h"]),
        # Normal D=epsilon*E for vertical and radial-horizontal sources.
        "normal-d:z_v": true_relative(
            lower.epsilon_relative * below_limit["z_v"],
            upper.epsilon_relative * above_limit["z_v"]),
        "normal-d:z_h": true_relative(
            lower.epsilon_relative * below_limit["z_h"],
            upper.epsilon_relative * above_limit["z_h"]),
    }

    worst_identical = max(errors.values())
    worst_appendix_f = max(appendix_f_errors.values())
    worst_boundary = max(boundary_errors.values())
    worst = max(worst_identical, worst_appendix_f, worst_boundary)
    target = mp.mpf("1e-9")
    if worst > target:
        raise AssertionError(("Appendix G identical-media limit", mp.nstr(worst, 12)))
    evidence = {
        "schema": "nec2c-nec4-m8-pairing-gate", "schema_version": 1,
        "equations": [f"G-{number}" for number in range(6, 21)],
        "target": mp.nstr(target, 8),
        "identical_media": list(media),
        "worst_identical_medium_error": mp.nstr(worst_identical, 16),
        "worst_appendix_f_21_error": mp.nstr(worst_appendix_f, 16),
        "worst_interface_condition_error": mp.nstr(worst_boundary, 16),
        "appendix_f_21_errors": {name: mp.nstr(value, 14)
                                  for name, value in appendix_f_errors.items()},
        "interface_condition_errors": {name: mp.nstr(value, 14)
                                       for name, value in boundary_errors.items()},
        "errors": {name: mp.nstr(value, 14) for name, value in errors.items()},
    }
    output = ROOT / "tests" / "nec4-similarity" / "m8-pairing-gate.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M8 Gate C1 Appendix F/G pairings passed: identical-medium "
          f"{mp.nstr(worst_identical, 6)}, direct-21 {mp.nstr(worst_appendix_f, 6)}, "
          f"interface {mp.nstr(worst_boundary, 6)}")


if __name__ == "__main__":
    main()
