#!/usr/bin/env python3
"""M8 Gate B2: Figure 5-3 contour, sheet tracking, and overlap identity."""

import argparse
import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import (  # noqa: E402
    COMPONENTS, bessel_contour, hankel_contour,
)
from nec4_m8_equations import (  # noqa: E402
    CrossGeometry, expanded_field_densities, hankel_field_densities,
    material, relative_error,
)


def compare(left, right):
    return {name: relative_error(left[name], right[name]) for name in COMPONENTS}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true",
                        help="include the multi-minute shallow-region qualification")
    arguments = parser.parse_args()
    mp.mp.dps = 80
    frequency = mp.mpf("7.1e6")
    upper = material(frequency, 1)
    lower = material(frequency, 13, conductivity=mp.mpf("0.005"))
    overlap = CrossGeometry(mp.mpf("2"), mp.mpf("-0.7"), mp.mpf("0.4"))

    # On the original real contour the negative-lambda Hankel contribution
    # plus the positive one must exactly recover the 2*J Bessel density.
    pointwise_errors = []
    for lambda_value in (mp.mpf("0.05"), mp.mpf("0.12"), mp.mpf("0.2"),
                         mp.mpf("0.6"), mp.mpf("2.0")):
        bessel = expanded_field_densities(
            lambda_value, overlap, lower, upper, frequency)
        positive = hankel_field_densities(
            lambda_value, overlap, lower, upper, frequency)
        negative = hankel_field_densities(
            -lambda_value - 1j * mp.mpf("1e-50"), overlap, lower, upper,
            frequency)
        pointwise_errors.extend(compare(
            {name: positive[name] + negative[name] for name in COMPONENTS},
            bessel).values())

    refined = hankel_contour(frequency, lower, upper, overlap,
                              precision_digits=58, cutoff_decay=85,
                              panel_subdivisions=4, quadrature_order=56)
    bessel = bessel_contour(frequency, lower, upper, overlap,
                             precision_digits=70, cutoff_decay=80,
                             panel_subdivisions=2, quadrature_order=56)
    refinement_errors = {}
    coarse = None
    if arguments.full:
        coarse = hankel_contour(frequency, lower, upper, overlap,
                                 precision_digits=58, cutoff_decay=120,
                                 panel_subdivisions=8, quadrature_order=64)
        refined = hankel_contour(frequency, lower, upper, overlap,
                                  precision_digits=60, cutoff_decay=150,
                                  panel_subdivisions=10, quadrature_order=68)
        refinement_errors = compare(refined.components, coarse.components)
    representation_errors = compare(refined.components, bessel.components)

    shallow_errors = {}
    shallow_representation = {}
    if arguments.full:
        # The shallow point is intentionally expensive.  The final two Hankel
        # levels and the independently stable Bessel path must all agree.
        shallow = CrossGeometry(mp.mpf("4"), mp.mpf("-0.35"), mp.mpf("0.25"))
        shallow_coarse = hankel_contour(
            frequency, lower, upper, shallow, precision_digits=58,
            cutoff_decay=135, panel_subdivisions=12, quadrature_order=68)
        shallow_refined = hankel_contour(
            frequency, lower, upper, shallow, precision_digits=60,
            cutoff_decay=150, panel_subdivisions=14, quadrature_order=72)
        shallow_bessel = bessel_contour(
            frequency, lower, upper, shallow, precision_digits=70,
            cutoff_decay=85, panel_subdivisions=4, quadrature_order=56,
            permit_outside_region=True)
        shallow_errors = compare(shallow_refined.components,
                                 shallow_coarse.components)
        shallow_representation = compare(shallow_refined.components,
                                         shallow_bessel.components)

    worst_pointwise = max(pointwise_errors)
    worst_representation = max(list(representation_errors.values()) +
                               list(shallow_representation.values()))
    refinement_values = list(refinement_errors.values()) + list(shallow_errors.values())
    worst_refinement = max(refinement_values) if refinement_values else mp.mpf(0)
    target = mp.mpf("1e-9")
    if max(worst_pointwise, worst_representation, worst_refinement) > target:
        raise AssertionError(("M8 Hankel gate", mp.nstr(worst_pointwise, 12),
                              mp.nstr(worst_representation, 12),
                              mp.nstr(worst_refinement, 12)))

    evidence = {
        "schema": "nec2c-nec4-m8-hankel-gate", "schema_version": 1,
        "mode": "full" if arguments.full else "compact",
        "target": mp.nstr(target, 8),
        "worst_pointwise_identity": mp.nstr(worst_pointwise, 16),
        "worst_bessel_hankel_overlap": mp.nstr(worst_representation, 16),
        "worst_refinement": mp.nstr(worst_refinement, 16),
        "overlap_errors": {name: mp.nstr(value, 14)
                           for name, value in representation_errors.items()},
        "overlap_refinement_errors": {name: mp.nstr(value, 14)
                                      for name, value in refinement_errors.items()},
        "shallow_refinement_errors": {name: mp.nstr(value, 14)
                                      for name, value in shallow_errors.items()},
        "shallow_representation_errors": {name: mp.nstr(value, 14)
                                          for name, value in shallow_representation.items()},
        "root_tracking": {
            "coarse_principal_sign_changes": (
                coarse.principal_root_sign_changes if coarse else None),
            "refined_principal_sign_changes": refined.principal_root_sign_changes,
            "refined_panels": refined.panel_count,
            "refined_order": refined.quadrature_order,
        },
    }
    output_name = ("m8-hankel-gate.json" if arguments.full else
                   "m8-hankel-compact-gate.json")
    output = ROOT / "tests" / "nec4-similarity" / output_name
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M8 Gate B2 Hankel contour passed: pointwise "
          f"{mp.nstr(worst_pointwise, 5)}, overlap {mp.nstr(worst_representation, 5)}, "
          f"refinement {mp.nstr(worst_refinement, 5)}")


if __name__ == "__main__":
    main()
