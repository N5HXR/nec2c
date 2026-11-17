#!/usr/bin/env python3
"""M8 Gate B1: Figure 5-2 contour against exact finite-medium fields."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, bessel_contour  # noqa: E402
from nec4_m8_equations import (  # noqa: E402
    CrossGeometry,
    homogeneous_components,
    material,
    relative_error,
)


def compare(left, right):
    return {name: relative_error(left[name], right[name]) for name in COMPONENTS}


def encoded(values):
    return {name: {"real": mp.nstr(mp.re(value), 30),
                   "imag": mp.nstr(mp.im(value), 30)}
            for name, value in values.items()}


def main():
    mp.mp.dps = 80
    frequency = mp.mpf("7.1e6")
    media = [
        ("vacuum", material(frequency, 1)),
        ("lossless-er4", material(frequency, 4)),
        ("lossy-er13", material(frequency, 13, conductivity=mp.mpf("0.005"))),
    ]
    geometries = [
        CrossGeometry(mp.mpf("1"), mp.mpf("-0.3"), mp.mpf("0.31")),
        CrossGeometry(mp.mpf("1"), mp.mpf("-0.7"), mp.mpf("0.4")),
        CrossGeometry(mp.mpf("0.125"), mp.mpf("-0.8"), mp.mpf("0.45")),
    ]
    evidence = {"schema": "nec2c-nec4-m8-bessel-gate",
                "schema_version": 1, "samples": []}
    analytic_errors = []
    refinement_errors = []

    for medium_name, medium in media:
        for index, geometry in enumerate(geometries):
            coarse = bessel_contour(frequency, medium, medium, geometry,
                                     precision_digits=55, cutoff_decay=65,
                                     panel_subdivisions=1, quadrature_order=40)
            refined = bessel_contour(frequency, medium, medium, geometry,
                                      precision_digits=70, cutoff_decay=80,
                                      panel_subdivisions=2, quadrature_order=56)
            analytic = homogeneous_components(
                frequency, medium, geometry.rho, geometry.source_z,
                geometry.observer_z)
            sample_analytic = compare(refined.components, analytic)
            sample_refinement = compare(refined.components, coarse.components)
            analytic_errors.extend(sample_analytic.values())
            refinement_errors.extend(sample_refinement.values())
            evidence["samples"].append({
                "medium": medium_name,
                "geometry": {"rho": mp.nstr(geometry.rho, 20),
                             "source_z": mp.nstr(geometry.source_z, 20),
                             "observer_z": mp.nstr(geometry.observer_z, 20)},
                "contour_break": {"real": mp.nstr(mp.re(refined.contour_break), 20),
                                  "imag": mp.nstr(mp.im(refined.contour_break), 20)},
                "panels": refined.panel_count,
                "analytic_errors": {name: mp.nstr(value, 12)
                                    for name, value in sample_analytic.items()},
                "refinement_errors": {name: mp.nstr(value, 12)
                                      for name, value in sample_refinement.items()},
                "refined_components": encoded(refined.components),
            })

    # Analytic rho=0 limits avoid evaluating J0'/rho.  Cylindrical symmetry
    # requires the mixed radial/vertical terms to vanish and rho_h=-phi_h.
    axis_geometry = CrossGeometry(mp.mpf(0), mp.mpf("-0.7"), mp.mpf("0.4"))
    axis = bessel_contour(
        frequency, media[1][1], media[1][1], axis_geometry,
        precision_digits=70, cutoff_decay=80,
        panel_subdivisions=2, quadrature_order=56)
    axis_analytic = homogeneous_components(
        frequency, media[1][1], axis_geometry.rho,
        axis_geometry.source_z, axis_geometry.observer_z)
    axis_errors = compare(axis.components, axis_analytic)
    axis_errors["rho_v_null"] = abs(axis.components["rho_v"])
    axis_errors["z_h_null"] = abs(axis.components["z_h"])
    axis_errors["transverse_symmetry"] = relative_error(
        axis.components["rho_h"], -axis.components["phi_h"])
    analytic_errors.extend(axis_errors.values())
    evidence["axis"] = {
        "errors": {name: mp.nstr(value, 12)
                   for name, value in axis_errors.items()},
        "components": encoded(axis.components),
    }

    # A true interface has no elementary closed form here.  It is retained as
    # a refinement-only witness so the material mismatch and lossy branch are
    # exercised before the independent Hankel representation is added.
    interface_geometry = CrossGeometry(mp.mpf("1"), mp.mpf("-0.7"), mp.mpf("0.4"))
    lower = media[2][1]
    upper = media[0][1]
    coarse = bessel_contour(frequency, lower, upper, interface_geometry,
                             precision_digits=55, cutoff_decay=65,
                             panel_subdivisions=1, quadrature_order=40)
    refined = bessel_contour(frequency, lower, upper, interface_geometry,
                              precision_digits=70, cutoff_decay=80,
                              panel_subdivisions=2, quadrature_order=56)
    interface_errors = compare(refined.components, coarse.components)
    refinement_errors.extend(interface_errors.values())
    evidence["interface_witness"] = {
        "lower": "lossy-er13", "upper": "vacuum",
        "refinement_errors": {name: mp.nstr(value, 12)
                              for name, value in interface_errors.items()},
        "refined_components": encoded(refined.components),
    }

    # The documented rule must be an enforced domain boundary, not a comment.
    try:
        bessel_contour(frequency, upper, upper,
                       CrossGeometry(mp.mpf("1"), mp.mpf("-0.2"),
                                     mp.mpf("0.2")))
    except ValueError as error:
        if "0.5*rho" not in str(error):
            raise
        evidence["selection_rule_rejected"] = True
    else:
        raise AssertionError("Bessel contour accepted the Hankel region")

    worst_analytic = max(analytic_errors)
    worst_refinement = max(refinement_errors)
    target = mp.mpf("1e-10")
    if worst_analytic > target or worst_refinement > target:
        raise AssertionError(("M8 Bessel gate", mp.nstr(worst_analytic, 12),
                              mp.nstr(worst_refinement, 12)))
    evidence["target"] = mp.nstr(target, 8)
    evidence["worst_analytic_error"] = mp.nstr(worst_analytic, 16)
    evidence["worst_refinement_error"] = mp.nstr(worst_refinement, 16)
    output = ROOT / "tests" / "nec4-similarity" / "m8-bessel-gate.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M8 Gate B1 Bessel contour passed: worst analytic error "
          f"{mp.nstr(worst_analytic, 6)}, refinement {mp.nstr(worst_refinement, 6)}")


if __name__ == "__main__":
    main()
