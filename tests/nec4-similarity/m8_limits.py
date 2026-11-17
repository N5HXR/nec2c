#!/usr/bin/env python3
"""M8 Gate C2: PEC approach and far-ray Fresnel limit."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, transformed_pairing  # noqa: E402
from nec4_m8_equations import (  # noqa: E402
    assemble_cartesian, homogeneous_components, material,
)


def true_relative(left, right):
    return abs(left - right) / max(abs(left), abs(right), mp.mpf("1e-80"))


def matrix_relative(left, right):
    numerator = max(abs(left[row, column] - right[row, column])
                    for row in range(3) for column in range(3))
    denominator = max(abs(right[row, column])
                      for row in range(3) for column in range(3))
    return numerator / max(denominator, mp.mpf("1e-80"))


def main():
    mp.mp.dps = 70
    frequency = mp.mpf("7.1e6")
    vacuum = material(frequency, 1)
    rho = mp.mpf(1)
    source_z = mp.mpf("0.7")
    observer_z = mp.mpf("0.4")
    direct = homogeneous_components(
        frequency, vacuum, rho, source_z, observer_z)
    image = homogeneous_components(
        frequency, vacuum, rho, -source_z, observer_z)
    pec = {name: direct[name] + image[name]
           if name in ("rho_v", "z_v") else direct[name] - image[name]
           for name in COMPONENTS}
    pec_errors = {}
    for conductivity in (mp.mpf(1), mp.mpf(100), mp.mpf("1e4")):
        lower = material(frequency, 1, conductivity=conductivity)
        exact = transformed_pairing(
            frequency, lower, vacuum, rho, source_z, observer_z,
            method="bessel", precision_digits=55, cutoff_decay=70,
            panel_subdivisions=2, quadrature_order=44)
        pec_errors[mp.nstr(conductivity, 8)] = max(
            true_relative(exact[name], pec[name]) for name in COMPONENTS)
    ordered_pec = list(pec_errors.values())
    if not all(right < left for left, right in zip(ordered_pec, ordered_pec[1:])):
        raise AssertionError(("PEC limit not monotone", pec_errors))
    if ordered_pec[-1] > mp.mpf("1e-4"):
        raise AssertionError(("PEC limit too far from image solution", ordered_pec[-1]))

    # Equation (5-35): apply TM (Rv) in the plane of incidence and TE (Rh)
    # normal to it.  The perfect-ground image orientation is already included
    # before multiplying by the finite-ground Fresnel coefficients.
    dielectric = material(frequency, 4)
    fresnel_errors = {}
    for distance_scale in (20, 50, 100, 200):
        scale = mp.mpf(distance_scale)
        radial = mp.mpf("1.2") * scale
        source = mp.mpf("0.4") * scale
        observer = mp.mpf("0.7") * scale
        exact_components = transformed_pairing(
            frequency, dielectric, vacuum, radial, source, observer,
            method="bessel", precision_digits=55, cutoff_decay=75,
            panel_subdivisions=2, quadrature_order=48)
        exact = assemble_cartesian(exact_components, mp.mpf(0))
        direct_dyadic = assemble_cartesian(homogeneous_components(
            frequency, vacuum, radial, source, observer), mp.mpf(0))
        image_dyadic = assemble_cartesian(homogeneous_components(
            frequency, vacuum, radial, -source, observer), mp.mpf(0))
        # PEC image: horizontal source columns reverse; vertical does not.
        image_dyadic[:, 0] = -image_dyadic[:, 0]
        image_dyadic[:, 1] = -image_dyadic[:, 1]
        image_distance = mp.sqrt(radial ** 2 + (observer + source) ** 2)
        sine = radial / image_distance
        cosine = (observer + source) / image_distance
        n1 = dielectric.wave_number
        n2 = vacuum.wave_number
        transmitted_cosine = mp.sqrt(1 - (n2 / n1 * sine) ** 2)
        rv = ((n1 * cosine - n2 * transmitted_cosine) /
              (n1 * cosine + n2 * transmitted_cosine))
        rh = ((n1 * transmitted_cosine - n2 * cosine) /
              (n1 * transmitted_cosine + n2 * cosine))
        reflected = mp.diag([rv, rh, rv]) * image_dyadic
        approximation = direct_dyadic + reflected
        fresnel_errors[str(distance_scale)] = matrix_relative(
            exact, approximation)
    ordered_fresnel = list(fresnel_errors.values())
    if not all(right < left for left, right in
               zip(ordered_fresnel, ordered_fresnel[1:])):
        raise AssertionError(("Fresnel residual not monotone", fresnel_errors))
    if ordered_fresnel[-1] > mp.mpf("0.02"):
        raise AssertionError(("far-ray residual", ordered_fresnel[-1]))

    evidence = {
        "schema": "nec2c-nec4-m8-limits-gate",
        "schema_version": 1,
        "pec_conductivity_s_per_m_to_error": {
            name: mp.nstr(value, 16) for name, value in pec_errors.items()},
        "fresnel_distance_scale_to_error": {
            name: mp.nstr(value, 16) for name, value in fresnel_errors.items()},
        "pec_final_error": mp.nstr(ordered_pec[-1], 16),
        "fresnel_final_error": mp.nstr(ordered_fresnel[-1], 16),
    }
    output = ROOT / "tests" / "nec4-similarity" / "m8-limits-gate.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M8 Gate C2 limits passed: PEC "
          f"{mp.nstr(ordered_pec[-1], 6)}, far-ray Fresnel "
          f"{mp.nstr(ordered_fresnel[-1], 6)}")


if __name__ == "__main__":
    main()
