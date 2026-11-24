#!/usr/bin/env python3
"""M9 Gate B3: qualify linear tables across materials and depth boundaries."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, appendix_h_bessel_contour, canonical_cross
from nec4_m8_equations import (CrossGeometry,
                               appendix_h_residual_singular_fields,
                               appendix_h_v0_fields, material)
from nec4_m9_tables import restore_value, storage_factor, table_scales, trilinear


FREQUENCY = mp.mpf("7.1e6")
MATERIALS = (
    ("lossless-er4", mp.mpf("4"), mp.mpf("0")),
    ("average-soil", mp.mpf("13"), mp.mpf("0.005")),
    ("seawater", mp.mpf("81"), mp.mpf("4")),
)


def reference(lower, upper, rho, z, depth, finite=False):
    geometry = CrossGeometry(rho, -depth, z)
    if finite:
        result = appendix_h_bessel_contour(
            FREQUENCY, lower, upper, geometry, precision_digits=44,
            cutoff_decay=50, panel_subdivisions=1, quadrature_order=28)
        return result.finite_remainder, result.components
    values = canonical_cross(
        FREQUENCY, lower, upper, geometry, method="bessel",
        precision_digits=42, cutoff_decay=48, panel_subdivisions=1,
        quadrature_order=28)
    return values, values


def add_table1_singular(lower, upper, rho, z, depth, finite):
    geometry = CrossGeometry(rho, -depth, z)
    v0 = appendix_h_v0_fields(FREQUENCY, lower, upper, geometry)
    residual = appendix_h_residual_singular_fields(
        FREQUENCY, lower, upper, geometry)
    return {name: finite[name] + v0[name] + residual[name]
            for name in COMPONENTS}


def interpolate(region, lower, upper, bounds, point):
    corners = {name: {} for name in COMPONENTS}
    # At a boundary face, corners on the zero-weight depth plane cannot affect
    # the published trilinear expression. Avoid four expensive oracle calls,
    # but populate those dictionary entries to exercise the same interpolator.
    active_depths = (0, 1)
    if point[2] == bounds[2][0]:
        active_depths = (0,)
    elif point[2] == bounds[2][1]:
        active_depths = (1,)
    for ir in (0, 1):
        for iz in (0, 1):
            for idepth in (0, 1):
                if idepth not in active_depths:
                    for name in COMPONENTS:
                        corners[name][(ir, iz, idepth)] = 0j
                    continue
                rho = mp.mpf(bounds[0][ir])
                z = mp.mpf(bounds[1][iz])
                depth = mp.mpf(bounds[2][idepth])
                stored_values, _ = reference(
                    lower, upper, rho, z, depth, finite=(region == 1))
                factor = storage_factor(
                    region, float(rho), float(z), -float(depth),
                    complex(lower.wave_number), complex(upper.wave_number))
                for name in COMPONENTS:
                    corners[name][(ir, iz, idepth)] = (
                        complex(stored_values[name]) * factor)
    restored = {}
    for name in COMPONENTS:
        stored = trilinear(bounds, corners[name], point)
        restored[name] = restore_value(
            stored, region, point[0], point[1], -point[2],
            complex(lower.wave_number), complex(upper.wave_number))
    if region == 1:
        restored = add_table1_singular(
            lower, upper, *(mp.mpf(value) for value in point), restored)
    return {name: complex(value) for name, value in restored.items()}


def errors(actual, expected):
    scale = max(abs(value) for value in expected.values())
    return {name: float(abs(actual[name] - complex(expected[name])) /
                        max(abs(expected[name]), scale * mp.mpf("1e-11")))
            for name in COMPONENTS}


def centered_case(region, lower, upper, scales):
    lg = mp.mpf(scales.lg)
    depth_factor = {1: mp.mpf("0.20"), 2: mp.mpf("0.70"),
                    3: mp.mpf("1.50")}[region]
    center = (mp.mpf("0.18") * lg, mp.mpf("0.22") * lg,
              depth_factor * lg)
    half = mp.mpf("0.0075") * lg
    bounds = tuple((float(value - half), float(value + half))
                   for value in center)
    point = tuple(float(value) for value in center)
    actual = interpolate(region, lower, upper, bounds, point)
    _, expected = reference(lower, upper, *center, finite=(region == 1))
    component_errors = errors(actual, expected)
    return {"center_over_Lg": [float(value / lg) for value in center],
            "half_width_over_Lg": float(half / lg),
            "maximum_null_aware_relative_error": max(component_errors.values()),
            "component_errors": component_errors}


def boundary_side(region, boundary, side, lower, upper, scales):
    lg = mp.mpf(scales.lg)
    half = mp.mpf("0.0075") * lg
    point_mp = (mp.mpf("0.18") * lg, mp.mpf("0.22") * lg, boundary)
    depth_bounds = ((boundary - 2 * half, boundary) if side == "below"
                    else (boundary, boundary + 2 * half))
    bounds = ((float(point_mp[0] - half), float(point_mp[0] + half)),
              (float(point_mp[1] - half), float(point_mp[1] + half)),
              tuple(float(value) for value in depth_bounds))
    point = tuple(float(value) for value in point_mp)
    return interpolate(region, lower, upper, bounds, point)


def boundary_case(left_region, right_region, boundary, lower, upper, scales):
    point = (mp.mpf("0.18") * scales.lg,
             mp.mpf("0.22") * scales.lg, boundary)
    _, expected = reference(lower, upper, *point,
                            finite=(left_region == 1))
    left = boundary_side(left_region, boundary, "below", lower, upper, scales)
    right = boundary_side(right_region, boundary, "above", lower, upper, scales)
    left_errors = errors(left, expected)
    right_errors = errors(right, expected)
    cross_errors = errors(left, right)
    return {"depth_over_Lg": float(boundary / scales.lg),
            "lower_numbered_owner": left_region,
            "left_maximum_error": max(left_errors.values()),
            "right_maximum_error": max(right_errors.values()),
            "cross_boundary_maximum_jump": max(cross_errors.values()),
            "left_component_errors": left_errors,
            "right_component_errors": right_errors,
            "cross_component_jumps": cross_errors}


def main():
    mp.mp.dps = 55
    upper = material(FREQUENCY, mp.mpf("1"))
    report = {"schema": "nec2c-nec4-m9-linear-qualification",
              "schema_version": 1, "frequency_hz": float(FREQUENCY),
              "materials": {}}
    maxima = []
    jumps = []
    for name, epsilon, conductivity in MATERIALS:
        lower = material(FREQUENCY, epsilon, conductivity=conductivity)
        wavelength0 = mp.mpf("299792458") / FREQUENCY
        scales = table_scales(complex(lower.wave_number),
                              complex(upper.wave_number), float(wavelength0))
        centered = {str(region): centered_case(
            region, lower, upper, scales) for region in (1, 2, 3)}
        boundaries = {
            "1-2": boundary_case(1, 2, mp.mpf("0.4") * scales.lg,
                                 lower, upper, scales),
            "2-3": boundary_case(2, 3, mp.mpf(scales.lg),
                                 lower, upper, scales),
        }
        report["materials"][name] = {
            "epsilon_r": float(epsilon),
            "conductivity_s_per_m": float(conductivity),
            "Lg_m": scales.lg, "Lmax_m": scales.lmax,
            "centered_cells": centered, "shared_boundaries": boundaries}
        maxima.extend(case["maximum_null_aware_relative_error"]
                      for case in centered.values())
        for case in boundaries.values():
            maxima.extend((case["left_maximum_error"],
                           case["right_maximum_error"]))
            jumps.append(case["cross_boundary_maximum_jump"])
    report["worst_interpolation_error"] = max(maxima)
    report["worst_shared_boundary_jump"] = max(jumps)
    if report["worst_interpolation_error"] > 5e-3:
        raise AssertionError(("Gate B3 interpolation envelope",
                              report["worst_interpolation_error"]))
    if report["worst_shared_boundary_jump"] > 8e-3:
        raise AssertionError(("Gate B3 shared-boundary jump",
                              report["worst_shared_boundary_jump"]))
    output = ROOT / "tests/nec4-similarity/m9-linear-qualification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M9 Gate B3 material and shared-depth-boundary qualification passed")
    print(f"worst interpolation {report['worst_interpolation_error']:.3e}; "
          f"boundary jump {report['worst_shared_boundary_jump']:.3e}")


if __name__ == "__main__":
    main()
