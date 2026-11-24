#!/usr/bin/env python3
"""Qualify independently selected regular grids for the NEC4 M9 linear tables."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, canonical_cross  # noqa: E402
from nec4_m8_equations import (CrossGeometry,  # noqa: E402
    appendix_h_residual_singular_fields, appendix_h_v0_fields, material)
from nec4_m9_tables import (restore_value, storage_factor, table_scales,  # noqa: E402
                            trilinear)


FREQUENCY = mp.mpf("7.1e6")
MATERIALS = (
    ("lossless-er4", "4", "0"),
    ("average-soil", "13", "0.005"),
    ("seawater", "81", "4"),
)
MAX_STEPS_OVER_LG = {
    1: {"rho": 0.02, "observer_z": 0.02, "source_depth": 0.02},
    2: {"rho": 0.03, "observer_z": 0.03, "source_depth": 0.015},
    3: {"rho": 0.03, "observer_z": 0.03, "source_depth": 0.03},
}
PROBES = {
    1: ((.083, .077, .073), (.553, .047, .033), (.897, .103, .023),
        (.203, .797, .347), (.803, .797, .197)),
    2: ((.103, .097, .453), (.803, .047, .453), (.897, .797, .903),
        (.403, .297, .697)),
    # Negative depth values are fractions measured back from Lmax.
    3: ((.103, .097, 1.053), (.897, .047, 1.053), (.897, .803, 1.503),
        (.403, .297, -0.083)),
}


def axis(start, stop, maximum_step):
    intervals = max(1, math.ceil((stop - start) / maximum_step))
    return {"start": start, "stop": stop, "intervals": intervals,
            "count": intervals + 1, "step": (stop - start) / intervals}


def layout(region, lmax_over_lg):
    maximum = MAX_STEPS_OVER_LG[region]
    depth_ranges = {1: (0., .4), 2: (.4, 1.), 3: (1., lmax_over_lg)}
    depth_start, depth_stop = depth_ranges[region]
    return {"rho": axis(0., 1., maximum["rho"]),
            "observer_z": axis(0., 1., maximum["observer_z"]),
            "source_depth": axis(depth_start, depth_stop,
                                 maximum["source_depth"])}


def cell_bounds(axis_record, coordinate):
    position = (coordinate - axis_record["start"]) / axis_record["step"]
    index = min(axis_record["intervals"] - 1, max(0, math.floor(position)))
    low = axis_record["start"] + index * axis_record["step"]
    return low, low + axis_record["step"]


def physical_fields(lower, upper, rho, observer_z, source_depth):
    geometry = CrossGeometry(rho, -source_depth, observer_z)
    # Grid qualification asks whether the published transformations and
    # trilinear interpolation preserve the physical field.  Use the Bessel
    # representation everywhere so this sweep does not conflate that question
    # with Hankel-contour cost or stability.  M8 independently demonstrated
    # Bessel/Hankel agreement in both overlap and shallow geometries.
    total = canonical_cross(
        FREQUENCY, lower, upper, geometry, method="bessel",
        precision_digits=44, cutoff_decay=55, panel_subdivisions=1,
        quadrature_order=30)
    production_contour = (
        "bessel" if observer_z + source_depth > mp.mpf("0.5") * rho
        else "hankel")
    return total, production_contour


def stored_fields(region, lower, upper, rho, observer_z, source_depth):
    total, contour = physical_fields(lower, upper, rho, observer_z, source_depth)
    if region != 1:
        return total, contour
    geometry = CrossGeometry(rho, -source_depth, observer_z)
    v0 = appendix_h_v0_fields(FREQUENCY, lower, upper, geometry)
    residual = appendix_h_residual_singular_fields(
        FREQUENCY, lower, upper, geometry)
    finite = {name: total[name] - v0[name] - residual[name]
              for name in COMPONENTS}
    return finite, contour


def task(payload):
    mp.mp.dps = 55
    material_name, epsilon_text, conductivity_text, region, probe = payload
    lower = material(FREQUENCY, mp.mpf(epsilon_text),
                     conductivity=mp.mpf(conductivity_text))
    upper = material(FREQUENCY, mp.mpf("1"))
    wavelength0 = mp.mpf("299792458") / FREQUENCY
    scales = table_scales(complex(lower.wave_number),
                          complex(upper.wave_number), float(wavelength0))
    lmax_over_lg = scales.lmax / scales.lg
    grid = layout(region, lmax_over_lg)
    probe = list(probe)
    if region == 3 and probe[2] < 0:
        probe[2] = lmax_over_lg + probe[2] * (lmax_over_lg - 1.)
    # Keep the far-depth probe inside the final cell rather than on Lmax.
    probe[2] = min(probe[2], lmax_over_lg - .25 *
                   grid["source_depth"]["step"])
    bounds_normalized = tuple(cell_bounds(grid[name], coordinate)
                              for name, coordinate in zip(
                                  ("rho", "observer_z", "source_depth"), probe))
    bounds = tuple(tuple(value * scales.lg for value in pair)
                   for pair in bounds_normalized)
    point = tuple(value * scales.lg for value in probe)
    corners = {name: {} for name in COMPONENTS}
    contours = set()
    for ir in (0, 1):
        for iz in (0, 1):
            for idepth in (0, 1):
                rho = mp.mpf(bounds[0][ir])
                observer_z = mp.mpf(bounds[1][iz])
                source_depth = mp.mpf(bounds[2][idepth])
                values, contour = stored_fields(
                    region, lower, upper, rho, observer_z, source_depth)
                contours.add(contour)
                factor = storage_factor(
                    region, float(rho), float(observer_z), -float(source_depth),
                    complex(lower.wave_number), complex(upper.wave_number))
                for name in COMPONENTS:
                    corners[name][(ir, iz, idepth)] = complex(values[name]) * factor
    interpolated = {}
    for name in COMPONENTS:
        stored = trilinear(bounds, corners[name], point)
        interpolated[name] = restore_value(
            stored, region, point[0], point[1], -point[2],
            complex(lower.wave_number), complex(upper.wave_number))
    if region == 1:
        geometry = CrossGeometry(*(mp.mpf(value) for value in
                                   (point[0], -point[2], point[1])))
        v0 = appendix_h_v0_fields(FREQUENCY, lower, upper, geometry)
        residual = appendix_h_residual_singular_fields(
            FREQUENCY, lower, upper, geometry)
        interpolated = {name: interpolated[name] + complex(v0[name] + residual[name])
                        for name in COMPONENTS}
    expected, target_contour = physical_fields(
        lower, upper, *(mp.mpf(value) for value in point))
    scale = max(abs(value) for value in expected.values())
    errors = {name: float(abs(interpolated[name] - complex(expected[name])) /
                          max(abs(expected[name]), scale * mp.mpf("1e-10")))
              for name in COMPONENTS}
    return {"material": material_name, "region": region,
            "probe_over_Lg": probe, "cell_bounds_over_Lg": bounds_normalized,
            "contours": sorted(contours | {target_contour}),
            "component_errors": errors,
            "maximum_null_aware_relative_error": max(errors.values()),
            "layout": grid}


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--jobs", type=int, default=12)
    result.add_argument("--material", action="append",
                        choices=tuple(item[0] for item in MATERIALS),
                        help="qualify only this material (repeatable)")
    result.add_argument("--region", action="append", type=int,
                        choices=(1, 2, 3),
                        help="qualify only this table region (repeatable)")
    result.add_argument("--output", type=Path,
                        default=ROOT / "tests/nec4-similarity/m9-table-layout.json")
    return result


def main():
    arguments = parser().parse_args()
    selected_materials = tuple(
        item for item in MATERIALS
        if not arguments.material or item[0] in arguments.material)
    selected_regions = tuple(
        region for region in PROBES
        if not arguments.region or region in arguments.region)
    payloads = [(name, epsilon, conductivity, region, probe)
                for name, epsilon, conductivity in selected_materials
                for region in selected_regions for probe in PROBES[region]]
    if arguments.jobs < 1:
        raise SystemExit("--jobs must be positive")
    with ProcessPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = {executor.submit(task, payload): payload for payload in payloads}
        cases = []
        for completed, future in enumerate(as_completed(futures), 1):
            case = future.result()
            cases.append(case)
            print(f"[{completed}/{len(futures)}] {case['material']} "
                  f"region {case['region']}: "
                  f"{case['maximum_null_aware_relative_error']:.6e}",
                  flush=True)
    cases.sort(key=lambda case: (case["material"], case["region"],
                                 case["probe_over_Lg"]))
    worst = max(case["maximum_null_aware_relative_error"] for case in cases)
    summaries = {}
    for material_name, epsilon, conductivity in selected_materials:
        material_cases = [case for case in cases if case["material"] == material_name]
        layouts = {str(region): next(case["layout"] for case in material_cases
                                     if case["region"] == region)
                   for region in (1, 2, 3)}
        total_values = 0
        for grid in layouts.values():
            total_values += (grid["rho"]["count"] *
                             grid["observer_z"]["count"] *
                             grid["source_depth"]["count"] * len(COMPONENTS))
        summaries[material_name] = {
            "epsilon_r": float(epsilon),
            "conductivity_s_per_m": float(conductivity),
            "layouts": layouts, "complex_value_count": total_values,
            "payload_bytes": total_values * 16,
            "worst_error": max(case["maximum_null_aware_relative_error"]
                               for case in material_cases)}
    document = {"schema": "nec2c-nec4-m9-table-layout-qualification",
                "schema_version": 1, "frequency_hz": float(FREQUENCY),
                "reference_evaluator": "validated-bessel-representation",
                "contour_equivalence_evidence":
                "tests/nec4-similarity/m8-hankel-gate.json",
                "target_maximum_error": .005,
                "maximum_steps_over_Lg": MAX_STEPS_OVER_LG,
                "worst_error": worst, "materials": summaries, "cases": cases}
    arguments.output.write_text(json.dumps(document, indent=2) + "\n",
                                encoding="utf-8")
    print(arguments.output)
    print(f"worst null-aware interpolation error {worst:.6e}")
    for name, summary in summaries.items():
        print(f"{name}: {summary['worst_error']:.6e}, "
              f"{summary['payload_bytes'] / (1024 * 1024):.2f} MiB")
    if worst > document["target_maximum_error"]:
        raise SystemExit("candidate table layout misses the 0.5% target")


if __name__ == "__main__":
    main()
