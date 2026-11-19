#!/usr/bin/env python3
"""Qualify the M9 original-NEC4 grazing warning against the M8 oracle."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material  # noqa: E402
from nec4_m9_asymptotic import (COMPONENTS,  # noqa: E402
                                first_order_cross,
                                first_order_uniform_cross,
                                second_order_cross_reference,
                                second_order_uniform_cross_reference)
from nec4_m9_dispatch import dispatch  # noqa: E402
from nec4_m9_model_fit import direct_worker  # noqa: E402


FREQUENCY = 7.1e6
MATERIALS = (
    ("lossless-er16", 16.0, 0.0),
    ("lossy-soil", 13.0, 0.005),
)
HEIGHTS = (0.02, 0.05, 0.10, 0.20, 0.40, 0.80)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tests/nec4-similarity" /
                        "m9-grazing-warning-gate.json")
    arguments = parser.parse_args()
    wavelength = 299792458.0 / FREQUENCY
    tasks = []
    metadata = []
    for name, epsilon_r, conductivity in MATERIALS:
        for height in HEIGHTS:
            point = (4.0 * wavelength, height * wavelength,
                     height * wavelength)
            tasks.append((point, FREQUENCY, epsilon_r, conductivity))
            metadata.append((name, epsilon_r, conductivity, height, point))
    with ProcessPoolExecutor(max_workers=arguments.jobs) as executor:
        references = list(executor.map(direct_worker, tasks))

    rows = []
    for metadata_row, reference in zip(metadata, references):
        name, epsilon_r, conductivity, height, point = metadata_row
        lower = material(FREQUENCY, epsilon_r, conductivity=conductivity)
        upper = material(FREQUENCY, 1)
        k_lower = complex(lower.wave_number)
        k_upper = complex(upper.wave_number)
        arguments_common = (FREQUENCY, k_lower, k_upper, *point)
        approximations = {
            "first_order": first_order_cross(*arguments_common),
            "first_order_uniform": first_order_uniform_cross(*arguments_common),
            "second_order": second_order_cross_reference(*arguments_common),
            "second_order_uniform": second_order_uniform_cross_reference(
                *arguments_common),
        }
        reference_norm = math.sqrt(sum(abs(value) ** 2 for value in reference))
        method_errors = {}
        for method, approximate in approximations.items():
            components = {
                component: float(abs(approximate[component] - reference[index]) /
                                 max(abs(reference[index]), 1e-300))
                for index, component in enumerate(COMPONENTS)
            }
            vector_error = math.sqrt(sum(
                abs(approximate[component] - reference[index]) ** 2
                for index, component in enumerate(COMPONENTS))) / max(
                    reference_norm, 1e-300)
            method_errors[method] = {
                "component_relative_error": components,
                "maximum_component_relative_error": max(components.values()),
                "vector_relative_error": float(vector_error),
            }
        decision = dispatch(point[0], point[1], -point[2],
                            k_lower, k_upper, wavelength)
        rows.append({
            "material": name,
            "epsilon_r": epsilon_r,
            "conductivity_s_per_m": conductivity,
            "rho_lambda0": 4.0,
            "observer_height_lambda0": height,
            "source_depth_lambda0": height,
            "direct_ray_angle_degrees": math.degrees(
                math.atan2(2 * height, 4.0)),
            "warning": decision.grazing_risk,
            "reference": {component: [float(reference[index].real),
                                      float(reference[index].imag)]
                          for index, component in enumerate(COMPONENTS)},
            "method_errors": method_errors,
        })
    report = {
        "schema": "nec2c-nec4-m9-grazing-warning-gate",
        "schema_version": 1,
        "frequency_hz": FREQUENCY,
        "warning_is_diagnostic_not_correction": True,
        "rows": rows,
    }
    arguments.output.write_text(json.dumps(report, indent=2) + "\n",
                                encoding="utf-8")
    for row in rows:
        print(f"{row['material']} h={row['observer_height_lambda0']:.2f} "
              f"angle={row['direct_ray_angle_degrees']:.2f}deg "
              f"warn={row['warning']} "
              f"uniform-vector="
              f"{row['method_errors']['second_order_uniform']['vector_relative_error']:.1%} "
              f"uniform-max="
              f"{row['method_errors']['second_order_uniform']['maximum_component_relative_error']:.1%}")
    print(arguments.output)


if __name__ == "__main__":
    main()
