#!/usr/bin/env python3
"""Generate M8 comparisons for the M9 grazing asymptotic paths."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material  # noqa: E402
from nec4_m9_asymptotic import (first_order_cross,  # noqa: E402
                                second_order_cross_reference,
                                second_order_uniform_cross_reference)
from nec4_m9_model_fit import direct_worker  # noqa: E402


FREQUENCY = 7.1e6
COMPONENTS = ("rho_v", "z_v", "rho_h", "phi_h", "z_h")


def pair(value):
    return [float(value.real), float(value.imag)]


def metrics(values, reference):
    return {name: float(abs(values[name] - reference[index]) /
                        max(abs(reference[index]), 1e-300))
            for index, name in enumerate(COMPONENTS)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tests/nec4-similarity" /
                        "m9-asymptotic-gate.json")
    arguments = parser.parse_args()
    wavelength = 299792458.0 / FREQUENCY
    cases = [
        {"name": "lossless-er16", "epsilon_r": 16.0, "conductivity": 0.0},
        {"name": "lossy-soil", "epsilon_r": 13.0, "conductivity": 0.005},
    ]
    point = (4 * wavelength, 0.02 * wavelength, 0.02 * wavelength)
    tasks = [(point, FREQUENCY, case["epsilon_r"], case["conductivity"])
             for case in cases]
    with ProcessPoolExecutor(max_workers=arguments.jobs) as executor:
        references = list(executor.map(direct_worker, tasks))
    results = []
    for case, reference in zip(cases, references):
        lower = material(FREQUENCY, case["epsilon_r"],
                         conductivity=case["conductivity"])
        upper = material(FREQUENCY, 1)
        arguments_common = (FREQUENCY, complex(lower.wave_number),
                            complex(upper.wave_number), *point)
        first = first_order_cross(*arguments_common)
        second = second_order_cross_reference(*arguments_common)
        uniform = second_order_uniform_cross_reference(*arguments_common)
        results.append({
            **case,
            "point_lambda0": {"rho": 4.0, "z": 0.02,
                              "source_depth": 0.02},
            "direct": {name: pair(reference[index])
                       for index, name in enumerate(COMPONENTS)},
            "first_order_error": metrics(first, reference),
            "second_order_error": metrics(second, reference),
            "second_order_uniform_error": metrics(uniform, reference),
        })
    report = {
        "schema": "nec2c-nec4-m9-asymptotic-gate",
        "schema_version": 1,
        "frequency_hz": FREQUENCY,
        "classification": "documented-gn2-grazing-interference-limit",
        "cases": results,
    }
    arguments.output.write_text(json.dumps(report, indent=2) + "\n",
                                encoding="utf-8")
    for result in results:
        print(result["name"], "second", max(result["second_order_error"].values()),
              "uniform", max(result["second_order_uniform_error"].values()))
    print(arguments.output)


if __name__ == "__main__":
    main()
