#!/usr/bin/env python3
"""M9 Gate D1: lock first-order saddle scale, phase, and limitations."""

import json
import os
from pathlib import Path
import sys


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material  # noqa: E402
from nec4_m9_asymptotic import (first_order_cross,  # noqa: E402
                                second_order_cross_reference)


def relative(actual, expected):
    return abs(actual - expected) / max(abs(expected), 1e-300)


def evaluate(filename, evaluator):
    report = json.loads((ROOT / "tests/nec4-similarity" / filename).read_text())
    frequency = report["frequency_hz"]
    lower = material(frequency, report["lower_medium"]["epsilon_r"],
                     conductivity=report["lower_medium"]["conductivity_s_per_m"])
    upper = material(frequency, 1)
    record = report["heldout_samples"][-1]
    approximate = evaluator(
        frequency, complex(lower.wave_number), complex(upper.wave_number),
        record["rho_m"], record["z_m"], record["source_depth_m"])
    reference = {name: complex(*value)
                 for name, value in record["components"].items()}
    return {name: relative(approximate[name], reference[name])
            for name in approximate}


def main():
    dielectric = evaluate("m9-model-fit-er16.json", first_order_cross)
    ordinary = ("rho_v", "rho_h", "phi_h", "z_h")
    if max(dielectric[name] for name in ordinary) > 0.20:
        raise AssertionError(("first-order far-ray scale/phase", dielectric))
    if dielectric["z_v"] < 0.5:
        raise AssertionError("grazing vertical-field limitation is no longer exposed")
    second_dielectric = evaluate("m9-model-fit-er16.json",
                                 second_order_cross_reference)
    second_lossy = evaluate("m9-model-fit.json",
                            second_order_cross_reference)
    if max(second_dielectric.values()) > 0.013:
        raise AssertionError(("second-order dielectric far ray", second_dielectric))
    if max(second_lossy.values()) > 0.013:
        raise AssertionError(("second-order lossy far ray", second_lossy))

    grazing = json.loads((ROOT / "tests/nec4-similarity" /
                          "m9-asymptotic-gate.json").read_text())
    if grazing["classification"] != "documented-gn2-grazing-interference-limit":
        raise AssertionError("grazing limitation classification")
    if not all(max(case["second_order_uniform_error"].values()) > 1
               for case in grazing["cases"]):
        raise AssertionError("the published GN2 grazing limitation is no longer exposed")
    print("M9 Gate D1/D2 first- and second-order SDP branches passed")
    print(("second order is <=1.3% at 2.8 lambda; the independently sampled "
           "4-lambda grazing-interference limit remains explicit"))


if __name__ == "__main__":
    main()
