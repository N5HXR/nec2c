#!/usr/bin/env python3
"""M9 Gate E: dispatcher boundaries and grazing diagnostic policy."""

import json
import os
from pathlib import Path
import sys


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material  # noqa: E402
from nec4_m9_dispatch import (ASYMPTOTIC, LINEAR_TABLE, MODEL_FIT,  # noqa: E402
                              dispatch, model_subregion)
from nec4_m9_tables import table_scales  # noqa: E402


def main():
    frequency = 7.1e6
    wavelength = 299792458.0 / frequency
    lower = complex(material(frequency, 13, conductivity=0.005).wave_number)
    upper = complex(material(frequency, 1).wave_number)
    scales = table_scales(lower, upper, wavelength)

    def route(rho, z, depth):
        return dispatch(rho * wavelength, z * wavelength, -depth * wavelength,
                        lower, upper, wavelength)

    # Linear-table shared boundaries belong to the lowest numbered table.
    near_edge = dispatch(scales.lg, scales.lg, -0.4 * scales.lg,
                         lower, upper, wavelength)
    assert near_edge.method == LINEAR_TABLE
    just_outside = dispatch(scales.lg * (1 + 1e-12), scales.lg,
                            -0.4 * scales.lg, lower, upper, wavelength)
    assert just_outside.method == MODEL_FIT

    # Figure 5-4 polygon ownership: region 1 wins exact shared boundaries.
    region = lambda rho, z: model_subregion(
        rho * wavelength, z * wavelength, wavelength, scales.lg)
    assert region(0.1, 0.1) == 1
    assert region(0.1, 0.099) == 2
    assert region(0.5, scales.lg / wavelength) == 1
    assert region(0.5, scales.lg / wavelength - 1e-6) == 2
    assert region(0.6, 0.6) == 1
    assert region(0.600001, 0.600001) == 1
    assert region(0.600001, 0.6) == 3

    # Documented outer extents are inclusive; the next representable point is
    # asymptotic.  Source depth has the same deterministic ownership.
    assert route(3.0, 2.0, 0.8 * scales.lg / wavelength).method == MODEL_FIT
    assert route(3.0 + 1e-12, 0.2, 0.01).method == ASYMPTOTIC
    assert route(1.0, 2.0 + 1e-12, 0.01).method == ASYMPTOTIC
    assert route(1.0, 0.2, 0.8 * scales.lg / wavelength + 1e-12).method == ASYMPTOTIC

    grazing = route(4.0, 0.02, 0.02)
    assert grazing.method == ASYMPTOTIC and grazing.grazing_risk
    assert route(4.0, 0.2, 0.02).grazing_risk
    assert not route(4.0, 0.8, 0.8).grazing_risk
    assert not route(1.0, 0.02, 0.02).grazing_risk  # model path
    # Segment length is intentionally absent from the routing contract.

    report = {
        "schema": "nec2c-nec4-m9-dispatch-gate",
        "schema_version": 1,
        "frequency_hz": frequency,
        "lg_lambda0": scales.lg / wavelength,
        "lmax_lambda0": scales.lmax / wavelength,
        "outer_extent_lambda0": {"rho": 3.0, "z": 2.0,
                                  "source_depth_lg": 0.8},
        "grazing_warning_envelope": {
            "maximum_direct_ray_angle_degrees": 15.0,
            "method": ASYMPTOTIC,
            "witness_angle_degrees": grazing.grazing_angle_degrees,
        },
        "boundary_policy": "lowest-numbered/nearer method owns exact boundary",
    }
    output = ROOT / "tests/nec4-similarity/m9-dispatch-gate.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M9 dispatcher boundaries and grazing warning policy passed")
    print(output)


if __name__ == "__main__":
    main()
