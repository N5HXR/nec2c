#!/usr/bin/env python3
"""M9 Gate B1: transformed linear-table convergence against the M8 oracle."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, appendix_h_bessel_contour
from nec4_m8_equations import CrossGeometry, material
from nec4_m9_tables import restore_value, storage_factor, table_scales, trilinear


FREQUENCY = mp.mpf("7.1e6")


def finite_remainder(lower, upper, rho, z, depth):
    result = appendix_h_bessel_contour(
        FREQUENCY, lower, upper, CrossGeometry(rho, -depth, z),
        precision_digits=50, cutoff_decay=60, panel_subdivisions=2,
        quadrature_order=36)
    return result.finite_remainder


def main():
    mp.mp.dps = 65
    lower = material(FREQUENCY, mp.mpf("13"),
                     conductivity=mp.mpf("0.005"))
    upper = material(FREQUENCY, mp.mpf("1"))
    wavelength0 = mp.mpf("299792458") / FREQUENCY
    scales = table_scales(complex(lower.wave_number),
                          complex(upper.wave_number), float(wavelength0))
    center = tuple(mp.mpf("0.12") * scales.lg for _ in range(3))
    direct = finite_remainder(lower, upper, *center)
    field_scale = max(abs(value) for value in direct.values())
    levels = []

    for half_fraction in ("0.06", "0.03", "0.015", "0.0075", "0.00375",
                          "0.001875", "0.0009375"):
        half = mp.mpf(half_fraction) * scales.lg
        bounds = tuple((float(coordinate - half), float(coordinate + half))
                       for coordinate in center)
        component_corners = {name: {} for name in COMPONENTS}
        for ir in (0, 1):
            for iz in (0, 1):
                for idepth in (0, 1):
                    rho = mp.mpf(bounds[0][ir])
                    z = mp.mpf(bounds[1][iz])
                    depth = mp.mpf(bounds[2][idepth])
                    values = finite_remainder(lower, upper, rho, z, depth)
                    factor = storage_factor(
                        1, float(rho), float(z), -float(depth),
                        complex(lower.wave_number), complex(upper.wave_number))
                    for name in COMPONENTS:
                        component_corners[name][(ir, iz, idepth)] = (
                            complex(values[name]) * factor)

        errors = {}
        for name in COMPONENTS:
            stored = trilinear(bounds, component_corners[name],
                               tuple(float(value) for value in center))
            restored = restore_value(
                stored, 1, float(center[0]), float(center[1]),
                -float(center[2]), complex(lower.wave_number),
                complex(upper.wave_number))
            absolute = abs(restored - complex(direct[name]))
            denominator = max(abs(direct[name]), field_scale * mp.mpf("1e-11"))
            errors[name] = float(absolute / denominator)
        levels.append({
            "half_width_over_Lg": float(mp.mpf(half_fraction)),
            "maximum_null_aware_relative_error": max(errors.values()),
            "component_errors": errors,
        })

    maxima = [level["maximum_null_aware_relative_error"] for level in levels]
    if any(newer >= older for older, newer in zip(maxima, maxima[1:])):
        raise AssertionError(("non-monotonic table convergence", maxima))
    ratios = [older / newer for older, newer in zip(maxima, maxima[1:])]
    if min(ratios) < 2.5:
        raise AssertionError(("weaker than expected linear-cell convergence", ratios))
    if maxima[-1] > 1e-4:
        raise AssertionError(("refined cell misses Gate B1 budget", maxima[-1]))

    report = {
        "schema": "nec2c-nec4-m9-linear-table-convergence",
        "schema_version": 1,
        "frequency_hz": float(FREQUENCY),
        "lower_medium": {"epsilon_r": 13.0, "conductivity_s_per_m": 0.005},
        "Lg_m": scales.lg,
        "center_over_Lg": [0.12, 0.12, 0.12],
        "levels": levels,
        "refinement_ratios": ratios,
    }
    output = ROOT / "tests" / "nec4-similarity" / "m9-linear-gate.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M9 Gate B1 transformed trilinear table converges against the M8 finite remainder")
    print("errors " + ", ".join(f"{value:.3e}" for value in maxima))


if __name__ == "__main__":
    main()
