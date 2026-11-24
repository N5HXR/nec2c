#!/usr/bin/env python3
"""M9 Gate B2: table-2/table-3 transformed interpolation convergence."""

import json
import os
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, canonical_cross
from nec4_m8_equations import CrossGeometry, material
from nec4_m9_tables import restore_value, storage_factor, table_scales, trilinear


FREQUENCY = mp.mpf("7.1e6")


def field(lower, upper, rho, z, depth):
    return canonical_cross(
        FREQUENCY, lower, upper, CrossGeometry(rho, -depth, z),
        method="bessel", precision_digits=48, cutoff_decay=55,
        panel_subdivisions=1, quadrature_order=32)


def sweep_region(region, center, lower, upper, scales):
    direct = field(lower, upper, *center)
    field_scale = max(abs(value) for value in direct.values())
    levels = []
    for half_fraction in ("0.06", "0.03", "0.015", "0.0075"):
        half = mp.mpf(half_fraction) * scales.lg
        bounds = tuple((float(coordinate - half), float(coordinate + half))
                       for coordinate in center)
        component_corners = {name: {} for name in COMPONENTS}
        for ir in (0, 1):
            for iz in (0, 1):
                for idepth in (0, 1):
                    rho, z, depth = (mp.mpf(bounds[0][ir]),
                                     mp.mpf(bounds[1][iz]),
                                     mp.mpf(bounds[2][idepth]))
                    values = field(lower, upper, rho, z, depth)
                    factor = storage_factor(
                        region, float(rho), float(z), -float(depth),
                        complex(lower.wave_number), complex(upper.wave_number))
                    for name in COMPONENTS:
                        component_corners[name][(ir, iz, idepth)] = (
                            complex(values[name]) * factor)
        errors = {}
        for name in COMPONENTS:
            stored = trilinear(bounds, component_corners[name],
                               tuple(float(value) for value in center))
            restored = restore_value(
                stored, region, float(center[0]), float(center[1]),
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
    maxima = [entry["maximum_null_aware_relative_error"] for entry in levels]
    if any(newer >= older for older, newer in zip(maxima, maxima[1:])):
        raise AssertionError((region, "non-monotonic convergence", maxima))
    ratios = [older / newer for older, newer in zip(maxima, maxima[1:])]
    if min(ratios) < 2.5:
        raise AssertionError((region, "weak convergence", ratios))
    return {"center_over_Lg": [float(value / scales.lg) for value in center],
            "levels": levels, "refinement_ratios": ratios}


def main():
    mp.mp.dps = 60
    lower = material(FREQUENCY, mp.mpf("13"),
                     conductivity=mp.mpf("0.005"))
    upper = material(FREQUENCY, mp.mpf("1"))
    wavelength0 = mp.mpf("299792458") / FREQUENCY
    scales = table_scales(complex(lower.wave_number),
                          complex(upper.wave_number), float(wavelength0))
    lg = mp.mpf(scales.lg)
    cases = {
        "2": sweep_region(2, (mp.mpf("0.18") * lg,
                              mp.mpf("0.18") * lg,
                              mp.mpf("0.70") * lg), lower, upper, scales),
        "3": sweep_region(3, (mp.mpf("0.18") * lg,
                              mp.mpf("0.18") * lg,
                              mp.mpf("1.50") * lg), lower, upper, scales),
    }
    for region, case in cases.items():
        if case["levels"][-1]["maximum_null_aware_relative_error"] > 5e-3:
            raise AssertionError((region, "refined cell error", case["levels"][-1]))
    report = {"schema": "nec2c-nec4-m9-outer-linear-convergence",
              "schema_version": 1, "frequency_hz": float(FREQUENCY),
              "Lg_m": scales.lg, "regions": cases}
    output = ROOT / "tests" / "nec4-similarity" / "m9-linear-outer-gate.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M9 Gate B2 table-2/table-3 transformations converge against M8")
    for region, case in cases.items():
        values = [entry["maximum_null_aware_relative_error"]
                  for entry in case["levels"]]
        print(f"table {region}: " + ", ".join(f"{value:.3e}" for value in values))


if __name__ == "__main__":
    main()
