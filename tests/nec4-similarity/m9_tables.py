#!/usr/bin/env python3
"""M9 Gate A: documented table coordinates and interpolation primitives."""

import json
import math
import os
from pathlib import Path
import sys


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m9_tables import (MODEL_RAY_TERM_COUNT,
                            MODEL_SUPPLEMENTAL_TERM_COUNT,
                            MODEL_TOTAL_TERM_COUNT,
                            model_sample_count, model_sample_lattices,
                            near_table_contains, restore_value,
                            source_depth_region, storage_factor, store_value,
                            table_scales, trilinear)


def close(actual, expected, tolerance=2e-14):
    scale = max(1.0, abs(expected))
    if abs(actual - expected) > tolerance * scale:
        raise AssertionError(f"{actual!r} != {expected!r}")


def main():
    scales = table_scales(4 + 3j, 1 + 0j, 2.0)
    close(scales.lg, 0.2 * 2.0)
    close(scales.lmax, 1.0 * 2.0)

    # Shared boundaries belong to the lower-numbered table, making lookup
    # deterministic without changing the guide's closed ranges.
    boundaries = [(0.0, 1), (-0.4 * scales.lg, 1),
                  (-scales.lg, 2), (-scales.lmax, 3),
                  (-math.nextafter(scales.lmax, math.inf), None)]
    for source_z, expected in boundaries:
        if source_depth_region(source_z, scales) != expected:
            raise AssertionError((source_z, source_depth_region(source_z, scales)))
    if not near_table_contains(scales.lg, scales.lg, -scales.lmax, scales):
        raise AssertionError("closed near-table extent")
    if near_table_contains(math.nextafter(scales.lg, math.inf), 0, 0, scales):
        raise AssertionError("rho outside near table")

    value = 1.25 - 0.75j
    transform_errors = []
    for region, geometry in ((1, (0.13, 0.08, -0.02)),
                             (2, (0.25, 0.15, -0.3)),
                             (3, (0.8, 0.2, -1.1))):
        arguments = (region, *geometry, 4 + 3j, 1 + 0j)
        stored = store_value(value, *arguments)
        restored = restore_value(stored, *arguments)
        close(restored, value)
        transform_errors.append(abs(restored - value))
        if storage_factor(*arguments) == 0:
            raise AssertionError("nonzero witness factor")

    bounds = ((0.2, 0.9), (0.1, 0.8), (0.3, 1.1))
    point = (0.37, 0.62, 0.91)

    def affine(rho, z, depth):
        return (2 - 3j) + (0.4 + 0.2j) * rho + (-0.7j) * z + 1.3 * depth

    corners = {}
    for ir in (0, 1):
        for iz in (0, 1):
            for idepth in (0, 1):
                corners[(ir, iz, idepth)] = affine(
                    bounds[0][ir], bounds[1][iz], bounds[2][idepth])
    interpolated = trilinear(bounds, corners, point)
    close(interpolated, affine(*point))

    counts = {str(region): model_sample_count(region, scales)
              for region in (1, 2, 3)}
    if counts != {"1": 32, "2": 76, "3": 88}:
        raise AssertionError(counts)
    if (MODEL_RAY_TERM_COUNT, MODEL_SUPPLEMENTAL_TERM_COUNT,
            MODEL_TOTAL_TERM_COUNT) != (36, 8, 44):
        raise AssertionError("equation 5-10 term count")
    lattice_summary = {}
    for region in (1, 2, 3):
        general, surface = model_sample_lattices(region, scales)
        points = general.points() + (() if surface is None else surface.points())
        lattice_summary[str(region)] = {
            "count": len(points),
            "unique_count": len(set(points)),
            "rho_range": [min(p[0] for p in points), max(p[0] for p in points)],
            "z_range": [min(p[1] for p in points), max(p[1] for p in points)],
            "depth_range": [min(p[2] for p in points), max(p[2] for p in points)],
        }

    report = {
        "schema": "nec2c-nec4-m9-table-equation-lock",
        "schema_version": 1,
        "medium_convention": "1=lower, 2=upper",
        "boundary_ownership": "lowest-numbered table",
        "scales": {"lambda0": scales.wavelength0, "Lg": scales.lg,
                   "Lmax": scales.lmax},
        "sample_counts": counts,
        "equation_5_10": {
            "ray_terms": MODEL_RAY_TERM_COUNT,
            "supplemental_terms": MODEL_SUPPLEMENTAL_TERM_COUNT,
            "total_terms": MODEL_TOTAL_TERM_COUNT,
            "region_1_documented_samples": counts["1"],
            "region_1_documented_supplemental_terms": 0,
            "rank_question": "36 printed basis terms but only 32 samples",
        },
        "sample_lattices": lattice_summary,
        "maximum_transform_roundtrip_error": max(transform_errors),
        "affine_interpolation_error": abs(interpolated - affine(*point)),
    }
    output = ROOT / "tests" / "nec4-similarity" / "m9-table-equation-lock.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M9 Gate A table regions, transforms, trilinear interpolation, and sample lattices passed")


if __name__ == "__main__":
    main()
