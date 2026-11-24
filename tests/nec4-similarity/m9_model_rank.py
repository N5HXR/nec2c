#!/usr/bin/env python3
"""M9 Gate C1: equation (5-10) structural rank and ambiguity audit."""

import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import qr


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material
from nec4_m9_tables import (MODEL_RAY_TERM_COUNT, MODEL_TOTAL_TERM_COUNT,
                            model_basis, model_sample_lattices, table_scales)


def design_matrix(region, scales, k_lower, k_upper, component=None):
    general, surface = model_sample_lattices(region, scales)
    points = general.points() + (() if surface is None else surface.points())
    rows = []
    labels = None
    for rho, z, depth in points:
        values, current_labels = model_basis(
            rho, z, depth, k_lower, k_upper,
            include_supplemental=region != 1, component=component)
        if labels is None:
            labels = current_labels
        elif labels != current_labels:
            raise AssertionError("basis labels changed between points")
        rows.append(values)
    return np.asarray(rows, dtype=np.complex128), labels


def audit(matrix, labels):
    norms = np.linalg.norm(matrix, axis=0)
    nonzero = norms > np.finfo(float).eps * max(1.0, norms.max())
    scaled = matrix[:, nonzero] / norms[nonzero]
    singular = np.linalg.svd(scaled, compute_uv=False)
    tolerance = max(scaled.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.count_nonzero(singular > tolerance))
    condition = float(singular[0] / singular[rank - 1]) if rank else float("inf")
    _, _, pivots = qr(scaled, mode="economic", pivoting=True)
    nonzero_indices = np.flatnonzero(nonzero)
    selected = [labels[int(nonzero_indices[index])] for index in pivots[:rank]]
    return {"rows": int(matrix.shape[0]), "columns": int(matrix.shape[1]),
            "nonzero_columns": int(np.count_nonzero(nonzero)), "rank": rank,
            "scaled_condition_at_rank": condition,
            "largest_singular_value": float(singular[0]),
            "smallest_retained_singular_value": float(singular[rank - 1]),
            "rank_revealing_columns": selected}


def main():
    frequency = 7.1e6
    lower = material(frequency, 13, conductivity=0.005)
    upper = material(frequency, 1)
    k_lower = complex(lower.wave_number)
    k_upper = complex(upper.wave_number)
    wavelength0 = 299792458.0 / frequency
    scales = table_scales(k_lower, k_upper, wavelength0)
    components = (None, "rho_v", "z_v", "rho_h", "phi_h", "z_h")
    result = {}
    for region in (1, 2, 3):
        result[str(region)] = {}
        for component in components:
            matrix, labels = design_matrix(region, scales, k_lower, k_upper,
                                           component=component)
            result[str(region)][component or "generic"] = audit(matrix, labels)

    generic1 = result["1"]["generic"]
    if generic1["rows"] != 32 or generic1["columns"] != MODEL_RAY_TERM_COUNT:
        raise AssertionError(generic1)
    if generic1["rank"] > generic1["rows"]:
        raise AssertionError("rank exceeds samples")
    for region, expected_rows in ((2, 76), (3, 88)):
        entry = result[str(region)]["generic"]
        if entry["rows"] != expected_rows or entry["columns"] != MODEL_TOTAL_TERM_COUNT:
            raise AssertionError(entry)
        if entry["rank"] < 30:
            raise AssertionError(("unexpected model collapse", region, entry))

    report = {
        "schema": "nec2c-nec4-m9-model-rank-audit",
        "schema_version": 1,
        "frequency_hz": frequency,
        "lower_medium": {"epsilon_r": 13.0, "conductivity_s_per_m": 0.005},
        "printed_equation_interpretation": {
            "ray_terms": MODEL_RAY_TERM_COUNT,
            "supplemental_terms": 8,
            "supplemental_phase": "exp(-j*(beta3*z+k1*R2))",
            "printed_trailing_half": "treated as typographical; literal form is dimensionally inconsistent",
        },
        "regions": result,
    }
    output = ROOT / "tests" / "nec4-similarity" / "m9-model-rank.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M9 Gate C1 equation (5-10) structural rank audit passed")
    for region in ("1", "2", "3"):
        entry = result[region]["generic"]
        print(f"region {region}: {entry['rows']}x{entry['columns']}, rank {entry['rank']}, "
              f"scaled condition {entry['scaled_condition_at_rank']:.3e}")


if __name__ == "__main__":
    main()
