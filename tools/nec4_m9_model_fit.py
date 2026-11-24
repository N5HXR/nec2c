#!/usr/bin/env python3
"""Fit documented NEC4 equation (5-10) to M8 direct-field samples."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from scipy import linalg


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, canonical_cross  # noqa: E402
from nec4_m8_equations import CrossGeometry, material  # noqa: E402
from nec4_m9_tables import (model_basis, model_sample_lattices,  # noqa: E402
                            table_scales)


FREQUENCY = 7.1e6
LOWER_ER = 13.0
LOWER_CONDUCTIVITY = 0.005


def complex_pair(value):
    return [float(value.real), float(value.imag)]


def sample_record(point, values):
    return {"rho_m": point[0], "z_m": point[1], "source_depth_m": point[2],
            "components": {name: complex_pair(value)
                           for name, value in zip(COMPONENTS, values)}}


def direct_worker(task):
    point, frequency, lower_er, lower_conductivity = task
    rho, z, depth = point
    mp.mp.dps = 60
    lower = material(frequency, lower_er, conductivity=lower_conductivity)
    upper = material(frequency, 1)
    wavelength0 = mp.mpf("299792458") / frequency
    # The canonical oracle accepts either endpoint on the interface, but not
    # both.  Use a symmetric one-sided limit only for that single corner.
    if z == 0 and depth == 0:
        epsilon = wavelength0 * mp.mpf("1e-8")
        source_z, observer_z = -epsilon, epsilon
    else:
        source_z, observer_z = -mp.mpf(depth), mp.mpf(z)
    values = canonical_cross(
        frequency, lower, upper,
        CrossGeometry(mp.mpf(rho), source_z, observer_z),
        method="auto", precision_digits=55, cutoff_decay=80,
        panel_subdivisions=3, quadrature_order=48)
    return tuple(complex(values[name]) for name in COMPONENTS)


def rank_revealing_fit(matrix, values):
    norms = np.linalg.norm(matrix, axis=0)
    nonzero = norms > np.finfo(float).eps * max(1.0, norms.max())
    scaled = matrix[:, nonzero] / norms[nonzero]
    _, triangular, pivots = linalg.qr(scaled, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(triangular))
    tolerance = max(scaled.shape) * np.finfo(float).eps * diagonal[0]
    rank = int(np.count_nonzero(diagonal > tolerance))
    selected_local = pivots[:rank]
    selected_global = np.flatnonzero(nonzero)[selected_local]
    reduced = scaled[:, selected_local]
    reduced_coefficients, _, qr_rank, _ = linalg.lstsq(
        reduced, values, cond=tolerance, lapack_driver="gelsy")
    if qr_rank != rank:
        raise RuntimeError(("QR rank mismatch", rank, qr_rank))
    coefficients = np.zeros(matrix.shape[1], dtype=np.complex128)
    coefficients[selected_global] = reduced_coefficients / norms[selected_global]
    prediction = matrix @ coefficients
    return coefficients, selected_global, prediction, rank


def metrics(prediction, reference):
    scale = max(abs(value) for value in reference)
    errors = []
    for actual, expected in zip(prediction, reference):
        errors.append(abs(actual - expected) / max(abs(expected), scale * 1e-11))
    return {"maximum_null_aware_relative_error": float(max(errors)),
            "rms_relative_to_peak": float(
                np.linalg.norm(prediction - reference) /
                max(np.linalg.norm(reference), np.finfo(float).tiny))}


def unique_in_order(points):
    return list(dict.fromkeys(tuple(float(value) for value in point)
                              for point in points))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--frequency", type=float, default=FREQUENCY)
    parser.add_argument("--lower-er", type=float, default=LOWER_ER)
    parser.add_argument("--lower-conductivity", type=float,
                        default=LOWER_CONDUCTIVITY)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tests" / "nec4-similarity" /
                        "m9-model-fit.json")
    arguments = parser.parse_args()
    lower = material(arguments.frequency, arguments.lower_er,
                     conductivity=arguments.lower_conductivity)
    upper = material(arguments.frequency, 1)
    k_lower, k_upper = complex(lower.wave_number), complex(upper.wave_number)
    wavelength0 = 299792458.0 / arguments.frequency
    scales = table_scales(k_lower, k_upper, wavelength0)
    general, surface = model_sample_lattices(3, scales)
    documented_points = list(general.points()) + list(surface.points())
    unique_training = unique_in_order(documented_points)
    heldout_dimensionless = (
        (0.72, 0.10, 0.10), (0.90, 0.35, 0.20),
        (1.10, 0.65, 0.40), (1.35, 0.25, 0.65),
        (1.60, 0.90, 0.15), (1.90, 1.20, 0.50),
        (2.20, 0.40, 0.70), (2.50, 1.50, 0.30),
        (2.80, 1.80, 0.75),
    )
    heldout_points = [(rho * wavelength0, z * wavelength0, depth * scales.lg)
                      for rho, z, depth in heldout_dimensionless]
    all_unique = unique_training + heldout_points
    print(f"evaluating {len(unique_training)} unique documented samples and "
          f"{len(heldout_points)} held-out points with {arguments.jobs} workers",
          flush=True)
    tasks = [(point, arguments.frequency, arguments.lower_er,
              arguments.lower_conductivity) for point in all_unique]
    with ProcessPoolExecutor(max_workers=arguments.jobs) as executor:
        all_values = list(executor.map(direct_worker, tasks, chunksize=1))
    value_map = {point: values for point, values in zip(all_unique, all_values)}
    training_values = [value_map[tuple(float(value) for value in point)]
                       for point in documented_points]
    heldout_values = [value_map[point] for point in heldout_points]

    results = {}
    for component_index, component in enumerate(COMPONENTS):
        rows, labels = [], None
        for rho, z, depth in documented_points:
            values, current_labels = model_basis(
                rho, z, depth, k_lower, k_upper,
                include_supplemental=True, component=component)
            rows.append(values)
            labels = current_labels if labels is None else labels
        matrix = np.asarray(rows, dtype=np.complex128)
        references = np.asarray([values[component_index]
                                 for values in training_values])
        coefficients, selected, fitted, rank = rank_revealing_fit(matrix, references)

        heldout_rows = np.asarray([
            model_basis(rho, z, depth, k_lower, k_upper,
                        include_supplemental=True, component=component)[0]
            for rho, z, depth in heldout_points], dtype=np.complex128)
        heldout_reference = np.asarray([values[component_index]
                                        for values in heldout_values])
        heldout_prediction = heldout_rows @ coefficients
        results[component] = {
            "rank": rank,
            "selected_columns": [labels[index] for index in selected],
            "training": metrics(fitted, references),
            "heldout": metrics(heldout_prediction, heldout_reference),
            "coefficients": {labels[index]: complex_pair(coefficients[index])
                             for index in selected},
        }

    worst_training = max(value["training"]["maximum_null_aware_relative_error"]
                         for value in results.values())
    worst_heldout = max(value["heldout"]["maximum_null_aware_relative_error"]
                        for value in results.values())
    report = {
        "schema": "nec2c-nec4-m9-model-fit",
        "schema_version": 1,
        "region": 3,
        "frequency_hz": arguments.frequency,
        "lower_medium": {"epsilon_r": arguments.lower_er,
                         "conductivity_s_per_m": arguments.lower_conductivity},
        "documented_sample_rows": len(documented_points),
        "unique_documented_points": len(unique_training),
        "heldout_points": len(heldout_points),
        "interface_limit_offset_lambda0": 1e-8,
        "documented_samples": [
            sample_record(tuple(float(value) for value in point), values)
            for point, values in zip(documented_points, training_values)],
        "heldout_samples": [sample_record(point, values)
                            for point, values in zip(heldout_points,
                                                     heldout_values)],
        "components": results,
        "worst_training_error": worst_training,
        "worst_heldout_error": worst_heldout,
    }
    arguments.output.write_text(json.dumps(report, indent=2) + "\n",
                                encoding="utf-8")
    print(f"worst training error {worst_training:.6e}")
    print(f"worst held-out error {worst_heldout:.6e}")
    print(arguments.output)


if __name__ == "__main__":
    main()
