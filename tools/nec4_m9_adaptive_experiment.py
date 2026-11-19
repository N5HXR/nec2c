#!/usr/bin/env python3
"""Screen interpolation candidates against persisted M9 direct-contour data.

This is an offline selection tool and never changes production routing.
Training points are the public-guide lattice; every reported error uses the
separately generated held-out direct-integral samples.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import RBFInterpolator


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material  # noqa: E402
from nec4_m9_tables import storage_factor, table_scales  # noqa: E402


COMPONENTS = ("rho_v", "z_v", "rho_h", "phi_h", "z_h")


def decode(record):
    point = np.asarray((record["rho_m"], record["z_m"],
                        record["source_depth_m"]), dtype=float)
    values = np.asarray([complex(*record["components"][name])
                         for name in COMPONENTS])
    return point, values


def common_factor(point, lower_k, upper_k, mode):
    rho, z, depth = point
    if mode == "raw":
        return 1 + 0j
    region = int(mode[1])
    return storage_factor(region, rho, z, -depth, lower_k, upper_k)


def real_matrix(values):
    return np.concatenate((values.real, values.imag), axis=1)


def complex_matrix(values):
    half = values.shape[1] // 2
    return values[:, :half] + 1j * values[:, half:]


def metrics(actual, expected):
    component_scale = np.max(np.abs(expected), axis=0)
    denominator = np.maximum(np.abs(expected), component_scale * 1e-11)
    component = np.abs(actual - expected) / np.maximum(denominator, 1e-300)
    vector = np.linalg.norm(actual - expected, axis=1) / np.maximum(
        np.linalg.norm(expected, axis=1), np.finfo(float).tiny)
    return {
        "maximum_null_aware_component_error": float(np.max(component)),
        "rms_vector_error": float(np.sqrt(np.mean(vector * vector))),
        "maximum_vector_error": float(np.max(vector)),
    }


def evaluate(source, transform, kernel, neighbors, smoothing, epsilon):
    frequency = source["frequency_hz"]
    medium = source["lower_medium"]
    lower = material(frequency, medium["epsilon_r"],
                     conductivity=medium["conductivity_s_per_m"])
    upper = material(frequency, 1)
    k1, k2 = complex(lower.wave_number), complex(upper.wave_number)
    scales = table_scales(k1, k2, 299792458.0 / frequency)
    training = [decode(row) for row in source["documented_samples"]]
    heldout = [decode(row) for row in source["heldout_samples"]]
    x_train = np.asarray([point for point, _ in training])
    x_test = np.asarray([point for point, _ in heldout])
    coordinate_scale = np.asarray((scales.wavelength0, scales.wavelength0,
                                   scales.lg))
    x_train = x_train / coordinate_scale
    x_test = x_test / coordinate_scale
    y_train = np.asarray([values * common_factor(point, k1, k2, transform)
                          for point, values in training])
    y_test = np.asarray([values for _, values in heldout])
    output_scale = np.maximum(np.linalg.norm(y_train, axis=0),
                              np.finfo(float).tiny)
    normalized = y_train / output_scale
    arguments = {"kernel": kernel, "smoothing": smoothing}
    if neighbors:
        arguments["neighbors"] = min(neighbors, len(training))
    if epsilon is not None:
        arguments["epsilon"] = epsilon
    interpolator = RBFInterpolator(x_train, real_matrix(normalized), **arguments)
    predicted = complex_matrix(interpolator(x_test)) * output_scale
    for index, (point, _) in enumerate(heldout):
        predicted[index] /= common_factor(point, k1, k2, transform)
    return metrics(predicted, y_test)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tests/nec4-similarity" /
                                "m9-adaptive-screen.json")
    args = parser.parse_args()
    paths = args.input or [ROOT / "tests/nec4-similarity" /
                           "m9-model-fit-er16.json",
                           ROOT / "tests/nec4-similarity" /
                           "m9-model-fit.json"]
    candidates = []
    for transform in ("raw", "f1", "f2", "f3"):
        for kernel in ("thin_plate_spline", "cubic", "quintic"):
            for neighbors in (24, 40, 64, None):
                for smoothing in (0.0, 1e-12, 1e-8, 1e-5):
                    candidates.append((transform, kernel, neighbors,
                                       smoothing, None))
        for kernel in ("multiquadric", "inverse_multiquadric",
                       "inverse_quadratic", "gaussian"):
            for epsilon in (0.3, 0.6, 1.0, 1.8, 3.0):
                for neighbors in (24, 40, 64, None):
                    candidates.append((transform, kernel, neighbors, 1e-12,
                                       epsilon))
    sources = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    reports = []
    for transform, kernel, neighbors, smoothing, epsilon in candidates:
        try:
            material_results = [evaluate(source, transform, kernel, neighbors,
                                         smoothing, epsilon)
                                for source in sources]
        except Exception:
            continue
        reports.append({
            "transform": transform, "kernel": kernel,
            "neighbors": neighbors, "smoothing": smoothing,
            "epsilon": epsilon,
            "score": max(item["rms_vector_error"] for item in material_results),
            "maximum_vector_error": max(item["maximum_vector_error"]
                                        for item in material_results),
            "materials": material_results,
        })
    reports.sort(key=lambda row: (row["score"], row["maximum_vector_error"]))
    report = {
        "schema": "nec2c-nec4-m9-adaptive-screen",
        "schema_version": 1,
        "training_policy": "documented samples only",
        "validation_policy": "persisted disjoint direct-contour heldout samples",
        "candidate_count": len(reports), "results": reports,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for row in reports[:12]:
        print(f"{row['transform']:3s} {row['kernel']:22s} "
              f"n={str(row['neighbors']):>4s} eps={str(row['epsilon']):>4s} "
              f"smooth={row['smoothing']:.0e}: "
              f"RMS {row['score']:.3%}, max {row['maximum_vector_error']:.3%}")
    print(args.output)


if __name__ == "__main__":
    main()
