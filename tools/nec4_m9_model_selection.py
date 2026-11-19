#!/usr/bin/env python3
"""Audit equation (5-10) phase readings and numerical rank.

This is deliberately a characterization tool, not a production dispatcher.
It reuses persisted M8 samples so the expensive direct integrals are not
recomputed while model-selection rules are studied.
"""

import argparse
import cmath
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy import linalg


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_equations import material  # noqa: E402
from nec4_m9_tables import _beta, saddle_lambda1  # noqa: E402


COMPONENTS = ("rho_v", "z_v", "rho_h", "phi_h", "z_h")
PHASE_MODES = ("ordinary-ray", "literal-phase-root",
               "literal-exponential-root")


def decode(value):
    return complex(value[0], value[1])


def basis(point, k_lower, k_upper, component, phase_mode):
    rho, z, depth = point
    source_z = -depth
    r1 = math.hypot(rho, z + abs(k_upper / k_lower) * depth)
    r2 = math.hypot(rho, depth)
    lambda1 = saddle_lambda1(rho, z, depth, k_lower, k_upper)
    beta1 = _beta(k_lower, lambda1)
    beta2 = _beta(k_upper, lambda1)
    phase1 = cmath.exp(-1j * (beta2 * z - beta1 * source_z + lambda1 * rho))
    radial = rho / r1
    values = []
    for ell in range(1, 4):
        for m in range(4):
            for n in range(3):
                value = (r1 ** -ell * (z / r1) ** m
                         * (source_z / r1) ** n * phase1)
                if component in ("rho_v", "z_h") or (component == "z_v" and ell == 1):
                    value *= radial
                values.append(value)
    if r2 == 0:
        supplemental = [0j] * 8
    else:
        lambda2 = k_lower * rho / r2
        argument = _beta(k_upper, lambda2) * z + k_lower * r2
        if phase_mode == "ordinary-ray":
            phase2 = cmath.exp(-1j * argument)
        elif phase_mode == "literal-phase-root":
            phase2 = cmath.exp(-1j * cmath.sqrt(argument))
        elif phase_mode == "literal-exponential-root":
            phase2 = cmath.sqrt(cmath.exp(-1j * argument))
        else:
            raise ValueError(phase_mode)
        f1 = r1 ** -1 * radial * (source_z / r1) * phase2
        f2 = r1 ** -1 * radial * (source_z / r1) ** 2 * phase2
        f3 = r1 ** -1 * r2 ** -1 * radial * phase2
        f4 = f3 * (source_z / r1)
        supplemental = [f1, f2, f3, f4,
                        f1 * z / r1, f2 * z / r1,
                        f3 * z / r1, f4 * z / r1]
        if component in ("rho_v", "z_h"):
            supplemental = [value * radial for value in supplemental]
    return values + supplemental


def relative_rms(predicted, reference):
    return float(np.linalg.norm(predicted - reference) /
                 max(np.linalg.norm(reference), np.finfo(float).tiny))


def qr_fit(matrix, reference, calibration_matrix, calibration_reference, rank):
    norms = np.linalg.norm(matrix, axis=0)
    keep = norms > np.finfo(float).eps * max(1.0, norms.max())
    kept_indices = np.flatnonzero(keep)
    scaled = matrix[:, keep] / norms[keep]
    _, _, pivots = linalg.qr(scaled, mode="economic", pivoting=True)
    selected = pivots[:min(rank, len(pivots))]
    coefficients, _, _, _ = linalg.lstsq(
        scaled[:, selected], reference, lapack_driver="gelsy")
    physical = np.zeros(matrix.shape[1], dtype=complex)
    physical[kept_indices[selected]] = coefficients / norms[kept_indices[selected]]
    return (relative_rms(matrix @ physical, reference),
            relative_rms(calibration_matrix @ physical, calibration_reference))


def records(report, key):
    result = []
    for record in report[key]:
        point = (record["rho_m"], record["z_m"], record["source_depth_m"])
        values = {name: decode(value) for name, value in record["components"].items()}
        result.append((point, values))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=ROOT / "tests/nec4-similarity/m9-model-fit-er16.json")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tests/nec4-similarity/m9-model-selection.json")
    arguments = parser.parse_args()
    source = json.loads(arguments.input.read_text(encoding="utf-8"))
    frequency = source["frequency_hz"]
    lower = material(frequency, source["lower_medium"]["epsilon_r"],
                     conductivity=source["lower_medium"]["conductivity_s_per_m"])
    upper = material(frequency, 1)
    k_lower, k_upper = complex(lower.wave_number), complex(upper.wave_number)
    training = records(source, "documented_samples")
    calibration = records(source, "heldout_samples")
    result = {}
    for phase_mode in PHASE_MODES:
        phase_result = {}
        for component in COMPONENTS:
            matrix = np.asarray([basis(point, k_lower, k_upper, component, phase_mode)
                                 for point, _ in training])
            reference = np.asarray([values[component] for _, values in training])
            calibration_matrix = np.asarray([
                basis(point, k_lower, k_upper, component, phase_mode)
                for point, _ in calibration])
            calibration_reference = np.asarray([
                values[component] for _, values in calibration])
            candidates = []
            for rank in range(4, matrix.shape[1] + 1):
                train_error, calibration_error = qr_fit(
                    matrix, reference, calibration_matrix,
                    calibration_reference, rank)
                candidates.append({"rank": rank,
                                   "training_rms": train_error,
                                   "calibration_rms": calibration_error})
            selected = min(candidates,
                           key=lambda item: item["training_rms"] +
                           item["calibration_rms"])
            phase_result[component] = {
                "selected": selected,
                "full_rank": candidates[-1],
            }
        result[phase_mode] = phase_result
    report = {
        "schema": "nec2c-nec4-m9-model-selection",
        "schema_version": 1,
        "source": str(arguments.input.relative_to(ROOT)),
        "frequency_hz": frequency,
        "lower_medium": source["lower_medium"],
        "training_rows": len(training),
        "calibration_points": len(calibration),
        "selection_metric": "minimum training_rms + calibration_rms",
        "phase_modes": result,
        "production_status": "characterization-only",
    }
    arguments.output.write_text(json.dumps(report, indent=2) + "\n",
                                encoding="utf-8")
    ordinary = result["ordinary-ray"]
    for component in COMPONENTS:
        selected = ordinary[component]["selected"]
        print(f"{component}: rank {selected['rank']}, "
              f"training {selected['training_rms']:.3%}, "
              f"calibration {selected['calibration_rms']:.3%}")
    print(arguments.output)


if __name__ == "__main__":
    main()
