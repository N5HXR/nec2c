#!/usr/bin/env python3
"""Independent M1 numerical references for the NEC4-similarity project.

This module intentionally imports no nec2c production code.  It implements
the public equations by direct quadrature and multiprecision arithmetic so
later production operators can be tested against an independent path.
"""

import argparse
import json
import math
import os
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")


PRODUCTION_RELATIVE_TARGET = 1.0e-11
REFERENCE_STABILITY_TARGET = 1.0e-13


def scaled_error(value, reference):
    return abs(value - reference) / max(abs(reference), 1.0)


def complex_pair(value):
    return [float(np.real(value)), float(np.imag(value))]


def basis_reference():
    """Conditioning of cos(x)-1 in the NEC4 short-segment basis."""
    rows = []
    stable_max = 0.0
    for x in np.logspace(-18, -1, 35):
        with mp.workdps(90):
            reference = mp.cos(mp.mpf(str(x))) - 1
        classic = math.cos(float(x)) - 1.0
        stable = -2.0 * math.sin(0.5 * float(x)) ** 2
        denominator = max(abs(float(reference)), np.finfo(float).tiny)
        classic_error = abs(classic - float(reference)) / denominator
        stable_error = abs(stable - float(reference)) / denominator
        stable_max = max(stable_max, stable_error)
        rows.append({"x": float(x), "classic_relative_error": classic_error,
                     "stable_relative_error": stable_error})
    return {"grid": rows, "stable_max_relative_error": stable_max}


def _quadrature_nodes(order, lower, upper):
    nodes, weights = leggauss(order)
    return (0.5 * (upper - lower) * nodes + 0.5 * (upper + lower),
            0.5 * (upper - lower) * weights)


def _composite_nodes(order, lower, upper, scales):
    """Gauss nodes split around the narrow near-singular region."""
    breaks = {float(lower), 0.0, float(upper)}
    limit = min(abs(lower), abs(upper))
    for scale in scales:
        for multiplier in (0.25, 1.0, 4.0):
            value = multiplier * scale
            if 0.0 < value < limit:
                breaks.add(value)
                breaks.add(-value)
    ordered = sorted(breaks)
    pieces = [_quadrature_nodes(order, a, b)
              for a, b in zip(ordered[:-1], ordered[1:]) if b > a]
    return (np.concatenate([piece[0] for piece in pieces]),
            np.concatenate([piece[1] for piece in pieces]))


def kernel_triplet(length_over_radius, rho_over_radius, order):
    """Equations (1)-(4) of the public NEC thin-wire accuracy report.

    Radius is the unit length.  Evaluation points avoid rho/a == 1, where
    specialized singular extraction is required and belongs to M4.
    """
    radius = 1.0
    length = float(length_over_radius)
    rho = float(rho_over_radius)
    wave_number = 0.05
    separation = max(abs(rho - radius), 1.0e-6)
    z, wz = _composite_nodes(order, -0.5 * length, 0.5 * length,
                             [separation, radius])
    phi, wp = _composite_nodes(
        order, -math.pi, math.pi,
        [separation / max(math.sqrt(rho), 1.0), 1.0])
    zz = z[:, None]
    pp = phi[None, :]
    distance = np.sqrt(rho * rho + radius * radius + zz * zz -
                       2.0 * radius * rho * np.cos(pp))
    exact = np.sum(wz[:, None] * wp[None, :] *
                   np.exp(-1j * wave_number * distance) / distance)
    axis_distance = np.sqrt(rho * rho + radius * radius + z * z)
    reduced = 2.0 * math.pi * np.sum(
        wz * np.exp(-1j * wave_number * axis_distance) / axis_distance)
    static_surface = np.sum(wz[:, None] * wp[None, :] / distance)
    extended = (2.0 * math.pi * np.sum(
        wz * (np.exp(-1j * wave_number * axis_distance) - 1.0) /
        axis_distance) + static_surface)
    return exact, extended, reduced


def kernel_reference():
    lengths = [2.0, 4.0, 8.0, 16.0, 32.0]
    radii = [0.0, 0.5, 0.9, 1.1, 2.0, 4.0]
    extended_error = []
    reduced_error = []
    convergence = []
    rows = []
    for length in lengths:
        extended_row = []
        reduced_row = []
        convergence_row = []
        for rho in radii:
            coarse = kernel_triplet(length, rho, 48)
            fine = kernel_triplet(length, rho, 80)
            reference, extended, reduced = fine
            conv = scaled_error(coarse[0], reference)
            ext_error = scaled_error(extended, reference)
            red_error = scaled_error(reduced, reference)
            extended_row.append(ext_error)
            reduced_row.append(red_error)
            convergence_row.append(conv)
            rows.append({
                "length_over_radius": length, "rho_over_radius": rho,
                "exact": complex_pair(reference),
                "extended": complex_pair(extended),
                "reduced": complex_pair(reduced),
                "quadrature_scaled_delta": conv,
                "extended_scaled_error": ext_error,
                "reduced_scaled_error": red_error,
            })
        extended_error.append(extended_row)
        reduced_error.append(reduced_row)
        convergence.append(convergence_row)
    return {
        "length_over_radius": lengths, "rho_over_radius": radii,
        "extended_scaled_error": extended_error,
        "reduced_scaled_error": reduced_error,
        "quadrature_scaled_delta": convergence,
        "maximum_quadrature_scaled_delta": float(np.max(convergence)),
        "samples": rows,
    }


def junction_reference():
    """Direct evaluation of the triangular charge potential in Eq. 3-44."""
    cases = []
    maximum_delta = mp.mpf("0")
    for angle_degrees, source_radius, match_s in (
            (180, "0.01", "0.05"), (90, "0.02", "0.08"),
            (45, "0.005", "0.12")):
        values = []
        for precision in (60, 90):
            with mp.workdps(precision):
                angle = mp.radians(angle_degrees)
                radius = mp.mpf(source_radius)
                match = mp.mpf(match_s)
                delta = mp.mpf("0.2")

                def integrand(source_s):
                    shape = delta - source_s
                    distance = mp.sqrt(match * match + source_s * source_s -
                                       2 * match * source_s * mp.cos(angle) +
                                       radius * radius)
                    return shape * mp.e**(-mp.mpc(0, 1) * mp.mpf("0.3") *
                                             distance) / distance
                values.append(mp.quad(integrand, [0, delta]))
        delta_value = scaled_error(complex(values[0]), complex(values[1]))
        maximum_delta = max(maximum_delta, mp.mpf(delta_value))
        cases.append({"angle_degrees": angle_degrees,
                      "source_radius": float(source_radius),
                      "match_distance": float(match_s),
                      "potential": complex_pair(complex(values[1])),
                      "precision_scaled_delta": delta_value})
    return {"cases": cases,
            "maximum_precision_scaled_delta": float(maximum_delta)}


def cap_reference():
    """Uniform flat-cap potential, with direct radial integration."""
    cases = []
    maximum_delta = mp.mpf("0")
    for radius_text, height_text in (("1", "0.2"), ("1", "1.5"),
                                     ("0.25", "0.03")):
        values = []
        analytics = []
        for precision in (60, 90):
            with mp.workdps(precision):
                radius = mp.mpf(radius_text)
                height = mp.mpf(height_text)
                numerical = 2 * mp.pi * mp.quad(
                    lambda radial: radial / mp.sqrt(radial * radial +
                                                     height * height),
                    [0, radius])
                analytic = 2 * mp.pi * (mp.sqrt(radius * radius +
                                                 height * height) - height)
                values.append(numerical)
                analytics.append(analytic)
        delta_value = scaled_error(float(values[0]), float(values[1]))
        identity_error = scaled_error(float(values[1]), float(analytics[1]))
        maximum_delta = max(maximum_delta, mp.mpf(delta_value),
                            mp.mpf(identity_error))
        cases.append({"radius": float(radius_text),
                      "axial_distance": float(height_text),
                      "potential": float(values[1]),
                      "precision_scaled_delta": delta_value,
                      "analytic_scaled_error": identity_error})
    return {"cases": cases,
            "maximum_scaled_error": float(maximum_delta)}


def sommerfeld_reference():
    """Direct Sommerfeld spectral identity for a scalar lossy Green function."""
    cases = []
    maximum_error = mp.mpf("0")
    for rho_text, height_text, alpha_text in (
            ("0.7", "1.3", "0.8"), ("2.0", "0.5", "0.4"),
            ("0.05", "2.0", "1.1"), ("1.5", "1.5", "0.2")):
        with mp.workdps(70):
            rho = mp.mpf(rho_text)
            height = mp.mpf(height_text)
            alpha = mp.mpf(alpha_text)

            def integrand(spectral):
                vertical = mp.sqrt(spectral * spectral + alpha * alpha)
                return (spectral / vertical * mp.besselj(0, spectral * rho) *
                        mp.exp(-vertical * height))
            direct = mp.quad(integrand, [0, mp.inf])
            distance = mp.sqrt(rho * rho + height * height)
            closed = mp.exp(-alpha * distance) / distance
            error = abs(direct - closed) / max(abs(closed), mp.mpf(1))
        maximum_error = max(maximum_error, error)
        cases.append({"rho": float(rho), "height_sum": float(height),
                      "attenuation": float(alpha), "direct": float(direct),
                      "closed": float(closed), "scaled_error": float(error)})
    return {"cases": cases, "maximum_scaled_error": float(maximum_error)}


def make_figures(output, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kernel = results["kernel"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for axis, key, title in (
            (axes[0], "extended_scaled_error", "Extended vs exact"),
            (axes[1], "reduced_scaled_error", "Reduced vs exact")):
        values = np.maximum(np.asarray(kernel[key]), 1.0e-16)
        image = axis.imshow(np.log10(values), origin="lower", aspect="auto",
                            extent=[min(kernel["rho_over_radius"]),
                                    max(kernel["rho_over_radius"]),
                                    math.log2(min(kernel["length_over_radius"])),
                                    math.log2(max(kernel["length_over_radius"]))],
                            cmap="viridis")
        axis.set_title(title)
        axis.set_xlabel(r"observation radius $\rho/a$")
        axis.set_ylabel(r"$\log_2(\Delta/a)$")
        fig.colorbar(image, ax=axis, label="log10 scaled error")
    fig.savefig(output / "m1-kernel-error-map.png", dpi=160)
    plt.close(fig)

    cases = results["sommerfeld"]["cases"]
    fig, axis = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
    labels = [f"{c['rho']:g}, {c['height_sum']:g}, {c['attenuation']:g}"
              for c in cases]
    axis.bar(labels, [max(c["scaled_error"], 1.0e-35) for c in cases])
    axis.set_yscale("log")
    axis.set_ylabel("direct-integral scaled error")
    axis.set_xlabel(r"$(\rho, z+z', \alpha)$")
    axis.set_title("Sommerfeld spectral quadrature vs closed identity")
    axis.tick_params(axis="x", rotation=25)
    fig.savefig(output / "m1-sommerfeld-error-map.png", dpi=160)
    plt.close(fig)


def evaluate():
    results = {
        "schema": "nec2c-nec4-m1-oracles", "schema_version": 1,
        "production_relative_target": PRODUCTION_RELATIVE_TARGET,
        "reference_stability_target": REFERENCE_STABILITY_TARGET,
        "provenance": {
            "basis_junction_cap": "nec4-theory.pdf sections 3.1, 3.4, 4.4 and appendix A",
            "kernels": "reports/accuracy.pdf equations (1)-(4)",
            "sommerfeld": "scalar Sommerfeld spectral identity used to qualify direct quadrature",
        },
        "basis": basis_reference(),
        "kernel": kernel_reference(),
        "junction": junction_reference(),
        "cap": cap_reference(),
        "sommerfeld": sommerfeld_reference(),
    }
    results["checks"] = {
        "kernel_quadrature_stable":
            results["kernel"]["maximum_quadrature_scaled_delta"] <=
            REFERENCE_STABILITY_TARGET,
        "junction_multiprecision_stable":
            results["junction"]["maximum_precision_scaled_delta"] <=
            REFERENCE_STABILITY_TARGET,
        "cap_multiprecision_stable":
            results["cap"]["maximum_scaled_error"] <=
            REFERENCE_STABILITY_TARGET,
        "sommerfeld_direct_identity":
            results["sommerfeld"]["maximum_scaled_error"] <=
            REFERENCE_STABILITY_TARGET,
    }
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    results = evaluate()
    destination = arguments.output_dir / "m1-oracles.json"
    destination.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    make_figures(arguments.output_dir, results)
    if arguments.self_test and not all(results["checks"].values()):
        failed = [name for name, passed in results["checks"].items()
                  if not passed]
        raise SystemExit("oracle self-test failed: " + ", ".join(failed))
    print("M1 independent oracle checks passed")
    print(f"kernel quadrature delta: {results['kernel']['maximum_quadrature_scaled_delta']:.3e}")
    print(f"Sommerfeld identity error: {results['sommerfeld']['maximum_scaled_error']:.3e}")
    print(destination)


if __name__ == "__main__":
    main()
