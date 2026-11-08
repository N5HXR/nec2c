#!/usr/bin/env python3
"""Analyze and plot the emerging M8 thickness/permittivity error surface."""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT/"tests/coated-wire/m8-openems-thickness.json"
OUTPUT = ROOT/"tests/coated-wire/m8-thickness-analysis.json"
FIGURE = ROOT/"tests/coated-wire/m8-thickness-error.png"


def main():
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source = json.loads(INPUT.read_text())
    cases = source["cases"]
    x = np.asarray([row["sqrt_er_k0_b"] for row in cases])
    y = np.asarray([row["error"]["frequency_percent"] for row in cases])
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope*x+intercept
    residual = y-predicted
    total = np.sum((y-y.mean())**2)
    r_squared = 1-np.sum(residual**2)/total
    crossing = (1.0-intercept)/slope
    ordered = sorted(cases, key=lambda row: row["sqrt_er_k0_b"])
    strict = [row for row in ordered
              if abs(row["error"]["frequency_percent"]) <= 1 and
              abs(row["error"]["resistance_percent"]) <= 1]
    failures = [row for row in ordered
                if abs(row["error"]["frequency_percent"]) > 1 or
                abs(row["error"]["resistance_percent"]) > 1]
    analysis = {
        "schema": "nec2c-coated-wire-m8-thickness-analysis",
        "case_count": len(cases),
        "observed_range": {
            "radius_ratio_b_over_a": [
                min(row["radius_ratio_b_over_a"] for row in cases),
                max(row["radius_ratio_b_over_a"] for row in cases)],
            "epsilon_r": [min(row["epsilon_r"] for row in cases),
                          max(row["epsilon_r"] for row in cases)],
            "sqrt_er_k0_b": [float(x.min()), float(x.max())],
            "frequency_error_percent": [float(y.min()), float(y.max())],
        },
        "affine_fit_frequency_error_vs_sqrt_er_k0_b": {
            "slope_percent_per_unit": float(slope),
            "intercept_percent": float(intercept),
            "r_squared": float(r_squared),
            "maximum_absolute_residual_percent": float(
                np.max(np.abs(residual))),
            "estimated_one_percent_crossing": float(crossing),
            "warning": ("The crossing is an interpolation aid, not a final "
                        "validity rule; mesh refinement, loss, and endpoint "
                        "sweeps remain."),
        },
        "strict_one_percent_boundary": {
            "largest_passing_sqrt_er_k0_b": max(
                row["sqrt_er_k0_b"] for row in strict),
            "first_failing_sqrt_er_k0_b": min(
                row["sqrt_er_k0_b"] for row in failures),
            "provisional_warning_threshold": 0.12,
            "criteria": ("absolute resonance-frequency and resistance "
                         "errors must both be no greater than one percent"),
        },
        "radiation_metrics": {
            "maximum_absolute_directivity_error_db": max(
                abs(row["error"]["directivity_db"]) for row in cases
                if row["error"]["directivity_db"] is not None),
            "maximum_pattern_error_db_above_minus_40_db": max(
                row["error"]["maximum_pattern_db_above_minus_40_db"]
                for row in cases
                if row["error"]["maximum_pattern_db_above_minus_40_db"]
                is not None),
        },
        "observed_knee": False,
        "interpretation": (
            "Across the current lossless flush-end grid the error rises "
            "smoothly. sqrt(Er)*k0*b collapses thickness and permittivity "
            "better than b/a alone; no abrupt physical failure is observed."),
    }
    OUTPUT.write_text(json.dumps(analysis, indent=2)+"\n")

    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    colors = {3.2: "tab:blue", 6.0: "tab:green", 10.0: "tab:red"}
    for er in sorted({row["epsilon_r"] for row in cases}):
        rows = sorted((row for row in cases if row["epsilon_r"] == er),
                      key=lambda row: row["radius_ratio_b_over_a"])
        axes[0].plot([row["radius_ratio_b_over_a"] for row in rows],
                     [row["error"]["frequency_percent"] for row in rows],
                     marker="o", label=f"Er={er:g}", color=colors.get(er))
        axes[2].plot([row["radius_ratio_b_over_a"] for row in rows],
                     [row["error"]["resistance_percent"] for row in rows],
                     marker="o", label=f"Er={er:g}", color=colors.get(er))
    for row in cases:
        axes[1].scatter(row["sqrt_er_k0_b"],
                        row["error"]["frequency_percent"],
                        color=colors.get(row["epsilon_r"]), s=38)
    fit_x = np.linspace(0, max(x)*1.08, 100)
    axes[1].plot(fit_x, slope*fit_x+intercept, color="0.25",
                 linestyle="--", label=f"affine fit, R²={r_squared:.3f}")
    axes[1].axvline(0.12, color="tab:purple", linestyle="-.",
                    label="0.12 warning threshold")
    for axis in axes[:2]:
        axis.axhline(1, color="tab:orange", linestyle=":",
                     label="1% reference band")
    axes[0].set(xlabel="Outer/conductor radius, b/a",
                ylabel="NEC resonance error (%)")
    axes[1].set(xlabel="Electrical radius in dielectric, sqrt(Er) k0 b",
                ylabel="NEC resonance error (%)")
    axes[2].set(xlabel="Outer/conductor radius, b/a",
                ylabel="NEC resistance error (%)")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=8)
    figure.suptitle("M8 coated-wire validity map (lossless, flush ends)")
    figure.tight_layout()
    figure.savefig(FIGURE, dpi=160)
    plt.close(figure)
    print(f"fit: error = {slope:.4f}*sqrt(Er)k0b + {intercept:.4f}, "
          f"R^2={r_squared:.4f}")
    print(f"estimated 1% crossing: {crossing:.4f}")
    print(OUTPUT)
    print(FIGURE)


if __name__ == "__main__":
    main()
