#!/usr/bin/env python3
"""Generate the NEC4 M8 direct-field convergence sweep and static map."""

import argparse
import csv
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, bessel_contour  # noqa: E402
from nec4_m8_equations import CrossGeometry, material, relative_error  # noqa: E402


def material_cases(frequency):
    vacuum = material(frequency, 1)
    return [
        ("vacuum/vacuum", vacuum, vacuum),
        ("er4/vacuum", material(frequency, 4), vacuum),
        ("low-loss/vacuum", material(frequency, 4,
                                      conductivity=mp.mpf("1e-5")), vacuum),
        ("soil/vacuum", material(frequency, 13,
                                  conductivity=mp.mpf("0.005")), vacuum),
        ("seawater/vacuum", material(frequency, 80,
                                      conductivity=mp.mpf("4")), vacuum),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=13)
    parser.add_argument("--output-json", type=Path,
                        default=ROOT / "tests" / "nec4-similarity" /
                        "m8-qualification.json")
    parser.add_argument("--output-csv", type=Path,
                        default=ROOT / "tests" / "nec4-similarity" /
                        "m8-qualification.csv")
    parser.add_argument("--output-figure", type=Path,
                        default=ROOT / "docs" / "nec4-m8-convergence.png")
    arguments = parser.parse_args()
    if arguments.points < 5:
        raise SystemExit("--points must be at least 5")

    mp.mp.dps = 80
    frequency = mp.mpf("7.1e6")
    vacuum = material(frequency, 1)
    k0 = mp.re(vacuum.wave_number)
    maximum_log = mp.log10(mp.mpf(200))
    logarithms = [mp.mpf(-5) + mp.mpf(index) * (maximum_log + 5) /
                  (arguments.points - 1)
                  for index in range(arguments.points)]
    krho_values = [mp.power(10, value) for value in logarithms]
    rows = []
    matrix = []

    for case_name, lower, upper in material_cases(frequency):
        case_errors = []
        for krho in krho_values:
            rho = krho / k0
            geometry = CrossGeometry(rho, -mp.mpf("0.65") * rho,
                                     mp.mpf("0.45") * rho)
            # The far-field saddle moves outward with k0*rho.  A fixed tail
            # cutoff that is ample in the near field can terminate before the
            # saddle contribution at k0*rho~100.  Keep two independently
            # converged levels beyond that transition.
            if mp.mpf("5e-4") <= krho <= mp.mpf("2e-2"):
                # Near-singular transition: the leading 1/R terms cancel and
                # the absolute tail error is most visible in z_v/rho_h.
                coarse_decay, refined_decay = 160, 200
                coarse_subdivisions, refined_subdivisions = 6, 8
                coarse_order, refined_order = 72, 76
            elif krho > 20:
                coarse_decay = int(mp.ceil(krho)) + 20
                refined_decay = int(mp.ceil(krho)) + 35
                coarse_subdivisions, refined_subdivisions = 4, 5
                coarse_order, refined_order = 60, 68
            elif krho > 5:
                coarse_decay, refined_decay = 100, 130
                coarse_subdivisions, refined_subdivisions = 3, 5
                coarse_order, refined_order = 60, 68
            else:
                coarse_decay, refined_decay = 100, 130
                coarse_subdivisions, refined_subdivisions = 3, 5
                coarse_order, refined_order = 60, 68
            coarse = bessel_contour(
                frequency, lower, upper, geometry,
                precision_digits=50, cutoff_decay=coarse_decay,
                panel_subdivisions=coarse_subdivisions,
                quadrature_order=coarse_order)
            refined = bessel_contour(
                frequency, lower, upper, geometry,
                precision_digits=65, cutoff_decay=refined_decay,
                panel_subdivisions=refined_subdivisions,
                quadrature_order=refined_order)
            errors = {name: relative_error(refined.components[name],
                                           coarse.components[name])
                      for name in COMPONENTS}
            worst_name = max(errors, key=errors.get)
            worst = errors[worst_name]
            case_errors.append(float(worst))
            rows.append({
                "material": case_name,
                "k0_rho": mp.nstr(krho, 16),
                "rho_m": mp.nstr(rho, 16),
                "worst_component": worst_name,
                "worst_refinement_error": mp.nstr(worst, 16),
                "coarse_cutoff_decay": coarse_decay,
                "refined_cutoff_decay": refined_decay,
                **{f"error_{name}": mp.nstr(value, 16)
                   for name, value in errors.items()},
            })
        matrix.append(case_errors)
        print(f"qualified {case_name}: worst {max(case_errors):.3e}", flush=True)

    worst_row = max(rows, key=lambda row: float(row["worst_refinement_error"]))
    target = mp.mpf("1e-9")
    if mp.mpf(worst_row["worst_refinement_error"]) > target:
        raise SystemExit("qualification target exceeded at " + json.dumps(worst_row))

    document = {
        "schema": "nec2c-nec4-m8-qualification",
        "schema_version": 1,
        "frequency_hz": mp.nstr(frequency, 16),
        "geometry": "source_z=-0.65*rho, observer_z=+0.45*rho",
        "range": {"minimum_k0_rho": "1e-5", "maximum_k0_rho": "200"},
        "target": mp.nstr(target, 8),
        "worst_point": worst_row,
        "rows": rows,
    }
    arguments.output_json.write_text(
        json.dumps(document, indent=2) + "\n", encoding="utf-8")
    with arguments.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure, axis = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)
    image = axis.imshow(
        [[float(max(-16, min(-6, mp.log10(max(value, 1e-30)))))
          for value in case] for case in matrix],
        aspect="auto", origin="lower", cmap="viridis",
        extent=[float(logarithms[0]), float(logarithms[-1]),
                -0.5, len(matrix) - 0.5], vmin=-16, vmax=-6)
    axis.set_yticks(range(len(matrix)))
    axis.set_yticklabels([name for name, _, _ in material_cases(frequency)])
    axis.set_xlabel(r"$\log_{10}(k_0\rho)$")
    axis.set_ylabel("lower / upper medium")
    axis.set_title("M8 direct Sommerfeld external refinement error")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(r"$\log_{10}$ worst component error")
    axis.axvline(mp.log10(mp.mpf("1e-5")), color="white", alpha=0.25)
    figure.savefig(arguments.output_figure, dpi=170)
    plt.close(figure)
    print(f"{arguments.output_json}\n{arguments.output_csv}\n{arguments.output_figure}")
    print("worst qualification error " + worst_row["worst_refinement_error"] +
          " at " + worst_row["material"] + " k0*rho=" + worst_row["k0_rho"])


if __name__ == "__main__":
    main()
