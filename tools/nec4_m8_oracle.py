#!/usr/bin/env python3
"""Diagnostic NEC4 M8 direct Sommerfeld electric-field oracle.

This is deliberately separate from the production solver.  It evaluates one
point or a short observer-height line with the public Section 5 / Appendix G
formulas and reports an external refinement difference.
"""

import argparse
import json
from pathlib import Path
import sys

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m8_direct import COMPONENTS, transformed_pairing  # noqa: E402
from nec4_m8_equations import assemble_cartesian, material, relative_error  # noqa: E402


def complex_record(value, digits=24):
    return {"real": mp.nstr(mp.re(value), digits),
            "imag": mp.nstr(mp.im(value), digits),
            "magnitude": mp.nstr(abs(value), digits),
            "phase_degrees": mp.nstr(mp.degrees(mp.arg(value)), digits)}


def pairing_name(source_z, observer_z):
    return ("1" if source_z < 0 else "2") + ("1" if observer_z < 0 else "2")


def evaluate(arguments, observer_z):
    lower = material(arguments.frequency, arguments.lower_er,
                     conductivity=arguments.lower_conductivity)
    upper = material(arguments.frequency, arguments.upper_er,
                     conductivity=arguments.upper_conductivity)
    common = dict(method=arguments.method)
    coarse = transformed_pairing(
        arguments.frequency, lower, upper, arguments.rho,
        arguments.source_z, observer_z,
        precision_digits=max(45, arguments.precision - 15),
        cutoff_decay=max(50, arguments.cutoff - 20),
        panel_subdivisions=max(1, arguments.subdivisions - 1),
        quadrature_order=max(28, arguments.order - 12), **common)
    refined = transformed_pairing(
        arguments.frequency, lower, upper, arguments.rho,
        arguments.source_z, observer_z,
        precision_digits=arguments.precision,
        cutoff_decay=arguments.cutoff,
        panel_subdivisions=arguments.subdivisions,
        quadrature_order=arguments.order, **common)
    errors = {name: relative_error(refined[name], coarse[name])
              for name in COMPONENTS}
    dyadic = assemble_cartesian(refined, mp.radians(arguments.phi_degrees))
    return {
        "pairing": pairing_name(arguments.source_z, observer_z),
        "geometry_m": {
            "rho": mp.nstr(arguments.rho, 20),
            "phi_degrees": mp.nstr(arguments.phi_degrees, 20),
            "source_z": mp.nstr(arguments.source_z, 20),
            "observer_z": mp.nstr(observer_z, 20),
        },
        "components": {name: complex_record(refined[name])
                       for name in COMPONENTS},
        "cartesian_electric_dyadic": [
            [complex_record(dyadic[row, column]) for column in range(3)]
            for row in range(3)],
        "refinement_error": {name: mp.nstr(value, 14)
                             for name, value in errors.items()},
        "worst_refinement_error": mp.nstr(max(errors.values()), 14),
    }


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("mode", choices=("point", "line"))
    result.add_argument("--frequency", type=mp.mpf, default=mp.mpf("7.1e6"))
    result.add_argument("--lower-er", type=mp.mpf, default=mp.mpf("13"))
    result.add_argument("--lower-conductivity", type=mp.mpf,
                        default=mp.mpf("0.005"))
    result.add_argument("--upper-er", type=mp.mpf, default=mp.mpf("1"))
    result.add_argument("--upper-conductivity", type=mp.mpf, default=mp.mpf("0"))
    result.add_argument("--rho", type=mp.mpf, default=mp.mpf("1"))
    result.add_argument("--phi-degrees", type=mp.mpf, default=mp.mpf("0"))
    result.add_argument("--source-z", type=mp.mpf, default=mp.mpf("-0.7"))
    result.add_argument("--observer-z", type=mp.mpf, default=mp.mpf("0.4"))
    result.add_argument("--line-start", type=mp.mpf, default=mp.mpf("-0.8"))
    result.add_argument("--line-stop", type=mp.mpf, default=mp.mpf("0.8"))
    result.add_argument("--line-count", type=int, default=5)
    result.add_argument("--method", choices=("auto", "bessel", "hankel"),
                        default="auto")
    result.add_argument("--precision", type=int, default=65)
    result.add_argument("--cutoff", type=int, default=80)
    result.add_argument("--subdivisions", type=int, default=3)
    result.add_argument("--order", type=int, default=48)
    result.add_argument("--output", type=Path)
    return result


def main():
    arguments = parser().parse_args()
    if arguments.line_count < 2:
        raise SystemExit("--line-count must be at least 2")
    if arguments.source_z == 0:
        raise SystemExit("source must be off the interface")
    mp.mp.dps = max(80, arguments.precision + 10)
    if arguments.mode == "point":
        if arguments.observer_z == 0:
            raise SystemExit("observer must be off the interface")
        points = [evaluate(arguments, arguments.observer_z)]
    else:
        points = []
        for index in range(arguments.line_count):
            observer_z = (arguments.line_start +
                          (arguments.line_stop - arguments.line_start) *
                          mp.mpf(index) / (arguments.line_count - 1))
            # The direct formulas are one-sided.  Preserve the line shape but
            # do not silently choose a side for a point exactly at z=0.
            if observer_z == 0:
                continue
            points.append(evaluate(arguments, observer_z))
    document = {
        "schema": "nec2c-nec4-m8-direct-oracle",
        "schema_version": 1,
        "status": "reference-only; not connected to production matrix fill",
        "frequency_hz": mp.nstr(arguments.frequency, 20),
        "materials": {
            "lower": {"epsilon_real": mp.nstr(arguments.lower_er, 20),
                      "conductivity_s_per_m": mp.nstr(
                          arguments.lower_conductivity, 20)},
            "upper": {"epsilon_real": mp.nstr(arguments.upper_er, 20),
                      "conductivity_s_per_m": mp.nstr(
                          arguments.upper_conductivity, 20)},
        },
        "integration": {
            "method": arguments.method,
            "precision_decimal_digits": arguments.precision,
            "cutoff_decay_lengths": arguments.cutoff,
            "panel_subdivisions": arguments.subdivisions,
            "quadrature_order": arguments.order,
            "external_error_estimate": "difference from independently coarser run",
        },
        "points": points,
    }
    text = json.dumps(document, indent=2) + "\n"
    if arguments.output:
        arguments.output.write_text(text, encoding="utf-8")
        print(arguments.output)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
