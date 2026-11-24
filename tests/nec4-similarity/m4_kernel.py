#!/usr/bin/env python3
"""Validate the NEC4 reduced wire kernel against an independent MP oracle."""

import json
import math
import pathlib
import subprocess

import mpmath as mp


ROOT = pathlib.Path(__file__).resolve().parents[2]
ETA_OVER_4PI = mp.mpf("29.97922085")
FIELDS = ("ezs", "ers", "ezc_minus_ezk", "erc", "ezk", "erk")


def oracle(case):
    mp.mp.dps = 80
    s = mp.mpf(str(case["length"]))
    z = mp.mpf(str(case["z"]))
    rho = mp.mpf(str(case["rho"]))
    k = mp.mpf(str(case.get("wave_number", 2 * math.pi)))
    delta = s / 2
    r1 = mp.sqrt(rho * rho + (z + delta) ** 2)
    r2 = mp.sqrt(rho * rho + (z - delta) ** 2)
    e1 = mp.exp(-1j * k * r1)
    e2 = mp.exp(-1j * k * r2)
    g1 = e1 / r1
    g2 = e2 / r2
    sd = mp.sin(k * delta)
    cd = mp.cos(k * delta)
    p = -1j * ETA_OVER_4PI
    q = -p

    integral = mp.quad(
        lambda zp: mp.exp(-1j * k * mp.sqrt(rho * rho + (z - zp) ** 2))
        / mp.sqrt(rho * rho + (z - zp) ** 2),
        [-delta, delta],
    )
    return {
        "ezs": q * (g2 - g1) * cd,
        "ers": p
        / rho
        * (((z - delta) * g2 - (z + delta) * g1) * cd
           - 1j * (e2 + e1) * sd),
        "ezc_minus_ezk": p * ((g2 + g1) * sd - k * integral),
        "erc": q
        / rho
        * (((z - delta) * g2 + (z + delta) * g1) * sd
           + 1j * (e2 - e1) * cd),
        "ezk": p * k * integral,
        "erk": mp.mpc(0),
    }


def complex_pair(pair):
    return mp.mpc(mp.mpf(str(pair[0])), mp.mpf(str(pair[1])))


def scaled_error(actual, expected):
    return float(abs(actual - expected) / max(mp.mpf(1), abs(expected)))


def sample(length, z, rho, wave_number=2 * math.pi):
    output = subprocess.check_output(
        [
            str(ROOT / "tests" / "nec4_kernel_test"),
            f"{length:.17g}",
            f"{z:.17g}",
            f"{rho:.17g}",
            f"{wave_number:.17g}",
        ],
        text=True,
    )
    return json.loads(output)


def case_error(case):
    expected = oracle(case)
    return max(
        scaled_error(complex_pair(case[field]), expected[field]) for field in FIELDS
    )


def main():
    output = subprocess.check_output(
        [str(ROOT / "tests" / "nec4_kernel_test")], text=True
    )
    cases = json.loads(output)
    worst = (0.0, "", "")
    for case in cases:
        expected = oracle(case)
        errors = {}
        for field, reference in expected.items():
            error = scaled_error(complex_pair(case[field]), reference)
            errors[field] = error
            worst = max(worst, (error, case["name"], field))
        print(
            f"{case['name']:16s} {case['branch']:14s} "
            f"max scaled error {max(errors.values()):.3e}"
        )
    print(
        f"M4 worst scaled field error {worst[0]:.3e} "
        f"({worst[1]} {worst[2]})"
    )
    if worst[0] > 1.0e-11:
        raise SystemExit("M4 reduced-kernel accuracy gate failed")

    # Exercise both sides of every production switch.  The geometry changes by
    # only one part per million, so the oracle comparison measures both branch
    # accuracy and the absence of a numerical discontinuity at the handoff.
    epsilon = 1.0e-6
    boundaries = [
        ("small/exact", lambda f: (2e-5, 8e-5 * f, 6e-5 * f),
         ("small-kR", "exact")),
        ("exact/axis", lambda f: (.1, .08, 3e-5 * f),
         ("near-axis", "exact")),
        ("small/axis", lambda f: (2e-5, 3e-5, 1.00001250023438e-7 * f),
         ("small-kR-axis", "small-kR")),
        ("exact/far", lambda f: (1e-7 / f, .04, .03),
         ("exact", "far")),
        ("far-small/far", lambda f: (1e-11, 8e-5 * f, 6e-5 * f),
         ("far-small-kR", "far")),
    ]
    boundary_worst = 0.0
    for name, geometry, expected_branches in boundaries:
        lower = sample(*geometry(1.0 - epsilon))
        upper = sample(*geometry(1.0 + epsilon))
        actual_branches = (lower["branch"], upper["branch"])
        if actual_branches != expected_branches:
            raise SystemExit(
                f"M4 {name} branch gate failed: {actual_branches} != "
                f"{expected_branches}"
            )
        lower_error = case_error(lower)
        upper_error = case_error(upper)
        error = max(lower_error, upper_error)
        boundary_worst = max(boundary_worst, error)
        print(
            f"{name:16s} {actual_branches[0]} -> {actual_branches[1]} "
            f"boundary error {error:.3e} ({lower_error:.3e}/{upper_error:.3e})"
        )
    print(f"M4 worst switch-boundary error {boundary_worst:.3e}")
    if boundary_worst > 1.0e-11:
        raise SystemExit("M4 switch-boundary accuracy gate failed")

    quadrature_cases = [
        ("series-radial", .1, 0., .75 * (1 + epsilon), 2 * math.pi),
        ("series-axial", .1, .75 * (1 + epsilon), .02, 2 * math.pi),
        ("adaptive-far-r", .1, 0., .75 * (1 + epsilon), 20.),
        ("adaptive-far-z", .1, .75 * (1 + epsilon), .02, 20.),
        ("adaptive-near-r", .1, 0., .05 * (1 + epsilon), 100.),
        ("adaptive-near", .1, .03, .02, 100.),
    ]
    quadrature_worst = 0.0
    for name, length, z, rho, wave_number in quadrature_cases:
        error = case_error(sample(length, z, rho, wave_number))
        quadrature_worst = max(quadrature_worst, error)
        print(f"{name:16s} quadrature error {error:.3e}")
    if quadrature_worst > 1.0e-11:
        raise SystemExit("M4 scalar-integral quadrature gate failed")


if __name__ == "__main__":
    main()
