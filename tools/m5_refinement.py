#!/usr/bin/env python3
"""Run reproducible M5 cuboid and sphere refinement studies."""

import argparse
import math
import os
import pathlib
import re
import subprocess
import tempfile

EPS0 = 8.8541878128e-12
ETA0 = math.sqrt(1.25663706212e-6 / EPS0)
C0 = 299792458.0


def scalar(label, report):
    match = re.search(re.escape(label) + r"\s*([-+0-9.E]+)", report)
    if not match:
        raise RuntimeError(f"missing output value: {label}")
    return float(match.group(1))


def dipole(report):
    match = re.search(
        r"PX\(RE,IM\).*?\n\s*([-+0-9.E]+)\s+([-+0-9.E]+)\s+"
        r"([-+0-9.E]+)\s+([-+0-9.E]+)\s+([-+0-9.E]+)\s+([-+0-9.E]+)",
        report, re.S,
    )
    if not match:
        raise RuntimeError("missing dipole row")
    values = [float(value) for value in match.groups()]
    return tuple(complex(values[index], values[index + 1]) for index in range(0, 6, 2))


def mie_cross_section(radius, epsilon_r, frequency_mhz):
    k0 = 2.0*math.pi*frequency_mhz*1.0e6/C0
    x = k0*radius
    refractive_index = math.sqrt(epsilon_r)

    def spherical_bessel(argument, limit):
        j = [math.sin(argument)/argument,
             math.sin(argument)/argument**2-math.cos(argument)/argument]
        y = [-math.cos(argument)/argument,
             -math.cos(argument)/argument**2-math.sin(argument)/argument]
        for order in range(1, limit):
            j.append((2*order+1)*j[order]/argument-j[order-1])
            y.append((2*order+1)*y[order]/argument-y[order-1])
        return j, y

    limit = max(3, int(x+4*x**(1/3)+3))
    jx, yx = spherical_bessel(x, limit)
    jmx, _ = spherical_bessel(refractive_index*x, limit)
    total = 0.0
    for order in range(1, limit+1):
        psi_x = x*jx[order]
        psi_m = refractive_index*x*jmx[order]
        dpsi_x = x*jx[order-1]-order*jx[order]
        dpsi_m = refractive_index*x*jmx[order-1]-order*jmx[order]
        xi = complex(x*jx[order], x*yx[order])
        dxi = complex(x*jx[order-1]-order*jx[order],
                      x*yx[order-1]-order*yx[order])
        a = ((refractive_index*psi_m*dpsi_x-psi_x*dpsi_m) /
             (refractive_index*psi_m*dxi-xi*dpsi_m))
        b = ((psi_m*dpsi_x-refractive_index*psi_x*dpsi_m) /
             (psi_m*dxi-refractive_index*xi*dpsi_m))
        total += (2*order+1)*(abs(a)**2+abs(b)**2)
    return 2.0*math.pi*total/(k0*k0)


def boxes_for_grid(radius, divisions, sphere):
    edge = 2.0*radius/divisions
    boxes = []
    for ix in range(divisions):
        for iy in range(divisions):
            for iz in range(divisions):
                center = [-radius+(index+0.5)*edge for index in (ix, iy, iz)]
                if sphere and sum(value*value for value in center) > radius*radius:
                    continue
                low = [(index-divisions/2)*edge for index in (ix, iy, iz)]
                high = [(index+1-divisions/2)*edge for index in (ix, iy, iz)]
                boxes.append("BOX 1 " + " ".join(f"{value:.14g}" for value in low+high))
    return edge, boxes


def solve(binary, root, name, boxes, frequency_mhz):
    resource = root/f"{name}.diel"
    deck = root/f"{name}.nec"
    output = root/f"{name}.out"
    resource.write_text(
        "nec2c-dielectric 1\nMATERIAL 1 ER 3\nREGION 1 MATERIAL 1\n"+
        "\n".join(boxes)+"\nMESH MAX_KH 1 MAX_CELLS 500\nEND\n"
    )
    deck.write_text(
        "CM M5 REFINEMENT STUDY\nCE\nGE\n"
        f"FR 0 1 0 0 {frequency_mhz:.14g}\n"
        "EX 1 1 1 0 90 0 90 0 0 0\n"
        f"DF \"{resource.name}\"\nEN\n"
    )
    result = subprocess.run([str(binary), "-i", str(deck), "-o", str(output)],
                            text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"{name} failed:\n{result.stderr}\n{output.read_text()}")
    report = output.read_text()
    return {
        "cells": int(scalar("ACTIVE CELLS:", report)),
        "dipole": abs(dipole(report)[1]),
        "scattered": scalar("SCATTERED POWER  (W):", report),
        "residual": scalar("RELATIVE RESIDUAL:", report),
        "imbalance": scalar("RELATIVE IMBALANCE:", report),
    }


def percent(value):
    return f"{100.0*value:.4f}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nec2c", default=os.environ.get("NEC2C", "./nec2c"))
    parser.add_argument("--maximum-division", type=int, default=5)
    arguments = parser.parse_args()
    binary = pathlib.Path(arguments.nec2c).resolve()
    radius = 0.02
    epsilon_r = 3.0
    with tempfile.TemporaryDirectory(prefix="nec2c-m5-refinement-") as temporary:
        root = pathlib.Path(temporary)

        print("CUBE DECOMPOSITION REFINEMENT (100 MHz)")
        print("n,cells,|p_y|_C_m,change_from_previous_pct,Psca_W,residual,power_imbalance")
        previous = None
        for divisions in range(1, min(arguments.maximum_division, 4)+1):
            _, boxes = boxes_for_grid(radius, divisions, False)
            result = solve(binary, root, f"cube-{divisions}", boxes, 100.0)
            change = float("nan") if previous is None else result["dipole"]/previous-1.0
            print(f"{divisions},{result['cells']},{result['dipole']:.9e},"
                  f"{percent(change) if previous is not None else 'NA'},"
                  f"{result['scattered']:.9e},{result['residual']:.3e},"
                  f"{result['imbalance']:.3e}")
            previous = result["dipole"]

        print("\nRAYLEIGH SPHERE BOUNDARY REFINEMENT (requested k0*a=0.08)")
        print("n,cells,volume_ratio,equivalent_radius_m,|p_y|_C_m,rayleigh_error_pct,mean_E_error_pct,residual,power_imbalance")
        low_frequency = 0.08*C0/(2.0*math.pi*radius)/1.0e6
        for divisions in range(2, arguments.maximum_division+1):
            edge, boxes = boxes_for_grid(radius, divisions, True)
            volume = len(boxes)*edge**3
            volume_ratio = volume/(4.0*math.pi*radius**3/3.0)
            equivalent_radius = (3.0*volume/(4.0*math.pi))**(1.0/3.0)
            result = solve(binary, root, f"rayleigh-{divisions}", boxes, low_frequency)
            expected = 4.0*math.pi*EPS0*equivalent_radius**3*(epsilon_r-1)/(epsilon_r+2)
            mean_field = result["dipole"]/(EPS0*(epsilon_r-1)*volume)
            exact_field = 3.0/(epsilon_r+2)
            print(f"{divisions},{result['cells']},{volume_ratio:.8f},{equivalent_radius:.9e},"
                  f"{result['dipole']:.9e},{percent(result['dipole']/expected-1.0)},"
                  f"{percent(mean_field/exact_field-1.0)},{result['residual']:.3e},"
                  f"{result['imbalance']:.3e}")

        print("\nDYNAMIC SPHERE REFINEMENT (requested k0*a=0.8)")
        print("n,cells,volume_ratio,equivalent_k0a,Psca_W,Mie_error_pct,residual,power_imbalance")
        dynamic_frequency = 0.8*C0/(2.0*math.pi*radius)/1.0e6
        k0 = 2.0*math.pi*dynamic_frequency*1.0e6/C0
        for divisions in range(3, arguments.maximum_division+1):
            edge, boxes = boxes_for_grid(radius, divisions, True)
            volume = len(boxes)*edge**3
            volume_ratio = volume/(4.0*math.pi*radius**3/3.0)
            equivalent_radius = (3.0*volume/(4.0*math.pi))**(1.0/3.0)
            result = solve(binary, root, f"dynamic-{divisions}", boxes, dynamic_frequency)
            exact_power = mie_cross_section(equivalent_radius, epsilon_r,
                                             dynamic_frequency)/(2.0*ETA0)
            print(f"{divisions},{result['cells']},{volume_ratio:.8f},{k0*equivalent_radius:.8f},"
                  f"{result['scattered']:.9e},{percent(result['scattered']/exact_power-1.0)},"
                  f"{result['residual']:.3e},{result['imbalance']:.3e}")


if __name__ == "__main__":
    main()
