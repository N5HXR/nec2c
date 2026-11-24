#!/usr/bin/env python3
"""Independent homogeneous-medium scalar reference for NEC4 M7.

This is deliberately not an antenna solver or a Sommerfeld evaluator.  It
provides the analytic material constants and outgoing scalar Green function
used to validate M7's lower-medium foundation before M8 adds the planar-
interface field integrals.
"""

import argparse
import cmath
import json
import math


C0 = 299_792_458.0
ETA0 = 376.730313668
EPS0 = 8.8541878128e-12


def complex_pair(value):
    return {"real": value.real, "imag": value.imag}


def evaluate(frequency_mhz, relative_permittivity, conductivity, epsilon_loss,
             distances):
    omega = 2.0 * math.pi * frequency_mhz * 1.0e6
    if epsilon_loss is not None:
        epsilon = complex(relative_permittivity, -abs(epsilon_loss))
        loss_kind = "epsilon-imaginary"
        loss_value = abs(epsilon_loss)
    else:
        epsilon = complex(relative_permittivity,
                          -conductivity / (omega * EPS0))
        loss_kind = "conductivity-s-per-m"
        loss_value = conductivity
    refractive_index = cmath.sqrt(epsilon)
    wave_number = omega / C0 * refractive_index
    impedance = ETA0 / refractive_index
    samples = []
    for distance in distances:
        if not distance > 0.0:
            raise ValueError("every distance must be positive")
        green = cmath.exp(-1j * wave_number * distance) / (
            4.0 * math.pi * distance)
        samples.append({"distance_m": distance,
                        "green_per_m": complex_pair(green),
                        "magnitude_per_m": abs(green),
                        "phase_degrees": math.degrees(cmath.phase(green))})
    return {
        "schema": "nec2c-nec4-m7-homogeneous-reference",
        "schema_version": 1,
        "convention": "exp(+j omega t), outgoing exp(-j k r)",
        "frequency_mhz": frequency_mhz,
        "relative_permittivity": relative_permittivity,
        "loss": {"kind": loss_kind, "value": loss_value},
        "epsilon_relative": complex_pair(epsilon),
        "refractive_index": complex_pair(refractive_index),
        "wave_number_per_m": complex_pair(wave_number),
        "intrinsic_impedance_ohm": complex_pair(impedance),
        "wavelength_m": 2.0 * math.pi / wave_number.real,
        "attenuation_np_per_m": -wave_number.imag,
        "samples": samples,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequency-mhz", type=float, default=7.1)
    parser.add_argument("--er", type=float, default=13.0)
    loss = parser.add_mutually_exclusive_group()
    loss.add_argument("--sigma", type=float, default=0.005,
                      help="conductivity in S/m (default 0.005)")
    loss.add_argument("--epsilon-loss", type=float,
                      help="positive magnitude of Im(epsilon_r)")
    parser.add_argument("--distance", type=float, action="append",
                        dest="distances", default=[])
    parser.add_argument("--output")
    arguments = parser.parse_args()
    if arguments.frequency_mhz <= 0.0 or arguments.er <= 0.0:
        parser.error("frequency and relative permittivity must be positive")
    sigma = 0.0 if arguments.epsilon_loss is not None else arguments.sigma
    if sigma < 0.0:
        parser.error("conductivity must be nonnegative")
    distances = arguments.distances or [0.1, 1.0, 10.0]
    document = evaluate(arguments.frequency_mhz, arguments.er, sigma,
                        arguments.epsilon_loss, distances)
    encoded = json.dumps(document, indent=2) + "\n"
    if arguments.output:
        with open(arguments.output, "w", encoding="utf-8") as stream:
            stream.write(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
