#!/usr/bin/env python3
"""Independent high-precision witnesses for the NEC4 flat end-cap field."""

import cmath
import math
from pathlib import Path
import subprocess

import mpmath as mp


ROOT = Path(__file__).resolve().parents[2]
mp.mp.dps = 50


def disk_field(radius, radial, axial, wave_number):
    radius, radial, axial, wave_number = map(
        mp.mpf, (radius, radial, axial, wave_number))

    def component(use_radial):
        def radial_integral(source_radius):
            def azimuth_integral(phi):
                distance = mp.sqrt(radial**2 + source_radius**2 + axial**2 -
                                   2*radial*source_radius*mp.cos(phi))
                direction = (radial-source_radius*mp.cos(phi)
                             if use_radial else axial)
                return (mp.exp(-1j*wave_number*distance) *
                        (1+1j*wave_number*distance) / distance**3 *
                        direction*source_radius)
            return mp.quad(azimuth_integral, [0, 2*mp.pi])
        return mp.quad(radial_integral, [0, radius])

    return complex(component(True)), complex(component(False))


def scaled(left, right):
    return abs(left-right)/(1+max(abs(left), abs(right)))


def main():
    near = disk_field(".01", ".023", ".017", str(2*mp.pi))
    expected_near = (
        complex(0.3013680894343398246, -0.00059528019888976065),
        complex(0.2442152899578017898, -0.00044007666094748538),
    )
    errors = [scaled(a, b) for a, b in zip(near, expected_near)]
    if max(errors) > 2e-13:
        raise AssertionError((near, expected_near, errors))

    radius, axial, wave_number = .01, .037, 2*math.pi
    edge = math.hypot(radius, axial)
    axis = 2*math.pi*cmath.exp(-1j*wave_number*axial)*(
        1-axial/edge*cmath.exp(-1j*wave_number*(edge-axial)))
    if scaled(axis, complex(0.22363515245144316, -0.0009557297096793316)) > 1e-13:
        raise AssertionError(axis)

    executable = ROOT / "tests/nec4_endcap_test"
    completed = subprocess.run([str(executable)], cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((completed.stdout, completed.stderr))
    print(f"M5 end-cap direct quadrature maximum scaled error {max(errors):.3e}")
    print(completed.stdout.strip())


if __name__ == "__main__":
    main()
