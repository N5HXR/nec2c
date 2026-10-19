#!/usr/bin/env python3
"""M5 end-to-end validation for the dielectric-only dense JVIE solver."""

import math
import os
import pathlib
import re
import subprocess
import tempfile
import csv

NEC2C = pathlib.Path(os.environ.get("NEC2C", "./nec2c")).resolve()
SOURCE = pathlib.Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
EPS0 = 8.8541878128e-12
ETA0 = math.sqrt(1.25663706212e-6 / EPS0)


def run(deck: pathlib.Path) -> str:
    with tempfile.NamedTemporaryFile(prefix="nec2c-m5-output-", suffix=".out") as temporary:
        output = pathlib.Path(temporary.name)
        completed = subprocess.run(
            [str(NEC2C), "-i", str(deck), "-o", str(output)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        result = output.read_text()
        if completed.returncode != 0:
            raise AssertionError(f"nec2c failed for {deck}: {completed.stderr}\n{result}")
        return result


def number_after(label: str, text: str) -> float:
    match = re.search(re.escape(label) + r"\s*([-+0-9.E]+)", text)
    if not match:
        raise AssertionError(f"missing {label!r}")
    return float(match.group(1))


def dipole(text: str):
    match = re.search(
        r"PX\(RE,IM\).*?\n\s*([-+0-9.E]+)\s+([-+0-9.E]+)\s+"
        r"([-+0-9.E]+)\s+([-+0-9.E]+)\s+([-+0-9.E]+)\s+([-+0-9.E]+)",
        text, re.S,
    )
    if not match:
        raise AssertionError("missing dipole row")
    values = [float(value) for value in match.groups()]
    return tuple(complex(values[i], values[i + 1]) for i in range(0, 6, 2))


def write_case(directory: pathlib.Path, name: str, resource: str, ex: str,
               frequency_mhz: float = 100.0, pattern: str = "") -> pathlib.Path:
    resource_path = directory / f"{name}.diel"
    deck_path = directory / f"{name}.nec"
    resource_path.write_text(resource)
    deck_path.write_text(
        "CM GENERATED M5 VALIDATION\nCE\nGE\n"
        f"FR 0 1 0 0 {frequency_mhz}\n{ex}\nDF \"{resource_path.name}\"\n"
        f"{pattern}EN\n"
    )
    return deck_path


def mie_cross_section(radius: float, epsilon_r: float, frequency_mhz: float) -> float:
    """Exact nonmagnetic-sphere Lorenz-Mie scattering cross section."""
    c0 = 299792458.0
    k = 2.0 * math.pi * frequency_mhz * 1.0e6 / c0
    x = k * radius
    m = math.sqrt(epsilon_r)

    def sph_jy(z: float, limit: int):
        j = [math.sin(z) / z]
        y = [-math.cos(z) / z]
        j.append(math.sin(z) / z**2 - math.cos(z) / z)
        y.append(-math.cos(z) / z**2 - math.sin(z) / z)
        for n in range(1, limit):
            j.append((2*n + 1) * j[n] / z - j[n - 1])
            y.append((2*n + 1) * y[n] / z - y[n - 1])
        return j, y

    limit = max(3, int(x + 4*x**(1/3) + 3))
    jx, yx = sph_jy(x, limit)
    jmx, _ = sph_jy(m*x, limit)
    total = 0.0
    for n in range(1, limit + 1):
        psi_x, psi_m = x*jx[n], m*x*jmx[n]
        dpsi_x = x*jx[n - 1] - n*jx[n]
        dpsi_m = m*x*jmx[n - 1] - n*jmx[n]
        xi = complex(x*jx[n], x*yx[n])
        dxi = complex(x*jx[n - 1] - n*jx[n], x*yx[n - 1] - n*yx[n])
        a = (m*psi_m*dpsi_x - psi_x*dpsi_m) / (m*psi_m*dxi - xi*dpsi_m)
        b = (psi_m*dpsi_x - m*psi_x*dpsi_m) / (psi_m*dxi - m*xi*dpsi_m)
        total += (2*n + 1) * (abs(a)**2 + abs(b)**2)
    return 2.0 * math.pi * total / (k*k)


def main():
    fixture = run(SOURCE / "tests/dielectric/m5-cube.nec")
    assert number_after("RELATIVE RESIDUAL:", fixture) <= 1.0e-10
    assert number_after("RELATIVE IMBALANCE:", fixture) <= 1.0e-4
    p = dipole(fixture)
    expected_cube = EPS0 * (0.02**3) * 3.0 * (3.0 - 1.0) / (3.0 + 2.0)
    assert abs(abs(p[1]) / expected_cube - 1.0) < 2.0e-3
    assert fixture.count("BISTATIC RCS") == 1
    pattern_rows = re.findall(
        r"^\s*(0\.000|45\.000|90\.000)\s+90\.000\s+"
        r"[-+0-9.E]+\s+[-+0-9.E]+\s+[-+0-9.E]+\s+[-+0-9.E]+\s+([-+0-9.E]+)$",
        fixture, re.M,
    )
    pattern = {float(angle): float(rcs) for angle, rcs in pattern_rows}
    assert abs(pattern[45.0]/pattern[0.0] - 0.5) < 2.0e-4
    assert pattern[90.0]/pattern[0.0] < 1.0e-20

    with tempfile.TemporaryDirectory(prefix="nec2c-m5-") as temporary:
        root = pathlib.Path(temporary)
        vacuum = write_case(
            root, "vacuum",
            "nec2c-dielectric 1\nMATERIAL 1 ER 1\nREGION 1 MATERIAL 1\n"
            "BOX 1 -0.01 -0.01 -0.01 0.01 0.01 0.01\nEND\n",
            "EX 1 1 1 0 90 0 90 0 0 0",
            pattern="RP 0 1 1 0000 45 0 0 0 0 0\n",
        )
        vacuum_output = run(vacuum)
        assert "VACUUM IDENTITY" in vacuum_output
        assert all(value == 0 for value in dipole(vacuum_output))
        assert number_after("SCATTERED POWER  (W):", vacuum_output) == 0.0

        rotated = write_case(
            root, "rotated",
            (SOURCE / "tests/dielectric/m5-cube.diel").read_text(),
            "EX 1 1 1 0 90 0 0 0 0 0",
        )
        rotated_p = dipole(run(rotated))
        assert abs(abs(rotated_p[2]) / abs(p[1]) - 1.0) < 1.0e-12

        # A 4-per-diameter, equal-cube staircase sphere.  The equivalent-volume
        # radius separates occupied-shape error from the field/operator error.
        radius = 0.02
        divisions = 4
        edge = 2.0 * radius / divisions
        boxes = []
        for ix in range(divisions):
            for iy in range(divisions):
                for iz in range(divisions):
                    center = [-radius + (index + 0.5)*edge for index in (ix, iy, iz)]
                    if sum(value*value for value in center) <= radius*radius:
                        lo = [(index - divisions/2)*edge for index in (ix, iy, iz)]
                        hi = [(index + 1 - divisions/2)*edge for index in (ix, iy, iz)]
                        boxes.append("BOX 1 " + " ".join(f"{value:.12g}" for value in lo + hi))
        volume = len(boxes) * edge**3
        equivalent_radius = (3.0*volume/(4.0*math.pi))**(1.0/3.0)
        sphere_resource = (
            "nec2c-dielectric 1\nMATERIAL 1 ER 3\nREGION 1 MATERIAL 1\n" +
            "\n".join(boxes) + "\nMESH MAX_KH 1.0 MAX_CELLS 100\nEND\n"
        )
        sphere = write_case(
            root, "sphere4", sphere_resource,
            "EX 1 1 1 0 90 0 90 0 0 0", frequency_mhz=200.0,
        )
        sphere_output = run(sphere)
        sphere_p = abs(dipole(sphere_output)[1])
        rayleigh_p = 4.0*math.pi*EPS0*equivalent_radius**3*(3.0-1.0)/(3.0+2.0)
        assert abs(sphere_p/rayleigh_p - 1.0) < 0.08
        mean_internal = sphere_p / (EPS0*(3.0-1.0)*volume)
        assert abs(mean_internal - 3.0/(3.0+2.0)) < 0.05
        mie_power = mie_cross_section(equivalent_radius, 3.0, 200.0)/(2.0*ETA0)
        computed_power = number_after("SCATTERED POWER  (W):", sphere_output)
        assert abs(computed_power/mie_power - 1.0) < 0.10
        assert number_after("RELATIVE RESIDUAL:", sphere_output) <= 1.0e-10
        assert number_after("RELATIVE IMBALANCE:", sphere_output) <= 2.0e-3

        dynamic_frequency = 0.8*299792458.0/(2.0*math.pi*radius)/1.0e6
        dynamic = write_case(
            root, "sphere4_dynamic", sphere_resource,
            "EX 1 1 1 0 90 0 90 0 0 0", frequency_mhz=dynamic_frequency,
            pattern="RP 0 5 1 0000 0 0 30 0 0 0\n",
        )
        dynamic_output = run(dynamic)
        exact_power = mie_cross_section(equivalent_radius, 3.0, dynamic_frequency)/(2.0*ETA0)
        numerical_power = number_after("SCATTERED POWER  (W):", dynamic_output)
        assert abs(numerical_power/exact_power - 1.0) < 0.02, (numerical_power, exact_power)
        assert number_after("RELATIVE IMBALANCE:", dynamic_output) <= 0.03
        # Independent VoxScatter DDA fixture, generated by
        # tools/voxscatter_m5_reference.m at upstream master (32 voxels,
        # requested-radius size parameter 0.8, epsilon_r=3).
        with (SOURCE / "tests/voxscatter-reference.csv").open() as reference_file:
            rows = (line for line in reference_file if not line.startswith("#"))
            reference = next(csv.DictReader(rows))
        assert int(reference["cells"]) == len(boxes)
        assert int(reference["solver_flag"]) == 0
        voxscatter_dipole = abs(complex(float(reference["pz_re"]),
                                        float(reference["pz_im"])))
        numerical_dipole = abs(dipole(dynamic_output)[1])/(EPS0*radius**3)
        assert abs(numerical_dipole/voxscatter_dipole - 1.0) < 0.03, (
            numerical_dipole, voxscatter_dipole,
        )

    print("M5 dielectric solver validation passed")


if __name__ == "__main__":
    main()
