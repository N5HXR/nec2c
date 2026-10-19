#!/usr/bin/env python3
"""Independent openEMS N2 dipole/cube impedance comparison for M6."""

import argparse
import csv
import pathlib
import tempfile

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0


def simulate(root: pathlib.Path, epsilon_r: float, mesh_mm: float):
    f0 = 300.0e6
    frequencies = np.linspace(240.0e6, 360.0e6, 5)
    simulation = root / ("vacuum" if epsilon_r == 1.0 else f"er{epsilon_r:g}")
    simulation.mkdir(parents=True, exist_ok=True)

    fdtd = openEMS(NrTS=80000, EndCriteria=1.0e-5, CellConstantMaterial=True)
    fdtd.SetGaussExcite(f0, 90.0e6)
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1.0e-3)

    dipole = csx.AddMetal("Dipole")
    dipole.AddCylinder([0, 0, -237.5], [0, 0, -5], radius=1.0, priority=1)
    dipole.AddCylinder([0, 0, 5], [0, 0, 237.5], radius=1.0, priority=1)
    port = fdtd.AddLumpedPort(
        port_nr=1, R=50.0, start=[-1, -1, -5], stop=[1, 1, 5],
        p_dir="z", excite=True, priority=100,
    )

    if epsilon_r != 1.0:
        cube = csx.AddMaterial("Cube", epsilon=epsilon_r)
        cube.AddBox([50, -50, -50], [150, 50, 50], priority=10)
    # Keep the Yee mesh identical between the vacuum and dielectric runs.
    for direction in "xyz":
        axis = "xyz".index(direction)
        bounds = ([50, 150], [-50, 50], [-50, 50])[axis]
        lines = np.unique(np.r_[bounds, np.arange(bounds[0], bounds[1], mesh_mm)])
        grid.AddLine(direction, lines)

    grid.AddLine("x", [-600, -2, -1, 0, 1, 2, 650])
    grid.AddLine("y", [-600, -2, -1, 0, 1, 2, 600])
    grid.AddLine("z", [-600, -237.5, -5, 0, 5, 237.5, 600])
    grid.SmoothMeshLines("all", min(20.0, 2.0 * mesh_mm), 1.35)

    csx.Write2XML(str(simulation / "n2.xml"))
    fdtd.Run(str(simulation), cleanup=True, verbose=1)
    port.CalcPort(str(simulation), frequencies)
    impedance = port.uf_tot / port.if_tot
    return frequencies, impedance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-mm", type=float, default=10.0)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--workdir", type=pathlib.Path)
    args = parser.parse_args()
    if args.mesh_mm <= 0:
        parser.error("--mesh-mm must be positive")

    temporary = None
    root = args.workdir
    if root is None:
        temporary = tempfile.TemporaryDirectory(prefix="nec2c-openems-m6-")
        root = pathlib.Path(temporary.name)
    root.mkdir(parents=True, exist_ok=True)
    frequencies, vacuum = simulate(root, 1.0, args.mesh_mm)
    _, dielectric = simulate(root, 3.0, args.mesh_mm)

    rows = []
    for frequency, base, loaded in zip(frequencies, vacuum, dielectric):
        rows.append({
            "frequency_mhz": frequency / 1.0e6,
            "vacuum_r_ohm": base.real, "vacuum_x_ohm": base.imag,
            "er3_r_ohm": loaded.real, "er3_x_ohm": loaded.imag,
            "delta_r_ohm": (loaded - base).real,
            "delta_x_ohm": (loaded - base).imag,
            "mesh_mm": args.mesh_mm,
            "openems_method": "EC-FDTD",
        })
    fields = list(rows[0])
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    else:
        writer = csv.DictWriter(__import__("sys").stdout, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    if temporary is not None:
        temporary.cleanup()


if __name__ == "__main__":
    main()
