#!/usr/bin/env python3
"""Independent openEMS comparison for the M7 mixed PEC/dielectric model."""

import argparse
import csv
import pathlib
import tempfile

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


def simulate(root: pathlib.Path, epsilon_r: float, mesh_mm: float):
    frequency = np.array([300.0e6])
    simulation = root / ("vacuum" if epsilon_r == 1.0 else f"er{epsilon_r:g}")
    simulation.mkdir(parents=True, exist_ok=True)
    fdtd = openEMS(NrTS=80000, EndCriteria=1.0e-5, CellConstantMaterial=True)
    fdtd.SetGaussExcite(300.0e6, 90.0e6)
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1.0e-3)

    metal = csx.AddMetal("Conductors")
    metal.AddCylinder([0, 0, -250], [0, 0, -5], radius=1.0, priority=1)
    metal.AddCylinder([0, 0, 5], [0, 0, 250], radius=1.0, priority=1)
    # Zero-thickness PEC sheet corresponding to the NEC patch surface x=0.2 m.
    metal.AddBox([200, -200, -200], [200, 200, 200], priority=5)
    port = fdtd.AddLumpedPort(
        port_nr=1, R=50.0, start=[-1, -1, -5], stop=[1, 1, 5],
        p_dir="z", excite=True, priority=100,
    )
    if epsilon_r != 1.0:
        cube = csx.AddMaterial("Cube", epsilon=epsilon_r)
        cube.AddBox([-150, 300, -50], [-50, 400, 50], priority=10)

    # Identical geometry-aware mesh in both runs.
    grid.AddLine("x", [-600, -150, -50, -2, -1, 0, 1, 2, 200, 600])
    grid.AddLine("y", [-600, -200, 0, 200, 300, 400, 600])
    grid.AddLine("z", [-600, -250, -200, -50, -5, 0, 5, 50, 200, 250, 600])
    for direction, bounds in zip("xyz", ((-150, -50), (300, 400), (-50, 50))):
        grid.AddLine(direction, np.arange(bounds[0], bounds[1] + 0.5 * mesh_mm, mesh_mm))
    grid.SmoothMeshLines("all", min(20.0, 2.0 * mesh_mm), 1.35)

    csx.Write2XML(str(simulation / "m7-mixed.xml"))
    fdtd.Run(str(simulation), cleanup=True, verbose=1)
    port.CalcPort(str(simulation), frequency)
    return (port.uf_tot / port.if_tot)[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-mm", type=float, default=15.0)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--workdir", type=pathlib.Path)
    args = parser.parse_args()
    output = args.output.resolve() if args.output else None
    temporary = None
    root = args.workdir
    if root is None:
        temporary = tempfile.TemporaryDirectory(prefix="nec2c-openems-m7-")
        root = pathlib.Path(temporary.name)
    root.mkdir(parents=True, exist_ok=True)
    vacuum = simulate(root, 1.0, args.mesh_mm)
    loaded = simulate(root, 3.0, args.mesh_mm)
    row = {
        "frequency_mhz": 300.0,
        "vacuum_r_ohm": vacuum.real, "vacuum_x_ohm": vacuum.imag,
        "er3_r_ohm": loaded.real, "er3_x_ohm": loaded.imag,
        "delta_r_ohm": (loaded - vacuum).real,
        "delta_x_ohm": (loaded - vacuum).imag,
        "mesh_mm": args.mesh_mm, "openems_method": "EC-FDTD",
    }
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(row))
            writer.writeheader(); writer.writerow(row)
    else:
        writer = csv.DictWriter(__import__("sys").stdout, fieldnames=list(row))
        writer.writeheader(); writer.writerow(row)
    if temporary is not None:
        temporary.cleanup()


if __name__ == "__main__":
    main()
