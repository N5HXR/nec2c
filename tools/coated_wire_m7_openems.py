#!/usr/bin/env python3
"""Cartesian OpenEMS references for coated-wire M7 general geometries."""

import argparse
import json
import math
import os
from pathlib import Path
import platform
import time

from coated_wire_bend_geometry import bend_case


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article", choices=("straight", "v", "yagi"),
                        required=True)
    parser.add_argument("--material", choices=("bare", "coated"),
                        default="coated")
    parser.add_argument("--mode", choices=("plan", "setup", "run", "postprocess"),
                        default="plan")
    parser.add_argument("--mesh-mm", type=float, default=0.5)
    parser.add_argument("--z-mesh-mm", type=float,
                        help="optional axial pitch for straight/Yagi articles")
    parser.add_argument("--bend-chords", type=int, default=0,
                        help="chords per rounded corner for the v article")
    parser.add_argument("--bend-radius-ratio", type=float, default=0.0,
                        help="bend radius divided by the 4 mm coating radius")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--timesteps", type=int, default=600000)
    parser.add_argument("--end-criteria", type=float, default=1e-4)
    parser.add_argument("--port-current", choices=("symmetric", "stock"),
                        default="symmetric")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--frequency-min-mhz", type=float, default=240.0)
    parser.add_argument("--frequency-max-mhz", type=float, default=360.0)
    parser.add_argument("--frequency-count", type=int, default=241)
    return parser.parse_args()


def span(lo, hi, pitch):
    import numpy as np
    return np.linspace(lo, hi, max(1, math.ceil((hi-lo)/pitch))+1)


def half_offset_span(extent, pitch):
    """Uniform lines at half-pitch offsets, including feed faces +/-p/2."""
    import numpy as np
    count = math.ceil(extent/pitch)
    return (np.arange(-count, count, dtype=float)+0.5)*pitch


def zero_centered_span(lo, hi, pitch):
    import numpy as np
    return np.arange(math.floor(lo/pitch), math.ceil(hi/pitch)+1,
                     dtype=float)*pitch


def smooth_local(grid, direction, centers, radius, pitch, outer, nf2ff):
    import numpy as np
    values = [-outer, -nf2ff, nf2ff, outer]
    for center in centers:
        values.extend(span(center-radius, center+radius, pitch))
    grid.AddLine(direction, np.asarray(values))
    grid.SmoothMeshLines(direction, 30.0, 1.35)


def resonance(rows):
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 == 0 or x0*x1 < 0:
            part = 0.0 if x0 == 0 else -x0/(x1-x0)
            return left["frequency_mhz"] + part*(right["frequency_mhz"]-
                                                   left["frequency_mhz"])
    return None


def build(a):
    import numpy as np
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    center = 0.5*(a.frequency_min_mhz+a.frequency_max_mhz)*1e6
    # Deliberately excite a broader band than the retained output sweep.  The
    # shorter pulse materially reduces Cartesian reference runtime while the
    # complete requested band remains well inside the Gaussian spectrum.
    bandwidth = 1.5*(a.frequency_max_mhz-a.frequency_min_mhz)*1e6
    fdtd = openEMS(NrTS=a.timesteps, EndCriteria=a.end_criteria)
    fdtd.SetGaussExcite(center, bandwidth)
    fdtd.SetBoundaryCond(["PML_8"]*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)

    z_pitch = a.z_mesh_mm if a.z_mesh_mm is not None else a.mesh_mm
    outer, near = 1300.0, 500.0
    if a.article == "v":
        bend = bend_case(a.bend_chords, a.bend_radius_ratio)
        points = [[1000*value for value in point]
                  for point in bend["points_m"]]
        grid.AddLine("x", np.r_[half_offset_span(72, a.mesh_mm), -outer, -near,
                                 near, outer])
        grid.AddLine("z", np.r_[zero_centered_span(-230, 230, z_pitch),
                                 -outer, -near,
                                 near, outer])
        smooth_local(grid, "y", [0], 12, a.mesh_mm, outer, near)
        grid.AddLine("y", [0])
        grid.AddLine("z", [0])
        grid.SmoothMeshLines("x", 30.0, 1.35)
        grid.SmoothMeshLines("z", 30.0, 1.35)
        gap_half = 0.5*a.mesh_mm
        arms = []
        for start, stop in zip(points, points[1:]):
            if (abs(start[2]) < 1e-9 and abs(stop[2]) < 1e-9 and
                    start[0] < 0 < stop[0]):
                arms.extend(((start, [-gap_half, 0, 0]),
                             ([gap_half, 0, 0], stop)))
            else:
                arms.append((start, stop))
        feed_start = [-gap_half, 0, 0]
        feed_stop, feed_direction = [gap_half, 0, 0], "x"
        feed_gap_mm = a.mesh_mm
        current_probes = []
    elif a.article == "straight":
        smooth_local(grid, "x", [0], 12, a.mesh_mm, outer, near)
        smooth_local(grid, "y", [0], 12, a.mesh_mm, outer, near)
        grid.AddLine("z", np.r_[half_offset_span(250, z_pitch),
                                 -outer, -near, near, outer])
        grid.SmoothMeshLines("z", 30.0, 1.35)
        gap_half = 0.5*a.mesh_mm
        half = 237.0
        arms = [([0, 0, -half], [0, 0, -gap_half]),
                ([0, 0, gap_half], [0, 0, half])]
        feed_start = [0, 0, -gap_half]
        feed_stop, feed_direction = [0, 0, gap_half], "z"
        feed_gap_mm = a.mesh_mm
        current_probes = []
    else:
        smooth_local(grid, "x", [-120, 0, 110], 12, a.mesh_mm,
                     outer, near)
        smooth_local(grid, "y", [0], 12, a.mesh_mm, outer, near)
        grid.AddLine("z", np.r_[span(-255, 255, z_pitch), -outer, -near,
                                 near, outer])
        grid.AddLine("x", [-120, 0, 110])
        grid.AddLine("y", [0])
        grid.AddLine("z", [-1, 0, 1])
        grid.SmoothMeshLines("z", 30.0, 1.35)
        arms = [([-120, 0, -245], [-120, 0, 245]),
                ([0, 0, -225], [0, 0, -1]),
                ([0, 0, 1], [0, 0, 225]),
                ([110, 0, -215], [110, 0, 215])]
        feed_start, feed_stop, feed_direction = [0, 0, -1], [0, 0, 1], "z"
        feed_gap_mm = 2.0
        current_probes = []

    coating = (csx.AddMaterial("coating", epsilon=3.2)
               if a.material == "coated" else None)
    metal = csx.AddMetal("conductor")
    for start, stop in arms:
        if coating is not None:
            coating.AddCylinder(start, stop, radius=4.0, priority=5)
        metal.AddCylinder(start, stop, radius=2.0, priority=10)

    port = fdtd.AddLumpedPort(1, 50.0, feed_start, feed_stop,
                              feed_direction, excite=True, priority=100)

    # Stock lumped-port current faces can snap asymmetrically when a terminal
    # lies on a primary Yee line.  Use an explicit, symmetric Ampere surface
    # just behind the negative terminal, as qualified by the M1/M5 harness.
    feed_current_name = "symmetric_port_it_1"
    if a.article == "v":
        feed_current = csx.AddProbe(feed_current_name, p_type=1, norm_dir=0)
        feed_current.AddBox([-1.5*a.mesh_mm, -6, -6],
                            [-1.5*a.mesh_mm, 6, 8], priority=100)
    else:
        feed_current = csx.AddProbe(feed_current_name, p_type=1, norm_dir=2)
        feed_current.AddBox([-6, -6, -1-z_pitch],
                            [6, 6, -1-z_pitch], priority=100)
    if a.port_current == "symmetric":
        port.I_filenames = [feed_current_name]

    if a.article == "yagi":
        # Axis-aligned Ampere surfaces provide representative current samples
        # on all three elements without changing the mesh.
        samples = [(-120, -120), (-120, 0), (-120, 120),
                   (0, -120), (0, -60), (0, 60), (0, 120),
                   (110, -120), (110, 0), (110, 120)]
        for index, (x, z) in enumerate(samples):
            name = f"wire_current_{index:02d}"
            probe = csx.AddProbe(name, p_type=1, norm_dir=2)
            probe.AddBox([x-6, -6, z], [x+6, 6, z], priority=100)
            current_probes.append((name, x, z))

    nf2ff = fdtd.CreateNF2FFBox(start=[-near, -near, -near],
                                stop=[near, near, near])
    return fdtd, csx, port, nf2ff, current_probes


def mesh_record(csx):
    import numpy as np
    lines = [np.asarray(csx.GetGrid().GetLines(i), dtype=float)
             for i in range(3)]
    return {
        "line_counts": [len(item) for item in lines],
        "cell_counts": [len(item)-1 for item in lines],
        "total_cells": math.prod(len(item)-1 for item in lines),
        "minimum_steps_mm": [float(np.min(np.diff(item))) for item in lines],
        "maximum_steps_mm": [float(np.max(np.diff(item))) for item in lines],
    }


def main():
    a = arguments()
    if (a.mesh_mm <= 0 or (a.z_mesh_mm is not None and a.z_mesh_mm <= 0)
            or a.threads < 1):
        raise SystemExit("mesh and thread counts must be positive")
    try:
        import numpy as np
        import CSXCAD
        import openEMS as openems_package
    except ImportError as error:
        raise SystemExit("run with the OpenEMS Python environment: " + str(error))

    fdtd, csx, port, nf2ff, probes = build(a)
    output = a.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    spec = {
        "schema": "nec2c-coated-wire-m7-openems",
        "article": a.article,
        "material": a.material,
        "coordinates": "cartesian",
        "geometry": "true cylinders with flat flush end faces",
        "conductor_radius_mm": 2.0,
        "coating_radius_mm": 4.0,
        "feed_gap_mm": 2.0 if a.article == "yagi" else a.mesh_mm,
        "port_current": a.port_current,
        "epsilon_r": 3.2,
        "mesh_pitch_mm": a.mesh_mm,
        "transverse_mesh_pitch_mm": a.mesh_mm,
        "z_mesh_pitch_mm": (a.z_mesh_mm if a.z_mesh_mm is not None
                            else a.mesh_mm),
        "mesh": mesh_record(csx),
        "sweep_mhz": [a.frequency_min_mhz, a.frequency_max_mhz,
                      a.frequency_count],
        "execution": {"timesteps": a.timesteps,
                      "end_criteria": a.end_criteria,
                      "threads": a.threads},
        "environment": {"python": platform.python_version(),
                        "openems_module": str(Path(openems_package.__file__).resolve()),
                        "csxcad_module": str(Path(CSXCAD.__file__).resolve())},
    }
    if a.article == "v":
        bend_spec = bend_case(a.bend_chords, a.bend_radius_ratio)
        spec["bend"] = {
            "chords_per_corner": a.bend_chords,
            "radius_ratio": a.bend_radius_ratio,
            "maximum_turn_deg": bend_spec["maximum_turn_deg"],
            "minimum_chord_ratio": bend_spec["minimum_chord_ratio"],
            "constant_centerline_length_m": 0.474,
            "points_m": bend_spec["points_m"],
        }
    if a.mode == "postprocess" and (output/"spec.json").exists():
        previous = json.loads((output/"spec.json").read_text())
        if "wall_seconds" in previous:
            spec["wall_seconds"] = previous["wall_seconds"]
    csx.Write2XML(str(output/"geometry.xml"))
    (output/"spec.json").write_text(json.dumps(spec, indent=2)+"\n")
    if a.mode == "plan":
        print(output/"spec.json")
        return
    if a.mode != "postprocess":
        original = Path.cwd()
        started = time.monotonic()
        try:
            result = fdtd.Run(str(output/"solver"), cleanup=True,
                              setup_only=a.mode == "setup", verbose=1,
                              numThreads=a.threads, exact_endcriteria=True,
                              dump_statistics=True)
        finally:
            os.chdir(original)
        if result not in (None, 0):
            raise RuntimeError(f"OpenEMS returned {result}")
        spec["wall_seconds"] = time.monotonic()-started
    if a.mode in ("run", "postprocess"):
        frequencies = np.linspace(a.frequency_min_mhz, a.frequency_max_mhz,
                                  a.frequency_count)*1e6
        port.CalcPort(str(output/"solver"), frequencies)
        impedance = port.uf_tot/port.if_tot
        rows = [{"frequency_mhz": float(f/1e6),
                 "impedance_ohm": {"re": float(z.real), "im": float(z.imag)}}
                for f, z in zip(frequencies, impedance)]
        (output/"impedance.json").write_text(json.dumps(rows, indent=2)+"\n")
        f0 = resonance(rows)
        if f0 is None:
            raise RuntimeError("resonance was not bracketed")
        current_rows = []
        if probes:
            from openEMS.ports import UI_data
            for name, x, z in probes:
                sample = UI_data(name, str(output/"solver"), f0*1e6)
                value = complex(sample.ui_f_val[0][0])
                current_rows.append({"x_mm": x, "z_mm": z,
                                     "current_a": {"re": value.real,
                                                   "im": value.imag},
                                     "magnitude_a": abs(value)})
            peak = max(row["magnitude_a"] for row in current_rows)
            for row in current_rows:
                row["relative_magnitude"] = row["magnitude_a"]/peak
            (output/"wire-current.json").write_text(
                json.dumps({"frequency_mhz": f0,
                            "samples": current_rows}, indent=2)+"\n")
        theta = np.linspace(0, 180, 37)
        phi = np.asarray([0.0, 90.0])
        far = nf2ff.CalcNF2FF(str(output/"solver"), f0*1e6, theta, phi,
                              outfile="nf2ff_resonance.h5", read_cached=False)
        field = np.asarray(far.E_norm[0])
        peak = float(np.max(field))
        pattern = []
        for ti, angle in enumerate(theta):
            for pi, azimuth in enumerate(phi):
                relative = float(field[ti, pi])/peak
                pattern.append({"theta_deg": float(angle),
                                "phi_deg": float(azimuth),
                                "relative_db": (-math.inf if not relative else
                                                20*math.log10(relative))})
        record = {"frequency_mhz": f0,
                  "directivity_dbi": float(10*math.log10(far.Dmax[0])),
                  "radiated_power_w": float(far.Prad[0]),
                  "samples": pattern}
        (output/"far-field.json").write_text(json.dumps(record, indent=2)+"\n")
    (output/"spec.json").write_text(json.dumps(spec, indent=2)+"\n")
    print(output)


if __name__ == "__main__":
    main()
