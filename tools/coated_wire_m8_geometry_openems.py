#!/usr/bin/env python3
"""Cartesian OpenEMS references for the M8 thick-coating geometry battery."""

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import platform
import time


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article", choices=("dogleg45", "staged90",
                                               "groundplane"), required=True)
    parser.add_argument("--material", choices=("bare", "coated"), required=True)
    parser.add_argument("--mode", choices=("plan", "setup", "run",
                                            "postprocess"), default="plan")
    parser.add_argument("--wire-radius-mm", type=float, default=2.0)
    parser.add_argument("--coating-radius-mm", type=float, default=6.0)
    parser.add_argument("--epsilon-r", type=float, default=3.2)
    parser.add_argument("--mesh-mm", type=float, default=1.5)
    parser.add_argument("--feed-gap-mm", type=float, default=2.0)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--timesteps", type=int, default=600000)
    parser.add_argument("--end-criteria", type=float, default=1e-4)
    parser.add_argument("--frequency-min-mhz", type=float, required=True)
    parser.add_argument("--frequency-max-mhz", type=float, required=True)
    parser.add_argument("--frequency-count", type=int, default=301)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def line(a, b):
    return [list(a), list(b)]


def article_geometry(name, coating_radius, gap):
    """Return millimetre centerlines, feed terminals, and validation facts."""
    if name == "dogleg45":
        center_half = 60.0
        outer_length = 225.0
        delta = outer_length/math.sqrt(2)
        points = [(-delta, 0, -center_half-delta),
                  (0, 0, -center_half), (0, 0, center_half),
                  (delta, 0, center_half+delta)]
        conductors = [line(points[0], points[1]),
                      line(points[1], (0, 0, -gap/2)),
                      line((0, 0, gap/2), points[2]),
                      line(points[2], points[3])]
        return {
            "points_mm": points, "conductors_mm": conductors,
            "feed_start_mm": [0, 0, -gap/2],
            "feed_stop_mm": [0, 0, gap/2], "feed_direction": "z",
            "maximum_turn_deg": 45.0,
            "minimum_chord_outer_radius_ratio": min(
                outer_length, 2*center_half)/coating_radius,
            "experimental_topology": False,
            "current_samples": [
                {"name": "center-negative-near", "point_mm": [0, 0, -10],
                 "direction": "z"},
                {"name": "center-negative", "point_mm": [0, 0, -30],
                 "direction": "z"},
                {"name": "center-negative-far", "point_mm": [0, 0, -50],
                 "direction": "z"},
                {"name": "center-positive-near", "point_mm": [0, 0, 10],
                 "direction": "z"},
                {"name": "center-positive", "point_mm": [0, 0, 30],
                 "direction": "z"},
                {"name": "center-positive-far", "point_mm": [0, 0, 50],
                 "direction": "z"},
            ],
        }
    if name == "staged90":
        radius = 4*coating_radius
        half_span = 110.0
        # Preserve about 570 mm total centerline length as b changes.
        height = (570-2*half_span-2*radius*(math.sqrt(2)-2))/2
        points = [(-half_span, 0, -height),
                  (-half_span, 0, -radius),
                  (-half_span+radius, 0, 0),
                  (half_span-radius, 0, 0),
                  (half_span, 0, radius),
                  (half_span, 0, height)]
        conductors = []
        for left, right in zip(points, points[1:]):
            if left[2] == 0 and right[2] == 0 and left[0] < 0 < right[0]:
                conductors.extend((line(left, (-gap/2, 0, 0)),
                                   line((gap/2, 0, 0), right)))
            else:
                conductors.append(line(left, right))
        chord_lengths = [math.dist(left, right)
                         for left, right in zip(points, points[1:])]
        return {
            "points_mm": points, "conductors_mm": conductors,
            "feed_start_mm": [-gap/2, 0, 0],
            "feed_stop_mm": [gap/2, 0, 0], "feed_direction": "x",
            "maximum_turn_deg": 45.0,
            "minimum_chord_outer_radius_ratio": (
                min(chord_lengths)/coating_radius),
            "experimental_topology": False,
            "current_samples": [
                {"name": "center-negative", "point_mm": [-50, 0, 0],
                 "direction": "x"},
                {"name": "center-positive", "point_mm": [50, 0, 0],
                 "direction": "x"},
                {"name": "outer-negative", "point_mm": [-110, 0, -100],
                 "direction": "z"},
                {"name": "outer-positive", "point_mm": [110, 0, 100],
                 "direction": "z"},
            ],
        }
    length = 145.0
    projection = length/math.sqrt(2)
    conductors = [line((0, 0, gap), (0, 0, gap+length)),
                  line((0, 0, 0), (projection, 0, -projection)),
                  line((0, 0, 0), (-projection, 0, -projection)),
                  line((0, 0, 0), (0, projection, -projection)),
                  line((0, 0, 0), (0, -projection, -projection))]
    return {
        "points_mm": None, "conductors_mm": conductors,
        "feed_start_mm": [0, 0, 0], "feed_stop_mm": [0, 0, gap],
        "feed_direction": "z", "maximum_turn_deg": 90.0,
        "minimum_chord_outer_radius_ratio": length/coating_radius,
        "experimental_topology": True,
        "radial_count": 4,
        "radial_slope_deg": -45.0,
        "current_samples": [
            {"name": "monopole-near", "point_mm": [0, 0, gap+35],
             "direction": "z"},
            {"name": "monopole-middle", "point_mm": [0, 0, gap+75],
             "direction": "z"},
            {"name": "monopole-far", "point_mm": [0, 0, gap+115],
             "direction": "z"},
        ],
    }


def uniform_lines(low, high, pitch, anchor=0.0):
    import numpy as np
    first = math.floor((low-anchor)/pitch)
    last = math.ceil((high-anchor)/pitch)
    return anchor+np.arange(first, last+1, dtype=float)*pitch


def build(a):
    import numpy as np
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    geometry = article_geometry(a.article, a.coating_radius_mm,
                                a.feed_gap_mm)
    center = 0.5*(a.frequency_min_mhz+a.frequency_max_mhz)*1e6
    bandwidth = 1.5*(a.frequency_max_mhz-a.frequency_min_mhz)*1e6
    fdtd = openEMS(NrTS=a.timesteps, EndCriteria=a.end_criteria)
    fdtd.SetGaussExcite(center, bandwidth)
    fdtd.SetBoundaryCond(["PML_8"]*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)

    values = [coordinate for segment in geometry["conductors_mm"]
              for point in segment for coordinate in (point,)]
    bounds = [[min(point[axis] for point in values),
               max(point[axis] for point in values)] for axis in range(3)]
    pad = max(3*a.coating_radius_mm, 6*a.mesh_mm)
    near = 400.0 if a.article != "groundplane" else 300.0
    outer = 1200.0 if a.article != "groundplane" else 900.0
    feed_axis = "xyz".index(geometry["feed_direction"])
    for axis, direction in enumerate("xyz"):
        low, high = bounds[axis]
        anchor = geometry["feed_start_mm"][axis] if axis == feed_axis else 0
        local = uniform_lines(low-pad, high+pad, a.mesh_mm, anchor)
        grid.AddLine(direction, np.r_[local, -outer, -near, near, outer])
        grid.SmoothMeshLines(direction, 30.0, 1.35)
    for point in (geometry["feed_start_mm"], geometry["feed_stop_mm"]):
        for axis, direction in enumerate("xyz"):
            grid.AddLine(direction, point[axis])

    coating = (csx.AddMaterial("coating", epsilon=a.epsilon_r)
               if a.material == "coated" else None)
    metal = csx.AddMetal("conductor")
    for start, stop in geometry["conductors_mm"]:
        if coating is not None:
            coating.AddCylinder(start, stop, radius=a.coating_radius_mm,
                                priority=5)
        metal.AddCylinder(start, stop, radius=a.wire_radius_mm, priority=10)

    port_start = list(geometry["feed_start_mm"])
    port_stop = list(geometry["feed_stop_mm"])
    if a.article == "groundplane":
        feed_axis = "xyz".index(geometry["feed_direction"])
        for axis in range(3):
            if axis != feed_axis:
                port_start[axis] -= a.wire_radius_mm
                port_stop[axis] += a.wire_radius_mm
    port = fdtd.AddLumpedPort(1, 50.0, port_start, port_stop,
                              geometry["feed_direction"], excite=True,
                              priority=100)
    # Use a symmetric Ampere surface one cell behind the negative terminal.
    extent = a.coating_radius_mm+2*a.mesh_mm
    current_name = "symmetric_port_it_1"
    direction_index = "xyz".index(geometry["feed_direction"])
    probe = csx.AddProbe(current_name, p_type=1, norm_dir=direction_index)
    start = list(geometry["feed_start_mm"])
    stop = list(geometry["feed_start_mm"])
    plane = (geometry["feed_stop_mm"][direction_index]+a.mesh_mm
             if a.article == "groundplane" else
             geometry["feed_start_mm"][direction_index]-a.mesh_mm)
    start[direction_index] = stop[direction_index] = plane
    for axis in range(3):
        if axis != direction_index:
            start[axis] -= extent
            stop[axis] += extent
    probe.AddBox(start, stop, priority=100)
    # The ground-plane article samples the positive element one cell beyond
    # the feed; two-arm articles use the qualified loop just behind the
    # negative terminal.  Both avoid the stock two-face cancellation.
    port.I_filenames = [current_name]

    current_probes = []
    for index, sample in enumerate(geometry["current_samples"]):
        name = f"wire_current_{index:02d}"
        direction_index = "xyz".index(sample["direction"])
        probe = csx.AddProbe(name, p_type=1, norm_dir=direction_index)
        start = list(sample["point_mm"])
        stop = list(sample["point_mm"])
        for axis in range(3):
            if axis != direction_index:
                start[axis] -= extent
                stop[axis] += extent
        probe.AddBox(start, stop, priority=100)
        current_probes.append((name, sample))

    nf2ff = fdtd.CreateNF2FFBox(start=[-near]*3, stop=[near]*3)
    return fdtd, csx, port, nf2ff, geometry, current_probes


def mesh_record(csx):
    import numpy as np
    lines = [np.asarray(csx.GetGrid().GetLines(axis), dtype=float)
             for axis in range(3)]
    return {
        "line_counts": [len(item) for item in lines],
        "cell_counts": [len(item)-1 for item in lines],
        "total_cells": math.prod(len(item)-1 for item in lines),
        "minimum_steps_mm": [float(np.min(np.diff(item))) for item in lines],
        "maximum_steps_mm": [float(np.max(np.diff(item))) for item in lines],
    }


def resonance(rows):
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 == 0 or x0*x1 < 0:
            part = 0 if x0 == 0 else -x0/(x1-x0)
            return {
                "frequency_mhz": left["frequency_mhz"]+part*(
                    right["frequency_mhz"]-left["frequency_mhz"]),
                "resistance_ohm": left["impedance_ohm"]["re"]+part*(
                    right["impedance_ohm"]["re"]-
                    left["impedance_ohm"]["re"]),
                "accepted_power_w": left["accepted_power_w"]+part*(
                    right["accepted_power_w"]-left["accepted_power_w"]),
            }
    raise RuntimeError("resonance not bracketed")


def main():
    a = arguments()
    if (a.wire_radius_mm <= 0 or
            a.coating_radius_mm < a.wire_radius_mm or
            a.epsilon_r < 1 or a.mesh_mm <= 0 or a.feed_gap_mm <= 0):
        raise SystemExit("invalid geometry, material, or mesh")
    try:
        import numpy as np
        import CSXCAD
        import openEMS as openems_package
    except ImportError as error:
        raise SystemExit("run with the OpenEMS Python environment: "+str(error))

    fdtd, csx, port, nf2ff, geometry, current_probes = build(a)
    output = a.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    spec = {
        "schema": "nec2c-coated-wire-m8-geometry-openems",
        "article": a.article, "material": a.material,
        "coordinates": "cartesian",
        "wire_radius_mm": a.wire_radius_mm,
        "coating_radius_mm": a.coating_radius_mm,
        "epsilon_r": a.epsilon_r, "feed_gap_mm": a.feed_gap_mm,
        "mesh_pitch_mm": a.mesh_mm, "geometry": geometry,
        "mesh": mesh_record(csx),
        "sweep_mhz": [a.frequency_min_mhz, a.frequency_max_mhz,
                      a.frequency_count],
        "execution": {"threads": a.threads, "timesteps": a.timesteps,
                      "end_criteria": a.end_criteria},
        "environment": {"python": platform.python_version(),
                        "openems_module": str(Path(
                            openems_package.__file__).resolve()),
                        "csxcad_module": str(Path(CSXCAD.__file__).resolve())},
    }
    if a.mode == "postprocess" and (output/"spec.json").exists():
        old = json.loads((output/"spec.json").read_text())
        if "wall_seconds" in old:
            spec["wall_seconds"] = old["wall_seconds"]
    csx.Write2XML(str(output/"geometry.xml"))
    (output/"spec.json").write_text(json.dumps(spec, indent=2)+"\n")
    if a.mode == "plan":
        print(json.dumps({"output": str(output), "mesh": spec["mesh"],
                          "geometry": geometry}, indent=2))
        return
    if a.mode != "postprocess":
        original = Path.cwd()
        started = time.monotonic()
        try:
            status = fdtd.Run(str(output/"solver"), cleanup=True,
                              setup_only=a.mode == "setup", verbose=1,
                              numThreads=a.threads, exact_endcriteria=True,
                              dump_statistics=True)
        finally:
            os.chdir(original)
        if status not in (None, 0):
            raise RuntimeError(f"OpenEMS returned {status}")
        spec["wall_seconds"] = time.monotonic()-started
    if a.mode in ("run", "postprocess"):
        frequencies = np.linspace(a.frequency_min_mhz,
                                  a.frequency_max_mhz,
                                  a.frequency_count)*1e6
        port.CalcPort(str(output/"solver"), frequencies)
        impedance = port.uf_tot/port.if_tot
        rows = [{"frequency_mhz": float(f/1e6),
                 "impedance_ohm": {"re": float(z.real),
                                   "im": float(z.imag)},
                 "accepted_power_w": float(power)}
                for f, z, power in zip(frequencies, impedance, port.P_acc)]
        (output/"impedance.json").write_text(json.dumps(rows, indent=2)+"\n")
        resonant = resonance(rows)
        from openEMS.ports import UI_data
        current_rows = []
        for name, sample in current_probes:
            data = UI_data(name, str(output/"solver"),
                           resonant["frequency_mhz"]*1e6)
            value = complex(data.ui_f_val[0][0])
            current_rows.append({**sample,
                                 "current_a": {"re": value.real,
                                               "im": value.imag},
                                 "magnitude_a": abs(value)})
        peak_current = max(row["magnitude_a"] for row in current_rows)
        for row in current_rows:
            row["relative_magnitude"] = row["magnitude_a"]/peak_current
        (output/"wire-current.json").write_text(json.dumps({
            "frequency_mhz": resonant["frequency_mhz"],
            "samples": current_rows,
        }, indent=2)+"\n")
        theta = np.linspace(0, 180, 37)
        phi = np.linspace(0, 360, 73)
        cache = output/("nf2ff_%.6fMHz.h5" %
                        resonant["frequency_mhz"])
        far = nf2ff.CalcNF2FF(str(output/"solver"),
                              resonant["frequency_mhz"]*1e6, theta, phi,
                              outfile=str(cache),
                              read_cached=False)
        field = np.asarray(far.E_norm[0])
        peak = float(np.max(field))
        samples = []
        for ti, angle in enumerate(theta):
            for pi, azimuth in enumerate(phi):
                relative = float(field[ti, pi])/peak
                samples.append({"theta_deg": float(angle),
                                "phi_deg": float(azimuth),
                                "relative_db": (-math.inf if not relative else
                                                20*math.log10(relative))})
        far_record = {
            "frequency_mhz": resonant["frequency_mhz"],
            "resistance_ohm": resonant["resistance_ohm"],
            "directivity_dbi": float(10*math.log10(far.Dmax[0])),
            "radiated_power_w": float(far.Prad[0]),
            "accepted_power_w": resonant["accepted_power_w"],
            "radiation_efficiency": (float(far.Prad[0]) /
                                     resonant["accepted_power_w"]),
            "samples": samples,
        }
        (output/"far-field.json").write_text(
            json.dumps(far_record, indent=2)+"\n")
    (output/"spec.json").write_text(json.dumps(spec, indent=2)+"\n")
    print(output)


if __name__ == "__main__":
    main()
