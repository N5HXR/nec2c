#!/usr/bin/env python3
"""Generate or run one independently reproducible coated-dipole OpenEMS case.

Run this with the Python interpreter shipped with OpenEMS, for example:
  /home/josh/opt/openEMS/venv/bin/python tools/coated_wire_openems.py --plan
"""

import argparse
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import time


C0 = 299792458.0


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("plan", "setup", "run", "postprocess"),
                   default="plan")
    p.add_argument("--coordinates", choices=("cylindrical", "cartesian"),
                   default="cylindrical")
    p.add_argument("--material", choices=("bare", "coated"), default="bare")
    p.add_argument("--negative-arm-material",
                   choices=("same", "bare", "coated"), default="same",
                   help="override material on the negative-z arm")
    p.add_argument("--feed", choices=("embedded", "split"), default="split")
    p.add_argument("--cylindrical-source",
                   choices=("axis", "ring-voltage"),
                   default="ring-voltage",
                   help="finite-radius source used by cylindrical cases")
    p.add_argument("--mesh-mm", type=float, default=0.5,
                   help="maximum axial pitch near the antenna")
    p.add_argument("--radial-mesh-mm", type=float, default=0.5,
                   help="local transverse/radial pitch")
    p.add_argument("--gap-mm", type=float, default=2.0)
    p.add_argument("--length-mm", type=float, default=475.0)
    p.add_argument("--conductor-radius-mm", type=float, default=2.0)
    p.add_argument("--coating-radius-mm", type=float, default=4.0)
    p.add_argument("--epsilon-r", type=float, default=3.2)
    p.add_argument("--tan-delta", type=float, default=0.0)
    p.add_argument("--negative-arm-coating-radius-mm", type=float,
                   help="negative-arm coating radius; defaults to the common value")
    p.add_argument("--negative-arm-epsilon-r", type=float,
                   help="negative-arm relative permittivity; defaults to the common value")
    p.add_argument("--negative-arm-tan-delta", type=float,
                   help="negative-arm loss tangent; defaults to the common value")
    p.add_argument("--coating-end",
                   choices=("flush", "flat-overhang", "hemisphere",
                            "recessed"),
                   default="flush",
                   help="free-end dielectric geometry; feed faces stay flush")
    p.add_argument("--domain-radius-mm", type=float, default=3200.0)
    p.add_argument("--domain-half-height-mm", type=float, default=3300.0)
    p.add_argument("--nf2ff-radius-mm", type=float, default=500.0)
    p.add_argument("--nf2ff-half-height-mm", type=float, default=600.0)
    p.add_argument("--air-mesh-mm", type=float, default=40.0)
    p.add_argument("--pml-cells", type=int, default=8)
    p.add_argument("--angular-cells", type=int, default=8)
    p.add_argument("--axisymmetric-nf2ff-samples", type=int, default=512,
                   help="postprocess m=0 near fields on this many azimuths; "
                        "zero disables symmetry enforcement")
    p.add_argument("--current-samples", type=int, default=0,
                   help="axial current probes per dipole arm; zero disables")
    p.add_argument("--end-criteria", type=float, default=1.0e-5)
    p.add_argument("--timesteps", type=int, default=720000)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--frequency-min-mhz", type=float, default=250.0)
    p.add_argument("--frequency-max-mhz", type=float, default=350.0)
    p.add_argument("--frequency-count", type=int, default=201)
    p.add_argument("--output-dir", type=Path)
    return p


def validate(a, p):
    negative_radius = (a.negative_arm_coating_radius_mm
                       if a.negative_arm_coating_radius_mm is not None
                       else a.coating_radius_mm)
    negative_er = (a.negative_arm_epsilon_r
                   if a.negative_arm_epsilon_r is not None else a.epsilon_r)
    negative_loss = (a.negative_arm_tan_delta
                     if a.negative_arm_tan_delta is not None else a.tan_delta)
    positive = ("mesh_mm", "radial_mesh_mm", "gap_mm", "length_mm",
                "conductor_radius_mm", "coating_radius_mm", "air_mesh_mm",
                "domain_radius_mm", "domain_half_height_mm",
                "nf2ff_radius_mm", "nf2ff_half_height_mm", "end_criteria")
    for name in positive:
        if getattr(a, name) <= 0:
            p.error(f"--{name.replace('_', '-')} must be positive")
    if a.coating_radius_mm < a.conductor_radius_mm:
        p.error("coating radius must not be smaller than conductor radius")
    if negative_radius < a.conductor_radius_mm:
        p.error("negative-arm coating radius must not be smaller than conductor radius")
    if a.epsilon_r < 1 or a.tan_delta < 0:
        p.error("material must be passive with epsilon-r >= 1")
    if negative_er < 1 or negative_loss < 0:
        p.error("negative-arm material must be passive with epsilon-r >= 1")
    if a.length_mm <= 2 * a.gap_mm:
        p.error("dipole is too short for its feed gap")
    if a.current_samples < 0:
        p.error("current samples must be nonnegative")
    if (a.nf2ff_radius_mm >= a.domain_radius_mm or
            a.nf2ff_half_height_mm >= a.domain_half_height_mm):
        p.error("NF2FF surface must be inside the simulation domain")
    if a.angular_cells < 8 or a.angular_cells % 4:
        p.error("--angular-cells must be a multiple of four and at least eight")
    if (a.axisymmetric_nf2ff_samples and
            (a.axisymmetric_nf2ff_samples < a.angular_cells or
             a.axisymmetric_nf2ff_samples % 4)):
        p.error("--axisymmetric-nf2ff-samples must be zero or a multiple of "
                "four no smaller than --angular-cells")
    if a.material != "coated" and a.coating_end != "flush":
        p.error("--coating-end variants require --material coated")
    if a.negative_arm_material != "same" and a.coating_end != "flush":
        p.error("asymmetric arms currently require --coating-end flush")
    if a.negative_arm_material != "same" and a.feed != "split":
        p.error("asymmetric arms require --feed split")
    if a.frequency_count < 3 or a.frequency_min_mhz >= a.frequency_max_mhz:
        p.error("invalid frequency sweep")
    if a.timesteps < 100 or a.threads < 1 or a.pml_cells < 4:
        p.error("invalid execution limit")


def evenly_spaced(lo, hi, pitch):
    import numpy as np
    count = max(1, int(math.ceil((hi - lo) / pitch)))
    return np.linspace(lo, hi, count + 1)


def axial_current_positions(gap, half, pitch, count):
    """Return symmetric sample planes that stay clear of feed and end cells."""
    import numpy as np
    if not count:
        return np.array([], dtype=float)
    positive = np.linspace(gap/2 + pitch, half - pitch, count)
    # Probe planes are observations, not mesh-refinement requests. Snap them
    # to the existing uniform antenna lattice so a convenient sample count
    # cannot introduce pathological sliver cells or change the FDTD solution.
    positive = np.round((positive + half) / pitch) * pitch - half
    positive = np.unique(positive)
    return np.concatenate((-positive[::-1], positive))


def build(a):
    import numpy as np
    from CSXCAD import ContinuousStructure
    from CSXCAD.CSRectGrid import CoordinateSystem
    from openEMS import openEMS

    cylindrical = a.coordinates == "cylindrical"
    fdtd = openEMS(NrTS=a.timesteps, EndCriteria=a.end_criteria,
                   CoordSystem=1 if cylindrical else 0)
    center_hz = 0.5 * (a.frequency_min_mhz + a.frequency_max_mhz) * 1e6
    half_band_hz = 0.65 * (a.frequency_max_mhz - a.frequency_min_mhz) * 1e6
    fdtd.SetGaussExcite(center_hz, half_band_hz)
    if cylindrical:
        # A complete 2*pi angular mesh plus matching alpha boundary types asks
        # the cylindrical operator for a closed rotation.  PMC here would make
        # the -pi/pi seam a physical symmetry boundary instead.
        pml = f"PML_{a.pml_cells}"
        fdtd.SetBoundaryCond(["PEC", pml, "PEC", "PEC", pml, pml])
        csx = ContinuousStructure(CoordSystem=CoordinateSystem.CYLINDRICAL)
    else:
        fdtd.SetBoundaryCond([f"PML_{a.pml_cells}"] * 6)
        csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)

    a0 = a.conductor_radius_mm
    negative_material = (a.material if a.negative_arm_material == "same"
                         else a.negative_arm_material)
    negative_b0 = (a.negative_arm_coating_radius_mm
                   if a.negative_arm_coating_radius_mm is not None
                   else a.coating_radius_mm)
    negative_er = (a.negative_arm_epsilon_r
                   if a.negative_arm_epsilon_r is not None else a.epsilon_r)
    negative_loss = (a.negative_arm_tan_delta
                     if a.negative_arm_tan_delta is not None else a.tan_delta)
    positive_b0 = a.coating_radius_mm
    b0 = max(a0, positive_b0 if a.material == "coated" else a0,
             negative_b0 if negative_material == "coated" else a0)
    gap = a.gap_mm
    half = a.length_mm / 2

    if cylindrical:
        radial_local = evenly_spaced(0, max(3*b0, 12.0), a.radial_mesh_mm)
        coating_radial = np.array([], dtype=float)
        coating_boundaries = []
        for arm_material, arm_radius in ((a.material, positive_b0),
                                         (negative_material, negative_b0)):
            if arm_material != "coated" or arm_radius <= a0:
                continue
            # As with a PCB substrate, a dielectric layer needs several cells
            # through its thickness.  Treat the requested radial pitch as a
            # maximum and guarantee four coating cells for the reference.
            coating_pitch = min(a.radial_mesh_mm, (arm_radius - a0) / 4)
            coating_boundaries.extend(evenly_spaced(a0, arm_radius,
                                                     coating_pitch))
        if coating_boundaries:
            coating_radial = np.asarray(coating_boundaries)
        # Keep the cylindrical Huygens surface on an actual Yee-grid plane.
        # Otherwise a nominal surface-radius study is confounded by the dump
        # machinery selecting a neighbouring radial plane.
        grid.AddLine("r", np.r_[radial_local, coating_radial,
                                 a.nf2ff_radius_mm,
                                 a.domain_radius_mm])
        # Stay one ulp inside the engine's independently rounded 2*pi limit.
        # The cylindrical operator still recognizes this as a closed rotation.
        alpha_max = math.pi - 1.0e-12
        grid.AddLine("a", np.linspace(-math.pi, alpha_max,
                                      a.angular_cells + 1))
        grid.SmoothMeshLines("r", a.air_mesh_mm, 1.35)
    else:
        local = evenly_spaced(-max(3*b0, 12.0), max(3*b0, 12.0),
                              a.radial_mesh_mm)
        for direction in "xy":
            grid.AddLine(direction, np.r_[-a.domain_radius_mm,
                                           -a.nf2ff_radius_mm, local,
                                           a.nf2ff_radius_mm,
                                           a.domain_radius_mm])
            grid.SmoothMeshLines(direction, a.air_mesh_mm, 1.35)

    # Exact geometry, feed, cap, and Huygens planes plus controlled axial
    # resolution.  The free-end variants extend/recess by one coating radius.
    cap_planes = [-(half + b0), -(half - b0), half - b0, half + b0]
    if a.coating_end == "hemisphere":
        cap_planes.extend(evenly_spaced(half, half + b0, a.mesh_mm))
        cap_planes.extend(evenly_spaced(-half - b0, -half, a.mesh_mm))
    current_planes = axial_current_positions(
        gap, half, a.mesh_mm, a.current_samples)
    axial = np.r_[evenly_spaced(-half, half, a.mesh_mm),
                   -a.nf2ff_half_height_mm,
                   -half, -gap/2, 0, gap/2, half,
                   cap_planes, current_planes,
                   a.nf2ff_half_height_mm,
                   -a.domain_half_height_mm, a.domain_half_height_mm]
    grid.AddLine("z", axial)
    grid.SmoothMeshLines("z", a.air_mesh_mm, 1.35)

    def box(prop, r0, r1, z0, z1, priority):
        if cylindrical:
            prop.AddBox([r0, -math.pi, z0], [r1, math.pi, z1],
                        priority=priority)
        else:
            # Cartesian round geometry is deliberately represented by true
            # cylinders so the refinement study measures staircasing.
            prop.AddCylinder([0, 0, z0], [0, 0, z1], radius=r1,
                             priority=priority)

    split = a.feed == "split"
    intervals = [(-half, -gap/2), (gap/2, half)] if split else [(-half, half)]
    arm_intervals = [
        ("negative", negative_material, negative_b0, negative_er,
         negative_loss, -half, -gap/2),
        ("positive", a.material, positive_b0, a.epsilon_r, a.tan_delta,
         gap/2, half),
    ] if split else [
        ("common", a.material, positive_b0, a.epsilon_r, a.tan_delta,
         -half, half),
    ]
    coatings = []
    for name, arm_material, arm_radius, arm_er, arm_loss, z0, z1 in arm_intervals:
        if arm_material == "coated":
            kappa = (2 * math.pi * center_hz * 8.8541878128e-12 *
                     arm_er * arm_loss)
            coating = csx.AddMaterial(f"coating_{name}", epsilon=arm_er,
                                      kappa=kappa)
            if cylindrical:
                coating.AddBox([a0, -math.pi, z0],
                               [arm_radius, math.pi, z1],
                               priority=5)
            else:
                box(coating, 0, arm_radius, z0, z1, 5)
            coatings.append(coating)

    if a.material == "coated" and a.negative_arm_material == "same":
        coating = coatings[0]
        # Only the two free antenna ends vary in this audit.  The dielectric
        # remains flush at the feed faces so cap effects are not conflated with
        # feed-gap filling.  Metal has higher priority than every cap primitive.
        if a.coating_end == "flat-overhang":
            box(coating, 0, b0, -half - b0, -half, 5)
            box(coating, 0, b0, half, half + b0, 5)
        elif a.coating_end == "hemisphere":
            for sign in (-1, 1):
                edges = evenly_spaced(0, b0, a.mesh_mm)
                for d0, d1 in zip(edges, edges[1:]):
                    # Cell-centre radius gives a convergent axisymmetric
                    # staircase representation of the outward hemisphere.
                    middle = 0.5 * (d0 + d1)
                    radius = math.sqrt(max(0.0, b0*b0 - middle*middle))
                    if sign > 0:
                        z0, z1 = half + d0, half + d1
                    else:
                        z0, z1 = -half - d1, -half - d0
                    box(coating, 0, radius, z0, z1, 5)
        elif a.coating_end == "recessed":
            # Override the annulus close to each free conductor end with air.
            # A distinct background material at higher priority is needed
            # because primitives cannot subtract from a CSX property.
            air_cut = csx.AddMaterial("coating_end_air", epsilon=1.0)
            box(air_cut, a0, b0, -half, -half + b0, 7)
            box(air_cut, a0, b0, half - b0, half, 7)

    metal = csx.AddMetal("conductor")
    for z0, z1 in intervals:
        box(metal, 0, a0, z0, z1, 10)

    if cylindrical:
        # At rho=0 all azimuths denote the same physical z-edge.  Extending a
        # lumped port over rho/alpha duplicates its excitation and leaves the
        # cylindrical current probe with no meaningful enclosing contour.
        # A zero-area transverse box selects the solver's dedicated axis edge.
        port_start = [0, 0, -gap/2]
        port_stop = [0, 0, gap/2]
    else:
        port_start = [0, 0, -gap/2]
        port_stop = [0, 0, gap/2]
    distributed_source = cylindrical and a.cylindrical_source != "axis"
    port = fdtd.AddLumpedPort(1, 50.0, port_start, port_stop, "z",
                              excite=not distributed_source, priority=100)
    if distributed_source:
        alpha = grid.GetLines("a")
        # Axisymmetric delta-gap/frill generator.  Every z-directed edge on
        # the conductor circumference sees the same E field, so the imposed
        # terminal voltage does not grow with angular resolution.  A trial
        # H_alpha drive was rejected: OpenEMS defines it as a magnetic source,
        # not an impressed axial electric current.
        source = csx.AddExcitation("ring_voltage_source", exc_type=0,
                                   exc_val=[0, 0, -1])
        source.SetWeightFunction(["0", "0", f"{1.0/gap:.17g}"])
        source.AddBox([a0, alpha[0], -gap/2],
                      [a0, alpha[-1], gap/2], priority=100)
    if cylindrical:
        # The generic lumped-port current box degenerates on the coordinate
        # axis.  Measure H_alpha just outside the coating and recover total
        # current from the axisymmetric Ampere contour, I=2*pi*r*H_alpha.
        # Keeping this explicit also makes the singular-axis treatment visible
        # in the audit record and gives us a radius-convergence check.
        radial = grid.GetLines("r")
        probe_index = next(index for index, value in enumerate(radial)
                           if value > b0 + 1e-12)
        probe_r = radial[probe_index]
        probe_dual_r = 0.5 * (radial[probe_index - 1] + probe_r)
        probe_z = -gap/2 - a.mesh_mm
        alpha = grid.GetLines("a")
        axis_current_name = "axis_port_it_1"
        current_probe = csx.AddProbe(axis_current_name, p_type=1,
                                     weight=a.angular_cells / 2, norm_dir=2)
        current_probe.AddBox([0, alpha[0], probe_z],
                             [probe_r, alpha[2], probe_z], priority=100)
        axis_h_name = "axis_port_h_1"
        h_probe = csx.AddProbe(axis_h_name, p_type=3)
        h_probe.AddBox([probe_r, 0, probe_z], [probe_r, 0, probe_z],
                       priority=100)
        port.axis_h_filename = axis_h_name
        port.axis_probe_radius_m = probe_dual_r * 1e-3
        port.axis_probe_z_mm = probe_z
        port.I_filenames = [axis_current_name]
        port.axial_current_probes = []
        if len(current_planes):
            for index, position in enumerate(current_planes):
                name = f"axis_current_z_{index:03d}"
                probe = csx.AddProbe(name, p_type=1,
                                     weight=a.angular_cells/2, norm_dir=2)
                probe.AddBox([0, alpha[0], float(position)],
                             [probe_r, alpha[2], float(position)],
                             priority=100)
                port.axial_current_probes.append((name, float(position)))
    else:
        # The stock lumped-port current face snaps asymmetrically when its
        # requested edge lies on a primary Yee line.  Place an explicit square
        # on known dual-grid midplanes so that refinement changes resolution,
        # not the physical measurement contour.
        probe_z = -gap/2 - a.mesh_mm
        probe_lo = -b0 + 0.5 * a.radial_mesh_mm
        probe_hi = b0 + 0.5 * a.radial_mesh_mm
        cart_current_name = "symmetric_port_it_1"
        current_probe = csx.AddProbe(cart_current_name, p_type=1, norm_dir=2)
        current_probe.AddBox([probe_lo, probe_lo, probe_z],
                             [probe_hi, probe_hi, probe_z], priority=100)
        port.I_filenames = [cart_current_name]

    radial_or_x = grid.GetLines(0)
    z_lines = grid.GetLines(2)
    outer_inner_pml = radial_or_x[-a.pml_cells - 1]
    lower_inner_pml = z_lines[a.pml_cells]
    upper_inner_pml = z_lines[-a.pml_cells - 1]
    if (a.nf2ff_radius_mm >= outer_inner_pml or
            -a.nf2ff_half_height_mm <= lower_inner_pml or
            a.nf2ff_half_height_mm >= upper_inner_pml):
        raise ValueError("fixed NF2FF surface overlaps the PML")

    if cylindrical:
        # A closed cylindrical integration surface consists of the outer-rho
        # wall and both z caps.  The axis has zero area and the alpha seam is
        # periodic; treating those disabled faces as PEC mirrors corrupts the
        # radiated-power sign in the generic auto-configuration.
        nf2ff = fdtd.CreateNF2FFBox(
            start=[0, grid.GetLines("a")[0], -a.nf2ff_half_height_mm],
            stop=[a.nf2ff_radius_mm, grid.GetLines("a")[-1],
                  a.nf2ff_half_height_mm],
            directions=[False, True, False, False, True, True],
            mirror=[0] * 6)
    else:
        nf2ff = fdtd.CreateNF2FFBox(
            start=[-a.nf2ff_radius_mm, -a.nf2ff_radius_mm,
                   -a.nf2ff_half_height_mm],
            stop=[a.nf2ff_radius_mm, a.nf2ff_radius_mm,
                  a.nf2ff_half_height_mm])
    return fdtd, csx, port, nf2ff


def mesh_record(csx, pml_cells):
    import numpy as np
    grid = csx.GetGrid()
    lines = [np.asarray(grid.GetLines(i), dtype=float) for i in range(3)]
    counts = [len(v) for v in lines]
    steps = [np.diff(v) for v in lines]
    return {
        "line_counts": counts,
        "cell_counts": [n - 1 for n in counts],
        "total_cells": int(math.prod(n - 1 for n in counts)),
        "minimum_steps": [float(v.min()) for v in steps],
        "maximum_steps": [float(v.max()) for v in steps],
        "maximum_adjacent_step_ratio": [float(max(
            np.max(v[1:] / v[:-1]), np.max(v[:-1] / v[1:])))
            if len(v) > 1 else 1.0 for v in steps],
        "pml_cells": pml_cells,
        "pml_inner_faces": [[float(v[pml_cells]),
                             float(v[-pml_cells - 1])] for v in lines],
    }


def specification(a, mesh):
    record = {
        "schema": "nec2c-coated-wire-openems-case",
        "schema_version": 1,
        "coordinates": a.coordinates,
        "material": a.material,
        "feed": a.feed,
        "cylindrical_source": (a.cylindrical_source
                               if a.coordinates == "cylindrical" else None),
        "geometry_mm": {
            "length": a.length_mm, "feed_gap": a.gap_mm,
            "conductor_radius": a.conductor_radius_mm,
            "coating_outer_radius": a.coating_radius_mm,
            "coating_free_end": a.coating_end,
        },
        "material_parameters": {
            "epsilon_r": 1.0 if a.material == "bare" else a.epsilon_r,
            "tan_delta": 0.0 if a.material == "bare" else a.tan_delta,
        },
        "mesh_controls_mm": {
            "antenna_axial": a.mesh_mm, "local_transverse": a.radial_mesh_mm,
            "air_maximum": a.air_mesh_mm,
            "domain_radius": a.domain_radius_mm,
            "domain_half_height": a.domain_half_height_mm,
            "nf2ff_radius": a.nf2ff_radius_mm,
            "nf2ff_half_height": a.nf2ff_half_height_mm,
            "angular_cells": a.angular_cells if a.coordinates == "cylindrical" else None,
            "axisymmetric_nf2ff_samples":
                (a.axisymmetric_nf2ff_samples
                 if a.coordinates == "cylindrical" else None),
            "current_samples_per_arm": a.current_samples,
            "coating_radial_cells":
                (max(4, int(math.ceil(
                    (a.coating_radius_mm - a.conductor_radius_mm) /
                    a.radial_mesh_mm)))
                 if (a.coordinates == "cylindrical" and
                     a.material == "coated" and
                     a.coating_radius_mm > a.conductor_radius_mm)
                 else None),
        },
        "sweep_mhz": [a.frequency_min_mhz, a.frequency_max_mhz,
                      a.frequency_count],
        "execution": {"timesteps": a.timesteps,
                      "end_criteria": a.end_criteria, "threads": a.threads},
        "mesh": mesh,
    }
    if a.negative_arm_material != "same":
        negative_is_coated = a.negative_arm_material == "coated"
        record["negative_arm"] = {
            "material": a.negative_arm_material,
            "coating_radius_mm": ((a.negative_arm_coating_radius_mm
                                    if a.negative_arm_coating_radius_mm is not None
                                    else a.coating_radius_mm)
                                   if negative_is_coated else None),
            "epsilon_r": ((a.negative_arm_epsilon_r
                            if a.negative_arm_epsilon_r is not None
                            else a.epsilon_r)
                           if negative_is_coated else 1.0),
            "tan_delta": ((a.negative_arm_tan_delta
                            if a.negative_arm_tan_delta is not None
                            else a.tan_delta)
                           if negative_is_coated else 0.0),
        }
    return record


def resonance_from_rows(rows):
    """Linearly interpolate the first reactance zero in the sweep."""
    if rows and max(math.hypot(row["impedance_ohm"]["re"],
                               row["impedance_ohm"]["im"])
                    for row in rows) < 1.0e-12:
        raise ValueError("degenerate zero terminal impedance; the excitation "
                         "is shorted or was not sampled")
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 == 0 or x0 * x1 < 0:
            fraction = 0.0 if x0 == 0 else -x0 / (x1 - x0)
            return (left["frequency_mhz"] + fraction *
                    (right["frequency_mhz"] - left["frequency_mhz"]))
    return None


@contextmanager
def axisymmetric_nf2ff_data(solver, samples):
    """Yield m=0-oversampled NF2FF data and its measured symmetry error.

    OpenEMS' generic cylindrical NF2FF quadrature consumes the physical alpha
    samples directly.  An axisymmetric solve needs few alpha cells for FDTD,
    but many quadrature points for an accurate ring integral.  Average the
    measured cylindrical field components, reject a non-axisymmetric dataset,
    and replicate only temporary HDF5 postprocessing files.
    """
    import h5py
    import numpy as np

    source = Path(solver)
    field_energy = 0.0
    asymmetric_energy = 0.0
    with tempfile.TemporaryDirectory(prefix="coated-wire-nf2ff-") as name:
        target = Path(name)
        alpha = np.linspace(-math.pi, math.pi, samples, endpoint=False)
        for path in sorted(source.glob("nf2ff_[EH]_[012].h5")):
            copied = target / path.name
            shutil.copy2(path, copied)
            with h5py.File(copied, "r+") as handle:
                mesh = handle["Mesh"]
                del mesh["alpha"]
                mesh.create_dataset("alpha", data=alpha)
                time_data = handle["FieldData/TD"]
                for dataset_name in list(time_data):
                    old = time_data[dataset_name]
                    data = old[...]
                    attributes = dict(old.attrs)
                    mean = data.mean(axis=2, keepdims=True)
                    # Accumulate one global energy ratio.  Normalizing each
                    # time slice separately lets roundoff in an almost-zero
                    # field component dominate the maximum even when its
                    # absolute contribution is immaterial.
                    field_energy += float(np.sum(data * data))
                    asymmetric_energy += float(np.sum((data - mean) ** 2))
                    expanded = np.repeat(mean, samples, axis=2)
                    del time_data[dataset_name]
                    new = time_data.create_dataset(dataset_name, data=expanded)
                    for key, value in attributes.items():
                        new.attrs[key] = value
        maximum = (math.sqrt(asymmetric_energy / field_energy)
                   if field_energy else 0.0)
        if maximum > 1.0e-5:
            raise RuntimeError(
                "axisymmetric NF2FF requested for non-axisymmetric fields: "
                f"relative RMS residual {maximum:.6g}")
        yield target, maximum


def main():
    p = parser()
    a = p.parse_args()
    validate(a, p)
    if a.mode == "plan":
        # Imports here also qualify the selected OpenEMS Python environment.
        pass
    try:
        import numpy as np
        import CSXCAD
        import openEMS as openems_package
    except ImportError as error:
        p.error("OpenEMS Python modules unavailable; run with its venv Python: " + str(error))

    fdtd, csx, port, nf2ff = build(a)
    mesh = mesh_record(csx, a.pml_cells)
    spec = specification(a, mesh)
    if hasattr(port, "axis_probe_z_mm"):
        spec["port_measurement"] = {
            "voltage_reference": "centered finite gap",
            "current_contour_z_mm": port.axis_probe_z_mm,
            "deembedding_distance_mm": 0.0,
            "deembedding_note": (
                "The adjacent Yee-cell contour is the lumped terminal-current "
                "measurement, not a transmission-line reference plane."),
        }
    spec["environment"] = {
        "python": platform.python_version(),
        "openems_module": str(Path(openems_package.__file__).resolve()),
        "csxcad_module": str(Path(CSXCAD.__file__).resolve()),
    }
    if a.output_dir is None:
        print(json.dumps(spec, indent=2, sort_keys=True))
        return
    output = a.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if a.mode == "postprocess" and (output / "spec.json").exists():
        previous = json.loads((output / "spec.json").read_text())
        if "wall_seconds" in previous:
            spec["wall_seconds"] = previous["wall_seconds"]
    (output / "spec.json").write_text(json.dumps(spec, indent=2,
                                                  sort_keys=True) + "\n")
    csx.Write2XML(str(output / "geometry.xml"))
    if a.mode == "plan":
        print(output / "spec.json")
        return
    if a.mode != "postprocess":
        original = Path.cwd()
        started = time.monotonic()
        try:
            result = fdtd.Run(str(output / "solver"), cleanup=True,
                              setup_only=a.mode == "setup", verbose=1,
                              numThreads=a.threads, exact_endcriteria=True,
                              dump_statistics=True)
        finally:
            os.chdir(original)
        if result not in (None, 0):
            raise RuntimeError(f"OpenEMS returned setup code {result}")
        spec["wall_seconds"] = time.monotonic() - started
    if a.mode in ("run", "postprocess"):
        if not (output / "solver").is_dir():
            p.error("postprocess requires an existing solver directory")
        frequencies = np.linspace(a.frequency_min_mhz, a.frequency_max_mhz,
                                  a.frequency_count) * 1e6
        port.CalcPort(str(output / "solver"), frequencies)
        impedance = port.uf_tot / port.if_tot
        accepted = port.P_acc
        rows = [{"frequency_mhz": float(f/1e6),
                 "impedance_ohm": {"re": float(z.real), "im": float(z.imag)},
                 "accepted_power_w": float(power)}
                for f, z, power in zip(frequencies, impedance, accepted)]
        (output / "impedance.json").write_text(json.dumps(rows, indent=2) + "\n")
        if np.max(np.abs(impedance)) < 1.0e-12:
            (output / "far-field.json").unlink(missing_ok=True)
            raise RuntimeError("degenerate zero terminal impedance: this feed "
                               "construction is shorted or unobservable")
        resonance_mhz = resonance_from_rows(rows)
        if resonance_mhz is not None:
            if getattr(port, "axial_current_probes", None):
                from openEMS.ports import UI_data
                current_rows = []
                peak = 0.0
                for name, position in port.axial_current_probes:
                    probe = UI_data(name, str(output / "solver"),
                                    resonance_mhz*1e6)
                    value = complex(probe.ui_f_val[0][0])
                    magnitude = abs(value)
                    peak = max(peak, magnitude)
                    current_rows.append({
                        "z_mm": position,
                        "current_a": {"re": value.real, "im": value.imag},
                        "magnitude_a": magnitude,
                    })
                for row in current_rows:
                    row["relative_magnitude"] = (
                        row["magnitude_a"]/peak if peak else 0.0)
                (output / "axial-current.json").write_text(
                    json.dumps({"frequency_mhz": resonance_mhz,
                                "samples": current_rows}, indent=2)+"\n")
            theta = np.linspace(0, 180, 37)
            phi = np.array([0.0, 90.0])
            solver = output / "solver"
            if (a.coordinates == "cylindrical" and
                    a.axisymmetric_nf2ff_samples):
                context = axisymmetric_nf2ff_data(
                    solver, a.axisymmetric_nf2ff_samples)
            else:
                @contextmanager
                def unchanged():
                    yield solver, None
                context = unchanged()
            with context as (far_solver, symmetry_residual):
                far = nf2ff.CalcNF2FF(str(far_solver),
                                      resonance_mhz * 1e6, theta, phi,
                                      outfile="nf2ff_resonance.h5",
                                      read_cached=False)
            norm = np.asarray(far.E_norm[0])
            peak = float(np.max(norm))
            far_rows = []
            for ti, angle in enumerate(theta):
                for pi, azimuth in enumerate(phi):
                    value = float(norm[ti, pi])
                    far_rows.append({
                        "theta_deg": float(angle),
                        "phi_deg": float(azimuth),
                        "relative_db": (-math.inf if value == 0 else
                                        20 * math.log10(value / peak)),
                    })
            directivity = float(far.Dmax[0])
            far_record = {
                "frequency_mhz": resonance_mhz,
                "directivity_db": (None if directivity <= 0 else
                                   float(10 * math.log10(directivity))),
                "radiated_power_w": float(far.Prad[0]),
                "axisymmetric_postprocess": ({
                    "azimuth_samples": a.axisymmetric_nf2ff_samples,
                    "maximum_relative_rms_residual": symmetry_residual,
                } if symmetry_residual is not None else None),
                "samples": far_rows,
            }
            (output / "far-field.json").write_text(
                json.dumps(far_record, indent=2) + "\n")
    (output / "spec.json").write_text(json.dumps(spec, indent=2,
                                                  sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
