#!/usr/bin/env python3
"""Generate and measure the NEC side of the M8 geometry battery."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT/"nec2c"


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article", choices=("dogleg45", "staged90",
                                               "groundplane"), required=True)
    parser.add_argument("--material", choices=("bare", "coated"), required=True)
    parser.add_argument("--wire-radius-mm", type=float, default=2.0)
    parser.add_argument("--coating-radius-mm", type=float, default=6.0)
    parser.add_argument("--epsilon-r", type=float, default=3.2)
    parser.add_argument("--frequency-min-mhz", type=float, required=True)
    parser.add_argument("--frequency-max-mhz", type=float, required=True)
    parser.add_argument("--frequency-count", type=int, default=321)
    parser.add_argument("--segment-mm", type=float, default=3.0)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def geometry(name, coating_radius_mm):
    if name == "dogleg45":
        d = 225/math.sqrt(2)
        points = [(-d, 0, -60-d), (0, 0, -60), (0, 0, 60),
                  (d, 0, 60+d)]
        samples = [("center-negative-near", (0, 0, -10)),
                   ("center-negative", (0, 0, -30)),
                   ("center-negative-far", (0, 0, -50)),
                   ("center-positive-near", (0, 0, 10)),
                   ("center-positive", (0, 0, 30)),
                   ("center-positive-far", (0, 0, 50))]
        return points, samples, False
    if name == "staged90":
        radius = 4*coating_radius_mm
        half_span = 110.0
        height = (570-2*half_span-2*radius*(math.sqrt(2)-2))/2
        points = [(-half_span, 0, -height), (-half_span, 0, -radius),
                  (-half_span+radius, 0, 0),
                  (half_span-radius, 0, 0),
                  (half_span, 0, radius), (half_span, 0, height)]
        samples = [("center-negative", (-50, 0, 0)),
                   ("center-positive", (50, 0, 0)),
                   ("outer-negative", (-110, 0, -100)),
                   ("outer-positive", (110, 0, 100))]
        return points, samples, False
    length = 145.0
    projection = length/math.sqrt(2)
    wires = [((0, 0, 0), (0, 0, length)),
             ((0, 0, 0), (projection, 0, -projection)),
             ((0, 0, 0), (-projection, 0, -projection)),
             ((0, 0, 0), (0, projection, -projection)),
             ((0, 0, 0), (0, -projection, -projection))]
    samples = [("monopole-near", (0, 0, 35)),
               ("monopole-middle", (0, 0, 75)),
               ("monopole-far", (0, 0, 115))]
    return wires, samples, True


def wire_rows(a):
    shape, samples, experimental = geometry(a.article, a.coating_radius_mm)
    pairs = (shape if experimental else list(zip(shape, shape[1:])))
    rows = []
    for tag, (left, right) in enumerate(pairs, 1):
        left = tuple(value/1000 for value in left)
        right = tuple(value/1000 for value in right)
        count = max(1, round(math.dist(left, right)/(a.segment_mm/1000)))
        # The central wire in each two-terminal article gets a center segment.
        crosses_origin = any(
            left[axis] < 0 < right[axis] or right[axis] < 0 < left[axis]
            for axis in range(3))
        if crosses_origin and count % 2 == 0:
            count += 1
        rows.append({"tag": tag, "count": count, "left": left,
                     "right": right})
    if experimental:
        source = (1, 1)
    else:
        source_row = next(row for row in rows if any(
            row["left"][axis] < 0 < row["right"][axis] or
            row["right"][axis] < 0 < row["left"][axis]
            for axis in range(3)))
        source = (source_row["tag"], (source_row["count"]+1)//2)
    return rows, samples, experimental, source


def deck(a, frequency=None, pattern=False):
    rows, samples, experimental, source = wire_rows(a)
    lines = [f"CM M8 GEOMETRY {a.article} {a.material}", "CE"]
    radius = a.wire_radius_mm/1000
    for row in rows:
        lines.append("GW %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g" %
                     (row["tag"], row["count"], *row["left"],
                      *row["right"], radius))
    lines.append("GE 0")
    if a.material == "coated":
        coating = a.coating_radius_mm/1000
        lines.extend("LD 7 %d 0 0 %.17g %.17g 0 0" %
                     (row["tag"], a.epsilon_r, coating) for row in rows)
    if frequency is None:
        step = ((a.frequency_max_mhz-a.frequency_min_mhz)/
                (a.frequency_count-1))
        lines.append("FR 0 %d 0 0 %.17g %.17g" %
                     (a.frequency_count, a.frequency_min_mhz, step))
    else:
        lines.append(f"FR 0 1 0 0 {frequency:.17g} 0")
    lines.append(f"EX 0 {source[0]} {source[1]} 0 1 0")
    lines.append("RP 0 37 73 1000 0 0 5 5 0 0" if pattern else "XQ")
    lines.extend(("EN", ""))
    return "\n".join(lines), rows, samples, experimental


def run(a, directory, name, frequency=None, pattern=False):
    text, geometry_rows, samples, experimental = deck(a, frequency, pattern)
    source = directory/f"{name}.nec"
    report = directory/f"{name}.out"
    result = directory/f"{name}.json"
    source.write_text(text, encoding="ascii")
    command = [str(BINARY), "--backend", "auto", "--jobs", "1",
               "--threads", "1", "--result-json", str(result)]
    if experimental:
        command.extend(("--coated-wire-junctions", "experimental"))
    command.extend(("-i", str(source), "-o", str(report)))
    completed = subprocess.run(command, cwd=ROOT, capture_output=True,
                               text=True)
    if completed.returncode:
        raise RuntimeError(completed.stdout+completed.stderr)
    return json.loads(result.read_text()), geometry_rows, samples


def complex_value(item):
    return complex(item["re"], item["im"])


def find_resonance(executions):
    for left, right in zip(executions, executions[1:]):
        zl = complex_value(left["feeds"][0]["impedance_ohm"])
        zr = complex_value(right["feeds"][0]["impedance_ohm"])
        if zl.imag == 0 or zl.imag*zr.imag < 0:
            part = 0 if zl.imag == 0 else -zl.imag/(zr.imag-zl.imag)
            return {"frequency_mhz": left["frequency_mhz"]+part*(
                        right["frequency_mhz"]-left["frequency_mhz"]),
                    "resistance_ohm": zl.real+part*(zr.real-zl.real)}
    raise RuntimeError("resonance not bracketed")


def directivity_and_pattern(samples):
    values = []
    for row in samples:
        et = complex_value(row["r_e_theta"])
        ep = complex_value(row["r_e_phi"])
        values.append(abs(et)**2+abs(ep)**2)
    peak = max(values)
    by_phi = {}
    for row, value in zip(samples, values):
        by_phi.setdefault(row["phi_deg"], []).append((row["theta_deg"], value))
    theta_integrals = []
    for phi, points in sorted(by_phi.items()):
        points.sort()
        integral = sum(math.radians(b[0]-a[0])*0.5*(
            a[1]*math.sin(math.radians(a[0]))+
            b[1]*math.sin(math.radians(b[0]))) for a, b in zip(points, points[1:]))
        theta_integrals.append((phi, integral))
    sphere = sum(math.radians(b[0]-a[0])*0.5*(a[1]+b[1])
                 for a, b in zip(theta_integrals, theta_integrals[1:]))
    pattern = [{"theta_deg": row["theta_deg"], "phi_deg": row["phi_deg"],
                "relative_db": -math.inf if value == 0 else
                10*math.log10(value/peak)}
               for row, value in zip(samples, values)]
    return 10*math.log10(4*math.pi*peak/sphere), pattern


def nearest_currents(execution, wanted):
    result = []
    for name, point_mm in wanted:
        row = min(execution["wire_currents"], key=lambda item: math.dist(
            [1000*x for x in item["center_m"]], point_mm))
        value = complex_value(row["current_a"])
        result.append({"name": name,
                       "point_mm": [1000*x for x in row["center_m"]],
                       "magnitude_a": abs(value)})
    peak = max(row["magnitude_a"] for row in result)
    for row in result:
        row["relative_magnitude"] = row["magnitude_a"]/peak
    return result


def main():
    a = arguments()
    with tempfile.TemporaryDirectory(prefix="nec2c-m8-geometry-") as name:
        directory = Path(name)
        sweep, _, _ = run(a, directory, "sweep")
        resonance = find_resonance(sweep["jobs"][0]["executions"])
        pattern_result, _, wanted = run(a, directory, "pattern",
                                         resonance["frequency_mhz"], True)
    execution = pattern_result["jobs"][0]["executions"][0]
    directivity, pattern = directivity_and_pattern(
        execution["far_fields"][0]["samples"])
    power = execution["power"]
    efficiency = (power["air_radiated_w"]/power["accepted_w"]
                  if power and power["air_radiated_w"] is not None else None)
    record = {"schema": "nec2c-coated-wire-m8-geometry-nec",
              "article": a.article, "material": a.material,
              "wire_radius_mm": a.wire_radius_mm,
              "coating_radius_mm": a.coating_radius_mm,
              "epsilon_r": a.epsilon_r,
              "segment_mm": a.segment_mm,
              "resonance": resonance,
              "directivity_dbi": directivity,
              "radiation_efficiency": efficiency,
              "pattern": pattern,
              "current": nearest_currents(execution, wanted)}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({"resonance": resonance,
                      "directivity_dbi": directivity,
                      "radiation_efficiency": efficiency}, indent=2))


if __name__ == "__main__":
    main()
