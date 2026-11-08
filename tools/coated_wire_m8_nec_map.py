#!/usr/bin/env python3
"""Generate the inexpensive NEC side of the M8 thickness validity map."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
C0 = 299792458.0
ETA0 = 376.730313668


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT/"nec2c")
    parser.add_argument("--output", type=Path,
                        default=ROOT/"tests/coated-wire/m8-nec-thickness.json")
    parser.add_argument("--outer-radii-mm", type=float, nargs="+",
                        default=(2, 2.5, 3, 4, 5, 6, 8, 10, 12))
    parser.add_argument("--epsilon-r", type=float, default=3.2)
    parser.add_argument("--tan-delta", type=float, default=0.0)
    parser.add_argument("--conductivity-s-per-m", type=float, default=0.0,
                        help=("frequency-independent dielectric conductivity; "
                              "use this to match an OpenEMS kappa sweep"))
    parser.add_argument("--evaluate-power", action="store_true",
                        help="run a radiation/power evaluation at each resonance")
    parser.add_argument("--threads", default="auto")
    return parser.parse_args()


def deck(radius_mm, epsilon_r, tan_delta, conductivity, *,
         frequency_mhz=180.0, count=401, step_mhz=0.5, pattern=False):
    return "\n".join((
        "CM M8 COATED-WIRE THICKNESS MAP",
        "CE",
        "GW 1 95 0 0 -0.2375 0 0 0.2375 0.002",
        "GE 0",
        f"LD 7 1 0 0 {epsilon_r:.17g} {radius_mm/1000:.17g} "
        f"{tan_delta:.17g} {conductivity:.17g}",
        f"FR 0 {count} 0 0 {frequency_mhz:.17g} {step_mhz:.17g}",
        "EX 0 1 48 0 1 0",
        ("RP 0 37 1 1000 0 0 5 0 0 0" if pattern else "XQ"),
        "EN",
        "",
    ))


def resonance(rows):
    for left, right in zip(rows, rows[1:]):
        zl = left["feeds"][0]["impedance_ohm"]
        zr = right["feeds"][0]["impedance_ohm"]
        x0, x1 = zl["im"], zr["im"]
        if x0 == 0 or x0*x1 < 0:
            part = 0.0 if x0 == 0 else -x0/(x1-x0)
            frequency = left["frequency_mhz"] + part*(
                right["frequency_mhz"]-left["frequency_mhz"])
            resistance = zl["re"] + part*(zr["re"]-zl["re"])
            return frequency, resistance
    raise RuntimeError("resonance is not bracketed")


def run(args, radius, directory):
    name = f"b-{radius:g}mm"
    source, report, result = (directory/f"{name}{suffix}"
                              for suffix in (".nec", ".out", ".json"))
    def execute(text):
        source.write_text(text, encoding="ascii")
        completed = subprocess.run(
            (str(args.binary), "--backend", "auto", "--jobs", "auto",
             "--threads", args.threads, "--result-json", str(result),
             "-i", str(source), "-o", str(report)), cwd=ROOT, text=True,
            capture_output=True)
        if completed.returncode:
            raise RuntimeError(
                f"{name} failed\n{completed.stdout}\n{completed.stderr}")
        return json.loads(result.read_text())

    parsed = execute(deck(radius, args.epsilon_r, args.tan_delta,
                          args.conductivity_s_per_m))
    frequency, resistance = resonance(parsed["jobs"][0]["executions"])
    radius_m = radius/1000
    record = {
        "outer_radius_mm": radius,
        "thickness_mm": radius-2.0,
        "radius_ratio_b_over_a": radius/2.0,
        "resonance_mhz": frequency,
        "resistance_ohm": resistance,
        "k0_b": 2*math.pi*(frequency*1e6/C0)*radius_m,
    }
    if args.evaluate_power:
        evaluated = execute(deck(
            radius, args.epsilon_r, args.tan_delta,
            args.conductivity_s_per_m, frequency_mhz=frequency, count=1,
            step_mhz=0.0, pattern=True))["jobs"][0]["executions"][0]
        power = evaluated["power"]
        field_samples = evaluated["far_fields"][0]["samples"]
        intensities = []
        for sample in field_samples:
            theta = sample["r_e_theta"]
            phi = sample["r_e_phi"]
            intensities.append(
                theta["re"]**2+theta["im"]**2+
                phi["re"]**2+phi["im"]**2)
        peak = max(intensities)
        pattern = [{"theta_deg": sample["theta_deg"],
                    "relative_db": (10*math.log10(value/peak)
                                    if value > 0 else None)}
                   for sample, value in zip(field_samples, intensities)]
        record["resonance_evaluation"] = {
            "impedance_ohm": evaluated["feeds"][0]["impedance_ohm"],
            "power": power,
            "radiation_efficiency": (power["air_radiated_w"] /
                                     power["accepted_w"]),
            "directivity_db": 10*math.log10(
                4*math.pi*peak/(2*ETA0*power["air_radiated_w"])),
            "pattern": pattern,
        }
    return record


def main():
    args = arguments()
    if args.tan_delta < 0 or args.conductivity_s_per_m < 0:
        raise SystemExit("loss parameters must be nonnegative")
    if args.tan_delta and args.conductivity_s_per_m:
        raise SystemExit("select tan delta or conductivity, not both")
    if any(radius < 2 for radius in args.outer_radii_mm):
        raise SystemExit("outer radii must be at least the 2 mm conductor radius")
    # Always retain the b=a null limit; it is the only valid zero-thickness
    # baseline when a caller requests a sparse stress grid.
    radii = sorted(set((2.0, *args.outer_radii_mm)))
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m8-") as temporary:
        cases = [run(args, radius, Path(temporary)) for radius in radii]
    bare = next(case["resonance_mhz"] for case in cases
                if case["outer_radius_mm"] == 2.0)
    for case in cases:
        case["resonance_shift_from_zero_thickness_percent"] = 100*(
            case["resonance_mhz"]/bare-1)
        print(f"b={case['outer_radius_mm']:g} mm: "
              f"{case['resonance_mhz']:.6f} MHz / "
              f"{case['resistance_ohm']:.6f} ohm, "
              f"shift {case['resonance_shift_from_zero_thickness_percent']:.3f}%")
    record = {
        "schema": "nec2c-coated-wire-m8-nec-thickness",
        "article": {"length_mm": 475, "conductor_radius_mm": 2,
                    "segments": 95, "epsilon_r": args.epsilon_r,
                    "tan_delta": args.tan_delta,
                    "conductivity_s_per_m": args.conductivity_s_per_m,
                    "free_end_geometry": "reduced-operator flush convention"},
        "cases": cases,
        "note": ("This is the reduced model's prediction, not independent "
                 "accuracy evidence. OpenEMS comparison supplies the error curve."),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    print(args.output)


if __name__ == "__main__":
    main()
