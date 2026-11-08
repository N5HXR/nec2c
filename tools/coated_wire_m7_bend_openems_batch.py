#!/usr/bin/env python3
"""Run a bounded OpenEMS bend battery and retain compact resonance records."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT/"tools/coated_wire_m7_openems.py"
RECORD = ROOT/"tests/coated-wire/m7-bend-openems.json"
CASES = {
    "pilot": ((1, 2), (1, 4), (1, 8), (2, 4), (3, 4), (4, 4), (6, 8)),
    "refine": ((1, 4), (2, 4), (4, 4)),
    "coating": ((1, 4), (1, 8), (2, 4), (4, 4)),
}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=tuple(CASES), required=True)
    parser.add_argument("--material", choices=("bare", "coated"), required=True)
    parser.add_argument("--mesh-mm", type=float, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    return parser.parse_args()


def resonance(rows):
    for left, right in zip(rows, rows[1:]):
        x0 = left["impedance_ohm"]["im"]
        x1 = right["impedance_ohm"]["im"]
        if x0 == 0 or x0*x1 < 0:
            part = 0.0 if x0 == 0 else -x0/(x1-x0)
            return {
                "frequency_mhz": left["frequency_mhz"] + part*(
                    right["frequency_mhz"]-left["frequency_mhz"]),
                "resistance_ohm": left["impedance_ohm"]["re"] + part*(
                    right["impedance_ohm"]["re"]-
                    left["impedance_ohm"]["re"]),
            }
    raise RuntimeError("resonance not bracketed")


def main():
    args = arguments()
    record = (json.loads(RECORD.read_text()) if RECORD.exists() else
              {"schema": "nec2c-coated-wire-m7-bend-openems", "runs": []})
    by_key = {row["key"]: row for row in record["runs"]}
    for chords, ratio in CASES[args.stage]:
        key = f"{args.material}-{args.mesh_mm:g}mm-n{chords}-r{ratio}b"
        output = Path("/tmp")/f"m7-bend-{key}"
        impedance_path = output/"impedance.json"
        command = [sys.executable, str(SOLVER), "--article", "v",
                   "--material", args.material, "--mode", "run",
                   "--mesh-mm", str(args.mesh_mm), "--threads",
                   str(args.threads), "--bend-chords", str(chords),
                   "--bend-radius-ratio", str(ratio), "--output-dir",
                   str(output)]
        started = time.monotonic()
        status = "ok"
        if not impedance_path.exists():
            output.mkdir(parents=True, exist_ok=True)
            with (output/"batch.log").open("w") as log:
                try:
                    subprocess.run(command, cwd=ROOT, stdout=log,
                                   stderr=subprocess.STDOUT, check=True,
                                   timeout=args.timeout_seconds)
                except subprocess.TimeoutExpired:
                    status = "timeout"
                except subprocess.CalledProcessError:
                    status = "error"
        if status == "ok":
            spec = json.loads((output/"spec.json").read_text())
            item = {
                "key": key, "status": status, "stage": args.stage,
                "material": args.material, "mesh_mm": args.mesh_mm,
                "chords_per_corner": chords, "radius_ratio": ratio,
                "maximum_turn_deg": spec["bend"]["maximum_turn_deg"],
                "minimum_chord_ratio": spec["bend"]["minimum_chord_ratio"],
                "mesh": spec["mesh"], "solver_wall_seconds":
                    spec.get("wall_seconds"),
                "batch_wall_seconds": time.monotonic()-started,
                "resonance": resonance(json.loads(
                    impedance_path.read_text())),
            }
        else:
            item = {"key": key, "status": status, "stage": args.stage,
                    "material": args.material, "mesh_mm": args.mesh_mm,
                    "chords_per_corner": chords, "radius_ratio": ratio,
                    "batch_wall_seconds": time.monotonic()-started}
        # Reload before every write so independent sequential/resumed stages
        # cannot discard records added by the sharp-convergence collector.
        if RECORD.exists():
            current = json.loads(RECORD.read_text())
            by_key.update({row["key"]: row for row in current["runs"]})
        by_key[key] = item
        record["runs"] = sorted(by_key.values(), key=lambda row: row["key"])
        RECORD.write_text(json.dumps(record, indent=2)+"\n")
        print(key, item.get("resonance", status), flush=True)


if __name__ == "__main__":
    main()
