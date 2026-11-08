#!/usr/bin/env python3
"""Generate the M9 NEC-only VSWR and radiation-pattern acceptance bundle."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def run(binary, source, directory, name):
    result = directory/f"{name}.json"
    report = directory/f"{name}.out"
    command = [str(binary), "--backend", "auto", "--jobs", "1",
               "--threads", "auto", "--result-json", str(result),
               "-i", str(source), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stdout+completed.stderr)
    return json.loads(result.read_text())


def vswr(z, reference=50.0):
    impedance = complex(z["re"], z["im"])
    reflection = abs((impedance-reference)/(impedance+reference))
    return (1+reflection)/(1-reflection) if reflection < 1 else math.inf


def material_summary(execution):
    rows = execution["material_evaluations"]
    if not rows:
        raise AssertionError("coated acceptance execution omitted material metadata")
    return {
        "maximum_validity_metric": max(
            row["validity_metric_sqrt_er_k0_b"] for row in rows),
        "statuses": sorted({row["validity_status"] for row in rows}),
    }


def vswr_summary(result):
    rows = []
    maximum_q = 0.0
    statuses = set()
    for execution in result["jobs"][0]["executions"]:
        feed = execution["feeds"][0]
        impedance = feed["impedance_ohm"]
        rows.append({"frequency_mhz": execution["frequency_mhz"],
                     "resistance_ohm": impedance["re"],
                     "reactance_ohm": impedance["im"],
                     "vswr_50_ohm": vswr(impedance)})
        material = material_summary(execution)
        maximum_q = max(maximum_q, material["maximum_validity_metric"])
        statuses.update(material["statuses"])
    best = min(rows, key=lambda row: row["vswr_50_ohm"])
    return {"reference_impedance_ohm": 50.0, "samples": rows,
            "minimum": best, "maximum_validity_metric": maximum_q,
            "validity_statuses": sorted(statuses)}


def pattern_summary(result):
    execution = result["jobs"][0]["executions"][0]
    samples = execution["far_fields"][0]["samples"]
    values = []
    peak = 0.0
    for sample in samples:
        et = sample["r_e_theta"]
        ep = sample["r_e_phi"]
        magnitude = math.sqrt(et["re"]**2+et["im"]**2+
                              ep["re"]**2+ep["im"]**2)
        peak = max(peak, magnitude)
        values.append((sample, magnitude))
    rows = [{"theta_deg": sample["theta_deg"],
             "phi_deg": sample["phi_deg"],
             "relative_db": (20*math.log10(magnitude/peak)
                             if magnitude and peak else -300.0)}
            for sample, magnitude in values]
    cut = {row["theta_deg"]: row["relative_db"] for row in rows
           if row["phi_deg"] == 0}
    symmetry = max(abs(value-cut[180-theta]) for theta, value in cut.items()
                   if value >= -40 and cut.get(180-theta, -300) >= -40)
    peak_row = max(rows, key=lambda row: row["relative_db"])
    material = material_summary(execution)
    return {"frequency_mhz": execution["frequency_mhz"],
            "samples": rows,
            "peak_direction": {"theta_deg": peak_row["theta_deg"],
                               "phi_deg": peak_row["phi_deg"]},
            "symmetry_maximum_db_above_minus_40": symmetry,
            **material}


def plot_vswr(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axis = plt.subplots(figsize=(7.5, 4.5))
    axis.plot([row["frequency_mhz"] for row in summary["samples"]],
              [row["vswr_50_ohm"] for row in summary["samples"]])
    best = summary["minimum"]
    axis.scatter([best["frequency_mhz"]], [best["vswr_50_ohm"]], color="red")
    axis.set(xlabel="Frequency (MHz)", ylabel="VSWR (50 ohm)",
             title="Lossy coated dipole")
    axis.set_ylim(1, min(10, max(row["vswr_50_ohm"]
                                for row in summary["samples"])))
    axis.grid(True, alpha=.3)
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_pattern(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axis = plt.subplots(subplot_kw={"projection": "polar"},
                                figsize=(6.5, 6.5))
    rows = sorted((row for row in summary["samples"] if row["phi_deg"] == 0),
                  key=lambda row: row["theta_deg"])
    angles = [math.radians(row["theta_deg"]) for row in rows]
    levels = [max(-30, row["relative_db"]) for row in rows]
    axis.plot(angles+[2*math.pi-angle for angle in reversed(angles)],
              levels+list(reversed(levels)), label="phi=0 deg")
    axis.set_ylim(-30, 0)
    axis.set_theta_zero_location("E")
    axis.set_theta_direction(-1)
    axis.set_title("Lossy coated dipole, relative field")
    axis.grid(True, alpha=.3)
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT/"nec2c")
    parser.add_argument("--vswr-deck", type=Path,
                        default=ROOT/"tests/coated-wire/m9-vswr.nec")
    parser.add_argument("--pattern-deck", type=Path,
                        default=ROOT/"tests/coated-wire/m9-pattern.nec")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure-directory", type=Path)
    arguments = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m9-acceptance-") as name:
        directory = Path(name)
        result = {
            "schema": "nec2c-coated-wire-m9-acceptance",
            "vswr": vswr_summary(run(arguments.binary.resolve(),
                                      arguments.vswr_deck.resolve(), directory,
                                      "vswr")),
            "pattern": pattern_summary(run(arguments.binary.resolve(),
                                            arguments.pattern_deck.resolve(),
                                            directory, "pattern")),
        }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2)+"\n")
    if arguments.figure_directory:
        arguments.figure_directory.mkdir(parents=True, exist_ok=True)
        plot_vswr(result["vswr"], arguments.figure_directory/"m9-vswr.png")
        plot_pattern(result["pattern"], arguments.figure_directory/"m9-pattern.png")
    print(f"VSWR minimum {result['vswr']['minimum']['vswr_50_ohm']:.3f} at "
          f"{result['vswr']['minimum']['frequency_mhz']:.3f} MHz")
    print("Dipole pattern symmetry defect "
          f"{result['pattern']['symmetry_maximum_db_above_minus_40']:.3e} dB")
    print(arguments.output)


if __name__ == "__main__":
    main()
