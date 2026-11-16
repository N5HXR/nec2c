#!/usr/bin/env python3
"""M4 reduced-kernel antenna convergence and deterministic fill checks."""

import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def execute(text, engine, directory, stem, threads=1, model_check="warn"):
    deck = directory / f"{stem}.nec"
    result = directory / f"{stem}-{engine}-{threads}.json"
    report = directory / f"{stem}-{engine}-{threads}.out"
    deck.write_text(text, encoding="ascii")
    command = [
        str(BINARY), "--wire-engine", engine, "--ground-engine", "nec2",
        "--backend", "builtin", "--jobs", "1", "--threads", str(threads),
        "--model-check", model_check,
        "--result-json", str(result), "-i", str(deck), "-o", str(report),
    ]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((command, completed.returncode,
                              completed.stdout, completed.stderr,
                              report.read_text(errors="replace")[-3000:]
                              if report.exists() else ""))
    return (json.loads(result.read_text(encoding="utf-8")),
            report.read_text(encoding="utf-8", errors="replace"))


def impedance(document):
    value = document["jobs"][0]["executions"][0]["feeds"][0]["impedance_ohm"]
    return complex(value["re"], value["im"])


def scaled(left, right):
    return abs(left-right) / max(1.0, abs(left), abs(right))


def numeric_delta(left, right):
    values = []
    if isinstance(left, dict) and isinstance(right, dict):
        for key in left.keys() & right.keys():
            if key not in {"elapsed_seconds", "threads", "requested", "used"}:
                values.extend(numeric_delta(left[key], right[key]))
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise AssertionError(("result length", len(left), len(right)))
        for a, b in zip(left, right):
            values.extend(numeric_delta(a, b))
    elif (isinstance(left, (int, float)) and not isinstance(left, bool) and
          isinstance(right, (int, float)) and not isinstance(right, bool)):
        values.append(abs(left-right) / max(1.0, abs(left), abs(right)))
    return values


def cards(lines, source_tag, source_segment, frequency=299.792458):
    return ("CM NEC4 SIMILARITY M4 REDUCED KERNEL\nCE\n" +
            "".join(lines) + "GE 0\n" +
            f"FR 0 1 0 0 {frequency:.12g} 0\n" +
            f"EX 0 {source_tag} {source_segment} 0 1 0\n" +
            "RP 0 5 8 1000 0 0 22.5 45 1000 0\nEN\n")


def dipole(count):
    return cards([f"GW 1 {count} 0 0 -0.24 0 0 0.24 0.005\n"],
                 1, (count+1)//2)


def loop(count):
    rows = [
        f"GW 1 {count} -0.1 -0.1 0 0.1 -0.1 0 0.003\n",
        f"GW 2 {count} 0.1 -0.1 0 0.1 0.1 0 0.003\n",
        f"GW 3 {count} 0.1 0.1 0 -0.1 0.1 0 0.003\n",
        f"GW 4 {count} -0.1 0.1 0 -0.1 -0.1 0 0.003\n",
    ]
    return cards(rows, 1, (count+1)//2)


def close_pair(count):
    rows = [
        f"GW 1 {count} -0.015 0 -0.24 -0.015 0 0.24 0.003\n",
        f"GW 2 {count} 0.015 0 -0.24 0.015 0 0.24 0.003\n",
    ]
    return cards(rows, 1, (count+1)//2)


def radius_step(count):
    rows = [
        f"GW 1 {count} 0 0 -0.24 0 0 0 0.002\n",
        f"GW 2 {count} 0 0 0 0 0 0.24 0.004\n",
    ]
    return cards(rows, 1, count)


def main():
    studies = {
        "thick-dipole": (dipole, [11, 21, 31]),
        "square-loop": (loop, [5, 9, 13]),
        "close-two-wire": (close_pair, [11, 21, 31]),
        "radius-step": (radius_step, [6, 10, 14]),
    }
    evidence = {"schema": "nec2c-nec4-m4-acceptance", "schema_version": 1,
                "studies": {}}
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m4-") as name:
        temporary = Path(name)
        for study, (builder, counts) in studies.items():
            values = {engine: [] for engine in ("nec2", "nec4")}
            validity = []
            for count in counts:
                for engine in values:
                    document, report = execute(builder(count), engine, temporary,
                                               f"{study}-{count}-{engine}")
                    value = impedance(document)
                    if not math.isfinite(value.real) or not math.isfinite(value.imag):
                        raise AssertionError((study, count, engine, value))
                    values[engine].append(value)
                    if engine == "nec4":
                        match = re.search(
                            r"NEC4 WIRE VALIDITY: (\d+) SEGMENTS BELOW L/A=2, "
                            r"(\d+) CLOSE PARALLEL PAIRS, (\d+) BURIED MATCH "
                            r"POINTS, (\d+) OVERLAPPING BENDS; MINIMUM L/A ([^\s]+)",
                            report)
                        validity.append(list(map(float, match.groups()))
                                        if match else [0, 0, 0, 0, math.inf])
            errors = {
                engine: [scaled(value, values[engine][-1])
                         for value in values[engine][:-1]]
                for engine in values
            }
            evidence["studies"][study] = {
                "counts": counts,
                "impedance_ohm": {
                    engine: [[value.real, value.imag] for value in sequence]
                    for engine, sequence in values.items()
                },
                "scaled_error_to_finest": errors,
                "validity": validity,
            }
            print(f"{study}: " + ", ".join(
                f"{count}={values['nec4'][index].real:.6f}"
                f"{values['nec4'][index].imag:+.6f}j"
                for index, count in enumerate(counts)))

        # Matrix assembly must preserve the same arithmetic under OpenMP.
        deck = close_pair(31)
        serial, _ = execute(deck, "nec4", temporary, "parallel-identity", 1)
        parallel, _ = execute(deck, "nec4", temporary, "parallel-identity", 4)
        deltas = numeric_delta(serial, parallel)
        maximum_parallel_delta = max(deltas, default=0.0)
        evidence["maximum_serial_parallel_scaled_delta"] = maximum_parallel_delta
        if maximum_parallel_delta > 1.0e-12:
            raise AssertionError(("serial/parallel delta", maximum_parallel_delta))

        # A deliberately over-refined thick wire exercises warn/strict policy.
        _, warning_report = execute(dipole(61), "nec4", temporary,
                                    "validity-warning")
        match = re.search(
            r"NEC4 WIRE VALIDITY: (\d+) SEGMENTS BELOW L/A=2, "
            r"(\d+) CLOSE PARALLEL PAIRS, (\d+) BURIED MATCH POINTS, "
            r"(\d+) OVERLAPPING BENDS; MINIMUM L/A ([^\s]+)",
            warning_report)
        if not match or int(match.group(1)) != 61:
            raise AssertionError(("validity warning", match, warning_report[-2000:]))
        strict_deck = temporary / "validity-strict.nec"
        strict_deck.write_text(dipole(61), encoding="ascii")
        strict = subprocess.run([
            str(BINARY), "--wire-engine", "nec4", "--ground-engine", "nec2",
            "--backend", "builtin", "--jobs", "1", "--threads", "1",
            "--model-check", "strict", "-i", str(strict_deck), "-o",
            str(temporary / "validity-strict.out")], cwd=ROOT, text=True,
            capture_output=True)
        if strict.returncode == 0:
            raise AssertionError("strict validity policy accepted L/A < 2")
        evidence["validity_policy"] = {
            "warn_short_thick_segments": int(match.group(1)),
            "warn_close_parallel_pairs": int(match.group(2)),
            "warn_buried_match_points": int(match.group(3)),
            "warn_overlapping_bends": int(match.group(4)),
            "minimum_length_radius_ratio": float(match.group(5)),
            "strict_rejected": True,
        }

    if os.environ.get("NEC4_M4_RECORD") == "1":
        path = ROOT / "tests" / "nec4-similarity" / "m4-acceptance.json"
        path.write_text(json.dumps(evidence, indent=2, sort_keys=True)+"\n",
                        encoding="utf-8")
    if os.environ.get("NEC4_M4_PLOTS") == "1":
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(2, 2, figsize=(10, 7),
                                     constrained_layout=True)
        for axis, (study, row) in zip(axes.flat, evidence["studies"].items()):
            counts = row["counts"][:-1]
            for engine, marker in (("nec2", "o"), ("nec4", "s")):
                axis.semilogy(counts, row["scaled_error_to_finest"][engine],
                              marker=marker, label=engine)
            axis.set_title(study)
            axis.set_xlabel("segments per wire or side")
            axis.set_ylabel("scaled impedance delta to finest")
            axis.grid(True, which="both", alpha=.3)
            axis.legend()
        figure.savefig(ROOT / "tests" / "nec4-similarity" /
                       "m4-model-convergence.png", dpi=160)
        plt.close(figure)
    print(f"M4 serial/parallel maximum scaled delta "
          f"{maximum_parallel_delta:.3e}")
    print("M4 antenna convergence corpus and deterministic fill passed")


if __name__ == "__main__":
    main()
