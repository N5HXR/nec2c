#!/usr/bin/env python3
"""Compare feed impedance across nec2c, momwire, and an external AEG build."""

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def run(command, cwd=None):
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(json.dumps({"command": command, "returncode": completed.returncode,
                                       "stdout": completed.stdout,
                                       "stderr": completed.stderr}, indent=2))
    return completed


def impedance_from_nec2c(binary, deck, directory, backend):
    result = directory / "nec2c-result.json"
    run([str(binary), "--backend", backend, "--result-json", str(result), "-i", str(deck),
         "-o", str(directory / "nec2c.out")], cwd=ROOT)
    value = json.loads(result.read_text(encoding="utf-8"))
    feed = value["jobs"][0]["executions"][0]["feeds"][0]
    impedance = feed["impedance_ohm"]
    return complex(impedance["re"], impedance["im"])


def impedance_from_aeg(binary, deck, directory):
    copied = directory / "aeg-dipole.nec"
    copied.write_bytes(deck.read_bytes())
    run([str(binary), "--segments=1000", str(copied)], cwd=directory)
    output = directory / "aeg-dipole.res"
    text = output.read_text(encoding="ascii", errors="replace")
    heading = text.find("ANTENNA INPUT PARAMETERS")
    if heading < 0:
        raise RuntimeError("AEG output lacks ANTENNA INPUT PARAMETERS")
    rows = text[heading:].splitlines()
    for line in rows:
        numbers = re.findall(FLOAT, line)
        if len(numbers) >= 11 and numbers[0].isdigit() and numbers[1].isdigit():
            return complex(float(numbers[6]), float(numbers[7]))
    raise RuntimeError("could not parse AEG feed impedance")


def impedance_from_momwire(root, deck):
    python = root / ".venv/bin/python"
    program = """
import json
import sys
from pathlib import Path
from momwire.deck import build_solver, parse
model = parse(Path(sys.argv[1]).read_text())
built = build_solver(model, basis="sinusoidal")
solution = built.solver.compute_port_solution()
port = built.ports.feed_ports[0]
z = 1.0 / solution.y[port, port]
print(json.dumps([float(z.real), float(z.imag)]))
"""
    completed = run([str(python), "-c", program, str(deck)], cwd=root)
    real, imaginary = json.loads(completed.stdout.strip().splitlines()[-1])
    return complex(real, imaginary)


def delta(reference, observed):
    absolute = abs(observed - reference)
    relative = absolute / max(abs(reference), sys.float_info.min)
    return {"absolute_ohms": absolute, "relative": relative}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nec2c", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--nec2c-backend", choices=("auto", "builtin", "lapack"),
                        default="auto")
    parser.add_argument("--momwire-root", type=Path, default=Path.home() / "src/momwire")
    parser.add_argument("--aegnec2", type=Path, required=True,
                        help="external AEG executable; never imported into this repository")
    parser.add_argument("--deck", type=Path, action="append", dest="decks",
                        help="compatible comparison deck; may be repeated")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--aeg-relative-tolerance", type=float, default=1e-3)
    parser.add_argument("--momwire-relative-tolerance", type=float, default=5e-3)
    arguments = parser.parse_args()

    decks = arguments.decks or [
        ROOT / "tests/performance/cross-solver-dipole.nec",
        ROOT / "tests/performance/cross-solver-yagi.nec",
    ]
    cases = []
    with tempfile.TemporaryDirectory(prefix="nec2c-cross-solver-") as temporary:
        base = Path(temporary)
        for index, deck_argument in enumerate(decks):
            deck = deck_argument.resolve()
            directory = base / str(index)
            directory.mkdir()
            values = {
                "nec2c": impedance_from_nec2c(arguments.nec2c.resolve(), deck, directory,
                                                arguments.nec2c_backend),
                "aegnec2": impedance_from_aeg(arguments.aegnec2.resolve(), deck, directory),
                "momwire-sinusoidal": impedance_from_momwire(
                    arguments.momwire_root.resolve(), deck),
            }
            reference = values["nec2c"]
            comparisons = {
                name: delta(reference, value) for name, value in values.items()
                if name != "nec2c"
            }
            case_passed = (
                comparisons["aegnec2"]["relative"] <= arguments.aeg_relative_tolerance and
                comparisons["momwire-sinusoidal"]["relative"] <= arguments.momwire_relative_tolerance)
            cases.append({
                "deck": str(deck),
                "values": {name: {"real": value.real, "imaginary": value.imag}
                           for name, value in values.items()},
                "comparisons_to_nec2c": comparisons,
                "passed": case_passed,
            })
    passed = all(case["passed"] for case in cases)
    report = {
        "schema": "nec2c-cross-solver-comparison",
        "schema_version": 1,
        "observable": "complex feed impedance, ohms",
        "nec2c_backend": arguments.nec2c_backend,
        "cases": cases,
        "tolerances": {
            "aegnec2_relative": arguments.aeg_relative_tolerance,
            "momwire_sinusoidal_relative": arguments.momwire_relative_tolerance,
        },
        "passed": passed,
        "notes": [
            "AEG is executed as a GPL external oracle; no AEG source or object code is copied.",
            "AEG's default result file prints impedance to six significant digits, which limits comparison precision.",
            "Momwire sinusoidal uses the NEC2 three-term basis but has an independent implementation.",
        ],
    }
    arguments.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                encoding="utf-8")
    for case in cases:
        print(Path(case["deck"]).name)
        for name, value in case["values"].items():
            print(f"  {name}: {value['real']:.12g} {value['imaginary']:+.12g}j ohm")
    print(f"comparison: {'PASS' if passed else 'FAIL'}")
    if not passed or any(not math.isfinite(part)
                         for case in cases for value in case["values"].values()
                         for part in (value["real"], value["imaginary"])):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
