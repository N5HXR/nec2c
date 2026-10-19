#!/usr/bin/env python3
"""Numerical regression runner for legacy (non-dielectric) NEC2C behavior."""

import argparse
import copy
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


FLOAT = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")

CASES = {
    "ex1-wire-near-field": {
        "deck": "Input/EX1.nec",
        "features": ["wire", "voltage-source", "near-electric-field"],
        "metrics": [
            ["input", "ANTENNA INPUT PARAMETERS", 0, 0, [4, 5, 6, 7, 10]],
            ["feed-current", "CURRENTS AND LOCATION", 0, 3, [6, 7]],
            ["near-e", "NEAR ELECTRIC FIELDS", 0, 1, [3, 4, 7, 8]],
        ],
    },
    "ex4-mixed-ground": {
        "deck": "Input/EX4.nec",
        "features": ["wire", "surface-patch", "mixed-geometry", "perfect-ground", "radiation-pattern"],
        "metrics": [
            ["input", "ANTENNA INPUT PARAMETERS", 0, 0, [4, 5, 6, 7, 10]],
            ["pattern", "RADIATION PATTERNS", 0, 5, [0, 1, 4, 7, 8, 9, 10]],
        ],
    },
    "ex5-transmission-line": {
        "deck": "Input/EX5.nec",
        "features": ["wire", "TL", "radiation-pattern"],
        "metrics": [
            ["network", "NETWORK DATA", 0, 0, [4, 5, 6, 7, 8, 9]],
            ["input", "ANTENNA INPUT PARAMETERS", 0, 0, [4, 5, 6, 7, 10]],
            ["pattern", "RADIATION PATTERNS", 0, 0, [0, 1, 4, 7, 8, 9, 10]],
        ],
    },
    "ex7-plane-wave": {
        "deck": "Input/EX7.nec",
        "features": ["wire", "plane-wave", "receiving-pattern"],
        "metrics": [
            ["pattern-a", "RADIATION PATTERNS", 0, 0, [0, 1, 4, 7, 8, 9, 10]],
            ["pattern-b", "RADIATION PATTERNS", 1, 0, [0, 1, 4, 7, 8, 9, 10]],
        ],
    },
    "patch-only": {
        "deck": "tests/decks/patch-only.nec",
        "features": ["surface-patch", "plane-wave", "radiation-pattern"],
        "metrics": [
            ["patch-current", "SURFACE PATCH CURRENTS", 0, 0, [4, 5, 6, 7]],
            ["pattern", "RADIATION PATTERNS", 0, 0, [0, 1, 4, 7, 8, 9, 10]],
        ],
    },
    "nt-network": {
        "deck": "tests/decks/nt-network.nec",
        "features": ["wire", "NT", "voltage-source", "radiation-pattern"],
        "metrics": [
            ["network", "NETWORK DATA", 0, 0, [4, 5, 6, 7, 8, 9]],
            ["input", "ANTENNA INPUT PARAMETERS", 0, 0, [4, 5, 6, 7, 10]],
            ["pattern", "RADIATION PATTERNS", 0, 0, [0, 1, 4, 7, 8, 9, 10]],
        ],
    },
}


def numeric_rows(text, heading, occurrence):
    matches = list(re.finditer(re.escape(heading), text))
    if occurrence >= len(matches):
        raise ValueError(f"missing occurrence {occurrence + 1} of section {heading!r}")
    rows = []
    for line in text[matches[occurrence].end():].splitlines()[1:120]:
        if re.match(r"^\s*[-+]?\d", line):
            rows.append([float(value) for value in FLOAT.findall(line)])
        if rows and "DATA CARD No:" in line:
            break
    return rows


def extract_metrics(text, case):
    result = {}
    for name, heading, occurrence, row_index, columns in case["metrics"]:
        rows = numeric_rows(text, heading, occurrence)
        if row_index >= len(rows):
            raise ValueError(f"{heading!r} has {len(rows)} numerical rows, need row {row_index}")
        row = rows[row_index]
        try:
            result[name] = [row[column] for column in columns]
        except IndexError as error:
            raise ValueError(f"short numerical row in {heading!r}: {row}") from error
    return result


def run_cases(binary, source_root):
    actual = {}
    with tempfile.TemporaryDirectory(prefix="nec2c-regression-") as temp_dir:
        for name, case in CASES.items():
            output = Path(temp_dir) / f"{name}.out"
            command = [str(binary), "-i", str(source_root / case["deck"]), "-o", str(output)]
            process = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if process.returncode != 0:
                raise RuntimeError(
                    f"{name}: NEC2C exited {process.returncode}\n{process.stdout}{process.stderr}"
                )
            actual[name] = extract_metrics(output.read_text(encoding="utf-8"), case)
    return actual


def compare(actual, expected, rel_tol, abs_tol):
    failures = []
    for case_name, metrics in expected.items():
        if case_name not in actual:
            failures.append(f"{case_name}: case was not run")
            continue
        for metric_name, reference in metrics.items():
            observed = actual[case_name].get(metric_name)
            if observed is None:
                failures.append(f"{case_name}/{metric_name}: metric missing")
                continue
            if len(observed) != len(reference):
                failures.append(f"{case_name}/{metric_name}: length {len(observed)} != {len(reference)}")
                continue
            for index, (got, want) in enumerate(zip(observed, reference)):
                if not math.isclose(got, want, rel_tol=rel_tol, abs_tol=abs_tol):
                    failures.append(
                        f"{case_name}/{metric_name}[{index}]: got {got:.9g}, "
                        f"expected {want:.9g} (rel={rel_tol:g}, abs={abs_tol:g})"
                    )
    return failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate", action="store_true", help="print a regenerated baseline as JSON")
    parser.add_argument("--binary", type=Path, default=Path(os.environ.get("NEC2C", "./nec2c")))
    parser.add_argument("--source-root", type=Path, default=Path(os.environ.get("TEST_SRCDIR", ".")))
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    binary = args.binary.resolve()
    actual = run_cases(binary, source_root)
    if args.generate:
        print(json.dumps(actual, indent=2, sort_keys=True))
        return 0

    baseline_path = source_root / "tests" / "baseline.json"
    document = json.loads(baseline_path.read_text(encoding="utf-8"))
    failures = compare(actual, document["cases"], document["relative_tolerance"], document["absolute_tolerance"])

    # Negative control: a one-percent perturbation must be diagnosed.
    perturbed = copy.deepcopy(document["cases"])
    perturbed["ex1-wire-near-field"]["input"][0] *= 1.01
    if not compare(actual, perturbed, document["relative_tolerance"], document["absolute_tolerance"]):
        failures.append("runner negative control did not detect a one-percent numerical change")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    feature_count = len({feature for case in CASES.values() for feature in case["features"]})
    metric_count = sum(len(metrics) for metrics in actual.values())
    print(f"PASS: {len(CASES)} legacy decks, {feature_count} feature classes, {metric_count} numerical samples")
    print("PASS: negative control rejected a one-percent reference perturbation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
