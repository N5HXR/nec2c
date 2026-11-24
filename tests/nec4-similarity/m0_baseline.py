#!/usr/bin/env python3
"""Verify the immutable pre-NEC4 numerical baseline."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
DEFAULT_BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
HERE = ROOT / "tests" / "nec4-similarity"
GOLDEN = HERE / "m0-baseline.json"
CASES = (
    ("free-space", ROOT / "tests" / "coated-wire" / "m0-bare-dipole.nec"),
    ("real-ground", ROOT / "tests" / "performance" /
     "m3-real-ground-sweep.nec"),
)


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def digest_json(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         allow_nan=False).encode("utf-8")
    return digest_bytes(encoded)


def physics_projection(document):
    """Select stable physical results, excluding timing/backend metadata."""
    projected_jobs = []
    for job in document["jobs"]:
        ground = {entry["ground_state_id"]: entry["ground"]
                  for entry in job["ground_states"]}
        executions = []
        for execution in job["executions"]:
            executions.append({
                "status": execution["status"],
                "frequency_mhz": execution["frequency_mhz"],
                "ground": ground[execution["ground_state_id"]],
                "excitation": execution["excitation"],
                "trigger": execution["trigger"],
                "material_evaluations": execution["material_evaluations"],
                "feeds": execution["feeds"],
                "wire_currents": execution["wire_currents"],
                "patch_currents": execution["patch_currents"],
                "volume_currents": execution["volume_currents"],
                "near_fields": execution["near_fields"],
                "far_fields": execution["far_fields"],
                "power": execution["power"],
            })
        projected_jobs.append({"dielectric": job["dielectric"],
                               "executions": executions})
    return projected_jobs


def normalized_report(text):
    text = re.sub(r"TOTAL RUN TIME: \d+ msec",
                  "TOTAL RUN TIME: X msec", text)
    text = re.sub(r"FILL: \d+ msec  FACTOR: \d+ msec",
                  "FILL: X msec  FACTOR: X msec", text)
    return text


def supports_engine_flags(binary):
    completed = subprocess.run([str(binary), "-h"], cwd=ROOT, text=True,
                               capture_output=True)
    return "--wire-engine" in (completed.stdout + completed.stderr)


def run_case(binary, name, deck, temporary, engine_flags):
    result = temporary / f"{name}.json"
    report = temporary / f"{name}.out"
    command = [str(binary), "--backend", "builtin", "--fft-backend",
               "builtin", "--jobs", "1", "--threads", "1"]
    command.extend(engine_flags)
    command.extend(["--result-json", str(result), "-i", str(deck),
                    "-o", str(report)])
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((completed.returncode, completed.stdout,
                              completed.stderr))

    document = json.loads(result.read_text(encoding="utf-8"))
    executions = [execution for job in document["jobs"]
                  for execution in job["executions"]]
    return {
        "deck": str(deck.relative_to(ROOT)),
        "deck_sha256": digest_bytes(deck.read_bytes()),
        "physics_sha256": digest_json(physics_projection(document)),
        "normalized_report_sha256": digest_bytes(
            normalized_report(report.read_text(encoding="utf-8")).encode(
                "utf-8")),
        "frequencies_mhz": [item["frequency_mhz"] for item in executions],
        "feed_impedance_ohm": [item["feeds"][0]["impedance_ohm"]
                                for item in executions],
        "wire_current_count": sum(len(item["wire_currents"])
                                  for item in executions),
        "far_field_sample_count": sum(
            len(request["samples"]) for item in executions
            for request in item["far_fields"]),
    }


def observe(binary):
    engine_flags = (["--wire-engine", "nec2", "--ground-engine", "nec2"]
                    if supports_engine_flags(binary) else [])
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m0-") as directory:
        temporary = Path(directory)
        return {
            "classic_selection": engine_flags or "pre-selector baseline",
            "cases": [run_case(binary, name, deck, temporary, engine_flags)
                      for name, deck in CASES],
        }


parser = argparse.ArgumentParser()
parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
parser.add_argument("--record", action="store_true",
                    help="print observed baseline instead of checking it")
arguments = parser.parse_args()
observed = observe(arguments.binary.resolve())

if arguments.record:
    print(json.dumps(observed, indent=2, sort_keys=True, allow_nan=False))
else:
    expected = json.loads(GOLDEN.read_text(encoding="utf-8"))
    if len(observed["cases"]) != len(expected["cases"]):
        raise AssertionError(
            f"case count: expected {len(expected['cases'])}, "
            f"got {len(observed['cases'])}")
    for observed_case, expected_case in zip(observed["cases"],
                                            expected["cases"]):
        for key in ("deck", "deck_sha256", "physics_sha256",
                    "normalized_report_sha256", "frequencies_mhz",
                    "feed_impedance_ohm", "wire_current_count",
                    "far_field_sample_count"):
            if observed_case[key] != expected_case[key]:
                raise AssertionError(
                    f"{observed_case['deck']} {key}: expected "
                    f"{expected_case[key]!r}, got {observed_case[key]!r}")
    print("NEC4-similarity M0 classic free-space and real-ground baselines "
          "are unchanged")
