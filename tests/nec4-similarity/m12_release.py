#!/usr/bin/env python3
"""M12 automatic-selection and explicit-classic release gates."""

import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def projection(document):
    jobs = []
    for job in document["jobs"]:
        executions = []
        for execution in job["executions"]:
            executions.append({key: execution[key] for key in (
                "status", "frequency_mhz", "ground_state_id", "excitation",
                "trigger", "material_evaluations", "feeds", "wire_currents",
                "patch_currents", "volume_currents", "power", "far_fields",
                "near_fields")})
        jobs.append({
            "dielectric": job["dielectric"],
            "coated_wires": job.get("coated_wires"),
            "ground_states": job["ground_states"],
            "executions": executions,
        })
    return jobs


def run(name, deck, directory, options):
    result = directory / f"{name}.json"
    report = directory / f"{name}.out"
    command = [str(BINARY), "--backend", "builtin", "--jobs", "1",
               "--threads", "1", "--result-json", str(result), *options,
               "-i", str(deck), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((command, completed.returncode,
                              completed.stdout, completed.stderr))
    document = json.loads(result.read_text(encoding="utf-8"))
    selection = [item for item in document["diagnostics"]
                 if item["code"] == "ENGINE_AUTO_SELECTION"]
    automatic = (("--wire-engine" in options and
                  options[options.index("--wire-engine") + 1] == "auto") or
                 ("--ground-engine" in options and
                  options[options.index("--ground-engine") + 1] == "auto"))
    if len(selection) != (1 if automatic else 0):
        raise AssertionError((name, "selection diagnostic", selection))
    engines = document["jobs"][0]["engines"]
    context = selection[0]["context"] if automatic else {
        "wire": engines["wire"], "ground": engines["ground"],
        "reason": engines["selection_reason"]}
    if engines["wire"] != context["wire"] or \
       engines["ground"] != context["ground"] or \
       engines["selection_reason"] != context["reason"]:
        raise AssertionError((name, "job engine metadata", engines,
                              context))
    return document, context, report.read_text(encoding="utf-8")


def pair(name, deck, directory, expected, explicit):
    automatic, selection, report = run(
        name + "-auto", deck, directory,
        ["--wire-engine", "auto", "--ground-engine", "auto"])
    oracle, _, _ = run(name + "-explicit", deck, directory, explicit)
    resolved = (selection["wire"]["resolved"],
                selection["ground"]["resolved"])
    if resolved != expected:
        raise AssertionError((name, resolved, expected, selection))
    if "ENGINE AUTO SELECTION:" not in report:
        raise AssertionError((name, "missing human selection report"))
    if projection(automatic) != projection(oracle):
        raise AssertionError((name, "automatic/explicit numerical mismatch"))
    return selection["reason"]


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4ish-m12-") as name:
        directory = Path(name)
        reasons = {}
        reasons["free_space"] = pair(
            "free", ROOT / "tests/coated-wire/m0-bare-dipole.nec",
            directory, ("nec4ish", "nec2"),
            ["--wire-engine", "nec4ish", "--ground-engine", "nec2"])
        reasons["real_ground"] = pair(
            "ground",
            ROOT / "tests/nec4-similarity/m10a-identical-buried.nec",
            directory, ("nec4ish", "nec4ish"),
            ["--wire-engine", "nec4ish", "--ground-engine", "nec4ish",
             "--nec4ish-ground-model", "grid"])
        reasons["patch"] = pair(
            "patch", ROOT / "tests/decks/patch-only.nec", directory,
            ("nec4ish", "nec2"),
            ["--wire-engine", "nec4ish", "--ground-engine", "nec2"])
        reasons["coated"] = pair(
            "coated", ROOT / "tests/coated-wire/m3-coated-dipole.nec",
            directory, ("nec2", "nec2"),
            ["--wire-engine", "nec2", "--ground-engine", "nec2"])

        dielectric, selection, report = run(
            "dielectric-auto", ROOT / "tests/dielectric/valid.nec", directory,
            ["--wire-engine", "auto", "--ground-engine", "auto"])
        if (selection["wire"]["resolved"],
                selection["ground"]["resolved"]) != ("nec2", "nec2"):
            raise AssertionError(("dielectric fallback", selection))
        if "ENGINE AUTO SELECTION:" not in report:
            raise AssertionError("missing dielectric fallback report")
        reasons["dielectric"] = selection["reason"]

        # Explicit classic selection must not add a human-facing line; this is
        # necessary for byte-stable normalized M0 reports.
        default, selection, report = run(
            "classic-default",
            ROOT / "tests/coated-wire/m0-bare-dipole.nec", directory, [])
        if "ENGINE AUTO SELECTION:" in report:
            raise AssertionError("explicit NEC2 report changed")
        if selection["wire"]["resolved"] != "nec2":
            raise AssertionError(("default classic", selection))
        explicit, _, _ = run(
            "classic-explicit",
            ROOT / "tests/coated-wire/m0-bare-dipole.nec", directory,
            ["--wire-engine", "nec2", "--ground-engine", "nec2"])
        if projection(default) != projection(explicit):
            raise AssertionError("default and explicit NEC2 differ")

    print("M12 classic default, opt-in automatic selection, and explicit NEC2 identity passed")
    print("  " + ", ".join(f"{key}={value}" for key, value in reasons.items()))


if __name__ == "__main__":
    main()
