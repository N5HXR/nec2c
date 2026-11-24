#!/usr/bin/env python3
"""M6 above-ground, combined-feature, field, power and concurrency gates."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def scrub(value):
    if isinstance(value, dict):
        return {key: scrub(item) for key, item in value.items()
                if key != "elapsed_seconds"}
    if isinstance(value, list):
        return [scrub(item) for item in value]
    return value


def numeric_delta(left, right):
    values = []
    if isinstance(left, dict) and isinstance(right, dict):
        for key in left.keys() & right.keys():
            if key not in {"elapsed_seconds", "threads", "requested", "used"}:
                values.extend(numeric_delta(left[key], right[key]))
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise AssertionError("result length mismatch")
        for a, b in zip(left, right):
            values.extend(numeric_delta(a, b))
    elif (isinstance(left, (int, float)) and not isinstance(left, bool) and
          isinstance(right, (int, float)) and not isinstance(right, bool)):
        values.append(abs(left-right)/max(1., abs(left), abs(right)))
    return values


def run(deck, directory, stem, engine="nec4", jobs=1, threads=1,
        expect_success=True):
    model = directory / f"{stem}.nec"
    report = directory / f"{stem}.out"
    result = directory / f"{stem}.json"
    model.write_text(deck, encoding="ascii")
    command = [str(BINARY), "--wire-engine", engine, "--ground-engine", "nec2",
               "--backend", "builtin", "--jobs", str(jobs), "--threads",
               str(threads), "--result-json", str(result), "-i", str(model),
               "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                               timeout=180)
    if expect_success and completed.returncode:
        raise AssertionError((command, completed.stdout, completed.stderr,
                              report.read_text(errors="replace")[-3000:]
                              if report.exists() else ""))
    if not expect_success:
        if completed.returncode == 0:
            raise AssertionError(f"unsupported combination succeeded: {stem}")
        return None, report.read_text(errors="replace"), completed.stderr
    return json.loads(result.read_text()), report.read_text(errors="replace"), \
        completed.stderr


def horizontal_dipole(ground, frequencies=3):
    return ("CM M6 HORIZONTAL DIPOLE ABOVE GROUND\nCE\n"
            "GW 1 41 -.24 0 .12 .24 0 .12 .001\nGE 0\n"
            f"FR 0 {frequencies} 0 0 180 10\n{ground}\n"
            "EX 0 1 21 0 1 0\n"
            "RP 0 10 8 1000 0 0 10 45 1000 0\nEN\n")


def pec_boundary():
    return ("CM M6 PEC CAP IMAGE BOUNDARY\nCE\n"
            "GW 1 41 0 0 .001 0 0 .251 .001\nGE 1\n"
            "FR 0 1 0 0 299.792458 0\nGN 1\nEX 0 1 1 0 1 0\n"
            "NE 0 5 1 1 .05 0 0 .05 0 0\nEN\n")


def reciprocal_pair(ground, source):
    return ("CM M6 GROUND RECIPROCITY\nCE\n"
            "GW 1 21 -.24 -.08 .15 .24 -.08 .15 .001\n"
            "GW 2 21 -.24 .08 .15 .24 .08 .15 .001\nGE 0\n"
            f"FR 0 1 0 0 180 0\n{ground}\nEX 0 {source} 11 0 1 0\nXQ\nEN\n")


def impedance(document):
    item = document["jobs"][0]["executions"][0]["feeds"][0]["impedance_ohm"]
    return complex(item["re"], item["im"])


def transfer(document, target):
    rows = document["jobs"][0]["executions"][0]["wire_currents"]
    item = rows[(target-1)*21+10]["current_a"]
    return complex(item["re"], item["im"])


def main():
    evidence = {"schema": "nec2c-nec4-m6-acceptance", "schema_version": 1}
    grounds = {"pec": "GN 1", "rca": "GN 0 0 0 0 13 .02",
               "sommerfeld": "GN 2 0 0 0 13 .02"}
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m6-") as name:
        directory = Path(name)

        boundary, _, _ = run(pec_boundary(), directory, "pec-boundary")
        samples = boundary["jobs"][0]["executions"][0]["near_fields"][0]["samples"]
        tangential = max(math.hypot(row["value"][axis]["re"],
                                    row["value"][axis]["im"])
                         for row in samples for axis in ("x", "y"))
        normal = max(math.hypot(row["value"]["z"]["re"],
                                row["value"]["z"]["im"]) for row in samples)
        boundary_defect = tangential/max(normal, 1e-300)
        if boundary_defect > 1e-12:
            raise AssertionError(("PEC tangential boundary", boundary_defect))
        evidence["pec_cap_image"] = {"tangential_relative_defect": boundary_defect}

        ground_results = {}
        for label, card in grounds.items():
            serial, report, _ = run(horizontal_dipole(card), directory,
                                    f"{label}-sweep", jobs=1)
            parallel, _, _ = run(horizontal_dipole(card), directory,
                                  f"{label}-sweep", jobs=2)
            if scrub(serial) != scrub(parallel):
                raise AssertionError(f"{label} serial/parallel mismatch")
            executions = serial["jobs"][0]["executions"]
            if len(executions) != 3 or any(not e["far_fields"] for e in executions):
                raise AssertionError(f"{label} field sweep incomplete")
            if any(e["power"]["accepted_w"] <= 0 for e in executions):
                raise AssertionError(f"{label} non-passive accepted power")
            old, _, _ = run(horizontal_dipole(card, 1), directory,
                            f"{label}-nec2", engine="nec2")
            new_z = impedance(serial)
            old_z = impedance(old)
            ground_results[label] = {
                "impedance_nec4_ohm": {"re": new_z.real, "im": new_z.imag},
                "impedance_nec2_ohm": {"re": old_z.real, "im": old_z.imag},
                "nec2_nec4_relative_delta": abs(new_z-old_z)/max(abs(old_z), 1.),
                "serial_parallel_exact": True,
                "sommerfeld_notice": "cap-specific Sommerfeld correction remains provisional"
                    in report,
            }

            transfers = []
            for source, target in ((1, 2), (2, 1)):
                result, _, _ = run(reciprocal_pair(card, source), directory,
                                   f"{label}-reciprocity-{source}")
                transfers.append(transfer(result, target))
            defect = abs(transfers[0]-transfers[1])/max(1., *map(abs, transfers))
            if defect > 2e-10:
                raise AssertionError((label, "reciprocity", defect, transfers))
            ground_results[label]["reciprocity_scaled_defect"] = defect
        if ground_results["pec"]["sommerfeld_notice"] or \
           ground_results["rca"]["sommerfeld_notice"] or \
           not ground_results["sommerfeld"]["sommerfeld_notice"]:
            raise AssertionError("Sommerfeld qualification notice is not scoped to GN 2")
        evidence["ground_modes"] = ground_results

        patch = (ROOT / "tests/performance/m6-connected-patch-657.nec").read_text()
        patch_result, patch_report, _ = run(patch, directory, "patch")
        execution = patch_result["jobs"][0]["executions"][0]
        if not execution["patch_currents"] or "surface patches enabled" not in patch_report:
            raise AssertionError("wire/patch combination was not exercised")
        evidence["wire_patch"] = {
            "wire_unknowns": len(execution["wire_currents"]),
            "patch_current_rows": len(execution["patch_currents"]),
            "impedance_ohm": execution["feeds"][0]["impedance_ohm"],
        }

        network = (ROOT / "tests/performance/m6-loaded-network-402.nec").read_text()
        serial, _, _ = run(network, directory, "network-serial", threads=1)
        parallel, _, _ = run(network, directory, "network-parallel", threads=8)
        delta = max(numeric_delta(serial, parallel), default=0.)
        if delta > 1e-12:
            raise AssertionError(("network serial/parallel", delta))
        evidence["loaded_network"] = {"serial_parallel_scaled_delta": delta}

        coated = (ROOT / "tests/coated-wire/m3-coated-dipole.nec").read_text()
        _, coated_report, _ = run(coated, directory, "unsupported-coated",
                                  expect_success=False)
        if "LD 7 coated wire has not been qualified" not in coated_report:
            raise AssertionError("coated-wire refusal missing")
        dielectric = (ROOT / "tests/dielectric/m5-cube.nec").read_text()
        _, dielectric_report, dielectric_error = run(
            dielectric, directory, "unsupported-dielectric", expect_success=False)
        if ("dielectric VIE has not been qualified" not in dielectric_report or
                "not yet supported" not in dielectric_error):
            raise AssertionError("dielectric refusal missing")
        evidence["combination_gates"] = {
            "coated_wire": "rejected explicitly",
            "dielectric_vie": "rejected explicitly",
        }

    output = ROOT / "tests/nec4-similarity/m6-acceptance.json"
    output.write_text(json.dumps(evidence, indent=2)+"\n", encoding="utf-8")
    print(f"M6 PEC cap-image tangential defect {boundary_defect:.3e}")
    for label, row in ground_results.items():
        print(f"M6 {label}: reciprocity {row['reciprocity_scaled_defect']:.3e}, "
              f"NEC2/NEC4 impedance delta {row['nec2_nec4_relative_delta']:.3%}, "
              "serial/parallel exact")
    print(f"M6 patch rows {evidence['wire_patch']['patch_current_rows']}; "
          f"loaded/network parallel delta {delta:.3e}; unsupported combinations gated")


if __name__ == "__main__":
    main()
