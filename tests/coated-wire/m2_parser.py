#!/usr/bin/env python3
"""M2 acceptance tests for LD 7 parsing, topology, and result metadata."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", "./nec2c")).resolve()


def geometry(lines=None, ge="GE"):
    if lines is None:
        lines = ["GW 1 11 0 0 -0.25 0 0 0.25 0.001"]
    return "CM M2 LD7 TEST\nCE\n" + "\n".join(lines) + f"\n{ge}\n"


def deck(ld_cards, lines=None, ge="GE", trailer=None):
    if trailer is None:
        trailer = "FR 0 1 0 0 300\nEX 0 1 6 0 1 0\nEN\n"
    return geometry(lines, ge) + "\n".join(ld_cards) + "\n" + trailer


def invoke(directory, text, validate=True):
    source = directory / "case.nec"
    report = directory / "case.out"
    result = directory / "case.json"
    source.write_text(text, encoding="ascii")
    command = [str(BINARY)]
    if validate:
        command.append("--validate-only")
    command += ["--result-json", str(result), "-i", str(source), "-o", str(report)]
    completed = subprocess.run(command, cwd=directory, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    parsed = json.loads(result.read_text(encoding="utf-8")) if result.exists() else None
    human = report.read_text(encoding="utf-8") if report.exists() else ""
    return completed.returncode, parsed, human + completed.stdout + completed.stderr


def expect_error(root, name, text, code):
    directory = root / name
    directory.mkdir()
    status, result, output = invoke(directory, text)
    observed = result["diagnostics"][0]["code"] if result and result["diagnostics"] else None
    if status == 0 or observed != code or result["completion"]["status"] != "validation-error":
        raise AssertionError(f"{name}: expected {code}, got status={status}, code={observed}\n{output}")


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m2-") as temporary:
        root = Path(temporary)

        # Whole-tag shorthand, explicit whole range, omitted optional loss fields,
        # and disconnected domains with different materials are all legal.
        valid = [
            ("whole", deck(["LD 7 1 0 0 2.8 0.003 0.01 0.0002"]), 1),
            ("explicit", deck(["LD 7 1 1 11 2.8 0.003"]), 1),
            ("two-disconnected", deck([
                "LD 7 1 0 0 2.1 0.002",
                "LD 7 2 0 0 4.0 0.004 0.02 0.001",
            ], lines=[
                "GW 1 5 0 0 -0.25 0 0 -0.05 0.001",
                "GW 2 5 0.1 0 0.05 0.1 0 0.25 0.001",
            ], trailer="EN\n"), 2),
            ("equal-junction", deck([
                "LD 7 1 0 0 2.8 0.003 0.01 0.0002",
                "LD 7 2 0 0 2.8 0.003 0.01 0.0002",
            ], lines=[
                "GW 1 5 0 0 -0.25 0 0 0 0.001",
                "GW 2 5 0 0 0 0 0 0.25 0.001",
            ], trailer="EN\n"), 2),
        ]
        for name, text, count in valid:
            directory = root / name
            directory.mkdir()
            status, result, output = invoke(directory, text)
            if status != 0:
                raise AssertionError(f"{name}: valid deck failed\n{output}")
            job = result["jobs"][0]
            assert job["validation"] == "coated-wire-preflight"
            assert job["coated_wires"]["operator_status"] == "validation-only"
            assert len(job["coated_wires"]["domains"]) == count
            resolved = [segment for item in job["coated_wires"]["domains"]
                        for segment in item["segments"]]
            assert len(resolved) == sum(item["segment_count"]
                                        for item in job["coated_wires"]["domains"])
            assert len(resolved) == len(set(resolved))
            assert not result["diagnostics"]

        whole = json.loads((root / "whole" / "case.json").read_text())["jobs"][0]
        domain = whole["coated_wires"]["domains"][0]
        assert domain == {
            "tag": 1, "tag_relative_range": [1, 11], "segment_count": 11,
            "segments": list(range(1, 12)),
            "conductor_radius_m": 0.001, "outer_radius_m": 0.003,
            "er": 2.8, "tan_delta": 0.01, "sigma_dc_s_per_m": 0.0002,
        }

        invalid = [
            ("zero-tag", deck(["LD 7 0 0 0 2 0.003"]), "COATED_WIRE_INVALID_TAG"),
            ("missing-tag", deck(["LD 7 9 0 0 2 0.003"]), "COATED_WIRE_TAG_NOT_FOUND"),
            ("partial", deck(["LD 7 1 2 10 2 0.003"]), "COATED_WIRE_PARTIAL_COVERAGE_UNSUPPORTED"),
            ("half-zero-range", deck(["LD 7 1 0 11 2 0.003"]), "COATED_WIRE_INVALID_RANGE"),
            ("reversed-range", deck(["LD 7 1 8 2 2 0.003"]), "COATED_WIRE_INVALID_RANGE"),
            ("er-low", deck(["LD 7 1 0 0 0.9 0.003"]), "COATED_WIRE_INVALID_ER"),
            ("er-nan", deck(["LD 7 1 0 0 nan 0.003"]), "COATED_WIRE_INVALID_ER"),
            ("radius-zero", deck(["LD 7 1 0 0 2 0"]), "COATED_WIRE_INVALID_OUTER_RADIUS"),
            ("radius-small", deck(["LD 7 1 0 0 2 0.0009"]), "COATED_WIRE_OUTER_RADIUS_TOO_SMALL"),
            ("loss-negative", deck(["LD 7 1 0 0 2 0.003 -0.1"]), "COATED_WIRE_INVALID_LOSS_TANGENT"),
            ("sigma-negative", deck(["LD 7 1 0 0 2 0.003 0 -1"]), "COATED_WIRE_INVALID_CONDUCTIVITY"),
            ("reserved-field", deck(["LD 7 1 0 0 2 0.003 0 0 1"]), "COATED_WIRE_RESERVED_FIELD_NONZERO"),
            ("duplicate", deck(["LD 7 1 0 0 2 0.003", "LD 7 1 0 0 2 0.003"]), "COATED_WIRE_DUPLICATE_TAG"),
            ("mixed-radius", deck(["LD 7 1 0 0 2 0.003"], lines=[
                "GW 1 5 0 0 -0.25 0 0 0 0.001",
                "GW 1 5 0 0 0 0 0 0.25 0.0015",
            ], trailer="EN\n"), "COATED_WIRE_MIXED_CONDUCTOR_RADIUS"),
            ("bare-junction", deck(["LD 7 1 0 0 2 0.003"], lines=[
                "GW 1 5 0 0 -0.25 0 0 0 0.001",
                "GW 2 5 0 0 0 0 0 0.25 0.001",
            ], trailer="EN\n"), "COATED_WIRE_BARE_JUNCTION_UNSUPPORTED"),
            ("unequal-junction", deck([
                "LD 7 1 0 0 2 0.003", "LD 7 2 0 0 3 0.003",
            ], lines=[
                "GW 1 5 0 0 -0.25 0 0 0 0.001",
                "GW 2 5 0 0 0 0 0 0.25 0.001",
            ], trailer="EN\n"), "COATED_WIRE_JUNCTION_MISMATCH"),
            ("ground-junction", deck(["LD 7 1 0 0 2 0.003"], lines=[
                "GW 1 5 0 0 0 0 0 0.25 0.001",
            ], ge="GE 1", trailer="EN\n"), "COATED_WIRE_GROUND_JUNCTION_UNSUPPORTED"),
        ]
        for name, text, code in invalid:
            expect_error(root, name, text, code)

        # M4 promotes M2's lossy parser witness through the same isolated-wire
        # topology gate.  EN alone performs no solve but must no longer reject
        # the loss fields.
        execution_dir = root / "execution-gate"
        execution_dir.mkdir()
        status, result, output = invoke(execution_dir, valid[0][1], validate=False)
        assert status == 0
        assert result["jobs"][0]["coated_wires"]["operator_status"] == "lossy-isolated"

        # Every legacy LD type still reaches the classic load implementation.
        legacy_cards = {
            0: "LD 0 1 1 1 10 1e-8 1e-12",
            1: "LD 1 1 1 1 1000 1e-8 1e-12",
            2: "LD 2 1 1 1 5.8e7",
            3: "LD 3 1 1 1 10 2",
            4: "LD 4 1 1 1 10 2",
            5: "LD 5 1 1 1 10 2",
        }
        for kind, card in legacy_cards.items():
            directory = root / f"legacy-{kind}"
            directory.mkdir()
            legacy = deck([card], trailer="FR 0 1 0 0 300\nEX 0 1 6 0 1 0\nXQ\nEN\n")
            status, result, output = invoke(directory, legacy, validate=False)
            if status != 0 or "STRUCTURE IMPEDANCE LOADING" not in output:
                raise AssertionError(f"legacy LD {kind} regressed\n{output}")
            assert "coated_wires" not in result["jobs"][0]

    print("M2 LD 7 parser: 4 valid topology/metadata cases passed")
    print("M2 LD 7 parser: 17 invalid card/topology diagnostics passed")
    print("M2 LD 7 loss promotion and LD 0-5 regression cases passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
