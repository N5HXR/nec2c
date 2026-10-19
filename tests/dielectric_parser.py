#!/usr/bin/env python3
"""Acceptance tests for the M2 DF card and dielectric resource parser."""

import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


BASIC = """NEC2C-DIELECTRIC 1
MATERIAL 1 ER 3.0
REGION 1 MATERIAL 1
BOX 1 0 0 0 1 1 1
MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2.0
END
"""


def invoke(binary, deck, cwd):
    output = cwd / ("nec2c-" + deck.stem + ".out")
    process = subprocess.run(
        [str(binary), "-i", str(deck), "-o", str(output)],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    text = output.read_text(encoding="utf-8") if output.exists() else ""
    return process.returncode, text + process.stdout + process.stderr


def write_case(root, resource_text=BASIC, df_cards=None, trailer="EN\n"):
    case = root / "case"
    resources = case / "resources"
    resources.mkdir(parents=True)
    resource = resources / "model.diel"
    resource.write_text(resource_text, encoding="utf-8")
    if df_cards is None:
        df_cards = ['DF "resources/model.diel"']
    deck = case / "model.nec"
    deck.write_text(
        "CE\nGW 1 1 100 0 0 100 0 1 0.001\nGE\nFR 0 1 0 0 1.0\n"
        + "\n".join(df_cards)
        + "\n"
        + trailer,
        encoding="utf-8",
    )
    return deck, resource


def expect_failure(binary, temp_root, name, resource_text, expected, df_cards=None):
    root = temp_root / name
    deck, _ = write_case(root, resource_text, df_cards)
    status, output = invoke(binary, deck, temp_root)
    if status == 0 or expected not in output:
        raise AssertionError(
            f"{name}: expected failure containing {expected!r}; status={status}\n{output}"
        )


def normalized_summary(output):
    match = re.search(
        r" {23}-------- DIELECTRIC RESOURCE --------.*?"
        r"STATUS: PARSED; READY FOR DIELECTRIC SOLVER",
        output,
        re.DOTALL,
    )
    if match is None:
        raise AssertionError("canonical dielectric summary is missing")
    return re.sub(
        r"(?m)^  RESOLVED RESOURCE: .+$",
        "  RESOLVED RESOURCE: <RESOURCE>",
        match.group(0),
    )


def main():
    binary = Path(os.environ.get("NEC2C", "./nec2c")).resolve()
    source_root = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
    fixture = source_root / "tests" / "dielectric" / "valid.nec"
    expected = (source_root / "tests" / "dielectric" / "expected-summary.txt").read_text(
        encoding="utf-8"
    ).rstrip("\n")

    # Run from an unrelated working directory to prove deck-relative resolution.
    with tempfile.TemporaryDirectory(prefix="nec2c-df-cwd-") as unrelated:
        status, output = invoke(binary, fixture, Path(unrelated))
    if status != 0:
        raise AssertionError(f"valid resource failed with status {status}\n{output}")
    if normalized_summary(output) != expected:
        raise AssertionError(
            "canonical summary mismatch\nEXPECTED:\n"
            + expected
            + "\nACTUAL:\n"
            + normalized_summary(output)
        )

    with tempfile.TemporaryDirectory(prefix="nec2c-df-tests-") as temporary:
        temp_root = Path(temporary)
        failures = [
            ("bad-record", BASIC.replace("BOX 1 0 0 0 1 1 1", "BOX 1 0 0"),
             "expected BOX <region-id>"),
            ("duplicate-id", BASIC.replace("REGION 1", "MATERIAL 1 ER 4\nREGION 1"),
             "duplicate material ID"),
            ("missing-end", BASIC.replace("END\n", ""), "required record is missing [value: END]"),
            ("non-finite", BASIC.replace("ER 3.0", "ER NaN"),
             "ER must be finite and greater than or equal to 1 [value: NaN]"),
            ("old-version", BASIC.replace("DIELECTRIC 1", "DIELECTRIC 2"),
             "unsupported format version; expected 1 [value: 2]"),
            ("unknown-keyword", BASIC.replace("MATERIAL 1", "SUBSTANCE 1"),
             "unknown record keyword [value: SUBSTANCE]"),
            ("unknown-material", BASIC.replace("REGION 1 MATERIAL 1", "REGION 1 MATERIAL 99"),
             "referenced material ID is not declared [value: 99]"),
            ("unknown-region", BASIC.replace("BOX 1", "BOX 99"),
             "referenced region ID is not declared [value: 99]"),
            ("inverted-box", BASIC.replace("BOX 1 0 0 0 1 1 1", "BOX 1 1 0 0 0 1 1"),
             "each minimum bound must be less than its maximum"),
            ("duplicate-mesh", BASIC.replace("END", "MESH MAX_KH 0.4\nEND"),
             "duplicate MESH record"),
            ("duplicate-mesh-option", BASIC.replace(
                "MESH MAX_KH 0.5", "MESH MAX_KH 0.5 MAX_KH 0.4"
             ), "duplicate mesh option [value: MAX_KH]"),
            ("after-end", BASIC + "MATERIAL 2 ER 4\n", "content is not permitted after END"),
        ]
        for name, resource_text, message in failures:
            expect_failure(binary, temp_root, name, resource_text, message)

        duplicate_root = temp_root / "duplicate-df"
        duplicate_deck, _ = write_case(
            duplicate_root,
            BASIC,
            ['DF "resources/model.diel"', 'DF "resources/model.diel"'],
        )
        status, output = invoke(binary, duplicate_deck, temp_root)
        if status == 0 or "A SECOND DF CARD IS NOT PERMITTED" not in output:
            raise AssertionError("second DF card was not rejected")

        quote_root = temp_root / "unterminated-df"
        quote_deck, _ = write_case(quote_root, BASIC, ['DF "resources/model.diel'])
        status, output = invoke(binary, quote_deck, temp_root)
        if status == 0 or "unterminated quoted resource path" not in output:
            raise AssertionError("unterminated DF quote was not rejected")

        long_root = temp_root / "overlong-df"
        long_deck, _ = write_case(long_root, BASIC, ["DF " + "x" * 140])
        status, output = invoke(binary, long_deck, temp_root)
        if status == 0 or "LINE EXCEEDS THE 132 CHARACTER LIMIT" not in output:
            raise AssertionError("overlong DF card was not rejected")

        missing_path_root = temp_root / "missing-df-path"
        missing_path_deck, _ = write_case(missing_path_root, BASIC, ["DF"])
        status, output = invoke(binary, missing_path_deck, temp_root)
        if status == 0 or "resource path is missing" not in output:
            raise AssertionError("missing DF resource path was not rejected")

        missing_file_root = temp_root / "missing-resource"
        missing_file_deck, missing_file = write_case(missing_file_root, BASIC)
        missing_file.unlink()
        status, output = invoke(binary, missing_file_deck, temp_root)
        if status == 0 or "cannot open dielectric resource" not in output:
            raise AssertionError("missing dielectric resource was not diagnosed")

        defaults_root = temp_root / "mesh-defaults"
        defaults_text = BASIC.replace(
            "MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2.0\n", ""
        )
        defaults_deck, _ = write_case(defaults_root, defaults_text)
        status, output = invoke(binary, defaults_deck, temp_root)
        if status != 0 or "MAX_ASPECT 2 (DEFAULTS)" not in output:
            raise AssertionError("omitted MESH record did not use documented defaults")

        absolute_root = temp_root / "absolute"
        absolute_deck, absolute_resource = write_case(absolute_root, BASIC)
        absolute_deck.write_text(
            f'CE\nGW 1 1 100 0 0 100 0 1 0.001\nGE\nFR 0 1 0 0 1.0\nDF "{absolute_resource}"\nEN\n',
            encoding="utf-8",
        )
        status, output = invoke(binary, absolute_deck, temp_root)
        if status != 0 or "PATH PORTABILITY: ABSOLUTE (NONPORTABLE)" not in output:
            raise AssertionError("absolute DF path was not accepted and reported")

        spaced_root = temp_root / "space-path"
        spaced_deck, spaced_resource = write_case(spaced_root, BASIC)
        renamed_resource = spaced_resource.with_name("model with space.diel")
        spaced_resource.rename(renamed_resource)
        spaced_deck.write_text(
            'CE\nGW 1 1 100 0 0 100 0 1 0.001\nGE\nFR 0 1 0 0 1.0\n'
            'DF "resources/model with space.diel"\nEN\n',
            encoding="utf-8",
        )
        status, output = invoke(binary, spaced_deck, temp_root)
        if status != 0 or str(renamed_resource) not in output:
            raise AssertionError("quoted deck-relative path containing spaces failed")

        lifecycle_root = temp_root / "nx-lifecycle"
        lifecycle_deck, _ = write_case(lifecycle_root, BASIC)
        lifecycle_deck.write_text(
            'CE\nGW 1 1 100 0 0 100 0 1 0.001\nGE\nFR 0 1 0 0 1.0\nDF "resources/model.diel"\nNX\n'
            'CE\nGW 1 1 100 0 0 100 0 1 0.001\nGE\nFR 0 1 0 0 1.0\nDF "resources/model.diel"\nEN\n',
            encoding="utf-8",
        )
        status, output = invoke(binary, lifecycle_deck, temp_root)
        if status != 0 or output.count("-------- DIELECTRIC RESOURCE --------") != 2:
            raise AssertionError("NX did not release dielectric state for the next job")

        gate_root = temp_root / "solver-gate"
        gate_deck, _ = write_case(gate_root, BASIC, trailer="XQ\nEN\n")
        status, output = invoke(binary, gate_deck, temp_root)
        if status == 0 or "WIRE MODE REQUIRES AN ORDINARY VOLTAGE-SOURCE EX CARD" not in output:
            raise AssertionError("dielectric solver gate was not enforced")

    print("PASS: canonical dielectric model summary")
    print("PASS: 23 DF/resource parsing, diagnostics, path, and lifecycle checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
