#!/usr/bin/env python3
"""M6 end-to-end validation for voltage-driven wire/dielectric coupling."""

import os
import pathlib
import re
import subprocess
import tempfile


NEC2C = pathlib.Path(os.environ.get("NEC2C", "./nec2c")).resolve()
SOURCE = pathlib.Path(os.environ.get("TEST_SRCDIR", ".")).resolve()


def run(deck: pathlib.Path, expect_success: bool = True) -> str:
    with tempfile.NamedTemporaryFile(prefix="nec2c-m6-output-", suffix=".out") as temporary:
        output = pathlib.Path(temporary.name)
        completed = subprocess.run(
            [str(NEC2C), "-i", str(deck), "-o", str(output)],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        result = output.read_text()
        if expect_success and completed.returncode != 0:
            raise AssertionError(f"nec2c failed for {deck}: {completed.stderr}\n{result}")
        if not expect_success and completed.returncode == 0:
            raise AssertionError(f"nec2c unexpectedly accepted {deck}\n{result}")
        return result


def coupled_impedance(text: str) -> complex:
    match = re.search(
        r"COUPLED ANTENNA INPUT PARAMETERS.*?\n\s*TAG SEG.*?\n"
        r"\s*\d+\s+\d+\s+[-+0-9.E]+\s+[-+0-9.E]+\s+"
        r"[-+0-9.E]+\s+[-+0-9.E]+\s+([-+0-9.E]+)\s+([-+0-9.E]+)",
        text, re.S,
    )
    if not match:
        raise AssertionError("missing coupled feed impedance")
    return complex(float(match.group(1)), float(match.group(2)))


def legacy_impedance(text: str) -> complex:
    match = re.search(
        r"ANTENNA INPUT PARAMETERS.*?\n\s*TAG.*?\n\s*No:.*?\n"
        r"\s*\d+\s+\d+\s+(?:[-+0-9.E]+\s+){4}([-+0-9.E]+)\s+([-+0-9.E]+)",
        text, re.S,
    )
    if not match:
        raise AssertionError("missing legacy feed impedance")
    return complex(float(match.group(1)), float(match.group(2)))


def residual(text: str) -> float:
    match = re.search(r"RELATIVE RESIDUAL:\s*([-+0-9.E]+)", text)
    if not match:
        raise AssertionError("missing coupled residual")
    return float(match.group(1))


def coupled_current(text: str, segment: int) -> complex:
    section = text.split("COUPLED WIRE CURRENTS", 1)[1].split(
        "COUPLED ANTENNA INPUT PARAMETERS", 1
    )[0]
    match = re.search(
        rf"^\s*{segment}\s+\d+\s+(?:[-+0-9.E]+\s+){{3}}"
        r"([-+0-9.E]+)\s+([-+0-9.E]+)$", section, re.M,
    )
    if not match:
        raise AssertionError(f"missing coupled current for segment {segment}")
    return complex(float(match.group(1)), float(match.group(2)))


def legacy_current(text: str, segment: int) -> complex:
    section = text.split("CURRENTS AND LOCATION", 1)[1]
    match = re.search(
        rf"^\s*{segment}\s+\d+\s+(?:[-+0-9.E]+\s+){{4}}"
        r"([-+0-9.E]+)\s+([-+0-9.E]+)", section, re.M,
    )
    if not match:
        raise AssertionError(f"missing legacy current for segment {segment}")
    return complex(float(match.group(1)), float(match.group(2)))


def main() -> None:
    root = SOURCE / "tests/dielectric"
    legacy_text = run(root / "m6-legacy.nec")
    legacy = legacy_impedance(legacy_text)
    vacuum_text = run(root / "m6-vacuum.nec")
    vacuum = coupled_impedance(vacuum_text)
    assert "VOLUME UNKNOWNS:            0" in vacuum_text
    assert abs(vacuum / legacy - 1.0) < 3.0e-5, (vacuum, legacy)
    assert residual(vacuum_text) <= 1.0e-10
    for segment in range(1, 12):
        reference = legacy_current(legacy_text, segment)
        actual = coupled_current(vacuum_text, segment)
        assert abs(actual - reference) <= max(8.0e-8, 5.0e-5*abs(reference))

    near_text = run(root / "m6-cube.nec")
    fine_text = run(root / "m6-cube-fine.nec")
    far_text = run(root / "m6-cube-far.nec")
    near = coupled_impedance(near_text)
    fine = coupled_impedance(fine_text)
    far = coupled_impedance(far_text)
    near_effect = abs(near - vacuum)
    far_effect = abs(far - vacuum)
    assert near_effect > 0.02
    assert far_effect < 0.8 * near_effect, (near_effect, far_effect)
    assert abs(fine - near) < 0.02 * near_effect, (near, fine)
    assert max(residual(near_text), residual(fine_text), residual(far_text)) <= 1.0e-10

    vacuum_21 = coupled_impedance(run(root / "m6-vacuum-21.nec"))
    loaded_21 = coupled_impedance(run(root / "m6-cube-21.nec"))
    refined_effect = loaded_21 - vacuum_21
    assert abs(refined_effect - (near - vacuum)) < 0.02 * near_effect

    with tempfile.TemporaryDirectory(prefix="nec2c-m6-gates-") as temporary:
        temp = pathlib.Path(temporary)
        resource = temp / "cube.diel"
        resource.write_text((root / "m6-cube.diel").read_text())
        deck = temp / "pattern.nec"
        deck.write_text(
            "CE\nGW 1 11 0 0 -0.24 0 0 0.24 0.001\nGE\n"
            "FR 0 1 0 0 299.792458\nEX 0 1 6 0 1 0\n"
            "DF \"cube.diel\"\nRP 0 1 1 0000 90 0 0 0 0 0\nEN\n"
        )
        rejected = run(deck, expect_success=False)
        assert "DIELECTRIC CONTRIBUTIONS TO RADIATION PATTERNS ARE SCHEDULED FOR M8" in rejected

        reciprocal_geometry = (
            "CE\nGW 1 11 0 0 -0.24 0 0 0.24 0.001\n"
            "GW 2 11 0.5 0 -0.24 0.5 0 0.24 0.001\nGE\n"
            "FR 0 1 0 0 299.792458\n"
        )
        first = temp / "reciprocal-first.nec"
        second = temp / "reciprocal-second.nec"
        first.write_text(reciprocal_geometry +
                         "EX 0 1 6 0 1 0\nDF \"cube.diel\"\nXQ\nEN\n")
        second.write_text(reciprocal_geometry +
                          "EX 0 2 6 0 1 0\nDF \"cube.diel\"\nXQ\nEN\n")
        transfer_21 = coupled_current(run(first), 17)
        transfer_12 = coupled_current(run(second), 6)
        reciprocity_error = abs(transfer_21 - transfer_12) / max(
            abs(transfer_21), abs(transfer_12)
        )
        # This is a full solved transfer-admittance check, including the
        # collocation/Galerkin mismatch, rather than the 1e-9 dyadic-kernel
        # identity checked by the lower-level VIE suite.
        assert reciprocity_error <= 1.0e-5, reciprocity_error

    print("M6 wire/dielectric coupling validation passed")


if __name__ == "__main__":
    main()
