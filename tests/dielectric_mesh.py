#!/usr/bin/env python3
"""Acceptance tests for the M3 deterministic dielectric mesher and cost report."""

import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def invoke(binary, deck, cwd):
    output_path = cwd / ("nec2c-mesh-" + deck.stem + ".out")
    process = subprocess.run(
        [str(binary), "-i", str(deck), "-o", str(output_path)],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    output = output_path.read_text(encoding="utf-8") if output_path.exists() else ""
    return process.returncode, output + process.stdout + process.stderr


def write_model(root, resource, fr="FR 0 1 0 0 10.0", wire=None, name="model"):
    root.mkdir(parents=True)
    resource_path = root / f"{name}.diel"
    resource_path.write_text(resource, encoding="utf-8")
    if wire is None:
        wire = "GW 1 1 100 0 0 100 0 0.1 0.001"
    deck = root / f"{name}.nec"
    deck.write_text(
        f"CE\n{wire}\nGE\n{fr}\nDF \"{resource_path.name}\"\nEN\n",
        encoding="utf-8",
    )
    return deck


def mesh_section(output):
    match = re.search(
        r" {23}---------- DIELECTRIC MESH ----------.*?"
        r"STATUS: MESH AND COST VALIDATED; READY FOR DIELECTRIC SOLVER",
        output,
        re.DOTALL,
    )
    if match is None:
        raise AssertionError("dielectric mesh report is missing")
    return match.group(0)


def normalized_report(output):
    report = mesh_section(output)
    return re.sub(
        r"(?m)^(  (?:DENSE MATRIX|ESTIMATED WORKSPACE) BYTES:) \d+$",
        r"\1 <BYTES>",
        report,
    )


def expect_failure(binary, root, name, resource, message, fr="FR 0 1 0 0 10.0", wire=None):
    deck = write_model(root / name, resource, fr=fr, wire=wire)
    status, output = invoke(binary, deck, root)
    if status == 0 or message not in output:
        raise AssertionError(
            f"{name}: expected failure containing {message!r}; status={status}\n{output}"
        )


def base_resource(boxes, mesh="MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2"):
    return (
        "NEC2C-DIELECTRIC 1\n"
        "MATERIAL 1 ER 3\n"
        "REGION 1 MATERIAL 1\n"
        + boxes
        + mesh
        + "\nEND\n"
    )


def number_after(label, output):
    match = re.search(re.escape(label) + r"\s*([-+0-9.eE]+)", output)
    if match is None:
        raise AssertionError(f"missing numeric report field {label!r}")
    return float(match.group(1))


def main():
    binary = Path(os.environ.get("NEC2C", "./nec2c")).resolve()
    source_root = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
    fixture = source_root / "tests" / "dielectric" / "valid.nec"
    expected = (source_root / "tests" / "dielectric" / "expected-mesh.txt").read_text(
        encoding="utf-8"
    ).rstrip("\n")

    with tempfile.TemporaryDirectory(prefix="nec2c-mesh-fixture-") as temporary:
        status, output = invoke(binary, fixture, Path(temporary))
    if status != 0 or normalized_report(output) != expected:
        raise AssertionError("canonical M3 report mismatch\n" + normalized_report(output))

    with tempfile.TemporaryDirectory(prefix="nec2c-mesh-tests-") as temporary:
        root = Path(temporary)

        ordered = """NEC2C-DIELECTRIC 1
MATERIAL 2 ER 4
REGION 20 MATERIAL 2
BOX 20 0 0 0 0.12 0.08 0.04
BOX 20 0.12 0 0 0.20 0.08 0.04
MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2
END
"""
        reordered = """NEC2C-DIELECTRIC 1
BOX 20 0.12 0 0 0.20 0.08 0.04
BOX 20 0 0 0 0.12 0.08 0.04
REGION 20 MATERIAL 2
MESH MAX_ASPECT 2 MAX_CELLS 1000 MAX_KH 0.5
MATERIAL 2 ER 4
END
"""
        deck_a = write_model(root / "order-a", ordered, fr="FR 0 3 0 0 50 25")
        deck_b = write_model(root / "order-b", reordered, fr="FR 0 3 0 0 50 25")
        status_a, output_a = invoke(binary, deck_a, root)
        status_b, output_b = invoke(binary, deck_b, root)
        if status_a != 0 or status_b != 0 or mesh_section(output_a) != mesh_section(output_b):
            raise AssertionError("equivalent reordered resources produced different meshes")
        if "MESH FREQUENCY: 100 MHz" not in output_a:
            raise AssertionError("linear FR sweep did not mesh at its highest frequency")

        log_deck = write_model(root / "frequency-log", ordered, fr="FR 1 3 0 0 25 2")
        status, output = invoke(binary, log_deck, root)
        if status != 0 or "MESH FREQUENCY: 100 MHz" not in output:
            raise AssertionError("multiplicative FR sweep did not mesh at its highest frequency")

        half_q = base_resource("BOX 1 0 0 0 0.1 0.1 0.1\n")
        deck = write_model(root / "q-half", half_q)
        status, output = invoke(binary, deck, root)
        if status != 0 or "ACCURACY WARNING" in output:
            raise AssertionError("MAX_KH=0.5 boundary was not accepted without warning")

        one_q = half_q.replace("MAX_KH 0.5", "MAX_KH 1.0")
        deck = write_model(root / "q-one", one_q)
        status, output = invoke(binary, deck, root)
        if status != 0 or "ACCURACY WARNING: 0.5 < MAX_KH <= 1.0" not in output:
            raise AssertionError("MAX_KH=1 boundary did not produce the accuracy warning")

        expect_failure(
            binary,
            root,
            "q-too-large",
            half_q.replace("MAX_KH 0.5", "MAX_KH 1.0001"),
            "MAX_KH greater than 1 is not permitted",
        )
        expect_failure(
            binary,
            root,
            "cell-limit",
            base_resource(
                "BOX 1 0 0 0 1 1 1\n",
                "MESH MAX_KH 0.1 MAX_CELLS 10 MAX_ASPECT 2",
            ),
            "required cells exceed MAX_CELLS before cell allocation",
            fr="FR 0 1 0 0 100",
        )
        expect_failure(
            binary,
            root,
            "overlap",
            base_resource("BOX 1 0 0 0 1 1 1\nBOX 1 0.5 0 0 1.5 1 1\n"),
            "physical box interiors overlap",
            fr="FR 0 1 0 0 1",
        )

        cross_material = """NEC2C-DIELECTRIC 1
MATERIAL 1 ER 2
MATERIAL 2 ER 3
REGION 1 MATERIAL 1
REGION 2 MATERIAL 2
BOX 1 0 0 0 0.1 0.1 0.1
BOX 2 0.1 0 0 0.2 0.1 0.1
MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2
END
"""
        deck = write_model(root / "cross-material-touch", cross_material)
        status, output = invoke(binary, deck, root)
        assert status == 0 and number_after("ACTIVE CELLS:", output) == 2, output

        aspect_model = base_resource(
            "BOX 1 0 0 0 0.2 0.02 0.02\n",
            "MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 1.5",
        )
        deck = write_model(root / "aspect", aspect_model, fr="FR 0 1 0 0 1")
        status, output = invoke(binary, deck, root)
        if status != 0 or number_after("MAXIMUM CELL ASPECT:", output) > 1.5 + 1e-12:
            raise AssertionError("aspect-ratio refinement did not enforce MAX_ASPECT")
        if number_after("VOLUME RELATIVE ERROR:", output) > 1e-10:
            raise AssertionError("meshing did not conserve physical volume")

        vacuum = """NEC2C-DIELECTRIC 1
MATERIAL 1 ER 1
REGION 1 MATERIAL 1
BOX 1 0 0 0 1 2 3
END
"""
        deck = write_model(root / "vacuum", vacuum)
        status, output = invoke(binary, deck, root)
        if status != 0 or "ACTIVE CELLS: 0" not in output or "VACUUM VOLUME OMITTED: 6" not in output:
            raise AssertionError("vacuum material cells were not omitted")

        clearance_model = base_resource("BOX 1 0 0 0 0.1 0.1 0.1\n")
        expect_failure(
            binary,
            root,
            "clearance",
            clearance_model,
            "dielectric-to-conductor clearance is insufficient",
            wire="GW 1 1 0.05 0.05 -0.1 0.05 0.05 0.2 0.001",
        )

        patch_root = root / "patch-clearance"
        patch_root.mkdir()
        patch_resource = patch_root / "patch.diel"
        patch_resource.write_text(clearance_model, encoding="utf-8")
        patch_deck = patch_root / "patch.nec"
        patch_deck.write_text(
            "CE\nSM 1 1 0 0 0 1 0 0\nSC 0 0 1 1 0 0 1 0\nGE\n"
            "FR 0 1 0 0 10\nDF \"patch.diel\"\nEN\n",
            encoding="utf-8",
        )
        status, patch_output = invoke(binary, patch_deck, root)
        if status == 0 or "dielectric-to-conductor clearance is insufficient" not in patch_output:
            raise AssertionError("conducting-patch clearance violation was not rejected")

        q_report = re.search(
            r"Q MIN/MEDIAN/P95/MAX:\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
            r"([-+0-9.eE]+)\s+([-+0-9.eE]+)",
            output_a,
        )
        if q_report is None or float(q_report.group(4)) > 0.5 + 1e-9:
            raise AssertionError("emitted cell exceeds MAX_KH")
        if "DENSE MATRIX BYTES:" not in output_a or "TOTAL UNKNOWNS:" not in output_a:
            raise AssertionError("cost estimate is incomplete")

    print("PASS: canonical deterministic mesh and cost report")
    print("PASS: 14 M3 ordering, frequency, geometry, resolution, volume, and clearance checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
