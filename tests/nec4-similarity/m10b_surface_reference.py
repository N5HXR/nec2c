#!/usr/bin/env python3
"""M10B published one-sided surface-wave reference and grazing witness."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from nec4_m9_asymptotic import (  # noqa: E402
    second_order_cross_reference,
    second_order_pole,
    second_order_zv_interface_reference,
)

BINARY = ROOT / "nec2c"
FREQUENCY_MHZ = 299.792458  # one free-space metre wavelength
EPS0 = 8.8541878128e-12
OMEGA = 2.0 * math.pi * FREQUENCY_MHZ * 1.0e6
PUBLISHED = {
    1.0: 3.088e-10,
    10.0: 2.872e-10,
    100.0: 2.763e-8,
    1000.0: 1.209e-5,
    "perfect": 1.814e-4,
}


def deck(loss: float | str) -> str:
    if loss == "perfect":
        ground = "GN 1\n"
    else:
        conductivity = float(loss) * OMEGA * EPS0
        ground = f"GN 2 0 0 0 10 {conductivity:.17g}\n"
    return "".join((
        "CM BURKE 1984 0.1-LAMBDA SURFACE-WAVE REFERENCE DIPOLE\n",
        "CE\n",
        # The report omits the reference-dipole radius.  A=0.0001 lambda is
        # stated for the related ground-stake models and is the documented
        # inference used here; the residual uncertainty is reported.
        "GW 1 31 0 0 0 0 0 .1 .0001\n",
        "GE 0\n", ground,
        f"FR 0 1 0 0 {FREQUENCY_MHZ:.9f} 0\n",
        "EX 0 1 16 0 1 0\n",
        "NE 0 1 1 1 1000 0 0 0 0 0\n",
        "EN\n",
    ))


def run(directory: pathlib.Path, label: str, loss: float | str,
        engine: str, threads: str, cache: pathlib.Path | None = None) -> dict:
    model = directory / f"{label}.nec"
    result = directory / f"{label}.json"
    output = directory / f"{label}.out"
    model.write_text(deck(loss), encoding="ascii")
    command = [str(BINARY), "--wire-engine", engine,
               "--ground-engine", engine, "--threads", threads,
               "--result-json", str(result)]
    if engine == "nec4":
        command += ["--nec4-ground-model", "grid", "--nec4-ground-cache",
                    str(cache)]
    command += [f"-i{model}", f"-o{output}"]
    subprocess.run(command, cwd=ROOT, check=True, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, text=True)
    return json.loads(result.read_text(encoding="utf-8"))


def execution(result: dict) -> dict:
    return result["jobs"][0]["executions"][0]


def normalized_vertical_field(result: dict) -> float:
    item = execution(result)
    field = item["near_fields"][0]["samples"][0]["value"]["z"]
    return (field["re"] ** 2 + field["im"] ** 2) / item["power"]["accepted_w"]


def feed(result: dict) -> dict:
    return execution(result)["feeds"][0]


def pair(value: complex) -> list[float]:
    return [value.real, value.imag]


def main() -> None:
    rows = []
    upper_k = OMEGA / 299792458.0
    lower_k = upper_k * complex(10.0, -1.0) ** 0.5
    saddle_generic = second_order_cross_reference(
        FREQUENCY_MHZ * 1.0e6, lower_k, upper_k, 1000.0, 0.05, 0.0,
        precision_digits=80)["z_v"]
    saddle_j10 = second_order_zv_interface_reference(
        FREQUENCY_MHZ * 1.0e6, lower_k, upper_k, 1000.0, 0.05,
        precision_digits=80)
    if abs(saddle_generic / saddle_j10 - 1.0) > 2.0e-8:
        raise AssertionError(("Appendix J-10 interface limit",
                              saddle_generic, saddle_j10))
    pole = second_order_pole(
        FREQUENCY_MHZ * 1.0e6, lower_k, upper_k, 1000.0, 0.05, 0.0)["z_v"]
    with tempfile.TemporaryDirectory(prefix="nec2c-m10b-") as name:
        directory = pathlib.Path(name)
        for loss, published in PUBLISHED.items():
            result = run(directory, f"classic-{loss}", loss, "nec2", "1")
            measured = normalized_vertical_field(result)
            error = abs(measured / published - 1.0)
            rows.append({"ground_imaginary_epsilon": loss,
                         "published_v2_per_m2_w": published,
                         "classic_v2_per_m2_w": measured,
                         "relative_error": error})
            if error > 0.16:
                raise AssertionError(("published surface-wave table", rows[-1]))
        if sum(row["relative_error"] <= 0.05 for row in rows) < 4:
            raise AssertionError(("published surface-wave majority", rows))

        cache = directory / "cache"
        auto = run(directory, "nec4-auto", 1.0, "nec4", "auto", cache)
        serial = run(directory, "nec4-serial", 1.0, "nec4", "1", cache)
        a, s = execution(auto), execution(serial)
        if a != s:
            raise AssertionError("NEC4 one-sided reference changed with workers")
        warnings = [item for item in auto["diagnostics"]
                    if item["code"] == "NEC4_M9_GRAZING_INTERFERENCE"]
        if len(warnings) != 1:
            raise AssertionError(("grazing warning", auto["diagnostics"]))
        preflight = [item for item in auto["diagnostics"]
                     if item["code"] == "NEC4_M7_MEDIUM_PREFLIGHT"]
        if len(preflight) != 1 or \
           preflight[0]["context"]["one_sided_interface_ends"] != 1:
            raise AssertionError(("one-sided endpoint", preflight))

        classic = run(directory, "classic-nec4-pair", 1.0, "nec2", "1")
        cf, nf = feed(classic), feed(auto)
        zi = cf["impedance_ohm"]
        zn = nf["impedance_ohm"]
        impedance_delta = abs(complex(zn["re"], zn["im"]) /
                              complex(zi["re"], zi["im"]) - 1.0)
        if impedance_delta > 0.001:
            raise AssertionError(("reference input impedance", impedance_delta))
        nec4_field = normalized_vertical_field(auto)
        published = PUBLISHED[1.0]
        published_ratio = nec4_field / published
        # The point-field adapter must retain the direct/image cancellation
        # over this 1000-wavelength path.  The original NEC4 grazing-risk
        # diagnostic remains conservative, but it is not permission to emit
        # a numerically wrong field.
        if abs(published_ratio - 1.0) > 0.05:
            raise AssertionError(("NEC4 published surface field",
                                  published_ratio))

        report = {
            "schema": "nec2c-nec4-m10b-surface-reference",
            "source": "Burke-King-Miller 1984 surface-wave report, Table 1",
            "reference_radius_lambda": 0.0001,
            "reference_radius_status": "inferred-report-omits-radius",
            "appendix_j_decomposition": {
                "rho_lambda": 1000.0,
                "observer_z_lambda": 0.05,
                "second_order_saddle_generic": pair(saddle_generic),
                "second_order_saddle_j10": pair(saddle_j10),
                "second_order_pole": pair(pole),
                "generic_j10_relative_delta": abs(
                    saddle_generic / saddle_j10 - 1.0),
            },
            "classic_external_table": rows,
            "nec4_one_sided": {
                "input_impedance_relative_delta_vs_classic": impedance_delta,
                "vertical_field_v2_per_m2_w": nec4_field,
                "published_vertical_field_v2_per_m2_w": published,
                "field_ratio_to_published": published_ratio,
                "warning": warnings[0],
                "worker_identity": "exact",
            },
        }
        target = ROOT / "tests/nec4-similarity/m10b-surface-reference.json"
        target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("M10B published surface-wave controls and one-sided interface field passed")


if __name__ == "__main__":
    main()
