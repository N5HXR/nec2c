#!/usr/bin/env python3
"""M10D Burke/King/Miller buried radial-wire ground-screen qualification."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
FREQUENCY_MHZ = 299.792458
EPS0 = 8.8541878128e-12
OMEGA = 2.0 * math.pi * FREQUENCY_MHZ * 1.0e6
REFERENCE = {1.0: 3.088e-10, 10.0: 2.872e-10,
             100.0: 2.763e-8, 1000.0: 1.209e-5}
RADII = (0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60)
COUNTS = (5, 10, 20, 50, 100)


def deck(radials: int, screen_radius: float, loss: float,
         maximum_pitch: float) -> str:
    junction_cell = 0.005
    outer_segments = max(1, math.ceil(
        (screen_radius - junction_cell) / maximum_pitch - 1e-12))
    conductivity = loss * OMEGA * EPS0
    return "".join((
        "CM M10D BURKE 1984 BURIED RADIAL GROUND SCREEN\nCE\n",
        f"GW 10 1 0 0 0 {junction_cell:.9g} 0 "
        "-0.00001 0.000001\n",
        f"GW 20 {outer_segments} {junction_cell:.9g} 0 -0.00001 "
        f"{screen_radius:.9g} 0 -0.00001 0.000001\n",
        f"GR 100 {radials}\n",
        "GW 1 25 0 0 0 0 0 0.25 0.000001\n",
        "GE 0\n",
        f"GN 2 0 0 0 10 {conductivity:.17g}\n",
        f"FR 0 1 0 0 {FREQUENCY_MHZ:.9f} 0\n",
        "EX 0 1 1 0 1 0\n",
        "NE 0 1 1 1 1000 0 0 0 0 0\n",
        "EN\n",
    ))


def reference_deck(loss: float) -> str:
    conductivity = loss * OMEGA * EPS0
    return "".join((
        "CM M10D BURKE 1984 0.1-LAMBDA RCE REFERENCE DIPOLE\nCE\n",
        # Burke does not state this radius.  Match the documented M10B
        # inference, and retain the absolute-field uncertainty separately
        # from the same-path RCE normalization.
        "GW 1 31 0 0 0 0 0 .1 .0001\n",
        "GE 0\n",
        f"GN 2 0 0 0 10 {conductivity:.17g}\n",
        f"FR 0 1 0 0 {FREQUENCY_MHZ:.9f} 0\n",
        "EX 0 1 16 0 1 0\n",
        "NE 0 1 1 1 1000 0 0 0 0 0\n",
        "EN\n",
    ))


def run_model(directory: Path, cache: Path, stem: str, source_text: str,
              threads: int) -> tuple[dict, str, float]:
    source = directory / f"{stem}.nec"
    result = directory / f"{stem}.json"
    output = directory / f"{stem}.out"
    source.write_text(source_text, encoding="ascii")
    command = [str(BINARY), "--wire-engine", "nec4",
               "--ground-engine", "nec4", "--model-check", "warn",
               "--threads", str(threads), "--nec4-ground-model", "grid",
               "--nec4-ground-cache", str(cache), "--result-json",
               str(result), "-i", str(source), "-o", str(output)]
    started = time.monotonic()
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True, timeout=7200)
    elapsed = time.monotonic() - started
    report = output.read_text(errors="replace") if output.exists() else ""
    if completed.returncode:
        raise AssertionError((stem, completed.returncode, completed.stderr,
                              report[-5000:]))
    return json.loads(result.read_text(encoding="utf-8")), report, elapsed


def normalized_vertical_field(document: dict) -> tuple[float, float]:
    execution = document["jobs"][0]["executions"][0]
    field = execution["near_fields"][0]["samples"][0]["value"]["z"]
    power = execution["power"]["accepted_w"]
    normalized = (field["re"] ** 2 + field["im"] ** 2) / power
    return normalized, power


def run_reference(directory: Path, cache: Path, loss: float,
                  threads: int = 16) -> dict:
    stem = f"reference-l{loss:g}"
    document, report, elapsed = run_model(
        directory, cache, stem, reference_deck(loss), threads)
    normalized, power = normalized_vertical_field(document)
    warnings = [item for item in document.get("diagnostics", [])
                if item.get("code") == "NEC4ISH_M9_GRAZING_INTERFERENCE"]
    if not (power > 0.0 and math.isfinite(normalized)):
        raise AssertionError(("non-finite reference result", stem))
    return {
        "ground_imaginary_epsilon": loss,
        "published_v2_per_m2_w": REFERENCE[loss],
        "computed_v2_per_m2_w": normalized,
        "computed_minus_published_db": 10.0 * math.log10(
            normalized / REFERENCE[loss]),
        "reference_radius_lambda": 0.0001,
        "reference_radius_status": "inferred-report-omits-radius",
        "accepted_power_w": power,
        "grazing_warning_count": len(warnings),
        "wall_seconds": elapsed,
    }


def run(directory: Path, cache: Path, radials: int, screen_radius: float,
        loss: float, maximum_pitch: float, reference_normalized: float,
        threads: int = 16) -> dict:
    stem = (f"n{radials}-r{screen_radius:g}-l{loss:g}-"
            f"p{maximum_pitch:g}")
    document, report, elapsed = run_model(
        directory, cache, stem,
        deck(radials, screen_radius, loss, maximum_pitch), threads)
    junction = re.search(
        r"(\d+) INTERFACE-PRESCRIBED, (\d+) FALLBACKS", report)
    if not junction or tuple(map(int, junction.groups())) != (1, 0):
        raise AssertionError(("radial-screen junction prescription", stem,
                              junction.groups() if junction else None))
    execution = document["jobs"][0]["executions"][0]
    normalized, power = normalized_vertical_field(document)
    impedance = execution["feeds"][0]["impedance_ohm"]
    if not (power > 0.0 and math.isfinite(normalized) and
            math.isfinite(impedance["re"]) and
            math.isfinite(impedance["im"])):
        raise AssertionError(("non-finite radial-screen result", stem))
    outer_segments = max(1, math.ceil(
        (screen_radius - 0.005) / maximum_pitch - 1e-12))
    segments = 1 + outer_segments
    return {
        "radial_count": radials,
        "screen_radius_lambda": screen_radius,
        "radial_segments_each": segments,
        "interface_transition_cell_length_lambda": 0.005,
        "actual_outer_radial_pitch_lambda":
            (screen_radius - 0.005) / outer_segments,
        # RCE is explicitly a same-location, same-ground comparison.  Using
        # Burke's published absolute reference here mixed two propagation
        # implementations and exposed their common grazing-path scale error
        # as a false radial-screen error.
        "rce_db": 10.0 * math.log10(normalized / reference_normalized),
        "vertical_field_v2_per_m2_w": normalized,
        "impedance_ohm": impedance,
        "accepted_power_w": power,
        "wall_seconds": elapsed,
    }


def plot_evidence(evidence: dict, plot_path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True,
                                sharey=True)
    for axis, loss in zip(axes.flat, REFERENCE):
        for curve in evidence["curves"]:
            if -curve["ground_epsilon_relative"]["im"] != loss:
                continue
            axis.plot([row["screen_radius_lambda"]
                       for row in curve["samples"]],
                      [row["rce_db"] for row in curve["samples"]],
                      marker="o", linewidth=1.2,
                      label=f"N={curve['radial_count']}")
        axis.set_title(f"epsilon = 10 - j{loss:g}")
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=8)
    for axis in axes[-1, :]:
        axis.set_xlabel("Screen radius / lambda")
    for axis in axes[:, 0]:
        axis.set_ylabel("RCE (dB)")
    figure.suptitle("M10D buried radial-wire ground screen")
    figure.tight_layout()
    figure.savefig(plot_path, dpi=160)
    plt.close(figure)


def collect(output_path: Path, plot_path: Path, maximum_pitch: float,
            threads: int) -> None:
    evidence = {
        "schema": "nec2c-nec4-m10d-radial-screen",
        "schema_version": 2,
        "source": "Burke-King-Miller 1984 surface-wave report, Figure 5 and appendix",
        "frequency_mhz": FREQUENCY_MHZ,
        "monopole_height_lambda": 0.25,
        "wire_radius_lambda": 1e-6,
        "screen_depth_lambda": 1e-5,
        "observation": {"rho_lambda": 1000.0, "z_lambda": 0.0,
                        "component": "Ez"},
        "maximum_radial_pitch_lambda": maximum_pitch,
        "curves": [],
        "reference_controls": [],
        "refinement_checks": [],
    }
    partial_path = output_path.with_suffix(".partial.json")

    def checkpoint() -> None:
        partial_path.write_text(json.dumps(evidence, indent=2) + "\n",
                                encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="nec2c-m10d-") as name:
        directory = Path(name)
        cache = ROOT / "tmp/m10d/cache"
        cache.mkdir(parents=True, exist_ok=True)
        for loss in REFERENCE:
            print(f"M10D ground epsilon 10-j{loss:g}", flush=True)
            reference = run_reference(directory, cache, loss, threads)
            evidence["reference_controls"].append(reference)
            checkpoint()
            reference_normalized = reference["computed_v2_per_m2_w"]
            print("  reference: "
                  f"{reference['computed_minus_published_db']:+.3f} dB "
                  "vs published absolute field", flush=True)
            for count in COUNTS:
                # NECGS used rotational symmetry.  The general solver does
                # not yet reduce the axis-fixed monopole plus radial orbits,
                # so retain the complete 5/10/20-wire curves and published
                # 50/100-wire ordering anchors without turning qualification
                # into an hour-long dense-matrix endurance test per point.
                radii = (RADII if count == 5 else
                         RADII[:6] if count == 10 else
                         RADII[:6] if count == 20 else
                         (0.10,))
                samples = []
                for radius in radii:
                    sample = run(directory, cache, count, radius, loss,
                                 maximum_pitch, reference_normalized, threads)
                    samples.append(sample)
                    checkpoint()
                    print(f"  N={count:3d} r={radius:.2f}: "
                          f"{sample['rce_db']:.3f} dB, "
                          f"{sample['wall_seconds']:.2f}s", flush=True)
                evidence["curves"].append({
                    "ground_epsilon_relative": {"re": 10.0, "im": -loss},
                    "radial_count": count,
                    "samples": samples,
                })

        refinement_cases = ((5, 0.20, 1.0), (20, 0.20, 10.0),
                            (100, 0.10, 10.0))
        by_key = {
            (-curve["ground_epsilon_relative"]["im"],
             curve["radial_count"], row["screen_radius_lambda"]): row
            for curve in evidence["curves"] for row in curve["samples"]
        }
        fine_pitch = maximum_pitch / 2.0
        for count, radius, loss in refinement_cases:
            coarse = by_key[(loss, count, radius)]
            fine = run(directory, cache, count, radius, loss, fine_pitch,
                       next(row["computed_v2_per_m2_w"]
                            for row in evidence["reference_controls"]
                            if row["ground_imaginary_epsilon"] == loss),
                       threads)
            delta = fine["rce_db"] - coarse["rce_db"]
            evidence["refinement_checks"].append({
                "ground_imaginary_epsilon": loss,
                "radial_count": count,
                "screen_radius_lambda": radius,
                "coarse_pitch_lambda": maximum_pitch,
                "fine_pitch_lambda": fine_pitch,
                "coarse_rce_db": coarse["rce_db"],
                "fine_rce_db": fine["rce_db"],
                "fine_minus_coarse_db": delta,
            })
            checkpoint()
            if abs(delta) > 0.50:
                raise AssertionError(("radial-screen refinement", count,
                                      radius, loss, delta))

    output_path.write_text(json.dumps(evidence, indent=2) + "\n",
                           encoding="utf-8")
    partial_path.unlink(missing_ok=True)
    plot_evidence(evidence, plot_path)


def verify(path: Path) -> None:
    evidence = json.loads(path.read_text(encoding="utf-8"))
    if evidence.get("schema") != "nec2c-nec4-m10d-radial-screen":
        raise AssertionError("M10D evidence schema")
    curves = evidence["curves"]
    references = evidence.get("reference_controls", [])
    if len(references) != len(REFERENCE):
        raise AssertionError("M10D same-path reference coverage")
    if any(not math.isfinite(row["computed_v2_per_m2_w"])
           for row in references):
        raise AssertionError("M10D non-finite same-path reference")
    if len(curves) != len(REFERENCE) * len(COUNTS):
        raise AssertionError("M10D curve coverage")
    for curve in curves:
        if not curve["samples"]:
            raise AssertionError("empty M10D curve")
        if any(not math.isfinite(row["rce_db"]) for row in curve["samples"]):
            raise AssertionError("non-finite M10D curve")
    if max(abs(row["fine_minus_coarse_db"])
           for row in evidence["refinement_checks"]) > 0.50:
        raise AssertionError("M10D retained refinement gate")
    print("M10D retained radial-screen evidence passed")


def acceptance(path: Path, threads: int) -> None:
    evidence = json.loads(path.read_text(encoding="utf-8"))
    references = {
        float(row["ground_imaginary_epsilon"]): row
        for row in evidence["reference_controls"]
    }
    anchors = {}
    for curve in evidence["curves"]:
        loss = -float(curve["ground_epsilon_relative"]["im"])
        if curve["radial_count"] != 20:
            continue
        for row in curve["samples"]:
            if row["screen_radius_lambda"] == 0.1:
                anchors[loss] = row

    print("\nM10D.1 retained normalization diagnosis")
    print(" ground       reference offset    old mixed RCE    same-path RCE")
    for loss in REFERENCE:
        reference = references[loss]
        anchor = anchors[loss]
        print(f" 10-j{loss:<5g}   "
              f"{reference['computed_minus_published_db']:+9.3f} dB"
              f"       {anchor['rce_db_published_reference_original']:+8.3f} dB"
              f"       {anchor['rce_db']:+8.3f} dB")

    loss = 1000.0
    retained_reference = references[loss]
    retained_anchor = anchors[loss]
    with tempfile.TemporaryDirectory(prefix="nec2c-m10d1-acceptance-") as name:
        directory = Path(name)
        cache = ROOT / "tmp/m10d/acceptance-cache"
        cache.mkdir(parents=True, exist_ok=True)
        print("\nFresh 10-j1000 reference solve...", flush=True)
        reference = run_reference(directory, cache, loss, threads)
        print("Fresh N=20, radius=0.1 lambda screen solve...", flush=True)
        screen = run(directory, cache, 20, 0.1, loss, 0.02,
                     reference["computed_v2_per_m2_w"], threads)

    reference_delta = (reference["computed_minus_published_db"] -
                       retained_reference["computed_minus_published_db"])
    rce_delta = screen["rce_db"] - retained_anchor["rce_db"]
    print("\nFresh high-conductivity witness")
    print(" reference absolute-field offset: "
          f"{reference['computed_minus_published_db']:+.6f} dB")
    print(f" screen normalized field:          "
          f"{screen['vertical_field_v2_per_m2_w']:.12e} V^2/m^2/W")
    print(f" same-path RCE:                    {screen['rce_db']:+.6f} dB")
    print(f" retained-data RCE:                "
          f"{retained_anchor['rce_db']:+.6f} dB")
    print(f" fresh minus retained:             {rce_delta:+.6f} dB")
    if abs(reference_delta) > 0.01:
        raise AssertionError(("fresh reference changed", reference_delta))
    if abs(rce_delta) > 0.01:
        raise AssertionError(("fresh same-path RCE changed", rce_delta))
    if not 0.5 <= screen["rce_db"] <= 0.9:
        raise AssertionError(("Figure 5 high-conductivity band",
                              screen["rce_db"]))
    print("M10D.1 acceptance passed: the mixed-path error is removed and "
          "the fresh Figure 5 witness is in the published band")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--acceptance", action="store_true")
    parser.add_argument("--plot-existing", action="store_true")
    parser.add_argument("--pitch", type=float, default=0.02)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--output", type=Path, default=ROOT /
                        "tests/nec4-similarity/m10d-radial-screen.json")
    parser.add_argument("--plot", type=Path, default=ROOT /
                        "tests/nec4-similarity/m10d-radial-screen.png")
    arguments = parser.parse_args()
    if arguments.collect:
        collect(arguments.output, arguments.plot, arguments.pitch,
                arguments.threads)
    verify(arguments.output)
    if arguments.plot_existing:
        plot_evidence(json.loads(arguments.output.read_text(encoding="utf-8")),
                      arguments.plot)
    if arguments.acceptance:
        acceptance(arguments.output, arguments.threads)


if __name__ == "__main__":
    main()
