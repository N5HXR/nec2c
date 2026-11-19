#!/usr/bin/env python3
"""M10C published Burke/King/Miller monopole-on-ground-stake trends."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
FREQUENCY_MHZ = 299.792458
EPS0 = 8.8541878128e-12
OMEGA = 2.0 * math.pi * FREQUENCY_MHZ * 1.0e6
REFERENCE = {1.0: 3.088e-10, 10.0: 2.872e-10,
             100.0: 2.763e-8, 1000.0: 1.209e-5}


def deck(height: float, depth: float, loss: float,
         segment_length: float) -> str:
    fixed = 0.005
    upper_segments = round((height - fixed) / segment_length)
    if not math.isclose(upper_segments * segment_length + fixed, height,
                        abs_tol=1e-12):
        raise ValueError("height must be an integral number of segments")
    wires = [f"GW 1 1 0 0 0 0 0 {fixed:.9g} .0001\n",
             f"GW 3 {upper_segments} 0 0 {fixed:.9g} 0 0 "
             f"{height:.9g} .0001\n"]
    if depth:
        lower_segments = round((depth - fixed) / segment_length)
        if depth < fixed or not math.isclose(
                lower_segments * segment_length + fixed, depth,
                            abs_tol=1e-12):
            raise ValueError("depth must be an integral number of segments")
        wires.append(f"GW 2 1 0 0 0 0 0 {-fixed:.9g} .0001\n")
        if lower_segments:
            wires.append(
                f"GW 4 {lower_segments} 0 0 {-fixed:.9g} 0 0 "
                f"{-depth:.9g} .0001\n")
    conductivity = loss * OMEGA * EPS0
    return "".join((
        "CM M10C BURKE 1984 MONOPOLE ON GROUND STAKE\nCE\n",
        *wires,
        "GE 0\n",
        f"GN 2 0 0 0 10 {conductivity:.17g}\n",
        f"FR 0 1 0 0 {FREQUENCY_MHZ:.9f} 0\n",
        "EX 0 1 1 0 1 0\n",
        "NE 0 1 1 1 1000 0 0 0 0 0\n",
        "EN\n",
    ))


def run(directory: Path, cache: Path, height: float, depth: float,
        loss: float, segment_length: float) -> dict:
    stem = (f"h{height:g}-d{depth:g}-l{loss:g}-s{segment_length:g}")
    source = directory / f"{stem}.nec"
    result = directory / f"{stem}.json"
    output = directory / f"{stem}.out"
    source.write_text(deck(height, depth, loss, segment_length),
                      encoding="ascii")
    command = [str(BINARY), "--wire-engine", "nec4",
               "--ground-engine", "nec4", "--threads", "1",
               "--nec4-ground-model", "grid", "--nec4-ground-cache",
               str(cache), "--result-json", str(result),
               "-i", str(source), "-o", str(output)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True, timeout=1200)
    if completed.returncode:
        raise AssertionError((stem, completed.returncode, completed.stderr,
                              output.read_text(errors="replace")[-4000:]
                              if output.exists() else ""))
    execution = json.loads(result.read_text(encoding="utf-8"))[
        "jobs"][0]["executions"][0]
    field = execution["near_fields"][0]["samples"][0]["value"]["z"]
    normalized = ((field["re"] ** 2 + field["im"] ** 2) /
                  execution["power"]["accepted_w"])
    impedance = execution["feeds"][0]["impedance_ohm"]
    return {
        "depth_lambda": depth,
        "vertical_field_v2_per_m2_w": normalized,
        "absolute_delta_from_published_reference_db":
            10.0 * math.log10(normalized / REFERENCE[loss]),
        "rce_db": 10.0 * math.log10(normalized / REFERENCE[loss]),
        "impedance_ohm": impedance,
    }


def main() -> None:
    depths = [index * 0.005 for index in range(41)]
    evidence = {
        "schema": "nec2c-nec4-m10c-ground-stake",
        "schema_version": 1,
        "source": "Burke-King-Miller 1984 surface-wave report, Figures 3-4",
        "frequency_mhz": FREQUENCY_MHZ,
        "wire_radius_lambda": 0.0001,
        "observation": {"rho_lambda": 1000.0, "z_lambda": 0.0,
                        "component": "Ez"},
        "segment_length_lambda": 0.005,
        "curves": [],
        "segmentation_checks": [],
    }
    with tempfile.TemporaryDirectory(prefix="nec2c-m10c-") as name:
        directory = Path(name)
        cache = ROOT / "tmp/m10c/cache"
        cache.mkdir(parents=True, exist_ok=True)
        for loss in REFERENCE:
            print(f"M10C ground epsilon 10-j{loss:g}", flush=True)
            # Populate the material/frequency cache before concurrent readers
            # enter it.  Each solver process remains numerically single-threaded.
            run(directory, cache, 0.1, 0.0, loss, 0.005)
            for height in (0.1, 0.25):
                with ThreadPoolExecutor(max_workers=32) as executor:
                    samples = list(executor.map(
                        lambda depth: run(directory, cache, height, depth,
                                          loss, 0.005), depths))
                peak = max(samples, key=lambda item: item["rce_db"])
                null = min(samples[1:], key=lambda item: item["rce_db"])
                refractive_index = complex(10.0, -loss) ** 0.5
                expected_peak = 0.25 / refractive_index.real
                evidence["curves"].append({
                    "ground_epsilon_relative": {"re": 10.0, "im": -loss},
                    "monopole_height_lambda": height,
                    "expected_quarter_ground_wavelength_lambda": expected_peak,
                    "measured_peak": peak,
                    "measured_minimum_excluding_zero": null,
                    "samples": samples,
                })
        by_loss = {}
        for curve in evidence["curves"]:
            loss = -curve["ground_epsilon_relative"]["im"]
            by_loss.setdefault(loss, []).append(curve)
            if abs(curve["measured_peak"]["depth_lambda"] -
                   curve["expected_quarter_ground_wavelength_lambda"]) > 0.03:
                raise AssertionError(("ground-stake resonance location", curve))
            points = [curve["measured_peak"]]
            if loss == 1.0:
                expected_null = 2.0 * curve[
                    "expected_quarter_ground_wavelength_lambda"]
                if abs(curve["measured_minimum_excluding_zero"][
                        "depth_lambda"] - expected_null) > 0.01:
                    raise AssertionError(("ground-stake half-wave null", curve))
                points.append(curve["measured_minimum_excluding_zero"])
            for point in points:
                fine = run(directory, cache,
                           curve["monopole_height_lambda"],
                           point["depth_lambda"], loss, 0.0025)
                delta = fine["rce_db"] - point["rce_db"]
                evidence["segmentation_checks"].append({
                    "ground_imaginary_epsilon": loss,
                    "monopole_height_lambda": curve[
                        "monopole_height_lambda"],
                    "depth_lambda": point["depth_lambda"],
                    "coarse_segment_length_lambda": 0.005,
                    "fine_segment_length_lambda": 0.0025,
                    "coarse_rce_db": point["rce_db"],
                    "fine_rce_db": fine["rce_db"],
                    "fine_minus_coarse_db": delta,
                })
                # The very sharp low-loss null is intentionally retained in
                # this check and is more segmentation-sensitive than a peak.
                if abs(delta) > 0.15:
                    raise AssertionError(("ground-stake segmentation", delta))
        for loss, pair in by_loss.items():
            pair.sort(key=lambda item: item["monopole_height_lambda"])
            differences = [high["rce_db"] - low["rce_db"]
                           for low, high in zip(pair[0]["samples"],
                                                pair[1]["samples"])]
            if min(differences) < 2.0:
                raise AssertionError(("published height ordering", loss,
                                      min(differences)))
            pair[0]["height_separation_db"] = {
                "minimum": min(differences),
                "maximum": max(differences),
            }
        target = ROOT / "tests/nec4-similarity/m10c-ground-stake.json"
        target.write_text(json.dumps(evidence, indent=2) + "\n",
                          encoding="utf-8")
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
        for axis, loss in zip(axes.flat, REFERENCE):
            for curve in by_loss[loss]:
                axis.plot([row["depth_lambda"] for row in curve["samples"]],
                          [row["rce_db"] for row in curve["samples"]],
                          marker=".", linewidth=1.2,
                          label=f"h={curve['monopole_height_lambda']}λ")
            axis.axvline(by_loss[loss][0][
                "expected_quarter_ground_wavelength_lambda"],
                color="0.55", linestyle="--", linewidth=0.8,
                label="¼ λ in ground")
            axis.set_title(f"εr = 10 − j{loss:g}")
            axis.grid(True, alpha=0.25)
            axis.legend(fontsize=8)
        for axis in axes[-1, :]:
            axis.set_xlabel("Stake depth d / λ₀")
        for axis in axes[:, 0]:
            axis.set_ylabel("RCE (dB)")
        figure.suptitle("M10C monopole-on-ground-stake reproduction")
        figure.tight_layout()
        figure.savefig(ROOT / "tests/nec4-similarity/"
                       "m10c-ground-stake.png", dpi=160)
        plt.close(figure)
    print("M10C published monopole/ground-stake sweep completed")


if __name__ == "__main__":
    main()
