#!/usr/bin/env python3
"""M10D closure: sensitive mesh convergence and appendix impedance anchors."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile

import m10d_radial_screen as m10d


ROOT = Path(os.environ.get(
    "TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
DEFAULT_OUTPUT = ROOT / "tests/nec4-similarity/m10d-closure.json"
DEFAULT_PLOT = ROOT / "tests/nec4-similarity/m10d-closure.png"
PITCHES = (0.02, 0.01, 0.005)
CONVERGENCE_CASES = (
    # Figure 5's low-loss sparse-screen standing-wave region.
    (1.0, 5, 0.30),
    (1.0, 5, 0.40),
    (1.0, 5, 0.50),
    # The original M10D refinement point nearest its 0.50 dB gate.
    (10.0, 20, 0.20),
)
APPENDIX_CASES = tuple(
    (loss, count, 0.10) for loss in (1.0, 10.0, 100.0)
    for count in (20, 50, 100)
)
PEC_MONOPOLE_OHM = complex(38.6, 22.2)
# Approximate NECGS symbol readings at screen radius 0.1 lambda.  Bounds are
# intentionally wider than line thickness and scan noise; these are an
# independent current-solution check, not data fitted to the Wait-Pope curves.
APPENDIX_DIGITIZED = {
    (1.0, 20): ((6.3, 1.3), (-16.0, 2.0)),
    (1.0, 50): ((5.5, 1.0), (-15.5, 1.5)),
    (1.0, 100): ((5.0, 1.0), (-15.5, 1.5)),
    (10.0, 20): ((15.0, 2.0), (-9.5, 1.5)),
    (10.0, 50): ((11.5, 1.5), (-11.0, 1.5)),
    (10.0, 100): ((10.0, 1.5), (-11.5, 1.5)),
    (100.0, 20): ((8.2, 0.8), (3.0, 1.0)),
    (100.0, 50): ((6.3, 0.7), (0.8, 0.8)),
    (100.0, 100): ((5.6, 0.6), (-0.8, 0.8)),
}


def complex_value(pair: dict) -> complex:
    return complex(pair["re"], pair["im"])


def retained_references() -> dict[float, float]:
    path = ROOT / "tests/nec4-similarity/m10d-radial-screen.json"
    evidence = json.loads(path.read_text(encoding="utf-8"))
    return {
        float(row["ground_imaginary_epsilon"]):
            float(row["computed_v2_per_m2_w"])
        for row in evidence["reference_controls"]
    }


def plot_evidence(evidence: dict, path: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    convergence, appendix = axes
    for case in evidence["convergence"]:
        loss = -case["ground_epsilon_relative"]["im"]
        label = (f"10-j{loss:g}, N={case['radial_count']}, "
                 f"r={case['screen_radius_lambda']:g}")
        rows = case["samples"]
        convergence.plot([row["radial_segments_each"] for row in rows],
                         [row["rce_db"] for row in rows], marker="o",
                         label=label)
    convergence.set_xlabel("Segments per radial")
    convergence.set_ylabel("Same-path RCE (dB)")
    convergence.set_title("Sensitive radial-mesh convergence")
    convergence.grid(True, alpha=0.25)
    convergence.legend(fontsize=8)

    colors = {1.0: "tab:blue", 10.0: "tab:orange", 100.0: "tab:green"}
    for component, marker in (("re", "o"), ("im", "s")):
        for row in evidence["appendix_impedance_anchors"]:
            loss = -row["ground_epsilon_relative"]["im"]
            target = row["digitized_necgs_symbol_ohm"][component]
            uncertainty = row["digitization_uncertainty_ohm"][component]
            measured = row["delta_impedance_ohm"][component]
            appendix.errorbar(target, measured, xerr=uncertainty,
                              color=colors[loss], marker=marker,
                              linestyle="none", capsize=2)
    limits = (-19.0, 18.0)
    appendix.plot(limits, limits, color="black", linewidth=1,
                  linestyle="--", label="exact agreement")
    appendix.set_xlim(limits)
    appendix.set_ylim(limits)
    appendix.set_aspect("equal", adjustable="box")
    appendix.set_xlabel("Digitized NECGS symbol (ohm)")
    appendix.set_ylabel("nec2c delta impedance (ohm)")
    appendix.set_title("Appendix A1-A3 anchors")
    appendix.grid(True, alpha=0.25)
    appendix.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def collect(output: Path, plot: Path, threads: int) -> None:
    references = retained_references()
    evidence = {
        "schema": "nec2c-nec4-m10d-closure",
        "schema_version": 1,
        "source": "Burke-King-Miller 1984 surface-wave report, Figure 5 and Appendix A",
        "pec_quarter_wave_monopole_impedance_ohm": {
            "re": PEC_MONOPOLE_OHM.real, "im": PEC_MONOPOLE_OHM.imag},
        "convergence": [],
        "appendix_impedance_anchors": [],
    }
    partial = output.with_suffix(".partial.json")

    def checkpoint() -> None:
        partial.write_text(json.dumps(evidence, indent=2) + "\n",
                           encoding="utf-8")

    cache = ROOT / "tmp/m10d/closure-cache"
    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nec2c-m10d-closure-") as name:
        directory = Path(name)
        completed: dict[tuple[float, int, float, float], dict] = {}

        def execute(loss: float, count: int, radius: float,
                    pitch: float) -> dict:
            key = (loss, count, radius, pitch)
            if key not in completed:
                completed[key] = m10d.run(
                    directory, cache, count, radius, loss, pitch,
                    references[loss], threads)
            return completed[key]

        for loss, count, radius in CONVERGENCE_CASES:
            print(f"closure convergence: 10-j{loss:g}, N={count}, "
                  f"r={radius:g}", flush=True)
            rows = []
            previous = None
            for pitch in PITCHES:
                row = execute(loss, count, radius, pitch)
                impedance = complex_value(row["impedance_ohm"])
                item = {
                    "maximum_pitch_lambda": pitch,
                    "radial_segments_each": row["radial_segments_each"],
                    "rce_db": row["rce_db"],
                    "impedance_ohm": row["impedance_ohm"],
                    "wall_seconds": row["wall_seconds"],
                }
                if previous is not None:
                    previous_impedance = complex_value(
                        previous["impedance_ohm"])
                    item["rce_change_from_previous_db"] = (
                        item["rce_db"] - previous["rce_db"])
                    item["impedance_change_from_previous_ohm"] = abs(
                        impedance - previous_impedance)
                rows.append(item)
                previous = item
                checkpoint()
                print(f"  pitch={pitch:g}: {row['rce_db']:+.6f} dB, "
                      f"Z={impedance.real:+.6f}{impedance.imag:+.6f}j, "
                      f"{row['wall_seconds']:.2f}s", flush=True)
            evidence["convergence"].append({
                "ground_epsilon_relative": {"re": 10.0, "im": -loss},
                "radial_count": count,
                "screen_radius_lambda": radius,
                "samples": rows,
            })
            checkpoint()

        for loss, count, radius in APPENDIX_CASES:
            print(f"appendix anchor: 10-j{loss:g}, N={count}, r={radius:g}",
                  flush=True)
            row = execute(loss, count, radius, 0.01)
            impedance = complex_value(row["impedance_ohm"])
            delta = impedance - PEC_MONOPOLE_OHM
            resistance, reactance = APPENDIX_DIGITIZED[(loss, count)]
            evidence["appendix_impedance_anchors"].append({
                "ground_epsilon_relative": {"re": 10.0, "im": -loss},
                "radial_count": count,
                "screen_radius_lambda": radius,
                "maximum_pitch_lambda": 0.01,
                "radial_segments_each": row["radial_segments_each"],
                "input_impedance_ohm": row["impedance_ohm"],
                "delta_impedance_ohm": {"re": delta.real, "im": delta.imag},
                "digitized_necgs_symbol_ohm": {
                    "re": resistance[0], "im": reactance[0]},
                "digitization_uncertainty_ohm": {
                    "re": resistance[1], "im": reactance[1]},
                "digitization_provenance": (
                    "Appendix A1-A3 scanned NEC symbols at radius 0.1 lambda; "
                    "bounds include line thickness and scan-reading error"),
                "rce_db": row["rce_db"],
                "wall_seconds": row["wall_seconds"],
            })
            checkpoint()
            print(f"  delta Z={delta.real:+.6f}{delta.imag:+.6f}j ohm",
                  flush=True)

    output.write_text(json.dumps(evidence, indent=2) + "\n",
                      encoding="utf-8")
    partial.unlink(missing_ok=True)
    plot_evidence(evidence, plot)


def verify(output: Path) -> None:
    evidence = json.loads(output.read_text(encoding="utf-8"))
    if evidence.get("schema") != "nec2c-nec4-m10d-closure":
        raise AssertionError("M10D closure evidence schema")
    if len(evidence["convergence"]) != len(CONVERGENCE_CASES):
        raise AssertionError("M10D closure convergence coverage")
    if len(evidence["appendix_impedance_anchors"]) != len(APPENDIX_CASES):
        raise AssertionError("M10D closure appendix coverage")
    for case in evidence["convergence"]:
        rows = case["samples"]
        if [row["maximum_pitch_lambda"] for row in rows] != list(PITCHES):
            raise AssertionError("M10D closure pitch sequence")
        if any(not math.isfinite(row["rce_db"]) for row in rows):
            raise AssertionError("M10D closure non-finite RCE")
        final = rows[-1]
        if abs(final["rce_change_from_previous_db"]) > 0.03:
            raise AssertionError(("M10D closure RCE convergence", case))
        if final["impedance_change_from_previous_ohm"] > 0.5:
            raise AssertionError(("M10D closure impedance convergence", case))
    for row in evidence["appendix_impedance_anchors"]:
        delta = row["delta_impedance_ohm"]
        if not (math.isfinite(delta["re"]) and math.isfinite(delta["im"])):
            raise AssertionError("M10D closure non-finite impedance")
        target = row["digitized_necgs_symbol_ohm"]
        uncertainty = row["digitization_uncertainty_ohm"]
        if abs(delta["re"] - target["re"]) > uncertainty["re"] or \
           abs(delta["im"] - target["im"]) > uncertainty["im"]:
            raise AssertionError(("M10D closure appendix symbol", row))
    print("M10D closure retained convergence and appendix anchors passed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--plot-existing", action="store_true")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    arguments = parser.parse_args()
    if arguments.collect:
        collect(arguments.output, arguments.plot, arguments.threads)
    verify(arguments.output)
    if arguments.plot_existing:
        plot_evidence(json.loads(arguments.output.read_text(encoding="utf-8")),
                      arguments.plot)


if __name__ == "__main__":
    main()
