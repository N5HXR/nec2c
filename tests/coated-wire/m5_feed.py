#!/usr/bin/env python3
"""M5 acceptance checks for the selected coated-wire feed model."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def main():
    acceptance = json.loads((ROOT / "tests/coated-wire/m5-acceptance.json").read_text(
        encoding="utf-8"))
    gap = acceptance["gap_convergence"]
    comparison = acceptance["nec2c_vs_physical_gap"]
    if abs(gap["resonance_percent"]) > 0.05:
        raise AssertionError(gap)
    if abs(gap["resistance_ohm"]) > 0.1:
        raise AssertionError(gap)
    if gap["maximum_relative_current_away_from_feed"] > 0.02:
        raise AssertionError(gap)
    if gap["maximum_pattern_db_above_minus_40"] > 0.05:
        raise AssertionError(gap)
    if abs(comparison["resonance_percent_vs_gap2"]) > 1.0:
        raise AssertionError(comparison)
    if abs(comparison["resistance_ohm_vs_gap2"]) > 2.0:
        raise AssertionError(comparison)
    if comparison["maximum_relative_current_vs_gap2_away_from_feed"] > 0.05:
        raise AssertionError(comparison)
    reference = acceptance["reference_plane"]
    if reference["deembedding_distance_mm"] != 0.0:
        raise AssertionError(reference)

    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m5-") as temporary:
        result = Path(temporary) / "result.json"
        report = Path(temporary) / "report.out"
        completed = subprocess.run([
            str(BINARY), "--backend", "builtin", "--jobs", "1",
            "--threads", "1", "--result-json", str(result),
            "-i", str(ROOT / "tests/coated-wire/m5-continuous-source.nec"),
            "-o", str(report),
        ], cwd=ROOT, text=True, capture_output=True)
        if completed.returncode:
            raise RuntimeError(completed.stdout+completed.stderr)
        parsed = json.loads(result.read_text(encoding="utf-8"))
    execution = parsed["jobs"][0]["executions"][0]
    feed = execution["feeds"][0]
    if feed["tag"] != 1 or feed["segment"] != 48:
        raise AssertionError(feed)
    currents = execution["wire_currents"]
    if len(currents) != 95:
        raise AssertionError("continuous article did not retain all 95 currents")
    left = currents[46]["current_a"]
    right = currents[48]["current_a"]
    if not (math.isclose(left["re"], right["re"], rel_tol=1e-12, abs_tol=1e-14)
            and math.isclose(left["im"], right["im"], rel_tol=1e-12,
                             abs_tol=1e-14)):
        raise AssertionError("center excitation broke dipole current symmetry")

    print("M5 1-2 mm finite-gap resonance, impedance, current, and pattern convergence passed")
    print(f"M5 gap delta: {gap['resonance_percent']:.6f}% resonance, "
          f"{gap['resistance_ohm']:+.6f} ohm")
    print(f"M5 NEC/physical-gap: {comparison['resonance_percent_vs_gap2']:+.4f}% resonance, "
          f"{comparison['resistance_ohm_vs_gap2']:+.6f} ohm")
    print("M5 continuous center-segment production source and zero-distance reference plane passed")


if __name__ == "__main__":
    main()
