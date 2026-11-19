#!/usr/bin/env python3
"""M9 Gate C2: lock the equation (5-10) interpretation and limits."""

import json
import os
from pathlib import Path


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()


def main():
    directory = ROOT / "tests" / "nec4-similarity"
    source = json.loads((directory / "m9-model-fit-er16.json").read_text(
        encoding="utf-8"))
    report = json.loads((directory / "m9-model-selection.json").read_text(
        encoding="utf-8"))
    if source["schema"] != "nec2c-nec4-m9-model-fit":
        raise AssertionError("schema")
    if report["schema"] != "nec2c-nec4-m9-model-selection":
        raise AssertionError("selection schema")
    if report["training_rows"] != 88:
        raise AssertionError("documented sample count")
    if report["calibration_points"] != 9:
        raise AssertionError("calibration sample count")
    ordinary = report["phase_modes"]["ordinary-ray"]
    if ordinary["z_v"]["full_rank"]["training_rms"] > 0.04:
        raise AssertionError("guide's dielectric Ez witness was not reproduced")
    for literal in ("literal-phase-root", "literal-exponential-root"):
        if (report["phase_modes"][literal]["z_v"]["full_rank"]["training_rms"]
                < 0.30):
            raise AssertionError("ambiguous literal phase reading was not rejected")
    selected = [component["selected"] for component in ordinary.values()]
    if max(value["training_rms"] for value in selected) > 0.13:
        raise AssertionError("rank-selected training envelope")
    if max(value["calibration_rms"] for value in selected) > 0.09:
        raise AssertionError("rank-selected calibration envelope")
    if not any(component["full_rank"]["calibration_rms"] > 1
               for component in ordinary.values()):
        raise AssertionError("full-rank instability is no longer exposed")
    if report["production_status"] != "characterization-only":
        raise AssertionError("an unvalidated model was marked production ready")
    print("M9 Gate C2 reproduces the guide witness and rejects unstable full-rank extrapolation")
    print("rank-selected ordinary-ray envelope: <=13% training, <=9% calibration RMS")


if __name__ == "__main__":
    main()
