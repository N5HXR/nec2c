#!/usr/bin/env python3
"""M8 Gate D2: shipped qualification artifact integrity."""

import csv
import json
import os
from pathlib import Path


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()


def main():
    directory = ROOT / "tests" / "nec4-similarity"
    evidence = json.loads((directory / "m8-qualification.json").read_text())
    rows = evidence["rows"]
    if evidence["schema"] != "nec2c-nec4-m8-qualification":
        raise AssertionError("qualification schema")
    if len(rows) != 65:
        raise AssertionError(("qualification row count", len(rows)))
    worst = max(float(row["worst_refinement_error"]) for row in rows)
    if worst > float(evidence["target"]):
        raise AssertionError(("qualification target", worst, evidence["target"]))
    materials = {row["material"] for row in rows}
    if len(materials) != 5:
        raise AssertionError(("qualification material count", materials))
    with (directory / "m8-qualification.csv").open(newline="") as stream:
        csv_rows = list(csv.DictReader(stream))
    if len(csv_rows) != len(rows):
        raise AssertionError("CSV/JSON row mismatch")
    figure = ROOT / "docs" / "nec4-m8-convergence.png"
    if figure.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError("qualification figure is not PNG")
    print(f"M8 Gate D2 qualification artifacts passed: {len(rows)} points, "
          f"5 media, worst {worst:.3e}")


if __name__ == "__main__":
    main()
