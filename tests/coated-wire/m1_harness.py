#!/usr/bin/env python3
"""Fast, dependency-free checks for the M1 OpenEMS reference harness."""

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "coated_wire_openems.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("coated_wire_openems", TOOL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    tool = load_tool()
    rows = [
        {"frequency_mhz": 290.0,
         "impedance_ohm": {"re": 70.0, "im": -10.0}},
        {"frequency_mhz": 300.0,
         "impedance_ohm": {"re": 80.0, "im": 30.0}},
    ]
    assert tool.resonance_from_rows(rows) == 292.5
    try:
        tool.resonance_from_rows([
            {"frequency_mhz": 290.0,
             "impedance_ohm": {"re": 0.0, "im": 0.0}},
            {"frequency_mhz": 300.0,
             "impedance_ohm": {"re": 0.0, "im": 0.0}},
        ])
    except ValueError as error:
        assert "degenerate zero terminal impedance" in str(error)
    else:
        raise AssertionError("shorted feed response was accepted")

    parser = tool.parser()
    args = parser.parse_args([
        "--coordinates", "cylindrical",
        "--cylindrical-source", "ring-voltage",
        "--material", "coated",
        "--coating-end", "hemisphere",
        "--axisymmetric-nf2ff-samples", "512",
    ])
    tool.validate(args, parser)
    assert args.conductor_radius_mm < args.coating_radius_mm
    assert args.coating_end == "hemisphere"
    assert args.axisymmetric_nf2ff_samples == 512

    source = TOOL.read_text()
    assert 'source.SetWeightFunction(["0", "0"' in source
    assert "directions=[False, True, False, False, True, True]" in source
    assert "(b0 - a0) / 4" in source
    assert "axisymmetric_nf2ff_data" in source
    print("M1 OpenEMS harness invariants passed")


if __name__ == "__main__":
    main()
