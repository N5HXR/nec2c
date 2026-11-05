#!/usr/bin/env python3
"""M3 acceptance tests for the lossless isolated coated-wire operator."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def dipole(segments=95, frequency=276.65, er=None, outer=0.004,
           source=None, sweep=False, pattern=False, extended=False):
    if source is None:
        source = (segments + 1) // 2
    lines = [
        "CM M3 LOSSLESS COATED-WIRE TEST", "CE",
        f"GW 1 {segments} 0 0 -0.2375 0 0 0.2375 0.002", "GE 0",
    ]
    if er is not None:
        lines.append(f"LD 7 1 0 0 {er:.17g} {outer:.17g} 0 0")
    if sweep:
        lines.append("FR 0 141 0 0 250 0.5")
    else:
        lines.append(f"FR 0 1 0 0 {frequency:.17g} 0")
    lines.append(f"EX 0 1 {source} 0 1 0")
    if extended:
        lines.append("EK")
    lines.append("RP 0 37 1 1000 0 0 5 0 0 0" if pattern else "XQ 0")
    lines.append("EN")
    return "\n".join(lines) + "\n"


def run(directory, name, text, backend="builtin", threads=1):
    source = directory / f"{name}.nec"
    report = directory / f"{name}.out"
    result = directory / f"{name}.json"
    source.write_text(text, encoding="ascii")
    command = [str(BINARY), "--backend", backend, "--jobs", "1",
               "--threads", str(threads), "--result-json", str(result),
               "-i", str(source), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    parsed = json.loads(result.read_text()) if result.exists() else None
    if completed.returncode != 0:
        code = None
        if parsed and parsed.get("diagnostics"):
            code = parsed["diagnostics"][0]["code"]
        raise RuntimeError(f"{name} failed ({completed.returncode}, {code})\n"
                           f"{completed.stdout}\n{completed.stderr}")
    return parsed


def expect_error(directory, name, text, code):
    source = directory / f"{name}.nec"
    report = directory / f"{name}.out"
    result = directory / f"{name}.json"
    source.write_text(text, encoding="ascii")
    completed = subprocess.run(
        [str(BINARY), "--result-json", str(result), "-i", str(source),
         "-o", str(report)], cwd=ROOT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True)
    parsed = json.loads(result.read_text()) if result.exists() else None
    observed = (parsed["diagnostics"][0]["code"] if parsed and
                parsed.get("diagnostics") else None)
    if completed.returncode == 0 or observed != code:
        raise AssertionError(f"{name}: expected {code}, got {observed}")


def without_timing(value):
    if isinstance(value, dict):
        return {key: without_timing(item) for key, item in value.items()
                if key != "elapsed_seconds"}
    if isinstance(value, list):
        return [without_timing(item) for item in value]
    return value


def numerically_equal(left, right, path="$"):
    if isinstance(left, bool) or left is None or isinstance(left, str):
        if left != right:
            raise AssertionError((path, left, right))
    elif isinstance(left, (int, float)):
        if not isinstance(right, (int, float)) or not math.isclose(
                left, right, rel_tol=1e-8, abs_tol=1e-10):
            raise AssertionError((path, left, right))
    elif isinstance(left, list):
        if not isinstance(right, list) or len(left) != len(right):
            raise AssertionError(path)
        for index, (a, b) in enumerate(zip(left, right)):
            numerically_equal(a, b, f"{path}[{index}]")
    else:
        if not isinstance(right, dict) or left.keys() != right.keys():
            raise AssertionError(path)
        for key in left:
            if path.endswith(".solver") and key in {
                    "requested", "used", "elapsed_seconds"}:
                continue
            numerically_equal(left[key], right[key], f"{path}.{key}")


def execution(result):
    return without_timing(result["jobs"][0]["executions"])


def feed(execution_row):
    return execution_row["feeds"][0]["impedance_ohm"]


def resonance(result):
    rows = []
    for item in result["jobs"][0]["executions"]:
        impedance = feed(item)
        rows.append((item["frequency_mhz"], impedance["re"], impedance["im"]))
    for left, right in zip(rows, rows[1:]):
        if left[2] <= 0 <= right[2] or right[2] <= 0 <= left[2]:
            fraction = -left[2] / (right[2] - left[2])
            return (left[0] + fraction*(right[0]-left[0]),
                    left[1] + fraction*(right[1]-left[1]))
    raise AssertionError("resonance was not bracketed")


def current(result, segment):
    row = result["jobs"][0]["executions"][0]["wire_currents"][segment-1]
    assert row["segment"] == segment
    value = row["current_a"]
    return complex(value["re"], value["im"])


def normalized_pattern(result):
    samples = result["jobs"][0]["executions"][0]["far_fields"][0]["samples"]
    magnitudes = {}
    for sample in samples:
        theta = sample["theta_deg"]
        et = sample["r_e_theta"]
        ep = sample["r_e_phi"]
        magnitudes[theta] = math.sqrt(et["re"]**2 + et["im"]**2 +
                                      ep["re"]**2 + ep["im"]**2)
    peak = max(magnitudes.values())
    return {theta: (20*math.log10(value/peak) if value else -math.inf)
            for theta, value in magnitudes.items()}


def sampled_directivity(result):
    pattern = normalized_pattern(result)
    values = [(math.radians(theta), 10**(db/10))
              for theta, db in sorted(pattern.items())]
    integral = sum((right_theta-left_theta)*0.5*
                   (left_power*math.sin(left_theta) +
                    right_power*math.sin(right_theta))
                   for (left_theta, left_power), (right_theta, right_power)
                   in zip(values, values[1:]))
    return 10*math.log10(2/integral)


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m3-") as temporary:
        root = Path(temporary)

        bare = run(root, "bare", dipole(er=None, pattern=True))
        er_one = run(root, "er-one", dipole(er=1.0, pattern=True))
        zero_thickness = run(root, "zero-thickness",
                             dipole(er=3.2, outer=0.002, pattern=True))
        assert execution(bare) == execution(er_one)
        assert execution(bare) == execution(zero_thickness)
        bare_extended = run(root, "bare-extended",
                            dipole(er=None, extended=True))
        er_one_extended = run(root, "er-one-extended",
                              dipole(er=1.0, extended=True))
        zero_extended = run(root, "zero-extended",
                            dipole(er=3.2, outer=0.002, extended=True))
        assert execution(bare_extended) == execution(er_one_extended)
        assert execution(bare_extended) == execution(zero_extended)

        thin_1 = run(root, "thin-1e-4",
                     dipole(er=3.2, outer=0.002*(1+1e-4)))
        thin_2 = run(root, "thin-1e-6",
                     dipole(er=3.2, outer=0.002*(1+1e-6)))
        bare_z = feed(bare["jobs"][0]["executions"][0])
        thin_1_z = feed(thin_1["jobs"][0]["executions"][0])
        thin_2_z = feed(thin_2["jobs"][0]["executions"][0])
        delta_1 = complex(thin_1_z["re"]-bare_z["re"],
                          thin_1_z["im"]-bare_z["im"])
        delta_2 = complex(thin_2_z["re"]-bare_z["re"],
                          thin_2_z["im"]-bare_z["im"])
        if not math.isclose(abs(delta_1/delta_2), 100.0, rel_tol=0.002):
            raise AssertionError("thin-coating limit is not first-order stable")

        coated_1 = run(root, "coated-1t", dipole(er=3.2, pattern=True),
                       threads=1)
        coated_4 = run(root, "coated-4t", dipole(er=3.2, pattern=True),
                       threads=4)
        coated_lapack = run(root, "coated-lapack", dipole(er=3.2, pattern=True),
                            backend="lapack", threads=1)
        assert execution(coated_1) == execution(coated_4)
        numerically_equal(coated_1["jobs"][0]["executions"],
                          coated_lapack["jobs"][0]["executions"])
        assert coated_1["jobs"][0]["coated_wires"]["operator_status"] == \
            "lossless-isolated"

        reciprocal_left = run(root, "reciprocal-left",
                              dipole(er=3.2, source=24))
        reciprocal_right = run(root, "reciprocal-right",
                               dipole(er=3.2, source=72))
        transfer_left = current(reciprocal_left, 72)
        transfer_right = current(reciprocal_right, 24)
        reciprocal_scale = max(abs(transfer_left), abs(transfer_right), 1e-30)
        if abs(transfer_left-transfer_right)/reciprocal_scale > 1e-12:
            raise AssertionError("coated operator violates reciprocity")

        bare_sweep = run(root, "bare-sweep", dipole(er=None, sweep=True))
        coated_sweep = run(root, "coated-sweep", dipole(er=3.2, sweep=True))
        coarse_sweep = run(root, "coarse-sweep",
                           dipole(segments=47, er=3.2, sweep=True))
        medium_sweep = run(root, "medium-sweep",
                           dipole(segments=71, er=3.2, sweep=True))
        bare_resonance, bare_resistance = resonance(bare_sweep)
        coated_resonance, coated_resistance = resonance(coated_sweep)
        coarse_resonance, coarse_resistance = resonance(coarse_sweep)
        medium_resonance, medium_resistance = resonance(medium_sweep)
        if abs(coarse_resonance-coated_resonance)/coated_resonance > 0.0025:
            raise AssertionError("47/95 segment resonance convergence exceeds 0.25%")
        if abs(coarse_resistance-coated_resistance) > 0.5:
            raise AssertionError("47/95 segment resistance convergence exceeds 0.5 ohm")
        if abs(bare_resonance-295.815987)/295.815987 > 0.01:
            raise AssertionError("bare resonance misses the OpenEMS reference")
        if abs(coated_resonance-275.599687)/275.599687 > 0.01:
            raise AssertionError("coated resonance misses the OpenEMS reference")
        if abs(bare_resistance-72.797845) > 2.0:
            raise AssertionError("bare resistance misses the OpenEMS reference")
        if abs(coated_resistance-64.083814) > 2.0:
            raise AssertionError("coated resistance misses the OpenEMS reference")

        reference_cases = json.loads((ROOT / "tests/coated-wire/"
            "m1-openems-caps.json").read_text())["cases"]
        reference = next(case for case in reference_cases
                         if case["coating_free_end"] == "flush" and
                         case["coating_radial_cells"] == 8)
        reference_samples = [sample for sample in reference["far_field"]["samples"]
                             if sample["phi_deg"] == 0]
        reference_peak = max(sample["relative_db"] for sample in reference_samples)
        reference_pattern = {sample["theta_deg"]:
                             sample["relative_db"]-reference_peak
                             for sample in reference_samples}
        nec_pattern = normalized_pattern(coated_1)
        pattern_error = max(abs(nec_pattern[theta]-value)
                            for theta, value in reference_pattern.items()
                            if value >= -40 and theta in nec_pattern)
        directivity = sampled_directivity(coated_1)
        if pattern_error > 0.15:
            raise AssertionError("coated pattern misses OpenEMS by over 0.15 dB")
        if abs(directivity-reference["far_field"]["directivity_db"]) > 0.1:
            raise AssertionError("coated directivity misses OpenEMS by over 0.1 dB")

        expect_error(root, "loss", dipole(er=3.2).replace(
            "3.2000000000000002 0.0040000000000000001 0 0",
            "3.2 0.004 0.01 0"), "COATED_WIRE_LOSS_UNAVAILABLE")
        bent = """CM M3 BEND REJECTION
CE
GW 1 5 0 0 -0.25 0 0 0 0.001
GW 1 5 0 0 0 0.25 0 0 0.001
GE 0
LD 7 1 0 0 2.5 0.002 0 0
FR 0 1 0 0 300
EX 0 1 3 0 1 0
XQ
EN
"""
        expect_error(root, "bend", bent, "COATED_WIRE_M3_GEOMETRY_UNSUPPORTED")
        grounded = """CM M3 GROUND REJECTION
CE
GW 1 11 0 0 0.1 0.5 0 0.1 0.001
GE 0
LD 7 1 0 0 2.5 0.002 0 0
GN 1
FR 0 1 0 0 300
EX 0 1 6 0 1 0
XQ
EN
"""
        expect_error(root, "ground", grounded,
                     "COATED_WIRE_M3_GROUND_UNSUPPORTED")

        shift = 100*(coated_resonance/bare_resonance-1)
        print("M3 null limits, backend/process identity, and reciprocity passed")
        print(f"M3 convergence: 47={coarse_resonance:.6f} MHz/{coarse_resistance:.6f} ohm, "
              f"71={medium_resonance:.6f} MHz/{medium_resistance:.6f} ohm, "
              f"95={coated_resonance:.6f} MHz/{coated_resistance:.6f} ohm")
        print(f"M3 resonance: bare {bare_resonance:.6f} MHz/{bare_resistance:.6f} ohm, "
              f"coated {coated_resonance:.6f} MHz/{coated_resistance:.6f} ohm, "
              f"shift {shift:.3f}%")
        print(f"M3 pattern: {pattern_error:.6f} dB maximum error above -40 dB, "
              f"{directivity:.6f} dBi sampled directivity")
        print("M3 loss, bend, and ground gates passed")


if __name__ == "__main__":
    main()
