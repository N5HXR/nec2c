#!/usr/bin/env python3
"""M7 general-geometry and mutual-coupling acceptance checks."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()


def finish(geometry, tags, source_tag, source_segment, frequency=285.0,
           sweep=False, coating=True, pattern=False, loss_tangent=0.0):
    lines = ["CM M7 GENERAL COATED-WIRE GEOMETRY", "CE", *geometry, "GE 0"]
    if coating:
        lines.extend(f"LD 7 {tag} 0 0 3.2 0.004 {loss_tangent:.17g} 0"
                     for tag in tags)
    lines.append("FR 0 121 0 0 240 1" if sweep else
                 f"FR 0 1 0 0 {frequency:.17g} 0")
    lines.append(f"EX 0 {source_tag} {source_segment} 0 1 0")
    lines.append("RP 0 37 1 1000 0 0 5 0 0 0" if pattern else "XQ")
    lines.extend(("EN", ""))
    return "\n".join(lines)


def v_deck(segments=25, source_tag=2, source_segment=None, **options):
    center_segments = round(8*segments/15)
    center_segments += 1-center_segments % 2
    if source_segment is None:
        source_segment = ((center_segments+1)//2 if source_tag == 2 else
                          (segments+1)//2)
    geometry = [
        f"GW 1 {segments} -0.05 0 -0.187 -0.05 0 0 0.002",
        f"GW 2 {center_segments} -0.05 0 0 0.05 0 0 0.002",
        f"GW 3 {segments} 0.05 0 0 0.05 0 0.187 0.002",
    ]
    return finish(geometry, (1, 2, 3), source_tag, source_segment, **options)


def chamfer_deck(offset=0.01131370849898476, **options):
    """Planar Z with four 45-degree vertices around two chamfers."""
    points = [(-0.05, 0, -0.19637258300203048),
              (-0.05, 0, -offset), (-0.05+offset, 0, 0),
              (0.05-offset, 0, 0), (0.05, 0, offset),
              (0.05, 0, 0.19637258300203048)]
    geometry, counts = [], []
    for tag, (left, right) in enumerate(zip(points, points[1:]), 1):
        count = max(1, round(math.dist(left, right)/0.005))
        if tag == 3 and count % 2 == 0:
            count += 1
        counts.append(count)
        geometry.append(
            "GW %d %d %.17g %.17g %.17g %.17g %.17g %.17g 0.002" %
            (tag, count, *left, *right))
    return finish(geometry, tuple(range(1, 6)), 3,
                  (counts[2]+1)//2, **options)


def hairpin_deck():
    geometry = [
        "GW 1 41 0 0 -0.2 0 0 0 0.002",
        "GW 2 41 0 0 0 0.1414213562373095 0 -0.1414213562373095 0.002",
    ]
    return finish(geometry, (1, 2), 1, 21)


def parallel_deck(source_tag=1):
    geometry = [
        "GW 1 41 -0.08 0 -0.23 -0.08 0 0.23 0.002",
        "GW 2 41 0.08 0 -0.23 0.08 0 0.23 0.002",
    ]
    return finish(geometry, (1, 2), source_tag, 21, frequency=290.0)


def yagi_deck(source_tag=2):
    geometry = [
        # Equal 10 mm segment lengths keep the reciprocity comparison from
        # measuring the classic NEC unequal-discretization defect.
        "GW 1 49 -0.12 0 -0.245 -0.12 0 0.245 0.002",
        "GW 2 45 0 0 -0.225 0 0 0.225 0.002",
        "GW 3 43 0.11 0 -0.215 0.11 0 0.215 0.002",
    ]
    return finish(geometry, (1, 2, 3), source_tag,
                  {1: 25, 2: 23, 3: 22}[source_tag],
                  frequency=300.0, pattern=True)


def loop_deck(source_tag=1, source_segment=5):
    geometry = [
        "GW 1 10 -0.12 -0.12 0 0.12 -0.12 0 0.002",
        "GW 2 10 0.12 -0.12 0 0.12 0.12 0 0.002",
        "GW 3 10 0.12 0.12 0 -0.12 0.12 0 0.002",
        "GW 4 10 -0.12 0.12 0 -0.12 -0.12 0 0.002",
    ]
    return finish(geometry, (1, 2, 3, 4), source_tag, source_segment,
                  frequency=300.0)


def spatial_deck(source_tag=1, source_segment=7):
    geometry = [
        # Three equal-length, non-coplanar legs.
        "GW 1 13 -0.20 0 0 0 0 0 0.002",
        "GW 2 13 0 0 0 0.16 0.12 0 0.002",
        "GW 3 13 0.16 0.12 0 0.16 0.12 0.20 0.002",
    ]
    return finish(geometry, (1, 2, 3), source_tag, source_segment,
                  frequency=310.0)


def run(directory, name, deck, backend="builtin", threads=1,
        experimental=False):
    source = directory / f"{name}.nec"
    report = directory / f"{name}.out"
    result = directory / f"{name}.json"
    source.write_text(deck, encoding="ascii")
    command = [
        str(BINARY), "--backend", backend, "--jobs", "1", "--threads",
        str(threads), "--result-json", str(result)]
    if experimental:
        command.extend(("--coated-wire-junctions", "experimental"))
    command.extend(("-i", str(source), "-o", str(report)))
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    parsed = json.loads(result.read_text()) if result.exists() else None
    if completed.returncode:
        diagnostic = (parsed.get("diagnostics", [{}])[0].get("code")
                      if parsed else None)
        raise RuntimeError(f"{name}: {diagnostic}\n{completed.stdout}\n"
                           f"{completed.stderr}")
    return parsed


def value(item):
    return complex(item["re"], item["im"])


def rows(result):
    return result["jobs"][0]["executions"]


def current(result, segment):
    return value(rows(result)[0]["wire_currents"][segment-1]["current_a"])


def reciprocity(directory, name, left_deck, right_deck, left_observation,
                right_observation, experimental=False):
    left = run(directory, f"{name}-left", left_deck,
               experimental=experimental)
    right = run(directory, f"{name}-right", right_deck,
                experimental=experimental)
    transfer_left = current(left, left_observation)
    transfer_right = current(right, right_observation)
    return abs(transfer_left-transfer_right)/max(
        abs(transfer_left), abs(transfer_right), 1.0e-30)


def resonance(result):
    for left, right in zip(rows(result), rows(result)[1:]):
        zl = value(left["feeds"][0]["impedance_ohm"])
        zr = value(right["feeds"][0]["impedance_ohm"])
        if zl.imag == 0 or zl.imag*zr.imag < 0:
            part = 0 if zl.imag == 0 else -zl.imag/(zr.imag-zl.imag)
            return (left["frequency_mhz"] +
                    part*(right["frequency_mhz"]-left["frequency_mhz"]),
                    zl.real+part*(zr.real-zl.real))
    raise AssertionError("resonance was not bracketed")


def normalized_execution(result):
    def clean(item):
        if isinstance(item, dict):
            return {key: clean(val) for key, val in item.items()
                    if key not in {"elapsed_seconds", "material_evaluations"}}
        if isinstance(item, list):
            return [clean(val) for val in item]
        return item
    return clean(rows(result))


def main():
    acceptance = json.loads(
        (ROOT / "tests/coated-wire/m7-acceptance.json").read_text(
            encoding="utf-8"))
    comparison = acceptance["comparison"]
    if acceptance["decision"]["status"] != "experimental-only":
        raise AssertionError(acceptance["decision"])
    if abs(acceptance["bare_convergence_1mm_to_half_mm_transverse"][
            "resonance_percent_coarse_to_fine"]) <= 1.0:
        raise AssertionError("retained negative reference unexpectedly converged")
    for name in ("bare", "coated"):
        if abs(comparison[name]["resistance_error_ohm"]) > 2.0:
            raise AssertionError(comparison[name])
    effect = comparison["coating_effect"]
    if abs(effect["shift_error_percentage_points"]) > 1.0:
        raise AssertionError(effect)
    if abs(effect["resistance_shift_error_ohm"]) > 2.0:
        raise AssertionError(effect)
    coated_reference = comparison["coated"]
    if abs(coated_reference["directivity_error_db"]) > 0.2:
        raise AssertionError(coated_reference)
    if coated_reference["pattern_maximum_db_above_minus_40"] > 0.5:
        raise AssertionError(coated_reference)
    if coated_reference["current_maximum_normalized_delta"] > 0.05:
        raise AssertionError(coated_reference)

    with tempfile.TemporaryDirectory(prefix="nec2c-coated-m7-") as temporary:
        root = Path(temporary)
        qualified = run(root, "qualified-bend", chamfer_deck())
        coated = qualified["jobs"][0]["coated_wires"]
        if coated["operator_status"] != "lossless-qualified-bend":
            raise AssertionError("qualified bend status was not stamped")
        policy = coated.get("bend_policy", {})
        if policy.get("qualified_junction_count") != 4:
            raise AssertionError("qualified bend junction count changed")
        if policy.get("observed_maximum_turn_deg", 100) > 45.000000001:
            raise AssertionError("qualified bend turn limit changed")
        if policy.get("observed_minimum_chord_outer_radius_ratio", 0) < 4:
            raise AssertionError("qualified bend chord limit changed")
        warning_codes = {row["code"] for row in
                         qualified.get("diagnostics", [])}
        if "COATED_WIRE_EXPERIMENTAL_GEOMETRY" in warning_codes:
            raise AssertionError("qualified bend was marked experimental")
        qualified_lossy = run(root, "qualified-bend-lossy",
                              chamfer_deck(loss_tangent=0.02))
        if (qualified_lossy["jobs"][0]["coated_wires"]["operator_status"] !=
                "lossy-qualified-bend"):
            raise AssertionError("lossy qualified bend status was not stamped")
        qualified_power = rows(qualified_lossy)[0]["power"]
        if qualified_power["dielectric_absorbed_w"] <= 0 or \
                qualified_power["balance_defect_relative"] > 0.003:
            raise AssertionError("lossy qualified bend power gate failed")
        try:
            run(root, "short-chamfer-rejection", chamfer_deck(0.008))
            raise AssertionError("strict mode accepted a short bend chord")
        except RuntimeError as error:
            if "COATED_WIRE_BEND_UNVALIDATED" not in str(error):
                raise
        try:
            run(root, "hairpin-rejection", hairpin_deck())
            raise AssertionError("strict mode accepted a 135-degree reversal")
        except RuntimeError as error:
            if "COATED_WIRE_BEND_UNVALIDATED" not in str(error):
                raise
        try:
            run(root, "v-strict-rejection", v_deck(pattern=True))
            raise AssertionError("strict mode accepted an unvalidated bend")
        except RuntimeError as error:
            if "COATED_WIRE_BEND_UNVALIDATED" not in str(error):
                raise
        try:
            run(root, "yagi-strict-rejection", yagi_deck(2))
            raise AssertionError("strict mode accepted unvalidated disconnected geometry")
        except RuntimeError as error:
            if "COATED_WIRE_GENERAL_GEOMETRY_UNVALIDATED" not in str(error):
                raise
        experimental_yagi = run(root, "yagi-experimental-status",
                                yagi_deck(2), experimental=True)
        if (experimental_yagi["jobs"][0]["coated_wires"]["operator_status"] !=
                "experimental-geometry"):
            raise AssertionError("experimental Yagi status was not stamped")
        warning_codes = {row["code"] for row in
                         experimental_yagi.get("diagnostics", [])}
        if "COATED_WIRE_EXPERIMENTAL_GEOMETRY" not in warning_codes:
            raise AssertionError("experimental Yagi warning is missing")

        bare_v = run(root, "v-bare", v_deck(coating=False, pattern=True),
                     experimental=True)
        er_one_v = run(root, "v-er-one", v_deck(pattern=True).replace(
            "3.2 0.004", "1 0.004"), experimental=True)
        if normalized_execution(bare_v) != normalized_execution(er_one_v):
            raise AssertionError("bent-wire Er=1 null limit changed")

        coated_v = run(root, "v-coated", v_deck(pattern=True),
                       experimental=True)
        coated_v_threads = run(root, "v-coated-4t", v_deck(pattern=True),
                                threads=4, experimental=True)
        coated_v_lapack = run(root, "v-coated-lapack", v_deck(pattern=True),
                               backend="lapack", experimental=True)
        if normalized_execution(coated_v) != normalized_execution(coated_v_threads):
            raise AssertionError("bent-wire threaded result changed")

        # Backend comparison is tolerant only at the expected factorization
        # roundoff scale.
        z_builtin = value(rows(coated_v)[0]["feeds"][0]["impedance_ohm"])
        z_lapack = value(rows(coated_v_lapack)[0]["feeds"][0]["impedance_ohm"])
        if abs(z_builtin-z_lapack)/max(abs(z_builtin), 1.0) > 1.0e-9:
            raise AssertionError("bent-wire backend result changed")

        v_defect = reciprocity(
            root, "v", v_deck(source_tag=1, source_segment=13),
            v_deck(source_tag=3, source_segment=13), 25+13+13, 13,
            experimental=True)
        parallel_defect = reciprocity(
            root, "parallel", parallel_deck(1), parallel_deck(2), 41+21, 21,
            experimental=True)
        yagi_defect = reciprocity(
            root, "yagi", yagi_deck(2), yagi_deck(3), 49+45+22, 49+23,
            experimental=True)
        loop_defect = reciprocity(
            root, "loop", loop_deck(1, 5), loop_deck(3, 5), 20+5, 5,
            experimental=True)
        spatial_defect = reciprocity(
            root, "spatial", spatial_deck(1, 7), spatial_deck(3, 7),
            2*13+7, 7, experimental=True)
        defects = (v_defect, parallel_defect, yagi_defect, loop_defect,
                   spatial_defect)

        # NEC's point-matched wire basis has a small, known reciprocity defect
        # for unequal element lengths and multi-leg junctions.  Require exact
        # reciprocity where the bare model has it, and otherwise require that
        # coating does not materially enlarge the bare discretization defect.
        def bare(deck):
            return "\n".join(line for line in deck.splitlines()
                             if not line.startswith("LD 7"))
        bare_defects = (
            reciprocity(root, "v-bare-recip", bare(v_deck(source_tag=1,
                        source_segment=13)), bare(v_deck(source_tag=3,
                        source_segment=13)), 25+13+13, 13,
                        experimental=True),
            reciprocity(root, "parallel-bare-recip", bare(parallel_deck(1)),
                        bare(parallel_deck(2)), 41+21, 21,
                        experimental=True),
            reciprocity(root, "yagi-bare-recip", bare(yagi_deck(2)),
                        bare(yagi_deck(3)), 49+45+22, 49+23,
                        experimental=True),
            reciprocity(root, "loop-bare-recip", bare(loop_deck(1, 5)),
                        bare(loop_deck(3, 5)), 20+5, 5,
                        experimental=True),
            reciprocity(root, "spatial-bare-recip", bare(spatial_deck(1, 7)),
                        bare(spatial_deck(3, 7)), 2*13+7, 7,
                        experimental=True),
        )
        for coated, uncoated in zip(defects, bare_defects):
            if coated > max(1.0e-7, 1.25*uncoated):
                raise AssertionError(
                    f"coating enlarged reciprocity defect: {coated} vs {uncoated}")

        # M4's loss and power formulation must remain closed when exercised
        # through a bend and through disconnected mutual coupling.
        lossy_results = [
            run(root, "v-lossy", v_deck(segments=51, loss_tangent=0.02),
                experimental=True),
            run(root, "parallel-lossy", parallel_deck(1).replace(
                "0.004 0 0", "0.004 0.02 0"), experimental=True),
        ]
        lossy_defects = []
        for result in lossy_results:
            power = rows(result)[0]["power"]
            if power["dielectric_absorbed_w"] <= 0:
                raise AssertionError("lossy general geometry is not passive")
            if power["balance_defect_relative"] > 0.003:
                raise AssertionError(power)
            lossy_defects.append(power["balance_defect_relative"])

        convergence = []
        for count in (15, 25, 35):
            result = run(root, f"v-sweep-{count}",
                         v_deck(segments=count, sweep=True), experimental=True)
            convergence.append((count, *resonance(result)))
        coarse, medium, fine = convergence
        if abs(coarse[1]-fine[1])/fine[1] > 0.01:
            raise AssertionError("V resonance convergence exceeds 1 percent")
        if abs(coarse[2]-fine[2]) > 2.0:
            raise AssertionError("V resistance convergence exceeds 2 ohms")

    print("M7 disconnected/Yagi geometry solved; experimental bend, loop, and 3-D cases solved")
    print("M7 X=45/Y=4 qualified bend accepted; short chord, 90-degree, and "
          "135-degree vertices rejected")
    print("M7 reciprocity defects: " + ", ".join(f"{item:.3e}" for item in defects))
    print("M7 bare reference defects: " +
          ", ".join(f"{item:.3e}" for item in bare_defects))
    print("M7 lossy power defects: " +
          ", ".join(f"{100*item:.4f}%" for item in lossy_defects))
    print("M7 experimental 90-degree-bend NEC convergence: " + ", ".join(
        f"{count}x2={frequency:.6f} MHz/{resistance:.6f} ohm"
        for count, frequency, resistance in convergence))
    print(f"M7 Yagi coating-effect errors: "
          f"{effect['shift_error_percentage_points']:+.4f} resonance point, "
          f"{effect['resistance_shift_error_ohm']:+.4f} ohm")
    print("M7 bare Cartesian reference failed convergence; general geometry remains experimental")
    print("M7 null-limit, backend, thread, pattern, and topology gates passed")


if __name__ == "__main__":
    main()
