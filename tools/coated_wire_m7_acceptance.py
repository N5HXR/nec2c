#!/usr/bin/env python3
"""Assemble M7 NEC/OpenEMS Yagi evidence and comparison plots."""

import argparse
import importlib.util
import json
import math
from pathlib import Path
import tempfile


def parser():
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--bare-openems", type=Path, required=True)
    value.add_argument("--coated-openems", type=Path, required=True)
    value.add_argument("--coarse-coated-openems", type=Path)
    value.add_argument("--coarse-bare-openems", type=Path)
    value.add_argument("--refined-bare-openems", type=Path)
    value.add_argument("--output-dir", type=Path, required=True)
    return value.parse_args()


def load_test(root):
    path = root/"tests/coated-wire/m7_geometry.py"
    spec = importlib.util.spec_from_file_location("m7_geometry", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def interpolate_resistance(rows, frequency):
    for left, right in zip(rows, rows[1:]):
        if left["frequency_mhz"] <= frequency <= right["frequency_mhz"]:
            part = ((frequency-left["frequency_mhz"])/
                    (right["frequency_mhz"]-left["frequency_mhz"]))
            return (left["impedance_ohm"]["re"] + part *
                    (right["impedance_ohm"]["re"]-
                     left["impedance_ohm"]["re"]))
    raise ValueError("frequency outside impedance sweep")


def openems_case(path):
    impedance = json.loads((path/"impedance.json").read_text())
    far = json.loads((path/"far-field.json").read_text())
    current = json.loads((path/"wire-current.json").read_text())
    spec = json.loads((path/"spec.json").read_text())
    frequency = far["frequency_mhz"]
    return {
        "frequency_mhz": frequency,
        "resistance_ohm": interpolate_resistance(impedance, frequency),
        "directivity_dbi": far["directivity_dbi"],
        "pattern": [row for row in far["samples"] if row["phi_deg"] == 0],
        "current": current["samples"],
        "mesh": spec["mesh"],
        "mesh_pitch_mm": spec["mesh_pitch_mm"],
        "wall_seconds": spec.get("wall_seconds"),
        "end_criteria": spec["execution"]["end_criteria"],
    }


def nec_case(module, directory, coated):
    deck = module.yagi_deck(2)
    if not coated:
        deck = "\n".join(line for line in deck.splitlines()
                         if not line.startswith("LD 7"))
    sweep = deck.replace("FR 0 1 0 0 300 0", "FR 0 241 0 0 240 0.5")
    sweep = sweep.replace("RP 0 37 1 1000 0 0 5 0 0 0", "XQ")
    result = module.run(directory, "nec-coated" if coated else "nec-bare",
                        sweep, experimental=True)
    frequency, resistance = module.resonance(result)
    pattern_deck = deck.replace("FR 0 1 0 0 300 0",
                                f"FR 0 1 0 0 {frequency:.17g} 0")
    pattern_deck = pattern_deck.replace(
        "RP 0 37 1 1000 0 0 5 0 0 0",
        "RP 0 37 73 1000 0 0 5 5 0 0")
    pattern_result = module.run(
        directory, "nec-coated-pattern" if coated else "nec-bare-pattern",
        pattern_deck, experimental=True)
    execution = module.rows(pattern_result)[0]
    samples = execution["far_fields"][0]["samples"]
    magnitudes = [math.hypot(
        math.hypot(row["r_e_theta"]["re"], row["r_e_theta"]["im"]),
        math.hypot(row["r_e_phi"]["re"], row["r_e_phi"]["im"]))
                  for row in samples]
    peak = max(magnitudes)
    pattern = [{"theta_deg": row["theta_deg"], "phi_deg": row["phi_deg"],
                "relative_db": (-math.inf if magnitude == 0 else
                                20*math.log10(magnitude/peak))}
               for row, magnitude in zip(samples, magnitudes)
               if row["phi_deg"] == 0]
    by_phi = {}
    for row, magnitude in zip(samples, magnitudes):
        by_phi.setdefault(row["phi_deg"], []).append(
            (row["theta_deg"], magnitude*magnitude))
    theta_integrals = []
    for phi, values in sorted(by_phi.items()):
        values.sort()
        integral = sum(
            math.radians(right[0]-left[0])*0.5*(
                left[1]*math.sin(math.radians(left[0])) +
                right[1]*math.sin(math.radians(right[0])))
            for left, right in zip(values, values[1:]))
        theta_integrals.append((phi, integral))
    sphere_integral = sum(
        math.radians(right[0]-left[0])*0.5*(left[1]+right[1])
        for left, right in zip(theta_integrals, theta_integrals[1:]))
    directivity = 10*math.log10(4*math.pi*peak*peak/sphere_integral)
    wanted = {(-120, -120), (-120, 0), (-120, 120),
              (0, -120), (0, -60), (0, 60), (0, 120),
              (110, -120), (110, 0), (110, 120)}
    currents = []
    for row in execution["wire_currents"]:
        x = round(1000*row["center_m"][0])
        z = round(1000*row["center_m"][2])
        if (x, z) in wanted:
            item = row["current_a"]
            currents.append({"x_mm": x, "z_mm": z,
                             "magnitude_a": math.hypot(item["re"], item["im"])})
    peak_current = max(row["magnitude_a"] for row in currents)
    for row in currents:
        row["relative_magnitude"] = row["magnitude_a"]/peak_current
    return {"frequency_mhz": frequency, "resistance_ohm": resistance,
            "directivity_dbi": directivity,
            "pattern": pattern, "current": currents}


def pattern_error(nec, reference):
    by_angle = {row["theta_deg"]: row["relative_db"] for row in nec}
    differences = [abs(by_angle[row["theta_deg"]]-row["relative_db"])
                   for row in reference
                   if row["relative_db"] >= -40 and row["theta_deg"] in by_angle]
    return max(differences)


def current_error(nec, reference):
    by_position = {(row["x_mm"], row["z_mm"]): row["relative_magnitude"]
                   for row in nec}
    return max(abs(by_position[(row["x_mm"], row["z_mm"])]-
                   row["relative_magnitude"]) for row in reference)


def convergence(coarse, fine):
    return {
        "resonance_percent_coarse_to_fine": 100*(
            coarse["frequency_mhz"]/fine["frequency_mhz"]-1),
        "resistance_ohm_coarse_to_fine": (
            coarse["resistance_ohm"]-fine["resistance_ohm"]),
        "pattern_maximum_db_above_minus_40": pattern_error(
            coarse["pattern"], fine["pattern"]),
        "current_maximum_normalized_delta": current_error(
            coarse["current"], fine["current"]),
    }


def plot(nec, reference, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    current_keys = [(row["x_mm"], row["z_mm"])
                    for row in reference["current"]]
    element_name = {-120: "Ref", 0: "Drv", 110: "Dir"}
    for label, case in (("NEC coated", nec), ("OpenEMS coated", reference)):
        axes[0].plot([row["theta_deg"] for row in case["pattern"]],
                     [max(-40, row["relative_db"]) for row in case["pattern"]],
                     label=label)
        by_key = {(row["x_mm"], row["z_mm"]): row["relative_magnitude"]
                  for row in case["current"]}
        axes[1].plot(range(len(current_keys)),
                     [by_key[key] for key in current_keys], "o-", label=label)
    axes[0].set(xlabel="Theta at phi=0 (degrees)", ylabel="Relative field (dB)",
                ylim=(-40, 1))
    axes[1].set(xlabel="Element and axial position (mm)",
                ylabel="Normalized current magnitude", ylim=(0, 1.05),
                xticks=range(len(current_keys)),
                xticklabels=[f"{element_name[x]}\n{z:+g}" for x, z in current_keys])
    for axis in axes:
        axis.grid(True, alpha=.25)
        axis.legend()
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def main():
    a = parser()
    root = Path(__file__).resolve().parents[1]
    module = load_test(root)
    with tempfile.TemporaryDirectory(prefix="nec2c-m7-accept-") as name:
        directory = Path(name)
        nec_bare = nec_case(module, directory, False)
        nec_coated = nec_case(module, directory, True)
    openems_bare = openems_case(a.bare_openems)
    openems_coated = openems_case(a.coated_openems)
    comparison = {
        "bare": {
            "resonance_error_percent": 100*(nec_bare["frequency_mhz"]/
                                               openems_bare["frequency_mhz"]-1),
            "resistance_error_ohm": (nec_bare["resistance_ohm"]-
                                       openems_bare["resistance_ohm"]),
        },
        "coated": {
            "resonance_error_percent": 100*(nec_coated["frequency_mhz"]/
                                               openems_coated["frequency_mhz"]-1),
            "resistance_error_ohm": (nec_coated["resistance_ohm"]-
                                       openems_coated["resistance_ohm"]),
            "directivity_error_db": (nec_coated["directivity_dbi"]-
                                       openems_coated["directivity_dbi"]),
            "pattern_maximum_db_above_minus_40": pattern_error(
                nec_coated["pattern"], openems_coated["pattern"]),
            "current_maximum_normalized_delta": current_error(
                nec_coated["current"], openems_coated["current"]),
        },
        "coating_effect": {
            "nec_resonance_shift_percent": 100*(
                nec_coated["frequency_mhz"]/nec_bare["frequency_mhz"]-1),
            "openems_resonance_shift_percent": 100*(
                openems_coated["frequency_mhz"]/
                openems_bare["frequency_mhz"]-1),
            "shift_error_percentage_points": 100*(
                nec_coated["frequency_mhz"]/nec_bare["frequency_mhz"]-
                openems_coated["frequency_mhz"]/
                openems_bare["frequency_mhz"]),
            "resistance_shift_error_ohm": (
                (nec_coated["resistance_ohm"]-nec_bare["resistance_ohm"])-
                (openems_coated["resistance_ohm"]-
                 openems_bare["resistance_ohm"])),
        },
    }
    record = {"schema": "nec2c-coated-wire-m7-acceptance",
              "nec": {"bare": nec_bare, "coated": nec_coated},
              "openems": {"bare": openems_bare, "coated": openems_coated},
              "comparison": comparison}
    if a.coarse_coated_openems:
        coarse = openems_case(a.coarse_coated_openems)
        record["openems"]["coated_coarse"] = coarse
        record["openems_convergence"] = convergence(coarse, openems_coated)
    if a.coarse_bare_openems:
        coarse_bare = openems_case(a.coarse_bare_openems)
        record["openems"]["bare_coarse"] = coarse_bare
        record["bare_convergence_2mm_to_1mm"] = convergence(
            coarse_bare, openems_bare)
    if a.refined_bare_openems:
        refined_bare = openems_case(a.refined_bare_openems)
        record["openems"]["bare_refined_transverse"] = refined_bare
        record["bare_convergence_1mm_to_half_mm_transverse"] = convergence(
            openems_bare, refined_bare)
        failed = abs(record["bare_convergence_1mm_to_half_mm_transverse"][
            "resonance_percent_coarse_to_fine"]) > 1.0
        record["decision"] = {
            "status": "experimental-only" if failed else "promotion-candidate",
            "reason": ("bare Cartesian resonance failed the 1 percent "
                       "reference-convergence gate" if failed else
                       "bare Cartesian reference met the convergence gate"),
        }
    a.output_dir.mkdir(parents=True, exist_ok=True)
    (a.output_dir/"m7-acceptance.json").write_text(
        json.dumps(record, indent=2, sort_keys=True)+"\n")
    plot(nec_coated, openems_coated, a.output_dir/"m7-yagi.png")
    print(json.dumps(comparison, indent=2))
    if "openems_convergence" in record:
        print("OpenEMS coarse/fine convergence:")
        print(json.dumps(record["openems_convergence"], indent=2))


if __name__ == "__main__":
    main()
