#!/usr/bin/env python3
"""M10A matrix normalization, passivity, and first convergence gates."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
HERE = ROOT / "tests" / "nec4-similarity"


def pair(value):
    return complex(value["re"], value["im"])


def run(source, directory, stem, model="grid", ground="nec4", threads="auto"):
    result = directory / f"{stem}.json"
    report = directory / f"{stem}.out"
    command = [str(BINARY), "--wire-engine", "nec4", "--ground-engine", ground,
               "--nec4-ground-model", model, "--threads", str(threads),
               "--result-json", str(result), "-i", str(source), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                               timeout=300)
    if completed.returncode:
        raise AssertionError((command, completed.returncode, completed.stderr,
                              report.read_text(errors="replace")[-4000:]
                              if report.exists() else ""))
    return json.loads(result.read_text())["jobs"][0]["executions"][0]


def deck(upper_rest, lower_rest):
    return ("CM M10A SEAWATER SEGMENTATION WITNESS\nCE\n"
            "GW 1 1 0 0 .1 0 0 0 .003\n"
            f"GW 3 {upper_rest} 0 0 5 0 0 .1 .003\n"
            "GW 2 1 0 0 0 0 0 -.1 .003\n"
            f"GW 4 {lower_rest} 0 0 -.1 0 0 -1 .003\n"
            "GE 0\nGN 2 0 0 0 80 4\nFR 0 1 0 0 14 0\n"
            "EX 0 1 1 0 1 0\nXQ\nEN\n")


def crossing_deck(ground):
    return ("CM M10A IDENTICAL-MEDIA CROSSING NULL\nCE\n"
            "GW 1 21 0 0 -4.75 0 0 0 .003\n"
            "GW 2 21 0 0 0 0 0 4.75 .003\n"
            "GE 0\n" + ground + "\nFR 0 1 0 0 15 0\n"
            "EX 0 2 1 0 1 0\nXQ\nEN\n")


def main():
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m10a-") as name:
        directory = Path(name)
        source = HERE / "m9-seawater-whip.nec"
        grid = run(source, directory, "grid")
        direct = run(source, directory, "direct", "direct")
        zg = pair(grid["feeds"][0]["impedance_ohm"])
        zd = pair(direct["feeds"][0]["impedance_ohm"])
        z_delta = abs(zg-zd) / max(abs(zd), 1e-300)
        ig = [pair(row["current_a"]) for row in grid["wire_currents"]]
        idirect = [pair(row["current_a"]) for row in direct["wire_currents"]]
        i_delta = math.sqrt(sum(abs(a-b)**2 for a, b in zip(ig, idirect)) /
                            sum(abs(value)**2 for value in idirect))
        if zg.real <= 0 or zd.real <= 0:
            raise AssertionError(("non-passive feed resistance", zg, zd))
        if z_delta > 5e-3 or i_delta > 5e-3:
            raise AssertionError(("grid/direct conditioned solve", z_delta, i_delta))

        # A ground plane separating identical media is physically invisible.
        # Compare the complete buried-wire matrix path with the ordinary
        # free-space NEC4 wire path, including feed impedance and every solved
        # segment current.  This also catches mixed historical/SI constants.
        free = run(HERE / "m10a-identical-free.nec", directory,
                   "identical-free", ground="nec2", threads=1)
        buried = run(HERE / "m10a-identical-buried.nec", directory,
                     "identical-buried", threads=1)
        z_free = pair(free["feeds"][0]["impedance_ohm"])
        z_buried = pair(buried["feeds"][0]["impedance_ohm"])
        identical_z_delta = abs(z_buried-z_free) / max(abs(z_free), 1e-300)
        i_free = [pair(row["current_a"]) for row in free["wire_currents"]]
        i_buried = [pair(row["current_a"]) for row in buried["wire_currents"]]
        identical_i_delta = math.sqrt(
            sum(abs(a-b)**2 for a, b in zip(i_free, i_buried)) /
            sum(abs(value)**2 for value in i_free))
        if identical_z_delta > 1e-9 or identical_i_delta > 1e-9:
            raise AssertionError(("identical-media null", identical_z_delta,
                                  identical_i_delta))

        crossing_free_path = directory / "crossing-free.nec"
        crossing_ground_path = directory / "crossing-ground.nec"
        crossing_free_path.write_text(crossing_deck("GN -1"), encoding="ascii")
        crossing_ground_path.write_text(
            crossing_deck("GN 2 0 0 0 1 0"), encoding="ascii")
        crossing_free = run(crossing_free_path, directory, "crossing-free",
                            ground="nec2", threads=1)
        crossing_ground = run(crossing_ground_path, directory,
                              "crossing-ground", threads=1)
        z_crossing_free = pair(crossing_free["feeds"][0]["impedance_ohm"])
        z_crossing_ground = pair(crossing_ground["feeds"][0]["impedance_ohm"])
        crossing_z_delta = abs(z_crossing_ground-z_crossing_free) / abs(z_crossing_free)
        i_crossing_free = [pair(row["current_a"])
                           for row in crossing_free["wire_currents"]]
        i_crossing_ground = [pair(row["current_a"])
                             for row in crossing_ground["wire_currents"]]
        crossing_i_delta = math.sqrt(
            sum(abs(a-b)**2 for a, b in zip(i_crossing_free,
                                             i_crossing_ground)) /
            sum(abs(value)**2 for value in i_crossing_free))
        if crossing_z_delta > 1e-9 or crossing_i_delta > 1e-9:
            raise AssertionError(("identical-media crossing null",
                                  crossing_z_delta, crossing_i_delta))

        # Thread count must not alter the ordered matrix solution.  The JSON
        # execution metadata may differ, so compare the numerical payload.
        serial = run(source, directory, "grid-serial", threads=1)
        z_serial = pair(serial["feeds"][0]["impedance_ohm"])
        i_serial = [pair(row["current_a"]) for row in serial["wire_currents"]]
        parallel_z_delta = abs(zg-z_serial) / max(abs(z_serial), 1e-300)
        parallel_i_delta = math.sqrt(
            sum(abs(a-b)**2 for a, b in zip(ig, i_serial)) /
            sum(abs(value)**2 for value in i_serial))
        if parallel_z_delta > 2e-13 or parallel_i_delta > 2e-13:
            raise AssertionError(("thread identity", parallel_z_delta,
                                  parallel_i_delta))

        sweep = []
        for upper, lower in ((49, 9), (98, 18), (147, 27)):
            path = directory / f"seawater-{upper}-{lower}.nec"
            path.write_text(deck(upper, lower), encoding="ascii")
            execution = run(path, directory, f"seawater-{upper}-{lower}")
            impedance = pair(execution["feeds"][0]["impedance_ohm"])
            if impedance.real <= 0:
                raise AssertionError(("segmentation passivity", upper, lower,
                                      impedance))
            sweep.append({"upper_rest_segments": upper,
                          "lower_rest_segments": lower,
                          "fixed_feed_and_interface_segments": 2,
                          "impedance_ohm": {"re": impedance.real,
                                             "im": impedance.imag}})

        z_penultimate = complex(
            sweep[-2]["impedance_ohm"]["re"],
            sweep[-2]["impedance_ohm"]["im"])
        z_finest = complex(sweep[-1]["impedance_ohm"]["re"],
                           sweep[-1]["impedance_ohm"]["im"])
        convergence_delta = abs(z_finest-z_penultimate) / abs(z_finest)
        convergence_re_delta = abs(z_finest.real-z_penultimate.real) / abs(z_finest.real)
        convergence_im_delta = abs(z_finest.imag-z_penultimate.imag) / abs(z_finest.imag)
        if convergence_re_delta > 1e-2 or convergence_im_delta > 1e-2:
            raise AssertionError(("segmentation convergence", convergence_re_delta,
                                  convergence_im_delta))

        evidence = {
            "schema": "nec2c-nec4-m10a-normalization",
            "schema_version": 1,
            "seawater_5m_whip_1m_stake": {
                "grid_impedance_ohm": {"re": zg.real, "im": zg.imag},
                "direct_impedance_ohm": {"re": zd.real, "im": zd.imag},
                "impedance_relative_delta": z_delta,
                "current_l2_relative_delta": i_delta,
                "passive": True,
                "thread_identity": {
                    "impedance_relative_delta": parallel_z_delta,
                    "current_l2_relative_delta": parallel_i_delta,
                },
                "segmentation_sweep": sweep,
                "two_finest_impedance_relative_delta": convergence_delta,
                "two_finest_resistance_relative_delta": convergence_re_delta,
                "two_finest_reactance_relative_delta": convergence_im_delta,
            },
            "identical_media_null": {
                "free_space_impedance_ohm": {"re": z_free.real,
                                               "im": z_free.imag},
                "buried_impedance_ohm": {"re": z_buried.real,
                                           "im": z_buried.imag},
                "impedance_relative_delta": identical_z_delta,
                "current_l2_relative_delta": identical_i_delta,
                "crossing_impedance_relative_delta": crossing_z_delta,
                "crossing_current_l2_relative_delta": crossing_i_delta,
            },
            "unit_gates": {
                "si_to_nec_maxwell_similarity": "tests/nec4_ground_matrix_test.c",
                "complex_medium_kernel_similarity": "tests/nec4_kernel_test.c",
            },
        }
        (HERE / "m10a-normalization-gate.json").write_text(
            json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
        print(f"M10A seawater grid/direct: Z {z_delta:.3e}, currents {i_delta:.3e}")
        print(f"M10A identical-media null: Z {identical_z_delta:.3e}, "
              f"currents {identical_i_delta:.3e}")
        print(f"M10A identical-media crossing: Z {crossing_z_delta:.3e}, "
              f"currents {crossing_i_delta:.3e}")
        print(f"M10A thread identity: Z {parallel_z_delta:.3e}, "
              f"currents {parallel_i_delta:.3e}")
        print("M10A seawater segmentation: " + ", ".join(
            f"{row['upper_rest_segments']}+{row['lower_rest_segments']}+2="
            f"{row['impedance_ohm']['re']:.6f}"
            f"{row['impedance_ohm']['im']:+.6f}j ohm" for row in sweep))


if __name__ == "__main__":
    main()
