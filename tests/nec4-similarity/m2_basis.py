#!/usr/bin/env python3
"""M2 representation equivalence and short-segment conditioning acceptance."""

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")

ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
FIXTURES = [
    ROOT / "tests" / "coated-wire" / "m0-bare-dipole.nec",
    ROOT / "tests" / "nec4-similarity" / "m2-loop.nec",
    ROOT / "tests" / "decks" / "nt-network.nec",
    ROOT / "tests" / "performance" / "m3-real-ground-sweep.nec",
]
SHORT_FREQUENCIES_MHZ = [300.0, 3.0, 0.03, 0.0003]


def run(deck, engine, directory, stem):
    result = directory / f"{stem}-{engine}.json"
    report = directory / f"{stem}-{engine}.out"
    command = [str(BINARY), "--wire-engine", engine, "--ground-engine",
               "nec2", "--backend", "builtin", "--fft-backend", "builtin",
               "--jobs", "1", "--threads", "1", "--result-json",
               str(result), "-i", str(deck), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((command, completed.returncode,
                              completed.stdout, completed.stderr))
    return json.loads(result.read_text(encoding="utf-8"))


def capabilities():
    completed = subprocess.run([str(BINARY), "--capabilities-json"],
                               cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise AssertionError((completed.args, completed.returncode,
                              completed.stdout, completed.stderr))
    return json.loads(completed.stdout)


def scaled_numeric_differences(left, right, path=""):
    values = []
    if isinstance(left, dict) and isinstance(right, dict):
        for key in left.keys() & right.keys():
            if key not in {"elapsed_seconds", "backend", "requested", "used"}:
                values.extend(scaled_numeric_differences(
                    left[key], right[key], f"{path}/{key}"))
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise AssertionError((path, len(left), len(right)))
        for index, (a, b) in enumerate(zip(left, right)):
            values.extend(scaled_numeric_differences(
                a, b, f"{path}/{index}"))
    elif (isinstance(left, (int, float)) and not isinstance(left, bool) and
          isinstance(right, (int, float)) and not isinstance(right, bool)):
        values.append((abs(left-right)/max(1.0, abs(left), abs(right)), path))
    return values


def feed_impedance(document):
    feed = document["jobs"][0]["executions"][0]["feeds"][0]
    value = feed["impedance_ohm"]
    return complex(value["re"], value["im"])


def currents(document):
    execution = document["jobs"][0]["executions"][0]
    return [abs(complex(item["current_a"]["re"], item["current_a"]["im"]))
            for item in execution["wire_currents"]]


def short_deck(path, frequency):
    path.write_text(
        "CM M2 ELECTRICALLY SHORT FIXED DIPOLE\nCE\n"
        "GW 1 21 0 0 -0.5 0 0 0.5 0.001\nGE 0\n"
        f"FR 0 1 0 0 {frequency:.12g} 0\n"
        "EX 0 1 11 0 1 0\nXQ\nEN\n", encoding="ascii")


def main():
    advertised = capabilities()
    wire_engines = advertised["features"]["wire_engines"]
    if "nec4" not in wire_engines:
        raise AssertionError(("NEC4 engine is not advertised", wire_engines))

    maximum_equivalence = 0.0
    equivalence = []
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m2-") as name:
        temporary = Path(name)
        dipole_documents = None
        loop_documents = None
        for fixture in FIXTURES:
            pair = {engine: run(fixture, engine, temporary,
                                fixture.stem.replace(".", "-"))
                    for engine in ("nec2", "nec4")}
            differences = scaled_numeric_differences(pair["nec2"], pair["nec4"])
            worst = max(differences, default=(0.0, ""))
            maximum_equivalence = max(maximum_equivalence, worst[0])
            equivalence.append({"deck": str(fixture.relative_to(ROOT)),
                                "maximum_scaled_delta": worst[0],
                                "maximum_delta_path": worst[1]})
            if fixture.name == "m0-bare-dipole.nec": dipole_documents = pair
            if fixture.name == "m2-loop.nec": loop_documents = pair

        short_results = {engine: [] for engine in ("nec2", "nec4")}
        for index, frequency in enumerate(SHORT_FREQUENCIES_MHZ):
            deck = temporary / f"short-{index}.nec"
            short_deck(deck, frequency)
            for engine in short_results:
                document = run(deck, engine, temporary, f"short-{index}")
                impedance = feed_impedance(document)
                short_results[engine].append({
                    "frequency_mhz": frequency,
                    "impedance_ohm": [impedance.real, impedance.imag],
                    "capacitive_scale": -impedance.imag*frequency,
                })

        reference = short_results["nec4"][1]["capacitive_scale"]
        nec2_drift = abs(short_results["nec2"][-1]["capacitive_scale"] /
                         reference-1.0)
        nec4_drift = abs(short_results["nec4"][-1]["capacitive_scale"] /
                         reference-1.0)
        if maximum_equivalence > 5.0e-9:
            raise AssertionError(("ordinary representation equivalence",
                                  maximum_equivalence, equivalence))
        if nec4_drift > 2.0e-3:
            raise AssertionError(("NEC4 short-segment drift", nec4_drift))
        if nec2_drift < 5.0*nec4_drift:
            raise AssertionError(("conditioning improvement", nec2_drift,
                                  nec4_drift))

        evidence = {
            "schema": "nec2c-nec4-m2-acceptance", "schema_version": 1,
            "ordinary_maximum_scaled_delta": maximum_equivalence,
            "ordinary_equivalence": equivalence,
            "short_frequency_results": short_results,
            "short_reference_frequency_mhz": SHORT_FREQUENCIES_MHZ[1],
            "nec2_deep_short_drift": nec2_drift,
            "nec4_deep_short_drift": nec4_drift,
            "conditioning_improvement": nec2_drift/nec4_drift,
        }
        destination = ROOT / "tests" / "nec4-similarity" / "m2-acceptance.json"
        if os.environ.get("NEC4_M2_RECORD") == "1":
            destination.write_text(json.dumps(evidence, indent=2,
                                              sort_keys=True)+"\n",
                                   encoding="utf-8")

        if os.environ.get("NEC4_M2_PLOTS") == "1":
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            output = ROOT / "tests" / "nec4-similarity"
            fig, axes = plt.subplots(1, 2, figsize=(10, 4),
                                     constrained_layout=True)
            for axis, pair, title in ((axes[0], dipole_documents, "Dipole"),
                                      (axes[1], loop_documents, "Square loop")):
                for engine, style in (("nec2", "-"), ("nec4", "--")):
                    values = currents(pair[engine])
                    axis.plot(range(1, len(values)+1), values, style,
                              label=engine)
                axis.set_title(title)
                axis.set_xlabel("ordered segment")
                axis.set_ylabel("current magnitude (A)")
                axis.legend()
            fig.savefig(output / "m2-current-equivalence.png", dpi=160)
            plt.close(fig)

            fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
            for engine, marker in (("nec2", "o"), ("nec4", "s")):
                # The 300 MHz point is an ordinary half-wave-scale antenna,
                # not part of the quasi-static -X*f plateau under test.
                rows = short_results[engine][1:]
                axis.plot([row["frequency_mhz"] for row in rows],
                          [row["capacitive_scale"]/reference for row in rows],
                          marker=marker, label=engine)
            axis.set_xscale("log")
            axis.set_xlabel("frequency (MHz)")
            axis.set_ylabel(r"normalized $-X f$")
            axis.set_title("Electrically short fixed dipole conditioning")
            axis.axhline(1.0, color="black", linewidth=0.8)
            axis.legend()
            fig.savefig(output / "m2-short-conditioning.png", dpi=160)
            plt.close(fig)

    print(f"M2 ordinary maximum scaled delta {maximum_equivalence:.3e}")
    print(f"M2 deep-short drift: nec2={nec2_drift:.3%}, "
          f"nec4={nec4_drift:.3%}, improvement={nec2_drift/nec4_drift:.1f}x")
    print("M2 representation equivalence and short-segment conditioning passed")


if __name__ == "__main__":
    main()
