#!/usr/bin/env python3
"""M3 local junction solve, convergence, reciprocity, and crash acceptance."""

import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import mpmath as mp


ROOT = Path(os.environ.get("TEST_SRCDIR", ".")).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
LOCAL_TEST = ROOT / "tests" / "junction_charge_test"


def scaled(a, b):
    return abs(a-b) / max(1.0, abs(a), abs(b))


def oracle_factors(length, radius, direction):
    mp.mp.dps = 80
    count = len(length)
    order = 2*count
    length = list(map(mp.mpf, length))
    radius = list(map(mp.mpf, radius))
    direction = [[mp.mpf(value) for value in row] for row in direction]

    def dot(left, right):
        return sum(a*b for a, b in zip(left, right))

    wave_number = 2*mp.pi

    def triangle(observer, source):
        def integrand(position):
            delta = [observer[d]-position*direction[source][d]
                     for d in range(3)]
            distance = mp.sqrt(dot(delta, delta)+radius[source]**2)
            return ((length[source]-position)*
                    mp.exp(-1j*wave_number*distance)/distance)
        return mp.quad(integrand, [0, length[source]])

    def tail(observer, source):
        axial = dot(observer, direction[source])
        rho = mp.sqrt(dot(observer, observer)+radius[source]**2-axial**2)

        def integrand(position):
            distance = mp.sqrt(rho**2+(position-axial)**2)
            return mp.exp(-1j*wave_number*distance)/distance

        finite = mp.mpc(0)
        if axial:
            finite = mp.quad(integrand, [0, axial])
        return finite-.5j*mp.pi*mp.hankel2(0, wave_number*rho)

    matrix = mp.matrix(order, order)
    rhs = mp.matrix(order, 1)
    row = 0
    for observer_wire in range(count):
        for match in (mp.mpf("0.5"), mp.mpf(1)):
            observer = [match*length[observer_wire]*value
                        for value in direction[observer_wire]]
            for source in range(count):
                matrix[row, 2*source] = triangle(observer, source)
                matrix[row, 2*source+1] = length[source]*tail(observer, source)
            rhs[row] = 1
            row += 1
    solution = mp.lu_solve(matrix, rhs)
    factors = [(solution[2*i]+solution[2*i+1])*length[i]
               for i in range(count)]
    mean = sum(factors)/count
    return [complex(value/mean) for value in factors]


def local_oracle_check():
    completed = subprocess.run([str(LOCAL_TEST)], cwd=ROOT, text=True,
                               capture_output=True, check=True)
    production = json.loads(completed.stdout)
    for name in ("step", "tee", "star"):
        production[name] = [complex(*value) for value in production[name]]
    cases = {
        "step": ([0.0125, 0.0125], [0.001, 0.005],
                 [[0, 0, -1], [0, 0, 1]]),
        "tee": ([0.02]*3, [0.001]*3,
                [[1, 0, 0], [-1, 0, 0], [0, 0, 1]]),
        "star": ([0.02]*4, [0.001, 0.0015, 0.002, 0.003],
                 [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]),
    }
    maximum = 0.0
    references = {}
    for name, arguments in cases.items():
        reference = oracle_factors(*arguments)
        references[name] = reference
        maximum = max(maximum, *(scaled(a, b) for a, b in
                                 zip(production[name], reference)))
    if maximum > 1.0e-11:
        raise AssertionError(("local junction oracle", maximum,
                              production, references))
    return maximum, production, references


def step_deck(segments, source_segment, frequency=299.792458):
    return ("CM M3 STEPPED-RADIUS CONVERGENCE\nCE\n"
            f"GW 1 {segments} 0 0 -0.25 0 0 0 0.0005\n"
            f"GW 2 {segments} 0 0 0 0 0 0.25 0.002\n"
            "GE 0\n"
            f"FR 0 1 0 0 {frequency:.12g} 0\n"
            f"EX 0 0 {source_segment} 0 1 0\nXQ\nEN\n")


def run_deck(text, engine, directory, stem):
    deck = directory / f"{stem}.nec"
    report = directory / f"{stem}-{engine}.out"
    result = directory / f"{stem}-{engine}.json"
    deck.write_text(text, encoding="ascii")
    command = [str(BINARY), "--wire-engine", engine, "--ground-engine",
               "nec2", "--backend", "builtin", "--jobs", "1",
               "--threads", "1", "--result-json", str(result),
               "-i", str(deck), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               capture_output=True)
    if completed.returncode:
        raise AssertionError((command, completed.returncode,
                              completed.stdout, completed.stderr))
    document = json.loads(result.read_text(encoding="utf-8"))
    return document, report.read_text(encoding="utf-8", errors="replace")


def feed_impedance(document):
    value = document["jobs"][0]["executions"][0]["feeds"][0]["impedance_ohm"]
    return complex(value["re"], value["im"])


def wire_current(document, segment):
    row = document["jobs"][0]["executions"][0]["wire_currents"][segment-1]
    value = row["current_a"]
    return complex(value["re"], value["im"])


def main():
    oracle_error, production_factors, reference_factors = local_oracle_check()
    # Keep the legacy reduced kernel in its pre-M4 length/radius envelope.
    # Odd counts place the source at the same physical point (-0.125 m).
    counts = [5, 9, 13, 17]
    convergence = {engine: [] for engine in ("nec2", "nec4")}
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m3-") as name:
        temporary = Path(name)
        for count in counts:
            source = (count+1)//2
            text = step_deck(count, source)
            for engine in convergence:
                document, report = run_deck(text, engine, temporary,
                                            f"step-{count}-{engine}")
                convergence[engine].append(feed_impedance(document))
                if engine == "nec4":
                    match = re.search(r"(\d+) FALLBACKS", report)
                    if not match or int(match.group(1)):
                        raise AssertionError(("junction fallback", count,
                                              report[-2000:]))

        # Reciprocal transfer admittance between corresponding segments on
        # opposite sides of the same unequal-radius junction.
        count = 20
        port_a, port_b = count//2, count+count//2+1
        reciprocity = {}
        for engine in ("nec2", "nec4"):
            first, _ = run_deck(step_deck(count, port_a), engine, temporary,
                                f"reciprocity-a-{engine}")
            second, _ = run_deck(step_deck(count, port_b), engine, temporary,
                                 f"reciprocity-b-{engine}")
            reciprocity[engine] = scaled(wire_current(first, port_b),
                                         wire_current(second, port_a))

        # The historical high-order center-junction crash must remain fixed.
        crash_result = temporary / "circular.json"
        crash_report = temporary / "circular.out"
        command = [str(BINARY), "--wire-engine", "nec4", "--ground-engine",
                   "nec2", "--backend", "builtin", "--jobs", "1",
                   "--threads", "1", "--result-json", str(crash_result),
                   "-i", str(ROOT / "circular-mesh-center-junction-segfault.nec"),
                   "-o", str(crash_report)]
        completed = subprocess.run(command, cwd=ROOT, text=True,
                                   capture_output=True)
        if completed.returncode or not crash_result.is_file():
            raise AssertionError(("circular junction regression", completed))
        crash_text = crash_report.read_text(encoding="utf-8", errors="replace")
        if not re.search(r"0 FALLBACKS", crash_text):
            raise AssertionError(("circular junction fallback", crash_text[-2000:]))

        topology_models = {}
        for fixture_name in ("m3-tee.nec", "m3-y.nec", "m3-star.nec"):
            fixture = ROOT / "tests" / "nec4-similarity" / fixture_name
            document, report = run_deck(fixture.read_text(encoding="ascii"),
                                        "nec4", temporary, fixture.stem)
            match = re.search(
                r"(\d+) JUNCTIONS, (\d+) LOCAL SOLVES, (\d+) SYMMETRIC, "
                r"(\d+) INTERFACE-PRESCRIBED, (\d+) FALLBACKS",
                report)
            if not match or int(match.group(5)):
                raise AssertionError(("topology fixture", fixture_name,
                                      report[-2000:]))
            topology_models[fixture_name] = {
                "junctions": int(match.group(1)),
                "local_solves": int(match.group(2)),
                "symmetric": int(match.group(3)),
                "feed_impedance_ohm": [feed_impedance(document).real,
                                        feed_impedance(document).imag],
            }

    errors = {engine: [scaled(value, values[-1]) for value in values]
              for engine, values in convergence.items()}
    if oracle_error > 1.0e-11:
        raise AssertionError(("oracle error", oracle_error))
    if production_factors["permutation_defect"] > 2.0e-12:
        raise AssertionError(("local permutation reciprocity",
                              production_factors["permutation_defect"]))
    # NEC's collocation voltage-gap source is not an exact reciprocal port
    # observable on unequal radii.  Retain it as a characterized end-to-end
    # check and require the new treatment to improve, rather than obscure, it.
    if (reciprocity["nec4"] > 1.0e-3 or
            reciprocity["nec4"] >= reciprocity["nec2"]):
        raise AssertionError(("source-model reciprocity", reciprocity))
    if not all(errors["nec4"][index] > errors["nec4"][index+1]
               for index in range(len(counts)-2)):
        raise AssertionError(("NEC4 step refinement did not converge",
                              errors["nec4"]))
    if not all(errors["nec4"][index] < errors["nec2"][index]
               for index in range(len(counts)-1)):
        raise AssertionError(("NEC4 step refinement did not improve",
                              errors))

    evidence = {
        "schema": "nec2c-nec4-m3-acceptance", "schema_version": 1,
        "local_oracle_maximum_scaled_error": oracle_error,
        "production_factors": {
            **production_factors,
            **{name: [[value.real, value.imag] for value in
                       production_factors[name]]
               for name in ("step", "tee", "star")},
        },
        "reference_factors": {
            name: [[value.real, value.imag] for value in values]
            for name, values in reference_factors.items()
        },
        "local_permutation_defect": production_factors["permutation_defect"],
        "source_model_reciprocity_defect": reciprocity,
        "segment_counts_per_arm": counts,
        "impedance_ohm": {engine: [[value.real, value.imag]
                                    for value in values]
                          for engine, values in convergence.items()},
        "scaled_error_to_nec4_17": errors,
        "circular_center_junction": "passed without fallback",
        "topology_models": topology_models,
    }
    if os.environ.get("NEC4_M3_RECORD") == "1":
        destination = ROOT / "tests" / "nec4-similarity" / "m3-acceptance.json"
        destination.write_text(json.dumps(evidence, indent=2,
                                          sort_keys=True)+"\n",
                               encoding="utf-8")
    if os.environ.get("NEC4_M3_PLOTS") == "1":
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/nec2c-matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
        for engine, marker in (("nec2", "o"), ("nec4", "s")):
            axis.semilogy(counts[:-1], errors[engine][:-1], marker=marker,
                          label=engine)
        axis.set_xlabel("segments per dipole arm")
        axis.set_ylabel("scaled complex-impedance error to 17/arm")
        axis.set_title("Stepped-radius junction convergence")
        axis.grid(True, which="both", alpha=.3)
        axis.legend()
        fig.savefig(ROOT / "tests" / "nec4-similarity" /
                    "m3-junction-convergence.png", dpi=160)
        plt.close(fig)
    print(f"M3 local junction oracle maximum scaled error {oracle_error:.3e}")
    print(f"M3 local permutation defect "
          f"{production_factors['permutation_defect']:.3e}")
    print(f"M3 source-model reciprocity defect: nec2={reciprocity['nec2']:.3e}, "
          f"nec4={reciprocity['nec4']:.3e}")
    print("M3 step convergence (NEC4 scaled to 17/arm): " +
          ", ".join(f"{count}={error:.3e}" for count, error in
                    zip(counts, errors["nec4"])))
    print("M3 junction solve, convergence, reciprocity, and crash gates passed")


if __name__ == "__main__":
    main()
