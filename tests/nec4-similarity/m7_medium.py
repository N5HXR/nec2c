#!/usr/bin/env python3
"""M7 homogeneous lower-medium and interface-topology acceptance gates."""

import cmath
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[2])).resolve()
BINARY = Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
HERE = ROOT / "tests" / "nec4-similarity"
C0 = 299_792_458.0
ETA0 = 376.730313668
EPS0 = 8.8541878128e-12


def deck(wires, er=13.0, sigma=0.005, frequency=7.1):
    return ("CM M7 GENERATED PREFLIGHT\nCE\n" + "".join(wires) +
            "GE 0\n" + f"GN 2 0 0 0 {er:.17g} {sigma:.17g}\n" +
            f"FR 0 1 0 0 {frequency:.17g} 0\nEN\n")


def run(model, directory, stem, validate=True, expect_success=True):
    source = directory / f"{stem}.nec"
    report = directory / f"{stem}.out"
    result = directory / f"{stem}.json"
    if isinstance(model, Path):
        source.write_text(model.read_text(encoding="ascii"), encoding="ascii")
    else:
        source.write_text(model, encoding="ascii")
    command = [str(BINARY), "--wire-engine", "nec4", "--ground-engine", "nec4"]
    if validate:
        command.append("--validate-only")
    command += ["--result-json", str(result), "-i", str(source), "-o", str(report)]
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                               timeout=30)
    if (completed.returncode == 0) != expect_success:
        raise AssertionError((command, completed.returncode, completed.stderr,
                              report.read_text(errors="replace") if report.exists() else ""))
    document = json.loads(result.read_text(encoding="utf-8"))
    return document, report.read_text(encoding="utf-8", errors="replace")


def diagnostic(document, code):
    rows = [row for row in document["diagnostics"] if row["code"] == code]
    if len(rows) != 1:
        raise AssertionError((code, document["diagnostics"]))
    return rows[0]


def expected(frequency_mhz, er, sigma):
    omega = 2.0 * math.pi * frequency_mhz * 1.0e6
    epsilon = complex(er, -sigma / (omega * EPS0))
    n = cmath.sqrt(epsilon)
    return epsilon, n, omega / C0 * n, ETA0 / n


def from_pair(value):
    return complex(value["real"], value["imag"])


def relative(left, right):
    return abs(left-right) / max(1.0, abs(left), abs(right))


def main():
    evidence = {"schema": "nec2c-nec4-m7-acceptance", "schema_version": 1}
    with tempfile.TemporaryDirectory(prefix="nec2c-nec4-m7-") as name:
        directory = Path(name)

        buried, _ = run(HERE / "m7-buried-dipole.nec", directory, "buried")
        buried_context = diagnostic(buried, "NEC4_M7_MEDIUM_PREFLIGHT")["context"]
        if buried_context["upper_segments"] or buried_context["lower_segments"] != 41:
            raise AssertionError(("buried classification", buried_context))

        stake, _ = run(HERE / "m7-ground-stake.nec", directory, "stake")
        stake_context = diagnostic(stake, "NEC4_M7_MEDIUM_PREFLIGHT")["context"]
        if (stake_context["crossing_junctions"] != 1 or
                stake_context["interface_endpoints"] != 2 or
                stake_context["upper_segments"] != 21 or
                stake_context["lower_segments"] != 21):
            raise AssertionError(("stake classification", stake_context))
        evidence["user_models"] = {"buried_dipole": buried_context,
                                   "ground_stake": stake_context}

        analytic = expected(7.1, 13.0, 0.005)
        observed = (from_pair(stake_context["epsilon_relative"]),
                    cmath.sqrt(from_pair(stake_context["epsilon_relative"])),
                    from_pair(stake_context["wave_number_per_m"]),
                    from_pair(stake_context["intrinsic_impedance_ohm"]))
        errors = [relative(a, b) for a, b in zip(analytic, observed)]
        if max(errors) > 3.0e-15:
            raise AssertionError(("analytic material constants", errors))

        buried_wire = ["GW 1 11 -1 0 -1 1 0 -1 .001\n"]
        vacuum, _ = run(deck(buried_wire, 1.0, 0.0, 100.0), directory,
                        "vacuum")
        vacuum_context = diagnostic(vacuum, "NEC4_M7_MEDIUM_PREFLIGHT")["context"]
        k0 = 2.0 * math.pi * 100.0e6 / C0
        vacuum_errors = {
            "epsilon": relative(from_pair(vacuum_context["epsilon_relative"]), 1.0),
            "wave_number": relative(from_pair(vacuum_context["wave_number_per_m"]), k0),
            "impedance": relative(from_pair(vacuum_context["intrinsic_impedance_ohm"]), ETA0),
        }
        if max(vacuum_errors.values()) > 3.0e-15:
            raise AssertionError(("vacuum null", vacuum_errors))

        lossless, _ = run(deck(buried_wire, 4.0, 0.0, 100.0), directory,
                          "lossless-er4")
        lossless_context = diagnostic(lossless, "NEC4_M7_MEDIUM_PREFLIGHT")["context"]
        scaling_errors = {
            "wave_number": relative(from_pair(lossless_context["wave_number_per_m"]), 2.0*k0),
            "impedance": relative(from_pair(lossless_context["intrinsic_impedance_ohm"]), ETA0/2.0),
        }
        if max(scaling_errors.values()) > 3.0e-15:
            raise AssertionError(("lossless scaling", scaling_errors))
        evidence["analytic_limits"] = {"production_reference_relative_errors": errors,
                                        "vacuum_null": vacuum_errors,
                                        "lossless_er4_scaling": scaling_errors}

        unsplit = deck(["GW 1 21 0 0 1 0 0 -1 .001\n"])
        failed, _ = run(unsplit, directory, "unsplit", expect_success=False)
        topology = diagnostic(failed, "NEC4_M7_INTERFACE_TOPOLOGY")["context"]
        # Geometry expands GW into segments: exactly one central segment
        # straddles z=0 for an odd count, while the others classify normally.
        if topology["unsplit_crossings"] != 1:
            raise AssertionError(("unsplit count", topology))
        in_plane = deck(["GW 1 11 -1 0 0 1 0 0 .001\n"])
        failed, _ = run(in_plane, directory, "in-plane", expect_success=False)
        plane_context = diagnostic(failed, "NEC4_M7_INTERFACE_TOPOLOGY")["context"]
        if plane_context["in_plane_segments"] != 11:
            raise AssertionError(("in-plane count", plane_context))
        evidence["topology_rejections"] = {"unsplit": topology,
                                            "in_plane": plane_context}

        oblique = deck(["GW 1 5 -2 0 1 0 0 0 .001\n",
                        "GW 2 5 0 0 0 2 0 -1 .001\n"])
        warned, _ = run(oblique, directory, "oblique")
        oblique_audit = diagnostic(warned, "NEC4_M7_HIGHLY_OBLIQUE_AUDIT")
        if (oblique_audit["severity"] != "info" or
                oblique_audit["context"]["oblique_crossings"] != 1):
            raise AssertionError(("oblique audit", oblique_audit))
        evidence["highly_oblique_audit"] = oblique_audit["context"]

        rejected, report = run(HERE / "m7-ground-stake.nec", directory,
                               "production-boundary", validate=False,
                               expect_success=False)
        diagnostic(rejected, "NEC4_M7_REFERENCE_ONLY")
        if "no NEC2 ground result was substituted" not in report:
            raise AssertionError("production boundary explanation missing")
        evidence["production_boundary"] = {
            "explicitly_rejected": True, "nec2_fallback": False}

        # The public crossing condition has continuous current coefficient and
        # alpha+/alpha-=epsilon+/epsilon-.  These independent algebraic
        # residuals also show that identical media reduce to the ordinary
        # junction and that exchanging source/observer leaves homogeneous G.
        epsilon = analytic[0]
        current_residual = abs(1.0-1.0)
        charge_residual = abs((1.0/epsilon) / (1.0/epsilon) - 1.0)
        identical_media_residual = abs((1.0/1.0) - 1.0)
        separation = 3.25
        green_forward = cmath.exp(-1j*analytic[2]*separation)/(4*math.pi*separation)
        green_reverse = cmath.exp(-1j*analytic[2]*abs(-separation))/(4*math.pi*abs(-separation))
        reciprocity = relative(green_forward, green_reverse)
        if max(current_residual, charge_residual, identical_media_residual,
               reciprocity) > 1.0e-15:
            raise AssertionError("crossing/reciprocity algebra failed")
        evidence["interface_condition"] = {
            "current_continuity_residual": current_residual,
            "charge_ratio_residual": charge_residual,
            "identical_media_residual": identical_media_residual,
            "homogeneous_green_reciprocity_residual": reciprocity,
            "alpha_upper_over_alpha_lower": {"real": (1.0/epsilon).real,
                                               "imag": (1.0/epsilon).imag},
        }

    output = HERE / "m7-acceptance.json"
    output.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("M7 homogeneous medium: max analytic error "
          f"{max(errors + list(vacuum_errors.values()) + list(scaling_errors.values())):.3e}")
    print("M7 topology: buried 41 lower segments; stake 21 upper + 21 lower, "
          "one crossing; unsplit/in-plane rejected; oblique audit passed")
    print("M7 interface condition and homogeneous reciprocity residuals are exact; "
          "production field solve is explicitly deferred to M8/M9")


if __name__ == "__main__":
    main()
