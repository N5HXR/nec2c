#!/usr/bin/env python3
"""Run and index the Milestone 9 release-candidate validation matrix.

The runner does not build the preserved M0 revision.  Pass an executable built
from commit 1166bc3bc5ba490c56fadf946070998a47ff8e00 with
``--reference-binary`` so the old implementation remains physically separate
from the optimized checkout.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_COMMIT = "1166bc3bc5ba490c56fadf946070998a47ff8e00"
SMALL_AUTO_CASES = {
    "classic-small", "classic-mixed", "dielectric-dense", "dielectric-fft",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_step(name, command, artifacts, steps):
    print(f"[{name}] {' '.join(map(str, command))}", flush=True)
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    step = {
        "name": name,
        "command": [str(item) for item in command],
        "elapsed_seconds": time.perf_counter() - started,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "artifacts": [],
    }
    for artifact in artifacts:
        if artifact.exists():
            step["artifacts"].append({
                "path": str(artifact.resolve()),
                "sha256": sha256(artifact),
            })
    steps.append(step)
    if completed.stdout:
        print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n")
    if completed.stderr:
        print(completed.stderr, file=sys.stderr,
              end="" if completed.stderr.endswith("\n") else "\n")
    return completed.returncode == 0


def benchmark_command(binary, output, repeat, warmup, controls):
    return [
        sys.executable, "tools/performance_benchmark.py",
        "--binary", str(binary), "--repeat", str(repeat),
        "--warmup", str(warmup), "--output", str(output), *controls,
    ]


def numerical_hashes(report_path):
    value = json.loads(report_path.read_text(encoding="utf-8"))
    return {
        case["id"]: sorted({run.get("physics_sha256", run["numerical_sha256"])
                            for run in case["runs"]})
        for case in value["cases"]
    }


def corpus_summary(report_path):
    value = json.loads(report_path.read_text(encoding="utf-8"))
    result = {}
    for case in value["cases"]:
        summary = case["summary"]
        cli_wall = summary["external_wall_seconds"]
        solver_wall = summary["wall_seconds"]
        rss = summary["process_tree_peak_rss_kib"]
        result[case["id"]] = {
            "cli_wall_seconds": cli_wall["median"],
            "solver_wall_seconds": solver_wall["median"] if solver_wall else None,
            "solver_wall_range": ([solver_wall["minimum"], solver_wall["maximum"]]
                                  if solver_wall else None),
            "peak_rss_kib": rss["maximum"] if rss else None,
            "stable_numerical_result": summary["stable_numerical_result"],
        }
    return result


def compare_corpora(auto_path, legacy_path, reference_path=None):
    auto = corpus_summary(auto_path)
    legacy = corpus_summary(legacy_path)
    reference = corpus_summary(reference_path) if reference_path else None
    comparisons = []
    for case_id in sorted(auto):
        item = {
            "case": case_id,
            "auto": auto[case_id],
            "current_legacy_controls": legacy[case_id],
            "legacy_over_auto_solver_speedup": (
                legacy[case_id]["solver_wall_seconds"] /
                auto[case_id]["solver_wall_seconds"]),
        }
        if reference is not None:
            item["m0_reference"] = reference[case_id]
            item["m0_over_auto_cli_speedup"] = (
                reference[case_id]["cli_wall_seconds"] /
                auto[case_id]["cli_wall_seconds"])
        comparisons.append(item)
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--reference-binary", type=Path,
                        help="executable built from the immutable M0 commit")
    parser.add_argument("--aegnec2", type=Path,
                        help="optional external AEG executable for cross-solver validation")
    parser.add_argument("--momwire-root", type=Path,
                        default=Path.home() / "src/momwire")
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--skip-check", action="store_true")
    parser.add_argument("--full", action="store_true",
                        help="also run the M4B-M8 phase/scaling benchmarks")
    arguments = parser.parse_args()
    if arguments.repeat < 1 or arguments.warmup < 0:
        parser.error("--repeat must be positive and --warmup nonnegative")

    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    binary = arguments.binary.resolve()
    if not binary.is_file():
        parser.error(f"current binary does not exist: {binary}")
    if arguments.reference_binary:
        arguments.reference_binary = arguments.reference_binary.resolve()
        if not arguments.reference_binary.is_file():
            parser.error(f"reference binary does not exist: {arguments.reference_binary}")

    steps = []
    passed = True
    if not arguments.skip_check:
        passed &= run_step("complete-test-suite", ["make", "check"], [], steps)

    auto = output_directory / "current-auto.json"
    passed &= run_step(
        "current-auto-corpus",
        benchmark_command(binary, auto, arguments.repeat, arguments.warmup,
                          ["--backend", "auto", "--fft-backend", "auto",
                           "--jobs", "auto", "--threads", "auto"]),
        [auto], steps)

    legacy = output_directory / "current-legacy-controls.json"
    passed &= run_step(
        "current-legacy-controls-corpus",
        benchmark_command(binary, legacy, arguments.repeat, arguments.warmup,
                          ["--backend", "builtin", "--fft-backend", "builtin",
                           "--jobs", "1", "--threads", "1"]),
        [legacy], steps)

    if arguments.reference_binary:
        reference = output_directory / "m0-reference.json"
        passed &= run_step(
            "preserved-m0-corpus",
            benchmark_command(arguments.reference_binary, reference,
                              arguments.repeat, arguments.warmup, ["--legacy-cli"]),
            [reference], steps)

    if arguments.aegnec2:
        cross = output_directory / "cross-solver.json"
        command = [
            sys.executable, "tools/cross_solver_compare.py",
            "--nec2c", str(binary), "--nec2c-backend", "auto",
            "--aegnec2", str(arguments.aegnec2.resolve()),
            "--momwire-root", str(arguments.momwire_root.resolve()),
            "--output", str(cross),
        ]
        passed &= run_step("cross-solver", command, [cross], steps)

    if arguments.full:
        phase_commands = [
            ("fft-backends", "tools/fft_backend_compare.py",
             ["--model", "tree", "--threads", "1,16,32",
              "--max-iterations", "300"], "fft-backends.json"),
            ("assembly", "tools/m5_assembly_benchmark.py",
             ["--threads", "1,2,4,8,16,32"], "m5-assembly.json"),
            ("fill", "tools/m6_fill_benchmark.py",
             ["--threads", "1,2,4,8,16,32"], "m6-fill.json"),
            ("coupling", "tools/m7_coupling_benchmark.py",
             ["--threads", "1,2,4,8,16,32"], "m7-coupling.json"),
            ("fields", "tools/m8_field_benchmark.py",
             ["--threads", "1,2,4,8,16,32"], "m8-fields.json"),
        ]
        for name, program, options, filename in phase_commands:
            artifact = output_directory / filename
            command = [sys.executable, program, *options,
                       "--warmup", str(arguments.warmup),
                       "--repeat", str(arguments.repeat),
                       "--output", str(artifact)]
            passed &= run_step(name, command, [artifact], steps)

    legacy_hashes = numerical_hashes(legacy) if legacy.exists() else None
    reference_hashes = (numerical_hashes(reference)
                        if arguments.reference_binary and reference.exists() else None)
    stable_legacy = (legacy_hashes is not None and
                     all(len(values) == 1 for values in legacy_hashes.values()))
    auto_hashes = numerical_hashes(auto) if auto.exists() else None
    stable_auto = (auto_hashes is not None and
                   all(len(values) == 1 for values in auto_hashes.values()))
    legacy_reference_exact = (
        legacy_hashes == reference_hashes if reference_hashes is not None else None
    )
    performance_comparison = (
        compare_corpora(auto, legacy,
                        reference if arguments.reference_binary else None)
        if auto.exists() and legacy.exists() and
        (not arguments.reference_binary or reference.exists()) else None
    )
    small_auto_gates = []
    for item in performance_comparison or []:
        if item["case"] not in SMALL_AUTO_CASES:
            continue
        auto_time = item["auto"]["solver_wall_seconds"]
        legacy_time = item["current_legacy_controls"]["solver_wall_seconds"]
        auto_range = item["auto"]["solver_wall_range"]
        legacy_range = item["current_legacy_controls"]["solver_wall_range"]
        overlap = auto_range[0] <= legacy_range[1] and legacy_range[0] <= auto_range[1]
        regression = auto_time / legacy_time - 1.0
        small_auto_gates.append({
            "case": item["case"],
            "auto_over_legacy_regression": regression,
            "sample_ranges_overlap": overlap,
            "passed": regression <= 0.05 or overlap,
        })
    small_auto_gate_passed = (len(small_auto_gates) == len(SMALL_AUTO_CASES) and
                              all(item["passed"] for item in small_auto_gates))
    report = {
        "schema": "nec2c-m9-release-candidate",
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "reference_commit": REFERENCE_COMMIT,
        "host": {
            "platform": platform.platform(),
            "logical_cpus": os.cpu_count(),
        },
        "current_binary": {"path": str(binary), "sha256": sha256(binary)},
        "reference_binary": ({
            "path": str(arguments.reference_binary),
            "sha256": sha256(arguments.reference_binary),
        } if arguments.reference_binary else None),
        "settings": {
            "repeat": arguments.repeat,
            "warmup": arguments.warmup,
            "full": arguments.full,
        },
        "required_evidence_present": {
            "complete_test_suite": not arguments.skip_check,
            "preserved_reference": arguments.reference_binary is not None,
            "cross_solver": arguments.aegnec2 is not None,
            "phase_matrix": arguments.full,
        },
        "legacy_reference_comparison": {
            "stable_automatic_defaults": stable_auto,
            "stable_current_legacy": stable_legacy,
            "exact_canonical_result_hashes": legacy_reference_exact,
            "note": ("Exact hashes apply to the explicitly selected built-in, "
                     "single-worker path. Optimized LAPACK/FFTW results are "
                     "validated by tolerance-based backend tests instead."),
        },
        "performance_comparison": performance_comparison,
        "small_auto_selection_gate": {
            "limit": "no more than 5% slower unless measured ranges overlap",
            "basis": "solver timing: automatic defaults versus current built-in single-worker controls",
            "passed": small_auto_gate_passed,
            "cases": small_auto_gates,
        },
        "m0_cli_startup_note": (
            "M0 comparisons use total CLI wall because M0 predates timing JSON. "
            "They disclose dynamic-library/bootstrap cost but do not replace the "
            "current auto-selection crossover gate."),
        "steps": steps,
        "commands_passed": passed,
        "acceptance_ready": passed and stable_auto and stable_legacy and
        small_auto_gate_passed and
        legacy_reference_exact is True and
        not arguments.skip_check and
        arguments.reference_binary is not None and arguments.aegnec2 is not None and
        arguments.full,
    }
    index = output_directory / "m9-release-candidate.json"
    temporary = index.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    temporary.replace(index)
    print(index)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
