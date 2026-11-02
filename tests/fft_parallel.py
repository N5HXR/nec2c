#!/usr/bin/env python3
"""Require bit-identical built-in FFT products across thread counts."""

import argparse
import os
import json
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = Path(os.environ.get("FFT_BENCHMARK", ROOT / "tools/m9_fft_benchmark"))
PROGRAM = Path(os.environ.get("NEC2C", ROOT / "nec2c"))
DECK = ROOT / "tests/dielectric/m3-frequency-sweep.nec"

parser = argparse.ArgumentParser()
parser.add_argument(
    "--require-parallel",
    action="store_true",
    help="fail when the benchmark was built without OpenMP",
)
arguments = parser.parse_args()


def run(threads):
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(threads)
    completed = subprocess.run(
        [str(BENCHMARK), "tree", "1e-8", "1", "builtin"],
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        check=False,
    )
    # One iteration intentionally reaches the documented iteration-limit exit.
    if completed.returncode != 1:
        raise AssertionError(
            f"{threads}-thread benchmark returned {completed.returncode}: "
            f"{completed.stderr}"
        )
    fields = dict(re.findall(r"([a-z_]+)=([^\s]+)", completed.stdout))
    required = {
        "matvec_checksum_real",
        "matvec_checksum_imag",
        "matvec_bitwise_hash",
        "residual",
        "threads",
    }
    if not required.issubset(fields):
        raise AssertionError(f"incomplete benchmark output: {completed.stdout}")
    return (
        int(fields["threads"]),
        tuple(fields[name] for name in sorted(required - {"threads"})),
    )


captured, reference = run(1)
if captured != 1:
    raise AssertionError(f"serial plan captured {captured} threads")
parallel_available = True
for count in (2, 4):
    captured, signature = run(count)
    if captured == 1:
        parallel_available = False
    elif captured != count:
        raise AssertionError(f"requested {count} threads, plan captured {captured}")
    if signature != reference:
        raise AssertionError(f"{count}-thread FFT product differs from serial")

if parallel_available:
    print("parallel built-in FFT is bit-identical at 1, 2, and 4 threads")
else:
    if arguments.require_parallel:
        raise SystemExit(
            "parallel FFT support is unavailable in this build; run "
            "`autoreconf -fi && ./configure --with-lapack=auto && make -j8`, "
            "and verify configure reports `parallel built-in FFT: enabled`"
        )
    print("OpenMP unavailable: built-in FFT serial fallback is bit-identical")


def run_cli(directory, label, jobs, threads):
    report = directory / f"{label}.out"
    result = directory / f"{label}.json"
    timing = directory / f"{label}-timing.json"
    subprocess.run(
        [
            str(PROGRAM),
            "-s", "fft",
            "--fft-backend", "builtin",
            "-r", "1e-8",
            "--jobs", jobs,
            "--threads", threads,
            "--result-json", str(result),
            "--timing-json", str(timing),
            "-i", str(DECK),
            "-o", str(report),
        ],
        cwd=ROOT,
        check=True,
        timeout=30,
    )
    result_document = json.loads(result.read_text())
    for job in result_document["jobs"]:
        for execution in job["executions"]:
            execution["solver"]["elapsed_seconds"] = 0
    return report.read_text(), result_document, json.loads(timing.read_text())


with tempfile.TemporaryDirectory(prefix="nec2c-m4a-") as temporary:
    directory = Path(temporary)
    serial_report, serial_result, serial_timing = run_cli(
        directory, "serial", "1", "1"
    )
    threaded_report, threaded_result, threaded_timing = run_cli(
        directory, "threaded", "1", "4"
    )
    automatic_report, automatic_result, automatic_timing = run_cli(
        directory, "automatic", "auto", "auto"
    )
    assert threaded_result == serial_result
    assert automatic_result == serial_result
    for report in (serial_report, threaded_report, automatic_report):
        frequencies = [
            float(value)
            for value in re.findall(r"FREQUENCY \(MHZ\):\s+([0-9.]+)", report)
        ]
        assert frequencies == [90.0, 100.0, 110.0, 120.0]
        assert "FFT ACCOUNTED APPLICATION BYTES:" in report
        assert "FFT THREAD BUDGET:" in report
    assert serial_timing["concurrency"]["jobs_used"] == 1
    assert serial_timing["concurrency"]["threads_per_job"] == 1
    assert threaded_timing["concurrency"]["jobs_used"] == 1
    assert threaded_timing["concurrency"]["threads_per_job"] == 4
    automatic = automatic_timing["concurrency"]
    assert automatic["jobs_used"] * automatic["threads_per_job"] <= automatic["allowed_cpus"]
    if not parallel_available:
        assert automatic["threads_per_job"] == 1

print("nec2c FFT CLI preserves ordered results and the CPU budget")
