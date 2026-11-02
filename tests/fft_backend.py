#!/usr/bin/env python3
"""Validate optional FFTW selection and numerical equivalence."""

import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(os.environ.get("TEST_SRCDIR", Path(__file__).resolve().parents[1]))
PROGRAM = Path(os.environ.get("NEC2C", ROOT / "nec2c"))
BENCHMARK = Path(os.environ.get("FFT_BENCHMARK", ROOT / "tools/m9_fft_benchmark"))


def invoke(command, expect=0, environment=None):
    completed = subprocess.run(command, cwd=ROOT, env=environment, text=True,
                               capture_output=True, timeout=60)
    if completed.returncode != expect:
        raise AssertionError((completed.returncode, completed.stdout, completed.stderr))
    return completed


def compare(left, right, path="$"):
    if isinstance(left, bool) or left is None or isinstance(left, str):
        assert left == right, (path, left, right)
    elif isinstance(left, (int, float)):
        assert isinstance(right, (int, float)), (path, left, right)
        if path.endswith(".accounted_storage_bytes"):
            assert abs(left - right) <= max(left, right) // 100, (path, left, right)
            return
        assert math.isclose(left, right, rel_tol=2e-9, abs_tol=2e-11), (
            path, left, right
        )
    elif isinstance(left, list):
        assert isinstance(right, list) and len(left) == len(right), path
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            compare(left_item, right_item, f"{path}[{index}]")
    else:
        assert isinstance(right, dict) and left.keys() == right.keys(), path
        for key in left:
            compare(left[key], right[key], f"{path}.{key}")


capabilities = json.loads(invoke([str(PROGRAM), "--capabilities-json"]).stdout)
available = "fftw" in capabilities["features"]["fft_backends"]
if not available:
    fallback_environment = os.environ.copy()
    fallback_environment["OMP_NUM_THREADS"] = "1"
    fallback = invoke([str(BENCHMARK), "skeleton", "1e-8", "1", "auto"],
                      expect=1, environment=fallback_environment)
    assert "backend=builtin" in fallback.stdout
    rejected = invoke([str(PROGRAM), "--fft-backend", "fftw", "-i",
                       str(ROOT / "tests/dielectric/m9-tree.nec"), "-o", os.devnull],
                      expect=1)
    assert "unavailable" in rejected.stderr
    print("FFTW unavailable; explicit-backend rejection passed")
    raise SystemExit(0)

environment = os.environ.copy()
environment["OMP_NUM_THREADS"] = "4"
observed = {}
for backend in ("builtin", "fftw", "auto"):
    run = invoke([str(BENCHMARK), "skeleton", "1e-8", "300", backend],
                 environment=environment)
    values = dict(re.findall(r"([a-z0-9_]+)=([^\s]+)", run.stdout))
    assert values["backend"] == ("fftw" if backend == "auto" else backend)
    observed[backend] = values
assert observed["auto"]["backend"] == "fftw"

for name in ("residual", "integrated_normalized_current_0",
             "integrated_normalized_current_1", "integrated_normalized_current_2",
             "extinction_normalized", "scattered16_normalized",
             "scattered32_normalized", "power_relative_defect"):
    left = complex(observed["builtin"][name].replace("j", "j"))
    right = complex(observed["fftw"][name].replace("j", "j"))
    assert math.isclose(left.real, right.real, rel_tol=2e-10, abs_tol=2e-12), name
    assert math.isclose(left.imag, right.imag, rel_tol=2e-10, abs_tol=2e-12), name

with tempfile.TemporaryDirectory(prefix="nec2c-fftw-") as temporary:
    directory = Path(temporary)
    results = {}
    for backend in ("builtin", "fftw"):
        output = directory / f"{backend}.out"
        result = directory / f"{backend}.json"
        invoke([str(PROGRAM), "-s", "fft", "--fft-backend", backend,
                "--jobs", "1", "--threads", "4", "--result-json", str(result),
                "-i", str(ROOT / "tests/dielectric/m9-tree.nec"), "-o", str(output)])
        assert f"FFT BACKEND: {backend}" in output.read_text()
        document = json.loads(result.read_text())
        for job in document["jobs"]:
            for execution in job["executions"]:
                execution["solver"]["elapsed_seconds"] = 0
        results[backend] = document
    compare(results["builtin"], results["fftw"])

    sweep = ROOT / "tests/dielectric/m3-frequency-sweep.nec"
    sweep_results = {}
    for backend in ("builtin", "fftw", "auto"):
        output = directory / f"sweep-{backend}.out"
        result = directory / f"sweep-{backend}.json"
        timing = directory / f"sweep-{backend}-timing.json"
        command = [str(PROGRAM), "-s", "fft"]
        if backend != "auto":
            command += ["--fft-backend", backend]
        command += ["--jobs", "4", "--threads", "2",
                    "--result-json", str(result), "--timing-json", str(timing),
                    "-i", str(sweep), "-o", str(output)]
        invoke(command)
        # This sweep has a one-cell padded grid, below the measured automatic
        # FFTW crossover.  Explicit FFTW remains available and authoritative.
        used = "builtin" if backend == "auto" else backend
        assert f"FFT BACKEND: {used}" in output.read_text()
        document = json.loads(result.read_text())
        for job in document["jobs"]:
            for execution in job["executions"]:
                execution["solver"]["elapsed_seconds"] = 0
        concurrency = json.loads(timing.read_text())["concurrency"]
        backend_timing = json.loads(timing.read_text())["fft_backend"]
        assert backend_timing == {"requested": backend, "used": used}, (
            backend, backend_timing
        )
        assert concurrency["jobs_used"] == 4
        assert concurrency["threads_per_job"] == 2
        sweep_results[backend] = document
    compare(sweep_results["builtin"], sweep_results["fftw"])
    compare(sweep_results["builtin"], sweep_results["auto"])

print("builtin/FFTW dense-oracle, full-solve, and CLI equivalence passed")
