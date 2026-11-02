#!/usr/bin/env python3
"""Run repeatable nec2c performance cases and emit one machine-readable report."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    temporary.replace(path)


def scrub_elapsed(value):
    """Remove timing-only result fields before calculating an accuracy hash."""
    if isinstance(value, dict):
        return {key: scrub_elapsed(item) for key, item in value.items()
                if key not in {"elapsed_seconds"}}
    if isinstance(value, list):
        return [scrub_elapsed(item) for item in value]
    return value


def numerical_hash(path):
    with path.open(encoding="utf-8") as stream:
        value = scrub_elapsed(json.load(stream))
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"),
                           allow_nan=False).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def git_value(*arguments):
    completed = subprocess.run(["git", *arguments], cwd=ROOT, text=True,
                               capture_output=True)
    return completed.stdout.strip() if completed.returncode == 0 else None


def command_output(command):
    completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    return completed.stdout.strip() if completed.returncode == 0 else None


def make_variable(name):
    try:
        pattern = re.compile(rf"^{re.escape(name)}\s*=\s*(.*)$")
        for line in (ROOT / "Makefile").read_text(encoding="utf-8").splitlines():
            match = pattern.match(line)
            if match:
                return match.group(1)
    except OSError:
        pass
    return None


def sysfs_value(path):
    try:
        return Path(path).read_text(encoding="ascii").strip()
    except OSError:
        return None


def cpu_model():
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or None


def summary(values):
    if not values:
        return None
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "maximum": max(values),
    }


def process_tree_rss_kib(root_pid):
    """Return aggregate resident memory for a Linux process tree."""
    total = 0
    pending = [root_pid]
    seen = set()
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        try:
            status = Path(f"/proc/{pid}/status").read_text(encoding="ascii")
            match = re.search(r"^VmRSS:\s+(\d+)\s+kB$", status, re.MULTILINE)
            if match:
                total += int(match.group(1))
            children = Path(f"/proc/{pid}/task/{pid}/children").read_text(
                encoding="ascii")
            pending.extend(int(value) for value in children.split())
        except (OSError, ValueError):
            continue
    return total


def summarize_runs(runs):
    phases = {}
    for run in runs:
        for phase in run["timing"]["phases"]:
            phases.setdefault(phase["name"], []).append(phase["seconds"])
    return {
        "wall_seconds": summary([run["timing"]["process"]["wall_seconds"] for run in runs]),
        "cpu_seconds": summary([run["timing"]["process"]["cpu_seconds"] for run in runs]),
        "peak_rss_kib": summary([run["timing"]["process"]["peak_rss_kib"] for run in runs
                                 if run["timing"]["process"]["peak_rss_kib"] is not None]),
        "external_wall_seconds": summary([run["external_wall_seconds"] for run in runs]),
        "process_tree_peak_rss_kib": summary([
            run["process_tree_peak_rss_kib"] for run in runs
            if run["process_tree_peak_rss_kib"] is not None]),
        "phases": {name: summary(values) for name, values in sorted(phases.items())},
        "stable_numerical_result": len({run["numerical_sha256"] for run in runs}) == 1,
        "dense_backends": sorted({
            run["timing"]["dense_backend"]["used"] for run in runs
            if run["timing"]["dense_backend"]["used"] is not None
        }),
    }


def execute(binary, case, root, directory, run_index, warmup, backend, jobs,
            threads, environment):
    deck = (root / case["deck"]).resolve()
    prefix = "warmup" if warmup else "run"
    stem = f"{case['id']}-{prefix}-{run_index}"
    output = directory / f"{stem}.out"
    result = directory / f"{stem}-result.json"
    timing = directory / f"{stem}-timing.json"
    command = [str(binary), "--backend", backend, "--jobs", jobs,
               "--threads", threads, *case.get("arguments", []),
               "--result-json", str(result), "--timing-json", str(timing),
               "-i", str(deck), "-o", str(output)]
    started = time.perf_counter()
    process = subprocess.Popen(command, cwd=root, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               env=environment)
    peak_rss = 0
    while process.poll() is None:
        if platform.system() == "Linux":
            peak_rss = max(peak_rss, process_tree_rss_kib(process.pid))
        time.sleep(0.005)
    stdout, stderr = process.communicate()
    external = time.perf_counter() - started
    if process.returncode != 0:
        raise RuntimeError(json.dumps({
            "case": case["id"], "command": command,
            "returncode": process.returncode,
            "stdout": stdout, "stderr": stderr,
        }, indent=2))
    with timing.open(encoding="utf-8") as stream:
        timing_value = json.load(stream)
    if timing_value["completion"]["status"] != "ok":
        raise RuntimeError(f"case {case['id']} emitted an unsuccessful timing artifact")
    return {
        "index": run_index,
        "command": command,
        "external_wall_seconds": external,
        "process_tree_peak_rss_kib": peak_rss or None,
        "timing": timing_value,
        "numerical_sha256": numerical_hash(result),
        "result_sha256": sha256(result),
        "output_sha256": sha256(output),
    }


def main():
    def auto_or_positive(value):
        if value == "auto":
            return value
        try:
            parsed = int(value)
        except ValueError as error:
            raise argparse.ArgumentTypeError("requires auto or a positive integer") from error
        if parsed < 1:
            raise argparse.ArgumentTypeError("requires auto or a positive integer")
        return str(parsed)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "nec2c")
    parser.add_argument("--manifest", type=Path,
                        default=ROOT / "tests/performance/models.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--backend", choices=("auto", "builtin", "lapack"),
                        default="auto")
    parser.add_argument("--jobs", type=auto_or_positive, default="auto",
                        help="frequency workers: auto or a positive integer")
    parser.add_argument("--threads", type=auto_or_positive, default="auto",
                        help="numerical-library threads per worker: auto or a positive integer")
    parser.add_argument("--case", action="append", dest="cases",
                        help="run only this case ID; may be repeated")
    arguments = parser.parse_args()
    if arguments.repeat < 1 or arguments.warmup < 0:
        parser.error("--repeat must be positive; --warmup cannot be negative")

    child_environment = os.environ.copy()
    if arguments.threads != "auto":
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "BLIS_NUM_THREADS"):
            child_environment[name] = str(arguments.threads)

    binary = arguments.binary.resolve()
    manifest_path = arguments.manifest.resolve()
    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    selected = [case for case in manifest["cases"]
                if not arguments.cases or case["id"] in arguments.cases]
    missing = set(arguments.cases or []) - {case["id"] for case in selected}
    if missing:
        parser.error("unknown case ID(s): " + ", ".join(sorted(missing)))

    report = {
        "schema": "nec2c-performance-report",
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "repository": {
            "path": str(ROOT),
            "commit": git_value("rev-parse", "HEAD"),
            "dirty": bool(git_value("status", "--porcelain")),
        },
        "host": {
            "platform": platform.platform(),
            "cpu_model": cpu_model(),
            "logical_cpus": os.cpu_count(),
            "cpu_affinity": (sorted(os.sched_getaffinity(0))
                             if hasattr(os, "sched_getaffinity") else None),
            "frequency_governor": sysfs_value(
                "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
            "frequency_driver": sysfs_value(
                "/sys/devices/system/cpu/cpu0/cpufreq/scaling_driver"),
            "python": platform.python_version(),
        },
        "build": {
            "cc": make_variable("CC"),
            "cflags": make_variable("CFLAGS"),
            "compiler_version": command_output(
                [(make_variable("CC") or "cc").split()[0], "--version"]),
        },
        "binary": {
            "path": str(binary),
            "sha256": sha256(binary),
            "linked_libraries": command_output(["ldd", str(binary)]),
        },
        "manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
        "settings": {"repeat": arguments.repeat, "warmup": arguments.warmup,
                     "backend": arguments.backend, "jobs": arguments.jobs,
                     "threads": arguments.threads},
        "execution_control": {
            "backend": arguments.backend,
            "thread_environment": {
                name: child_environment.get(name) for name in
                ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS", "FFTW_NUM_THREADS")
            },
        },
        "cases": [],
    }
    with tempfile.TemporaryDirectory(prefix="nec2c-benchmark-") as temporary:
        directory = Path(temporary)
        for case in selected:
            deck = (ROOT / case["deck"]).resolve()
            warmup_runs = [execute(binary, case, ROOT, directory, index, True,
                                   arguments.backend, arguments.jobs,
                                   arguments.threads, child_environment)
                           for index in range(arguments.warmup)]
            runs = [execute(binary, case, ROOT, directory, index, False,
                            arguments.backend, arguments.jobs,
                            arguments.threads, child_environment)
                    for index in range(arguments.repeat)]
            report["cases"].append({
                "id": case["id"],
                "description": case["description"],
                "provenance": case["provenance"],
                "deck": str(deck),
                "deck_sha256": sha256(deck),
                "arguments": case.get("arguments", []),
                "warmup_runs": warmup_runs,
                "runs": runs,
                "summary": summarize_runs(runs),
            })
    atomic_json(arguments.output.resolve(), report)
    for case in report["cases"]:
        median = case["summary"]["wall_seconds"]["median"]
        stable = "stable" if case["summary"]["stable_numerical_result"] else "CHANGED"
        print(f"{case['id']}: median {median:.6f}s, numerical result {stable}")
    print(arguments.output.resolve())


if __name__ == "__main__":
    main()
