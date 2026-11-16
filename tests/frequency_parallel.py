#!/usr/bin/env python3
import json
import os
import pathlib
import re
import resource
import subprocess
import tempfile

ROOT = pathlib.Path(os.environ.get("TEST_SRCDIR",
    pathlib.Path(__file__).resolve().parents[1])).resolve()
PROGRAM = pathlib.Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
DECKS = {
    "dielectric": ROOT / "tests/dielectric/m3-frequency-sweep.nec",
    "classic": ROOT / "tests/dielectric/m6-n2-legacy-sweep.nec",
    "real-ground": ROOT / "tests/performance/m3-real-ground-sweep.nec",
}


def run(directory: pathlib.Path, label: str, deck: pathlib.Path, jobs: str,
        fail_index=None, open_file_limit=None):
    report = directory / f"{label}-jobs-{jobs}.out"
    result = directory / f"{label}-jobs-{jobs}.json"
    timing = directory / f"{label}-jobs-{jobs}-timing.json"
    environment = os.environ.copy()
    if fail_index is not None:
        environment["NEC2C_TEST_CLASSIC_WORKER_FAIL_INDEX"] = str(fail_index)
    def set_open_file_limit():
        resource.setrlimit(resource.RLIMIT_NOFILE,
                           (open_file_limit, open_file_limit))
    completed = subprocess.run(
        [str(PROGRAM), "--jobs", jobs, "--threads", "1", "--backend", "builtin",
         "--fft-backend", "auto",
         "-i", str(deck), "-o", str(report), "--result-json", str(result),
         "--timing-json", str(timing)],
        cwd=ROOT, check=True, env=environment, text=True, capture_output=True,
        preexec_fn=set_open_file_limit if open_file_limit is not None else None,
    )
    return (report.read_text(), json.loads(result.read_text()),
            json.loads(timing.read_text()), completed.stderr)


def normalize_report(text: str) -> str:
    text = re.sub(r"FILL: \d+ msec  FACTOR: \d+ msec",
                  "FILL: TIME msec  FACTOR: TIME msec", text)
    return re.sub(r"TOTAL RUN TIME: \d+ msec", "TOTAL RUN TIME: TIME msec", text)


def normalize_result(document):
    for job in document["jobs"]:
        for execution in job["executions"]:
            solver = execution.get("solver")
            if solver:
                solver["elapsed_seconds"] = 0
    return document


with tempfile.TemporaryDirectory(prefix="nec2c-m3-") as temporary:
    directory = pathlib.Path(temporary)
    baselines = {}
    expected = {
        "dielectric": [90, 100, 110, 120],
        "classic": [240, 270, 300, 330, 360],
        "real-ground": [80, 85, 90, 95, 100, 105, 110, 115, 120],
    }
    for label, deck in DECKS.items():
        serial_report, serial_result, _, _ = run(directory, label, deck, "1")
        baselines[label] = (serial_report, serial_result)
        parallel_report, parallel_result, parallel_timing, _ = \
            run(directory, label, deck, "2")
        assert normalize_report(serial_report) == normalize_report(parallel_report)
        assert normalize_result(serial_result) == normalize_result(parallel_result)
        assert parallel_result["jobs"][0]["ground_states"] == \
            serial_result["jobs"][0]["ground_states"]
        assert parallel_timing["concurrency"]["jobs_used"] == 2
        executions = parallel_result["jobs"][0]["executions"]
        assert [row["frequency_mhz"] for row in executions] == expected[label]

    retry_report, retry_result, retry_timing, retry_stderr = run(
        directory, "real-ground-retry", DECKS["real-ground"], "2", 2)
    real_ground_report, real_ground_result = baselines["real-ground"]
    assert normalize_report(retry_report) == normalize_report(real_ground_report)
    assert normalize_result(retry_result) == normalize_result(real_ground_result)
    assert retry_timing["concurrency"]["jobs_used"] == 1
    assert "worker 2 at 85 MHz exited with status 124" in retry_stderr
    assert "retrying the frequency sweep serially" in retry_stderr

    long_deck = directory / "long-sweep.nec"
    long_deck.write_text(DECKS["classic"].read_text().replace(
        "FR 0 5 0 0 240.0 30.0", "FR 0 350 0 0 240.0 0.1"))
    _, long_result, long_timing, long_stderr = run(
        directory, "descriptor-pressure", long_deck, "2",
        open_file_limit=128)
    assert len(long_result["jobs"][0]["executions"]) == 350
    assert long_timing["concurrency"]["jobs_used"] == 2
    assert "Too many open files" not in long_stderr

print("classic, real-ground, and dielectric frequency sweeps preserve ordered results; worker failure retries serially and long sweeps keep bounded descriptors")
