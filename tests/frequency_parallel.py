#!/usr/bin/env python3
import json
import os
import pathlib
import re
import subprocess
import tempfile

ROOT = pathlib.Path(os.environ.get("TEST_SRCDIR",
    pathlib.Path(__file__).resolve().parents[1])).resolve()
PROGRAM = pathlib.Path(os.environ.get("NEC2C", ROOT / "nec2c")).resolve()
DECKS = {
    "dielectric": ROOT / "tests/dielectric/m3-frequency-sweep.nec",
    "classic": ROOT / "tests/dielectric/m6-n2-legacy-sweep.nec",
}


def run(directory: pathlib.Path, label: str, deck: pathlib.Path, jobs: str):
    report = directory / f"{label}-jobs-{jobs}.out"
    result = directory / f"{label}-jobs-{jobs}.json"
    subprocess.run(
        [str(PROGRAM), "--jobs", jobs, "--threads", "1", "--backend", "builtin",
         "--fft-backend", "auto",
         "-i", str(deck), "-o", str(report), "--result-json", str(result)],
        cwd=ROOT, check=True,
    )
    return report.read_text(), json.loads(result.read_text())


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
    expected = {
        "dielectric": [90, 100, 110, 120],
        "classic": [240, 270, 300, 330, 360],
    }
    for label, deck in DECKS.items():
        serial_report, serial_result = run(directory, label, deck, "1")
        parallel_report, parallel_result = run(directory, label, deck, "2")
        assert normalize_report(serial_report) == normalize_report(parallel_report)
        assert normalize_result(serial_result) == normalize_result(parallel_result)
        executions = parallel_result["jobs"][0]["executions"]
        assert [row["frequency_mhz"] for row in executions] == expected[label]
