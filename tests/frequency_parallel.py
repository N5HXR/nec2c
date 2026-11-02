#!/usr/bin/env python3
import json
import pathlib
import re
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROGRAM = ROOT / "nec2c"
DECK = ROOT / "tests/dielectric/m3-frequency-sweep.nec"


def run(directory: pathlib.Path, jobs: str):
    report = directory / f"jobs-{jobs}.out"
    result = directory / f"jobs-{jobs}.json"
    subprocess.run(
        [str(PROGRAM), "--jobs", jobs, "--threads", "1", "--backend", "builtin",
         "-i", str(DECK), "-o", str(report), "--result-json", str(result)],
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
    serial_report, serial_result = run(directory, "1")
    parallel_report, parallel_result = run(directory, "2")
    assert normalize_report(serial_report) == normalize_report(parallel_report)
    assert normalize_result(serial_result) == normalize_result(parallel_result)
    executions = parallel_result["jobs"][0]["executions"]
    assert [row["frequency_mhz"] for row in executions] == [90, 100, 110, 120]
