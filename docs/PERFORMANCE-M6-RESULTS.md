# Milestone 6: parallel classic NEC matrix fill

Status: implementation complete; awaiting user acceptance.

## Delivered behavior

Classic NEC matrix fill now uses process-isolated, disjoint column workers on
Linux.  Each worker executes the existing interaction kernels and retains the
legacy source order within every matrix entry.  Process isolation gives each
worker private copies of the historical mutable field, interpolation, and
connection scratch without rewriting the validated numerical kernels.

`--threads auto` remains the default and permits up to 16 fill workers when
the matrix has at least 256 reduced columns.  This automatic threshold avoids
fork overhead on small fills.  `--threads N` selects an explicit worker count
from 128 columns upward, while `--threads 1` selects the preserved serial
`cmset` path.
The existing `--jobs` planner divides the CPU budget between independent
frequency jobs and fill workers, preventing their product from exceeding the
available CPU count.  Automatic LAPACK factorization remains capped at four
threads; an explicit thread count continues to override that policy.

Workers write an anonymous shared mapping with exclusive column ownership.
The parent copies it into the legacy matrix only after every exact child PID
exits successfully.  A mapping, allocation, fork, wait, or worker diagnostic
failure discards the temporary matrix and reruns the original serial fill.
The planner accounts for the temporary matrix by charging twice the normal
matrix storage.  Unsupported platforms use the serial fallback.

This design follows the project's accepted Linux/POSIX assumption.  Forking
is allowed only before a threaded numerical runtime has been activated in the
current process.  Built-in factorization does not activate one, so repeated
built-in fills remain eligible.  A later fill after LAPACK activation in the
same process falls back to serial.  Ordinary multi-frequency runs already
isolate frequencies in pristine worker processes, so their per-job fill
worker counts still share the global CPU budget.  Linux parent-death signaling
also terminates nested fill workers if their owning NEC process is cancelled.

## Numerical coverage

The focused regression compares normalized final result JSON from one worker
and four workers and requires byte-identical hashes.  Its fixtures cover:

* free-space wire interactions;
* perfect and finite real ground;
* extended thin-wire kernels;
* wire/patch, patch/wire, and patch/patch interactions;
* a wire endpoint's secondary connected-patch correction;
* lumped loading and two-port networks; and
* fourfold rotational symmetry and symmetry-mode combination; and
* repeated fills in a three-frequency built-in serial-job sweep.

The original serial implementation remains in the source as the one-worker
path, numerical oracle, and automatic error fallback.  The AEG GPL source was
used only for architectural observation; no source or expression was copied.
The implemented process-isolation specification was developed independently.

The full regression suite passes 33 of 33 groups, and the new focused group
passes all interaction families plus repeated built-in fills.  No existing
numerical baseline changed.  The focused group also passes a combined
AddressSanitizer and UndefinedBehaviorSanitizer build at one and four workers;
leak detection is disabled because it cannot initialize under the host's
tracing mode.  A complete distribution build, test, install, uninstall, and
archive check passes.

## Local performance evidence

The development host is the 32-logical-CPU Intel Core i9-13900K acceptance
machine.  The table contains one warmup and five measured repetitions.  Times
are medians of the `matrix_fill` phase; numerical result hashes are identical
at every worker count.

| Workers | Wire 801 time | Speedup | Mixed 651 time | Speedup |
|---:|---:|---:|---:|---:|
| 1 | 0.086005 s | 1.00x | 0.040952 s | 1.00x |
| 2 | 0.061290 s | 1.40x | 0.032284 s | 1.27x |
| 4 | 0.042538 s | 2.02x | 0.024426 s | 1.68x |
| 8 | 0.026606 s | 3.23x | 0.013743 s | 2.98x |
| 16 | 0.019111 s | 4.50x | 0.008746 s | 4.68x |
| 32 | 0.014155 s | 6.08x | 0.005661 s | 7.23x |

The automatic 16-worker point exceeds the milestone's 2x fill target on both
fixtures.  The extra 32-worker gain remains available explicitly but increases
process CPU substantially and slightly regresses whole-process wall time on
these short fixtures.  Automatic selection therefore retains the measured
16-worker ceiling.  The reported maximum single-process RSS rose from 18 MiB
serial to approximately 26 MiB for the wire fixture and 19 MiB for the mixed
fixture.  This metric is not aggregate process-tree RSS; the planner separately
accounts for the parent matrix plus one shared temporary matrix.

## User acceptance

First run the deterministic focused regression:

```sh
python3 tests/classic_fill_parallel.py
```

Expected output:

```text
classic matrix fill is exact across interaction families and repeated builtin fills
```

Then run the complete scaling study:

```sh
python3 tools/m6_fill_benchmark.py \
  --threads 1,2,4,8,16,32 \
  --warmup 2 --repeat 7 \
  --output /tmp/m6-fill.json
```

The benchmark fails unless every repeat and worker count has the same
normalized numerical hash and both fixtures achieve at least a 2x fill
speedup.  The JSON file retains all samples, wall and CPU times, peak RSS,
affinity, environment policy, and validation decisions.  Do not add an
`OMP_PROC_BIND` setting for this host unless testing an affinity policy
explicitly.

For an application-level check, run an exported Antenna Turtle model twice:

```sh
./nec2c --jobs 1 --threads 1 --backend auto \
  --result-json /tmp/m6-model-serial.json \
  --timing-json /tmp/m6-model-serial-timing.json \
  -i MODEL.nec -o /tmp/m6-model-serial.out

./nec2c --jobs 1 --threads auto --backend auto \
  --result-json /tmp/m6-model-auto.json \
  --timing-json /tmp/m6-model-auto-timing.json \
  -i MODEL.nec -o /tmp/m6-model-auto.out
```

The normalized numerical results should agree.  On this host the automatic
timing document should report 16 threads per job for a single sufficiently
large model.  No additional software or system library is required for this
milestone.
