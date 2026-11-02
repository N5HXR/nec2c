# Milestone 7: parallel conductor-volume coupling

Status: accepted by the user on 2026-11-02.

## Delivered behavior

The coupled wire/patch and dielectric solver now freezes every conductor
unknown into explicit NEC wire-basis and Cartesian patch-current coefficients
before evaluating conductor-to-volume interactions.  Linux workers install
only a process-private copy of that state and own disjoint source/cell result
blocks in an anonymous shared mapping.  The parent publishes the completed
coupling block only after all exact worker PIDs exit successfully.  Allocation,
fork, wait, or worker-evaluation failure discards the partial result and reruns
the preserved one-worker evaluation.

`--threads auto` remains the default and permits up to 16 coupling workers.
`--threads N` selects an explicit process count, including 32 on the acceptance
host, while `--jobs 1 --threads 1` selects the serial reference path.  Frequency
jobs and nested coupling workers continue to share one CPU budget. This
milestone follows the accepted Linux execution contract; other platforms are
not part of M7 acceptance.

Coupling workers finish before OpenMP, BLAS, or FFTW threading is activated.
In particular, FFT system and FFTW plan creation is now deliberately delayed
until after the fork-based coupling phase.  Automatic direct/LAPACK solves keep
their four-thread numerical ceiling; explicit thread counts still override it.
The existing FFT and matrix backends remain selectable from the command line.

The volume-pair phase was moved after coupling so its OpenMP runtime cannot be
initialized before a fork. Its arithmetic and output ownership are unchanged.
The new `dielectric_wire_to_volume_assembly` timer isolates the optimized
phase from FFT planning, reverse coupling, and iteration. The text report also
records the actual coupling worker count and labels a serial error fallback;
the acceptance harness rejects a requested parallel run that did not fork.

## Numerical and safety coverage

The focused regression compares normalized final result JSON at 1, 2, and 4
configured workers and requires byte-identical hashes for:

* mixed wire and disconnected conducting-patch geometry in free space;
* a wire endpoint connected to the center of a conducting patch;
* a coupled wire/dielectric model over perfect ground; and
* a coupled wire/dielectric model over Sommerfeld real ground.

The automatic FFT backend also produces an identical coupled result at one and
four workers after the FFT planning reorder. The complete project suite passes
34 of 34 test groups with
no numerical baseline changes. A separate AddressSanitizer and
UndefinedBehaviorSanitizer distribution build passes the focused test; leak
detection is disabled because it cannot initialize under the host's tracing
mode. A complete distribution build, test, install, uninstall, and archive
check also passes. Momwire and AEG do not implement the same
coupled dielectric formulation, so the frozen built-in serial result is the
applicable M7 oracle; the M1 classic-wire cross-solver thresholds remain
unchanged.

The implementation is an independent process-isolation design. No AEG source
or expression was copied. The original implementation remains reachable in
the immutable branch reference for measurement and source-level comparison.

## Development performance evidence

The following three-repeat development medians were measured on the 32-logical
CPU host with one warmup. Numerical hashes were identical at every worker
count. The benchmark used the production `--fft-backend auto` policy.

| Workers | Coupling time | Coupling speedup | Complete solve | Solve speedup |
|---:|---:|---:|---:|---:|
| 1 | 0.388438 s | 1.00x | 0.536102 s | 1.00x |
| 2 | 0.317821 s | 1.22x | 0.428256 s | 1.25x |
| 4 | 0.139147 s | 2.79x | 0.235641 s | 2.28x |
| 8 | 0.096837 s | 4.01x | 0.183544 s | 2.92x |
| 16 | 0.059728 s | 6.50x | 0.165712 s | 3.24x |
| 32 | 0.052666 s | 7.38x | 0.188263 s | 2.85x |

This exceeds both milestone gates: at least 2x for the isolated coupling phase
and at least 1.2x for the complete coupled solve. The automatic ceiling remains
16 because 32 improved the isolated phase but regressed the complete solve in
the production-backend study; explicit `--threads 32` remains available.

Peak RSS in the JSON is the maximum RSS of one process, not aggregate
process-tree memory. The coupling state is quadratic in conductor unknowns and
the temporary shared coupling block is proportional to conductor unknowns
times volume unknowns. Milestone 9 will revisit copies, memory accounting, and
the automatic crossover with the full release-candidate workload.

As with the accepted M6 fill path, a `--jobs 1` multi-frequency run becomes
conservatively ineligible for additional forks after its first threaded
numerical solve. Default multi-frequency execution uses pristine frequency-job
processes and retains coupling parallelism. This conservative serial-job
fallback avoids forking a process after BLAS, OpenMP, or FFTW activation.

## User acceptance

The user returned a PASS on 2026-11-02. The clean acceptance run passed the
four-regime deterministic regression and confirmed exact normalized results at
1, 2, and 4 workers with the automatic FFT path. Two warmups and seven measured
repetitions produced these medians:

| Workers | Coupling time | Coupling speedup | Complete solve | Solve speedup |
|---:|---:|---:|---:|---:|
| 1 | 0.390809 s | 1.00x | 0.538325 s | 1.00x |
| 2 | 0.301898 s | 1.29x | 0.418166 s | 1.29x |
| 4 | 0.172332 s | 2.27x | 0.265102 s | 2.03x |
| 8 | 0.094385 s | 4.14x | 0.179263 s | 3.00x |
| 16 | 0.060990 s | 6.41x | 0.169977 s | 3.17x |
| 32 | 0.046760 s | 8.36x | 0.179302 s | 3.00x |

The accepted 16-worker automatic point exceeds both milestone gates. The
32-worker point accelerates the isolated coupling phase further but is slightly
slower end to end, validating the 16-worker automatic ceiling. The raw record
is `/tmp/m7-coupling.json` on the acceptance host.

The accepted validation procedure is retained below.

Run the focused deterministic regression first:

```sh
python3 tests/coupling_parallel.py
```

Expected output:

```text
free-space, connected patch, PEC, and real-ground coupling is bit-identical at 1, 2, and 4 workers; automatic FFT is exact
```

Then run the seven-repeat scaling study:

```sh
python3 tools/m7_coupling_benchmark.py \
  --threads 1,2,4,8,16,32 \
  --warmup 2 --repeat 7 \
  --output /tmp/m7-coupling.json
```

The benchmark fails unless all repeats and worker counts have the same
normalized numerical hash, the 16-worker acceptance point reaches at least 2x
for coupling and 1.2x for the complete solve, and actual-worker telemetry shows
no silent fallback. The JSON retains every sample, the resolved FFT backend,
phase and process times, peak single-process RSS, CPU affinity, environment,
and the validation decision.

For an application-level comparison, replace `MODEL.nec` with a coupled
wire/dielectric model containing a `DF` card:

```sh
./nec2c --jobs 1 --threads 1 -s fft --fft-backend auto \
  --result-json /tmp/m7-model-serial.json \
  --timing-json /tmp/m7-model-serial-timing.json \
  -i MODEL.nec -o /tmp/m7-model-serial.out

./nec2c --jobs 1 --threads auto -s fft --fft-backend auto \
  --result-json /tmp/m7-model-auto.json \
  --timing-json /tmp/m7-model-auto-timing.json \
  -i MODEL.nec -o /tmp/m7-model-auto.out
```

The normalized results should agree, and the automatic timing document should
report up to 16 threads for a single sufficiently large model. No additional
software or system library is required for this milestone.
