# Milestone 5: parallel dielectric assembly

Status: implemented; awaiting user acceptance.

## Delivered behavior

OpenMP now assembles independent dielectric matrix blocks concurrently while
preserving each individual entry's arithmetic order.  The parallel regions
cover dense volume-volume cell pairs and the volume-source blocks tested by
wire and patch rows.  The legacy wire-source field path remains serial because
it mutates NEC current and ground state; making that path reentrant remains
Milestone 7 work.

Each unordered volume pair owns both directed 3 by 3 blocks.  Each flattened
wire/cell or patch/cell task owns a distinct output block.  Workers publish
failures into private indexed status elements, and the serial code reports the
lowest canonical failing pair after the OpenMP barrier.  This avoids races in
diagnostics and makes error selection deterministic.

`--threads auto` is still the default.  OpenMP builds now allow up to 16
threads for assembly-bound builtin and dense-GMRES dielectric solves.  The
accepted four-thread LAPACK cap and 16-thread FFT cap remain unchanged.
`--threads 1` retains serial execution.  Builds without OpenMP retain the
serial loops and automatically use one thread for these paths.

Timing JSON exposes these phases:

* `dielectric_volume_assembly`
* `dielectric_volume_to_wire_assembly`
* `dielectric_volume_to_patch_assembly`

## Performance

The development host is an Intel Core i9-13900K with 32 logical CPUs.  Warm
medians below use one warmup and three measured repetitions.  The 216-cell
volume fixture has 648 volume unknowns.  The coupled fixture has the same 648
volume unknowns plus 11 wire and 8 patch unknowns.

| Threads | Volume assembly | Speedup | Coupled row assembly | Speedup |
|---:|---:|---:|---:|---:|
| 1 | 18.021853 s | 1.00x | 0.089256 s | 1.00x |
| 2 | 9.020507 s | 2.00x | 0.046094 s | 1.94x |
| 4 | 4.509922 s | 4.00x | 0.032448 s | 2.75x |
| 8 | 2.291416 s | 7.86x | 0.018345 s | 4.87x |
| 16 | 1.512987 s | 11.91x | 0.008485 s | 10.52x |
| 32 | 1.067452 s | 16.88x | 0.006313 s | 14.14x |

The automatic ceiling is 16 rather than 32 because the extra logical-thread
gain is smaller and leaves no CPU headroom.  Explicit `--threads 32` remains
available.

Do not force `OMP_PROC_BIND=close` for the acceptance run.  On this asymmetric
P-core/E-core host, the local OpenMP runtime's `cores` placement reduced the
effective run to roughly two cores.  NEC2C leaves affinity policy to the Linux
scheduler unless the caller explicitly overrides it.

## Numerical and operational validation

The benchmark requires bit-identical normalized result JSON across every
thread count and repeat.  Both workloads pass at 1, 2, 4, 8, 16, and 32
threads.  The focused regression independently verifies exact results at 1,
2, and 4 threads for volume-only and mixed wire/patch systems and verifies the
automatic CPU budget.

All 31 pre-existing regression groups pass.  The targeted dielectric suite
covers dense, GMRES, FFT, PEC, real ground, loss, wire, patch, and mixed
wire/patch cases.  A distribution build configured without OpenMP, LAPACK, or
FFTW compiles cleanly and passes the new serial-fallback regression.  A final
read-only race audit found no overlapping output ownership or parallel access
to the known mutable legacy field path.  The focused parallel regression also
passes combined AddressSanitizer and UndefinedBehaviorSanitizer instrumentation
at up to four threads; leak detection is disabled because LeakSanitizer cannot
initialize under this host's tracing mode.

The distribution archive rebuilds out of tree, runs all 32 regression groups,
installs, uninstalls, and passes `make distcheck`.

## User acceptance

Rebuild from the feature branch, then run:

```sh
python3 tests/dielectric_assembly_parallel.py

python3 tools/m5_assembly_benchmark.py \
  --threads 1,2,4,8,16,32 \
  --warmup 2 --repeat 7 \
  --output /tmp/m5-assembly.json
```

Expected regression output:

```text
dielectric assembly is bit-identical at 1, 2, and 4 threads
```

The benchmark passes only if both cases remain numerically identical and each
case reaches at least 2x assembly speedup.  It prints the median time and
speedup at every thread count.  `/tmp/m5-assembly.json` retains all raw samples,
CPU time, wall time, RSS, phase timings, numerical hashes, and validation
decisions.

For a manual default-policy check, compare the timing documents from these two
runs.  Both should produce the same numerical result; the automatic run should
report `threads_per_job` of 16 on this host.

```sh
./nec2c -s gmres --jobs 1 --threads 1 \
  --result-json /tmp/m5-serial-result.json \
  --timing-json /tmp/m5-serial-timing.json \
  -i tests/performance/m5-volume-216.nec -o /tmp/m5-serial.out

./nec2c -s gmres --jobs 1 --threads auto \
  --result-json /tmp/m5-auto-result.json \
  --timing-json /tmp/m5-auto-timing.json \
  -i tests/performance/m5-volume-216.nec -o /tmp/m5-auto.out
```
