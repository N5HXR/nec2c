# Milestone 3: automatic frequency concurrency

Status: ready for user acceptance on 2026-11-02.

## Delivered behavior

`--jobs auto` and `--threads auto` are the defaults.  The explicit compatibility
path is `--jobs 1 --threads 1`; adding `--backend builtin` selects the complete
classic numerical path.

On Linux, independent frequencies execute in forked worker processes.  Each
worker owns a copy-on-write solver state and private report and result-JSON
fragments.  The parent merges those fragments strictly in frequency order.
For a classic NEC sweep the parent computes the final frequency itself, leaving
the matrix, currents, geometry, ground, and network state in the same final
state expected by subsequent cards.  Dielectric-only and coupled
conductor/dielectric sweeps use the same ordered worker contract.

Automatic allocation uses the process CPU-affinity mask, frequency count, a
measured four-thread ceiling for LAPACK, and an available-memory estimate that
reserves the greater of 20 percent or 1 GiB.  Non-threaded backends receive one
thread automatically.  A one-time early re-exec establishes a one-thread BLAS
environment before library constructors run; worker-local runtime controls then
apply the selected thread count.  This prevents hidden OpenBLAS/OpenMP
oversubscription while preserving authoritative CLI settings.

Timing JSON now records requested and selected concurrency, allowed CPUs, and
whether memory constrained the plan.  The benchmark runner also samples
aggregate Linux process-tree RSS rather than reporting only one process.

## Performance

Host: Intel Core i9-13900K, 32 allowed logical CPUs, Linux 6.8, GCC 12.3,
LAPACKE/OpenBLAS from Ubuntu.  Values are warm medians.

| Model | Repetitions | Serial | Automatic | Speedup | Median aggregate RSS |
|---|---:|---:|---:|---:|---:|
| Four-frequency 801-segment classic dipole | 7 | 0.516 s | 0.247 s | 2.09x | 19.2 -> 56.3 MiB |
| Five-frequency coupled wire/dielectric sweep | 3 | 13.929 s | 2.863 s | 4.87x | 9.2 -> 36.5 MiB |

Both benchmark configurations produced stable numerical results across all
repetitions.  Automatic memory growth is intentional: independent frequencies
own independent dense matrices.  It is bounded by the allocator.

## Numerical validation

Forced serial LAPACK (`--jobs 1 --threads 1`) was compared with the automatic
LAPACK result at every frequency:

| Model | Quantity | Maximum relative complex error |
|---|---|---:|
| Classic 801-segment sweep | feed impedance | 2.31e-11 |
| Classic 801-segment sweep | wire current | 5.24e-11 |
| Coupled dielectric sweep | feed impedance | 1.06e-14 |
| Coupled dielectric sweep | wire current | 2.60e-14 |
| Coupled dielectric sweep | volume current | 2.17e-14 |

The regression fixture additionally requires byte-identical ordered report
content after timing fields are normalized and structurally identical result
JSON after solver elapsed time is removed.  It covers classic and dielectric
frequency sweeps.  Linear and multiplicative stepping, impedance summaries,
finite real ground, and a radiation-pattern command followed by later-card
reuse were also compared against forced serial execution.

The unchanged cross-solver harness passed against external AEG and momwire
executables for the dipole and Yagi fixtures.  AEG remained an external oracle;
no AEG source was copied, and its earlier inspection was isolated as a
clean-room behavioral specification rather than implementation material.

## Automated verification

* 29 of 29 normal tests pass.
* `make distcheck` passes.
* The focused classic/dielectric parallel regression passes under combined
  AddressSanitizer and UndefinedBehaviorSanitizer.
* A build configured without LAPACKE passes the frequency-parallel regression.
* A four-CPU `taskset` restriction selects four jobs and one thread per job,
  confirming that automatic allocation honors affinity.

No new system package is required for this milestone.

## Conservative serial fallbacks

Correctness takes precedence over concurrency.  Classic execution remains
serial for matrix-reuse entry paths, a later sweep after the parent has already
entered a numerical runtime, receiving-pattern plane-wave excitation, CP
coupling, plot-file output, and negative-conductivity ground notation.  A
frequency-dependent real-ground sweep is parallel for the normal text report,
but currently stays serial when result JSON is requested because schema-v1
ground-state fragments need a separate ordered merge.  Unsupported platforms
use the serial path; Linux/POSIX is the accepted Milestone 3 target.

## Manual acceptance

Build once, then run the automated serial/parallel equivalence fixture:

```sh
make -j8
python3 tests/frequency_parallel.py
```

Benchmark the classic sweep using the new default and the compatibility path:

```sh
python3 tools/performance_benchmark.py \
  --case classic-frequency-sweep-801 --backend lapack \
  --jobs 1 --threads 1 --repeat 7 --warmup 2 \
  --output /tmp/m3-classic-serial.json

python3 tools/performance_benchmark.py \
  --case classic-frequency-sweep-801 --backend lapack \
  --jobs auto --threads auto --repeat 7 --warmup 2 \
  --output /tmp/m3-classic-auto.json
```

For a heavier end-to-end check, substitute
`--case dielectric-frequency-sweep --repeat 3 --warmup 1`.  The input models
are `tests/performance/scaling-dipole-801-sweep.nec` and
`tests/dielectric/m6-n2.nec`.  Their ordinary `.out` reports should retain
ascending frequency order.
