# Milestone 4A: parallel built-in FFT

Status: ready for user acceptance on 2026-11-02.

## Delivered behavior

The built-in radix-2 FFT now batches tensor components and executes independent
one-dimensional lines with OpenMP.  A barrier remains between the z, y, and x
phases, and the arithmetic inside each line is unchanged.  Source scatter,
spectral tensor products, destination gather, preconditioning, and disjoint
conductor/volume rows also use the same fixed thread budget when large enough.

Each plan captures an immutable thread count when it is created.  `--threads 1`
therefore selects the exact serial arithmetic path.  Automatic execution uses
up to 16 threads for a single FFT solve, the measured optimum for a full GMRES
solve on the acceptance host.  The Milestone 3 allocator divides CPUs between
frequency jobs first: under an eight-CPU affinity mask, the four-frequency
fixture selects four jobs with two threads each.

OpenMP is detected at configure time.  A compiler without OpenMP builds and
passes the same tests using the serial fallback.  The ordinary report includes
the FFT plan's captured thread count.  No new system package is required.

## Performance

Host: Intel Core i9-13900K, 32 allowed logical CPUs, Linux 6.8, GCC 12.3.
Values are warm medians on the 24 x 24 x 48 padded production grid with 9,432
occupied cells and 28,296 unknowns.

| Threads | Matvec | Matvec speedup | One-iteration solve phase |
|---:|---:|---:|---:|
| 1 | 0.127471 s | 1.00x | 0.362531 s |
| 2 | 0.070241 s | 1.81x | 0.188641 s |
| 4 | 0.041126 s | 3.10x | 0.102662 s |
| 8 | 0.038473 s | 3.31x | 0.076352 s |
| 16 | 0.026281 s | 4.85x | 0.056639 s |
| 32 | 0.031543 s | 4.04x | 0.063462 s |

The full 26-iteration GMRES solve improved from a 3.429238-second median at one
thread to 0.677552 seconds at 16 threads, a 5.06x solve-phase speedup.  Median
peak RSS changed from 126,464 KiB to 124,928 KiB.  The 216-cell skeleton also
improved rather than regressed: its matvec changed from 0.010814 seconds at one
thread to 0.002436 seconds at 16 threads.  The explicit serial branch is within
0.3 percent of the pre-M4 production-grid baseline (0.127094 seconds), well
inside the approved five-percent small-workload limit.

## Numerical validation

The serial and parallel implementations produce bit-identical hashes over
every operator-product byte, as well as identical residuals, iteration counts,
integrated currents, far-axis fields,
extinction, scattering, and power defect at 1, 2, 4, and 16 threads.  This is
expected because scheduling changes only which disjoint line runs first; it
does not change any scalar summation or butterfly order.

The existing dense-oracle FFT, PEC-image, heterogeneous material, coupled-wire,
frequency-parallel, and result-schema tests remain unchanged and pass.  The new
production-grid regression executes the serial, two-thread, and four-thread
paths and requires exact hexadecimal product checksums.

## Verification

* 30 of 30 normal tests pass.
* The production-grid benchmark is repeatable across seven measured runs at
  1, 2, 4, 8, 16, and 32 threads.
* A configure run with OpenMP deliberately unavailable builds successfully,
  passes the serial-fallback FFT regression, and reports one automatic FFT
  thread through the real CLI timing contract.
* Combined AddressSanitizer and UndefinedBehaviorSanitizer execution of the
  four-thread production grid reports no memory or undefined-behavior error.
  LeakSanitizer is unavailable under the execution sandbox's tracing mode.
* GCC ThreadSanitizer cannot initialize on this host (`unexpected memory
  mapping`) with either PIE or non-PIE executables.  Race coverage therefore
  uses disjoint-output code inspection plus exact repeated results at every
  measured thread count; this environment limitation is not reported as a
  successful TSAN run.
* A four-frequency FFT sweep restricted to eight CPUs selects four jobs x two
  threads and completes successfully.

The AEG clean-room review found no FFT implementation to port.  Momwire's
SciPy behavior reinforced batched transforms, explicit fallback reporting, and
dense-equivalence validation; no source was copied from either project.

## Manual acceptance

Generated Autotools files are intentionally not versioned.  After checking out
or pulling this milestone, regenerate and reconfigure before building; a plain
`make` against an older Makefile will retain its old compiler flags.

Build and run the exact serial/parallel regression:

```sh
autoreconf -fi
./configure --with-lapack=auto
make -j8
python3 tests/fft_parallel.py --require-parallel
```

The configure output must include
`parallel built-in FFT: enabled with -fopenmp`.  The acceptance test now fails
with a rebuild instruction if it detects a serial-only benchmark binary.

Reproduce the scaling evidence:

```sh
python3 tools/fft_benchmark.py \
  --model tree --threads 1,2,4,8,16,32 \
  --warmup 2 --repeat 7 \
  --output /tmp/m4a-fft-scaling.json
```

Compare complete GMRES solves:

```sh
python3 tools/fft_benchmark.py \
  --model tree --threads 1,16 --max-iterations 300 \
  --warmup 1 --repeat 3 \
  --output /tmp/m4a-fft-full-solve.json
```

The committed evidence is in `tests/performance/m4a-fft-scaling.json` and
`tests/performance/m4a-fft-full-solve.json`.  Milestone 4B will add the optional
FFTW backend only after this built-in implementation is accepted, then apply
the approved 20-percent automatic-selection margin.
