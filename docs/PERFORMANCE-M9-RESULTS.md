# Performance Milestone 9 results

Status: accepted by the user on 2026-11-02.

Milestone 9 turns the previously accepted performance work into a portable,
measurable release candidate.  It does not merge the performance branch or
replace the user acceptance gate.

## Retained changes

* The portable default compiler optimization is `-O3`.  User-supplied
  `CFLAGS` remain authoritative.  Fast-math and `-march=native` are excluded:
  the former weakens the numerical contract and the latter changed canonical
  hashes in testing as well as limiting binary portability.
* Built-in batched solves reuse an existing scratch vector.  The common
  single-symmetry LAPACK path no longer allocates unused scratch.  Neither
  change alters arithmetic order.
* Automatic worker selection avoids parallel startup on tiny coupled-volume
  systems and on classic radiation work below the existing parallel-work
  threshold.  Explicit `--threads N` remains authoritative.
* Automatic FFT selection uses the built-in implementation below a 512-point
  padded grid and FFTW at or above it.  Explicit `--fft-backend` remains
  authoritative.  CPU-pinned, interleaved measurements showed built-in ahead
  by 3.3% for one cell and 1-2% for eight cells, parity at 27 cells, and FFTW
  pulling ahead thereafter.
* Parallel frequency workers flush inherited JSON streams before `fork()`,
  preventing duplicated buffered fragments in result documents.
* `tools/performance_benchmark.py` can select the FFT backend, run the
  historical command-line interface, report external CPU time and process-tree
  RSS, and compare backend-neutral physics hashes.
* `tools/m9_release_candidate.py` collects the complete correctness,
  performance, numerical-oracle, scaling, crossover, and memory evidence in a
  single indexed JSON bundle.

## Developer validation

The release candidate completed all 36 Automake tests and `make distcheck` in
an out-of-tree `-O3` build.  Controlled `-O2`/`-O3` tests retained exact
canonical hashes; their interleaved medians were 1.01x for the classic
801-frequency sweep and 1.21x for the dielectric FFT workload.

The current built-in, single-worker corpus retained exact backend-neutral
physics hashes against an executable built from the immutable M0 commit
`1166bc3bc5ba490c56fadf946070998a47ff8e00`.  Cross-solver impedance checks
also passed their frozen limits against AEG and momwire.  AEG was used only as
an external numerical oracle; no GPL source or object code was copied into the
project.

Small-model crossover tests confirmed that automatic execution suppresses
unproductive worker and FFTW startup.  The release gate compares interquartile
ranges rather than full sample ranges, preventing a lone scheduler outlier
from excusing a median regression.  Total command startup relative to M0 remains
separately visible because the current executable dynamically loads OpenMP,
LAPACK, and FFTW; it is not substituted for the solver crossover gate.

No additional software is required on the validated host.  OpenMP, LAPACKE,
and FFTW discovered by the existing build are sufficient.

## Acceptance

The complete reproducible command, required gates, and four manual scenarios
are in `docs/PERFORMANCE-M9-ACCEPTANCE.md`.  The top-level
`acceptance_ready` value in its generated JSON bundle must be true before the
milestone can be accepted.

The full user-run release bundle reported `acceptance_ready: true`, with all
commands passing, complete evidence present, stable automatic and legacy
results, and exact built-in single-worker physics hashes against M0.  A final
CPU-pinned 15-run check after tuning the small-grid FFT crossover measured:

| Case | Automatic median | Built-in single-worker median |
| --- | ---: | ---: |
| dielectric dense | 0.007283 s | 0.007287 s |
| dielectric FFT | 0.007355 s | 0.007330 s |

Both numerical results were stable.  The dense difference was effectively
zero, and the FFT automatic path was approximately 0.34% slower, comfortably
inside the 5% acceptance limit and normal measurement noise.  The user
returned PASS after reviewing these results.
