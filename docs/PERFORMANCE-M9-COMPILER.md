# Milestone 9 compiler optimization evaluation

NEC2C now uses `-Wall -O3` when the caller does not supply `CFLAGS`.  An
explicit `CFLAGS` value remains authoritative, including an intentionally empty
value.  This is a portable optimization level rather than a CPU-specific
binary contract: the default deliberately excludes `-march=native`,
`-mtune=native`, and all fast-math options.

## Evidence

Three isolated builds from the same commit were compared on the acceptance
host with GCC 12.3.0:

1. `-Wall -O2` (the previous default);
2. `-Wall -O3` (the new default); and
3. `-Wall -O3 -march=native -mtune=native` (evaluation only).

The portable `-O3` build produced the exact same canonical numerical hashes as
`-O2` for the 801-segment classic solve, its frequency sweep, and the
dielectric FFT solve.  It materially improved the frequency-sweep and FFT
medians in the exploratory run, although short phase timings were noisy enough
that the compiler change should not claim a universal fixed speedup.  The full
test suite is the release gate, not those exploratory timings.

| Case (18 interleaved measured runs) | `-O2` median | `-O3` median | Ratio |
| --- | ---: | ---: | ---: |
| Classic 801 frequency sweep | 0.492144 s | 0.488509 s | 1.01x |
| Dielectric FFT | 0.010108 s | 0.008364 s | 1.21x |

Each compiler variant ran once per round, with the order reversed on alternating
rounds and one numerical thread used to reduce heterogeneous-core scheduling
noise.  Short multithreaded and fill-only probes moved in both directions by
workload, so `-O3` is not assigned a universal fixed speedup.  It is adopted
because it exposes portable optimizer passes, preserves the tested numerical
artifacts, and has meaningful wins where optimization dominates startup.

The native build changed all three canonical hashes.  This is consistent with
the compiler using target-specific fused operations on the test host, and its
instructions would not run on every supported x86-64 CPU.  Its incremental
performance was also inconsistent.  NEC2C therefore does not enable native
targeting or runtime CPU dispatch in this release.  Users building for a known
machine may still opt in through `CFLAGS`, accepting the portability and
last-bit reproducibility tradeoffs.

Fast-math was not benchmarked as a candidate.  It permits transformations that
violate the IEEE behavior assumed by complex arithmetic and convergence/error
checks, so it is outside the accuracy-preserving performance scope.

## Reproduction

Regenerate the build system, then confirm the default and override paths:

```sh
autoreconf -fi
./configure --with-lapack=auto --with-fftw=auto
make -j8
make check

make distclean
./configure CFLAGS='-Wall -O2' --with-lapack=auto --with-fftw=auto
make -j8
make check
```

The first configure run must report compilation commands containing `-O3`.
The second must contain `-O2` and no injected `-O3`.  Release validation also
runs `make distcheck` from the default build.
