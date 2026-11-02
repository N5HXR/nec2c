# Milestone 4B: optional threaded FFTW backend

Status: ready for user acceptance on 2026-11-02.

## Delivered behavior

`--fft-backend auto|builtin|fftw` selects the implementation of the existing
dielectric FFT convolution.  FFTW is detected through its header, core
library, and threads library; this host provides FFTW 3.3.8.  `auto` uses FFTW
when it is compiled in and otherwise falls back to the accepted parallel
built-in transform.  Explicit `fftw` fails clearly when unavailable, while
explicit `builtin` preserves the Milestone 4A implementation.

Each FFTW plan captures the Milestone 3 per-worker thread budget.  Three
in-place batched plans are retained for the six-component kernel transform and
the three-component forward and inverse applications.  Planning uses
`FFTW_ESTIMATE`; the existing power-of-two padding, grid layout, signed kernel,
and scatter/multiply/gather operations are unchanged.  Forward transforms are
unscaled and the inverse is explicitly divided by the padded grid volume,
matching the built-in convention.

FFTW initialization is lazy.  All plans are destroyed before its process-local
thread resources are released, so the preflight plan is gone before Milestone
3 forks frequency workers.  No FFTW activity occurs during option parsing or
capability reporting.

## Performance

Host: Intel Core i9-13900K, 32 allowed logical CPUs, Linux 6.8, GCC 12.3,
FFTW 3.3.8.  Values are warm medians on the 24 x 24 x 48 padded production
grid with 9,432 occupied cells and 28,296 unknowns.

| Threads | Built-in matvec | FFTW matvec | FFTW advantage |
|---:|---:|---:|---:|
| 1 | 0.127503 s | 0.031706 s | 4.02x |
| 2 | 0.069749 s | 0.025572 s | 2.73x |
| 4 | 0.041214 s | 0.021053 s | 1.96x |
| 8 | 0.038951 s | 0.019160 s | 2.03x |
| 16 | 0.029033 s | 0.017500 s | 1.66x |
| 32 | 0.028796 s | 0.036245 s | 0.79x |

The automatic thread cap remains 16; FFTW's 32-thread regression therefore
does not affect default execution.  Interleaved backend runs avoid
thermal/order bias.  On complete 26-iteration production solves, FFTW improved
the 16-thread median from 0.574810 to 0.374051 seconds, a 1.54x solve speedup
over M4A; the matvec improved 1.46x.  At one thread it improved the solve from
3.452000 to 0.755782 seconds, or 4.57x.  Including the common kernel-build
phase, complete 16-thread invocation time improved 1.17x.  Peak RSS increased
from 124,928 to 127,232 KiB.

On the 216-cell skeleton at 16 threads, interleaved medians show effectively
equal cold matvec latency (0.002398 versus 0.002406 seconds, a 0.33-percent
difference).  The complete solve improved 1.77x and build-plus-solve latency
improved 1.09x.  This supports using FFTW by default without a size crossover
on the measured workload set.

## Numerical and operational validation

The independent dense operator remains the primary oracle.  Both backends pass
the one-cell, zero-frequency, lossy, sparse, anisotropic-grid, far-cell,
PEC-image, heterogeneous-material, and coupled-system checks.  Full GMRES
residuals, integrated complex currents, far-axis fields, extinction,
scattering, and power defect agree within `2e-9` relative and `2e-11`
absolute.  CLI result documents agree at the same tolerance.

Within each backend and thread count the production result is repeatable.
Different FFT implementations are compared by tolerance rather than byte hash
because their butterfly order is allowed to differ.  A four-frequency run at
four jobs x two threads produces ordered, equivalent builtin/FFTW results and
stays within the declared CPU budget.

The build matrix also verifies:

* `--with-fftw=no` compiles and passes the built-in parallel tests;
* explicit `--fft-backend fftw` is rejected when FFTW is absent; and
* forced `--with-fftw=yes` fails configure when the dependency probe fails;
* an FFTW build without compiler OpenMP gives explicit FFTW its requested
  thread budget while keeping explicit built-in execution and scheduling
  serial; and
* the pthread dependency is discovered explicitly for static/non-glibc link
  environments as well as the supported Linux shared-library build.

All 31 normal regression groups pass.  The distribution archive rebuilds
out-of-tree, repeats all 31 tests, installs, uninstalls, and passes
`make distcheck`.  A four-thread production FFTW solve passes combined
AddressSanitizer and UndefinedBehaviorSanitizer instrumentation.  LeakSanitizer
cannot initialize under this host's tracing mode, so leak detection was
disabled for that run.  ThreadSanitizer has the same `unexpected memory
mapping` host limitation recorded in Milestone 4A; thread safety is covered by
the non-reentrant per-plan ownership, FFTW's documented execution contract,
forked-frequency tests, and repeated numerical results.

FFTW exposes no plan-memory byte count.  NEC2C conservatively retains the
built-in root-table allowance in its 512 MiB account and separately measures
whole-process peak RSS.  FFTW's measured production RSS increment is under
2 MiB at the default thread count.  Library-private metadata is the only
portion not directly itemized by `FFT ACCOUNTED APPLICATION BYTES`.  The
public API and report deliberately call this application-owned accounting
rather than claiming an opaque-library hard cap.

Normal timing JSON records `fft_backend.requested` and `fft_backend.used`, so
saved measurements establish whether `auto` used FFTW or the built-in
fallback.  A backend change within one run is represented as `mixed`.  This
addition is timing schema version 2; the version 1 schema remains unchanged
for validation of retained Milestone 1 artifacts.

The AEG clean-room comparison found no FFT implementation.  Momwire uses
SciPy's FFT rather than FFTW; its transferable behaviors were cached kernel
spectra, explicit fallback diagnostics, dense/matvec/end-to-end validation,
and a normalized forward/multiply/inverse contract.  No source was copied.
FFTW is GPL-licensed and is compatible with this GPLv3 program as an optional
linked backend.

## Manual acceptance

Regenerate the build system so the optional dependency probe is current:

```sh
autoreconf -fi
./configure --with-lapack=auto --with-fftw=auto
make -j8
```

Configure must report `optional threaded FFTW backend: yes`.  Run the backend,
dense-oracle, full-solve, CLI, forked-frequency, and M4A fallback tests:

```sh
python3 tests/fft_backend.py
python3 tests/fft_parallel.py --require-parallel
```

Reproduce the production full-solve comparison:

```sh
python3 tools/fft_benchmark.py \
  --model tree --backend builtin --threads 1,16 \
  --max-iterations 300 --warmup 1 --repeat 3 \
  --output /tmp/m4b-builtin.json

python3 tools/fft_benchmark.py \
  --model tree --backend fftw --threads 1,16 \
  --max-iterations 300 --warmup 1 --repeat 3 \
  --output /tmp/m4b-fftw.json

python3 tools/fft_backend_compare.py \
  --model tree --threads 1,16 --max-iterations 300 \
  --warmup 1 --repeat 5 --output /tmp/m4b-interleaved.json
```

The report generated by a normal `-s fft` run must say `FFT BACKEND: fftw`.
Adding `--fft-backend builtin` must instead report `FFT BACKEND: builtin`.
The committed scaling, full-solve, interleaved, and small-workload records are under
`tests/performance/m4b-*`.
