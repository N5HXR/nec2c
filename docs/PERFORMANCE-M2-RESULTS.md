# Milestone 2: optional LAPACK dense backend

Status: implementation complete; awaiting user acceptance.

Date: 2026-11-02

Branch: `codex/performance`

Preserved implementation reference: `1166bc3bc5ba490c56fadf946070998a47ff8e00`

## Delivered behavior

Milestone 2 adds one dense factor/solve abstraction shared by classic NEC,
dense VIE, and dense coupled conductor-volume systems. The command-line
interface is:

```text
--backend auto|builtin|lapack
```

`auto` is the default. It selects the original built-in implementation for
orders below 128 and LAPACK at order 128 or above when LAPACKE is available.
`builtin` always selects the unchanged scalar implementation. An explicit
`lapack` request fails clearly when support was not compiled in.

The LAPACK implementation uses the typed LAPACKE C interface. It factors NEC's
existing transposed storage directly and uses a transpose solve, avoiding an
extra matrix copy. Multi-right-hand-side classic solves are submitted as one
batched `zgetrs` operation. Dense dielectric and coupled systems use the same
backend selection.

Selected and requested backends are visible in text output, result JSON,
timing JSON, and `--capabilities-json`. The benchmark runner accepts
`--backend` and `--threads`, with the latter setting the child BLAS/OpenMP
environment for controlled experiments.

## Build dependency and fallback

The validated host has Ubuntu `liblapacke-dev` and `liblapacke` 3.10.0-2ubuntu1.
Configure checks both `lapacke.h` and `LAPACKE_zgetrf`, then runs a nonsymmetric
complex-double factor/solve probe. No untyped Fortran ABI is used.

These configurations were validated:

- normal auto-detected LAPACKE build, linked through `liblapacke` to OpenBLAS;
- `./configure --with-lapack=no`, with no LAPACK link dependency; and
- explicit `--with-lapack=yes`, which fails configure if compatible LAPACKE is
  unavailable.

## Performance results

The retained report is `tests/performance/m2-dense-backend.json`. It uses 2
warmups and 11 measured repetitions per backend, interleaves built-in and
LAPACK runs at each size, and fixes the crossover study to one BLAS thread.
The host remained on its `powersave` governor, so medians rather than individual
runs are authoritative.

| Matrix order | Auto backend | Built-in factor | LAPACK factor | Factor speedup |
| ---: | --- | ---: | ---: | ---: |
| 31 | built-in | 0.000056 s | 0.000204 s | 0.27x |
| 63 | built-in | 0.000376 s | 0.000382 s | 0.99x |
| 95 | built-in | 0.001206 s | 0.000808 s | 1.49x |
| 127 | built-in | 0.002793 s | 0.001551 s | 1.80x |
| 129 | LAPACK | 0.002951 s | 0.001631 s | 1.81x |
| 191 | LAPACK | 0.009183 s | 0.004425 s | 2.07x |
| 255 | LAPACK | 0.016407 s | 0.003830 s | 4.28x |
| 401 | LAPACK | 0.016863 s | 0.007191 s | 2.35x |
| 801 | LAPACK | 0.136205 s | 0.052406 s | **2.60x** |

This meets the declared target of at least 2x on the 801-order matrix. Auto
retains the built-in path through the small/near-crossover cases where LAPACK
startup cost can dominate.

The 801-order LAPACK factorization showed these thread results:

| Threads | Factor time | Factor speedup vs 1 | Total wall time |
| ---: | ---: | ---: | ---: |
| 1 | 0.052438 s | 1.00x | 0.128710 s |
| 2 | 0.030709 s | 1.71x | 0.109877 s |
| 4 | 0.015445 s | 3.40x | **0.087930 s** |
| 8 | 0.013754 s | 3.81x | 0.091058 s |
| 16 | 0.010216 s | 5.13x | 0.132862 s |
| 32 | 0.044707 s | 1.17x | 0.166027 s |

More BLAS workers do not monotonically improve end-to-end time because matrix
fill remains serial and thread startup/coordination matters. Milestone 3 will
therefore expose coherent thread and job controls rather than making a
Milestone 2 thread-count policy.

## Numerical validation

`tests/backend_equivalence.py` compares complete machine-readable results for
nine representative models: classic examples, patches, networks, the Antenna
Turtle loop-Yagi, the 801-segment scaling model, and a dense dielectric model.
It uses relative tolerance `1e-8` and absolute tolerance `1e-10`; all cases
pass. The maximum relative feed-impedance difference in the crossover sweep is
`4.50e-11`.

The explicit LAPACK path was also rerun against both external numerical
oracles. `tests/performance/m2-cross-solver.json` passes the frozen Milestone 1
tolerances:

| Model | NEC2C LAPACK impedance (ohm) | AEG relative delta | momwire relative delta |
| --- | --- | ---: | ---: |
| dipole | 67.634642159 - j31.287356736 | 4.94e-4 | 1.24e-3 |
| Yagi | 58.020328022 + j1.162785893 | 6.98e-4 | 1.70e-3 |

AEG was executed only as an external GPL oracle; no AEG source or object code
was copied. Momwire remains an independent MIT-licensed oracle.

## Automated verification

- Normal GCC `-O2` build: all 28 Automake tests pass.
- `make distcheck`: out-of-tree build, all 28 tests, install/uninstall, and
  distribution completeness pass.
- LAPACK-disabled GCC `-O2` build: focused backend, JSON, and regression tests
  pass; no LAPACK dependency is linked.
- ASan+UBSan build: the five backend/JSON/regression tests pass with halt-on-
  error enabled. Leak detection was disabled because LeakSanitizer cannot run
  under the host application's ptrace environment.
- The LAPACK unit test covers nonsymmetric complex matrices, NEC's strided
  transposed storage, row-major VIE storage, pivot conversion, and strided
  batched right-hand sides.

No in-tree parallel code is introduced in this milestone; concurrency is
inside the externally validated BLAS library. Race-enabled in-tree testing
begins with Milestone 3.

## Manual acceptance scenarios

Build and confirm discovery:

```sh
autoreconf -fi
./configure
make -j8
./nec2c --capabilities-json | python3 -m json.tool
```

The capability document should list `auto`, `builtin`, and `lapack`, with
`default_dense_backend` equal to `auto` and `lapack_auto_min_order` equal to
128.

Run the large model through all three choices with one BLAS thread:

```sh
OPENBLAS_NUM_THREADS=1 ./nec2c --backend builtin \
  --timing-json /tmp/m2-builtin.json \
  -i tests/performance/scaling-dipole-801.nec -o /tmp/m2-builtin.out

OPENBLAS_NUM_THREADS=1 ./nec2c --backend lapack \
  --timing-json /tmp/m2-lapack.json \
  -i tests/performance/scaling-dipole-801.nec -o /tmp/m2-lapack.out

OPENBLAS_NUM_THREADS=1 ./nec2c --backend auto \
  --timing-json /tmp/m2-auto.json \
  -i tests/performance/scaling-dipole-801.nec -o /tmp/m2-auto.out

grep 'DENSE BACKEND\|FILL:' /tmp/m2-builtin.out /tmp/m2-lapack.out /tmp/m2-auto.out
```

Expected: explicit choices report their selected backend; auto reports LAPACK;
LAPACK factorization is materially faster; feed impedance remains about
`66.232419 - j53.177760` ohm.

Confirm the small-case fallback:

```sh
./nec2c --backend auto -i Input/EX1.nec -o /tmp/m2-small.out
grep 'DENSE BACKEND' /tmp/m2-small.out
```

Expected: `builtin (REQUESTED auto)`.

For visual validation, run the accepted Antenna Turtle model with `auto` and
`lapack`, open both output patterns in the normal NEC2C viewer, and confirm the
same approximately 14.5 dBi forward lobe:

```sh
./nec2c --backend auto -i tests/performance/loop-yagi-438.nec -o /tmp/m2-yagi-auto.out
./nec2c --backend lapack -i tests/performance/loop-yagi-438.nec -o /tmp/m2-yagi-lapack.out
```

Optionally reproduce the complete crossover/thread report:

```sh
python3 tools/dense_backend_benchmark.py --repeat 7 --warmup 2 \
  --output /tmp/m2-dense-backend.json
```

Return `PASS`, `FAIL`, or `PASS WITH CHANGES`. Milestone 3 does not begin until
that verdict is received.
