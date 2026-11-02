# NEC2C performance roadmap

Status: approved for implementation on 2026-11-02.

This roadmap governs performance work on the `codex/performance` branch.  The
immutable implementation reference is commit
`1166bc3bc5ba490c56fadf946070998a47ff8e00`.  Optimized implementations remain
subject to explicit user acceptance at every milestone.

## Governing decisions

* A validated faster backend becomes the default through automatic selection.
* Command-line switches must retain access to the original scalar or serial
  implementation.
* Configuration is exposed through command-line options.  The intended common
  controls are `--backend auto|builtin|lapack`,
  `--fft-backend auto|builtin|fftw`, `--threads N`, and `--jobs N`.  Exact
  spellings are reviewed when the corresponding milestone begins.
* Existing numerical baselines are not regenerated merely because an optimized
  backend differs.
* New system software is installed only by the user.  Before a dependency is
  needed, development stops and reports the package, version, discovery name,
  and a verification command.
* Source from aegnec2 is not copied.  AEG is an external behavioral and
  numerical oracle.  Momwire is MIT licensed, but borrowed source must remain
  isolated and retain its required notice.

## Acceptance protocol

Each milestone ends with:

1. a design and change summary;
2. automated numerical, sanitizer, and race results appropriate to the change;
3. phase-level performance and peak-memory results;
4. comparison with the built-in implementation and, where compatible,
   momwire and aegnec2;
5. known limitations and regressions; and
6. exact models and commands for manual validation.

Development does not begin on the next milestone until the user returns
`PASS`, `FAIL`, or `PASS WITH CHANGES`.  Failed milestones are corrected and
presented again.

## Measurement contract

Performance reports include monotonic wall time, process CPU time, peak RSS,
backend, thread count, compiler and flags, linked libraries, CPU affinity, and
cold and warm repetitions.  Relevant phases include setup, ground tables,
matrix or kernel fill, factorization, solve, FFT planning and application,
conductor-volume coupling, and field or power output.

Thread studies use 1, 2, 4, 8, 16, and 32 threads unless the hardware or phase
makes a point inapplicable.  Automatic selection must not regress small target
workloads by more than 5% outside measured noise.  Performance is never gained
by relaxing an accuracy tolerance without separate approval.

For the same discretized equations, the built-in implementation is the primary
regression oracle.  All existing tests must pass unchanged.  Comparisons use a
scaled residual and complex currents, impedance, fields, and power rather than
requiring byte identity from different floating-point reduction or pivot
orders.  Cross-solver momwire and AEG thresholds are frozen from the Milestone
1 evidence because their formulations and output conventions are not identical
in every mode.

## Milestones

### M0: baseline freeze and branch

Record the immutable reference, create the performance branch, perform a clean
archived build, run the complete suite, retain representative measurements,
and document manual reference scenarios.  No numerical implementation changes.

### M1: benchmark and cross-solver harness

Add monotonic phase timers and machine-readable benchmark records.  Establish
cold and warm runners, the Antenna Turtle fixture corpus, compatible momwire
models, compatible AEG decks, and frozen numerical tolerances.  Timing overhead
must be below 2% when disabled and 5% when enabled.

### M2: optional LAPACK dense backend

Introduce a factor/solve abstraction for classic NEC, dense VIE, and dense
coupled systems.  Preserve the built-in solver, support batched right-hand
sides, report the selected backend, and choose the optimized backend by default
only above validated crossover sizes.  Target at least a 2x factorization
speedup on a declared large matrix without a small-case regression.

### M3: concurrency contract and frequency jobs

Implement deterministic process-parallel frequency execution and one coherent
thread/job control model.  Prevent nested BLAS, OpenMP, and FFT oversubscription.
Per-frequency numerical output must match separate serial executions.

### M4A: parallel built-in FFT

Parallelize independent FFT-line batches and tensor/grid work while preserving
each line's arithmetic order.  Retain the serial path and size threshold.
Target at least 1.5x on the declared large grid.

### M4B: optional threaded FFTW backend

Add FFTW as a selectable backend and retain the built-in FFT as reference and
fallback.  FFTW becomes automatic only if it beats the parallel built-in path
by an agreed material margin, initially 20%, on production-scale targets.

### M5: parallel dielectric assembly

Parallelize cell-pair and reentrant volume, wire-row, and patch-row block
assembly using disjoint output ownership and deterministic failure reporting.
Target at least 2x on the selected assembly workload.

### M6: reentrant classic NEC fill

Move mutable interaction scratch into worker-local contexts and parallelize
source columns or tiles without changing per-entry summation order.  Cover
free-space, perfect and real ground, extended kernels, patches, symmetry,
networks, and loading.  Target at least 2x fill speedup on a large model.

### M7: parallel conductor-volume coupling

Remove temporary global field/current mutations, pass explicit source state,
and parallelize disjoint source/cell blocks.  Require a material end-to-end
improvement on the coupled dielectric target as well as a 2x coupling-phase
goal.

### M8: reuse, batching, caching, and field outputs

Reuse factorizations and solved currents, batch right-hand sides, add exact-key
in-process setup caches, and evaluate independent field points in parallel with
deterministic serialized reductions and output.  A persistent cache requires a
separate design review.

### M9: compiler, memory, and release candidate

Evaluate conservative optimization and CPU dispatch without fast-math, remove
avoidable copies and allocations, tune crossover thresholds and worker waiting,
run the complete correctness and performance matrix, and prepare the final
manual acceptance bundle before merge.

## Production acceptance fixtures

Generated NEC resources are retained with their Antenna Turtle source path and
source hash so performance tests do not require Racket at execution time.  The
primary proposed sources are `dipole.rkt`, `yagi.rkt`, `loop-yagi.rkt`,
`folded-dipole-with-transformer.rkt`, `nec2c-car-roof-monopole.rkt`,
`nec2c-dielectric-quarter-panel.rkt`, and the HF/VHF/UHF scenarios in
`n5hxr.rkt`.  Milestone 1 freezes the final small, medium, and production sets.

The repeated-element/block-Toeplitz array solver, generic H-matrix methods,
GPU execution, and distributed execution remain separate research tracks and
are not allowed to expand these milestone scopes implicitly.
