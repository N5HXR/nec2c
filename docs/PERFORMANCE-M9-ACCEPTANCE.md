# Performance Milestone 9 acceptance bundle

Status: implementation and developer validation complete; awaiting user
acceptance.

This is the final gate for the `codex/performance` branch.  It measures the
optimized defaults, the selectable built-in single-worker behavior, and an
executable built from the immutable pre-optimization commit.  It also reruns
the independent AEG and momwire numerical witnesses.  AEG remains an external
GPL oracle; no AEG source or object code enters this repository.

## Build the preserved reference

Keep the reference build outside the optimized checkout.  The runner expects
commit `1166bc3bc5ba490c56fadf946070998a47ff8e00`, which is the revision frozen
in M0.  One reproducible preparation is:

```sh
mkdir -p /tmp/nec2c-m9-reference
git archive 1166bc3bc5ba490c56fadf946070998a47ff8e00 |
  tar -x -C /tmp/nec2c-m9-reference
cd /tmp/nec2c-m9-reference
./autogen.sh
./configure
make -j8
```

Return to the optimized checkout and ensure its normal clean build is current.
The final run uses existing LAPACKE, OpenMP and FFTW discovery; it does not
install or download anything.

## Automated release-candidate matrix

Run the complete matrix with two warmups and seven measured repetitions:

```sh
python3 tools/m9_release_candidate.py \
  --reference-binary /tmp/nec2c-m9-reference/nec2c \
  --aegnec2 /path/to/aegnec2 \
  --momwire-root "$HOME/src/momwire" \
  --full --warmup 2 --repeat 7 \
  --output-directory /tmp/nec2c-m9-acceptance
```

The runner performs and indexes:

* the complete Automake correctness suite;
* the nine-case performance corpus under automatic defaults;
* the same corpus under `--backend builtin --fft-backend builtin --jobs 1
  --threads 1`;
* the corpus with the M0 executable and its historical command-line interface;
* the AEG/momwire impedance witnesses and their frozen 0.1%/0.5% limits;
* FFTW versus built-in FFT at 1, 16, and 32 threads; and
* the M5 through M8 phase curves at 1, 2, 4, 8, 16, and 32 workers.

`m9-release-candidate.json` records every command, return code, captured output,
artifact hash, and whether all mandatory evidence is present.  The current
built-in single-worker corpus must have exact canonical result hashes matching
the M0 executable.  Optimized LAPACK and FFTW paths are instead governed by the
tolerance-based backend tests because their legal floating-point order and
pivot choices differ.

The top-level `acceptance_ready` field must be true.  Each performance corpus
case must report a stable numerical result.  The automatic defaults must not
regress a small corpus case by more than 5% outside measurement noise.  The
runner enforces that gate on solver timing against the current explicitly
selected built-in single-worker path; overlapping interquartile ranges count
as noise.  Minimum/maximum overlap is deliberately insufficient, so one
scheduling outlier cannot conceal a bad median.  It separately discloses total
CLI wall time against M0, including the
startup/shutdown cost of today's dynamically linked OpenMP, LAPACK and FFTW
libraries.  M0 predates timing JSON, so that CLI comparison is not silently
substituted for the backend crossover gate.  Larger
cases are reviewed against both M0 and the current built-in single-worker path;
phase speedups must remain consistent with the previously accepted milestone
results.  Peak process-tree RSS is compared as an absolute value and ratio;
any unexplained material increase is a release blocker.

## Manual models and user-visible controls

These checks deliberately exercise the workflows most likely to reveal a bad
automatic crossover or an ordering problem.

1. Run and view `tests/performance/loop-yagi-438.nec`.  Confirm the forward
   lobe is unchanged and the peak remains approximately 14.47 dBi.  Compare
   automatic execution with `--backend builtin --fft-backend builtin --jobs 1
   --threads 1`.
2. Run `tests/performance/scaling-dipole-801-sweep.nec` with automatic defaults
   and then with the same legacy controls.  Confirm frequency jobs retain deck
   order and feed impedances agree.
3. Run `tests/dielectric/m6-cube.nec` with `-s dense` and `-s fft`.  Confirm the
   reported residuals satisfy the deck tolerance and the requested backend is
   named in the timing artifact.
4. Run a representative Antenna Turtle model exported with
   `(export-nec2c "PATH")`, adding deterministic `FR`/`RP` analysis cards if
   necessary.  Confirm the optimized and legacy-control patterns agree and the
   automatic run uses the expected CPU resources.

The user returns `PASS`, `FAIL`, or `PASS WITH CHANGES` after reviewing the
JSON bundle and these model results.  No merge to the main development branch
occurs before that verdict.
