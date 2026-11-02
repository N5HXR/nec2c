# Performance Milestone 1 results

Status: implementation complete; awaiting user acceptance.

Date: 2026-11-02  
Branch: `codex/performance`  
Preserved reference: `1166bc3bc5ba490c56fadf946070998a47ff8e00`

Milestone 1 adds measurement and validation infrastructure only. It does not
change a solver algorithm, tolerance, arithmetic precision, or backend.

## Delivered interface

`nec2c` now accepts:

```text
--timing-json FILE
```

The option writes a separate, atomic timing artifact described by
`docs/nec2c-timing-v1.schema.json`. It deliberately does not change result
schema version 1. `--capabilities-json` advertises timing schema version 1.
Timing and result destinations must be different files, and `-` is rejected as
a timing destination so normal solver output cannot be mixed with measurement
data.

The artifact records monotonic wall time, process CPU time, peak RSS, completion
status, call counts, and applicable coarse phases. The current phases are
`ground_table`, `matrix_fill`, `matrix_factor`, `matrix_solve`, `near_field`,
`far_field`, `dielectric_mesh`, and `dielectric_solve`. Phase boundaries do not
change arithmetic order. When timing is not requested, the clock reader returns
without issuing a system clock call.

## Reproducible corpus and runners

`tools/performance_benchmark.py` runs warmups and measured repetitions from
`tests/performance/models.json`. It records the executable and deck hashes,
host identity, individual timing artifacts, median/minimum/maximum values, and
a canonical numerical-result hash with timing-only fields removed.

The frozen corpus contains small and mixed legacy examples, the repository's
19-element Yagi, dense and FFT dielectric cases, a deterministic 801-segment
scaling dipole, and the user-approved Antenna Turtle loop-Yagi. The retained
loop-Yagi deck has 324 `GW` records and SHA-256:

```text
96fcdf72c16095ba617297072d746a34c9b98fc670ea7fa9c111ebb787d620fb
```

`tools/prepare_nec_fixture.py` reproducibly adds 438 MHz and a 5-degree
full-sphere `RP` grid to an Antenna Turtle geometry export. It embeds the source
export's SHA-256 in the fixture without modifying the Antenna Turtle project.

`tools/timing_overhead.py` interleaves reference, timing-disabled, and
timing-enabled executions to reduce temperature and scheduling bias.

## Reference measurements

The retained seven-run warm report is `tests/performance/m1-reference.json`.
Every repeated case produced the same canonical numerical hash.

| Case | Median wall | Peak RSS | Dominant measured phases |
| --- | ---: | ---: | --- |
| classic small | 0.000698 s | 17,408 KiB | near field 0.000134 s |
| classic mixed | 0.000873 s | 17,408 KiB | far field 0.000203 s |
| 19-element Yagi | 0.033825 s | 17,664 KiB | fill 0.013690 s; factor 0.011952 s |
| Antenna Turtle loop-Yagi | 0.052093 s | 23,424 KiB | far field 0.025257 s; fill 0.011939 s; factor 0.008693 s |
| scaling dipole, 801 segments | 0.207503 s | 23,432 KiB | factor 0.131076 s; fill 0.072135 s |
| dielectric dense | 0.035134 s | 23,432 KiB | dielectric solve 0.033701 s |
| dielectric FFT | 0.025962 s | 23,432 KiB | dielectric solve 0.024904 s |

These short cases define repeatable baselines, not performance claims. The host
reported its `powersave` CPU-frequency governor, and the shortest dielectric
times vary materially between runs. Larger production fixtures will be added
when a milestone needs to establish a backend crossover or scaling curve.

Against the archived Milestone 0 executable, 25 interleaved repetitions of the
801-segment case measured:

| Mode | Median overhead versus M0 | Limit |
| --- | ---: | ---: |
| timing disabled | -1.611% | at most +2% |
| timing enabled | -0.910% | at most +5% |

Negative values are measurement noise, not an optimization. Both satisfy the
contract, and all 75 normalized solver outputs were identical. Raw samples are
retained in `tests/performance/m1-overhead.json`.

## Cross-solver numerical evidence

`tools/cross_solver_compare.py` treats AEG as an external GPL executable and
imports no AEG source or object code. Momwire is loaded from its own virtual
environment and uses its independent NEC2 sinusoidal implementation. Two
compatible free-space cases compare complex feed impedance:

| Model | nec2c | AEG | momwire sinusoidal |
| --- | ---: | ---: | ---: |
| 21-segment dipole | 67.634642 - j31.287357 Ω | 67.6397 - j31.2509 Ω | 67.622162 - j31.378545 Ω |
| two-element Yagi | 58.020328 + j1.162786 Ω | 58.0400 + j1.19817 Ω | 57.979488 + j1.073235 Ω |

The worst relative complex difference is 0.0698% for AEG and 0.1696% for
momwire. Frozen acceptance limits are 0.1% for AEG and 0.5% for momwire. AEG's
ordinary result file prints impedance to six significant digits, limiting the
available comparison precision. Momwire shares the three-term sinusoidal basis
but has independent geometry and integration code, so its threshold is
intentionally broader. The complete evidence is retained in
`tests/performance/m1-cross-solver.json`.

AEG required no system installation. Its old CMake project had to use the
installed matching GCC 11 C and Fortran pair; GCC 12 C objects and gfortran 11
objects could not share link-time-optimization bytecode. It then selected the
already-installed OpenBLAS/LAPACK. Momwire's existing virtual environment
passed `momwire-nec2c --selftest`.

## Automated verification

The normal GCC `-O2` build passes all 26 Automake tests, including the new
timing artifact contract test:

```text
TOTAL: 26  PASS: 26  FAIL: 0  ERROR: 0
```

An isolated GCC build with AddressSanitizer and UndefinedBehaviorSanitizer also
passes all 26 tests. Leak detection was disabled because LeakSanitizer cannot
run under this environment's ptrace supervision. No race test is applicable in
this milestone because it adds no threads or shared concurrent execution.

## Manual acceptance scenarios

Build and run the complete suite:

```sh
./autogen.sh
./configure
make -j8
make check
```

Inspect the loop-Yagi phase breakdown:

```sh
./nec2c --timing-json /tmp/m1-loop-timing.json \
  --result-json /tmp/m1-loop-result.json \
  -i tests/performance/loop-yagi-438.nec -o /tmp/m1-loop.out
python3 -m json.tool /tmp/m1-loop-timing.json
```

Run the complete retained performance corpus:

```sh
python3 tools/performance_benchmark.py --repeat 7 --warmup 2 \
  --output /tmp/m1-benchmark.json
```

Regenerate the loop-Yagi fixture after running the modified Antenna Turtle
example that exports `/tmp/loop-yagi.nec`:

```sh
python3 tools/prepare_nec_fixture.py /tmp/loop-yagi.nec \
  /tmp/loop-yagi-438.nec --frequency-mhz 438
sha256sum /tmp/loop-yagi-438.nec
```

The expected generated SHA-256 is the fixture hash above. For manual visual
validation, open `tests/performance/loop-yagi-438.nec` in the normal NEC2C
pattern viewer and compare the forward-lobe shape with the accepted screenshot.
The retained 5-degree grid reports a 14.47 dBi peak, consistent with the
approximately 14.5 dBi screenshot.

## Known limitations and next gate

The benchmark report captures the requested phase data but does not yet expose
backend and thread fields because Milestone 1 has no alternative backend or
thread-control interface. Those fields become mandatory when Milestone 2 adds
the dense backend selection. The dielectric phase currently measures the full
solve; FFT planning/application and conductor-volume coupling are split in the
milestones that modify those paths.

No Milestone 2 implementation begins until the user returns `PASS`, `FAIL`, or
`PASS WITH CHANGES` for this milestone.
