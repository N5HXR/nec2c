# NEC4 similarity M0 baseline

Status: accepted by the user on 2026-11-16.

Date: 2026-11-16  
Branch: `codex/nec4-similarity`  
Reference commit: `9f7f108456a28b3f961403d9bad2769ad6c55fcf`

## Milestone result

M0 changes project documentation and test infrastructure only. It makes no production numerical, algorithm, input, output, or default-selection change.

The accepted source was exported with `git archive` into a clean temporary directory, configured by `./autogen.sh`, built with `make -j16`, and tested with `make check -j16`. Untracked workspace documents and reference PDFs were not included in the build.

The clean suite result was:

```text
TOTAL: 46
PASS:  46
SKIP:   0
XFAIL:  0
FAIL:   0
XPASS:  0
ERROR:  0
```

After adding the semantic baseline runner, the feature-branch suite also passed `47/47`; the additional test is precisely the new classic free-space/real-ground preservation check.

## Preservation record

```text
reference commit:
9f7f108456a28b3f961403d9bad2769ad6c55fcf

git archive SHA-256:
cddf30dc8dc86e27f0a8d35760623d16036460a6efc3ec99fe07d1659c82a47d

clean reference executable SHA-256:
709c2e7f5daecc4bffca5b5a00fa79eca0fb172fefd01c6afd982e50f01a0057
```

The immutable commit is the authoritative source reference. The executable hash identifies the exact M0 build on this host; executable byte identity is not expected across compiler or library changes. `tests/nec4-similarity/m0_baseline.py` instead verifies a stable physical projection of classic results and normalized traditional reports.

The runner automatically uses:

```text
--wire-engine nec2 --ground-engine nec2
```

once those selectors exist. Before M1 introduces them, it runs the pre-selector classic implementation. This keeps the baseline useful after NEC4 becomes an available or default engine.

## Host and toolchain

- Ubuntu 22.04 userland; Linux `6.8.0-138-generic`, x86-64.
- Intel Core i9-13900K, 24 physical cores, 32 online logical CPUs, one NUMA node, 36 MiB L3.
- GCC 12.3.0, GNU Make 4.3, Python 3.10.12.
- Build flags: `-fopenmp -Wall -O3`.
- LAPACKE 3.10.0, OpenBLAS pthread 0.3.20, FFTW 3.3.8, GNU OpenMP 12.3.
- Configured optional capabilities: LAPACK dense backend, threaded FFTW, and OpenMP built-in FFT.

The reference executable dynamically links FFTW, LAPACKE, OpenBLAS/BLAS/LAPACK, GNU OpenMP, libm, libc, libgfortran, libquadmath, and the ELF loader. The exact `ldd` paths are host-specific and are not treated as numerical inputs.

No new library is required by M0.

## Frozen semantic witnesses

The golden file is `tests/nec4-similarity/m0-baseline.json`. Its physics hash includes complex feeds, currents, fields, power, material evaluations, excitations, triggers, and resolved ground state. It excludes elapsed time, solver/backend identity, executable path, and other nonphysical metadata.

### Free-space driven dipole and radiation grid

Deck: `tests/coated-wire/m0-bare-dipole.nec`

- Frequencies: 290, 299.896229, and 309.792458 MHz.
- Total wire-current samples: 303.
- Total complex far-field samples: 2,109.
- Physics SHA-256: `2bc268a240ae7b109264900497a13423ffae5f6c3454533b7f9ac995f04e38f3`.
- Normalized traditional-report SHA-256: `c4299cd2105827185d4a5d533cacb27085f23bfdb0bb69d37939987e91fb1d69`.

Feed impedances are:

| MHz | Resistance (ohm) | Reactance (ohm) |
| ---: | ---: | ---: |
| 290 | 74.82978708754783 | 12.398158167288537 |
| 299.896229 | 83.8410060490305 | 48.26302302343913 |
| 309.792458 | 93.992160363995 | 84.46686220252845 |

### Real-ground frequency sweep

Deck: `tests/performance/m3-real-ground-sweep.nec`

- Frequencies: 80 through 120 MHz in 5 MHz steps.
- Total wire-current samples: 3,609.
- Physics SHA-256: `df018166a78a9f50bc7732645bb453863e70f67ae66d2226d74b4d751400bf66`.
- Normalized traditional-report SHA-256: `a2d0f488173ded8d7becc179d7be161ecb8e1203c065f20a395b8b556587abf2`.

Feed impedances are:

| MHz | Resistance (ohm) | Reactance (ohm) |
| ---: | ---: | ---: |
| 80 | 99.7320371578203 | 128.35678540200385 |
| 85 | 126.50752446140723 | 216.26345742195673 |
| 90 | 162.43624132528004 | 310.2930701407098 |
| 95 | 211.46634032364585 | 412.67708860814383 |
| 100 | 279.7745714795854 | 525.7620172711789 |
| 105 | 377.4101447312206 | 651.5164129904587 |
| 110 | 521.1912349555377 | 789.7457461051658 |
| 115 | 739.4654137387879 | 932.4490239370026 |
| 120 | 1076.6450209716754 | 1047.2373208994954 |

These are preservation witnesses, not claims that the classic answer is the physically correct NEC4 target.

## Public reference-model inventory for M1

The following families are sufficiently important to reconstruct or formalize in M1:

| Source | Candidate witnesses | Intended use |
| --- | --- | --- |
| 1992 NEC4 theory | unequal-radius junctions, junction-potential solve, short-segment field limits, flat end caps, buried/cross-interface field regions | Equation-level high-precision oracles and branch tests |
| 2008 kernel accuracy report | thick half-wave dipole, square loop, close parallel-wire line | RTWK/ETWK/exact-kernel distinction and validity limits |
| 2011 GN3 report | 5 m whip with 0.2 m seawater ground stake; short vertical dipole near lossless water | Original-ground versus adaptive-ground comparisons and difficult H/interference regions |
| 2019 NEC5 validation manual | dipole convergence, small loops, horizontal/long wires over ground, radial screens | Independent published tolerance bands and known NEC4 limits |
| Existing repository | Antenna Turtle exports, performance real-ground sweeps, coated-wire and dielectric regressions | Backward compatibility, realistic workloads, VSWR/radiation-grid coverage |

M1 will record full geometry and excitation provenance before promoting any plot or reported number to an acceptance witness. Incomplete figures remain qualitative evidence.

## Manual acceptance

Build the feature branch normally:

```sh
./autogen.sh
make -j16
```

Run the frozen semantic witnesses:

```sh
python3 tests/nec4-similarity/m0_baseline.py
```

Expected output:

```text
NEC4-similarity M0 classic free-space and real-ground baselines are unchanged
```

Optionally run the complete suite:

```sh
make check -j16
```

The branch adds the M0 baseline test, so the expected post-M0 suite total is 47 rather than the clean reference build's 46.

For direct inspection of the two reports:

```sh
./nec2c --backend builtin --fft-backend builtin --jobs 1 --threads 1 \
  -i tests/coated-wire/m0-bare-dipole.nec -o /tmp/nec4-m0-free-space.out \
  --result-json /tmp/nec4-m0-free-space.json

./nec2c --backend builtin --fft-backend builtin --jobs 1 --threads 1 \
  -i tests/performance/m3-real-ground-sweep.nec \
  -o /tmp/nec4-m0-real-ground.out \
  --result-json /tmp/nec4-m0-real-ground.json
```

## Acceptance gate

M0 passes when the user confirms that:

1. the classic source and behavior are preserved adequately;
2. the free-space and real-ground witnesses execute successfully;
3. the public reference-model inventory is an appropriate foundation for M1; and
4. the feature may proceed to independent numerical oracles and the operator boundary.

No M1 implementation begins before that verdict.

The user returned `PASS` on 2026-11-16 and authorized M1.
