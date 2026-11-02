# Performance Milestone 0 baseline

Status: implementation complete; awaiting user acceptance.

Date: 2026-11-02  
Branch: `codex/performance`  
Reference commit: `1166bc3bc5ba490c56fadf946070998a47ff8e00`

Milestone 0 changes documentation and branch state only.  It makes no source,
algorithm, numerical, build-flag, or expected-result change.

## Preservation and build

The workspace was already configured in place, so a conventional out-of-tree
configure correctly refused to reuse it.  The clean validation instead used a
`git archive` of the exact reference commit in a temporary directory, followed
by the normal sequence:

```sh
./autogen.sh
./configure
make -j8
make check
```

This excludes untracked workspace files and avoids modifying the existing
configured source tree.  The untracked workspace `bugs.md` was not read,
modified, added, or removed.

Both the existing workspace executable and the executable from the clean
archived build have SHA-256:

```text
c63fa74641506f18743083e89244dbf5a171f62ccb16ef7933e51279757c5aaa
```

The baseline links only `libm`, `libc`, and the ELF loader.  It was compiled by
GCC 12.3.0 with the configured baseline flags `-Wall -O2`.

## Host record

* Linux 6.8.0-138-generic, x86-64
* Intel Core i9-13900K
* 24 physical cores, 32 online logical CPUs, one NUMA node
* 36 MiB L3 cache
* GNU Make 4.3
* Python 3.10.12
* Racket CS 8.2

Milestone 1 will capture frequency policy, affinity, library versions, repeated
samples, and monotonic phase times.  The single-run times below are smoke
measurements, not optimization claims.

## Automated result

The clean archived build passed the complete Automake suite:

```text
TOTAL: 25
PASS:  25
SKIP:  0
XFAIL: 0
FAIL:  0
XPASS: 0
ERROR: 0
```

This includes legacy NEC regression, dielectric parser and meshing, dense and
iterative solvers, wire and patch coupling, result JSON, cuboid kernels, FFT,
perfect and real ground, loss, material tables, touching materials, and the
retained tree report.

## Representative smoke measurements

| Model | Wall | User CPU | Peak RSS | Numerical note |
| --- | ---: | ---: | ---: | --- |
| `Input/EX1.nec` | <0.01 s | <0.01 s | 2,560 KiB | two legacy matrix fills |
| `Input/EX4.nec` | <0.01 s | <0.01 s | 2,816 KiB | mixed geometry/perfect ground |
| `Input/19el_432_tonna.nec` | 0.04 s | 0.04 s | 5,120 KiB | reported fill/factor 20/20 ms |
| `tests/dielectric/m6-cube.nec`, dense | <0.01 s | <0.01 s | 2,816 KiB | relative residual 7.6991e-15 |
| same model, FFT GMRES | 0.01 s | 0.01 s | 2,816 KiB | true relative residual 2.15045e-14 |

Raw output hashes from this run are retained here to identify the exact manual
acceptance evidence; Milestone 1 replaces these coarse timings with a
reproducible benchmark record.

```text
40830465aabbfe02f5e71fcbce3de6ae436eac4729bd4f8ba3cf05916c39b72f  EX1 output
2d4c56972824fdf5b0e2d83c08c885b010b13a70f670ad5ef2b1f970d493ff69  EX4 output
3217135c668683d211ef9091f7c4b929b6b6b061853c88f0efdf6bffe06e8c40  19-element output
f43d3c99bd738c22c7977c86121e174dcd2b36770168020f681df78488b30493  dense dielectric text
940b1d617b1e8fe9c70055d54d1373caa0fb0f2da4c22a8d79cc001af9f09e0c  dense dielectric JSON
eee8ed7230fb907e546ba804d92eed3b14b09fa5932ce427bf82f90e1f7ef43a  FFT dielectric text
8a0d2ff8eb93298a57d2e9a658e91cf26eeacb61b82b37a396b5ebca6949f976  FFT dielectric JSON
```

## Antenna Turtle fixture confirmation

The user-provided `examples/loop-yagi.rkt` export was exercised without
modification.  It produced `/tmp/loop-yagi.nec` with 324 `GW` records and one
driven source.  The generated resource intentionally contains geometry and its
source but no frequency or requested field grid.  Milestone 1 will apply
versioned deterministic analysis-card overlays rather than alter the exported
geometry.  The user's rendered pattern, with a reported peak near 14.5 dBi, is
the initial manual pattern reference.

## Manual acceptance scenarios

From the nec2c repository, rebuild and run the complete baseline:

```sh
./autogen.sh
./configure
make -j8
make check
```

Run a medium legacy wire model:

```sh
./nec2c -i Input/19el_432_tonna.nec -o /tmp/m0-19el.out
grep -A2 'MATRIX TIMING' /tmp/m0-19el.out
```

Compare the small dense and FFT dielectric paths:

```sh
./nec2c -s dense -i tests/dielectric/m6-cube.nec \
  -o /tmp/m0-dense.out --result-json /tmp/m0-dense.json
./nec2c -s fft -i tests/dielectric/m6-cube.nec \
  -o /tmp/m0-fft.out --result-json /tmp/m0-fft.json
grep 'RELATIVE RESIDUAL' /tmp/m0-dense.out /tmp/m0-fft.out
```

From the Antenna Turtle repository, reproduce the accepted loop-Yagi export:

```sh
racket examples/loop-yagi.rkt
grep -c '^GW' /tmp/loop-yagi.nec
```

The expected wire-record count is 324.  The existing pattern screenshot may be
used as the Milestone 0 manual visual reference; precise sampled pattern
comparisons begin in Milestone 1.

## Exit criteria

Milestone 0 is accepted when the user confirms that the baseline suite and
manual scenarios are credible, the selected branch/reference preserve the old
implementation adequately, and the fixture direction is appropriate.  No
Milestone 1 work begins before that verdict.
