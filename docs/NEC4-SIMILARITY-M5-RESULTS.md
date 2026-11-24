# NEC4-similarity M5 results

Status: complete pending user acceptance.

Date: 2026-11-16

## Outcome

M5 closes the documented free-space NEC4 wire formulation with flat end caps
at free wire ends, voltage-source gap caps, and impedance-load gap caps. The
classic NEC2 path remains unchanged. The NEC4 path now carries cap fields
through matrix assembly, source excitation, load stamping, and reported near
electric fields.

This is the first candidate free-space NEC4-similarity preview. It does not
change the default `auto` engine and it does not include the M6-M10 ground
work.

## Mathematical implementation

The implementation follows the public NEC4 theory guide:

- section 4.4 equations (4-35) through (4-37) for constant charge on a flat
  free-end disk;
- Appendix E equations (E-1) through (E-3) for far/off-axis and exact-axis
  disk fields;
- equations (4-38) through (4-40) for voltage-source gap capacitance and cap
  fields;
- the load-gap rule immediately following (4-40), applied only to the
  constant-current coefficient at the loaded segment center; and
- equations (3-12), (3-13), and (4-36), giving the NEC4 free-end parameter
  `X_i = k_s a/2` (represented by `pi*a` in wavelength-normalized geometry).

Production evaluation uses the exact on-axis expression, independently
converged product Gauss quadrature in the near off-axis region, and the
Appendix E series only beyond 1,000 cap radii. The conservative far switch
keeps the measured scaled error below `1e-11`.

The theory guide explicitly says that NEC4 did not include annular charge
surfaces at radius steps. M5 therefore reports radius-step junctions but does
not invent an undocumented cap correction for them.

## Numerical evidence

`tests/nec4-similarity/m5_endcap_oracle.py` performs independent 50-digit
direct integration over the disk. `tests/nec4_endcap_test.c` checks the
production axis, adaptive quadrature, and far-series branches.

- exact-axis and near-quadrature witnesses agree at approximately machine
  precision;
- the retained far-series witness has scaled error `6.39e-12`;
- a two-dipole source-port swap has reciprocity defect `8.00e-18`;
- the 402-segment loaded/networked model is identical at one and eight fill
  workers to the retained numeric comparison (`0.0` maximum scaled delta);
- its accepted, radiated, structure-loss, and network-loss powers close to
  the precision reported in JSON (`0.0` measured relative defect).

The recreated source-driven quarter-wave monopole has radius `0.01 lambda`,
length `0.25 lambda`, and the guide's 80 segments (`L/a = 0.3125`). At that
point the normalized alternating-current ripple is:

| engine | ripple | input impedance (ohm) |
| --- | ---: | ---: |
| NEC2 | `1.525` | `12.0005 - j26.4057` |
| NEC4 | `0.09770` | `99.4284 + j39.6860` |

Thus the documented NEC4 treatment reduces this deliberately pathological
ripple by about 15.6 times. The retained sweep also includes a finer
`L/a = 0.2083` point, where ripple rises to `0.9916`. This is a clear
stability knee for the simplified constant-charge cap model. The existing
diagnostic remains deliberately conservative below `L/a = 2`; M5 records
`L/a = 0.3125` as the demonstrated extreme and makes no accuracy claim below
it.

## Source and compatibility audit

- Ordinary `EX 0` voltage sources receive the documented two-cap correction.
- Impedance loads receive the documented current-dependent cap field.
- Networks and maximum-coupling consumers use the corrected solved currents.
- Applied plane waves require no source-gap cap and retain the free-end caps.
- `EX 5` current-slope-discontinuity sources remain supported through their
  legacy construction and are explicitly counted as `LEGACY EX5 SOURCES`;
  M5 does not falsely claim the documented `EX 0` gap-cap correction for them.
- Source segments more than 2:1 different in length from a connected segment
  produce a warning and are rejected by `--model-check strict`.
- Radius-step junctions are reported but not rejected solely for lacking the
  annular correction that NEC4 itself omitted.

No card syntax or result field was removed or reinterpreted. Explicit
`--wire-engine nec2` remains the original path.

## Regression, sanitizers, and performance

- optimized suite: 57/57 passed;
- M4 kernel and model corpus: passed after the additive diagnostic format was
  updated;
- ASan/UBSan: end-cap unit oracle and complete M5 model corpus passed with
  LeakSanitizer disabled because it cannot run under the app host;
- accepted M4 versus M5 median time on the 402-segment loaded/network model:
  approximately `0.14 s` versus `0.15 s` at one worker, a measured 7.1%
  increase and below the 10% ordinary-model budget;
- peak RSS remained in the same `9,472-9,728 KiB` range.

For a thin 95-segment half-wave dipole, M5 moves the accepted M4 result from
`78.1463 + j14.2262 ohm` to `73.0037 + j4.4077 ohm`. This is expected cap and
source-gap physics, not a classic-path regression; the result moves toward the
ordinary thin half-wave reference behavior while preserving exact power
closure.

Retained evidence:

- `tests/nec4-similarity/m5-acceptance.json`
- `tests/nec4-similarity/m5-cap-convergence.png`
- `tests/nec4-similarity/m5-cap-monopole.nec`
- `tests/nec4-similarity/m5_endcap_oracle.py`
- `tests/nec4-similarity/m5_models.py`
- `tests/nec4_endcap_test.c`

No additional system software or library is required.

## Manual acceptance

From the repository root:

```sh
make tests/nec4_endcap_test
python3 tests/nec4-similarity/m5_endcap_oracle.py
python3 tests/nec4-similarity/m5_models.py
```

Run the published stress geometry directly:

```sh
./nec2c --wire-engine nec4ish --model-check warn \
  -i tests/nec4-similarity/m5-cap-monopole.nec \
  -o /tmp/m5-cap-monopole.out
```

The current table should be smooth near the source and should retain only a
small residual ripple at the free tip. For a loaded/networked case:

```sh
./nec2c --wire-engine nec4ish --jobs 8 --threads 1 \
  -i tests/performance/m6-loaded-network-402.nec \
  -o /tmp/m5-loaded-network.out
```

The power table should report approximately `8.9541e-3 W` input,
`5.2305e-4 W` radiated, `5.8411e-5 W` structure loss, and `8.3727e-3 W`
network loss.
