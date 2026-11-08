# Coated-wire evidence

This directory holds the versioned evidence for the coated-wire project.  M0
contains only a bare-wire reference; it does not exercise a coating operator.

`m0-bare-dipole.nec` is the immutable compatibility witness.  It deliberately
uses only established NEC cards and requests three frequencies plus an ordered
far-field grid.  `m0-baseline.json` records semantic hashes and selected
observables from the unmodified solver at the M0 reference commit.  Timing is
excluded.

Run the check from the repository root:

```sh
python3 tests/coated-wire/m0_baseline.py
```

The check invokes the built-in dense backend with one job and one thread, then
compares the complete machine-result document after removing provenance fields
that necessarily depend on the checkout and output location.  It also verifies
the traditional report after normalizing timing lines.  A mismatch is a
compatibility failure; do not refresh the golden file merely because a later
implementation differs.

To inspect a fresh result manually:

```sh
./nec2c --backend builtin --jobs 1 --threads 1 \
  --result-json /tmp/coated-wire-m0.json \
  -i tests/coated-wire/m0-bare-dipole.nec \
  -o /tmp/coated-wire-m0.out
```

Run the dependency-free M1 harness checks with:

```sh
python3 tests/coated-wire/m1_harness.py
```

The full M1 simulations require the separately installed OpenEMS Python
environment. The milestone report records exact commands and compact versioned
results; the multi-gigabyte transient field dumps are deliberately not stored
in this repository. See [the M1 OpenEMS qualification report](../../docs/COATED-WIRE-M1-OPENEMS.md)
for the geometry convention, convergence criteria, and manual acceptance run.

The `m1-openems-caps*` artifacts contain the endpoint and mesh-refinement
results. The `m1-openems-audits*` artifacts contain the feed, gap, and transient
decay checks. JSON is the machine-readable record; PNG files are review plots.

M2 adds parser and topology validation without enabling a numerical coated-wire
solve. Run its complete acceptance test with:

```sh
python3 tests/coated-wire/m2_parser.py
```

For manual inspection, `m2-valid.nec` resolves one whole 11-segment tag:

```sh
./nec2c --validate-only --result-json /tmp/coated-m2.json \
  -i tests/coated-wire/m2-valid.nec -o /tmp/coated-m2.out
python3 -m json.tool /tmp/coated-m2.json
```

`m2-junction-mismatch.nec` deliberately joins two unequal coatings. The same
command should exit nonzero and emit `COATED_WIRE_JUNCTION_MISMATCH` in the
machine-readable diagnostic. Ordinary execution of either coated deck is
accepted only when it fits the M3 lossless isolated straight-wire envelope.
Lossy or unsupported geometry is rejected rather than silently returning a
bare-wire answer.

Run the complete M3 numerical acceptance test with:

```sh
python3 tests/coated-wire/m3_operator.py
```

The supplied manual pair and resonance sweeps are:

```sh
./nec2c --result-json /tmp/m3-bare.json \
  -i tests/coated-wire/m3-bare-dipole.nec -o /tmp/m3-bare.out
./nec2c --result-json /tmp/m3-coated.json \
  -i tests/coated-wire/m3-coated-dipole.nec -o /tmp/m3-coated.out

./nec2c --result-json /tmp/m3-bare-sweep.json \
  -i tests/coated-wire/m3-bare-sweep.nec -o /tmp/m3-bare-sweep.out
./nec2c --result-json /tmp/m3-coated-sweep.json \
  -i tests/coated-wire/m3-coated-sweep.nec -o /tmp/m3-coated-sweep.out
```

Regenerate the compact acceptance record and overlays with:

```sh
python3 tools/coated_wire_m3_acceptance.py \
  --output-dir tests/coated-wire
```

See `docs/COATED-WIRE-M3-RESULTS.md` for the numerical comparison and M3
scope boundary.

Run the M4 loss, material-evaluation, and power checks with:

```sh
python3 tests/coated-wire/m4_loss.py
```

The manual low-loss and lossy 50-ohm VSWR sweeps are:

```sh
./nec2c --result-json /tmp/m4-low-loss.json \
  -i tests/coated-wire/m4-low-loss-sweep.nec -o /tmp/m4-low-loss.out
./nec2c --result-json /tmp/m4-lossy.json \
  -i tests/coated-wire/m4-lossy-sweep.nec -o /tmp/m4-lossy.out
```

Regenerate the M4 acceptance record and figures with:

```sh
python3 tools/coated_wire_m4_acceptance.py \
  --output-dir tests/coated-wire
```

See `docs/COATED-WIRE-M4-RESULTS.md` for the per-frequency material semantics,
power derivation, OpenEMS comparison, and acceptance metrics.

Run the M5 feed-model decision checks with:

```sh
python3 tests/coated-wire/m5_feed.py
```

The manual witness uses the certified continuous center-segment source:

```sh
./nec2c --backend builtin --jobs 1 --threads 1 \
  --result-json /tmp/m5-continuous.json \
  -i tests/coated-wire/m5-continuous-source.nec \
  -o /tmp/m5-continuous.out
```

See `docs/COATED-WIRE-M5-RESULTS.md` for the 1--2 mm physical-gap convergence,
zero-distance reference-plane decision, current and pattern overlays, and the
scope boundary carried into M6.

Run the M6 junction-policy checks with:

```sh
python3 tests/coated-wire/m6_junctions.py
```

The production witness splits one physically uniform coated wire across two
tags. Both tags must designate exactly equal conductor radius, coating radius,
permittivity, loss tangent, and conductivity:

```sh
./nec2c --result-json /tmp/m6-uniform.json \
  -i tests/coated-wire/m6-uniform-junction.nec \
  -o /tmp/m6-uniform.out
```

Coated/bare and unequal-coating contacts remain errors in the default `strict`
policy. The research-only candidate can be exercised explicitly; its result is
stamped `experimental-junction` and includes a
`COATED_WIRE_EXPERIMENTAL_JUNCTION` warning:

```sh
./nec2c --coated-wire-junctions experimental \
  --result-json /tmp/m6-experimental.json \
  -i tests/coated-wire/m6-coated-bare-experimental.nec \
  -o /tmp/m6-experimental.out
```

This escape hatch is not a claim of validity. See
`docs/COATED-WIRE-M6-RESULTS.md` for the reciprocity evidence, asymmetric
OpenEMS references, supported-topology matrix, and M6 acceptance decision.

## M7 general geometry and mutual coupling

Run the automated geometry, null-limit, reciprocity, power, backend and
retained OpenEMS evidence checks with:

```sh
python3 tests/coated-wire/m7_geometry.py
```

The experimental manual witness is a uniformly coated, three-element Yagi. Its
elements are separate straight wires, so it exercises mutual coupling without
introducing a coating junction or sharp bend. The bare Cartesian reference
failed the 1% convergence gate, so the article requires explicit opt-in:

```sh
./nec2c --coated-wire-junctions experimental \
  --result-json /tmp/m7-yagi.json \
  -i tests/coated-wire/m7-coated-yagi.nec \
  -o /tmp/m7-yagi.out
```

Non-collinear connected coated wires inside the qualified bend envelope are
supported by default.
The initial coarse axis-aligned reference was outside the M7 budget, but later
fine OpenEMS and MATLAB references corrected that result and the bend addendum
qualified a conservative `X=45 degrees, Y=4 coating radii` envelope. Strict
mode rejects bends outside that envelope with `COATED_WIRE_BEND_UNVALIDATED`.
Run the accepted witness with:

```sh
./nec2c --result-json /tmp/m7-qualified-bend.json \
  -i tests/coated-wire/m7-qualified-bend.nec \
  -o /tmp/m7-qualified-bend.out
```

The accepted witness needs no experimental flag. Its
`coated_wires.operator_status` is `lossless-qualified-bend`, and its
`bend_policy` object reports the observed turns and chord/radius ratio. The
older sharp 90-degree witness remains research-only:

```sh
./nec2c --coated-wire-junctions experimental \
  --result-json /tmp/m7-bend.json \
  -i tests/coated-wire/m7-coated-v.nec \
  -o /tmp/m7-bend.out
```

See `docs/COATED-WIRE-M7-RESULTS.md` for the independent Cartesian OpenEMS
comparison, convergence record, topology support matrix and acceptance
decision.

## M8 thickness/permittivity envelope

Run the retained validity, loss, mesh-refinement, endpoint and warning-policy
checks with:

```sh
python3 tests/coated-wire/m8_validity.py
```

M8 uses `sqrt(ER)*k0*b`, where `b` is coating outer radius, rather than a fixed
thickness/conductor-radius ratio.  NEC2C warns and continues above 0.12.  The
largest tested point meeting 1% in both resonance and resistance is 0.1091;
the first failing point is 0.1276.  Each JSON material evaluation records the
metric and status, and the job records the warning policy.

The convention is flush dielectric and metal at a free wire end.  A thick,
extended-flat printed cap shifted the reference resonance by 2.178%, which is
reported as physical endpoint sensitivity rather than reduced-operator error.
The straight-wire envelope was also transferred to a 45-degree dogleg at
`b=6` and `12 mm`, a 90-degree route staged through two 45-degree transitions,
and an experimental ground plane with four 45-degree down-sloped radials.
Paired bare/coated
Cartesian OpenEMS references isolate coating error from ordinary geometry and
feed error.  The worst refined coating-effect errors were 0.229 percentage
point in resonance shift, 0.855 ohm in resistance change, 0.226 dB in pattern
change, and 0.0018 in normalized-current change.  The 45-degree ground-plane
reference has a roughly plus or minus 0.5-point resonance-shift uncertainty
from non-monotonic 2, 1.5, and 1 mm Cartesian-grid results.  The ground-plane
junction remains experimental; this evidence does not override the M6
topology policy.
See `docs/COATED-WIRE-M8-RESULTS.md` for the full table and reproduction notes.

For manual acceptance, compare the quiet in-envelope article with the
warn-and-continue article:

```sh
./nec2c --result-json /tmp/m8-within.json \
  -i tests/coated-wire/m8-within-envelope.nec \
  -o /tmp/m8-within.out
./nec2c --result-json /tmp/m8-warning.json \
  -i tests/coated-wire/m8-warning-envelope.nec \
  -o /tmp/m8-warning.out
```

The first result should contain `within-reference-envelope` and no M8 warning.
The second should complete normally, contain
`warning-outside-reference-envelope`, and emit exactly one
`COATED_WIRE_VALIDITY_ENVELOPE_EXCEEDED` diagnostic.

## M9 integration and release

Run the final compatibility, matrix-order, backend-shape, performance-record,
promoted-evidence, VSWR and pattern gates with:

```sh
python3 tests/coated-wire/m0_baseline.py
python3 tests/coated-wire/m9_release.py
```

Regenerate the compact end-to-end acceptance record and figures with:

```sh
python3 tools/coated_wire_m9_acceptance.py \
  --output /tmp/m9-acceptance.json \
  --figure-directory /tmp/m9-figures
```

The release VSWR witness has a 1.284944 minimum at 275 MHz.  The release
pattern is a promoted straight coated dipole, not the experimental disconnected
multi-element geometry, and is mirror-symmetric to `1.016e-9 dB` above -40 dB.
The retained performance record confirms equal bare/coated matrix order, no
material memory increase, and approximately 2.8--3.5 times total runtime for
the meaningful 401--1601-segment single-frequency cases.  See
`docs/COATED-WIRE-M9-RESULTS.md` and
`docs/COATED-WIRE-ANTENNA-TURTLE-SPEC.md` for the release decision and the
separate-project integration handoff.

### M7 bend-envelope addendum

The practical bend study keeps the centerline length fixed at 474 mm while
varying the number of chords per 90-degree corner and the bend radius relative
to the 4 mm coating radius.  Regenerate all 31 NEC cases and the representative
MATLAB input with:

```sh
python3 tools/coated_wire_m7_bend_study.py \
  --run-nec --write-matlab-config
matlab -batch \
  "run('/home/josh/homebrew/nec2c/tools/coated_wire_bend_matlab.m')"
```

The MATLAB command is deliberately stable; all runtime choices are read from
`tests/coated-wire/m7-bend-matlab-config.json`.  This permits a narrowly scoped
persistent execution approval without authorizing arbitrary MATLAB commands.

OpenEMS pilots and selected refinements use the same generated centerlines:

```sh
/home/josh/opt/openEMS/venv/bin/python3 \
  tools/coated_wire_m7_bend_openems_batch.py \
  --stage pilot --material bare --mesh-mm 2 --threads 16
```

Rebuild the compact comparison record and figure, then check the retained
evidence with:

```sh
python3 tools/coated_wire_m7_bend_analyze.py
python3 tests/coated-wire/m7_bend_envelope.py
```

The retained OpenEMS record includes sharp-bend convergence, a 61.85-million-
cell bare reference for the one-chord `R/b=4` candidate, and matched 2 mm and
1 mm bare/coated pairs. The engineering bands in the analysis JSON distinguish
high-confidence, "basically right," and qualitative agreement. They are
reporting bands, not automatic promotion gates; the candidate met the
high-confidence resonance, resistance, and coating-shift bands, but support
is enforced by the default topology validator.
