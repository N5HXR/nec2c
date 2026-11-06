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
