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
deliberately gated until M3 rather than silently returning a bare-wire answer.
