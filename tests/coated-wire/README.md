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
