# Coated-wire M9 integration and release results

M9 packages the coated-wire work for review without changing the numerical
operator accepted in M3--M8.  It certifies the legacy path, records the cost of
the new operator, supplies retained end-to-end examples, and defines the
handoff contract for a future Antenna Turtle implementation.

## Compatibility and integration

The immutable M0 bare-wire witness remains byte-for-byte equivalent after its
documented provenance and timing normalization.  A paired 401-segment bare and
coated article confirms that both jobs have exactly 401 wire unknowns, no patch
or volume unknowns, and the same traditional report matrix order.  The coated
result has one coating domain; the bare result acquires no coated-wire metadata
or material evaluations.  Built-in and LAPACK back ends return the same coated
result shape and ordering.

No existing NEC card changed meaning.  `LD 7` is the additive coated-wire
designation defined in M2, and old decks need no edits.  The existing
`nec2c-results` schema remains version 1.  Coated jobs add optional material and
validity metadata while uncoated jobs retain their prior shape.

Implementing the Antenna Turtle DSL or exporter is deliberately out of scope
for this repository.  `COATED-WIRE-ANTENNA-TURTLE-SPEC.md` proposes a lexical
`with-wire-coating` form, whole-edge/tag semantics, frequency-dependent
material handling, strict junction defaults, and `LD 7` export behavior for
that separate project.

## Performance characterization

The retained benchmark interleaves bare and coated straight-dipole runs, uses
the automatic dense backend with 16 threads, and records medians after three
warmups and 15 measured repetitions.

| Segments | Bare wall (s) | Coated wall (s) | Fill ratio | Factor ratio | Wall ratio | RSS ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 95 | 0.00187 | 0.01738 | 8.96 | 1.00 | 9.28 | 1.00 |
| 401 | 0.00954 | 0.02840 | 4.39 | 1.02 | 2.98 | 1.05 |
| 801 | 0.02827 | 0.10010 | 5.04 | 1.01 | 3.54 | 0.99 |
| 1601 | 0.11872 | 0.32908 | 4.64 | 1.00 | 2.77 | 1.00 |

The 95-segment ratio is dominated by fixed process and timing costs and is
retained as a small-case guard, not as a capacity-planning number.  At 401--1601
segments, the adaptive two-radius coating quadrature makes matrix fill roughly
4--5 times as expensive and total execution roughly 2.8--3.5 times as
expensive.  Factorization is effectively unchanged and peak memory has no
material increase because coating adds no unknowns.  The benchmark is a
single-frequency operator characterization; frequency sweeps repeat the fill
cost for each frequency in the usual way and can still use the project's
frequency-level parallelism.

## End-to-end acceptance articles

`m9-vswr.nec` is a lossy coated dipole swept from 250 through 320 MHz.  Its
minimum 50-ohm VSWR is 1.284944 at 275 MHz, with
`63.493607 - j4.304019` ohms.  Its largest validity metric is 0.047988, inside
the M8 reference envelope.

`m9-pattern.nec` is the same promoted straight topology at 276.5 MHz.  The
retained elevation cut has its maximum at 90 degrees and a maximum mirrored
sample difference of `1.016e-9 dB` above -40 dB.  Its validity metric is
0.041465, also inside the reference envelope.

The release pattern intentionally does not use the disconnected coated Yagi.
M7 studied that geometry, but disconnected multi-element coating behavior
remains experimental.  The release witness uses only a topology promoted by
the accepted validation program.

## Supported envelope at release

- A coating is concentric, in contact with the conductor, constant over every
  tagged segment, and electrically thin according to the M8 metric.
- Lossless relative permittivity and frequency-adjusted loss tangent are
  supported through `LD 7`; every frequency is evaluated independently.
- Straight wires and the M7-qualified gradual bend construction are supported.
- Identically coated junctions are supported.  Coated/bare and unequal-coating
  junctions require the experimental policy and remain outside promoted use.
- Values above `sqrt(ER)*k0*b = 0.12` warn and continue; they are outside the
  validated reduced-operator envelope.
- Radomes, air gaps, partial segment coverage, arbitrary thick sleeves,
  endpoint-cap selection, and multilayer coatings remain future work.

The promoted evidence bundle comprises the M1 cylindrical/Cartesian OpenEMS
reference and audits, M3--M5 straight and feed validation, M6 junction study,
M7 Cartesian bend comparisons, and M8 thickness, endpoint, dogleg, staged bend,
and 45-degree-radial ground-plane studies.  `m9_release.py` verifies that those
versioned records remain present.

## Reproduction and acceptance

From the repository root:

```sh
python3 tests/coated-wire/m0_baseline.py
python3 tests/coated-wire/m9_release.py
python3 tools/coated_wire_m9_acceptance.py \
  --output /tmp/m9-acceptance.json \
  --figure-directory /tmp/m9-figures
```

The first two commands must pass.  The final command regenerates the VSWR and
pattern summaries and figures for human review.  The committed reference
figures are `tests/coated-wire/m9-vswr.png` and
`tests/coated-wire/m9-pattern.png`.

M9 was accepted on 2026-11-07 after the user reviewed these results, the
figures, and the Antenna Turtle handoff contract.  This completes the scoped
coated-wire project; commit and merge remain separate repository actions.
