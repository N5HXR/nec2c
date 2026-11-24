# NEC4 Similarity M11B Results: Adaptive Total-E Tables

## Outcome

M11B is complete and ready for user acceptance. It adds an isolated adaptive
total-electric-field representation; production routing is deliberately
unchanged until M11E.

Every interpolated value accepted by the new representation met the frozen
`1e-4` difficult-domain gate at independent Halton witnesses. The worst result
was `5.140e-5`, in the seawater interface-line case.

For average soil, low-loss water, and seawater, all 12 material/region cases
accelerated every independent witness with no direct fallback. Adaptive table
sizes were small:

- switch and axis volumes usually required one cubic cell;
- grazing sheets required one to five cells;
- interface lines required nine to thirteen cells.

Lossless water exposed the expected total-field interference limitation:

- the interface line converged completely in 515 one-dimensional cells;
- the small-rho/large-z volume converged in one cell;
- the 3-wavelength switch volume reached its 512-cell budget and rejected
  `89.8%` of independent witnesses;
- the grazing sheet reached its 512-cell budget and rejected `90.6%`.

Rejection is safe and explicit: unresolved leaves return no interpolated
value. The future runtime must use the direct contour for those points. M11C
will investigate separated interfering contributions rather than increasing
the total-field cell budget without bound.

## Implementation

`nec4_ground_adaptive.c` implements deterministic hierarchical tensor-cubic
cells:

- fit nodes are `0`, `1/3`, `2/3`, and `1` on each active axis;
- error witnesses are excluded from the fit at `1/6`, `1/2`, and `5/6`;
- cells refine only when withheld scaled error exceeds `5e-5`;
- zero-width dimensions become exact lower-dimensional tables without
  duplicate samples;
- one-worker and four-worker builds produce identical table statistics and
  bit-identical evaluations;
- depth and cell budgets prevent unbounded preprocessing;
- unresolved or over-budget leaves force direct fallback.

`nec4_ground_adaptive_e.c` stores total E after removing the accepted F2 ray
phase and inverse-square spreading factor. Exact-interface samples use the
accepted Appendix-H interface limit. Existing near-region quasistatic tables
remain authoritative and unchanged.

The persistent cache is atomic, versioned, checksummed, and keyed by:

- adaptive algorithm and representation versions;
- frequency and both complex wave numbers;
- direct-contour panel, quadrature, and cutoff settings;
- bounds, tolerance, depth limit, cell budget, and tree payload.

Identity mismatch, trailing corruption, and checksum mismatch are rejected.

## Retained qualification

The retained matrix contains four materials and four difficult regions, with
128 independently located points per combination:

- average soil at 7.1 MHz;
- low-loss water at 5 MHz;
- seawater at 5 MHz;
- lossless water at one-metre wavelength;
- the current 3-wavelength switch volume;
- the exact interface line;
- a low-angle grazing sheet;
- the small-rho/large-z volume.

Evidence:

- `tests/nec4-similarity/m11b-adaptive.json`
- `tests/nec4-similarity/m11b-adaptive.png`
- `tests/nec4-similarity/m11b_adaptive.py`

## Acceptance commands

```sh
make -j16 tests/nec4_ground_adaptive_test tools/nec4_m11b_qualify
./tests/nec4_ground_adaptive_test
python3 tests/nec4-similarity/m11b_adaptive.py
```

Expected output includes:

```text
NEC4 adaptive cubic cells passed
M11B retained adaptive total-E accuracy and fallback gates passed
  worst accepted independent error: 5.140e-05
  lossless fallback regions: switch-volume=89.8%, grazing-sheet=90.6%
```

Regenerating the full retained data takes roughly three minutes on the current
32-thread host:

```sh
python3 tests/nec4-similarity/m11b_adaptive.py --collect --workers 16
```

## M11C entry gate

M11C may begin after user acceptance of:

1. the adaptive cubic representation and `1e-4` external accuracy gate;
2. the deterministic direct-fallback behavior at unresolved leaves;
3. a bounded build rather than unbounded refinement;
4. treating the high lossless switch/grazing rejection rate as the measured
   reason to develop an interference-safe representation next.
