# NEC4 Similarity M11C Results: Interference-Safe E Tables

## Outcome

M11C is complete and ready for user acceptance.  It closes both lossless-water
coverage holes measured in M11B without weakening the frozen accuracy target:

| Region | M11B 512-cell fallback | M11C cells | M11C fallback | Worst independent scaled E error |
|---|---:|---:|---:|---:|
| grazing sheet, 3-5 wavelengths | 233/256 | 2,341 | 0/256 | `4.086e-5` |
| three-dimensional switch volume | 227/256 | 2,115 | 0/256 | `4.741e-5` |

The cell gate remains `5e-5`, below the external difficult-domain gate of
`1e-4`.  Each result uses 256 deterministic Halton witnesses excluded from the
fit.  The tables stopped naturally below the hard 4,096-cell budget.

These are genuine destructive-interference witnesses.  An exact partition of
the direct contour measured maximum individual-contribution-to-total ratios of
`116.3` on the grazing sheet and `279.0` in the switch volume.

Production routing is still unchanged.  M11E will choose an accepted ordinary
table first, then an accepted directional table, and finally the direct contour
if neither representation passes.

## Published-mathematics boundary

Section 5.1.3 of the 1992 theory guide publishes the two-ray phase model and
the original NEC4 least-squares basis.  The 2011 GN3 report says that NEC4.2
sometimes interpolates two branch-cut integrals separately, but it does not
publish those integrals, their contour deformation, or selection thresholds.

M11C therefore did not invent GN3 formulas.  It tested three independently
specified candidates:

1. **Exact moving-contour partition.**  The Figure 5-3a Hankel contour was
   partitioned at `P2`.  Ordered recomposition agreed with the unsplit direct
   result within `2.206e-13`.  This is retained as a cancellation diagnostic.
   It is not a stable interpolation representation because `P1`, `P2`, and
   `P3` move with geometry: one contribution converged in five cells while the
   other exhausted 512 cells and accelerated none of the independent points.
2. **Published Appendix-J ray plus numerical remainder.**  The direct field
   was represented exactly as the accepted public asymptotic field plus a
   direct numerical remainder.  The remainder retained the difficult phase
   variation, exhausted 512 cells, and accelerated none of the grazing
   witnesses.
3. **Directional total-field cubic cells.**  The accepted M11B transformed
   total field and withheld error metric were retained.  A rejected cell is
   bisected only along the coordinate with the largest center-line withheld
   error.  This bounded k-d refinement passed both difficult regions.

The selected method follows M11's stop rule: public mathematics did not support
claiming GN3 branch-cut equivalence, so direct/adaptive total-field coverage was
extended instead.

## Implementation

`nec4_ground_direct_cross_hankel_split()` exposes the exact ordered `P2`
spectral partition.  Its type and comments explicitly state that it is a
diagnostic contour contribution, not NEC4.2 GN3 equivalence.

`nec4_ground_adaptive_build_anisotropic()` adds opt-in deterministic k-d
refinement:

- fit nodes remain `0`, `1/3`, `2/3`, and `1`;
- withheld witnesses remain `1/6`, `1/2`, and `5/6`;
- accuracy and direct-fallback behavior are unchanged;
- the split axis is chosen from center-line withheld errors;
- only one axis is bisected per rejected cell;
- one- and multi-worker sampling remains deterministic;
- the hard depth and cell budgets remain mandatory.

The M11B isotropic builder and its cache format remain unchanged.  Directional
cache persistence is deliberately rejected until M11E versions the production
representation and routing identity.

## Retained evidence

- `tests/nec4-similarity/m11c-interference.json`
- `tests/nec4-similarity/m11c-interference.png`
- `tests/nec4-similarity/m11c_interference.py`
- `tools/nec4_m11c_qualify.c`

Cold 16-worker builds on the current host took about 76 seconds for the grazing
sheet and 183 seconds for the switch volume.  These are table-construction
costs; cached startup and production matrix timing belong to M11E.

## Acceptance commands

```sh
make -j16 tests/nec4_ground_adaptive_test tests/nec4_ground_direct_test \
  tools/nec4_m11c_qualify
./tests/nec4_ground_adaptive_test
./tests/nec4_ground_direct_test
python3 tests/nec4-similarity/m11c_interference.py
```

Expected summary:

```text
NEC4 exact spectral split passed: recomposition 2.206e-13
M11C bounded directional interference tables passed
  grazing-sheet: 2341 cells, 4.086e-05 worst error, 116.3x cancellation
  switch-volume: 2115 cells, 4.741e-05 worst error, 279.0x cancellation
```

The full cold qualification can be regenerated with:

```sh
python3 tests/nec4-similarity/m11c_interference.py --collect --workers 16
```

It takes roughly five minutes on the current host.

## Acceptance decision

Acceptance means approving:

1. bounded directional total-field refinement as M11C's documented
   interference-safe representation;
2. retaining the exact contour split only as a diagnostic, with no GN3 claim;
3. zero tolerance relaxation and direct fallback for any unresolved cell;
4. deferring persistent directional caches and production routing to M11E;
5. proceeding to M11D direct magnetic-field derivation.
