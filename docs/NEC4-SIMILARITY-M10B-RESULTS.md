# NEC4 Similarity M10B Results

## Outcome

The one-sided non-contacting interface dipole now runs to completion in
NEC4-similarity mode.  Its matrix/input solution agrees closely with the
classic exact-ground control, and its 1000-wavelength interface observation
agrees with the published value while reliably raising the conservative
grazing warning.

## Published Table 1 control

Using a `1e-4 lambda` radius inferred from the report's ground-stake models and
31 segments, the classic exact-ground path produced these relative errors:

| Ground | Relative error |
| --- | ---: |
| `10-j1` | 15.248% |
| `10-j10` | 0.303% |
| `10-j100` | 4.164% |
| `10-j1000` | 0.618% |
| Perfect | 1.379% |

The first, lowest-loss case is the most sensitive to the report's omitted
reference geometry/discretization detail.  A joint radius/segmentation sweep
changed this value by more than ten percentage points without producing one
geometry-independent optimum.  The other four values are within 5 percent,
including the perfect-ground analytical limit.  The report's independently
quoted input impedance also does not identify the missing radius uniquely, so
the reference must not be tuned per material.

## NEC4-similarity result

- Input impedance relative delta from the classic solution: `4.557e-4`.
- One-worker versus automatic-worker result: byte-for-byte identical at the
  parsed execution level.
- NEC4 field ratio to the published reference: `0.96325` (3.675% low).
- Diagnostic count: exactly one `NEC4_M9_GRAZING_INTERFERENCE` warning.

The interface-height sweep from `z=1e-2 lambda` through `1e-8 lambda` and the
exact `z=0` limit is smooth.  From 100 through 1000 wavelengths the NEC4 result
remains 0.823--0.836 of the NEC2 result; at 1000 wavelengths NEC4 is closer to
the published value.  The warning remains because this geometry enters the
documented risk region, not because an inaccurate result is accepted.

## Defects found and fixed

1. M7 classified a physically valid single free endpoint on the interface as
   an unsupported junction.
2. Deferred NE/NH evaluation reused the last matrix source-segment index and
   interpreted `ij` as a matrix observer segment.  This could reject the field
   point, choose the wrong observer medium, or choose the wrong source medium.
3. A grazing interaction first encountered after matrix fill could miss the
   warning because warning publication happened too early.
4. While separating the point path, an initial implementation changed the
   classic near-field `ij` selector; the published control exposed it and the
   historical selector is now explicitly preserved.
5. The NEC4 arbitrary-point adapter used `legacy direct + (NEC4 complete -
   analytic homogeneous)`.  The two direct implementations differed only
   slightly, but their residual decayed like a space wave and overwhelmed the
   true surface field after hundreds of wavelengths.  Direct use of the
   complete NEC4 segment field restores the cancellation and changes the
   published-reference ratio from roughly 1824 to 0.963.

## Six-step investigation

1. Joint radius/segmentation sweeps showed that the report's omitted reference
   geometry explains the legacy low-loss sensitivity; no tuning was adopted.
2. The upper-side `z -> 0` sweep is finite and smooth.
3. Direct, reflected, saddle, and pole terms were separated.  The isolated
   Appendix-J saddle and pole terms are large but the complete antenna field
   is well conditioned when assembled without cross-implementation
   cancellation.
4. A 1--1000 wavelength distance sweep localized the former defect to the
   long-range point-field adapter, not the matrix solution.
5. Table 1's normalization was re-derived from the report as
   `|E_z|^2 / P_input`; the existing test uses exactly that quantity.
6. The original NEC4 second-order equations were retained.  No empirical
   tuning factor and no NEC4.2 replacement was introduced.
