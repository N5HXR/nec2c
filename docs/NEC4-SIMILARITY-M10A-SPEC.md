# NEC4 similarity M10A normalization specification

Status: implementation complete on 2026-11-19; pending user acceptance.

## Purpose

M10A qualifies the absolute normalization of the complete wire/interface
operator before realistic buried-antenna results are accepted.  Agreement
between the M8 direct contour and M9 accelerated tables is necessary but not
sufficient because both paths enter the same matrix adapter.

The public NEC-4 Theory Guide fixes the conventions used here:

- time dependence is `exp(+j omega t)` and outgoing waves are `exp(-j k R)`;
- lower medium is medium 1 and upper medium is medium 2;
- equations (5-1) and (5-2) give electric field for unit current moment
  `I l = 1 A m`;
- equation (3-3) uses `A + B sin(k_s s) + C[cos(k_s s)-1]`;
- section 3.2 selects `k_s` as the wave number of the medium containing the
  wire in ordinary cases; and
- NEC's solved wire coefficient is current per metre in the internal
  free-space-wavelength normalization and is multiplied by `lambda0` when
  reported as amperes.

## Gate A - SI field to NEC matrix units

The M8/M9 dyadic returns `E/(I l)` for a unit electric-current element.  A
physical integration over a source segment therefore returns `E/I` in
ohm/metre.  NEC's matrix coefficient is in ohms because its unknown is current
per metre.  The sole conversion at this boundary is

`Z_matrix = lambda0 * (E/I)`.

The conversion must be explicit in the matrix adapter and must not be hidden
in either the direct or accelerated field evaluator.

Acceptance:

- Maxwell similarity: doubling frequency while halving every coordinate and
  radius doubles physical `E/I` and leaves `lambda0 E/I` invariant;
- identical-media corrections remain null for every medium pairing and all
  three basis terms; and
- classic NEC2 and free-space NEC4 paths do not change.

Gate A passed its unit-level test on 2026-11-19.  The audit found and corrected
a missing `lambda0` factor shared by both M8 and M9 production adapters.  On
the 14 MHz seawater whip witness the correction changed the previously
unqualified `14.585-j2554.222 ohm` result to `38.124-j17.879 ohm`.  This is a
plausibility improvement, not yet an accepted physical oracle.

The complete Gate A production check now includes both a wholly buried wire
and a split wire crossing an invisible vacuum/vacuum interface.  Feed
impedance and the complete current vector reduce to the ordinary free-space
NEC4 path within `2.3e-11` relative error.

## Gate B - medium-specific current basis and direct kernel

Every source segment must use the complex normalized basis wave number

`k_s lambda0 = k_medium lambda0`.

Upper-medium segments therefore retain `2 pi`.  Lower-medium segments use the
frequency-dependent complex wave number prepared by M7.  This value must be
used consistently by:

- `sbf`, `tbf`, and `trio` coefficient construction;
- current reconstruction and charge differentiation;
- constant, sine, and cosine source-field integration;
- free-end and voltage-gap cap parameters;
- matrix loading/source terms; and
- current, charge, power, and machine-result output.

No implementation may combine a lower-medium complex basis with a
free-space-basis direct kernel and call the difference an interface
correction.  Same-medium singular terms and finite-radius matching require an
explicitly qualified lower-medium kernel.

Acceptance:

- a wholly buried antenna in identical media reduces to a geometrically
  scaled free-space antenna using the medium wavelength and impedance;
- lossless `epsilon_r=4` gives the expected wavelength and intrinsic-impedance
  scaling;
- lossy homogeneous results converge under segmentation and satisfy
  reciprocity and passivity;
- an electrically similar pair of complete decks has invariant impedance when
  geometry is scaled inversely with frequency and conductivity is scaled with
  frequency to hold complex relative permittivity constant; and
- source and observer exchange agrees after the documented current-moment and
  junction transformations.

## Gate C - interface limits and production diagnostics

Cases:

- vacuum/vacuum and identical non-vacuum media;
- `epsilon_r -> 1`, `sigma -> 0`;
- lossless dielectric and high-conductivity limits;
- PEC approach for an upper wire;
- source and observer approaching the interface from either side;
- normal and oblique split-wire crossings; and
- quasistatic short-dipole and long-distance radiating limits.

Each case records impedance, current, matrix residual, reciprocity defect,
power balance, and sensitivity to segmentation and source placement.

## Gate D - seawater witness explanation

The 5 m whip/1 m stake is retained as a stress case, not as a truth source.
The final report will separate:

- the effect of the corrected SI/NEC matrix normalization;
- the effect of the lower-medium complex current basis;
- direct-versus-grid interpolation error;
- segmentation convergence in both media; and
- comparison with an independent analytic, published, MATLAB, or OpenEMS
  witness of declared uncertainty.

M10A is complete only when the final value is explained by these components,
not merely because it appears plausible.

Implementation results and the remaining oracle boundary are recorded in
`NEC4-SIMILARITY-M10A-RESULTS.md`.
