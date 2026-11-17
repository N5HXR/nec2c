# NEC4 similarity M8 execution specification

Status: accepted by the user on 2026-11-17. All equation, contour,
singularity, pairing, limit, interface, and artifact gates pass.

## Purpose and boundary

M8 will build a slow, high-accuracy, production-independent electric-field
oracle for a point electric-current element near a planar interface. It will
cover every source/observer pairing:

- lower source, upper observer (`12`);
- upper source, lower observer (`21`);
- lower source and observer (`11`); and
- upper source and observer (`22`).

The oracle is validation infrastructure for the accelerated M9 ground engine.
It is not a matrix-fill implementation, is not selected automatically, and
must not make `--ground-engine nec4` perform a production solve in M8.

The superscript `H` in the theory guide denotes a *horizontal electric-dipole
source*. It does not denote magnetic field. The published equations used here
provide the complete electric-field response for vertical and horizontal
electric-current elements. A magnetic-field oracle is outside M8 unless a
public, independently auditable derivation is added.

## Public equation contract

The implementation will transcribe and test these parts of the NEC-4 Theory
Guide independently of the production C solver:

- Section 5.1.1, equations (5-1) through (5-7): cylindrical electric-field
  components and assembly for an arbitrarily oriented source;
- Section 5.1.2, Figures 5-2 and 5-3: Bessel- and Hankel-form integration
  contours, branch cuts, contour-selection rule, and adaptive quadrature;
- Appendix F: explicit field equations for all four medium pairings;
- Appendix G, equations (G-1) through (G-20): reciprocity, boundary-condition,
  direct-field, and image-field conversions; and
- Appendix H: singularity subtraction and cancellation-resistant integrands.

Every implemented formula will carry its equation identifier in source comments
and in the transcription tests. No constants or switching rules will be inferred
from NEC binaries or licensed source.

## Deliverable

The primary deliverable will be a Python reference tool using arbitrary-precision
`mpmath`, with JSON output containing:

- material and frequency inputs;
- source and observer coordinates and medium pairing;
- the cylindrical canonical components and Cartesian 3x3 electric dyadic;
- working precision, contour, subdivisions/evaluations, and termination reason;
- absolute and relative error estimates;
- results of independent refinement/representation comparisons; and
- explicit diagnostics for unsupported or unconverged points.

Python is intentional: a separate language, numerical library, and data path
reduce common-mode failures with the production C implementation. Sample points
may be distributed among processes, but each integral remains deterministic and
serial within one process.

## Internal implementation gates

These are engineering gates inside M8, not new user-acceptance milestones.

### Gate A - equation and convention lock

- Fix the `exp(-j k R)` time convention, complex square-root sheets, medium
  numbering, cylindrical axes, and current-moment normalization.
- Transcribe the `12` Bessel-form integrands for the five canonical electric
  components.
- Test algebraic identities among derivative and expanded forms before relying
  on quadrature.
- Reject ambiguous points exactly on the interface until one-sided limits are
  requested explicitly.

Exit condition: analytic identities and hand-audited witness values agree at
least to 30 decimal digits under arbitrary precision.

### Gate B - independent direct integration paths

- Implement the Figure 5-2 Bessel contour for
  `(|z| + |z'|) > 0.5 rho`.
- Implement the Figure 5-3 Hankel contour for the complementary region,
  including explicit branch-cut treatment and the `lambda = k1 + j s^2`
  substitution where required.
- Add Appendix H singularity subtraction and large-`lambda`
  cancellation-resistant evaluation.
- Provide a deliberately slower independent cross-check using a different
  integral representation and independently refined contour parameters.

Exit condition: both representations agree across the nonsingular qualification
grid within the measurement contract below, without sharing contour segmentation
or tail termination decisions.

### Gate C - all medium pairings and arbitrary orientation

- Use Appendix G transformations as the primary production-independent path
  from the canonical `12` evaluator to `21`, `11`, and `22`.
- Independently evaluate the Appendix F formulas at selected points so the
  transforms are not self-validating.
- Assemble the cylindrical components into the Cartesian 3x3 electric dyadic
  using equation (5-1), including the `rho = 0` symmetry limits.

Exit condition: reciprocity, interface conditions, and direct Appendix F versus
Appendix G comparisons pass across every pairing and source orientation.

### Gate D - qualification artifacts and diagnostic interface

- Generate convergence and residual maps over the complete grid.
- Save machine-readable JSON plus concise CSV summaries and static figures.
- Add a diagnostic point/sweep CLI; do not connect it to matrix fill.
- Document failures honestly rather than returning a value with an optimistic
  error estimate.

Exit condition: all automated gates pass and the manual commands below are
ready for user acceptance.

## Measurement contract

### Numerical convergence

For non-null components:

- default target: relative difference at most `1e-9` between independent
  representations and between the final two precision/refinement levels;
- preferred result: `1e-10` or better; and
- no value is accepted solely on the quadrature routine's internal estimate.

For exact or near-null components, use a scale-aware absolute residual of at
most `1e-11` times the largest non-null component for the same source/observer
configuration. Reports will show both absolute and relative metrics so field
nulls cannot create meaningless percentage errors.

The normal reference evaluation will compare at least 60- and 90-decimal-digit
runs. Difficult points may escalate precision. A point that cannot meet the
contract is marked unconverged and retained in the report.

### Physics identities

- Electric dyadic reciprocity:
  `G_E(r, r') = transpose(G_E(r', r))` for reciprocal media.
- At a source-free interface, one-sided limits preserve tangential `E` and
  normal `D` to `1e-8` relative or the null-aware absolute criterion.
- Exchanging identical upper and lower media reproduces the homogeneous
  finite-medium electric dyadic to `1e-10` relative.
- Increasing lower-medium conductivity toward the PEC limit reproduces the
  appropriate image solution above the interface and suppresses transmitted
  field; convergence of the limit is reported rather than judged at one
  arbitrary conductivity.
- Far-field results converge toward the plane-wave Fresnel transmitted or
  reflected ray, with the residual plotted versus distance.
- Quasistatic limits are checked only where a documented analytic expression
  applies; they do not substitute for full-wave validation elsewhere.

## Qualification grid

The grid is dimensionless where possible and deliberately includes the regions
called difficult by the theory guide.

- media: identical media, lossless `epsilon_r=4`, low-loss dielectric,
  representative dry/wet soil, and high-conductivity seawater/PEC approach;
- radial separation: `k0*rho` from approximately `1e-5` through `200`, with
  denser sampling near zero, one wavelength, grazing transitions, and the
  documented numerical range boundaries;
- source and observer heights: zero approached from either side, shallow,
  wavelength-scale, and deep points;
- geometry: `rho=0`, near-coincident points, grazing paths, far-ray paths, and
  cancellation-dominated combinations;
- orientation: vertical, radial-horizontal, azimuthal-horizontal, and three
  oblique Cartesian source directions; and
- pairings: `12`, `21`, `11`, and `22` at every applicable geometry.

The default acceptance run should finish in roughly 30 minutes on the current
machine. A comprehensive high-precision run may use all cores for independent
points but has a two-hour wall-clock cap. Exceeding the cap creates a documented
coverage gap rather than silently relaxing tolerance.

## M8 user acceptance

The acceptance package will provide:

1. field-component maps for magnitude, phase, convergence error, reciprocity
   residual, and interface-condition residual;
2. a table of worst points, including every unconverged point;
3. one command for a selected direct point in each medium pairing;
4. one command for a short line sweep crossing the interface; and
5. a concise explanation of any region deferred to M9 or later work.

M8 passes only after the automated measurement contract passes and the user
reviews the maps, runs the point/sweep examples, and returns `PASS`, `FAIL`, or
`PASS WITH CHANGES`.

## Explicit non-goals

- no production matrix-fill integration;
- no accelerated lookup tables or asymptotic switching (M9);
- no end-to-end buried antenna solve (M10);
- no claim of magnetic-field validation from horizontal-source superscripts;
- no NEC4 binary, licensed source, or reverse-engineered tuning data; and
- no automatic fallback to the classic NEC2 ground field under the NEC4 name.
