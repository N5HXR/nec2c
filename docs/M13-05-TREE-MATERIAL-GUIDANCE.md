# M13-05: tree-material assumptions and milestone closure

Status: COMPLETE (2026-10-21). This is the final increment of M13
for the agreed scalar-material and GN 2 wires-first scope. It adds reporting
and a small native sensitivity witness; it does not introduce a new field model.
The independently refined OpenEMS evidence in M13-01 through M13-04 remains
applicable. Another OpenEMS tree run would compare two discretizations of the
same assumed material; it would not establish that the assumed tree material
represents a particular specimen.

## Material meaning and required provenance

The unknown is cell-average **volume polarization current**, not a dielectric
surface-current model. All current materials are passive, isotropic and
nonmagnetic. With exp(+i omega t),

    epsilon_relative = ER - i * (EI + sigma / (omega * epsilon0)).

Conductivity is in S/m, EI is dimensionless, and frequency is in Hz in this
formula. The resource syntax permits either EI or sigma, not both. If converting
measured EI to conductivity at a witness frequency, do not add a second estimate
for that same loss. Resource version 2
supports constants; version 3 supports linearly interpolated frequency tables
and rejects extrapolation. Constant sigma implies frequency-dependent imaginary
permittivity; constant EI is a different material law. Passivity alone does not
establish that an interpolated table is a complete causal broadband model.

A tree result must travel with the following information:

* geometry dimensions, coordinate units, trunk/branch/crown construction,
  unresolved structure and whether the crown represents a homogenized volume;
* each region's material identity, scalar approximation, law (ER/EI/sigma or
  table), units, frequency coverage and the origin of every value;
* whether values were measured, fitted, obtained from a named dataset or merely
  assumed, including available specimen/condition information and any known
  temperature, moisture or orientation limitations of the source;
* ground model and parameters, conductor/ground clearances, mesh, solver,
  residual and power accounting; and
* scenarios considered, numerical refinement evidence, and excluded effects.

Do not label an assumed ER/sigma pair as validated wood, bark or foliage.
Material tables allow users to supply measured frequency dependence; they do
not infer moisture, grain direction or species from a geometry label. Treating
all tree regions as the same scalar material is an explicit idealization.
A homogeneous crown is a filled effective dielectric volume, not a resolved
collection of leaves. No effective-medium mixing law is selected automatically.
Different material regions can share boundaries under the M13-03 policy.

## Sensitivity and uncertainty

`tools/m13_05_tree_report.py` attaches assumptions to numerical results and
reports ranges across supplied material scenarios. The input format is shown
in `tests/m13-05-tree-input.json`. It requires named scenarios with a common
frequency and mesh pitch, valid passive material values, finite impedances and
absorption, residual <=1e-8 and power defect <0.3%. Ground and geometry assumptions
are shared metadata; the producer must ensure all scenarios use the same actual
geometry, solver configuration and ground. The tool does not reconstruct or
validate a model from an impedance alone.

The report includes the nominal case, every supplied scenario and sampled
ranges of ER, sigma, impedance components and dielectric absorption. These are
**sampled scenario envelopes**. They are not confidence intervals or guaranteed
bounds between samples; no monotonic material-to-impedance response is assumed.
Without an evidence-based probability model, do not attach a confidence level.
The example uses one-at-a-time variations, not a joint ER/sigma uncertainty grid.
Interaction effects and unsampled extrema are therefore untested.

Use a separate refined solve with the nominal material to assess numerical
mesh sensitivity. The reporter rejects mixing different mesh pitches into the
material envelope and reports a supplied refinement separately. A small linear
residual establishes linear-system convergence. A small power defect checks
energy consistency. Neither bounds material/geometry error, and a two-mesh
comparison is not a certified continuum-error estimate. The established
max(2 ohm, 2% of each reference component) comparison rule belongs to benchmark
acceptance; it is not a universal tree-prediction error bar.

For a production result, report uncertainties in geometry, homogenization,
material provenance and ground assumptions even if no quantitative interval is
available. State unknown uncertainty as unknown. If credible material bounds are
available, sample them using the same mesh, frequency and ground, then refine
representative nominal or difficult cases separately. Report changes to antenna
impedance, fields or absorption as scenario sensitivity, with units and fixed
source normalization. The initial reporter covers driven impedance and dielectric
absorption; field-pattern uncertainty still needs its own declared sampling.

## Supported production paths and resources

Free-space and PEC tree models can use the regular-grid FFT path. Touching
scalar materials and PEC-ground contact are supported with existing constraints.
The historical 9,432-cell tree records in M9/M11 establish particular lossless
geometry/resource checkpoints. The accepted PEC-tree target is 40 minutes per
frequency on the constrained VM; these timings are observations for that fixture,
not predictions for arbitrary lossy/material-contrast scenarios. The present
small sensitivity witness is not a full-size lossy-tree resource certification.

GN 2 real ground supports wires and strictly positive dielectric/ground and
conductor clearances with dense or assembled GMRES. It has no FFT acceleration.
M13-04 validates small fixtures up to 64 active volume cells. Do not transfer the
PEC 9,432-cell timing/memory claim to real ground; estimate and measure the actual
model before committing a large sweep. Zero substrate conductivity with ER>1
still leaves an interface and transmitted power. Lower-half-space flux is
separate from dielectric absorption; only their appropriate sum enters the
combined driven energy balance.

GN 0, real-ground patches, burial/ground contact, screens/cliff, ground FFT,
conductor/dielectric contact or penetration, anisotropy and magnetic media are
outside this completed scope. M12 remains deferred. These exclusions are explicit
boundaries, not additional M13 acceptance gates.

## Reproduction

The synthetic witness uses a rectilinear trunk, two branches and crown over PEC,
a positive-gap 11-segment wire and 100 MHz excitation. Five material scenarios
use ER=2/3/4 and sigma=.005/.01/.02 S/m, with nominal ER=3/sigma=.01. These values
are arbitrary sensitivity witnesses, **not measured vegetation parameters**.
All scenarios share a 0.1 m cube mesh; the nominal refinement uses 0.05 m cubes.

    python3 tools/m13_05_tree_experiment.py --workdir /tmp/m13-05-native --output-dir /tmp/m13-05-evidence
    python3 tools/m13_05_tree_report.py tests/m13-05-tree-input.json --output /tmp/m13-05-report.json

The experiment writes the actual solver decks, material files and raw output
under its work directory. Retained reports include the solver binary hash.
`tests/m13_tree_report.py` checks envelope arithmetic, invalid metadata/material
and solver refusal, mesh separation and exact retained-report reproduction.

## M13 exit-criterion evidence

| Exit criterion | Retained evidence |
| --- | --- |
| Passive signs and dielectric/ground loss power | M13-01 analytical Joule checks; M13-04 substrate-aware independent flux and driven/plane-wave power |
| Zero-loss and appropriate ground limits | M13-01/03 zero-loss limits; M13-04 PEC, matched-background and quasistatic/Fresnel checks; M13-05 localized-source remote-ground E/H limit |
| Explicit real-ground subset and reciprocity | M13-04 GN 2 positive-gap wire scope, reciprocal ports and deterministic unsupported-mode refusals |
| Frequency tables and extrapolation refusal | M13-02 parser, interpolation, sweep range and independent witness tests |
| Layered/concentric material analytical cases | M13-03 independent interface and coated-sphere references |
| Independent lossy solver comparisons/refinement | M13-01/02/03/04 returned OpenEMS decay/material/mesh/domain audits and native component comparisons |
| Tree-material assumptions and uncertainty accompany results | This guidance, material-report tool, five-scenario synthetic tree and separate nominal refinement |

M13-05 closes the reporting criterion: all six native witness/refinement solves
and all 24 regression groups pass. M13 is COMPLETE for the agreed scope.
No additional empirical tree-material validation is claimed.

The final limit review adds an explicit localized-source remote-ground check to
`tests/vie_halfspace_test.c`: source and receiver keep their 0.1 m direct
separation while both heights increase from 0.5 m to 4 m over ER=4-i*.4 ground
at k0=.01 rad/m. The reflected electric and magnetic dyads must decrease by
factors greater than 100 and 40 respectively. This tests the low-frequency
remote-interface limit; it does not assert that plane-wave illumination or
far-field interference becomes identical to free space after translation.

## Final witness results

The five 156-cell material scenarios produce sampled impedance ranges
R=9.499956..9.894527 ohm and X=-824.553969..-824.228466 ohm at 100 MHz.
The nominal 1,248-cell refinement changes R/X by 0.016868/0.022102 ohm.
Largest power defect is 0.272187%, within the 0.3% gate; largest true residual
is 9.62e-11. These observations apply only to this synthetic geometry and
sampling plan. In particular, they do not bound realistic tree-material error.
The remote-ground electric/magnetic ratios are 0.00202262/0.0157424. All 24
regression groups pass, including reporting refusals and retained-report
reproduction. `tests/m13-05-acceptance.json` records closure of every M13 exit
criterion and the unchanged production solver binary hash.
