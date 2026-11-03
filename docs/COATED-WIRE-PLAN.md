# Coated-wire project plan and M0 contract

Status: M0 and M1 accepted; M2 in progress.

Date: 2026-11-03  
Branch: `codex/coated-wire`  
Bare-wire reference commit: `60df8fc722fef936229a4a358ca70ed918462da0`

## Goal and boundary

Add a numerically validated, passive, concentric dielectric coating to the
existing thin-wire formulation without changing any uncoated NEC deck or its
results.  The first physical witness is a center-driven dipole.  The target is
insulated wire and solid printed housings in direct contact with the conductor.
A hollow radome or an air gap between conductor and housing is a later feature.

The coating is represented by a modified thin-wire operator, not by cuboidal
`DF` cells.  The intended formulation introduces no additional unknowns, so an
N-segment coated model retains the matrix order of the corresponding bare-wire
model.  This is a hypothesis to validate, not permission to hide disagreement
with a full-wave reference.

## Frozen user interface

The coating designation is the 4nec2-compatible `LD 7` extension:

```text
LD 7 TAG FIRST LAST ER OUTER_RADIUS TAN_DELTA SIGMA_DC
```

The fields are:

| Field | Meaning |
| --- | --- |
| `TAG` | Positive wire tag identifying one coating domain. |
| `FIRST`, `LAST` | Standard tag-relative segment selection. `0 0` means the complete tag. Partial coverage is rejected in the first release. |
| `ER` | Real relative permittivity, at least 1. |
| `OUTER_RADIUS` | Overall conductor-plus-coating radius in metres; it must be at least the selected `GW` conductor radius. |
| `TAN_DELTA` | Nonnegative dielectric loss tangent, optional and zero when omitted. |
| `SIGMA_DC` | Additional nonnegative bulk conductivity in S/m, optional and zero when omitted. |

The first two floating fields intentionally retain 4nec2's meaning and units.
The existing `GW` radius is the conductor radius.  User-facing APIs may call
`2*GW_RADIUS` the conductor diameter and `2*OUTER_RADIUS` the coating diameter;
the NEC card remains radius-based for compatibility.  Unlike 4nec2's conversion
to a distributed inductive load, nec2c will use the coating operator specified
below.

For the `exp(+i omega t)` convention, the material evaluated at each frequency
is

```text
epsilon_r_complex(f) = ER - i * (ER*TAN_DELTA + SIGMA_DC/(omega*epsilon_0)).
```

Thus a fixed loss tangent remains fixed over an `FR` sweep, while a fixed
conductivity contributes a loss tangent proportional to `1/f`.  The two loss
terms are additive.  Arbitrary tabulated or dispersive material laws are not in
the first card contract; they can later be added with an external resource
without reinterpreting these fields.

### Selection and topology rules

* `TAG` must be nonzero and must exist.
* The selected final segments must have one conductor radius.
* `0 0`, or an explicit range equal to the complete tag-relative range, is
  accepted.  Any subset is rejected as partial coverage.
* A tag is one coating domain.  A producer must split and retag geometry where
  coating parameters change.  Reusing a tag on unrelated `GW` primitives is
  discouraged because the post-geometry NEC representation does not retain
  enough provenance to distinguish them safely.
* Duplicate `LD 7` designations of the same tag are rejected rather than using
  order-dependent overrides.
* `ER=1`, zero loss and any legal outer radius must reduce to the bare-wire
  result within the milestone's explicit numerical tolerance.
* `OUTER_RADIUS=GW_RADIUS` must likewise reduce to the bare-wire result,
  including when `ER` is greater than one.

The initial release requires identical coating parameters at a physical wire
junction.  Unequal coated/coated and coated/bare contacts are rejected until
their junction basis and charge continuity are independently validated.  The
feed is treated separately: the validation program will compare both (a) the
ordinary NEC center-segment source on a continuously coated wire and (b) two
coated arms separated by a finite gap and driven by a differential/lumped port.
The gap will be swept and the port plane de-embedded.  The comparison will
determine which forms can be supported and what warning or constraint users
need.  A future asymmetric hot-coated/cold-bare experiment is permitted only as
an explicit validation case, not silently accepted production topology.

## Numerical operator specification

Let `a` be conductor radius, `b` coating outer radius, `epsilon_c` the complex
coating permittivity, and `I(s)` the existing wire-current expansion.  The
starting model follows the equivalent-polarization-current family of Richmond
and Newman and the two-potential specialization of Popovic, Djordjevic and
Kircanski:

```text
surface charge:       q_s(s) = -(1/(i*omega)) dI(s)/ds
coating polarization: P_r(r,s) proportional to
                      ((epsilon_c-epsilon_0)/epsilon_c) q_s(s)/(2*pi*r)
                      for a <= r <= b
```

After radial integration, the coating changes the charge/scalar-potential part
of each wire-wire matrix interaction through separate inner-radius (`a`) and
outer-radius (`b`) kernels.  It does not become an `LD 2` impedance, does not
alter the vector-potential current term as an empirical velocity factor, and
does not create independent volume unknowns.  The exact discrete correction,
self term and endpoint terms must be derived in M2 against the published
equations and the existing `sbf`/`trio` basis; no formula is to be inferred from
the prose in this document alone.

The likely implementation seam is `cmww`: it already evaluates each wire
source against every tested wire segment after `trio` expands the connected
basis.  Coating metadata must remain separate from legacy `zload`; `LD 7` must
not pass through the RLC conversion or matrix-load subtraction used for
`LD 0`--`LD 5`.  All excitation, current reporting, far-field and near-field
paths must be audited for the same physical current convention.  Radiation is
not to be modified merely to force agreement.

This starting operator has a known risk: later quasi-static analysis identifies
an axial coating field as well as the commonly assumed radial field.  Endpoint,
feed-gap, thick-coating and high-permittivity sweeps therefore define the
validated envelope rather than an optimistic thinness slogan.

## Diagnostics and machine results

Existing decks retain result schema version 1 and byte-for-byte numerical
content after timing/path normalization.  A coated job may add optional fields
only; existing required fields and enumerated solver names do not change.

Planned optional metadata includes the evaluated conductor radius, outer
radius, `ER`, `TAN_DELTA`, `SIGMA_DC`, and complex permittivity for every coated
tag and frequency.  Coating loss contributes to the existing
`power.dielectric_absorbed_w`; it must not be mislabeled as conductor or network
loss.

Invalid geometry or non-passive/nonfinite input is an error.  A model outside
the independently validated envelope emits a machine-readable warning and
continues, as requested.  The initial warning code is
`COATED_WIRE_OUTSIDE_VALIDATED_ENVELOPE`, with context containing the measured
dimensionless ratios and accepted envelope.  Users may run with scissors, but
the result will say where the scissors were used.

## Research basis and evidence hierarchy

1. J. H. Richmond and E. H. Newman, “Dielectric coated wire antennas,”
   *Radio Science* 11(1), 13–20 (1976),
   [doi:10.1029/RS011i001p00013](https://doi.org/10.1029/RS011i001p00013).
   This is the primary equivalent-polarization-current basis and explicitly
   retains the bare-wire matrix size.
2. B. D. Popovic, A. R. Djordjevic and N. M. Kircanski, “Simple method for
   analysis of dielectric-coated wire antennas,” *The Radio and Electronic
   Engineer* 51(3), 141–145 (1981),
   [open issue PDF](https://www.worldradiohistory.com/UK/British-Institution-of-Radio-Engineers/1981/TREE-1981-03.pdf).
   It supplies the two-radius, two-potential derivation and reports its strongest
   agreement for coatings no thicker than roughly one conductor radius and
   relative permittivity below 10.  That historical range is a starting test
   region, not the final product limit.
3. A. I. Mowete, A. Ogunsola and L. Sandrolini, “A quasi-static theory for
   dielectric-coated thin-wire antenna structures,” *PIER Letters* 20, 45–54
   (2011), [doi:10.2528/PIERL10121104](https://doi.org/10.2528/PIERL10121104).
   Its axial-field finding is treated as an explicit model-risk test.
4. OpenEMS FDTD is the independent numerical corroborator.  Analytical
   coaxial/capacitance and electrically thin limits are checked first; agreement
   between two numerical programs alone is not sufficient.

## Milestones and acceptance gates

Every milestone ends with a result summary, versioned artifacts, manual test
instructions and a user PASS/FAIL verdict.  Development stops at that gate.

### M0 — contract and immutable baseline (this milestone)

Goal: freeze scope, syntax, conventions, topology policy, evidence hierarchy
and a bare-wire reference before implementation.

Metrics:

* No production source or numerical behavior changes.
* Full existing test suite passes.
* The M0 bare dipole's canonical result and normalized traditional report match
  their committed SHA-256 values exactly.
* The feature branch points back to the recorded reference commit.

User acceptance: review this contract and run the baseline check.

### M1 — OpenEMS cylindrical reference qualification

Goal: establish an independently converged reference workflow before judging
the NEC operator.

Cases: bare dipole and lossless coated dipole, each with both a continuous
center excitation analogue and a two-arm finite-gap differential port where
the solvers permit it.  Sweep radial/axial mesh, PML distance, domain size,
time-step/end criteria, feed gap and port de-embedding.  Compare cylindrical
and Cartesian OpenEMS meshes and explicitly audit the cylindrical-axis
singularity.  The normative coated geometry has a solid, concentric dielectric
annulus with flat ends flush with the conductor ends and feed-gap faces.  An
endpoint sensitivity audit compares that convention with a flat dielectric
overhang of one outer coating radius, a hemispherical dielectric termination,
and bare conductor protruding beyond the dielectric.  These variants are
reference geometries, not additional NEC card syntax.

The flat variants must place their boundaries on mesh planes.  Because a
cylindrical Yee grid still staircases a hemisphere in the `rho,z` meridian,
qualify the hemispherical result with its own axial and radial shape-refinement
pair rather than inferring its uncertainty from either flat geometry.

Resolve the coating itself with at least four radial cells, analogous to the
through-thickness requirement for an FDTD PCB substrate, and confirm the
selected result with eight cells through the coating.  For later thick sleeves,
cell placement should follow the expected radial gradient rather than relying
only on a uniform whole-domain pitch.

Metrics: successive-refinement change below 0.25% in resonance frequency,
0.5 ohm in each feed-impedance component near resonance, and 0.05 dB in sampled
far-field magnitude.  If those targets cannot be reached, report the achieved
uncertainty and do not promote the reference.  Report resonance, impedance and
pattern deltas caused by each endpoint treatment separately from numerical
mesh uncertainty; retain flat/flush as the reduced model's convention unless
the audit demonstrates that an endpoint correction is necessary.

User acceptance: review convergence plots and run one supplied OpenEMS case.

### M2 — parser, data model and discrete derivation

Goal: implement `LD 7` parsing/validation and a reviewed derivation mapped to
the existing NEC basis, without yet changing solved results.

Metrics: exhaustive valid/invalid card tests; no `LD 0`--`LD 5` regression;
coating lookup covers every selected segment exactly once; documented self,
mutual, endpoint and junction terms; `--validate-only` emits stable metadata and
diagnostics.

User acceptance: review sample decks and their machine-readable diagnostics.

### M3 — lossless isolated-dipole operator

Goal: solve a single, uniformly coated, lossless center-fed dipole with the
existing matrix size.

Metrics: both null limits reproduce the bare-wire result to 1e-12 relative or
the existing floating-point floor; matrix symmetry/reciprocity within 1e-12;
converged current, impedance, resonance and pattern comparisons against M1 and
published thin-coating cases.  Initial promotion target is within the combined
OpenEMS numerical uncertainty plus 1% in resonance shift and 2 ohms per
impedance component.

User acceptance: inspect overlay plots and run the provided bare/coated pair.

### M4 — loss and frequency sweeps

Goal: add `TAN_DELTA` and `SIGMA_DC` with per-frequency evaluation and power
accounting.

Metrics: zero-loss identity; equivalence of tangent and conductivity inputs at
one chosen frequency to 1e-12; passive nonnegative absorbed power; power-balance
defect below 0.3%; sweep results agree with OpenEMS within the M3 envelope.

User acceptance: run a VSWR sweep for one low-loss and one lossy insulation.

### M5 — feed-model decision

Goal: compare continuous center-segment excitation against the finite-gap
differential port and decide the supported production forms.

Metrics: gap/excitation convergence tables, de-embedded impedance, current away
from the feed, resonance and pattern.  Differences must be bounded over the
validated envelope; otherwise only the demonstrated form is certified.

User acceptance: choose the final warning/constraint based on the plots and
construct one dipole model through Antenna Turtle or a NEC deck.

### M6 — junction policy experiments

Goal: test equal coated junctions, coated/bare junctions, unequal coatings and
asymmetric hot/cold arms.

Metrics: charge/current continuity, segmentation convergence, reciprocity,
power balance and OpenEMS comparison.  Equal coatings are the first promotion
candidate.  Other cases remain rejected unless they meet the M3/M4 accuracy
budgets.

User acceptance: approve the supported topology matrix and warning/error UX.

### M7 — thickness/permittivity validity map

Goal: quantify how far the reduced operator can extend into 6–12 mm solid
housings while all coatings remain electrically thin.

Research note for the point where the reduced radial-field operator ceases to
meet its error budget: do not assume that `t/a = 1` is a sharp physical
boundary.  Test whether increasing thickness makes axial dielectric field and
radial field variation material, then verify two possible extensions against
the later coated-wire literature:

1. A coupled thin-wire/VIE formulation in which the conductor retains its 1-D
   current basis while the dielectric sleeve carries volume-polarization
   unknowns.  Concentric `rho,z` shells may remain economical for the straight,
   axisymmetric validation article; arbitrary NEC wire orientation, bends,
   junctions and mutual coupling do not generally preserve a global `m=0`
   reduction.
2. A multilayer cylindrical Green-function or modified-kernel formulation that
   analytically satisfies the transverse layered boundary conditions and keeps
   unknowns on the wire.  Evaluate numerical stability and cost of its
   Bessel/Hankel or spectral integrals, and treat finite endpoints, feed gaps
   and coating transitions as separate validation problems rather than
   consequences of the infinite-cylinder kernel.

   As a specific M7 research candidate, formulate each concentric layer as a
   two-port radial field-transfer operator and cascade those operators from the
   PEC core to the outer boundary.  Track the dimensional convention carefully:
   conventional surface impedance is `E_z/H_phi` (ohms), whereas `E_z/I_z`
   carries ohms per metre when `I_z = 2*pi*rho*H_phi`.  The layer transformation
   will generally depend on axial spectral wavenumber as well as frequency, so
   it is a nonlocal longitudinal operator after inverse transformation.  Verify
   whether it can be incorporated efficiently into the full interaction kernel;
   do not presume that replacing only segment self-impedances captures mutual
   coupling, endpoints, bends, junctions, or coating transitions.  Compare the
   resulting one-, two-, and three-layer kernels with OpenEMS before claiming
   that the cascade removes the thick-coating limit.

These are design hypotheses, not yet accepted implementation requirements.
M2's literature review must identify the governing papers and applicability
limits before either path is selected.

Metrics: sweep `b/a`, `k0*b`, `ER`, loss and distance from ends/feed.  Publish
error surfaces for resonance, complex impedance, efficiency and pattern.  The
validated envelope is the largest connected region meeting the M3/M4 budgets;
outside it the requested warning-and-continue policy applies.  Include
flat/flush and representative printed-housing endpoint treatments so the
published envelope distinguishes reduced-operator error from physical end-cap
capacitance.  Archive enough geometric detail to reproduce a future printed
article and VNA/range comparison.

User acceptance: review the envelope and manually exercise a representative
printed-housing model inside and outside it.

### M8 — integration and release

Goal: complete documentation, Antenna Turtle-facing examples, performance
checks and backward-compatibility certification.

Metrics: complete suite and M0 golden pass; uncoated end-to-end output remains
unchanged; coated matrix order equals bare order; coated runtime/memory overhead
is characterized; CLI/result schema consumers accept both coated and uncoated
jobs; all promoted cases have archived OpenEMS evidence.

User acceptance: final example/VSWR/pattern review, then explicit merge verdict.

## Dependencies and administration

M0 needs no new package.  OpenEMS and AppCSXCAD are present under
`/home/josh/opt/openEMS`, while their Python modules are not currently on the
default Python module path.  M1 begins by testing the installed launchers and
locating the packaged Python/Octave interface.  If a system package, path change
or library installation is actually required, work stops and the exact request
is handed to the user; the agent will not perform system administration.

## M0 manual acceptance

From the repository root:

```sh
git branch --show-current
python3 tests/coated-wire/m0_baseline.py
make check
```

Expected branch: `codex/coated-wire`.  Expected baseline message:

```text
coated-wire M0 bare-wire baseline passed
```

M1 begins only after the user reviews this document and returns a PASS.
