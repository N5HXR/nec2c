# Dielectric volume support

**Revision 1 COMPLETE (2026-10-21):** M1–M11 and M13 are accepted;
M12 is intentionally deferred. Start with the
[current user guide](DIELECTRIC-REVISION-1.md),
[canonical examples](../examples/dielectric-r1/README.md) and
[deferred-work register](DIELECTRIC-DEFERRED-WORK.md).
The initial proof-of-concept requirements below are historical design stages;
the later milestone records and revision-1 guide state the implemented scope.

## Status

This document specifies a staged extension to the N5HXR `nec2c` fork for
lossless dielectric volumes.  It is a design and validation contract; the milestone table and validation
notes below state the implemented scope and remaining gates.

The first implementation is a numerically conservative reference solver:

- free-space background;
- real, isotropic relative permittivity;
- nonmagnetic materials;
- dielectric regions represented by unions of non-overlapping, axis-aligned
  cuboids;
- three piecewise-constant polarization-current unknowns per numerical volume
  cell;
- Galerkin volume-integral-equation discretization; and
- the existing dense direct matrix solution.

The design deliberately favors a small, validated problem over an initially
fast one.  FFT acceleration, dielectric/ground contact, conductor/dielectric
contact, loss, dispersion, and very large models are later phases.  A model
without a dielectric extension must retain the existing NEC-2 numerical path
and behavior.

## Motivation and production target

Antenna Turtle already emits NEC2C wire models and invokes `nec2c` as an
external executable.  It also has Manifold-backed CSG and openEMS backends.
The desired extension lets the same workflow model parasitic dielectric
effects on feed impedance, derived S-parameters, and radiation patterns.

The representative production case is an HF antenna near an idealized tree
or building.  A tree may initially use `epsilon_r = 3.0`.  This is a scale and
workflow target, not a claim that a real tree is accurately represented by a
lossless homogeneous material.  Real vegetation will eventually require loss
and heterogeneous-material validation.

## Goals

1. Couple dielectric polarization to the existing wire EFIE and conducting
   surface MFIE equations in one self-consistent MoM system.
2. Preserve existing deck compatibility and results when no dielectric file
   is selected.
3. Use a formulation with published provenance and independent analytical and
   numerical validation.
4. Let Antenna Turtle describe physical dielectric regions without embedding
   another CSG engine in NEC2C.
5. Express mesh quality and complexity limits in electrical as well as raw
   cell-count terms.
6. Establish a dense reference implementation against which future FFT or
   hierarchical solvers can be checked.
7. Leave a documented path to finite ground, ground contact, radomes,
   insulation, and embedded conductors.

## Non-goals for the proof of concept

- Loss tangent, conductivity, complex permittivity, or dispersive materials.
- Magnetic media or anisotropic permittivity.
- Dielectric cells touching or containing wires or conducting patches.
- Dielectric cells touching the ground plane.
- Sommerfeld, reflection-coefficient, cliff, radial-screen, or perfect-ground
  calculations when a dielectric model is active.
- NEC rotational or reflection symmetry reduction when a dielectric model is
  active.
- Production-scale voxel counts.
- Direct calculation of broadband S-parameter matrices in one excitation.
  Antenna Turtle may continue deriving single-port quantities from repeated
  standard NEC excitations.
- PMCHWT or another dielectric surface-integral-equation backend.

## Selected method and provenance

The selected method is an electric-current volume integral equation (JVIE)
discretized with Galerkin testing and piecewise-constant Cartesian basis
functions on cuboidal cells.

This choice follows the stable current-based VIE family described by
Polimeridis, Villena, Daniel, and White.  Their formulation uses equivalent
volumetric currents, Galerkin testing on uniform voxels, accurate singular
integration, and an FFT-accelerated matrix-vector product.  VoxScatter exposes
both DDA and this Galerkin family and reports better conditioning for the
Galerkin formulation at larger refractive indices.

DDA remains a useful independent comparison, especially through DDSCAT, but
is not the production formulation.  Point-dipole polarizability prescriptions
and surface error can otherwise become part of the solver definition.  PMCHWT
remains a possible future backend for clean homogeneous closed surfaces, but
it requires both electric and magnetic equivalent surface currents, conforming
interfaces, singular surface integration, and separate conditioning work.

The implementation must be independently derived from published equations
and openly available reference implementations.  It must not copy restricted
NEC-3 or NEC-4 material.

## Electromagnetic convention

The normative convention for the new derivation is time dependence

```text
exp(+j omega t)
```

with the outgoing free-space scalar Green function

```text
g(r,r') = exp(-j k0 |r-r'|) / (4 pi |r-r'|).
```

This matches the `exp(-j k r)` kernels visible in the existing near-field
routines, but it must be confirmed with executable source/observation tests
before any coupled block is accepted.  A sign copied from a source using the
opposite time convention is not acceptable.

For a nonmagnetic dielectric in free space, define

```text
chi(r) = epsilon_r(r) - 1
Jv(r)  = j omega epsilon0 chi(r) E(r).
```

Let `L_E[J]` denote the electric field radiated in the background medium by an
electric volume current under the convention above.  The total field inside
the dielectric is

```text
E = E_background + L_E[Jv],
```

and the current equation is therefore

```text
[1 / (j omega epsilon0 chi)] Jv - L_E[Jv] = E_background.
```

Cells with `chi = 0` are omitted rather than introducing a singular material
coefficient.

The implementation shall express `L_E` through vector- and scalar-potential
terms, or an algebraically equivalent weak form, rather than directly sampling
the strongly singular dyadic kernel.  For a piecewise-constant cell current,
the charge contribution can be represented by its normal current on cell
faces.  The final derivation must show:

- the time convention and every sign;
- SI dimensions of every unknown, matrix block, and right-hand side;
- treatment of the principal value and self term;
- reduction of singular volume-volume terms to regular or weakly singular
  integrals suitable for cubature; and
- the scaling used to join SI volume equations to NEC2C's wavelength-scaled
  conductor equations.

That derivation is a required implementation artifact.  The equations above
define the continuum problem but do not substitute for the discrete self-term
derivation.

## Coupled block system

Let `Ic` contain the existing wire and conducting-patch unknowns and let `Jv`
contain three Cartesian volume-current coefficients per dielectric cell.  The
augmented system has the form

```text
[ Zcc  Zcv ] [ Ic ] = [ bc ]
[ Zvc  Zvv ] [ Jv ]   [ bv ]
```

where:

- `Zcc` is the existing NEC2C matrix and must not be reformulated for the
  proof of concept;
- `Zcv` tests fields from dielectric volume currents on wire EFIE and patch
  MFIE observations;
- `Zvc` tests fields from wire and patch currents in dielectric cells; and
- `Zvv` is the Galerkin JVIE material-plus-interaction block.

For a plane wave, both `bc` and `bv` contain incident-field terms with signs
consistent with their tested equations.  For a voltage-driven antenna, `bv`
has no independent impressed field; the antenna excites the dielectric through
`Zvc`.  The dielectric reacts on the antenna through `Zcv` in the same solve.

The precise stored signs depend on the existing NEC convention that places
the negative impressed conductor field on its right-hand side.  Unit tests
must establish the block signs from physical limiting cases rather than from
matrix symmetry alone.

Reciprocity provides a strong diagnostic, but the assembled blocks need not be
literal unscaled transposes because the wire, patch, and volume bases and tests
use different normalizations.  The derivation shall state the diagonal scaling
under which the reciprocal relation is expected.

## Unknowns, testing, and units

Each cuboidal numerical cell has three constant Cartesian current-density
coefficients:

```text
Jv(r) = Jx x_hat + Jy y_hat + Jz z_hat, r inside the cell.
```

The coefficients have SI units A/m^2.  The initial volume equations use
volume-averaged Galerkin testing so their right-hand sides remain in V/m:

```text
<f,g>_cell = (1 / Vcell) integral_cell f(r) dot g(r) dV.
```

Using averages, rather than raw volume integrals, also avoids making row
magnitudes depend trivially on cell volume.  Any alternative normalization is
permitted only if it is documented and produces the same physical currents.

NEC2C stores wire and patch geometry in wavelengths during the solve, while
the dielectric input and material laws are SI.  The volume module should use
SI internally.  Narrow adapter functions shall perform all conversions at the
conductor/volume boundary; wavelength factors must not be scattered through
the new code.

The following dimensions remain distinct:

```text
legacy_neq = number_of_wires + 2 * number_of_pec_patches
volume_neq = 3 * number_of_volume_cells
total_neq  = legacy_neq + volume_neq
```

`data.np2m` must retain its existing meaning.  Dielectric solutions must be
stored separately so `cabc()` and patch-current conversion never interpret
volume coefficients as conductor currents.

## Geometry responsibility

Antenna Turtle owns:

- CSG operations and material-region assignment;
- conversion of the Manifold result to non-overlapping, axis-aligned physical
  boxes;
- deterministic region IDs; and
- writing the dielectric resource retained with the NEC deck.

NEC2C owns:

- parsing and validating the resource;
- checking box overlap, bounds, and material references;
- checking separation from conductors, patches, and ground;
- subdividing physical boxes into numerical cells; and
- reporting the resulting electrical resolution and cost.

NEC2C will not perform Boolean union.  Boxes may share faces, edges, or corners
within the same material region, but their interiors must not overlap.  The
proof of concept rejects touching regions with different material IDs.

The distinction between physical boxes and numerical cells is intentional.
The retained Antenna Turtle geometry remains stable while NEC2C selects a mesh
for the highest requested frequency.

## `DF` card

`DF` is an N5HXR NEC2C extension card in the control-card section after `GE`.
It refers to one dielectric resource:

```text
GE 0
DF "dielectric/tree.diel"
FR ...
EX ...
RP ...
EN
```

Requirements:

- zero or one `DF` card per NEC job in the first implementation;
- paths may be quoted and may contain spaces;
- an overlong or unterminated card is an error and must never be silently
  truncated; the first implementation may retain the existing NEC input-line
  limit if the diagnostic states that limit explicitly;
- a relative path is resolved against the directory containing the NEC input
  deck, not the process working directory;
- an absolute path is accepted but reported as nonportable;
- nested includes are not supported;
- the resource is parsed during card processing, while numerical subdivision
  is deferred until the complete frequency request is known;
- a changed `DF` invalidates the electromagnetic matrix;
- `NX` releases all dielectric state before the next job; and
- a deck without `DF` never enters dielectric parsing, allocation, or solve
  code.

The current numeric `readmn()` path cannot parse a filename.  `DF` therefore
requires a narrowly scoped free-form card parser before numeric control-card
decoding.  Other card parsing must remain unchanged.

## Dielectric resource grammar

The initial format is line-oriented ASCII and dependency-free.  Keywords are
case-insensitive; identifiers are positive decimal integers; numbers use C
locale decimal syntax.  `#` begins a comment outside a quoted string.

```text
NEC2C-DIELECTRIC 1

MATERIAL 1 ER 3.0

REGION 101 MATERIAL 1
BOX 101 -0.50 -0.50 0.10  0.50 0.50 10.00
BOX 101 -2.00 -2.00 8.00  2.00 2.00 15.00

MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2.0
END
```

Normative records:

```text
NEC2C-DIELECTRIC <format-version>
MATERIAL <material-id> ER <real-relative-permittivity>
REGION <region-id> MATERIAL <material-id>
BOX <region-id> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax>
MESH MAX_KH <q> MAX_CELLS <count> MAX_ASPECT <ratio>
END
```

Rules:

- coordinates are meters in the NEC global Cartesian frame;
- every minimum bound must be strictly less than its maximum bound;
- `epsilon_r` must initially be finite and greater than or equal to 1;
- `epsilon_r = 1` is permitted for validation but its cells are omitted from
  the physical solve;
- every box references a declared region and every region references a
  declared material;
- IDs need not be contiguous and file order does not define physical meaning;
- duplicate IDs and unknown keywords are errors;
- zero or one `MESH` record is allowed; omitted values use the documented
  defaults, and duplicate options in one record are errors;
- all non-comment content after `END` is an error;
- `GS` in the NEC deck does not scale the external resource; both files must
  already describe the same final geometry in meters; and
- multiple declared materials are parsed from the start, although an early
  solver milestone may reject a job that uses more than one non-vacuum
  material.

The parser reports resource filename, source line, record type, and offending
value on every error.  Antenna Turtle may use any filename suffix; `.diel` is
the documented convention.

## Meshing policy

Numerical subdivision is based on the highest frequency that will use the
matrix.  Define for each material and cell

```text
q = k0 sqrt(epsilon_r) hmax
  = 2 pi hmax / lambda_material,
```

where `hmax` is the cell's longest edge.

Initial policy:

- `q <= 0.5`: recommended default for impedance and pattern work;
- `0.5 < q <= 1.0`: permitted only with an accuracy warning;
- `q > 1.0`: error unless an explicit diagnostic override is added in a later
  format version;
- `MAX_ASPECT` defaults to `2.0`; and
- each physical box is divided into equal intervals per axis so every child
  meets both `MAX_KH` and `MAX_ASPECT`.

`q = 0.5` corresponds to about 12.6 cells per internal wavelength and `q = 1`
to about 6.3.  These are conservative starting guards motivated by DDSCAT's
published `|m| k d` criteria.  They are not assumed to be the final JVIE error
model.  Layer-1 convergence tests will determine whether the defaults should
change.

Wavelength resolution is necessary but not sufficient.  Reports also include:

- cells across each physical box dimension;
- smallest physical feature relative to local cell size;
- region electrical diameter `k0 sqrt(epsilon_r) D`;
- actual occupied volume versus requested physical volume; and
- warnings for features represented by only one cell across.

The proof of concept defaults to `MAX_CELLS 1000` and recommends no more than
500.  The effective limit also includes conductor unknowns and estimated
factorization memory.  The program calculates cost before matrix allocation
and fails early with an actionable report.

For complex-double dense storage, raw matrix memory is approximately

```text
16 * total_neq^2 bytes,
```

before pivots, right-hand sides, scratch arrays, and factorization overhead.
Overflow-safe `size_t` arithmetic is mandatory for every dimension and byte
calculation.

## Cell interaction evaluation

The reference solver must prioritize reproducibility over speed.

Interactions are classified as:

1. identical cell;
2. face-, edge-, or corner-adjacent cells;
3. geometrically near cells; and
4. well-separated cells.

Identical and adjacent cases may not use center-point sampling.  Their
weakly/strongly singular terms require an analytical reduction or a published
singularity-removal transformation.  Near cases use adaptive cubature with a
reported tolerance.  Well-separated cases may use fixed-order tensor-product
quadrature only after comparison with a higher-accuracy oracle.

The preferred implementation follows the current-based Galerkin reduction of
volume-volume terms to smoother surface-surface integrals.  A slow,
high-precision evaluator shall be retained in tests as an oracle for randomly
selected cell pairs.

Required numerical tests include:

- invariance under translation and Cartesian-axis permutation;
- reciprocal cell-pair interactions after documented normalization;
- continuity as separated cells approach face contact;
- agreement of fast and oracle quadrature over a sweep of size, separation,
  aspect ratio, and frequency; and
- stable static-limit behavior as `k0` approaches zero without evaluating a
  singular `1/omega` expression naively.

## Existing solver insertion points

The present solver constructs `n + 2m` conductor equations in `cmset()` and
then calls the dense factorization and symmetry-aware solve routines.  The
extension should be isolated in new modules, for example:

```text
dielectric.h
dielectric_io.c
dielectric_mesh.c
vie_operator.c
vie_quadrature.c
dielectric_coupling.c
dielectric_fields.c
```

Suggested internal interfaces:

```c
bool dielectric_load(const char *deck_path,
                     const char *resource_path,
                     dielectric_model *model,
                     dielectric_diagnostic *diagnostic);

bool dielectric_mesh(const dielectric_model *model,
                     const frequency_plan *frequencies,
                     dielectric_cells *cells,
                     dielectric_diagnostic *diagnostic);

void vie_apply(const vie_operator *op,
               const complex double *input,
               complex double *output);
```

The dense reference implementation may materialize `vie_operator` as a matrix,
but conductor code must call through a representation-independent boundary so
an FFT or hierarchical operator can replace it later.

Before extending matrix dimensions, audit and refactor every use of:

- `data.np2m` and `data.np3m`;
- `netcx.neq` and `netcx.npeq`;
- pivot and scratch-array allocation;
- `factrs()`, `solve()`, and `solves()`;
- `netwk()` and coupled-network right-hand sides;
- `cabc()` and current output;
- near- and far-field dispatch; and
- symmetry-mode transforms.

Generic dense linear algebra routines must allocate scratch space from their
explicit matrix order, not from the legacy conductor count.

## Compatibility and feature gates

When `DF` is absent:

- `legacy_neq == total_neq`;
- matrix allocation and fill follow the existing branch;
- no dielectric files are opened;
- no output headings are added; and
- retained regression decks should be byte-for-byte identical where platform
  formatting permits, with parsed numerical equality as the portable gate.

When `DF` is present in the proof of concept, reject:

- any non-free-space `GN` configuration;
- active NEC symmetry reduction (`GX` or `GR`-derived `data.ipsym`);
- dielectric/conductor intersection or insufficient configured clearance;
- unsupported network combinations until `netwk()` has been audited for the
  augmented equation count; and
- model dimensions exceeding safe memory or cell limits.

Rejecting symmetry reduction does not reject symmetric geometry.  Antenna
Turtle may emit every reflected or rotated object explicitly.

## Conductor clearance

The proof of concept requires a positive gap between every dielectric cell and
wire centerline/surface or conducting patch.  The initial conservative rule is

```text
minimum gap >= max(local volume-cell longest edge,
                   neighboring wire segment length),
```

with the exact default subject to conductor/volume quadrature experiments.
The report prints both the measured minimum and required minimum.  A later
milestone replaces this coarse restriction with accurate near-singular
coupling integration.

M3 computes wire clearance from the exact segment-to-axis-aligned-cell-box
distance, less wire radius.  The legacy processed-patch representation retains
patch center and area but not exact polygon extents. M7 deliberately retains
an equivalent-area-radius proxy under its separated-object constraint; exact
polygon-aware clearance is deferred to M12. This proxy is an explicit
preflight limitation, not a certified intersection test for elongated patches.

A radome separated from its antenna is therefore an earlier supported case
than insulation touching a wire.

Free space has no distinguished `z = 0` plane, so the proof of concept does
not reject negative dielectric coordinates.  Once a ground mode is enabled,
the pre-contact phase requires every dielectric cell to remain above the
configured ground interface by its documented positive clearance.

## Excitation, solution, and output

Initial supported excitations:

- ordinary voltage source on a wire segment; and
- incident plane wave.

The dielectric affects the same solved conductor currents used by existing
feed impedance, network, near-field, and pattern calculations.  Volume-current
contributions must be added to near and far fields before radiation-pattern
validation is considered complete.

Internal dielectric fields and polarization currents need not be printed in
normal output.  They must remain programmatically available to tests because
the analytical internal field of a sphere is an important independent check.

Add one opt-in dielectric report section containing:

- resource format version and resolved path;
- material and region tables;
- physical-box and numerical-cell counts;
- frequency used for meshing;
- minimum, maximum, and percentile `q`;
- cell aspect ratios and feature warnings;
- minimum conductor and ground clearances;
- conductor, volume, and total unknown counts;
- matrix and workspace memory estimates;
- fill, factor, and solve timings; and
- normalized residual and validation warnings.

Classic output sections and labels remain unchanged.

## Validation strategy

Validation is layered so a conductor-coupling error is not mistaken for a VIE
operator error.

### Layer 0: legacy regression

Create a machine-readable regression runner over representative existing
decks:

- free-space wire antenna;
- finite/perfect/Sommerfeld ground where currently supported;
- patch-only and mixed wire/patch models;
- `NT` and `TL` networks;
- near electric and magnetic fields; and
- radiation patterns.

Capture selected complex currents, feed impedances, powers, and pattern
samples rather than relying only on whole-text diffs.  The no-`DF` branch is a
release blocker if these change unexpectedly.

### Layer 1: analytical and algebraic dielectric tests

#### V1: vacuum identity

Run the same geometry with no cells and with cells declared at `epsilon_r=1`.
After vacuum cells are omitted, matrix order and results must match the legacy
case.  Target: roundoff-level numerical equality.

#### V2: Rayleigh sphere

For `k0 a <= 0.1`, compare induced dipole moment against

```text
alpha = 4 pi epsilon0 a^3 (epsilon_r - 1) / (epsilon_r + 2).
```

Use `epsilon_r = 2, 3, 4, 9` and at least four boundary refinements.  Record
requested radius, occupied cuboid volume, and equivalent-volume radius so
voxelized-boundary error is distinguishable from operator error.

#### V3: sphere internal field

Compare volume-averaged internal field against

```text
Einside = 3 E0 / (epsilon_r + 2)
```

in the Rayleigh/static regime.

#### V4: Lorenz-Mie sphere

At several size parameters through approximately `k0 a = 1`, compare:

- total scattering cross section;
- complex far field at selected angles; and
- normalized angular-pattern RMS error.

Refine the cuboid representation until both geometry and field discretization
show a clear convergence trend.

#### V5: lossless power balance

Absorbed dielectric power is zero.  Extinction/scattering for plane-wave tests
and accepted/radiated power for driven tests must close within the combined
linear-solve, quadrature, and mesh error budget.

#### V6: reciprocity and invariance

Test reciprocal source/observer pairs, translations, Cartesian-axis
permutations, and equivalent physical decompositions of the same rectangular
region.

Preliminary Layer-1 gates:

- normalized dense-system residual no worse than `1e-10` for well-conditioned
  reference problems;
- reciprocal interaction discrepancy no worse than `1e-9` after documented
  basis/test scaling;
- roundoff-level vacuum identity;
- monotone or explainably asymptotic refinement toward analytical results;
  and
- no more than 2% error in agreed Mie integral quantities at a documented
  feasible reference resolution.

The thresholds may be tightened after the first convergence study.  They may
not be loosened merely to accept a single unresolved model.

### Layer 2: independent numerical comparison

Use openEMS because its EC-FDTD discretization is independent of the integral
equation implementation.  Retain generator inputs, openEMS version, mesh,
boundary conditions, convergence settings, and raw results with each fixture.

#### N1: dielectric cube scattering

- dimensions `0.2 lambda0 x 0.1 lambda0 x 0.1 lambda0`;
- `epsilon_r = 3` and `4`;
- multiple incidence directions and polarizations; and
- compare selected complex fields and angular patterns.

#### N2: separated dielectric-loaded dipole

- thin center-fed dipole approximately `0.475 lambda0` long;
- dielectric cube side `0.1 lambda0`;
- lateral center offsets `0.1`, `0.2`, and `0.4 lambda0` where the PoC
  clearance rule permits;
- `epsilon_r = 2, 3, 4, 9`; and
- sweep `0.8 f0` through `1.2 f0`.

Compare complex feed impedance, resonance shift, accepted power, and principal
pattern cuts after independent refinement in both solvers.

#### N3: idealized free-space tree

- `epsilon_r = 3`;
- axis-aligned box decomposition of trunk, branches, and canopy;
- several antenna separations and cell limits; and
- report convergence, memory, fill/factor time, impedance shift, and pattern
  distortion.

This is initially a production-scale and stability benchmark rather than an
analytical truth model.

Preliminary Layer-2 targets after both solvers have converged:

- resistance and reactance within the greater of 2% or 2 ohms;
- principal resonance frequency within 1%; and
- principal normalized pattern cuts within an agreed RMS tolerance, initially
  0.25 dB away from deep nulls.

Disagreement must be investigated through independent mesh refinement and
power balance; openEMS is a comparator, not an oracle.

## Tracked milestones

There are fourteen tracked milestones, `M0` through `M13`.  Each milestone
must produce a usable or independently testable result.  A milestone is not
complete merely because its code path exists.

Use these status values in the table below:

- `NOT STARTED`: no implementation work has been accepted;
- `IN PROGRESS`: work exists but at least one exit criterion is unmet;
- `DEFERRED`: intentionally outside the current revision; retained for future work;
- `BLOCKED`: progress requires a recorded external decision or prerequisite;
- `IN REVIEW`: all exit criteria appear met and evidence is under review; and
- `COMPLETE`: every exit criterion and artifact has been accepted.

Design changes discovered during a milestone are expected.  Record them in
this document or a linked decision note.  A change may revise future exit
criteria.  It may revise a completed milestone only when new evidence shows
that the earlier result is invalid; in that case, explicitly reopen it rather
than silently editing its meaning.

| ID | Concrete result | Status |
| --- | --- | --- |
| M0 | Reviewed design and validation contract | COMPLETE |
| M1 | Reproducible legacy regression baseline | COMPLETE |
| M2 | `DF` and dielectric-resource parser | COMPLETE |
| M3 | Deterministic mesher and cost report | COMPLETE |
| M4 | Certified cuboid interaction kernel | COMPLETE |
| M5 | Validated dielectric-only dense JVIE solver | COMPLETE (2026-10-18) |
| M6 | Voltage-driven wire/dielectric coupling | COMPLETE (2026-10-19) |
| M7 | Conducting-patch/dielectric coupling | COMPLETE (2026-10-19) |
| M8 | Validated free-space dielectric release | IN PROGRESS |
| M9 | Production-scale free-space tree solver | COMPLETE |
| M10 | Dielectric above perfect ground with a positive gap | COMPLETE (2026-10-20) |
| M11 | Dielectric/perfect-ground contact | COMPLETE (2026-10-21) |
| M12 | Touching and embedded conductors | DEFERRED (later revision) |
| M13 | Loss, real ground and advanced material behavior | IN PROGRESS (scalar loss; preliminary real-ground work) |

The status table is the project-level progress tracker.  Test reports and
benchmark artifacts provide the evidence behind it.

Milestone scope decisions (2026-10-20): the current proof of concept remains
lossless through M11, using real, lossless dielectric permittivity and ideal
PEC conductors/ground. M10 uses `GN 1`; M11 extends dielectric contact with
that perfect ground. M12's wire/patch contact with, penetration into and
embedding in dielectric are deferred to a later revision. The current revision
retains documented positive wire/patch-to-dielectric clearance; M11's PEC-ground
interface is a separate contact boundary. M12's identifier and proposed criteria
are retained for tracking, not as a prerequisite for the next stage.

Real-ground approximations, conductivity and substrate-loss validation belong
to M13, alongside dielectric loss and advanced materials. M13 can proceed on
separated conductor/dielectric geometries without completing M12; any additional
contact scope must be selected explicitly. Existing GN 2 code and evidence are
retained as preliminary M13 work. Neither finishing real-ground work nor
obtaining further real-ground OpenEMS results blocks M10 or M11. These decisions
supersede the broader real-ground scope in historical M10 progress notes below.
M10's PEC validation gates remain in force. The planned current-revision path
is M10, M11, then M13, with M12 deferred.

### M0: reviewed design and validation contract

**Result:** A design that fixes the initial formulation, interfaces, scope,
validation layers, and milestone gates while clearly identifying unresolved
numerical decisions.

**Included:**

- this specification;
- references and provenance;
- selected JVIE family and initial constraints;
- `DF` and dielectric-resource design;
- Layer-0, Layer-1, and Layer-2 validation plans; and
- explicit deferred decisions.

**Exit criteria:**

- NEC2C and Antenna Turtle maintainers agree that the initial UI and division
  of geometry responsibilities are implementable;
- no known contradiction remains between PoC scope and a later production
  tree model; and
- unresolved mathematical choices are assigned to a later milestone rather
  than hidden as implementation details.

**Artifacts:** This document and any accepted review amendments.

### M1: reproducible legacy regression baseline

**Result:** One command builds NEC2C, runs representative unmodified decks,
and compares stable numerical outputs against reviewed fixtures.

**Included:** wire-only, patch-only, mixed wire/patch, voltage-source,
plane-wave, near-field, radiation-pattern, ground, `NT`, and `TL` coverage.

**Excluded:** All `DF` parsing and dielectric calculations.

**Exit criteria:**

- the suite passes on the primary development platform from a clean build;
- selected complex currents, impedances, powers, and pattern samples are
  compared numerically with documented tolerances;
- at least one expected-failure fixture proves the runner detects a numerical
  regression; and
- fixture provenance and regeneration instructions are committed.

**Artifacts:** Test runner, legacy fixtures, tolerance policy, and baseline
report.

### M2: `DF` and dielectric-resource parser

**Result:** A normal NEC deck can select a deck-relative dielectric resource,
and NEC2C can parse and report its materials, regions, and physical boxes
without attempting an electromagnetic dielectric solve.

**Included:** quoted paths, path resolution, format version, all initial
records, detailed source-line diagnostics, cleanup at `NX`, and no-`DF`
compatibility.

**Excluded:** numerical subdivision, matrix-size estimation, and dielectric
fields.

**Exit criteria:**

- valid parser fixtures reproduce a canonical model summary;
- malformed, duplicate, missing, non-finite, and unsupported records fail with
  deterministic diagnostics;
- a second `DF` and an overlong/truncated card fail explicitly;
- relative paths work independently of process working directory; and
- the complete M1 suite remains unchanged.

**Artifacts:** Parser implementation, format fixtures, parser tests, and
resource-format documentation suitable for Antenna Turtle.

### M3: deterministic mesher and cost report

**Result:** NEC2C converts physical boxes into deterministic numerical cells
using the highest requested frequency and reports whether the requested model
is computationally admissible.

**Included:** `MAX_KH`, `MAX_CELLS`, `MAX_ASPECT`, material wavelength,
overflow-safe unknown/memory estimates, box-overlap validation, and
conductor-clearance measurement.

**Excluded:** Green functions, matrix fill, and solution.

**Exit criteria:**

- repeated runs and reordered equivalent input records produce the same
  canonical cell ordering;
- overlapping box interiors and touching different-material regions fail with
  deterministic geometry diagnostics;
- every emitted cell satisfies the requested electrical size and aspect-ratio
  rules within documented floating-point tolerance;
- occupied numerical volume equals physical box-union volume within tolerance;
- excessive models fail before matrix allocation with an actionable report;
- boundary cases at `q = 0.5` and `q = 1.0` have tested warning/error behavior;
  and
- the M1 and M2 suites pass.

**Artifacts:** Mesher, cost estimator, geometry fixtures, and report examples.

### M4: certified cuboid interaction kernel

**Result:** A standalone numerical library evaluates every Cartesian component
of the Galerkin interaction between two cuboidal volume-current basis/test
functions, including self and adjacent cases.

**Included:** documented discrete equations, time convention, SI units,
singularity treatment, high-accuracy oracle, and far/near evaluation policy.

**Excluded:** global matrix assembly, plane-wave solving, and conductor
coupling.

**Exit criteria:**

- the derivation is independently reviewed against the cited JVIE literature;
- self, face-, edge-, corner-adjacent, near, and far cases pass oracle
  comparisons over documented size/frequency/aspect-ratio grids;
- translation, Cartesian permutation, reciprocity, and static-limit tests
  pass their specified tolerances;
- no center-point approximation is used for self or adjacent interactions;
  and
- numerical tolerances and failure behavior are reproducible.

**Artifacts:** Derivation note, interaction implementation, oracle, generated
reference data, and test report.

### M5: validated dielectric-only dense JVIE solver

**Result:** A `DF`-only free-space model under plane-wave excitation produces
internal polarization currents, total dipole moment, and scattered far fields
using the dense reference solver.

**Included:** volume matrix assembly, material term, plane-wave RHS, direct
factorization, residual reporting, and test-only internal-field access.

**Excluded:** wires, conducting patches, ground, and accelerated solution.

**Exit criteria:**

- vacuum identity V1 passes;
- Rayleigh polarizability V2 and sphere internal-field V3 show convergence;
- Lorenz-Mie sphere V4 meets the preliminary integral and pattern tolerances;
- lossless power balance V5 and invariance/reciprocity V6 pass;
- refinement data distinguishes shape error from operator error; and
- DDSCAT or VoxScatter comparison fixtures agree within their documented
  discretization limits.

**Artifacts:** Dense JVIE solver, analytical reference calculator or reviewed
fixtures, convergence plots/data, and Layer-1 validation report.

### M6: voltage-driven wire/dielectric coupling

**Result:** A conventional voltage-driven NEC wire antenna and separated
dielectric volume are solved self-consistently, producing dielectric-modified
wire currents and feed impedance.

**Included:** generic total matrix dimensions, wire-to-volume and
volume-to-wire blocks, voltage-source RHS, separate current storage, and
coupled reciprocity tests.

**Excluded:** conducting patches, dielectric radiation contribution to final
patterns, NEC ground, and conductor contact.

**Exit criteria:**

- the M1 suite passes after the linear-algebra dimension refactor;
- removing the dielectric or setting `epsilon_r = 1` restores the legacy
  wire result;
- the dielectric effect approaches zero with increasing antenna separation;
- block reciprocity and full-system residual meet their specified tolerances;
- feed impedance converges under both wire and dielectric refinement; and
- the N2 dipole/cube impedance and resonance shift agree with independently
  refined openEMS results within preliminary Layer-2 limits.

**Artifacts:** Coupled matrix implementation, wire/cube fixtures, refinement
data, and openEMS comparison report.

### M7: conducting-patch/dielectric coupling

**Result:** Existing conducting surface patches participate in the same
self-consistent dielectric solve as wires.

**Included:** patch-to-volume electric fields, volume-to-patch magnetic-field
testing consistent with the existing MFIE, mixed wire/patch models, and
patch-block reciprocity/scaling tests.

**Excluded:** dielectric touching a patch and ground.

**Exit criteria:**

- patch-only and mixed wire/patch vacuum identities pass;
- all four conductor/volume block directions have isolated tests;
- patch results converge under surface and volume refinement;
- at least one separated PEC/dielectric model agrees with openEMS; and
- M1 through M6 remain passing.

**Artifacts:** Patch coupling implementation, isolated block fixtures,
convergence data, and comparator model.

### M8: validated free-space dielectric release

**Result:** The dependency-light dense implementation is a complete supported
free-space feature for models within its published size and clearance limits.

**Included:** volume contributions to near and far fields, lossless power
accounting, supported `NT`/`TL` audit, stable dielectric reporting, Antenna
Turtle export, retained artifacts, N1 and N2 comparisons, and user guidance.

**Excluded:** production-scale trees, ground, and touching conductors.

**Exit criteria:**

- Layer-0 and Layer-1 suites pass;
- N1 cube scattering and N2 loaded dipole pass the Layer-2 gates;
- every accepted excitation/output combination is tested and every rejected
  combination has a clear diagnostic;
- Antenna Turtle can generate, execute, retain, and parse a complete model;
- no-`DF` decks remain compatible; and
- limitations and convergence procedure are documented for users.

**Artifacts:** Versioned free-space feature, Antenna Turtle integration,
validation bundle, examples, and release notes.

#### M8 implementation progress (2026-10-19)

The first M8 increment provides combined conductor/volume complex far fields,
exterior scattered electric and magnetic near fields, and driven lossless
power diagnostics. `tests/dielectric_fields.py` runs without openEMS, NumPy,
or another electromagnetic solver and is part of `make check`.

Field conventions and current limits:

- `RP` mode 0 reports complex SI `r E`, with the outgoing
  `exp(-j k r) / r` factor omitted. The conductor contribution is NEC's
  `ffld()` result multiplied by `data.wlam`. The volume contribution uses
  the exact constant-cell Fourier integral, including its sinc factors.
  Complex amplitudes are summed before intensity is calculated.
- `NE` and `NH` report the scattered field from all solved currents, without
  adding the incident plane wave. For a driven source these are its radiated
  electric or magnetic fields. Cartesian and legacy spherical sampling use
  the existing card coordinates. Each observation must be separated from
  every active dielectric cell by at least that cell's longest edge. Internal,
  interface, and closer observations are rejected; the M5 cell-averaged
  internal-field diagnostic remains available separately.
- Near-field cubature reuses the separated tenth-order electric and magnetic
  kernels used by conductor coupling. The electric field is
  `(k0^2 I + grad grad) integral(G J) / (j omega epsilon0)` and the magnetic
  field is `integral(grad(G) cross J)` under `exp(+j omega t)`.
- Ordinary voltage-driven models report accepted power
  `0.5 Re(sum(V conj(I)))`, radiated power integrated over the complete sphere,
  their relative defect, and an angular-integration refinement diagnostic
  comparing 974 and 3896 directions. No conservation correction is applied.
  The retained M6 vacuum and loaded dipoles both have approximately 0.29%
  power defect, with angular-integration change below `4e-7`.
- The initial analytical gates are `1e-4` relative error against the exact
  constant-cube far-field form factor, `4e-4` against electric-dipole exterior
  E/H fields, `6e-4` for near-to-far consistency, `3e-4` against the rounded
  legacy wire/patch/mixed vacuum field output, and `0.5%` driven power defect
  on the retained M6 examples. These are fixture-specific gates, not accuracy
  claims for all admissible models.

The `NT`/`TL` audit finds that legacy `netwk()` allocates its response vector
from `data.np3m` and invokes the legacy conductor factorization/solve path.
The coupled solver instead owns a separate augmented dense matrix. Enabling
networks requires port-response solves against that augmented matrix and
consistent network power accounting; changing the legacy equation count is
insufficient. `DF` with either network card remains explicitly rejected.
No-`DF` network behavior remains covered by the legacy regression suite.

M8 remains **IN PROGRESS**. Remaining release work includes the complete
excitation/output acceptance matrix, broader coupled/patch and plane-wave
power/refinement validation, network support or an accepted release-scope
decision, N1 and expanded N2 independent comparisons, Antenna Turtle
end-to-end integration, and finalized reporting, examples, and release notes.
M6/M7 acceptance remains unchanged: component-wise dielectric-induced
impedance differences use the greater of 2% or 2 ohms. Analytical validation
is prioritized locally; additional openEMS execution can be performed on a
separate machine when the independent-comparison stage is reached.

### M9: production-scale free-space tree solver

**Result:** The idealized free-space tree benchmark runs at an agreed
production cell count without changing the validated dense-reference result
on overlapping sizes.

**Included:** iterative operator interface, selected acceleration method,
preconditioning, convergence diagnostics, and N3 tree benchmarks.

**Excluded:** ground and realistic tree loss.

**Exit criteria:**

- a recorded decision compares regular-grid FFT, sparse-grid cost, and at
  least one hierarchical alternative against the tree geometry;
- accelerated and dense solutions agree within a tolerance substantially
  below discretization error on every overlapping fixture;
- iterative residual correlates with error in feed impedance and pattern;
- the agreed tree-size model fits its memory and runtime budget; and
- solver nonconvergence fails explicitly rather than returning plausible
  output.

**Artifacts:** Acceleration decision note, implementation, dense comparison
suite, N3 performance report, and production sizing guidance.

#### M9 first increment (2026-10-19)

`vie_iterative.h` defines a representation-independent operator callback and
an optional inverse-preconditioner callback. `vie_iterative.c` implements
right-preconditioned restarted complex GMRES with two modified Gram-Schmidt
passes and complex unitary Givens rotations. The implementation follows the
algorithm described in Netlib's *Templates for the Solution of Linear
Systems*, independently implemented here:
<https://www.netlib.org/linalg/old_html_templates/subsection2.6.3.4.html>.

Acceptance requires a recomputed, unpreconditioned residual:

```text
||b - A x||_2 <= max(absolute_tolerance, relative_tolerance * ||b||_2).
```

An estimated Arnoldi residual alone never certifies convergence. The initial
options are restart 30, at most 300 Arnoldi iterations, relative tolerance
`1e-12`, and zero absolute tolerance. The API permits other options. The
operator returns failures explicitly; invalid options, nonfinite arithmetic,
workspace overflow/allocation failure, Krylov breakdown, and iteration limits
are distinguished. Reports expose status, iteration count, restart count,
operator-application count, and true residual. A failed solve propagates as a
failed NEC2C job before solved-current/field output.

For a dielectric deck, the experimental command is:

```sh
./nec2c -s gmres -i tests/dielectric/m6-cube.nec -o /tmp/m9.out
```

`-s dense` is the default. Without `DF`, both settings retain the legacy solve.
The initial GMRES adapter retains the existing dense matrix and uses its
complex diagonal as a right Jacobi preconditioner. It is explicitly reported
as **DENSE-OPERATOR GMRES (EXPERIMENTAL; NO ACCELERATION)**. It does not raise
cell limits or reduce quadratic matrix storage. This increment tests the
iterative solve independently before changing interaction representation.

`tests/vie_iterative_test.c` covers a matrix-free asymmetric complex stencil,
restarts, known-solution accuracy, complex Jacobi preconditioning, exact
initial guesses, zero RHS, absolute tolerance, genuine/happy breakdown,
iteration limits, operator/preconditioner failure, nonfinite data, invalid
options, RHS aliasing, and allocation-size overflow. The physical comparison
in `tests/dielectric_iterative.py` checks dense versus GMRES currents and
complex principal cuts for volume-only, wire, patch, and mixed fixtures.
Compact patch/mixed comparisons use a single 50 mm dielectric cell to keep
this acceptance test inexpensive. It compares reported observables to
`1e-8` relative / `1e-12` absolute and requires GMRES true residual `<=1e-12`.
Power diagnostics are compared when available. These are solver-equivalence
gates on identical meshes, not new physical-accuracy claims.

#### N3 sizing evidence and acceleration decision status

`tools/m9_tree_sizing.py` inventories disjoint, lattice-aligned boxes without
allocating a dense matrix or enumerating all occupied cells. Its provisional
free-space tree has a 0.5 m square, 8 m trunk, two horizontal branches, and a
6 m by 6 m by 4 m box canopy. A separate skeleton omits the canopy. The
reference inventory uses `epsilon_r=3` at 28 MHz, with cell pitches 0.25,
0.125, and 0.0625 m. These are sizing scenarios, not an agreed production
accuracy, runtime, or memory target. They are regenerated with:

```sh
python3 tools/m9_tree_sizing.py --output tests/m9-tree-sizing.csv
```

| Full-tree pitch | Active cells | Two dense matrices | FFT arrays | GMRES(30) workspace |
| --- | ---: | ---: | ---: | ---: |
| 0.25 m | 9,432 | 25.6 GB | 100.7 MB | 28.5 MB |
| 0.125 m | 75,456 | 1.64 TB | 805.3 MB | 228.2 MB |
| 0.0625 m | 603,648 | 104.9 TB | 6.44 GB | 1.83 GB |

Units in this table are decimal. Counts use volume unknowns only; conductor
unknowns and their coupling blocks are additional. The candidate FFT inventory assumes
power-of-two padding of at least `2*naxis-1`, six symmetric interaction tensor
spectra and six vector grids (twelve complex arrays). It excludes FFT-library
workspace, maps, geometry, conductor coupling, and allocator overhead. GMRES
workspace is counted separately. It is not a peak-memory or runtime claim.
The two-matrix column counts only the existing matrix and residual/reference
copy; other dense allocations are additional.

The dense/FFT/sparse/hierarchical comparison is currently:

- Regular-grid FFT: a candidate for the filled-canopy model. It preserves the
  uniform-cell translation structure and would retain certified singular/near
  interactions in the convolution kernel. The existing per-box mesher does
  not generally guarantee a common lattice, so lattice validation or an
  explicit uniform-mesh policy is required before enabling this backend.
- Sparse geometry in a regular FFT embedding: the full-tree grid is 34.1%
  occupied; the skeleton is only 3.99% occupied. FFT arrays still cover the
  padded bounding box, rather than just active cells. At 0.0625 m the skeleton
  has 13,824 active cells but needs the candidate 805.3 MB FFT arrays.
  Restricting storage to occupied voxels does not by itself implement a sparse
  FFT or remove convolution cost. Directed pair counts and translation-table
  bounds are retained in the CSV for comparison with direct/cache approaches.
- Hierarchical alternative: the inventory counts occupied octree nodes and
  leaves with a 32-cell leaf capacity, without replacing empty space by
  dielectric. The finest skeleton has 1,995 nodes / 1,728 leaves; the finest
  full tree has 86,343 nodes / 75,456 leaves. These geometric counts motivate
  investigation of an H-matrix or multipole operator, but do not predict
  numerical ranks, accuracy, memory, or speed. Those require measured
  compression/error and kernel-generation benchmarks.

On 2026-10-19 the user selected the current 2 GB RAM machine and the small
9,432-cell tree as the first experimental target. Use the existing 0.25 m,
28 MHz, epsilon_r=3 scenario initially, then measure what is practical before
choosing a runtime budget or attempting larger cases. Total machine RAM is
the ceiling, not an allocation allowance; initial experiments should cap
solver-owned allocations at 512 MiB to leave room for the OS and other work.

The current implementation experiment is an isolated uniform-cube FFT volume
operator behind the existing GMRES callback. Validate convolution indexing,
self/adjacent interactions, far-interaction quadrature, and true residuals
against the certified dense operator on small meshes before running the full
tree. Retain measured peak RSS, kernel-build time, matvec/solve time, iteration
counts, and convergence status. The first tree experiment may use plane-wave
excitation to isolate volume acceleration; antenna coupling and observable
error checks remain M9 release gates. The isolated FFT operator and benchmark are now implemented in
`vie_fft.c` and `tools/m9_fft_benchmark.c`. The previous dense-GMRES prototype
and all ten regression programs passed before this increment. The FFT
experiment and its independent validation are recorded in
[M9-FFT-EXPERIMENT.md](M9-FFT-EXPERIMENT.md). The later production integration described below now connects this backend
to the NEC2C deck solver; validation and milestone closure are recorded below.

The first 9,432-cell FFT tree solve converged in 26 iterations to true
relative residual `9.83e-9`, with peak RSS 126,384 KiB (about 123 MiB).
Kernel build and solve took 40.5 s and 238.2 s on this machine. Solver-owned
accounted allocations, including GMRES and caller arrays, totaled
130,866,512 bytes (124.8 MiB), below the initial 512 MiB cap. These are
exploratory measurements with concurrent development work, not an agreed
production runtime guarantee. The 216-cell skeleton also converged, and its
`1e-8`/`1e-10` tolerance comparison changed sampled complex far fields by
`8.30e-10` relative to the largest sample. See the experiment report for raw
results, independent dense checks, power checks and limits of these claims.
A tighter full-tree run converged in 33 iterations to `9.80e-11` true residual,
with 129,840 KiB peak RSS; its lossless power defect was `2.84e-11` relative.

Regular-grid FFT is selected for the initial uniform-cube tree scope (see
the integration note below). At the isolated-experiment checkpoint, remaining M9 gates
were production integration, suitable preconditioning, dense comparisons,
observable tolerance sweeps and an integrated N3 benchmark. These are now
implemented and validated as described below. The runtime budget was agreed on
2026-10-20, completing M9. M8 retains its separate gates.

OpenEMS is not a mandatory M9 exit criterion. The primary M9 correctness check
is agreement with the dense operator on identical meshes, well below mesh
error, together with analytical and conservation tests. An independent
openEMS tree comparison can add physical-model evidence after acceleration
is established; it does not replace operator-equivalence tests or close the
remaining M8 N1/N2 independent-comparison gates.

#### M9 production integration and validation

`nec2c -s fft` now applies the volume block through the uniform-cube FFT
backend and stores only the NEC conductor matrix and the two cross blocks.
The plane-wave and driven wire/patch/mixed paths share this representation;
they do not allocate a quadratic volume matrix. Dense remains the default.
The adapter rejects incompatible lattices/materials and enforces the initial
512 MiB accounted solver cap. `-r` and `-n` expose GMRES residual and iteration
controls; failed solves stop before solved observables are printed.

The small-tree backend decision, supported mesh policy, implementation,
validation commands and N3 antenna fixture are recorded in
[M9-FFT-INTEGRATION.md](M9-FFT-INTEGRATION.md). Regular-grid FFT is selected
for this uniform-cube tree scope; sparse embedding and an occupied-octree
H-matrix alternative are compared there without claiming unmeasured
hierarchical ranks or performance. The initial end-to-end volume, wire,
patch and mixed dense comparisons and hybrid sanitizer checks passed.
The eight-cell overlap comparisons now pass at requested residual `1e-13`,
retaining `1e-8` relative / `1e-12` absolute comparison thresholds. All twelve
regression programs pass. Wire and mixed tolerance sweeps at `1e-10` meet
`1e-4` ohm feed-impedance and `1e-6` scaled complex-pattern goals against
dense references. A translated-contact quadrature inconsistency discovered
by the patch comparison was corrected and given independent regression tests.

The integrated N3 tree plus 11-segment dipole converged in 45 iterations to
`7.75e-11` true residual at requested `1e-10`, using peak RSS 140,064 KiB
(136.8 MiB) and accounted solver storage 141,975,080 bytes. Total wall time
was 1,224.76 s (20 min 25 s), including construction, solving and output;
CPU time was 86.46 s. Results, inputs and sampled complex patterns are retained
with the integration note. On 2026-10-20 the user agreed to a 30-minute-per-frequency
runtime acceptance budget for this target at requested residual `1e-10`, within
the initial 512 MiB accounted solver cap. This is a benchmark acceptance target,
not an enforced timeout or a limit on verification runs; the user permits
verification runs to take as long as needed. Final kernel/FFT/hybrid
AddressSanitizer and UndefinedBehaviorSanitizer
checks and the distribution build passed. The new gap-scaling and translated-
contact tests passed against the retained independent kernel suite. All
M9 exit criteria are satisfied, and M9 is complete. A multi-frequency N3
benchmark is optional additional coverage, not an M9 exit requirement. The
recorded resource acceptance applies to the specified 28 MHz case, not an
unmeasured frequency range. M8's separate gates are unchanged.

### M10: dielectric above perfect ground with a positive gap

**Result:** Wire, patch, and lossless dielectric models separated from a planar
perfect-ground interface by a documented positive clearance use one consistent
PEC-image background formulation.

**Included:** `GN 1`; all conductor/volume and volume/volume PEC reflected terms,
ground-aware excitation and near/far fields, and lossless power/reciprocity checks.

**Excluded:** real ground (`GN 0/2`, assigned to M13), dielectric contact with
ground (M11), cliff ground, and radial screens. Existing preliminary GN 2
support remains available but is outside M10 acceptance.

**Exit criteria:**

- the supported PEC `GN 1` subset and clearance policy are stated unambiguously;
- removing dielectric recovers the existing PEC ground result;
- moving perfect ground far away approaches the free-space dielectric result;
- independent PEC analytical/image cases pass;
- at least two independently refined openEMS PEC ground-gap models agree; and
- lossless power and reciprocity checks include the PEC background consistently.

**Artifacts:** PEC ground formulation note, implementation, canonical fixtures,
openEMS comparisons, and revised user constraints. Real-ground evidence is
retained in the [M13 handoff](M13-LOSS-AND-REAL-GROUND.md).

#### M10 acceptance completed (2026-10-20)

The revised PEC-only scope is complete. Independent image/analytical tests,
legacy-vacuum recovery, the coupled remote-ground limit, port reciprocity and
lossless power checks pass. All 18 regression groups pass. Volume refinement
from 27 through 64 to 125 cells changes dielectric-induced feed impedance by
0.00620 then 0.00288 ohm; wire refinement changes it by 0.0446 ohm.

The two independently refined OpenEMS PEC discretizations agree within
0.2631 ohm across 280/300/320 MHz. NEC differs from the finer OpenEMS delta Z
by 0.1894/0.3982/0.8218 ohm, below the agreed 2 ohm allowance. Additional
larger-domain and shorter-feed pairs all reach -60 dB decay; their maximum
delta-Z changes are 0.00574 and 0.00630 ohm respectively. The acceptance is
for dielectric-induced impedance changes, not absolute finite-feed versus
delta-gap feed impedance. It applies to the documented positive-clearance
PEC subset; dielectric/ground contact remains M11 and real ground remains M13.

Evidence: [acceptance summary](../tests/m10-pec-acceptance.json),
[sensitivity audit](../tests/m10-pec-sensitivity-review.json), and
[formulation record](M10-GROUND-FORMULATION.md). The progress entries below
are historical checkpoints; their pending M10 gates are closed by this record.

#### M10 first increment (2026-10-20)

The original scope included existing NEC2 perfect and real-ground
abstractions, with PEC analytical checks first. The later scope decision above
assigns real-ground work to M13. The formulation, GN mapping,
legacy Sommerfeld/patch compatibility issue and remaining integration gates
are recorded in [M10-GROUND-FORMULATION.md](M10-GROUND-FORMULATION.md).
`vie_ground.c` implements a positive-gap PEC cuboid volume kernel using
signed electric-current images. Independent static/dynamic dipole quadrature,
PEC boundary, weighted reciprocity, self-image, remote-ground and invalid-input
tests pass. At this first checkpoint, production ground decks remained rejected
until coupling, excitation, fields and power used the same background. This increment does
not close M10 or its independent OpenEMS comparison requirements.

#### M10 PEC deck integration (2026-10-20)

The dense and dense-operator GMRES paths now support positive-gap `GN 1`
models with wires, patches and dielectric volumes. The conductor block keeps
NEC's PEC state, both cross blocks and the volume block include current
images, plane-wave excitation includes the reflected background, and near/
far fields include the same images. Radiated/scattered power integrates the
upper hemisphere; volume extinction uses the total background electric field.
The initial scope requires patch centroid height greater than `sqrt(area)`,
wire endpoint/radius clearance and complete dielectric boxes above z=0,
upper-hemisphere incident/pattern directions and nonnegative near-field z.
Elliptical EX types 2/3 require a positive explicit F6 axial ratio. At this
PEC integration checkpoint, `-s fft`, real ground, screens, cliff ground and
contact were rejected; subsequent GN 2 work is preliminary M13 support.

`tests/dielectric_ground.py` compares volume, wire and patch solutions with
explicit mirrored free-space models, checks vacuum PEC compatibility for
wire/patch/mixed models, dense/GMRES agreement, ground boundary conditions,
power, invalid requests and frequency-sweep state restoration. The runnable
fixture is `tests/dielectric/m10-pec-wire.nec`. See the
[formulation record](M10-GROUND-FORMULATION.md) for verification details,
limits of the initial patch-clearance abstraction and outstanding gates.
M10 remains in progress for PEC power/refinement evidence and independently
refined PEC OpenEMS comparisons. Real-ground consistency and implementation
are now M13 work and do not gate this PEC milestone.

#### M10 finite-ground reference operator and laptop staging (2026-10-20)

The isolated `vie_halfspace` module now evaluates reflected electric and
magnetic point/cuboid fields for a passive dielectric/conducting half-space,
with exact pulse-cell spectral averages and explicit quadrature failure.
Fresnel, quasistatic, PEC image, finite-loss independent axial reference,
lossless branch and reciprocity checks pass. The subsequent GN 2 wire
integration below connects this module to the deck solver. GN 0 and
real-ground patches remain separate follow-up work.

A portable, externally time-limited OpenEMS runner and SSH instructions are
provided in [M10-OPENEMS-LAPTOP.md](M10-OPENEMS-LAPTOP.md). Development remains
on the rxtx cloud VM; OpenEMS is on the user's laptop. Begin with geometry
preparation and a bounded timing pilot to assess battery cost. The original
pilot has now returned (1.22394 ps, 5,000 steps in 20.29 s, clipped source);
its review and the user-authorized 5 mm radius revision are recorded in the
laptop report. The revised deterministic mesh bounds grading to 2:1 and
starts with another 30-second pilot. No passing independent ground comparison
is claimed.
See the [formulation record](M10-GROUND-FORMULATION.md) for limitations and
remaining gates. M10 remains in progress.

### M11: dielectric/perfect-ground contact

**Result:** A lossless dielectric region may touch the PEC ground interface,
so idealized tree and building models can be grounded without an artificial
air gap while retaining the lossless proof of concept.

**Included:** contact-cell self/near interactions and the air/dielectric/PEC
junction for perfect ground (`GN 1`).

**Excluded:** real ground and material loss (M13), and buried dielectric volume
unless added through a separate reviewed scope change.

**Exit criteria:**

- positive-gap models converge to the contact solution as the gap closes;
- contact-specific quadrature passes an independent high-accuracy oracle;
- canonical contact configurations agree with openEMS under mesh refinement;
- no nonphysical power generation or discontinuous impedance jump appears at
  the transition to contact; and
- the grounded idealized-tree model completes within the agreed accuracy and
  resource budget.

**Artifacts:** Contact derivation, quadrature tests, contact fixtures, and
grounded-tree validation report.

The following development checkpoints are historical; the completion record
below closes their pending gates.

#### M11 first increment (2026-10-20)

The PEC image kernel now permits exact face contact at zero requested minimum
clearance, reusing the existing touching-cuboid quadrature. Independent
face-oracle, resolved-gap approach and weighted-reciprocity tests precede
any relaxation of deck-level restrictions. Exploratory 0.1 mm and 0.01 mm
gaps exposed quadrature nonconvergence; resolving that near-contact regime
is required before claiming continuous contact support. The remaining M11 work,
including the grounded-tree backend/resource gate, is tracked in
[M11-PEC-CONTACT.md](M11-PEC-CONTACT.md). M11 is not complete.

#### M11 contact/FFT integration (2026-10-20, in progress)

The closest-corner quadrature extension resolves the originally failing tiny
gaps without relaxing tolerances. Dynamic kernel and solved-current tests
approach contact through 1 micrometre; an independent analytic Coulomb-face
reference checks static contact and near-contact tensors. PEC dielectric
contact is now accepted by the deck layer; real-ground contact, buried volume,
and conductor/dielectric contact retain their restrictions.

PEC FFT uses a reflected-source embedding with only physical cells as unknowns.
Small volume, wire, patch and mixed contact cases agree with dense solves and
pass power checks. The 216-cell grounded skeleton completes with true residual
5.43e-11 and 28,064 KiB peak RSS; the full 9,432-cell tree remains under test.
All four returned OpenEMS pillar verification runs reach -60 dB decay. The
25 mm NEC mesh has a 2.03 ohm complex-vector discrepancy at 320 MHz against
fine OpenEMS, but its separate resistance/reactance errors (1.06/1.73 ohm)
both pass the documented component-wise greater-of-2%-or-2-ohm criterion.
The vector magnitude is retained as an additional diagnostic, not a stricter
replacement acceptance gate. Finer volume/wire and full-tree checks remain
pending. M11 is not complete.

#### M11 validation checkpoint and remaining runtime decision (2026-10-20)

Contact-kernel, explicit-image, boundary, reciprocity, power, native-refinement,
and independent OpenEMS component-wise impedance gates pass. All 19 regression
groups pass. The full 9,432-cell grounded tree converges to true residual
9.07e-11 in 45 iterations, with 234.8 MiB peak RSS and a 0.2808% power defect.
Its elapsed time is 34.6 minutes (155.42 s process CPU), exceeding the nominal
30-minute target. A proposed 40-minute target for grounded trees versus further
optimization awaits user selection. The 512 MiB cap and accuracy gates remain
unchanged. M11 remains in progress until this runtime-target issue is resolved.
Evidence is retained in [m11-acceptance.json](../tests/m11-acceptance.json) and
[m11-pec-tree.json](../tests/m11-pec-tree.json).

#### M11 acceptance completed (2026-10-21)

The user accepted a 40-minute per-frequency grounded-tree target and deferred
further optimization to future revisions after experience with the models.
The measured 2075.89 s (34.60 minutes) passes this 2400 s target. The M9
free-space 30-minute target, 512 MiB accounted solver cap, and all accuracy
requirements remain unchanged.

All M11 exit criteria are satisfied: independent contact quadrature and gap
convergence, PEC boundaries and reciprocity, power, refined OpenEMS
component-wise impedance comparisons, native refinement, and the full grounded
tree accuracy/resource check. All 19 regression groups pass. M11 is COMPLETE;
M12 remains deferred, and M13 is the next planned milestone. The retained
acceptance and tree JSON records include the agreed runtime decision.

### M12: touching and embedded conductors (deferred to a later revision)

**Revision decision (2026-10-20):** intentionally deferred at the user's request.
Conductor/dielectric contact and penetration are not needed for the current
simulation scenarios. The proposal below is retained for a future revision;
it is not a current-revision exit gate or a prerequisite for M13. Until then,
wire/patch conductors and dielectric volumes must retain supported positive
clearance; M11's ideal PEC-ground contact is a separate boundary case.

**Result:** Supported wires or PEC surfaces may touch, be enclosed by, or be
embedded in dielectric material, enabling insulation, housings, and radomes
without an artificial clearance.

**Included:** near-singular conductor/volume quadrature, conductor exclusion
geometry, touching-interface conditions, and separated-to-contact convergence
with lossless dielectrics, ideal conductors and free space or PEC ground.

**Dependency note:** dielectric loss and real ground are M13 capabilities.
Combining them with future M12 contact requires additional validation.

**Exit criteria:**

- a separated radome case converges before contact cases are enabled;
- dielectric approach-to-wire/patch tests converge without arbitrary result
  jumps;
- insulation and enclosure benchmarks agree with analytical cases where
  available and openEMS elsewhere;
- conductor current, dielectric current, and total power remain consistent;
  and
- invalid overlap topology is rejected deterministically.

**Artifacts:** Contact formulation, topology rules, embedded-model fixtures,
openEMS comparisons, and Antenna Turtle examples.

### M13: loss, real ground and advanced material behavior

**Status: COMPLETE (2026-10-21)** for passive scalar materials and GN 2
positive-gap wire/volume real ground. All 24 regression groups pass; the
M13-05 exit-criterion review closes the agreed scope.

**Result:** The solver extends the lossless proof of concept through M11 to real
ground and the material behavior required for physically useful vegetation
and building models beyond a real constant `epsilon_r`.

**Included:** complex permittivity or conductivity and consistent loss accounting
in dielectric regions and the ground background; documented real-ground GN
and conductor/contact subsets; frequency-dependent material tables; and multiple
touching dielectric materials. Anisotropy and magnetic media require separate
scope decisions.

**Dependency boundary:** deferred M12 is not a prerequisite. Initial loss and
real-ground scenarios retain positive conductor/dielectric clearance. If M11
PEC-ground contact is extended to real ground, that contact scope and its
verification must be stated explicitly rather than inherited automatically.

**Accepted scope:** GN 2 positive-gap wire/volume coupling, reflected E/H
kernels, substrate-power accounting, analytical limits and independently refined
lossy OpenEMS comparisons. GN 0, real-ground patches, ground contact/burial,
screens/cliff and ground FFT remain outside this completed scope. Scalar loss,
frequency tables, touching materials and tree-material reporting are completed
by the increments below; anisotropy and magnetic media remain excluded.
See [M13-LOSS-AND-REAL-GROUND.md](M13-LOSS-AND-REAL-GROUND.md).

**M13-01 COMPLETE (2026-10-21):** constant isotropic complex
permittivity/conductivity in free space and over PEC. Anisotropy is excluded
from the current scope. Independent OpenEMS decay, material, mesh/domain and impedance checks pass,
as do native analytical/refinement checks and all 20 regression groups.
Later M13 increments remain open. See
[M13-01-SCALAR-LOSS.md](M13-01-SCALAR-LOSS.md).

**M13-02 COMPLETE (2026-10-21):** frequency-dependent scalar material tables,
linear frequency interpolation, explicit range checks and a conservative
fixed-mesh bound. Native checks, all 21 regressions and independent OpenEMS
witness acceptance pass. Maximum component error is 1.33655 ohm; all 12
verification runs reach -60 dB decay. M13 overall remains in progress. See
[M13-02-MATERIAL-TABLES.md](M13-02-MATERIAL-TABLES.md).

**M13-03 COMPLETE (2026-10-21):** touching scalar dielectric materials and
cell-local complex permittivity in FFT. Independent analytical coated-sphere
and interface checks, all 22 regressions, 11 native benchmark/refinement
comparisons and 12 OpenEMS verification runs pass. Largest component error is
1.97370 ohm; the joint refined hardest case gives 1.69712 ohm. M13 remains open
for real-ground work. See
[M13-03-TOUCHING-MATERIALS.md](M13-03-TOUCHING-MATERIALS.md).

**M13-04 COMPLETE (2026-10-21):** combined scalar material loss and GN 2
wire/volume real-ground acceptance with positive ground/conductor gaps. All 23
regression groups, six native comparisons and 12 independent OpenEMS verification
runs pass. Worst component discrepancy is 0.829761 ohm. GN 0, real-ground patches
and contact/burial remain excluded. The material/uncertainty guidance review
is completed by M13-05 below. See
[M13-04-REAL-GROUND-WIRES.md](M13-04-REAL-GROUND-WIRES.md).

**M13-05 COMPLETE (2026-10-21):** production-tree material provenance,
geometry/ground assumptions, sampled material-sensitivity reporting and separate
nominal mesh refinement. Five 156-cell material witnesses, a 1,248-cell refinement,
the explicit localized-source remote-ground E/H limit and all 24 regression
groups pass. No further OpenEMS run was needed. This closes M13 for the agreed
scope; it does not claim empirical vegetation material validation or full-size
lossy/real-ground tree performance. See
[M13-05-TREE-MATERIAL-GUIDANCE.md](M13-05-TREE-MATERIAL-GUIDANCE.md).

**Exit criteria:**

- passive sign conventions are verified against absorbed-power calculations;
- zero-loss dielectric limits and the appropriate PEC, matched-background and
  remote-ground limits reproduce the implemented lossless reference results
  through M11, and M12 only if subsequently implemented (zero substrate
  conductivity alone does not turn finite-permittivity ground into PEC);
- the supported real-ground GN, conductor and contact subsets are explicit,
  preserve supported legacy limits and pass ground-aware reciprocity/power checks;
- dispersive sweeps use documented interpolation and reject silent
  extrapolation;
- layered or concentric multi-material analytical cases pass;
- lossy openEMS comparisons pass independent refinement; and
- tree material assumptions and their uncertainty are reported with results.

**Artifacts:** Versioned material-format extension, dielectric/ground loss-power
implementation, real-ground formulation and supported-mode policy, analytical
and openEMS fixtures, and updated production-tree guidance.

## Deferred design questions

M4 resolved the cuboid singular integration and oracle question in
`VIE-CUBOID-KERNEL.md`. M6 resolved the initial joint NEC/JVIE scaling in
`DIELECTRIC-WIRE-M6.md`. The remaining questions must be resolved before the
applicable later milestone ships:

1. What minimum conductor clearance is justified by measured quadrature error
   rather than the conservative initial rule?
2. At production tree scale, is regular-grid FFT, sparse FFT, or a hierarchical
   method the best acceleration strategy?
3. M13-04 resolves the initial planar positive-gap Sommerfeld representation
   and its validation subset; contact/burial and large real-ground volumes
   require later accuracy/resource work. See M13-05 for the accepted limits.

Each answer should be recorded as a short decision note with equations,
references, numerical evidence, and the validation fixtures it affects.

## References

- A. G. Polimeridis, J. F. Villena, L. Daniel, and J. K. White, “Stable
  FFT-JVIE solvers for fast analysis of highly inhomogeneous dielectric
  objects,” *Journal of Computational Physics*, 269, 280–296, 2014,
  <https://doi.org/10.1016/j.jcp.2014.03.026>.
- M. M. Botha, “Solving the volume integral equations of electromagnetic
  scattering,” *Journal of Computational Physics*, 218(1), 141–158, 2006,
  <https://doi.org/10.1016/j.jcp.2006.02.004>.
- B. T. Draine and P. J. Flatau, “User Guide for the Discrete Dipole
  Approximation Code DDSCAT 7.3,” 2013,
  <https://arxiv.org/abs/1305.6497>.
- B. T. Draine and P. J. Flatau, “Discrete-dipole approximation for scattering
  calculations,” *Journal of the Optical Society of America A*, 11(4),
  1491–1499, 1994, <https://doi.org/10.1364/JOSAA.11.001491>.
- S. P. Groth, VoxScatter, Galerkin VIE and DDA reference implementation,
  <https://github.com/samuelpgroth/VoxScatter>.
- A. G. Polimeridis et al., MARIE, open-source integral-equation suite,
  <https://github.com/thanospol/MARIE>.
- P. R. Wiecha, “pyGDM—A Python toolkit for full-field electro-dynamical
  simulations and evolutionary optimization of nanostructures,” *Computer
  Physics Communications*, 233, 167–192, 2018,
  <https://doi.org/10.1016/j.cpc.2018.06.017>.
- openEMS numerical-method documentation,
  <https://docs.openems.de/en/latest/concepts/numerical_method.html>.
- openEMS discrete-material documentation,
  <https://docs.openems.de/en/latest/concepts/disc_material.html>.

#### Preliminary M13 GN 2 work, developed during M10 (2026-10-20)

The user selected GN 2 wires first, with GN 0 and real-ground patches as
separate follow-up work. The dense and GMRES paths now combine native NEC
Sommerfeld/Norton wire interactions with spectral finite-cell volume and
wire/volume reflected interactions. Ground-aware excitation, near/far fields
and an independent downward interface-flux power calculation are connected.
Driven-source and plane-wave power checks include energy delivered into the
substrate. Negative-conductivity GN semantics and frequency-sweep state are
preserved. The runnable fixture is `tests/dielectric/m10-real-wire.nec`;
[the formulation record](M10-GROUND-FORMULATION.md) documents approximations,
verification and limitations. Under the revised lossless milestone boundary,
this is preliminary M13 evidence. M10 remains in progress only for its PEC
verification gates; GN 0 and real-ground patch work do not block M10 or M11; M12 is separately deferred.

The completed 64-cell, 300 MHz GN 2 wire benchmark differs from the completed
coarse OpenEMS dielectric impedance perturbation by 0.27573 ohm, below the
agreed 2-ohm checkpoint allowance. Its power defect is 0.0876% and true GMRES
residual 2.30e-11. Results are in `tests/m10-nec-real-300.json`. This single
frequency/mesh comparison does not replace the pending fine OpenEMS result
or resolve source, absolute feed-impedance and boundary convergence gates.
