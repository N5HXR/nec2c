# Dielectric volume support

## Status

This document specifies a staged extension to the N5HXR `nec2c` fork for
lossless dielectric volumes.  It is a design and validation contract, not a
claim that dielectric support is already implemented.

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
| M8 | Validated free-space dielectric release | NOT STARTED |
| M9 | Production-scale free-space tree solver | NOT STARTED |
| M10 | Dielectric above ground with a positive gap | NOT STARTED |
| M11 | Dielectric/ground contact | NOT STARTED |
| M12 | Touching and embedded conductors | NOT STARTED |
| M13 | Loss and advanced material behavior | NOT STARTED |

The status table is the project-level progress tracker.  Test reports and
benchmark artifacts provide the evidence behind it.

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

### M10: dielectric above ground with a positive gap

**Result:** Wire, patch, and dielectric models separated from a planar ground
interface by a documented positive clearance use one consistent
layered-background formulation.

**Included:** all conductor/volume and volume/volume reflected or Sommerfeld
terms needed by the selected initial ground mode, plus ground-aware far fields.

**Excluded:** dielectric contact with ground, cliff ground, and radial screens
unless explicitly selected by the milestone decision note.

**Exit criteria:**

- the supported `GN` subset is stated unambiguously;
- removing dielectric recovers the existing supported ground result;
- moving ground far away approaches the free-space dielectric result;
- layered analytical or published canonical cases pass;
- at least two independently refined openEMS ground-gap models agree; and
- power and reciprocity checks include the ground background consistently.

**Artifacts:** Ground formulation note, implementation, canonical fixtures,
openEMS comparisons, and revised user constraints.

### M11: dielectric/ground contact

**Result:** A dielectric region may touch the supported ground interface, so
tree and building models can be physically grounded without an artificial air
gap.

**Included:** contact-cell self/near interactions and the air/dielectric/ground
material junction for the supported ground model.

**Excluded:** buried dielectric volume unless added through a separate reviewed
scope change.

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

### M12: touching and embedded conductors

**Result:** Supported wires or PEC surfaces may touch, be enclosed by, or be
embedded in dielectric material, enabling insulation, housings, and radomes
without an artificial clearance.

**Included:** near-singular conductor/volume quadrature, conductor exclusion
geometry, touching-interface conditions, and separated-to-contact convergence.

**Excluded:** material loss unless M13 is also complete.

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

### M13: loss and advanced material behavior

**Result:** The solver supports the material behavior required for physically
useful vegetation and building models beyond a real constant `epsilon_r`.

**Included in order:** complex permittivity or conductivity, loss accounting,
frequency-dependent material tables, and multiple touching dielectric
materials.  Anisotropy and magnetic media require separate scope decisions.

**Exit criteria:**

- passive sign conventions are verified against absorbed-power calculations;
- zero-loss limits reproduce M8 through M12 results;
- dispersive sweeps use documented interpolation and reject silent
  extrapolation;
- layered or concentric multi-material analytical cases pass;
- lossy openEMS comparisons pass independent refinement; and
- tree material assumptions and their uncertainty are reported with results.

**Artifacts:** Versioned material-format extension, loss/power implementation,
analytical and openEMS fixtures, and updated production-tree guidance.

## Deferred design questions

M4 resolved the cuboid singular integration and oracle question in
`VIE-CUBOID-KERNEL.md`. M6 resolved the initial joint NEC/JVIE scaling in
`DIELECTRIC-WIRE-M6.md`. The remaining questions must be resolved before the
applicable later milestone ships:

1. What minimum conductor clearance is justified by measured quadrature error
   rather than the conservative initial rule?
2. At production tree scale, is regular-grid FFT, sparse FFT, or a hierarchical
   method the best acceleration strategy?
3. Which layered-medium Green-function representation is accurate and
   efficient enough for volume cells at and near Sommerfeld ground?

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
