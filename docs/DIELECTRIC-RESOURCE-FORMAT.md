# NEC2C dielectric resource format, versions 1, 2 and 3

Revision 1 is complete; see the [user guide](DIELECTRIC-REVISION-1.md) and
[canonical examples](../examples/dielectric-r1/README.md) for current execution
and ground/solver scope.

This is the producer-facing contract for Antenna Turtle and other geometry
generators. Milestones M2/M3 parse, validate, and mesh this format and report
the projected dense-solver cost. Milestone M5 also solves free-space,
dielectric-only plane-wave models; see `DIELECTRIC-SOLVER-M5.md`.

## Selecting a resource

Place one `DF` card in the NEC control-card section after `GE`:

```text
GE 0
DF "dielectric/tree model.diel"
EN
```

A relative path is interpreted relative to the NEC deck, regardless of the
process working directory. Double quotes are required when the path contains
spaces. Absolute paths work but are reported as nonportable. There may be at
most one `DF` per NEC job; `NX` starts a new job and releases the prior model.
The complete card, including `DF`, must fit within NEC2C's 132-character input
line limit.

## File grammar

The resource is line-oriented text. Keywords are case-insensitive, whitespace
separates fields, and `#` begins a comment outside a quoted token. Identifiers
are positive decimal integers. Numbers use C-locale syntax and must be finite.

```text
NEC2C-DIELECTRIC 1

MATERIAL 1 ER 3.0
REGION 101 MATERIAL 1
BOX 101 -0.5 -0.5 0.1  0.5 0.5 10.0

MESH MAX_KH 0.5 MAX_CELLS 1000 MAX_ASPECT 2.0
END
```

The records are:

```text
NEC2C-DIELECTRIC <version>
MATERIAL <material-id> ER <relative-permittivity>
REGION <region-id> MATERIAL <material-id>
BOX <region-id> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax>
MESH [MAX_KH <q>] [MAX_CELLS <count>] [MAX_ASPECT <ratio>]
END
```

Exactly one version header must be the first content record, and versions `1`, `2`
and `3` are supported. At least one `MATERIAL`, `REGION`, and `BOX` is
required. `END` is mandatory, and only blank lines or comments may follow it.
Record order after the header has no semantic meaning, so forward references
are allowed. IDs need not be contiguous but must be unique within their record
type.

`MATERIAL` requires a real relative permittivity `ER >= 1`. Version 2 optionally
adds `EI <nonnegative-loss>` or `SIGMA <nonnegative-S/m>` after ER; these are
mutually exclusive. The exp(+i omega t) convention gives relative permittivity
`ER - i*(EI + SIGMA/(omega*epsilon0))`. Omitted loss is zero; version 1 rejects
loss fields. ER=1 represents vacuum only when loss is zero. See
[M13-01-SCALAR-LOSS.md](M13-01-SCALAR-LOSS.md) for the implementation and
acceptance status. Each `REGION`
must reference a declared material. Each `BOX` must reference a declared
region, uses meters in NEC's global Cartesian coordinates, and requires
`xmin < xmax`, `ymin < ymax`, and `zmin < zmax`.

There may be zero or one `MESH` record. Options may appear in any order and
may be omitted. An option may not occur twice. Version-1 defaults are:

| Option | Default | Parser constraint |
| --- | ---: | --- |
| `MAX_KH` | `0.5` | finite and greater than zero |
| `MAX_CELLS` | `1000` | positive decimal integer representable by `size_t` |
| `MAX_ASPECT` | `2.0` | finite and at least one |

The parser sorts materials and regions by ID and boxes by region and bounds,
producing the same canonical summary and cell order for equivalent record
orderings.

## M3 meshing and admissibility rules

Meshing uses the highest frequency in the NEC `FR` sweep. For every active
cell, M3 enforces

```text
q = k0 * sqrt(abs(epsilon_relative)) * longest_cell_edge <= MAX_KH.
```

`MAX_KH <= 0.5` is the recommended range. Values through `1.0` are accepted
with an accuracy warning, while larger values are rejected. Each physical box
is divided uniformly on each axis; the longest-axis division count is then
increased until `MAX_ASPECT` is met. Cells are emitted in canonical box order,
then increasing x, y, and z subdivision index. Constant `ER=1` boxes with zero
loss remain in physical volume and region reports but do not produce volume
unknowns. `ER=1` with positive EI or conductivity remains an active material.

Overlapping box interiors are rejected. Boxes of equal or different scalar
materials may share faces, partial faces, edges or corners (M13-03). This
extends versions 1/2/3 without new syntax. Use identical coordinates for a
shared boundary; no overlap tolerance or coordinate snapping is applied.
See [M13-03-TOUCHING-MATERIALS.md](M13-03-TOUCHING-MATERIALS.md) for acceptance
status and FFT lattice restrictions.
The report verifies numerical-volume conservation, reports per-box division
counts and per-region bounding-box diameter, and prints the electrical region
diameter `k0 * sqrt(ER) * D`.

Every dielectric cell contributes three Cartesian volume-current unknowns.
The preflight report combines those with one unknown per wire segment and two
per conducting patch, then estimates dense complex matrix bytes as
`16 * total_unknowns^2`. All counts and byte calculations are checked for
overflow, and `MAX_CELLS` is enforced before the cell array is allocated.

Wire clearance uses exact segment-to-axis-aligned-box distance and subtracts
the wire radius. The required gap is the larger of the cell's longest edge and
the wire segment length. NEC2C's processed patch representation retains patch
center and area but not the original polygon extents, so M3 uses an
equivalent-area radius as a patch-clearance proxy and the larger of cell edge
and square root of patch area as its required gap. M7 retains this conservative
proxy because touching conductors remain out of scope; exact polygon-aware
clearance is deferred to the M12 contact work.

## Current execution behavior

`DF` followed by `EN` without excitation validates, meshes, and prints the
resource and cost reports. With excitation, `RP`, `XQ`, `NE`/`NH` and `EN` execute
supported dielectric solves/output requests. Geometry-free models use plane-wave
EX types 1–3; conductor models accept ordinary wire voltage sources or plane-wave
EX types 1–3. Successive field requests can execute the same frequency again.

Revision 1 supports free space, PEC (`GN 1`, including dielectric-ground contact),
and GN 2 real ground with wires and strictly positive ground clearances. Dense
is the default; assembled GMRES and common-cubic-lattice FFT-backed GMRES are
selected with `-s gmres` or `-s fft`. FFT supports cell-local scalar materials in
free space/PEC and is excluded for GN 2. Conductor/dielectric contact, LD/NT/TL,
NEC symmetry and other explicitly excluded modes fail rather than returning a
silently uncoupled result. See the [revision-1 guide](DIELECTRIC-REVISION-1.md)
for excitation, ground, field-output, clearance and solver constraints. Earlier
M5/M6/M7 documents retain historical formulation and validation details.

## Version 3 frequency tables (M13-02)

Version 3 retains constant materials and adds `MATERIAL <id> TABLE` with at
least two `SAMPLE <id> FREQ_MHZ <MHz> ER <ER> EI|SIGMA <loss>` rows. ER and
the selected loss parameter interpolate linearly in frequency; one table must
use a consistent loss representation. Out-of-range sweeps are rejected before
solving. Table samples may be unordered and forward-referenced, but must have
unique positive frequencies, passive finite values, and nonzero susceptibility.
A conservative table-wide material bound sizes the fixed mesh. See
[M13-02-MATERIAL-TABLES.md](M13-02-MATERIAL-TABLES.md) for the complete contract
and current acceptance status.
