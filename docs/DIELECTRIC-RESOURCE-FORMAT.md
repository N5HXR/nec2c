# NEC2C dielectric resource format, version 1

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

Exactly one version header must be the first content record, and version `1`
is currently supported. At least one `MATERIAL`, `REGION`, and `BOX` is
required. `END` is mandatory, and only blank lines or comments may follow it.
Record order after the header has no semantic meaning, so forward references
are allowed. IDs need not be contiguous but must be unique within their record
type.

`MATERIAL` requires a real relative permittivity `ER >= 1`. A value of one is
valid and will represent vacuum in later solver milestones. Each `REGION`
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
q = k0 * sqrt(ER) * longest_cell_edge <= MAX_KH.
```

`MAX_KH <= 0.5` is the recommended range. Values through `1.0` are accepted
with an accuracy warning, while larger values are rejected. Each physical box
is divided uniformly on each axis; the longest-axis division count is then
increased until `MAX_ASPECT` is met. Cells are emitted in canonical box order,
then increasing x, y, and z subdivision index. `ER=1` boxes remain in physical
volume and region reports but do not produce volume unknowns.

Overlapping box interiors are rejected. Same-material boxes may touch, but
touching boxes with different material IDs are rejected in format version 1.
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
resource and cost reports.  A geometry-free model with an `EX` plane wave
invokes M5's dense free-space volume solver at `RP`, `XQ`, or `EN`.  A
free-space conductor model with either an ordinary wire voltage source or a
type 1--3 plane-wave `EX` card invokes the self-consistent wire/patch/volume
solver at `XQ` or `EN`; see `DIELECTRIC-WIRE-M6.md` and
`DIELECTRIC-PATCH-M7.md`. Other explicitly excluded features fail rather than
producing a silently uncoupled legacy result.
