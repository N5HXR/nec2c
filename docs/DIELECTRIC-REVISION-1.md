# Dielectric volume support: revision 1

Revision 1 is COMPLETE (2026-10-21): M1–M11 and M13 are accepted; M12 was
intentionally deferred. This is a **feature revision**, not a renumbering of
upstream nec2c 1.3. Resource-format versions 1–3 are separate file-grammar
versions, all supported by this feature revision. The supported model and its
exclusions are stated below.

NEC2C now couples volumetric dielectric polarization currents to NEC wire and
conducting-patch currents in one self-consistent frequency-domain system.
This is a volume/surface integral-equation extension to NEC2. Dielectrics use
three cell-average Cartesian **volume** current components per cell. Existing
NEC conductor representations remain in use. A deck without `DF` follows the
legacy NEC path.

## Start here

From a Git checkout with a C/C++ toolchain, Python >=3.6, make, Autoconf
and Automake:

```sh
autoreconf -fi
./configure
make
python3 examples/dielectric-r1/run.py --output-dir /tmp/nec2c-r1-examples
```

The solver links libm and includes its FFT implementation; it does not require
OpenEMS, Node.js, NumPy or an external FFT package. The existing configure
script checks C/C++ compilers and Python >=3.6. Python 3 is used by tests,
example runners and optional reporting tools. OpenEMS is an optional independent
validation environment, not a runtime dependency. From a generated source
archive, start with `./configure` and `make`. Unix builds are exercised here;
the inherited Windows/MSYS2 build path has not received a revision-1 validation
run. See the root README for that platform's existing build instructions.

The example runner needs only Python's standard library and runs four small
models. It writes raw `.out`/`.log` files plus `summary.json` to the output
directory. It checks completion, residual and power diagnostics; it does not
require copying or editing input files. Relative `DF` paths resolve beside the
NEC deck, so these commands also work from another working directory when the
binary and deck paths are absolute.

## Supported scope

| Capability | Free space | PEC (`GN 1`) | Real ground (`GN 2`) |
| --- | --- | --- | --- |
| Wire/dielectric coupling | Yes | Yes, positive wire-ground clearance | Yes, wires with positive ground clearance |
| Conducting-patch/dielectric coupling | Yes | Yes, conservative patch clearance | Excluded |
| Dielectric-only plane waves | Yes | Yes | Yes, positive ground clearance |
| Lossless or lossy passive scalar materials | Yes | Yes | Yes |
| Frequency tables and touching dielectric materials | Yes | Yes | Yes |
| Dielectric contact with ground | Not applicable | Yes, flat z=0 contact | Excluded; all boxes must have zmin>0 |
| Dense / assembled GMRES | Yes | Yes | Yes |
| Regular-grid FFT-backed GMRES | Common cubic lattice | Common cubic lattice and PEC image alignment | Excluded |

All material laws are isotropic and nonmagnetic, with ER>=1 and passive loss.
Dielectric interiors cannot overlap. Touching dielectric regions are supported;
conductor/dielectric contact, penetration and burial are not. FFT supports
cell-local material variation, including touching materials, when the mesh
satisfies its lattice rules. The ground is an infinite planar background, not
an automatically meshed soil volume. No GN 0/Fresnel approximation, radial screen
or cliff-ground dielectric coupling is provided.

Ordinary voltage sources (`EX 0`) and plane-wave types 1–3 are supported within
the established geometry/excitation restrictions. Wire loading (`LD`), networks
(`NT`), transmission lines (`TL`), discontinuity sources and NEC symmetry
reduction are not integrated with active dielectric solving. Expand geometry
explicitly where needed. Unsupported dielectric modes fail explicitly.
These restrictions concern active `DF` jobs, not legacy no-DF decks.
In particular, dielectric/ground loss support does not add lossy conductors:
LD-based finite-conductivity or lumped conductor loading remains excluded.
Conductors retain the established NEC thin-wire/patch approximations.

## Defining a model

Keep conductors, frequencies, excitation and output requests in the NEC deck.
Add one `DF` card after `GE` to select an external dielectric resource:

```text
CM SMALL LOSSY DIELECTRIC NEXT TO A DIPOLE
CE
GW 1 11 0 0 .36 0 0 .84 .001
GE
FR 0 1 0 0 300 0
EX 0 1 6 0 1 0
DF "object.diel"
RP 0 3 2 0000 0 0 45 90 0 0
EN
```

`object.diel` describes physical boxes, not explicit numerical voxel records:

```text
NEC2C-DIELECTRIC 2
MATERIAL 1 ER 3 SIGMA .02
REGION 1 MATERIAL 1
BOX 1 .10 -.01 .58 .12 .01 .60
MESH MAX_KH .5 MAX_ASPECT 1 MAX_CELLS 100
END
```

Coordinates and cell dimensions are in meters; NEC `FR` is in MHz; conductivity
is S/m. Header version 1 is lossless, version 2 adds constant EI or conductivity,
and version 3 adds frequency tables. For exp(+i omega t),
`epsilon_relative = ER - i*(EI + sigma/(omega*epsilon0))`, with **either EI or
sigma** in the resource, not both. ER=1 with loss is active material; zero-loss
ER=1 contributes no volume unknowns. Table interpolation is linear in frequency
and the supplied loss parameter; extrapolation is rejected. There is no automatic
moisture/species model or effective-medium mixing law.

See [the resource contract](DIELECTRIC-RESOURCE-FORMAT.md) for complete grammar,
frequency-table examples, `NX` lifetime, fixed-mesh sweep behavior and diagnostics.
The NEC card limit is 132 characters; path limits are inherited from NEC2C.

## Choosing a solver and mesh

```sh
./nec2c -s dense -i model.nec -o model.out
./nec2c -s gmres -r 1e-10 -n 300 -i model.nec -o model.out
./nec2c -s fft -r 1e-10 -n 300 -i model.nec -o model.out
```

`dense` is the default reference LU solver. `gmres` uses an assembled matrix;
it is not a matrix-free large-volume path. `fft` accelerates the volume block
using equal cubic cells on a common Cartesian lattice; it retains conductor
coupling/storage costs. There is no silent fallback when an FFT lattice, budget
or convergence check fails. GMRES defaults are a 1e-12 relative target, restart
30 and 300 maximum iterations; `-r` and `-n` control the target and iteration cap.

Each physical box is subdivided independently under the material/frequency
mesh bound. `MAX_KH` bounds `k0*sqrt(abs(epsilon_relative))*longest_edge`;
it is not an exact pitch request. Values <=.5 are the normal starting range;
>.5 through 1 warn, and >1 is rejected. `MAX_ASPECT` controls cell shape, and
`MAX_CELLS` is a hard allocation guard. Frequency tables use a conservative
whole-table bound to keep a fixed mesh across a sweep.

For FFT, arrange box dimensions and boundaries so their actual subdivisions
produce equal cubic cells at common lattice sites. PEC image embedding also
requires bottom clearance on the half-cell lattice. Review the emitted divisions
and active-cell count; do not infer them solely from the physical box count.
The FFT path enforces an initial 512 MiB **accounted solver-storage** budget,
which is not a cap on total process RSS. Sparse objects can still require a
large padded rectangular FFT grid. Dense storage grows quadratically in total
unknowns, and the dense factorization grows cubically. Mesh preflight prints
its projected dense cost even when the selected solve is accelerated.

Use `DF` with `EN` and no excitation for a mesh/cost preflight. For conductor
clearance the required gap is at least the larger of cell longest edge and
segment length (or sqrt(patch area)); wire radius is included. Patch clearance
uses an equivalent-area proxy, not exact polygon geometry. Refining cells alone
may not permit a narrow gap if the conductor discretization remains coarse.

## Reading results

The output reports conductor currents, volume polarization-current density
(A/m^2), input impedances for driven feeds, numerical residuals, and dielectric
absorption. `RP` mode 0 reports complex range-normalized far electric fields
`r*E_theta`/`r*E_phi` in volts. Plane-wave results are scattered fields and include
RCS where supported. Coupled fields include conductor and dielectric radiation;
they are not automatically the traditional NEC gain table. `NE`/`NH` return
exterior scattered fields (V/m or A/m), not incident-plus-scattered plane-wave
fields. Their observation points must remain at least one longest cell edge
outside all dielectric cells. Internal cell-averaged fields are test diagnostics,
not a general internal/interface near-field evaluator. Ground patterns are
restricted to the upper hemisphere and GN 2 observations to z>0. Ground
plane-wave incidence uses theta in 0..90 degrees; elliptical ground excitation
requires an explicit positive axial ratio.

`RP`, `NE` and `NH` can each execute a solve. The free-space example deliberately
shows all three requests and therefore has three current/power blocks for one
frequency. Read the section labels and requests rather than treating every block
as a new frequency. For a long model, request only the outputs you need.

Driven GN 2 power separates accepted power, air radiation, lower-half-space
flux and dielectric absorption. The lower flux is eventually absorbed in lossy
infinite ground, but represents transmission for lossless ground. Zero ground
conductivity at ER>1 does not produce free space or PEC. Plane-wave balances
exclude the unperturbed incident wave's substrate absorption.

A small residual confirms the linear solve; a small power defect checks energy
consistency. Neither certifies physical material data or mesh accuracy. Refine
volume and conductor discretizations separately and compare the observables you
need. The validated component-wise comparison allowance
`max(2 ohm, 2% of the reference impedance-perturbation component)` is a benchmark
acceptance rule, not a general error bar on arbitrary models.

## Reference models and validation

[Canonical examples](../examples/dielectric-r1/README.md) cover driven free-space
coupling and exterior fields, touching/table-valued materials on PEC, a lossy
cube above GN 2 ground, and dielectric-only plane-wave scattering. Their small
size makes them suitable for first use and output-parser development.

Larger reference work is retained in the M9/M11 9,432-cell tree fixtures and
M13-05's 156/1,248-cell lossy-tree sensitivity experiment. The accepted historical
PEC-tree runtime target is 40 minutes per frequency on the constrained VM.
This is fixture-specific; a full-sized lossy or GN 2 tree has not been resource
certified. GN 2 has no FFT acceleration, and its current combined-loss physical
benchmarks use up to 64 cells.

Run `make check` for the 24 regression groups and `make distcheck` for a fresh
source-distribution build/test/install check. OpenEMS reproduction scripts and
numbered acceptance records live in `tools/` and `tests/`. M13-01/02/03/04 retain
independent returned material, decay, mesh/domain and impedance checks. M13-05
adds material provenance/sensitivity reporting and an explicit remote-ground
limit; it does not claim measured vegetation parameters.

## Documentation map

* [Deferred features and simplifying decisions](DIELECTRIC-DEFERRED-WORK.md)
* [Revision-1 release notes and announcement](DIELECTRIC-REVISION-1-RELEASE.md)
* [Material provenance and uncertainty guidance](M13-05-TREE-MATERIAL-GUIDANCE.md)
* [Design, equations and milestone history](DIELECTRIC-VOLUME-DESIGN.md)
* [Cuboid kernel and singular/contact quadrature](VIE-CUBOID-KERNEL.md)
* [FFT integration history](M9-FFT-INTEGRATION.md)
* [Ground formulation and scope](M10-GROUND-FORMULATION.md)

Milestone documents retain their original intermediate assumptions and timing
records. This revision guide and the completed M13 scope supersede earlier
statements that loss, heterogeneous FFT or PEC ground were still future work.
