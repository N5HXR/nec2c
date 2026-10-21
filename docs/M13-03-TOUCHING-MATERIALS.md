# M13-03: touching scalar dielectric materials

Status: COMPLETE (2026-10-21). Native analytical, solver-equivalence and
refinement checks, all 22 regression groups and independent OpenEMS acceptance
pass. M13 as a whole remains in progress for real-ground work.

## Scope and formulation

Support multiple passive isotropic dielectric materials with disjoint interiors
and shared box faces, partial faces, edges or corners. Constants and version-3
tables keep their material identity. No new resource syntax is needed: versions
1/2/3 use the existing MATERIAL/REGION/BOX records. Overlapping interiors are
still rejected strictly, including tiny overlaps. Coordinates are not snapped:
use identical boundary coordinates to specify a touching interface. A zero-volume
shared boundary belongs to neither material's volume measure and is not counted
twice. Vacuum boxes remain omitted from the unknowns.

The unknown remains cell-average volume polarization current. For testing cell i,
with chi_i = epsilon_i - 1, the volume block is

    (A J)_i = epsilon_i J_i - chi_i sum_j N_ij J_j
    rhs_i = i omega epsilon0 chi_i E_inc,i.

The existing distributional cuboid kernel includes the self term and interface
charge from the discontinuous cellwise polarization. There is no additional
surface-current approximation and no new interface-current unknown. Do not impose
continuity of polarization current: its normal component can jump across a
material interface. Maxwell interface conditions are recovered under refinement.
The independent small-system reference integrates the electrostatic surface
kernel in closed form, including the finite incident-wave cell average.

Dense and assembled GMRES paths already use chi_i separately. The mesher's old
unlike-material contact rejection is removed. FFT now stores scalar complex
permittivity for every cell; the geometry-only background N kernel is still
translation invariant. The FFT computes N J, followed by local epsilon/chi
multiplication. Its self-block preconditioner uses each cell's own material.
PEC image cells copy their physical source material; only physical currents
are unknowns. The 512 MiB cap includes the added material storage.

FFT still requires equal cubes on a common lattice and the existing PEC image
embedding rules. The per-box mesher retains material-dependent subdivisions;
split boxes into a common cubic grid explicitly when selecting FFT for unlike
materials. Dense/GMRES also support unequal cuboids and nonconforming contacts.
Electrical cell-size limits alone do not certify high-contrast or thin-layer
interface accuracy; verify volume refinement for the intended geometry.
This step does not extend real-ground support, conductor/dielectric contact,
penetration, anisotropy or magnetic materials.

## Native verification

`tests/m13_materials.py` checks independent closed-form electrostatic two-cell
face interactions for normal and tangential polarization; local Joule absorption;
face/edge/corner/partial-face topology; strict overlap refusal; equal-material
limits; material relabeling; touching EI/SIGMA tables, including conductive ER=1;
free-space/PEC dense/GMRES/FFT agreement; and coupled-wire currents/power. Patch and mixed-conductor checks retain the
existing free-space/PEC conductor paths.
C FFT tests compare heterogeneous lossy free-space and PEC products to direct
cuboid integration, solve them with GMRES and check invalid material updates.

A second independent reference uses the exact electrostatic coated-sphere
polarizability. For core epsilon c, shell epsilon s, outer radius b and f=(a/b)^3:

    alpha/epsilon0 = 4 pi b^3 * [(s-1)(c+2s)+f(c-s)(1+2s)]
                                  / [(s+2)(c+2s)+2f(s-1)(c-s)].

This follows by matching potential and normal epsilon*E at both spherical
boundaries for the core, shell and exterior dipole potentials.
At exp(+i omega t), absorbed power in the electrostatic limit is
`-omega epsilon0 Im(alpha/epsilon0)/2` for unit incident field. Test core
ER=4, EI=.4 and shell ER=2, EI=.2 at 1 MHz; radii a=18 mm, b=30 mm.
Voxel-center inclusion approximates both curved boundaries using explicit
nonoverlapping cubes. Diameter grids 8/12/16 give 280/912/2176 active cells.
The preselected gates are <8% polarization/absorption error on each grid and
<3% on the finest, with finest polarization closer than coarse. Staircase
geometry can give nonmonotone absorption error; no smooth convergence order
is claimed. Reproduction:

    python3 tests/m13_materials.py
    python3 tools/m13_03_acceptance.py --workdir /tmp/m13-03-native

The antenna benchmark explicitly partitions the pillar into 25 mm cubes to
preserve one lattice across materials; use `--pitches .0125 --frequencies 320`
for volume refinement, or `--segments 21 --frequencies 320` for wire refinement.
Every solve must meet residual <=1e-10, nonnegative absorption and <0.3% global
power defect. The native benchmark caches exact binary/deck/material identities.

## OpenEMS experiment M13-03

Self-contained package: `/tmp/m13-03-openems.tar.gz`. On the Mac with the existing
OpenEMS Python environment:

    tar -xzf m13-03-openems.tar.gz
    cd m13-03-openems
    bash run-m13-03-materials.sh

The included README is sufficient; no commands from a prior chat are needed.
Return `m13-03-results.tar.gz` to `/tmp` on rxtx. The driver source is
`tools/run-m13-03-materials.sh`; incomplete results are packed on exit and exact
completed runs are cached on restart.

The PEC pillar/wire/feed/source/domain fixture is the accepted M13-01/02 setup.
Pillar x=50..150, y=-50..50, z=0..650 mm; wire radius 5 mm; conductor clearance
45 mm; finite OpenEMS feed gap 10 mm. Split-x uses x=100 mm and split-z z=325 mm.
Near/lower material ER=2, sigma=.01 S/m; far/upper ER=4, sigma=.03 S/m. These
orientations exercise predominantly tangential and normal polarization across
the interface. No dielectric/conductor contact is introduced.

For each orientation: one -40 dB pilot, then -60 dB vacuum/loaded pairs at
20 mm bulk/5 mm wire, 15/2.5 mm, and 15/2.5 mm with an extra 150 mm side/top
margin. The interface is an explicit mesh anchor even in vacuum. Twelve
verification runs plus two pilots run sequentially, one CPU thread. Largest
planned mesh: 1,401,744 cells. Shared pair wall caps are 900/1200/1500 seconds;
pilot cap 180 seconds. These are ceilings, not measured runtime predictions.

Both material laws are constant: all 280/300/320 MHz rows are valid, unlike
M13-02's frequency-designated witnesses. Cache identities and offline reviewer
grouping include material_split so orientations cannot be cross-paired.

Acceptance requires explicit final verification decay <=-60 dB before step cap,
no clipped source, XML material/interface audit, mesh/domain convergence, and
native induced delta R/delta X errors <=max(2 ohm, 2% of each OpenEMS reference
component). Complex-vector error is diagnostic only. Returned evidence gates
acceptance of this step; development continues without waiting for it.

## Conductor power fixture refinement

The initial new mixed wire/patch check used an 11-segment wire and exceeded
the unchanged 0.3% power-defect gate. A homogeneous-material control gives
0.37402% in free space and 0.39041% over PEC, essentially identical to the
heterogeneous result: this is conductor discretization error, not an interface
effect. Using the same 21-segment central-feed wire as M13-01 reduces the
heterogeneous defect to 0.16080% and 0.17741%. Patch-only cases already pass.
The corrected fixture passes dense/FFT agreement and the original power gate;
all 22 regression groups pass with an unchanged solver binary.

## Independent acceptance closure (2026-10-21)

All 12 OpenEMS verification runs explicitly reach -60 dB before the step cap;
final decay ranges from -64.88 to -60.15 dB. No source is clipped. Both scalar
materials and their exact interface coordinates match the intended XML geometry.
All 14 meshes, including pilots, have adjacent step ratios <=2 within rounding,
minimum steps >=2.5 mm, and fewer than two million cells. Maximum OpenEMS mesh
refinement change is 0.40850 ohm; domain enlargement changes induced impedance
by at most 0.04449 ohm.

All 11 native benchmark/refinement pairs meet residual <=1e-10, nonnegative
absorption and <0.3% power defect, and pass the agreed component-wise OpenEMS
impedance gate. All 280/300/320 MHz rows compare constant material laws.

The largest component discrepancy is 1.97370 ohm: stacked materials at 320 MHz,
12.5 mm volume pitch and 11 wire segments. This is inside the 2-ohm floor, with
a narrow margin. The joint 12.5 mm/21-segment refinement gives induced
Z=1.78358-j3.54676 ohm, against OpenEMS 1.60803-j5.24388 ohm: resistance error
0.17556 ohm and reactance error 1.69712 ohm. Its power defect is 0.011884%,
true residual 9.6934e-11, solver CPU 70.22 s, wall time 907.72 s, peak child RSS
33,308 KiB and accounted solver storage 31,252,640 bytes. The 21-segment wire
provides more comparison margin for this fixture; tolerances are unchanged.
Complex-vector discrepancy remains diagnostic only.

At 320 MHz, halving volume pitch changes induced impedance by 0.10467 ohm
(side-by-side) and 0.12080 ohm (stacked). Refining the wire from 11 to 21
segments at 25 mm pitch changes it by 0.14686 and 0.30438 ohm, respectively.
The analytical coated-sphere finest-grid errors are 2.1423% in polarizability
and 2.5457% in absorption. These results establish the documented scalar
material-interface scope, without extending conductor-contact or real-ground
acceptance.

Retained native evidence, refinement measurements and regression provenance:
[m13-03-acceptance.json](../tests/m13-03-acceptance.json). Returned archive hash,
DFT-reproduced impedance pairs, explicit decay, material/interface XML and mesh
audits are in [m13-03-openems-review.json](../tests/m13-03-openems-review.json)
and [m13-03-openems-pairs.csv](../tests/m13-03-openems-pairs.csv).
Source distribution and whitespace checks pass. No additional OpenEMS runs
are required to close M13-03; no later increment is started by this closure.
