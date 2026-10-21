# M13-04: scalar dielectric loss over GN 2 ground

Status: COMPLETE (2026-10-21). M13-01/02/03 remain complete. This increment
validates the retained GN 2 implementation with the scalar material extensions;
it does not yet close M13 or accept new real-ground conductor/contact modes.

## Scope

Nonmagnetic passive infinite planar ground, NEC `GN 2` Sommerfeld/Norton wires,
and positive-clearance scalar dielectric volumes. Support driven wires and
volume-only or wire/volume plane-wave excitation, constants and version-3 tables,
and disjoint or touching dielectric materials. Dielectric materials may touch
each other; all boxes, including vacuum boxes, must retain strictly positive
ground clearance. Conductors retain the existing dielectric clearance policy.

Retain the explicitly selected GN 2/wires-first boundary: GN 0, real-ground
patches, screens, cliff ground, contact/burial and ground FFT remain rejected.
Use dense or assembled GMRES. PEC-ground contact remains the separate accepted
M11/M13-03 scope. Anisotropy and magnetic media remain excluded.

The existing derivation is in
[M10-GROUND-FORMULATION.md](M10-GROUND-FORMULATION.md), especially the GN 2 wire
integration and independent interface-flux formulas. M13-04 supplies current
combined-loss acceptance evidence without replacing that history.

## Physics and analytical gates

For exp(+i omega t), volume epsilon is ER-i*(EI+sigma/(omega*epsilon0)). GN's
negative conductivity convention specifies a fixed imaginary relative
permittivity; positive GN conductivity is evaluated at each frequency. Every
solve/sweep must rebuild the appropriate ground state without mutating its law.

Driven power must satisfy

    P_accepted = P_air + P_lower_halfspace + P_absorbed_dielectric.

Plane-wave scattered-power checks use extinction, scattered air radiation and
scattered lower-half-space flux; the unperturbed incident-wave substrate
absorption is excluded. Lower-half-space flux is independently integrated with
wire/volume interference retained. For lossy infinite ground it is eventually
absorbed; for lossless ground it is transmitted. Zero ground conductivity with
ER>1 is therefore not a free-space or PEC limit.

The new `tests/m13_real_ground.py` covers:

* an independent quasistatic finite-cube image reference, using closed-form
  Coulomb face integrals with electrostatic coefficient (epsilon_g-1)/(epsilon_g+1),
  together with exact Fresnel plane-wave illumination and cell averaging;
* local Joule absorption, matched-ground recovery of free-space lossy currents,
  lossless-ground transmission, and dielectric EI/sigma/zero-loss limits;
* touching lossy materials, table interpolation, conductive/fixed-imaginary
  ground sweep state, dense/GMRES agreement and independent power channels;
* driven-port reciprocity with refined wires; and
* deterministic refusals of GN 0, real-ground patches, contact and FFT.

Existing half-space tests already cover Fresnel/Brewster limits, PEC-image,
matched and quasistatic limits, weighted kernel reciprocity, independent axial
Sommerfeld integrals, E/H curl and dipole work versus air/interface flux.
They remain part of the regression evidence. The new static finite-cube test
uses a 10 mm cube centered 10 mm above ER=4-i*.4 ground, ER=3-i*.4 material,
1 MHz normal incidence, and a 5e-6 relative-current gate. It checks both the
self-consistent polarization and its Joule dissipation.

## Native antenna experiments

`tools/m13_04_acceptance.py` runs cached, binary-hash-pinned counterparts of the
OpenEMS wire/cube fixture. Each run independently records accepted power,
upper radiation, lower flux, dielectric absorption, flux quadrature change,
true linear residual, CPU/wall time and peak child RSS. Gates: residual <=1e-10,
nonnegative absorption/radiation, positive lower flux, and power defect <0.3%.

The initial current-binary set is:

* lossy ER=3/sigma=.02 cube at 280/300/320 MHz with 3x3x3 cells, 11 wire segments;
* lossless ER=3 control at 300 MHz on the same mesh;
* lossy 300 MHz volume refinement to 4x4x4 cells; and
* lossy 300 MHz wire refinement to 21 segments with 3x3x3 cells.

The existing 64-cell lossless 300 MHz checkpoint is retained in
`tests/m10-nec-real-300.json`; it remains provenance, not a substitute for new
lossy verification. These experiments deliberately bound dense real-ground
work to <=64 cells on the 2 GB VM. The half-space operator is much more costly
than the free-space/PEC FFT path; time the actual runs before selecting larger
production fixtures. No unsupported ground FFT acceleration is implied.

Reproduce with:

    python3 tests/m13_real_ground.py
    python3 tools/m13_04_acceptance.py --workdir /tmp/m13-04-native
    python3 tools/m13_04_acceptance.py --workdir /tmp/m13-04-native-lossless --frequencies 300 --sigma 0
    python3 tools/m13_04_acceptance.py --workdir /tmp/m13-04-native-refined --frequencies 300 --divisions 4
    python3 tools/m13_04_acceptance.py --workdir /tmp/m13-04-native-wire --frequencies 300 --segments 21

## OpenEMS experiment M13-04

The self-contained `/tmp/m13-04-openems.tar.gz` contains its own README and
`run-m13-04-real-ground.sh`. On the Mac with the existing OpenEMS Python:

    tar -xzf m13-04-openems.tar.gz
    cd m13-04-openems
    bash run-m13-04-real-ground.sh

Return `m13-04-results.tar.gz` to `/tmp` on rxtx. The source driver is
`tools/run-m13-04-real-ground.sh`. It packs partial results on exit and caches
exact normally completed configurations on restart.

The cube is x=50..150, y=-50..50, z=550..650 mm; 5 mm wire radius, 475 mm wire
span centered at 600 mm, 10 mm finite OpenEMS feed gap. Soil ER=13/sigma=.02 S/m
occupies z<0. The cube has 550 mm ground clearance and 45 mm conductor clearance.
Control sets use cube ER=3/sigma=0 and ER=3/sigma=.02. Both laws are constant:
all 280/300/320 MHz rows are valid. No material is buried or touches ground.

Each control has one -40 dB pilot and -60 dB vacuum/loaded pairs on 20/5 mm,
15/2.5 mm, and 15/2.5 mm with 150 mm added side/top/soil-depth margin. The largest
planned grid has 1,966,500 cells, minimum step >=2.5 mm and <=2:1 grading.
Twelve verification runs plus two pilots run sequentially. The pair wall caps
are 1200/1800/2400 seconds and pilot cap 240 seconds; these are ceilings, not
measured runtime predictions. The larger soil domain and its transient decay
can add work; the actual engine time step must be read from the returned logs.

Review the controls separately: their vacuum jobs have identical physics but
are intentionally stored in distinct control directories. The wrapper prevents
cross-pairing or duplicate-job confusion:

    python3 tools/m13_04_review.py /tmp/m13-04-openems-review/raw --output-dir /tmp/m13-04-openems-review/audit

That wrapper reproduces native-timestamp port DFTs but does not grant physical
acceptance. Final gates require explicit verification decay <=-60 dB before
step cap, no clipped source, soil/cube material and geometry XML audits,
independent mesh/domain convergence, and native induced delta-R/delta-X errors
<=max(2 ohm, 2% of each corresponding OpenEMS reference component). Complex-vector
error is diagnostic only. Returned results gate step acceptance; development
continues while those runs are external.

Retained status and native evidence will be recorded in
[m13-04-acceptance.json](../tests/m13-04-acceptance.json).

## Returned acceptance evidence (2026-10-21)

All 23 regression groups, six native antenna pairs and all 12 returned OpenEMS verification runs pass
the physical acceptance gates. The 18 OpenEMS spectral pairs cover both material
controls at 280/300/320 MHz across three mesh/domain configurations. Each verify
run explicitly reached <=-60 dB before its step cap, with no clipped source.
The XML audits confirm isotropic ER=13/sigma=.02 soil, the selected cube loss,
5 mm wire radius, 10 mm source gap and all fixture bounds. Fourteen grids,
including pilots, have minimum step 2.5 mm, <=2:1 grading and at most
1,966,500 cells.

Largest native/OpenEMS component error is 0.829761 ohm. Largest OpenEMS mesh
change is 0.280831 ohm; largest domain change is 0.029949 ohm. At 300 MHz,
27-to-64-cell native volume refinement changes delta-R/delta-X by
0.000136/0.004346 ohm; 11-to-21 wire segments changes them by
0.050306/0.004057 ohm. Maximum native power defect is 0.264241%, below the
0.3% gate. The independent quasistatic current error is 5.81e-9 relative.
These are fixture-specific acceptance results, not a universal accuracy bound.

Evidence is retained in `tests/m13-04-openems-review.json`,
`tests/m13-04-openems-pairs.csv` and `tests/m13-04-acceptance.json`.
Reproduce the physical audit after extracting the returned archive:

    python3 tools/m13_04_audit.py /tmp/m13-04-openems-review/raw --archive /tmp/m13-04-results.tar.gz --native /tmp/m13-04-native /tmp/m13-04-native-lossless /tmp/m13-04-native-refined /tmp/m13-04-native-wire --output tests/m13-04-openems-review.json

The audit fails on incomplete decay, wrong materials/geometry, inadequate mesh
quality or failed component comparisons. It is separate from the regression
suite and requires the external archive and completed native benchmark sets.
No production solver changes were needed for this increment. GN 0, real-ground
patches, contact/burial, screens/cliff and ground FFT retain their exclusions.
