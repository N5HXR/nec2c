# Dielectric revision 1: deferred work and simplifying decisions

Revision 1 is complete for the scope in [the user guide](DIELECTRIC-REVISION-1.md).
This register distinguishes unimplemented features, deliberate model choices,
and validation/resource limits. These are opportunities for later revisions;
none is silently included in the revision-1 acceptance claim. Items are not
ranked commitments or estimates.

## Deferred capabilities

| ID | Deferred item | Revision-1 boundary and future work |
| --- | --- | --- |
| D01 | **M12: conductor/dielectric contact, penetration and embedding** | Deliberately deferred by the user. Current wires and patches need positive conservative clearance. Needs contact/interface formulation, topology handling, singular integration and independent validation; includes insulation, embedded wires and contacting enclosures. |
| D02 | Anisotropic permittivity | Scalar materials only. Tensor material laws, orientation/grain metadata, coupling and analytical/independent tests are needed for directional materials. |
| D03 | Magnetic materials | mu_r=1. Magnetic constitutive response/unknowns and validation need a separate design. |
| D04 | GN 0 dielectric coupling | Fresnel/reflection-coefficient ground approximation is rejected. GN 2 wires-first was explicitly selected; GN 0 needs its own consistent excitation, interaction and power treatment. |
| D05 | Real-ground conducting patches | GN 2 coupling accepts wires, not NEC patches. Reconcile the patch ground approximation with the chosen volume half-space operator before extending it. |
| D06 | Real-ground contact, buried/cross-interface volumes and conductors | All GN 2 dielectric boxes, including vacuum boxes, must have zmin>0; wires stay above ground. Requires a cross-interface model and dedicated interface/contact quadrature. PEC dielectric contact is already supported and is a different boundary problem. |
| D07 | Radial screens, cliff and nonplanar/layered ground | Infinite planar single half-space only; screens/cliff are rejected. Additional ground geometry/layers need compatible kernels, fields and flux accounting. |
| D08 | Ground FFT / scalable GN 2 trees | GN 2 uses dense or assembled GMRES. A reflected-kernel accelerator and realistic resource/accuracy tests are future work. The free-space/PEC full-tree timing cannot be transferred to real ground. |
| D09 | Wire loading, networks and transmission lines with DF | Active dielectric jobs reject LD, NT and TL, including LD-based finite-conductivity conductor loss. Their operators and source/power bookkeeping need integration with the joint solve. |
| D10 | Additional sources and NEC output modes | Ordinary voltage-source EX 0 and plane-wave EX 1–3 are supported. Discontinuity sources and unsupported excitations/RP modes need implementation; arbitrary legacy output/gain semantics are not implied by complex-field output. |
| D11 | NEC symmetry reduction with DF | Expanded explicit geometry is used. Symmetry-aware volume/conductor coupling and lattice handling remain future work. |
| D12 | General internal/interface near fields | NE/NH use exterior separated observations. Singular internal/boundary E/H evaluation needs a distinct formulation. Cell-averaged internal test diagnostics are not that feature. |
| D13 | Native arbitrary dielectric CAD/CSG and conforming boundaries | Resource geometry is a union of axis-aligned nonoverlapping boxes. External producers perform geometry conversion; arbitrary curves/surfaces, cut cells and conforming mesh support need additional geometry/kernel work. |
| D14 | Alternative dielectric formulations | No PMCHWT/surface-only dielectric backend, magnetic equivalent-current model, or point-dipole replacement. These require separate conditioning, singular integration and validation decisions. |
| D15 | Wider constitutive laws | ER>=1 and passive loss only. Gain, negative/sub-unity permittivity, nonlinear/time-dependent materials, automatic causal broadband fits and extrapolation are not supported. |
| D16 | Empirical vegetation/material library | No built-in species, moisture, grain, temperature or effective-medium calibration. The tree material values in examples are declared assumptions; measured datasets and calibrated models remain separate work. |
| D17 | Broadband network/S-parameter workflow | No automatic multiport sweep/matrix extraction. External tools may derive quantities from supported repeated excitations; this does not add NT/TL support to DF jobs. |
| D18 | Product integrations and platform certification | The DF contract is available for Antenna Turtle and other producers. A complete GUI/CSG export UX is external work. Revision-1 native Windows/MSYS2/macOS certification and cross-platform numerical evidence were not performed here. |

| D19 | Parallel, GPU and distributed dielectric backends | The implemented CPU paths were validated on the constrained VM. No dedicated GPU/distributed dielectric backend is provided; parallel kernel builds/matvecs and broader hardware benchmarking are future optimization work. |

## Simplifying choices and improvement opportunities

| ID | Current decision | Implication / possible improvement |
| --- | --- | --- |
| S01 | Three piecewise-constant Cartesian volume-current components per cuboid with Galerkin testing | Resolves volumetric polarization rather than replacing a whole dielectric object with its centroid. Thin/curved structure needs enough cells; higher-order bases, adaptive or conforming discretizations could improve cost/accuracy. |
| S02 | Per-box uniform subdivision | MAX_KH and aspect limits control electrical size and cell shape, not a direct user pitch. Different boxes may yield incompatible FFT pitches. A shared-grid producer, explicit pitch contract or adaptive mesh could improve usability. |
| S03 | Equal cubic common-lattice cells for FFT | Heterogeneous scalar values are supported, but arbitrary nonuniform meshes are not. A padded rectangular embedding can waste memory for sparse trees; sparse FFT, hierarchical methods or other accelerators remain alternatives. |
| S04 | Regular-grid FFT selected over an initial octree/H-matrix/tree accelerator | Chosen using the small-tree experiments on a 2 GB VM. This is an algorithm choice for the accepted case, not proof that it is optimal at every size or occupancy. |
| S05 | Initial 512 MiB accounted FFT solver-storage budget | Protects the chosen VM workload, not total RSS. Other allocations and conductor coupling remain. Configurable budgets, more complete memory accounting and scalable conductor blocks could extend capacity. |
| S06 | Dense LU reference and assembled-matrix GMRES | GMRES alone does not remove matrix assembly/storage. Matrix-free general cuboid/ground operators and improved preconditioning could broaden scale. GMRES restart is fixed at 30; target/iteration cap are configurable. |
| S07 | Conservative conductor gaps | Required gap uses max(cell longest edge, wire segment length or sqrt(patch area)); wire radius is included. May reject physically separated geometries that could be accurate with better quadrature. Do not loosen without evidence; M12 is a separate contact feature. |
| S08 | Equivalent-area patch-clearance proxy | NEC stores centroid, area and tangents rather than original polygon extents. Exact polygon-aware clearance/intersection requires preserving/reconstructing more geometry. Elongated patches are a known proxy limitation. |
| S09 | Strict disjoint interiors and exact shared boundaries | No coordinate snapping or overlap tolerance. Producers must emit identical coordinates to touch; tiny gaps are not certified as contact and tiny overlaps are rejected. Robust CSG/tolerance handling is future producer work. |
| S10 | Axis-aligned material interfaces | Touching scalar-material faces/edges/corners are supported. Curved and very thin interfaces are represented geometrically by boxes and require refinement; there is no universal small-gap accuracy guarantee. |
| S11 | Linear frequency interpolation with no extrapolation | Straightforward, explicit table contract. Not automatically a globally causal broadband constitutive fit. Higher-quality passive/causal fits or user-selected interpolation need a separate validated contract. |
| S12 | Whole-table mesh bound | Fixed geometry/unknowns across a sweep, potentially overmeshing a narrow requested band. Tighter range-specific certified bounds or remeshing could save work. |
| S13 | Nonmagnetic isotropic scalar material assumptions | Grain direction and unresolved mixtures are not inferred. Per-region metadata, tensor materials and measured-data provenance could improve realism. |
| S14 | Homogeneous box crown/trunk idealizations in tree witnesses | No automatic leaves/bark/roots or subvoxel mixing law. Sampled sensitivity envelopes are not measured uncertainty or confidence intervals; empirical geometry/material validation is separate. |
| S15 | Explicit fixed-source normalization and independent power quadrature | Residual and energy checks complement observable refinement; neither certifies physical material data or every local field. Adaptive angular/spectral quadrature and stronger a-posteriori estimates could help. |
| S16 | Bounded Sommerfeld quadrature and a 16 MiB optional reflected-wire-basis cache | Difficult grazing/high-contrast cases may explicitly fail to converge. Better asymptotics/quadrature/caching could improve speed and coverage without masking failure. |
| S17 | Exterior near fields, upper-hemisphere ground patterns and complex rE output | Useful validated outputs with explicit conventions; arbitrary internal samples and all traditional NEC gain/table modes require additional work. Repeated RP/NE/NH requests may repeat a solve; solution reuse is an optimization opportunity. |
| S18 | Benchmark component allowance of max(2 ohm, 2% of reference component) | An agreed acceptance test on induced impedance changes, not a universal absolute-impedance, field or material-error guarantee. Additional observables and tighter reference studies can extend validation. |
| S19 | Fixture-specific runtime budgets | Historical free-space and PEC tree targets were accepted separately; PEC target is 40 min/frequency. Optimization was deferred by the user. Full-size lossy trees, GN 2 trees, multi-frequency performance and larger conductor systems need their own measurements. |
| S20 | Text DF/NEC interface and existing path/card limits | Portable external resources and no extra runtime dependencies, but no binary geometry cache, machine-readable full solver API or automatic model repair. Tooling can improve reproducibility and diagnostics. |
| S21 | Existing NEC conductor representation | Thin-wire and ideal conducting-patch approximations are retained. Dielectric/ground loss does not introduce conductor loss because LD is excluded. Conductor segmentation, close-wire gaps and electrically large patch concerns remain subject to NEC modeling limits; a volumetric conductor or improved conductor basis is separate work. |

## Evidence boundaries

The 24 regression groups, analytical limits and independently refined OpenEMS
M13-01/02/03/04 fixtures establish the accepted scope. A 9,432-cell lossless PEC
tree and a small lossy-tree sensitivity study are separate checkpoints, not a
measured real-tree model or a full-size lossy GN 2 certification. The OpenEMS
finite feed gap, explicit wire radius, mesh and finite PML domain differ from
the NEC idealizations; their mesh/domain/source choices are retained in the
numbered audit records. Acceptance does not mean those modeling differences
vanish.

The production-tree report deliberately keeps numerical refinement, physical
material uncertainty and scenario sampling separate. Extend this register and
add independent evidence whenever a deferred item is brought into scope.
