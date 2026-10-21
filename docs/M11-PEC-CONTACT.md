# M11: lossless dielectric contact with perfect ground

Status: COMPLETE (2026-10-21). M10 PEC positive-gap validation is complete.
M11 permits dielectric/PEC-plane contact only; conductor/dielectric contact
remains deferred M12, and real-ground contact is not implied (M13).

## Formulation and first increment

For a constant Cartesian polarization current on a cuboid in z>=0, the PEC
background remains the free-space direct operator plus the mirrored source
with current transformation Q=diag(-1,-1,+1). When the bottom face is at z=0,
the cell and its image share a face but their interiors do not overlap.
Therefore the direct block retains its singular self treatment and the image
block uses the existing face-touching cuboid quadrature. Contact does not
justify a point-dipole substitute or an arbitrary positive gap. Source current
need not have vanishing normal component at PEC; the tangential total electric
field is the PEC boundary condition.

The low-level `vie_cuboid_pec_n` interface now permits zero clearance when
`minimum_gap=0`. A requested positive minimum still excludes contact, and
buried or straddling cells remain rejected. No fuzzy tolerance silently moves
buried volume onto the ground. The initial increment retained the deck-level gate. The subsequent integration
now permits PEC dielectric contact in dense, GMRES, and the supported FFT
subset; buried cells, real-ground contact, and conductor/dielectric contact
remain rejected.

The extended kernel tests compare the reflected contact block against the
independently generated high-order M4 face-contact oracle, rotated so its
contact normal is z. They also check direct-self/reflected-face classification,
continuous convergence from positive gaps, unequal-cell weighted reciprocity,
and continued minimum-gap/buried-geometry rejection. The extended ground
kernel test executable passes all these checks. These are the first
kernel checks, not complete M11 acceptance.

## Development and acceptance record

1. Validate static/dynamic self, face, edge and corner image interactions,
   including independent contact-oracle convergence beyond one canonical pair.
2. Enable validated PEC-only deck contact with mesher-aware handling of a
   bottom face at zero; keep GN 2 positive-gap and conductor clearance rules.
3. Check solved currents, feed impedance and power as a positive gap closes;
   verify PEC boundaries and image equivalence for contact models.
4. Compare canonical contact cases against independently refined OpenEMS
   results. Stage runnable scripts for the laptop when the model is ready.
5. Complete grounded-tree accuracy/resource validation. A reflected-source
   FFT embedding is now implemented and passes small dense-reference tests.
   The 216-cell grounded skeleton completes in 89.59 s wall time with true
   residual 5.43e-11 and 28,064 KiB peak RSS. The full 9,432-cell tree now converges in 45 iterations at 9.07e-11 true
   residual, with 240,416 KiB peak RSS and 244,753,200 accounted solver bytes.
   Its observed 34.6-minute elapsed time passes the 40-minute grounded-tree
   target accepted by the user on 2026-10-21. Further optimization is deferred.

The grounded-tree production target is 40 minutes per frequency, with the
512 MiB accounted solver cap and accuracy requirements unchanged. The original
M9 free-space target remains 30 minutes. Verification runs may run longer.

## Initial near-gap failure (resolved by the mapping below)

For a 0.3 x 0.4 x 0.2 m cuboid at k=2.1 rad/m, exact PEC contact passes
the independent face oracle. An exploratory approach sequence at bottom gaps
1e-3, 1e-4 and 1e-5 m produced explicit quadrature failure at the latter two
gaps under the existing defaults (relative tolerance 2e-8, absolute tolerance
2e-11, maximum order 64). No returned matrix from a failed call is accepted.
The retained regression approach sequence covers the resolved 1e-2, 3e-3,
1e-3 m range; it does not certify the unresolved smaller gaps.

This is a near-singular quadrature issue, not evidence that the contact
operator diverges. The production reduction uses Duffy sectors when all
varying difference-coordinate intervals have an endpoint exactly at zero;
a small nonzero interval endpoint can instead trigger tensor quadrature.
A candidate next step is a closest-corner Duffy mapping or local subdivision
that resolves the small offset without changing the geometry or silently
relaxing tolerances. It needs independent oracle and existing-kernel regression
coverage before use. Tiny-gap convergence remains an explicit M11 gate.


## Near-gap mapping and independent static oracle

The reduction now uses Duffy sectors about the closest interval corner when
its distance from the origin is small compared with the interval extent.
This is a partition of the actual integration box: offsets are preserved,
not snapped to zero. The existing radial power map uses the closest-corner
separation, including varying interval offsets, to resolve small nonzero gaps.
The original tolerances and maximum Gauss order remain unchanged.

Dynamic PEC kernel tests now converge at bottom gaps 10 mm, 3 mm, 1 mm,
0.1 mm, 0.01 mm and 1 micrometre. Integration tests of solved currents also
converge through 1 micrometre and preserve the exact-contact image solution
for linear and both elliptical incident polarizations. Wire, closed-patch
and mixed contact models pass dense/FFT current and lossless power checks.
The driven one-cell feed approaches contact with errors 0.000350, 0.0000404,
and 0.00000412 ohm at gaps 1 mm, 0.1 mm and 0.01 mm.

`tools/m11_static_reference.py` independently integrates parallel rectangular
Coulomb face potentials using analytic primitives of 1/r, u/r, v/r and uv/r.
It has no numerical quadrature or production imports. Its retained eight
static contact/near-gap tensors are in `tests/m11-static-reference.csv`.
The extended C ground suite checks these against the production kernel.

User guidance (2026-10-20): exact contact and separated geometry may have
different numerical support limits. A documented minimum reliable relative
gap or explicit accuracy warning is acceptable if extreme gaps remain
problematic. No universal gap bound is established by the tests above;
quadrature failure still rejects a solve rather than silently accepting an
inaccurate matrix. Do not infer support for arbitrary sub-cell gaps solely
from exact-contact support.

## PEC FFT representation

Embed real cells above their reflected cells on one doubled-height lattice.
Before each convolution, fill image sources with Q J, Q=diag(-1,-1,+1).
Read the operator result only on the real cells; image currents are not
additional unknowns. The free-space self-block preconditioner is retained
as an approximation; true residual tests certify the physical system.
Conductor and conductor/volume blocks keep their existing PEC formulation.

This path requires the existing equal-cube, single-real-permittivity lattice
subset, with bottom gap an integer multiple of half the cube pitch. Contact
is the zero-gap member of that subset. Real ground FFT remains unsupported.
Physical nonnegative cell bounds are checked before accommodating the usual
floating-point cube/lattice roundoff. Temporary image vectors and doubled
lattice maps are counted in the unchanged 512 MiB joint solver budget.

## Laptop and VM pillar benchmarks

The M10 wire is retained (475 mm outer span, 5 mm radius, center z=600 mm),
while the epsilon-3 dielectric cross-section x=50..150 mm, y=-50..50 mm now
extends from z=0 to 650 mm. The physical wire/dielectric gap remains 45 mm.
OpenEMS compares 20 mm/5 mm-spacing and 15 mm/2.5 mm-spacing grids at
280/300/320 MHz. The portable runner's `--cube-bottom-mm 0` selects this
PEC-only fixture and adds that physics to job identities and offline pairing.
The loaded pilot is selected with `--pilot-epsilon-r 3`.

The staged archive `/tmp/m11-openems-bundle.tar.gz` includes
`run-m11-pec-contact.sh`, a self-contained README and checksums. The driver
is also retained as `tools/run-m11-pec-contact.sh` in the source distribution. One script
runs the loaded pilot and both mesh pairs, packages results, and prints the
return command. No command transcription from session history is required.

`tools/m11_pec_acceptance.py` stages cached native FFT counterparts with binary
hash and exact input provenance. Default 300 MHz refinements use 25 mm cubes
(416 cells) and 12.5 mm cubes (3,328 cells). A 50 mm trial was rejected by the
unchanged 45 mm wire/dielectric clearance versus 50 mm required cell clearance;
it is not an accepted refinement checkpoint. Larger conductor-clearance
requirements are not relaxed by dielectric/PEC-plane contact support.


## Returned physical and resource evidence (2026-10-20)

All four OpenEMS verification jobs reached -60 dB decay without source clipping
or the maximum-step ceiling. The coarse/fine mesh differences in complex delta Z
are 0.237/0.319/0.596 ohm at 280/300/320 MHz. Retained independent port-trace
DFT and decay evidence are in `tests/m11-openems-contact-review.json` and its
paired CSV.

The agreed physical gate is **component-wise**: delta R and delta X each differ
by no more than max(2 ohm, 2% of the corresponding OpenEMS delta component).
The 25 mm NEC comparison passes at all three frequencies. At 320 MHz the
resistance/reactance errors are 1.056/1.734 ohm, while the combined vector
magnitude is 2.031 ohm. The latter is a retained diagnostic, not a replacement
gate; the initial progress statement calling it a failure was corrected.

At 300 MHz, refining volume cells from 416 (25 mm) to 3,328 (12.5 mm) changes
delta Z by 0.121 ohm. Refining the wire from 11 to 21 segments at 25 mm volume
pitch changes delta Z by 0.189 ohm. The 21-segment 320 MHz comparison has
0.965/1.404 ohm resistance/reactance discrepancies. Its driven power defect is
0.0160%; the 300 MHz defect is 0.0224% versus the 11-segment 0.225% inherited
wire-basis defect. Refinement is not claimed to monotonically approach a
finite-feed OpenEMS model from the different NEC delta-gap source model.

All 19 aggregate regression groups pass. The expanded contact integration group
also passes explicit PEC-plane tangential-E/normal-H checks and physical port
reciprocity (relative discrepancy 7.88e-5, consistent with the native NEC
reference). ASan/UBSan checks of the PEC FFT adapter pass with LeakSanitizer
disabled because the execution environment traces processes. The source
archive builds successfully.

The complete grounded tree has 9,432 cells and 28,307 joint unknowns. It uses
the same 28 MHz M9 geometry with PEC ground and upper-hemisphere patterns.
The true residual is 9.07e-11 in 45 iterations. Accepted power is 6.200078 mW
and upper-hemisphere radiation 6.182667 mW; relative power defect is 0.2808%,
comparable to M9's inherited wire discretization defect. The 974/3896-ray
angular quadrature change is 0.0585%. Peak RSS is 234.8 MiB, within the 512 MiB
solver budget and the 2 GB VM's resources. Metadata and sampled complex fields
are retained in `tests/m11-pec-tree.json`.

Observed runtime is 2075.89 s (34.60 minutes), with 155.42 s of process CPU time.
This passes the 40-minute grounded-tree target explicitly accepted by the user
on 2026-10-21. This verification overlapped other checkpoint/test jobs; its
timing is an observed constrained-VM result, not an isolated hardware throughput
guarantee. Further optimization is deferred to future revisions after experience
with the models. **All M11 exit criteria are satisfied; M11 is complete.**

Wires and patches retain the M10 conservative positive clearance from ground,
as well as the existing positive conductor/dielectric clearance. This milestone
enables dielectric/PEC-plane contact; it does not extend the supported wire or
patch contact topology.
