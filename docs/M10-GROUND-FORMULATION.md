# M10 ground formulation and implementation record

Status: M10 PEC-only validation is COMPLETE (2026-10-20); see the acceptance
and returned sensitivity records at the end of this document. Scope was revised
at the user's request on 2026-10-20: real ground and its loss/verification work now
belong to M13. The implemented GN 2 path and evidence are retained as preliminary
M13 work; they are not M10 acceptance requirements. See the
[M13 handoff](M13-LOSS-AND-REAL-GROUND.md).

M12's conductor/dielectric contact, penetration and embedding are separately
deferred to a later revision. The current revision retains positive clearance
between conductors and dielectric regions; M13 does not depend on M12.

## Current M10 scope

Support a horizontal, infinite PEC ground at z=0 (`GN 1`), with all conductors
and complete lossless dielectric cells separated from it by a documented
positive gap. Dielectric/perfect-ground contact belongs to M11. Real ground
belongs to M13; buried cells, radial screens and cliff ground remain outside
M10. Free-space behavior and the no-DF legacy path must be retained.

Earlier requests to support perfect and finite-conductivity ground together
are superseded by this milestone boundary. Historical real-ground work below
is retained for provenance; its pending checks do not block M10 or M11.

The existing GN flags are:

| Flag | Existing NEC2 meaning | Milestone scope/status |
| --- | --- | --- |
| -1 | Free space | Existing M9 support |
| 1 | Perfect electric conductor (PEC) | Integrated with dense and GMRES, within the subset below |
| 0 | Finite-ground Fresnel reflection approximation | M13; integration and validity policy pending |
| 2 | Finite-ground Sommerfeld/Norton | Preliminary M13 dense/GMRES integration for wires and volumes; real-ground patches excluded |

These distinctions are specified by the [NEC2 GN documentation](https://nec2.org/part_3/cards/gn.html).
The [NEC2 modeling guidelines](https://www.nec2.org/part_3/cards/secii.html)
state that Sommerfeld/Norton is applied to wire interactions, while patch
interactions use the reflection approximation. This is also visible in the
local implementation: `fields.c` electric wire fields add Sommerfeld terms,
whereas patch fields use reflection coefficients. A single exact layered
background for every block cannot simply inherit this mixed approximation.
The first real-ground integration supports GN 2 wires and volumes. It retains
NEC's native Sommerfeld/Norton conductor block and adds spectral layered-medium
volume and cross blocks. Its vacuum limit preserves legacy GN 2 currents and
feed impedance for finite-contrast ground. An exactly matched lower half-space
uses free-space conductor evaluation, avoiding a degenerate native Sommerfeld
limit. Real-ground patches and GN 0 remain rejected: the user selected
GN 2 wires first, with those approximations handled in a separate follow-up.
GN 0 is not silently replaced by an exact Sommerfeld solution.

## Perfect-ground volume operator

With the existing exp(+i omega t) convention, let

- S = diag(1,1,-1) reflect a position through z=0;
- Q = diag(-1,-1,1) reflect an electric current;
- N0(t,s) be the existing cuboid-averaged free-space electric N operator.

For source cuboid s, its reflected cuboid s* has center S center(s) and the
same dimensions. Then

    Npec(t,s) = N0(t,s) + N0(t,s*) Q.

Q acts on source components, hence matrix columns. Horizontal electric
current images reverse sign; vertical images retain sign. There is no extra
factor of two and no additional singular self term: a cell above ground is
separated from its own image. The direct block retains the certified cuboid
self/near treatment. The volume-system equation remains

    Aij = epsilon_r(i) delta_ij I - (epsilon_r(i)-1) Npec(i,j).

`vie_ground.c/.h` implements this operation as `vie_cuboid_pec_n`. It checks
complete-cell bounds and a caller minimum gap in metres, and leaves the
result unchanged on failure. At the M10 checkpoint strictly positive clearance
was required. The M11 extension permits PEC face contact when the caller
minimum is zero; see [M11-PEC-CONTACT.md](M11-PEC-CONTACT.md). Both underlying kernel reports are
available; fixed far-order reports must not be interpreted as adaptive error
bounds. Kernel wavelength and quadrature restrictions remain unchanged.

## PEC integration and supported deck subset

The production implementation now uses the same PEC background throughout:

1. Preserve ground state when constructing the NEC conductor block, replacing
   the M9 free-space reset and setting NEC's PEC image factor explicitly.
2. Add the PEC image block to volume/volume interactions, including diagonal
   blocks, and use NEC's PEC electric fields for conductor-to-volume terms.
3. Add electric-current images to both electric volume-to-wire and magnetic
   volume-to-patch fields. Magnetic image fields come from curl(G Jimage);
   they are not obtained by applying an electric observation-row sign rule.
4. Add the reflected incident plane wave to volume excitation, consistent
   with NEC conductor excitation and the electric boundary condition.
5. Add image radiation and image near fields. Report upper-half-space fields
   only and integrate power over the upper hemisphere. Account for the total
   background excitation in the plane-wave extinction calculation.
6. Validate conductor clearance, including wire endpoints/radius and the
   declared patch geometry; validate observation and pattern domains.
7. Test vacuum identities, reciprocity and power for wire, patch and mixed
   models. Compare against explicitly mirrored free-space reference models.

The implemented deck subset is:

- `GN 1` with no radial screen and only planar `RP 0` (or no RP);
- `-s dense` and `-s gmres`, with existing GMRES residual/iteration controls;
- ordinary voltage sources and plane waves, including elliptical types 2/3
  with an explicit positive axial ratio in EX field F6;
- every declared dielectric box, including vacuum boxes, wholly above z=0;
- every wire endpoint above ground by more than its wire radius;
- every patch centroid above ground by more than `sqrt(area)`, the initial
  conservative support-size policy for the NEC centroid/area abstraction;
  this does not reconstruct the original polygon boundary from NEC's stored
  patch data or establish clearance for arbitrarily elongated polygons;
- incident and far-field theta in [0,90] degrees, and near-field observations
  at z >= 0 with the existing one-cell-edge dielectric observation clearance.

`dielectric_ground_validate` runs before matrix allocation or geometry
normalization. Unsupported combinations fail explicitly before a solution
banner. The explicit elliptical F6 requirement avoids interpreting NEC's
zero axial ratio differently between conductor and volume excitation.
Free-space excitation conventions are unchanged. Far-field complex output
for volume-only decks now uses 12 decimal places to make image comparisons
meaningful below the previous printed precision. Scattered fields exclude
the incident/reflected background wave; the plane-wave extinction integral
uses the complete background electric field.

A runnable small fixture is `tests/dielectric/m10-pec-wire.nec`:

```sh
./nec2c -s dense -i tests/dielectric/m10-pec-wire.nec -o /tmp/m10-pec.out
./nec2c -s gmres -r 1e-10 -i tests/dielectric/m10-pec-wire.nec -o /tmp/m10-pec-gmres.out
```

The existing free-space FFT plan must not be used unchanged for PEC. The
reflected block depends on z_test + z_source, rather than their difference.
A future image embedding can scatter signed mirror currents into an extended
lattice and gather only physical-cell rows. Recompute padded dimensions and
memory accounting before allocation. The PEC deck path currently rejects `-s fft` explicitly, without
silently solving a free-space model. Image embedding remains a later
optimization; the dense PEC path supplies its reference.

## Verification completed in the first increment

`tests/vie_ground_test.c` checks:

- static and dynamic separated cuboid blocks against independent six-dimensional
  Gauss integration of the closed-form electric dipole dyad, including all
  Cartesian couplings;
- volume-weighted transpose reciprocity for unequal cuboids;
- tangential electric field approaching zero linearly as a testing cuboid
  approaches z=0, while the normal field remains finite;
- the direct singular self block plus its nonsingular image correction;
- the static free-space limit as ground separation grows;
- rejection of contact, buried test cells, insufficient clearance, invalid
  input and preservation of the output on failure.

The targeted Automake run passed all three selected programs: the new ground
kernel test, the existing independent cuboid-kernel suite and the legacy NEC
regression suite. The new test/module also passed strict C11 warnings and
AddressSanitizer/UndefinedBehaviorSanitizer checks (LeakSanitizer disabled for
the sandbox). The source distribution build passed and includes these files.

Those initial kernel tests alone do not validate conductor coupling or fields.
The second increment adds `tests/dielectric_ground.py` with explicit mirrored
free-space models and integration checks:

- linear and both elliptical handedness volume solutions equal the sum of two
  free-space plane-wave solves on physical plus reflected dielectric cells;
- lower image currents have reversed horizontal and retained vertical signs;
- volume far and electric/magnetic near fields match the explicit images;
- tangential E and normal H vanish at the PEC boundary;
- driven wire/volume and patch/volume currents and patterns match explicitly
  mirrored conductor and dielectric models (mirror patch winding is reversed);
- vacuum DF recovers NEC's PEC wire, patch and mixed currents/feed results;
- loaded wire, patch and mixed dense/GMRES solutions and fields agree;
- a two-frequency driven sweep matches independent single-frequency solves;
- volume-only plane-wave power and driven-wire upper-hemisphere power checks
  pass, with finite angular quadrature and the inherited conductor
  discretization error retained; this is not an independent physical oracle;
- unsupported modes, FFT, contact/clearance and below-ground output fail.

Mirror comparisons use the free-space operator rather than the PEC image
wrapper. Currents/patterns use scaled tolerances of 2e-7 (3e-7 for the wire
image model) plus 1e-12 absolute, and dense/GMRES overlap uses 2e-8 or tighter.
Legacy patch comparisons round to NEC's printed `%.2E` Cartesian precision;
legacy wire diagnostics are also compared at their printed precision.

## Second-increment validation record

All fourteen regression programs pass. The full run initially found one
stale free-space-only RP diagnostic assertion; after updating it and rebuilding
with the final input guards, the two affected Python suites passed through
the Automake driver and the aggregate report records 14/14 pass, no errors.
The additional PEC elliptical patch vacuum/dense/GMRES check and explicit F6
rejection are included in that final run. The independent existing kernel,
FFT, hybrid, wire, patch and legacy suites remain passing. The entire NEC2C
binary was rebuilt with AddressSanitizer and UndefinedBehaviorSanitizer for
the PEC integration suite; LeakSanitizer is disabled for the sandbox. Strict
C11 warnings/syntax checks pass for the new ground validation module.

The runnable fixture with `-s gmres -r 1e-10` converged to true relative
residual `1.595494884557e-12`. Accepted power was `6.634552638757e-3 W`,
upper-hemisphere radiated power `6.615389605082e-3 W`, relative power defect
`2.888368623924e-3` (0.289%), and angular quadrature change
`1.675115023806e-6`. The defect is comparable to the inherited coupled
free-space discretization defect; it is not claimed as independent physical
validation or as complete M10 power evidence.

## Finite-ground operator increment (2026-10-20)

`vie_halfspace.c/.h` implements a reflected angular-spectrum operator for a
passive, nonmagnetic, homogeneous lower half-space. It is built by the isolated
test target. Initially it was not connected to the production deck solver;
the GN 2 wire integration below now connects it consistently to the volume
and cross blocks, excitation, fields and power. GN 0 and real-ground patches
remain rejected.

For exp(+i omega t), epsilon is `er - i loss`, `er >= 1`, `loss >= 0`.
The outgoing vacuum vertical wave number is `k0*w`, with
`w=sqrt(1-u^2)` for propagating waves and `w=-i*sqrt(u^2-1)` for evanescent
waves. The substrate root uses the same outgoing/decaying branch. Both TE
and TM Fresnel coefficients enter the reflected electric dyad:

    rs = (w - wg)/(w + wg)
    rp = (epsilon*w - wg)/(epsilon*w + wg)
    A = rs s s^T + rp p_up p_down^T

The reflected Weyl integral has scalar weight
`-i*k0^3/(8*pi^2) * (u/w) du dphi`, lateral phase and vertical
`exp(-i*k0*w*(z_test+z_source))`. H comes from the curl of this same spectrum,
not an independently chosen magnetic approximation. The spectral framework
is described in [the layered dyadic Green-function derivation](https://arxiv.org/abs/1610.04479);
the implementation uses the project's opposite time convention explicitly.

Point-current, finite cuboid and mixed point/cuboid source/test APIs are available.
The mixed kernels integrate a finite source cell when testing at a point,
and average the field of a unit point current moment when testing a cell.
These are the building blocks for volume-to-wire fields and wire-to-volume
fields; integrating the NEC wire basis remains separate work. Electric
reciprocity includes the source-cell volume factor; magnetic fields are
computed from the spectral curl rather than inferred by transposition. Cuboid
averages use exact pulse-cell Fourier factors, rather than replacing a cell
with a centroid dipole. Evanescent vertical averaging evaluates endpoint
exponentials so large complex sinc factors do not overflow. The direct
free-space self/near operator remains separate and unchanged.

Quadrature splits propagating/evanescent spectra, scales the evanescent
variable by total bottom clearance, splits at the substrate's lossless branch
and decay scales, and compares successive Gauss/trapezoid orders. Failure at
the maximum order is explicit and leaves outputs unchanged. The evanescent
scaled cutoff is 64. Success is a successive-quadrature check, **not a rigorous
tail/error certificate**, and electrically large lateral separation or very
small clearance can require more than the initial maximum order. This
reference implementation is not a production performance claim.

Verification in `tests/vie_halfspace_test.c` includes normal, Brewster,
grazing and electrostatic Fresnel limits; finite-loss E/H values from an
independent scalar axial Sommerfeld oracle; lossless substrate branch
handling; finite-volume reciprocity; matched-medium zero reflection;
quasistatic loss sign; PEC point E/H and cuboid-image limits; invalid input
and explicit unresolved-integral failure. The independent scalar reference
is reproducible with `python3 tools/m10_halfspace_reference.py`, using
analytic azimuth integration and adaptive Simpson radial quadrature with no
NumPy/SciPy dependency. The PEC point maximum electric/magnetic errors were
`1.86e-12` / `3.73e-13`, and the finite-cuboid electric image error `1.06e-16`
on the tested fixtures. Mixed finite-source E/H also match independent
real-space PEC image integration to `2.33e-14` / `2.79e-15`, and the mixed
lossy electric blocks pass volume-weighted reciprocity. These limits do not replace finite-ground end-to-end
validation.

The two new test programs pass through Automake, giving an aggregate 16/16
passing report with the previous suites retained; unchanged production suites
were not needlessly rerun. The new C kernel/tests also pass strict C11
warnings and AddressSanitizer/UndefinedBehaviorSanitizer with halt-on-error
(LeakSanitizer disabled for the sandbox).

## Laptop OpenEMS staging

OpenEMS is on the user's laptop, while development runs on rxtx. The portable
`tools/openems_m10_ground.py` and [laptop instructions](M10-OPENEMS-LAPTOP.md)
provide a standard-library-only plan, preparation without time stepping, a
one-thread, wall-time-limited timing pilot and a shared wall-budget controller that
terminates worker and solver descendants on timeout. The timeout behavior
is tested without OpenEMS. Model setup and simulation still need execution
in the user's existing OpenEMS environment. The original timing pilot has
now returned and is reviewed below; revised-model runtimes have also returned
in the second batch, while physical accuracy remains under review. Pilot output is not accuracy
evidence, and verification CSVs remain explicitly unreviewed.

## Remaining M10 exit gates

Only PEC evidence gates M10. Finish review of the two independently refined
PEC OpenEMS ground-gap models, check NEC volume/wire refinement, and resolve
how source/feed and boundary-domain errors affect the agreed dielectric-induced
impedance comparison. The available 300 MHz PEC checkpoint differs from the
fine OpenEMS perturbation by 0.3982 ohm, but does not by itself close all
refinement and source/domain questions.

Consolidate PEC free-space/remote-ground limits and coupled reciprocity/power
evidence, including the remaining PEC patch/plane-wave power cases, and record
the supported clearance/backend subset and acceptance evidence. M10 remains
in progress until these PEC gates are satisfied; changing scope is not itself
completion. Image FFT acceleration is a separate capability, not a new M10
gate for the documented dense/GMRES subset.

GN 0, GN 2, substrate absorption, real-ground patches and further real-ground
OpenEMS results are preliminary/remaining M13 work. No real-ground run is
required to complete M10 or M11. See the [M13 handoff](M13-LOSS-AND-REAL-GROUND.md).

### Returned pilot review and revised antenna (2026-10-20)

The original OpenEMS pilot returned from the laptop is timing-only evidence:
1.22394 ps engine timestep, 5,000 steps in 20.29 s, clipped Gaussian excitation
and no energy-decay convergence. Its 1 mm wire/feed anchors and poorly graded
fixed-anchor intervals explain the timestep burden; this was not an
unexplained numerical sliver. The user authorized a thicker benchmark wire:
revision 2 uses 5 mm radius in OpenEMS and generated NEC counterparts, an
explicit 2:1 graded grid, separately controlled transverse refinement, shorter
broadband excitation and a full-source pilot step ceiling under the existing
external battery budget. See [the laptop report](M10-OPENEMS-LAPTOP.md) for
measured original costs, predicted revised mesh sizes and the next bounded
pilot. This does not constitute a passing independent M10 comparison.

### Second-batch verification review (2026-10-20)

Both PEC mesh pairs and the coarse real-ground pair returned normally without
source clipping or step-limit warnings. The fine real vacuum job finished;
the loaded job was killed by the shared wall budget near -58 dB, so its raw-
trace DFT is retained only as a provisional estimate. The independent offline
DFT reproduces completed native port impedances to 1.6e-13 ohm. Review logic
keeps incomplete pairs provisional and tests the magnetic half-step phase.

The PEC dielectric perturbation changes by at most 0.2631 ohm across the two
OpenEMS meshes. A 64-cell, 300 MHz NEC PEC checkpoint differs from the fine
OpenEMS dielectric perturbation by about 0.3982 ohm. However, absolute feed
impedances remain mesh-sensitive by tens of ohms, and NEC wire refinement
alone does not explain their difference. Source/wire discretization and
boundary-domain errors still need investigation. No additional battery-
powered refinement is requested at this point; the single missing loaded
run can be restarted later using the cached vacuum result. The detailed
[laptop review](M10-OPENEMS-LAPTOP.md) records timings, numbers, limitations
and exact resubmission instructions. M10 remains in progress.

### Preliminary M13 GN 2 wire integration, developed during M10 (2026-10-20)

The user selected GN 2 with wires for this integration; GN 0 and real-ground
patches remain separate follow-up work. Positive-gap GN 2 wires and dielectric
volumes now use dense or dense-operator GMRES. The standalone dielectric
plane-wave solver also supports this half-space. Screens, cliff ground,
contact/buried geometry, real-ground patches and ground FFT remain rejected.
Real-ground near-field observations require strictly positive height.

The conductor block retains NEC's Sommerfeld/Norton wire formulation, including
its existing interpolation and thin-wire approximations. Its material is
initialized at every frequency from GN parameters. A negative GN conductivity
specifies fixed imaginary relative permittivity; it is not converted in place
and reused incorrectly at later frequencies. NEC's Sommerfeld interpolation
cache is invalidated when the complex ground permittivity changes, even if
an observation reuses a previous interpolation stencil. Free-space-only field
evaluation also suppresses the post-loop Sommerfeld correction explicitly;
changing only the image-loop count was insufficient in the legacy routine.
This prevents counting the reflected field twice in spectral cross blocks
and electric near fields.

The volume block combines the established singular free-space cuboid kernel
with the finite-cell reflected spectrum. Volume-to-wire fields use point-tested
finite source cells. Wire-to-volume reflected fields integrate NEC's actual
constant/sine/cosine wire current basis with eight-point segment quadrature
against a cell-averaged point-source spectrum. The free-space wire contribution
retains native NEC's finite-radius fields and six-point-per-axis testing-cell
average. Reflected constant/sine/cosine basis fields are cached per cell and
segment before applying NEC's current-basis coefficients, avoiding repeated
Sommerfeld evaluations for different conductor unknowns. This optional cache
is capped at 16 MiB and falls back to uncached evaluation if allocation fails
or the cap would be exceeded. It uses 101,376 bytes for 64 cells and 11 wires.
It is freed after each solve, including failures; frequency sweeps rebuild it.
These fixed spatial rules inherit the existing need for wire/cell
refinement; radial/angular spectral convergence is checked separately.
Both real-ground electric and magnetic conductor near fields combine native
free-space fields with reflected spectral wire fields; the legacy approximate
magnetic reflection is not reused. Plane-wave reflection and upper-hemisphere
volume radiation use the same TE/TM Fresnel coefficients as the Green dyad.

Power delivered into the lower half-space is computed independently through
Parseval's theorem applied to downward interface Poynting flux. It includes
propagating and evanescent spectra and forms the combined wire-plus-volume
source amplitude before squaring, retaining interference. Cell transforms
and transforms of the NEC sinusoidal wire basis are analytic; stable endpoint
exponentials avoid evanescent overflow. Fresnel flux factors are simplified
algebraically to avoid cancellation near grazing incidence and the PEC limit.
A finite angular quadrature can still fail for an extremely narrow grazing
transmission layer; failure is explicit, not reported as zero ground power.

For lossy infinite ground this delivered power is eventually absorbed in the
substrate; for lossless ground it is transmitted power. Output therefore uses
`LOWER-HALF-SPACE POWER` rather than calling every such term absorption.
Driven-source power checks compare accepted power with air radiation plus
lower-half-space power. Plane-wave checks compare extinction with scattered
air radiation plus scattered power into the substrate. The substrate's
unperturbed incident-wave absorption is not part of this scattered-power check.

The independent vertical-dipole reference compares local source work with
upper-hemisphere radiation and interface flux; the tested finite-loss balance
agrees within 2e-9 W (about 3e-11 relative). Matched media transmit exactly half
of the corresponding free-space vertical-dipole radiation. These checks are
independent of the coupled linear-system residual and of accepted-minus-air
power. They do not certify every spatial quadrature or close M10's independent
physical verification gates.

Runnable fixture: `tests/dielectric/m10-real-wire.nec`.
Integration checks: `tests/dielectric_real_ground.py`.

For reproducibility, let `D(u,phi)` be the complete downward source transform,
`w=sqrt(1-u^2)` with the outgoing branch, and `v=sqrt(epsilon-u^2)` on the
substrate branch. With `Ds=s^T D` and `Dp=p_down^T D`, the simplified interface
flux implemented in `vie_halfspace_power` is

    P_lower = eta0*k0^2/(8*pi^2) * integral u du dphi [
        Re(v)/|w+v|^2 * |Ds|^2
        + Re(v*conj(epsilon))/|epsilon*w+v|^2 * |Dp|^2 ]

The transform convention is
`D = integral J(r) exp(+i*k0*u*rhohat.r_xy-i*k0*w*z) dV`.
For wires, replace `J dV` with `I(l) tangent dl`. The propagating substitution
is `u=sin(theta)`; the evanescent substitution is
`u=sqrt(1+(t/(k0*clearance))^2)`. No accepted-power difference enters this
integral. Extremely small clearance, large lateral extent or sharp grazing
transmission can exceed the supported quadrature order; such cases fail rather
than silently treating the unaccounted power as ground absorption.

Verification after integration: the aggregate Automake report is 17/17 passing.
Affected legacy, free-space and PEC solver/field/wire/patch suites were rerun;
unchanged parser/mesh and standalone iterative/FFT kernel reports were retained.
The new GN 2 integration suite checks native-vacuum compatibility, dense/GMRES
agreement, EX 1/2/3 power, matched media, lossless transmission, conductivity
and fixed-imaginary-permittivity sweeps, and Maxwell's near-field curl relation.
The small driven fixture's power defect is `2.675274284704e-4` (0.0268%);
coupled EX 1/2/3 defects are below `1.710e-3` (0.171%). The total scattered
E/H curl check passes at relative tolerance `3e-4`, including the precision
limits of NEC's inherited physical constants. These are discretized-model
checks, not a global physical accuracy bound.

Strict C11 warnings, AddressSanitizer and UndefinedBehaviorSanitizer checks pass
for the new half-space/flux test. The full instrumented executable also passes
GN 2 driven E/H observations, elliptical plane-wave power and matched-medium
compatibility, with halt-on-error and sandbox-incompatible LeakSanitizer disabled.
The source distribution builds successfully.

The reflected-basis cache agrees with the uncached executable for driven and
elliptical plane-wave currents, feed impedance and far fields (relative checks
at 1e-10, feed at 1e-11). The one-cell/11-segment fixtures used 2.12/2.15 CPU
seconds without caching and 1.02/1.03 seconds with caching on this VM; these
small-fixture timings are not a large-model performance guarantee.

The optional cache's allocation, indexing and cleanup also pass sanitizer
checks with the current wire translation unit instrumented. Its fallback
operator was compared with the saved uncached executable; no numerical
change is claimed from this assembly optimization.

After adding the cache, the GN 2 suite and aggregate report were rerun and
remain 17/17 passing. A separate two-cell case also passes cache indexing and
cleanup checks with the wire translation unit instrumented by ASan/UBSan.
The larger 64-cell, 300 MHz cloud NEC counterpart has completed; the reviewed
checkpoint below is supporting evidence, not full M10 physical acceptance.

### Completed preliminary M13 GN 2 physical checkpoint (2026-10-20)

The 64-cell, 300 MHz counterpart to the 5 mm-radius laptop benchmark completed
with true GMRES residual `2.298610538446e-11`. The wire-basis cache used
101,376 bytes; solver CPU time was 84.48 s (not wall time). NEC's vacuum feed
impedance was `79.08469476773 + j21.08074897213` ohm; loading the cube gave
`78.71405222046 + j19.89813918016` ohm, a perturbation of
`-0.37064254727 - j1.18260979197` ohm. The completed coarse OpenEMS pair gave
`-0.36782103506 - j1.45832265030` ohm. Their complex perturbation difference
is **0.27573 ohm**, below the agreed 2-ohm allowance at this checkpoint.

The loaded model's accepted power was 5.970569 mW, air radiation 2.599191 mW
and lower-half-space power 3.376606 mW, with relative defect `8.755894e-4`
(0.0876%). Metadata and results are recorded in
[`tests/m10-nec-real-300.json`](../tests/m10-nec-real-300.json). This is one
frequency, one VIE mesh and one completed real-ground OpenEMS mesh. The fine
loaded OpenEMS result remains pending, and absolute feed/source and boundary
convergence are unresolved. This is preliminary M13 evidence; real-ground
acceptance does not gate the revised PEC-only M10.

### PEC closure checks, 2026-10-20 (complete)

`tools/m10_pec_acceptance.py` stages reproducible GMRES checkpoints with a
pinned-binary SHA-256 manifest and input-aware completed-result reuse.
`tests/m10-pec-acceptance.json` retains the physical checks and comparison
results. The frequency sweep uses the same 5 mm radius, 475 mm span, extended
wire kernel and positive-gap epsilon-3 cube as the laptop benchmark.

At 280/300/320 MHz the magnitude of the difference between NEC and fine-grid
OpenEMS dielectric-induced impedance changes is respectively
0.1894/0.3982/0.8218 ohm. All are below the agreed 2 ohm allowance. The 300 MHz
wire refinement from 11 to 21 segments changes delta Z by 0.0446 ohm and
reduces the loaded driven power defect from 0.222% to 0.0188%. These compare
**dielectric-induced changes**, not the absolute impedances of the NEC
ideal delta-gap source and OpenEMS finite lumped feed.

`tests/m10_pec_validation.py` adds the following physical integration gates:

* Transfer-admittance reciprocity between unequal, offset wire ports:
  relative discrepancy 7.82e-5 for vacuum and 7.57e-5 with dielectric. This
  physical port check inherits NEC collocation/basis error; it is distinct
  from the stricter dyadic-kernel reciprocity identity.
* Translating the driven wire and dielectric away from the ground plane
  reduces the feed-impedance difference from free space from 2.875 ohm to
  0.0450 ohm at +4 m and 0.00345 ohm at +16 m. Driven excitation avoids the
  nondecaying reflected incident wave of a plane-wave remote-ground test.
* Linear and both elliptical plane-wave polarizations conserve power to
  within 0.157% for tested wire, closed-patch and mixed geometries, with and
  without dielectric. A refined closed-patch case has a 0.127% defect.
  Patch extinction integrates native Cartesian sheet currents with the same
  centroid/area approximation as NEC radiation; wire extinction integrates
  its constant/sine/cosine basis analytically. All sources and PEC-reflected
  incidence contribute before comparison with total scattered power.

Physical patch power checks use outward-facing **closed surfaces**. An open
sheet fixture remains useful for algebraic image/legacy-compatibility tests,
but does not supply a physical power validation of NEC's magnetic-field
surface equation. Patch discretization defects need not decrease
monotonically; the power checks bound them independently at two resolutions.

The volume refinement from 27 to 64 to 125 cells changes delta Z by
0.00620 ohm and then 0.00288 ohm. The 125-cell solve has a true GMRES
residual of 7.84e-11 and a 0.222% driven power defect. The aggregate
`make check -j2` passes all 18 test groups, including the new PEC physical
validation. The returned laptop sensitivity evidence below closes the last
PEC gates; the revised M10 scope is complete.

### Returned PEC feed/domain sensitivity: accepted (2026-10-20)

All four new OpenEMS jobs completed normally at the requested -60 dB energy
criterion, well before the 120000-step ceiling, without source clipping.
The larger-domain vacuum/loaded jobs took 23.97/23.93 engine seconds; the
shorter-feed jobs took 16.86/16.86 seconds. Port-trace DFT independently
reproduces the returned impedances using each voltage/current trace's own
timestamps. Periodic progress lines precede the final stopping criterion;
the audit records the explicit final decay rather than treating the last
periodic energy as the final energy.

| Frequency (MHz) | Domain delta-Z change (ohm) | Feed delta-Z change (ohm) |
| --- | --- | --- |
| 280 | 0.00282 | 0.00222 |
| 300 | 0.00397 | 0.00284 |
| 320 | 0.00574 | 0.00629 |

Changes are complex magnitudes against the original 20 mm/5 mm-spacing
baseline, with only the specified domain or feed choice varied (and mesh
anchors following the source geometry). NEC's delta-Z discrepancy against
these variants is at most 0.555 ohm. Combined with the original independent
mesh refinement and the fine-grid NEC comparison (maximum 0.822 ohm), these
results satisfy the agreed 2 ohm physical acceptance gate. They demonstrate
sampled sensitivity, not a rigorous error bound for arbitrary geometries.

The retained audit is `tests/m10-pec-sensitivity-review.json`, with raw archive
SHA-256, explicit decay results, mesh metrics, and per-frequency comparisons;
`tests/m10-pec-sensitivity-pairs.csv` is the compact numerical table. The
original archive remains `/tmp/m10-pec-sensitivity-results.tar.gz` on rxtx.
The current M10 PEC scope is complete. Earlier pending-gate statements in
this document describe historical checkpoints. M11 PEC-ground contact and
M13 real-ground/loss validation remain separate work.
