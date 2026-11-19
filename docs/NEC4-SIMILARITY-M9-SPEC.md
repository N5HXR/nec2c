# NEC4 similarity M9 execution specification

Status: implementation started on 2026-11-17 after user acceptance of M8.

## Purpose and boundary

M9 turns the M8 direct Sommerfeld oracle into an accelerated ground-field
evaluator for matrix fill.  The target is the documented NEC-3/NEC-4
Sommerfeld/asymptotic method, not a reconstruction from a licensed binary or
source tree.

The implementation remains split into independently testable layers:

1. direct samples from M8;
2. near-region transformed tables and linear interpolation;
3. outer-region model-based least-squares interpolation;
4. first- and second-order asymptotics, including the surface-wave term;
5. a deterministic dispatcher and an exact-key persistent cache; and
6. production C integration and parallel table generation.

No accelerated value is accepted solely because it is self-consistent.  Every
layer is compared with M8 before it can feed the next layer.

## Public equation contract

The equation lock is based on Section 5.1.3, Section 5.2, and Appendix J of the
public NEC-4 Theory Guide:

- equations (5-9a-e): singular and finite remainders;
- the three source-depth table regions and transforms `F1`, `F2`, and `F3`;
- the published trilinear interpolation equation;
- equation (5-10), its 36 ray-basis terms and eight supplemental terms;
- Figure 5-4 and Table 1 sample lattices;
- equations (5-11) through (5-23): steepest-descent and first-order fields;
- the modified saddle-point surface-wave treatment; and
- Appendix J second-order corrections and their small-rho interpolation.

The guide's symbols are retained in derivation comments.  Code-facing names
also state the medium and physical units so that `k1`, `k2`, `z`, and `z'`
cannot silently change meaning.

## Internal implementation gates

### Gate A - coordinates, regions, transforms, and interpolation primitives

- Lock medium 1 as the lower medium and medium 2 as the upper medium.
- Lock the three source-depth ranges, `Lg`, `Lmax`, and the near-table extent.
- Implement reversible `F1`/`F2`/`F3` transforms.
- Implement the published trilinear interpolation formula.
- Generate the Figure 5-4/Table 1 general and surface sample lattices.

Exit condition: boundary ownership is deterministic, every transform round
trips to floating-point precision, affine complex fields interpolate exactly,
and sample counts are 32, 76, and 88 for regions 1, 2, and 3.

Gate A passed on 2026-11-17.  Transform round-trip error was
`1.11e-16`, affine interpolation error was zero, and the three published
sample counts were reproduced exactly.

### Gate B - direct-sampled linear tables

- Generate all five canonical components from M8 on each near-region grid.
- Apply the equation (5-9) singular subtraction only in table 1.
- Store transformed values and reconstruct physical fields after interpolation.
- Sweep every cell interior, face, edge, and shared source-depth boundary.

Exit condition: null-aware error versus M8 meets a frozen per-component budget;
there are no ownership gaps or discontinuities at 0.4 `Lg` and `Lg`.

Gate B1 passed on 2026-11-17 for a representative lossy-ground table-1 cell.
Seven successive cell halvings reduced the worst null-aware error from
`3.462e-1` to `8.228e-5`; every observed convergence ratio was at least 4.000.
Gate B2 independently exercised the table-2 and table-3 phase transforms:
their four-level sweeps reached `2.491e-4` and `1.944e-4`, respectively, with
second-order convergence.

Gate B3 passed on 2026-11-17 across lossless `epsilon_r=4`, average lossy
soil, and seawater. It exercises all three storage transformations, restores
the analytic singular field around table 1, and approaches the `0.4 Lg` and
`Lg` source-depth boundaries independently from both sides. The worst
null-aware interpolation error was `5.917e-4`; the worst independently
restored boundary jump was `4.941e-4`. Boundary ownership remains with the
lower-numbered table exactly at a shared coordinate.

### Gate C - model-based interpolation

- Implement the complete equation (5-10) basis without pruning terms.
- Solve the overdetermined complex system by QR decomposition.
- Fit the documented general and surface samples independently per component.
- Characterize interpolation and extrapolation through the documented outer
  extent of 3 lambda0 in rho, 2 lambda0 in z, and 0.8 `Lg` in source depth.

Exit condition: rank, conditioning, coefficient residuals, and held-out M8
errors are reported for every component and subregion.

The public text contains an ambiguity that must be resolved by this gate rather
than guessed: equation (5-10) prints 36 `A_lmn` terms and eight `B_l` terms,
while Table 1 supplies only 32 region-1 samples and says the eight `F_l`
functions are dropped there.  This conflicts with the preceding statement that
the fit is overdetermined.  The implementation will expose matrix rank and test
candidate term restrictions or rank-revealing solutions against held-out M8
points.  The selected rule will be labeled as an independently validated
implementation choice, not an undocumented NEC-4 fact.

Gate C characterization completed on 2026-11-17.  The ordinary second-ray
phase, `exp[-j(beta3*z+k1*R2)]`, reproduced the guide's lossless-dielectric
`E_z^V` witness with `2.741%` RMS training error.  Treating the scanned
trailing `1/2` as either a square root of the phase or of the exponential gave
more than `36%` RMS and was rejected.  The nominal full-rank solve was not
stable between samples: at least one component exceeded 100% RMS on the
calibration points.  Rank-revealing QR selection produced component ranks
25--28 with a measured envelope of `12.627%` RMS on the documented lattice
and `8.702%` RMS on the initial calibration set.

This is a successful reconstruction and limitation measurement of the legacy
model, not a production-accuracy acceptance.  The 1984 public paper confirms
that deployed variants used separately tuned 32- and 28-term regional models,
but does not identify the individual retained terms.  The later public GN3
report explicitly replaced this model-based parameter estimation with
adaptive cubic interpolation.  M9 will retain the documented model as a
diagnostic/reference option and will not silently claim that an independently
selected rank is the original NEC-4 tuning.

### Gate D - asymptotic evaluator

- Implement first-order steepest descent independently of the table path.
- Add the modified saddle-point pole/surface-wave contribution.
- Implement Appendix J second-order corrections and the documented small-rho
  interpolation of higher-order terms.
- Retain a first-order-only diagnostic mode.

Exit condition: both paths reproduce M8 over expanding distance and angle
grids, including grazing, total-reflection, small-rho, and low-loss cases.

Gate D1/D2 equation and branch work completed on 2026-11-17.  The generic
first-order result exposed two scan/implementation hazards: equation (5-22)
contains `sqrt(rho/lambda_s)` in the denominator, and the complex `Gamma`
roots must be continued from the outgoing real-axis branch rather than chosen
independently by instantaneous real part.  The latter otherwise jumps sheets
as a lossy saddle moves just below the real lambda axis.

The Appendix-J second-order reference evaluates the universal `Q0+Q2`
formula with high-precision derivatives of the printed two-term Hankel
amplitudes.  At `rho=2.8 lambda0`, it agrees with M8 to `1.189%` or better for
all five lossless-dielectric components and `1.151%` or better for all five
lossy-soil components.  The best errors are below `0.04%`.

A separate two-medium direct run at `rho=4 lambda0`, `z=|z'|=0.02 lambda0`
retains errors above 100% for at least one component even with the published
second-order uniform pole correction.  This is the documented original-GN2
grazing limitation: the first-order lower-medium ray vanishes and its
interference is lost.  The public GN3 report explicitly identifies this
failure and replaces it with separated branch-cut interpolation and a
higher-order interface approximation.

The accepted M9 policy is fidelity to the documented original NEC4 method:
the evaluator will return that legacy result rather than silently substitute
an M8 direct integral or a later GN3-inspired correction.  When the dispatcher
can identify an interaction in the characterized grazing-risk envelope, it
will emit a deduplicated human-readable warning and a machine-readable
diagnostic.  The warning concerns source and observation points whose direct
ray is within 15 degrees of the interface in the legacy asymptotic region.
The 15-degree boundary is an independently qualified conservative diagnostic,
not an undocumented NEC4 switch: M8 comparisons found severe errors through
11.31 degrees and substantial improvement by 21.80 degrees. A short wire
segment by itself is not a trigger;
`0.02 lambda0` in the qualification witness is the source and observation
height/depth, not the segment length.  M11 owns any accuracy-first NEC4.2-like
remedy.

### Gate E - dispatcher and switch audit

- Select direct-table, model, or asymptotic evaluation from documented extents.
- Evaluate both neighboring methods in overlap bands used only for auditing.
- Freeze switch ownership and report value and derivative jumps.
- Detect the characterized original-NEC4 grazing-risk envelope without
  changing the selected numerical path.
- Deduplicate its warning per job and expose the condition in machine results.

Exit condition: no unexplained discontinuity and no accuracy regression hidden
by a switch; the grazing witness warns and continues with the documented NEC4
result, while short segments away from the interface do not warn merely due to
their length.

Gate E routing and diagnostic policy passed on 2026-11-17. Figure 5-4 boundary
ownership is deterministic, the linear/model/asymptotic extents have direct
unit coverage in Python and C, and an expanded M8 sweep places the conservative
grazing warning at a 15-degree direct-ray angle. The ordinary second-order
saddle remained severely wrong through 11.31 degrees and improved to 8.9--11.2
percent vector error at 21.80 degrees. The sweep also demonstrated that the
uniform pole correction must be limited to its intended near-interface regime;
applying it unconditionally at deeper points is invalid.

### Gate F - cache, production integration, and parallel generation

- Key caches by normalized material properties, frequency, extents, precision,
  equation version, and byte format.
- Validate checksum, dimensions, endianness, and version before loading.
- Require fresh and cached evaluations to be bit-identical.
- Integrate the accepted evaluator into the NEC4 ground matrix path.
- Parallelize deterministic table samples only after serial acceptance.

Exit condition: the complete M9 acceptance grid and user witness commands are
ready, with classic NEC2 behavior still selectable from the command line.

Gate F cache infrastructure passed its first slice on 2026-11-17. The binary
cache uses an exact material/frequency/axis-layout/precision/equation key,
validates grid dimensions and their payload product, format, host endianness,
floating-point width, payload size, trailing data, and checksum, writes by
atomic replacement on Linux, and reproduces serialized complex values bit for
bit. Invalid keys, truncation, incompatible axes, and corrupt payloads are
rejected. Production field integration remains open: the existing NEC2 `somnec`
contour and complex-Bessel machinery is reusable, but its integrand assumes an
upper-medium source and observer. It must be connected to the accepted M8
cross-interface equations rather than reused as an NEC2 reflection-field
substitute.

The second Gate F slice adds a reentrant production-C Bessel-contour evaluator
for the documented lower-source/upper-observer equations. The legacy complex
`J0` implementation was moved without changing its arithmetic into a shared,
thread-local module. Against independent high-precision M8 witnesses for a
lossless dielectric, average lossy soil, and seawater, all five canonical
components agree within `6.1e-6` null-aware relative error. The classic
regression corpus remains unchanged. The C path also implements both
Appendix-H subtractions and analytic add-backs for table 1; its finite
remainder and reconstructed total field meet the same independent error band.
Four-thread repeated evaluations are bit-identical to serial evaluation, so
the accepted direct sampler is ready for deterministic parallel table
generation.

The direct-table layout was independently qualified rather than inferred from
an undocumented NEC-4 tuning constant: the public guide specifies the table
regions and transformations but not the grid increments. The frozen maximum
steps, in units of `Lg`, are `0.02/0.02/0.02` in rho, observer height, and
source depth for table 1; `0.03/0.03/0.015` for table 2; and
`0.03/0.03/0.03` for table 3. The finer table-2 depth axis is measured rather
than aesthetic: a uniform `0.03 Lg` layout missed the gate at a shallow,
long-range average-soil probe, while refining only depth reduced the error
without multiplying the whole table by eight.

Across 39 strictly off-grid probes spanning lossless `epsilon_r=4`, average
soil, seawater, all three depth regions, shallow/Hankel geometries, and the
deep-table limit, the worst null-aware error is `4.090014e-3` against a
`5e-3` gate. Raw five-component complex payloads are 15.38 MiB, 18.28 MiB,
and 51.74 MiB for those materials. Layout qualification uses the Bessel
representation everywhere so it measures interpolation rather than contour
cost; the independent M8 gate already establishes Bessel/Hankel equivalence
in overlap and shallow cases. The C Hankel contour, C asymptotic evaluator,
physical axis limits, and segment integration remain before the production
matrix path can be enabled.

The third Gate F slice implements the production table representation around
that frozen layout. It allocates each region with the exact qualified axis
spacing, gives every sample a fixed contiguous five-component slot, fills the
grid with a reentrant callback under static OpenMP scheduling, performs the
published eight-corner interpolation, and reads or writes the complete table
through the exact cache key. A complex affine witness interpolates to
roundoff; one- and four-worker generation are byte-identical; the average-soil
three-table allocation is exactly 19,167,680 raw bytes; and a cache reload is
bit-identical.

The fourth Gate F slice closes the complete table-1 interface plane. Appendix
H requires the analytic `V0` field and the residual `1/R` terms to be removed
before the interface value is formed; Appendix I further shows why simply
making the Figure 5-3 rays vertical is not equivalent when both endpoints are
on the interface: the source saddle becomes a branch-cut integral from `k1`,
and a second branch-cut contribution can enter at `z=0`. The raw exact-height
Hankel entry points therefore reject `z=z'=0`. Production C instead forms the
finite Appendix-H `E''` value by first-order Richardson extrapolation from
three positive-height samples and independently checks it along a symmetric
and a source-only approach.

Every non-origin point on the frozen interface lattice was exercised: 50 rho
values from `0.02 Lg` through `Lg` for lossless `epsilon_r=4`, average soil,
and seawater. Across all 150 points the worst refinement estimate was
`4.391e-6` of the local five-component field scale and the worst independent
approach-path difference was `4.079e-6`. The ordinary Bessel/Hankel overlap
remains `6.220e-5`, and the high-precision M8 Hankel gate remains at
`1.694e-13` overlap. The transformed coincident origin remains the analytic
zero already specified by Section 5.1.3 rather than an extrapolated physical
field. The square-root continuation is initialized with Appendix H's
large-`|lambda|` rule, `gamma ~ S lambda`; this removed otherwise discrete
sheet changes in narrow radial bands.

The fifth Gate F slice connects Appendix H to the production table container.
`nec4_ground_table_fill_direct` now selects the accepted Bessel contour,
Hankel contour, or physical interface limit at each grid point, stores table
1 as the finite `E''` remainder multiplied by `F1`, and stores tables 2 and 3
as their physical fields multiplied by `F2` or `F3`. The corresponding table
evaluator interpolates the stored quantity, removes the transform, and for
table 1 adds the residual `1/R` field followed by the analytic `V0` field.
The coincident transformed origin remains zero and physical evaluation at
that singular point is rejected.

A reduced production table using the frozen `0.02 Lg` table-1 spacing is
byte-identical with one and four generation workers. Interface storage agrees
with the independently qualified interface evaluator to roundoff, grid-node
reconstruction agrees with direct Appendix-H evaluation, and a fully off-grid
cell-center physical field differs from direct evaluation by `5.119e-5`, well
inside the `5e-3` table target. Direct generation also verifies that the
supplied wave numbers match the material and free-space wavelength encoded in
the cache key, preventing a numerically valid table from being cached under
the wrong material identity. Production matrix integration remains disabled.

The callback boundary was deliberate in this slice: the then-current direct
integrator was qualified away from the cylindrical axis, while a complete
production grid includes rho=0 and the coincident table-1 origin. The next
slice derives those limits rather than hiding them with an epsilon offset.

The sixth Gate F slice derives those cylindrical limits from
`J0(x)=1-x^2/4+O(x^4)` and applies the corresponding `J0'`, `J0''`, and
`J0'/rho` limits directly in both the ordinary and Appendix-H integrands. It
also evaluates the `rho -> 0` limits of every cancellation-resistant residual
term analytically. Against evaluations at rho=`1e-3` and `5e-4` metres, the
two components that vanish on axis scale by 0.500 while all surviving direct
and finite-remainder components converge quadratically by approximately
4.000. Separate witnesses put the source exactly on the interface and the
observer exactly on the interface; the ordinary and subtracted/reconstructed
paths agree there within `4.289e-11` null-aware relative error.

The coincident physical field at rho=z=source-depth=0 remains correctly
rejected as singular. Its table-1 stored value is different: after the
published singular subtraction and `F1` transformation its analytic limit is
zero. The parallel table fill writes that zero directly and a test callback
that deliberately rejects the origin proves that no hidden epsilon-offset
sample is requested.

The seventh Gate F slice connects the accepted canonical field to the wire
matrix boundary without enabling an incomplete user-facing solve.  A new
matrix adapter applies Appendix G equations (G-6) through (G-20) to construct
all four source/observer pairings (`12`, `21`, `11`, and `22`), rotates the
five canonical vertical/horizontal-source components into the global
Cartesian dyadic of equation (5-1), and integrates the constant, sine, and
cosine current bases over the source segment.  The adapter can consume the
three production tables directly and gives the lower-numbered table ownership
at shared source-depth boundaries.

The matrix decomposition deliberately retains the accepted finite-radius
NEC4 direct-wire kernel.  The new segment integral supplies
`complete-interface field - upper homogeneous field`, so upper-medium direct
self terms are not replaced by a singular filament Green function.  Lower
and crossing pairings receive the corresponding material/transmission
correction through the same expression.  Matching points on a source axis
are displaced to the wire surface consistently with the established wire
kernel.  A thread-local direct-only request avoids mutating shared ground
state during parallel matrix fill.

Against independent high-precision M8 fixtures for lossless
`epsilon_r=4`, the Cartesian matrix adapter's worst null-aware error across
all five canonical components and all four pairings is `7.998e-5`.  The
identical-media null is `2.119e-7`, and splitting a nonidentical-medium source
segment into two halves changes the reconstructed constant/sine/cosine
integrals by at most `2.332e-16`.  The call is now present at
`wire_operator_efld`, but the matrix context is intentionally left unset and
the M7 reference-only execution guard remains active.  Enabling the CLI path
still requires per-frequency table/cache construction plus accepted
model/asymptotic coverage outside the three linear tables; no NEC2 ground
field is silently substituted.

The eighth Gate F slice supplies that frequency-scoped table and context
lifecycle. It derives `Lg=min(|k2/k1|,0.3)*lambda0` and
`Lmax=min(15|k2/k1|,1)*lambda0` from the active complex material, allocates
all three accepted table layouts, and either loads each exact-key cache entry
or generates and atomically stores it. Cache filenames use a stable hash of
the complete serialized key only as a lookup aid; the cache header still
checks every frequency, material, scale, axis, precision, equation-version,
and payload field exactly, so a hash collision cannot accept the wrong table.
Nested private cache directories are created with mode `0700`.

Once all three tables are ready, the runtime publishes one immutable matrix
context whose canonical callback owns the accepted shared-depth boundary
rule and whose material/wavenumber values are the same values used to build
the cache identity. Re-preparing the runtime destroys the prior frequency's
tables before attaching the next context. The medium state can likewise be
updated per frequency, including the complex interface charge factor, without
reclassifying unchanged geometry.

The lifecycle test allocates `17,011,680` bytes for its lossy witness. Its
first run records zero hits and generates three tables, an exact repeat loads
all three tables bit-for-bit, and a one-percent frequency change records zero
hits and generates a separate three-table bank. A cached matrix-context
evaluation is bit-identical to the fresh context. The runtime is not yet
installed by the CLI frequency loop: doing so before the accepted outer
model/asymptotic callback exists would turn legitimate interactions outside
the linear-table box into a late matrix-fill failure.

The ninth Gate F slice ports the documented outer asymptotic field to
production C without numerical differentiation. Truncated complex Taylor
jets carry the phase through fourth derivative and each component amplitude
through second derivative, reproducing the universal Appendix-J `Q0+Q2`
correction while preserving the continuously selected outgoing square-root
sheet. The first-order saddle, second-order saddle, first-order pole, and
second-order pole are separately callable and separately tested before being
composed into the uniform field. No complex-error-function library is
required: the modified `Q` function uses a convergent long-double series near
pole/saddle coalescence and a scaled asymptotic expansion farther away.

Against independent high-precision lossless and lossy fixtures, all four
terms agree within `5.444e-10` relative error; that common floor is set by the
rounded wave numbers serialized in the portable C fixtures. Appendix J's
near-vertical-incidence rule is implemented literally: below
`rho_A=0.02|z-z'|`, only the higher-order correction is interpolated between
the analytic cylindrical-axis limit and the ordinary second-order value at
`rho_A`, with odd/even component symmetry and the local phase restored after
interpolation. The analytic axis witnesses pass and the two branches join at
`rho_A` to roundoff.

This closes the asymptotic equation implementation but deliberately does not
enable the CLI yet. The public guide does not identify the deployed 32- and
28-term selections for equation (5-10), and the independently reconstructed
rank-revealing model reaches `12.627%` RMS in its documented region. The
remaining production decision must therefore be explicit: accept that
characterized legacy-like fit, use the accurate but slower direct contour in
the model region, or undertake a later adaptive-interpolation method. None of
those choices will be hidden behind the name of the original NEC-4 tuned
model.

The tenth Gate F slice makes that model-region choice executable and
observable rather than implicit. The frequency runtime now owns the production
callback: linear-table routes evaluate the cached table bank, asymptotic
routes evaluate the new second-order saddle plus pole path, and equation-(5-10)
model routes use the accepted direct Bessel or Hankel contour. The latter is a
correctness-first replacement for the unreconstructable tuned fit, not a claim
to possess NEC-4's missing coefficients, and the dispatch record marks every
such evaluation with `used_direct_fallback`. Runtime tests exercise all three
outcomes, and fresh/cached table evaluations remain bit-identical after the
callback change. This completes the numerical coverage needed before the
frequency-loop and user-diagnostic integration; model-region speed remains a
measured optimization target rather than a reason to accept a less accurate
field.

The follow-on optimization will use the direct contour as an oracle and will
not inherit the undocumented equation-(5-10) term selection. Candidate fields
will first apply the accepted `F1`/`F2`/`F3` phase and singularity transforms,
then compare regular transformed grids, local polynomial cells, and an
adaptive hierarchical grid. Training and withheld samples remain disjoint.
Selection reports must include null-aware component error, five-component
vector error, switch and cell-boundary jumps, cache size and construction
time, scalar and parallel evaluation throughput, and worst-case refinement
depth. The first target is the existing `5e-3` interpolation gate over the
entire documented model region, with no boundary jump larger than that error.

The direct path remains the default until a candidate passes on lossless
dielectric, average soil, seawater, low-loss, grazing, and near-axis witness
sets. A promoted representation will have a separately versioned cache key
and a command-line `direct` escape hatch. Failure to meet the requested error
at runtime falls back to the oracle rather than returning an extrapolated
field.

The eleventh Gate F slice completes that experiment and promotes an
independently measured model-region representation. A screen of 512 global
radial-basis candidates on the historical sparse lattice failed decisively:
the best retained `47.336%` RMS vector error. Regular trilinear grids were then
tested with raw fields and all three published storage transformations. `F2`
was consistently superior. The experiment corrected an initial moving-target
validation design and now uses the same 4096-point deterministic Halton set at
every resolution; no validation point becomes a nested-grid node.

On that fixed set the level-3 `F2` grid showed second-order convergence. A
conservative direct-contour strip for `z < 0.16 lambda0` isolates the steep
near-interface behavior. Outside that strip, lossless `epsilon_r=16`, lossless
`epsilon_r=4`, low-loss soil, and seawater produced `0.008%--0.026%` RMS vector
error and `0.284%--0.433%` worst vector error. The worst absolute error of any
individual component, normalized by the local five-component field peak, was
`0.240%--0.337%`. Average lossy soil already passed on level 2 at `0.040%` RMS
and `0.488%` worst vector error; production nevertheless uses level 3 for one
material-independent rule.

A center-certified adaptive refinement was rejected rather than promoted.
Even a `0.1%` center threshold missed an off-center `1.046%` error, proving
that the proposed estimator could not certify a leaf. The accepted grid is
therefore deliberately simple: `97 x 60 x 33` nodes over rho `0..3 lambda0`,
z `5/32..2 lambda0`, and source depth `0..0.8 Lg`, occupying 15,364,800 bytes
per frequency. Ordinary-material oracle builds took about 10--11 seconds on
16 workers; the seawater stress build took 124 seconds. Exact-key caching
amortizes repeated material/frequency use.

The production runtime now allocates, fills, caches, and evaluates this fourth
table alongside the three direct near tables. A model-region request at or
above the strip uses the constant-time `F2` lookup; a request below it retains
the accepted direct Bessel/Hankel contour and marks the fallback. The runtime
cache payload for its test material is 32,376,480 bytes across all four tables,
and fresh versus cached evaluations remain bit-identical. The direct path
remains the required escape hatch for later CLI integration.

The twelfth Gate F slice completes that CLI integration. Geometry is still
classified once per NEC job, while complex material constants, interface
charge factors, and the complete table bank are refreshed for each frequency
before matrix assembly. `--ground-engine nec4` is no longer reference-only.
It attaches the immutable M9 matrix context to the existing NEC4 wire
operator; `--nec4-ground-model grid` is the default and
`--nec4-ground-model direct` retains the contour oracle for comparison.
`--nec4-ground-cache DIR` selects a cache, `none` disables it, and the default
is `$XDG_CACHE_HOME/nec2c/nec4-ground` or the equivalent user cache directory.

The production witness uses two horizontal wires in opposite media so that
their mutual terms traverse the model region rather than only the near tables.
The level-3 result differs from direct contour by `1.123e-8` in feed impedance
and `1.417e-8` in the L2 norm of all complex wire currents. On this deliberately
small witness, cached model-region lookup reduced matrix fill from about
`1.04 s` to `0.002 s` (roughly `500x` for that phase); this is a routing benchmark,
not an end-to-end antenna speed claim. A cold run generates four tables, an
exact repeat loads all four and preserves all reported numerical results
bit-for-bit, and direct mode reuses the same three near-table cache entries.

The runtime records whether any asymptotic interaction enters the documented
`<=15 degree` grazing-risk region. The CLI emits one deduplicated warning and
machine diagnostic after matrix fill; short segments alone do not trigger it.
The `z < 0.16 lambda0` direct strip is likewise reported in the human output.
No NEC2 interface correction is silently substituted by the M9 matrix path.

The manual stress witness is a 5 m whip connected at the interface to a 1 m
stake in seawater (`epsilon_r=80`, `sigma=4 S/m`) at 14 MHz. Its cold four-table
build took `170.5 s` on 16 workers and its warm load took `0.127 s`, with exact
machine-result identity. The accelerated solve gave
`14.5853-j2554.2219 ohm`, versus `14.5822-j2554.2228 ohm` in direct mode:
`1.257e-6` relative impedance difference and `6.794e-6` L2 current difference.
Matrix fill was `1.893 s` accelerated versus `6.740 s` direct (`3.56x`). The
smaller speedup than the separated-wire witness is expected and desirable:
this antenna exercises the correctness-first direct interface strip heavily.

## Measurement contract

Reports use both relative error and a scale-aware absolute error for field
nulls.  Initial engineering targets are:

- near-table interpolation: 5e-3 null-aware relative across the frozen
  production layout, while centered convergence diagnostics retain tighter
  `1e-4`-class targets where supported;
- model-based interpolation: 5e-3 relative, matching the guide's stated
  several-percent historical envelope only where the documented basis cannot
  meet the tighter target;
- asymptotic path: 1e-3 beyond its accepted switch distance; and
- switch jumps: no larger than the worse accepted error of the adjacent paths.

These are starting gates, not claims about achieved accuracy.  The limits will
be tightened when M8 sweeps show that the documented method supports it.

## User acceptance

The final M9 handoff will include:

- error and switch-continuity figures over the full M8 qualification grid;
- cache fresh/load/corruption tests;
- a 5 m seawater whip/ground-stake model; and
- a near-interface field sweep that crosses every production method boundary.
