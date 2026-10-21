# Legacy numerical regression suite

This suite freezes representative NEC2C behavior before dielectric-volume
support is introduced. It intentionally contains no `DF` cards.

## Run it

From a configured source tree, one command builds NEC2C and runs the suite:

```sh
make regression
```

The same test is also part of `make check`. A clean-tree verification is:

```sh
./autogen.sh
make check
```

Python 3.6 or newer is required only for the test runner, not for NEC2C.

`make check` also runs the M2 dielectric-resource parser, M3 deterministic
mesher, M4 cuboid-kernel, M5 dielectric-solver, and fast M6 coupled-wire
acceptance suites. The
parser's canonical fixture and generated
negative cases cover `DF` quoting and path resolution, every version-1 record,
deterministic diagnostics, `NX` cleanup, long-card rejection, and the explicit
execution gate. See
`docs/DIELECTRIC-RESOURCE-FORMAT.md` for the producer-facing contract.

## Coverage and provenance

| Case | Provenance | Legacy paths covered |
| --- | --- | --- |
| `EX1.nec` | Pre-existing NEC example 1 in `Input/` | wire, voltage source, complex current and impedance, power, near electric field |
| `EX4.nec` | Pre-existing NEC example 4 in `Input/` | mixed wire/patch geometry, perfect ground, radiation pattern |
| `EX5.nec` | Pre-existing NEC example 5 in `Input/` | crossed transmission-line (`TL`) network, radiation pattern |
| `EX7.nec` | Pre-existing NEC example 7 in `Input/` | plane-wave excitation and receiving patterns |
| `patch-only.nec` | Added for M1; minimal rectangular patch mesh derived from the `SM`/`SC` syntax already used by repository examples | patch-only solve, plane wave, patch current, radiation pattern |
| `nt-network.nec` | Added for M1; minimal two-wire model | reciprocal admittance-matrix (`NT`) network, driven input, radiation pattern |

The baseline was recorded on 2026-10-18 from the unmodified NEC2C solver at
the start of dielectric milestone M1, using GCC on x86-64 Linux. The runner
uses temporary output files and leaves the source tree unchanged.

## Numerical contract

`baseline.json` stores values exactly as printed by NEC2C. The runner selects
complex current and impedance components, input power, network parameters,
near-field components, surface current components, and radiation-pattern
samples. It compares each scalar with:

```text
abs(actual - expected) <= max(5e-8, 5e-4 * max(abs(actual), abs(expected)))
```

The relative tolerance is slightly wider than NEC2C's four-to-five-digit
report precision. The absolute floor protects printed zeros and very small
cross-polarized components. A change outside this envelope needs review; do
not automatically replace the baseline.

Every normal run also perturbs one expected current by one percent and
requires the comparator to reject it. This is the suite's expected-failure
control and proves that the comparison path is active.

## Review and regenerate a baseline

First capture proposed values without overwriting the reviewed fixture:

```sh
python3 tests/regression.py --generate > /tmp/nec2c-baseline.proposed.json
diff -u tests/baseline.json /tmp/nec2c-baseline.proposed.json
```

`--generate` prints only the `cases` object, so preserve the format and
tolerance fields when accepting reviewed values. For every accepted change,
record the solver change that caused it, inspect the corresponding full NEC2C
report, and run `make regression` after editing `baseline.json`.

## Baseline report

M1 baseline result on 2026-10-18: six decks passed, covering eleven named feature
classes and fifteen selected numerical samples. The one-percent negative
control was rejected as expected. Both the normal in-tree build and a clean
build from the generated distribution archive passed.

## M2 parser report

M2 acceptance on 2026-10-18 passed the canonical-summary fixture and twenty-three
checks covering valid records, forward references, canonical ordering,
malformed and missing records, duplicate IDs and records, non-finite values,
unsupported versions and keywords, bad references and bounds, second and
overlong `DF` cards, quoted paths with spaces, deck-relative and absolute path
behavior, `NX` cleanup, and the parser-only solver gate. The M1 numerical
suite remained unchanged. The generated distribution archive built and passed
both suites from a clean directory.

## M3 mesher report

M3 acceptance on 2026-10-18 passed the canonical mesh/cost fixture and fourteen
checks covering canonical ordering, linear and multiplicative highest-frequency
selection, both `MAX_KH` boundaries and rejection above one, pre-allocation
`MAX_CELLS` failure, overlap rejection and cross-material contact acceptance, aspect-ratio
refinement, volume conservation, vacuum omission, wire and patch clearance,
maximum emitted-cell electrical size, and cost fields. Per-region physical and
electrical diameters are frozen in the canonical report. The unchanged M1 and
M2 suites also passed, including from a clean build of the generated source
distribution archive.

## M4 cuboid-kernel report

M4 acceptance on 2026-10-18 passed all nine Cartesian components for dynamic
self, face-, edge-, and corner-adjacent cuboids and unequal near and far
cuboids against the checked high-order reference tensors. Additional tests
cover the analytical static cube value `diag(2/3, 2/3, 2/3)`, the static
rectangular-cuboid trace identity, translation and Cartesian permutation,
volume-scaled reciprocity, electrical-scale invariance, direct `k0=0`
evaluation, near-contact continuity at two positive gaps, and deterministic
failure outside the certified electrical-size and geometry domain.

The reference generator compares Gauss orders 36 and 52 before emitting data.
Its singular surface integrals were separately checked against DIRECTFN for
same-surface, common-edge, and common-vertex quadrilaterals. See
`docs/VIE-CUBOID-KERNEL.md` for equations, units, tolerances, provenance, and
the recorded comparison values. DIRECTFN and NumPy are validation-time tools;
the NEC2C build has no new runtime dependency.

## M5 dielectric-only solver report

M5 acceptance on 2026-10-18 exercises the public NEC deck path for a vacuum
identity, an analytical small cube, a 32-cell staircase sphere, axis
permutation, complex residuals, lossless power closure, dipole radiation, and
`RP` far-field output.  The dynamic sphere is checked against an independent
Lorenz-Mie calculator at requested-radius size parameter 0.8 and against a
pinned 32-dipole VoxScatter result generated by
`tools/voxscatter_m5_reference.m`.  See
`docs/DIELECTRIC-SOLVER-M5.md` for equations, tolerances, provenance, and the
intentional dense/free-space scope. The complete five-suite test set also
passes from a clean build of the generated source distribution.

## M6 wire/dielectric coupling implementation report

The M6 implementation on 2026-10-19 covers the complete-vector vacuum identity,
combined-system residual, decreasing effect with separation, dielectric and
wire refinement, reciprocal transfer admittance in an asymmetric two-port
model, and explicit rejection of unsupported pattern output. The canonical
64-cell N2 dipole/cube comparison is retained as a slower manual fixture. It
passes the 300 MHz and resonance-shift comparisons, but its 360 MHz reactance
shift differs from openEMS by about 1.42 ohm, within the accepted 2-ohm
cross-solver tolerance. M6 is therefore `COMPLETE`. See
`docs/DIELECTRIC-WIRE-M6.md` and
`docs/M6-OPENEMS-REPORT.md` for the formulation, scope, raw-data links, and
performance note.

## M7 conducting-patch/dielectric coupling implementation report

M7 adds patch-only and mixed wire/patch coupling through NEC's existing MFIE
surface basis. The fast suite covers both vacuum identities, all four
conductor/volume block directions, residuals, separation, surface refinement,
volume refinement, and the retained 10 mm openEMS impedance-perturbation
comparison. The openEMS comparison is within the accepted 2-ohm component-wise
tolerance; its opposite-sign sub-ohm reactance shift remains documented as a
refinement caveat. See `docs/DIELECTRIC-PATCH-M7.md`.

## M8 analytical field increment

`make check` includes `tests/dielectric_fields.py`. It needs only Python's
standard library. It checks the finite-cube far-field form factor against the
reported dipole moment, analytical electric-dipole E/H fields at 1, 10, and
1000 m, near-to-far consistency, vacuum zero fields, observation-clearance
rejection, deferred frequency sweeps, legacy wire/patch/mixed vacuum field
identity, and driven accepted/radiated power. The M6 RP rejection test now
checks successful combined-field output. See the M8 progress section in
`docs/DIELECTRIC-VOLUME-DESIGN.md` for conventions, tolerances, limitations,
and remaining release gates. This increment does not complete M8.

Validation on 2026-10-19: all eight `make check` programs passed, including
the complete M1–M7 suites and the new analytical field suite. `make dist`
succeeded, and the archive includes the M8 tests and design/progress document.

## M9 iterative-solver foundation

`make check` includes `tests/vie_iterative_test` and
`tests/dielectric_iterative.py`. The C test checks a matrix-free complex
operator, preconditioned GMRES, restarts, known solutions, residuals, and
explicit failure paths. The Python test compares dense and experimental
GMRES currents, feed impedances, complex patterns, and available power
reports on volume-only, wire, patch, and mixed fixtures. It also validates
the tree-sizing inventory on known full/sparse grids and invalid geometry.

Use `./nec2c -s gmres -i model.nec -o model.out` to opt into the experimental
DF solver. It still uses dense storage and diagonal Jacobi preconditioning;
it is not yet an accelerated tree solver. The default remains `-s dense`.
`tools/m9_tree_sizing.py` regenerates `tests/m9-tree-sizing.csv` without an
external EM solver. See the M9 progress section in the design document for
its assumptions, estimates, solver-equivalence tolerances, and remaining gates.

M9 first-increment validation on 2026-10-19: all ten `make check`
programs passed with no failures or skips. The two affected Python suites
were rerun successfully after using Python-3.6-compatible subprocess calls.
The new iterative C module also passed strict C99 compiler checks with
`-Wall -Wextra -Werror -pedantic` and AddressSanitizer/UndefinedBehaviorSanitizer.
`make dist` succeeded and the archive includes the implementation, tests,
sizing tool, and retained CSV. M9 remains in progress.

### M9 experimental FFT operator

`make tests/vie_fft_test tools/m9_fft_benchmark` builds the isolated uniform-
cube volume backend. Run `./tests/vie_fft_test` for certified-dense product
and known-solution checks, including holes, boundary impulses, repeat calls,
static and dynamic fields, the certified electrical cell-size limit, and
invalid/memory-budget inputs. It is also included in `make check`.

`./tools/m9_fft_benchmark skeleton 1e-8 150` and
`./tools/m9_fft_benchmark tree 1e-8 150` measure the 216- and 9,432-cell
free-space cases under a 512 MiB allocation cap. The second numeric argument
is the maximum Arnoldi iteration count. A failed solve exits nonzero.
The benchmark reports normalized complex currents/pattern sums and independently
refined extinction/scattering power checks; see
[the experiment report](../docs/M9-FFT-EXPERIMENT.md) for their units and limits.
`nec2c -s gmres` retains dense storage; `nec2c -s fft` now selects the
production uniform-cube FFT integration. See below.

### M9 integrated FFT solver

`./nec2c -s fft -r 1e-10 -n 300 -i model.nec -o model.out` selects the
accelerated volume block with validated NEC conductor/cross blocks. It
requires equal cubic cells on a common lattice and one active permittivity,
with a 512 MiB accounted solver cap and no silent dense fallback.

`tests/dielectric_fft_test` verifies hybrid block layout, known-solution
accuracy, preconditioning, empty volumes, invalid lattices/materials and
budget rejection. `tests/dielectric_iterative.py` now compares dense, dense
GMRES and FFT solutions for volume, wire, patch and mixed eight-cell fixtures,
and checks no-DF identity, CLI validation and explicit nonconvergence.

`python3 tools/m9_tolerance_sweep.py --output /tmp/m9-tolerance.csv` measures
feed-impedance and complex-pattern error against dense references as residual
tolerance is tightened. `python3 tools/m9_nec_benchmark.py --model tree
--output-dir /tmp/m9-n3` generates the 9,432-cell production antenna benchmark
and retains inputs, output, elapsed time and peak RSS. See
[the integration note](../docs/M9-FFT-INTEGRATION.md) for scope and evidence.

## M10 perfect-ground kernel

`make check` includes `tests/vie_ground_test`, the isolated positive-gap PEC
volume operator tests. They use independent finite-volume dipole quadrature,
reciprocity, boundary and remote-ground limits, self-image checks and invalid
input rejection. This does not enable dielectric GN decks; see
[the formulation record](../docs/M10-GROUND-FORMULATION.md) for integration
requirements and the real-ground compatibility decision.

## M10 PEC deck integration

`tests/dielectric_ground.py` exercises positive-gap `GN 1` decks on dense and
GMRES backends. It compares currents and complex fields against explicit
mirrored free-space conductor/volume models, checks vacuum compatibility,
PEC boundary conditions, power and frequency-sweep state restoration, and
rejects unsupported ground/backend/clearance requests.

A runnable driven example is:

```sh
./nec2c -s gmres -r 1e-10 -i tests/dielectric/m10-pec-wire.nec -o /tmp/m10-pec.out
```

Only planar perfect ground is enabled in this increment. Real-ground modes,
FFT images and contact remain outside the implemented subset. See the
[M10 formulation record](../docs/M10-GROUND-FORMULATION.md) for its exact
scope and remaining milestone gates.

## M10 finite-ground reference and battery-limited OpenEMS runner

`tests/vie_halfspace_test` checks the isolated reflected E/H spectral kernel
against Fresnel/quasistatic/PEC limits, an independent finite-loss axial
reference, finite-volume reciprocity and explicit failure cases. Reproduce
its scalar oracle with `python3 tools/m10_halfspace_reference.py`. Mixed
point/cuboid checks cover source-volume normalization, lossy electric
reciprocity and independent real-space PEC electric/magnetic image integration.
GN 0/2 deck integration remains disabled.

`tests/m10_openems_runner.py` validates the portable runner's default plan,
argument guards and termination of a worker plus simulated child solver. It
requires no OpenEMS installation. Follow
[M10-OPENEMS-LAPTOP.md](../docs/M10-OPENEMS-LAPTOP.md) to stage bounded jobs
through SSH to the user's laptop. Pilot measurements are cost diagnostics;
verification output requires decay, mesh/domain convergence and physical
comparison review before acceptance.

`tools/m10_review_results.py` independently recomputes port impedances using
each trace's own timestamps, checks completed native JSON results, and labels
incomplete reconstructed pairs provisional. Returned normalized evidence is
in `m10-openems-ground-pairs.csv` and `m10-openems-ground-review.json`; the
single-frequency NEC PEC checkpoint is in `m10-nec-pec-300.json`. Neither
successful simulation completion nor a stable perturbation alone demonstrates
absolute feed/source/domain convergence. See the laptop report for review.

`tests/dielectric_real_ground.py` covers the GN 2 wire/volume integration:
vacuum compatibility with native NEC, loaded dense/GMRES agreement, linear
and elliptical plane-wave extinction/air/substrate power, a matched half-space,
lossless transmission, conductivity and negative-sigma frequency sweeps, and
Maxwell's curl relation between independently requested electric/magnetic
near fields. `tests/dielectric/m10-real-wire.nec` is a small runnable fixture.
The half-space C test additionally checks interface flux against independent
local dipole work and upper-hemisphere radiation; it explicitly rejects an
unresolved narrow grazing transmission limit.

`tests/m10-nec-real-300.json` records the completed 64-cell GN 2 wire/volume
checkpoint at 300 MHz and its comparison with the completed coarse OpenEMS
pair: complex dielectric-impedance-perturbation difference 0.27573 ohm and
NEC power defect 0.0876%. The metadata explicitly excludes full M10 acceptance.

Milestone boundary updated 2026-10-20: M10 accepts PEC-only ground and M11
continues lossless PEC-ground contact. The GN 2 kernels, integration tests and
real-ground checkpoint data above are retained as preliminary M13 work. Their
M10-prefixed filenames record provenance; further real-ground verification
does not block M10/M11. M12 conductor/dielectric contact and penetration are
deferred to a later revision; M13 initially retains positive conductor/
dielectric clearance. See [the handoff](../docs/M13-LOSS-AND-REAL-GROUND.md).

`m10_pec_validation.py` checks PEC port reciprocity with unequal wires,
the driven remote-ground limit, and plane-wave extinction versus scattered
power for wire, closed-patch, and mixed conductors. Physical patch power
checks use closed surfaces with outward normals. `m10-pec-acceptance.json`
retains the sampled-frequency impedance comparison and refinement evidence;
its status explicitly lists remaining acceptance work. Reproduce the larger
VM checkpoints with `tools/m10_pec_acceptance.py --group volume|frequency|wire
--workdir /tmp/your-checkpoints --binary ./nec2c`; choose one group per command.
Each work directory records its binary hash to prevent cross-build reuse.

M10 PEC acceptance completed on 2026-10-20. The returned four-job feed/domain
audit is retained in `m10-pec-sensitivity-review.json`, including archive hash,
normal -60 dB stopping results, mesh metrics, and sensitivity comparisons.
`m10-pec-sensitivity-pairs.csv` contains the six frequency/variant rows.
Domain/feed changes in delta Z are below 0.0063 ohm. These are sampled
sensitivity results; the 2 ohm acceptance criterion concerns dielectric-induced
impedance changes, not absolute NEC/OpenEMS feed impedance.


M11 contact tests (`m11_pec_contact.py`) compare exact contact with explicit
free-space images for all three plane-wave polarizations, verify solved-current
and feed continuity through narrow positive gaps, and compare dense/FFT
solutions and power for wire, closed-patch and mixed conductors. GN 2 contact
and buried dielectric remain rejected. `vie_ground_test` now also checks
`m11-static-reference.csv`, generated independently by the closed-form
Coulomb-face oracle in `tools/m11_static_reference.py`.

The PEC FFT adapter embeds mirrored source currents with reversed tangential
components and retains physical unknowns only. `dielectric_fft_test` compares
that operator against dense PEC matrices at contact and a positive gap and
checks its lattice restrictions. The free-space preconditioner is approximate;
acceptance uses the true residual of the PEC system.

Returned contact OpenEMS evidence is retained in `m11-openems-contact-review.json`
and `m11-openems-contact-pairs.csv`. Decay and the 25 mm NEC physical comparison
pass. The documented gate is component-wise delta R/delta X within the greater
of 2% or 2 ohm; complex-vector discrepancy is an additional diagnostic. Native
refinement and full-tree acceptance pass. M11 is complete (2026-10-21): the
34.60-minute grounded-tree run passes the user-approved 40-minute target, with
234.8 MiB peak RSS. See `m11-acceptance.json` and `m11-pec-tree.json`.

### M13-01 scalar dielectric loss

`m13_loss.py` verifies passive scalar complex permittivity and conductivity,
analytic small-cube polarization and Joule absorption, zero-loss compatibility,
frequency sweeps, conductive ER=1, invalid input rejection, and lossy
dense/GMRES/FFT current and energy agreement in free space and over PEC.
`tools/run-m13-01-loss.sh` stages independent laptop OpenEMS comparisons;
returned decay, mesh/domain refinement and native comparisons pass. M13-01
is complete (2026-10-21), with 20/20 regression groups passing. Retained evidence
is in `m13-01-acceptance.json`, `m13-01-openems-review.json` and
`m13-01-openems-pairs.csv`. See `docs/M13-01-SCALAR-LOSS.md`.

### M13-02 frequency tables

`m13_tables.py` checks interpolation against constant-material references and
an independent small-cube solution, EI/SIGMA distinction, absorption, solver
agreement, sweep order/geometric steps, whole-sweep extrapolation rejection,
interior material peaks in meshing, canonical table validation, coupled power
and NX lifecycle. Independent witnesses are staged with
`tools/run-m13-02-tables.sh`; see `docs/M13-02-MATERIAL-TABLES.md`.

M13-02 is complete: native analytical/refinement checks, all 21 regression
groups and independent OpenEMS witness acceptance pass. All 12 verification
runs reach -60 dB decay; maximum component error is 1.33655 ohm. Evidence is
retained in `m13-02-acceptance.json`, `m13-02-openems-review.json` and
`m13-02-openems-pairs.csv`. Only designated witness rows compare to the table;
the 300 MHz knot reuses the accepted M13-01 reference.

### M13-03 touching materials

`m13_materials.py` checks exact electrostatic contact blocks, local absorption,
material identities/tables, topology and overlap rules, solver equivalence,
coupled wire power and an independent coated-sphere polarizability/absorption
reference on three voxel refinements. The FFT C tests include heterogeneous
lossy free-space/PEC products and material-specific preconditioning.
Native antenna counterparts are in `tools/m13_03_acceptance.py`; independent
OpenEMS jobs are staged with `tools/run-m13-03-materials.sh`. See
`docs/M13-03-TOUCHING-MATERIALS.md` and `m13-03-acceptance.json` for scope and evidence.
M13-03 is complete: all 22 groups, analytical/refinement checks and independent
OpenEMS acceptance pass. All 11 native comparisons pass; largest component
error is 1.97370 ohm, with 1.69712 ohm in the joint refined hardest case.
Returned pair data and audits are retained in `m13-03-openems-pairs.csv` and
`m13-03-openems-review.json`.

### M13-04 GN 2 combined loss

`m13_real_ground.py` checks an independent quasistatic finite-cube image and
Fresnel drive, local absorption, ground/material limits, touching frequency tables,
sweep state, dense/GMRES power and reciprocity, and explicit unsupported-mode
refusals. Native antenna/refinement cases use `tools/m13_04_acceptance.py`;
OpenEMS lossless/lossy controls use `tools/run-m13-04-real-ground.sh` and are
reviewed separately by `tools/m13_04_review.py`. Status and scope are in
`docs/M13-04-REAL-GROUND-WIRES.md` and `m13-04-acceptance.json`.

M13-04 returned GN 2 controls pass all 12 verification-run physical audits and
six native comparisons (worst component error 0.829761 ohm). Retained evidence:
`m13-04-openems-review.json`, `m13-04-openems-pairs.csv` and
`m13-04-acceptance.json`. `tools/m13_04_audit.py` reproduces decay, material,
geometry, grid quality, refinement and native component checks from the archive
and native benchmark directories; see the M13-04 document for its command.

M13-05 completes the agreed M13 scope with production-tree material provenance
and sampled-sensitivity reporting. `m13_tree_report.py` tests finite/passive input,
convergence refusal, common-mesh material scenarios, separate nominal refinement
and retained-report reproduction. `vie_halfspace_test` includes a localized-source
remote-ground E/H limit. All 24 groups pass. Reproduce the synthetic native tree
witness using `tools/m13_05_tree_experiment.py`; report assumptions/ranges using
`tools/m13_05_tree_report.py`. See `docs/M13-05-TREE-MATERIAL-GUIDANCE.md` and
`m13-05-acceptance.json` for the exit-criterion mapping and explicit exclusions.

Revision-1 release entry points are `docs/DIELECTRIC-REVISION-1.md` and
`docs/DIELECTRIC-DEFERRED-WORK.md`. The canonical examples in
`examples/dielectric-r1/` include a standalone runner and numerical reference
summary. They supplement the 24 regression groups with an executable user
workflow; `make distcheck` checks the complete source distribution.
