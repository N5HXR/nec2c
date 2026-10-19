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
`MAX_CELLS` failure, overlap and cross-material-contact rejection, aspect-ratio
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
