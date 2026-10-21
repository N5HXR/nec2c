# M9 FFT integration and backend decision

> Historical milestone record. The [revision-1 guide](DIELECTRIC-REVISION-1.md)
> states current supported scope; later milestones extend some restrictions below.
> In particular, M13-03 adds cell-local heterogeneous scalar materials on the
> common FFT lattice; the original uniform-material restriction is superseded.


The formulation is a hybrid wire/surface–volume integral-equation extension
of NEC2: thin-wire and conducting-patch equations are coupled to dielectric
polarization-current volume equations. “VSIE extension” is an appropriate
short description, with the wire representation stated when precision matters.
FFT acceleration changes the volume operator's application, not the physical
formulation.

## Supported production path

`nec2c -s fft` uses FFT-backed GMRES for DF models with one uniform scalar permittivity
and equal cubic cells on a common lattice. The M9 release supported free
space; the M11 extension also supports PEC `GN 1` through a reflected-source
embedding, including dielectric contact. PEC bottom clearance must be on the
half-cell lattice; real-ground FFT remains unsupported. M11 physical/tree
acceptance is complete (2026-10-21), with a separately agreed 40-minute
grounded-tree target; see [M11-PEC-CONTACT.md](M11-PEC-CONTACT.md). The dense solver
remains the default; `-s gmres` remains the dense-reference iterative adapter.
The M13-01 increment extends the same operator to a passive complex scalar
permittivity, with independent lossy acceptance completed in
[M13-01-SCALAR-LOSS.md](M13-01-SCALAR-LOSS.md).
A no-DF deck retains its legacy solver for all three selections.

```
./nec2c -s fft -r 1e-10 -n 300 -i model.nec -o model.out
```

`-r` sets the true relative GMRES residual target (default `1e-12`), and `-n`
sets its maximum Arnoldi iterations (default 300, accepted range 1..10000).
GMRES restart remains 30. Nonconvergence, invalid geometry and budget failures
return a failed job before solved currents/patterns are printed. There is no
silent dense fallback.

The adapter verifies every cell's dimensions, material and lattice coordinates
before building the FFT plan. It rejects unequal cuboids, differing active
permittivities, nonlattice centers, duplicate sites and oversized embeddings.
Existing per-box meshing is unchanged. Select `MAX_KH` so the intended box
lengths divide into the same cube pitch; it is an electrical size bound, not
a request to snap or alter box boundaries. For the 28 MHz, epsilon_r=3 tree,
`MAX_KH=0.2541077043967435`, `MAX_ASPECT=1` and `MAX_CELLS=10000` produce
9,432 cubes of pitch 0.25 m. Other frequencies/geometries need revalidation.

## Joint operator

The system is stored as

```
[ C   B ] [ conductor coefficients ]
[ D   V ] [ volume currents        ]
```

`C` retains the validated dense NEC conductor matrix, including mixed
wire/patch rows. `B` and `D` retain the existing conductor/volume quadrature
and scaling. `V` is applied by FFT; it is never assembled into a dense matrix
on this path. Cross-block storage is `2*n_conductor*n_volume` complex values,
so many conductors can still be expensive. This is not an acceleration of
NEC's conductor matrix or of coupling quadrature construction.

The preconditioner uses the complex diagonals of `C` and the volume self
block. Both plane-wave and voltage-source right-hand sides, solved-current
conversion, combined patterns, exterior fields and power diagnostics retain
the existing formulas. The true unpreconditioned joint residual certifies
GMRES acceptance; observable comparisons are also required because a small
residual alone does not guarantee a small impedance error.

The 512 MiB cap includes requested plan/coupling allocations, conservatively
reserved GMRES workspace, caller solver arrays, temporary lattice map and
cell geometry. The printed `FFT ACCOUNTED SOLVER BYTES` is this conservative
inventory; allocator overhead and unrelated legacy/global process state are
not a total process-memory guarantee. Peak RSS is measured separately.

Exactly zero wire-current basis coefficients are skipped in near-electric
field evaluation. This removes unused kernel evaluations during unit-source
coupling assembly without approximating the remaining interactions.

## Backend decision for the selected tree

Choose regular-grid FFT for the initial 9,432-cell, uniform-cube production
scope. The filled canopy gives 34.1% occupancy of a 24x24x48 physical grid;
padding to 64x64x128 costs approximately 101 MB for twelve complex arrays.
The volume experiment already converged using about 127 MiB peak RSS,
compared with the 25.6 GB sizing estimate for two dense matrices.

A sparse embedding retains only active unknowns, but a conventional FFT still
covers the padded box. The skeleton has only 3.99% occupancy and therefore
wastes proportionally more FFT work. A sparse FFT is not implemented, and
zero-padding does not confer sparse-FFT complexity.

An occupied octree with capacity 32 has 1,345 nodes and 1,172 leaves for the
small full tree (see `tests/m9-tree-sizing.csv`). An H-matrix on this tree
could compress well-separated interaction blocks and avoid empty-grid work,
but these node counts do not establish ranks, accuracy or achieved memory.
It would introduce additional approximation/error controls and implementation
work. Since the selected FFT case fits the memory cap with substantial room,
there is no present resource requirement that justifies that added complexity.
A hierarchical backend remains a candidate for sparse, nonuniform or larger
models; it is not claimed inferior on unmeasured numerical ranks/performance.

The 0.125 m full tree requires about 1.03 GB for the current FFT arrays plus
GMRES alone and is outside this initial 512 MiB scope. Do not extrapolate the
small-tree result to arbitrary volume size or conductor count.

## Reproduction and evidence

```
make check
python3 tools/m9_tolerance_sweep.py --output tests/m9-tolerance-sweep.csv
python3 tools/m9_nec_benchmark.py --model tree --output-dir /tmp/m9-n3
```

The integration suite requests `rtol=1e-13` for strict per-component overlap
comparisons (relative `1e-8`, absolute `1e-12`); this is stricter than the
normal default and retains the comparison thresholds. It compares
volume-only, wire, patch and mixed fixtures
against dense solutions on the identical mesh. Eight-cell fixtures exercise
cross-cell convolution as well as conductor coupling. It checks currents,
complex principal cuts, available power diagnostics, no-DF legacy identity,
invalid solver controls and explicit nonconvergence. A separate hybrid unit
test checks both cross-block layouts, preconditioning, a known solution,
empty dielectric, lattice/material rejection and memory limits.

The tolerance tool compares an 11-segment wire and a mixed wire/patch fixture
at targets `1e-4`, `1e-6`, `1e-8`, `1e-10` against dense solutions. It retains
feed impedance error in ohms and complex pattern error scaled to the largest
reference field, including samples near pattern nulls. At `rtol<=1e-10` it
requires impedance error below `1e-4` ohm and scaled pattern error below
`1e-6`, much smaller than the separate 2 ohm physical accuracy allowance.

The N3 script generates the agreed tree plus an 11-segment, 5 m dipole at
x=-5 m, y=0, z=3.5..8.5 m, radius 2 mm, driven by 1 V at its center at 28 MHz.
The antenna is separated from the dielectric and is not grounded. Inputs,
full output, exit status, elapsed wall time and Linux child peak RSS are
retained. This tests the production path and coupling cost; it does not
certify antenna discretization or realistic vegetation material assumptions.

The prior volume-only measurements and analytical/power checks remain in
[M9-FFT-EXPERIMENT.md](M9-FFT-EXPERIMENT.md). OpenEMS is not required for these
operator-equivalence gates. M8 independent physical validation and later
milestones remain separate.

The initial wire sweep met those observable goals at `1e-8`, but the mixed
fixture at the same residual target had 0.000162 ohm impedance error and
`2.68e-6` scaled pattern error. The tighter `1e-10` point is therefore the
acceptance checkpoint; the failed looser point is retained as evidence of
fixture-dependent observable sensitivity, not removed or treated as a
validated observable at that looser target.
The default remains `1e-12`; requested residual alone is not a universal
observable-error certificate.

## Measured tolerance results

Both fixtures use eight dielectric cells. Errors are against dense solutions
on exactly those meshes, not against a physically exact antenna model.

| Fixture | Requested residual | True residual | Iterations | Impedance error (ohm) | Scaled complex-pattern error |
| --- | ---: | ---: | ---: | ---: | ---: |
| m6-cube.nec | 1e-04 | 2.73e-05 | 7 | 4.2e-05 | 6.04e-07 |
| m6-cube.nec | 1e-06 | 1.92e-07 | 9 | 1.28e-07 | 4.55e-09 |
| m6-cube.nec | 1e-08 | 3.61e-09 | 13 | 8.42e-10 | 1.72e-11 |
| m6-cube.nec | 1e-10 | 5.03e-11 | 15 | 0 | 1.84e-13 |
| m7-mixed-cube.nec | 1e-04 | 7.9e-05 | 6 | 0.425 | 0.00703 |
| m7-mixed-cube.nec | 1e-06 | 4.76e-07 | 8 | 0.00491 | 8.11e-05 |
| m7-mixed-cube.nec | 1e-08 | 7.23e-09 | 12 | 0.000162 | 2.68e-06 |
| m7-mixed-cube.nec | 1e-10 | 4.89e-11 | 15 | 3.63e-07 | 5.99e-09 |

At `1e-10`, the mixed fixture meets the retained `1e-4` ohm / `1e-6` pattern
comparison goals with substantial margin. The original `1e-8` checkpoint was
insufficient for that stricter pattern goal, although its impedance error
was already much smaller than the separate 2 ohm physical allowance. The
accuracy goals were retained and the residual target tightened. Raw values
and iteration counts are in `tests/m9-tolerance-sweep.csv`.

For a fixed discretized operator, `A*x=b` and computed `x_hat`, the residual
`r=b-A*x_hat` gives `x_hat-x=-A^{-1}r`. A sampled complex field is a linear
functional of the currents, whereas feed impedance `V/I` is sensitive to
current error when `I` is small. This is why the true residual and observable
errors are checked separately; no condition number was measured here.

## Integrated N3 measurement

On 2026-10-19 the production path completed the 9,432-cell tree plus the
11-segment dipole at `rtol=1e-10`: 45 Arnoldi iterations, true relative
residual `7.7532e-11`, peak RSS 140,064 KiB (136.8 MiB), and accounted solver
storage 141,975,080 bytes (135.4 MiB), below the 512 MiB cap. It returned exit
status zero. Computed feed impedance was `67.27269449915 - j34.57141258790`
ohm. These are benchmark results, not independently certified physical
impedance values for the idealized tree.

Total wall time was 1,224.76 s (20 min 25 s); solver CPU time was 86.46 s.
Other tests/development jobs ran concurrently. The total includes mesh and
kernel construction, conductor/cross-block assembly, iteration and output/
power integration. Separate phase times were not instrumented in this first
integrated run; the isolated volume phase measurements are in the prior
experiment note. Conductor/cross construction is additional work that the
volume FFT does not accelerate.

The lossless power defect was `2.8175e-3` (0.282%), comparable to the inherited
conductor/volume power-discretization defect; it is not the machine-precision
volume-only power check and does not close M8's separate validation gates.

Inputs are retained as `tests/dielectric/m9-tree.nec` and `.diel`; run this
large fixture with `-s fft`. Metadata and complex principal cuts are in
`tests/m9-n3-fft.json` and `tests/m9-n3-pattern.csv`. The full first-run output
is in `/tmp/m9-n3-integrated/model.out` on the development machine; the
benchmark tool regenerates and retains complete outputs in its output folder.
The strict multi-cell comparisons passed after the contact correction below
and tightening the comparison residual target. Runtime-budget agreement
remains open at this checkpoint.

## Contact-quadrature correction discovered by integration

The eight-cell translated patch fixture initially missed the strict current
comparison. An isolated touching-cell pair differed from its canonical
translation by `5.69e-8` under the dense kernel's `5e-7 / 5e-10` quadrature
tolerances. Floating-point subtraction produced a machine-level gap while
the topology classifier identified contact. The integration code treated
that gap as nonzero, selecting a different near-singular map and order 48
instead of 8. This was numerical geometry sensitivity, not FFT indexing.

Face-integration breakpoints and fixed separations now resolve coincidence
using a local geometry scale and the existing 64-machine-epsilon tolerance
constant. Tiny duplicate/zero breakpoints are coalesced. The tolerance scales
with local face coordinates rather than a fixed physical gap, so meaningful
positive gaps are retained when electrically equivalent geometry is scaled.
The translated pair now differs by `1.11e-16`, using order 8 in both cases.

New kernel regressions cover all seven touching offsets, a one-ulp outward
perturbation and a resolved positive gap under `1e-10` electrical scaling.
The independent kernel reference suite, contact continuity and other physical
regressions remain required. This corrects the dense reference as well as
FFT kernel generation; the comparison thresholds were not weakened.

For the N3 benchmark, quarter-metre cell displacements and half-widths are
exact binary fractions. All 13 canonical self/near blocks were audited
against the previous kernel and are exactly unchanged. Adaptive far volume
quadrature is unchanged. Thus the integrated N3 measurement remains
applicable. The audit is retained in `tests/m9-contact-kernel-audit.txt`.

The tolerance CSV was regenerated after this correction; the pre-correction
sweep is retained separately as `tests/m9-tolerance-before-contact-fix.csv`.
A zero impedance error in the CSV means agreement at the printed reference
precision, not proof of mathematically exact currents.

## Validation completion and agreed budget

All twelve regression programs pass, including the stricter multi-cell
comparison rerun at requested residual `1e-13`. The full kernel reference
suite, all touching-offset/one-ulp tests and the resolved-gap electrical-scale
test pass. Kernel, FFT and hybrid tests also passed AddressSanitizer and
UndefinedBehaviorSanitizer; LeakSanitizer was disabled after it failed in the
sandbox environment. Strict C99 warnings/syntax checks passed for the new
adapter/tests, and the source distribution build passed.

The backend decision, implementation, dense equivalence, observable tolerance
sweeps, explicit failure behavior and integrated N3 resource benchmark are
complete. On 2026-10-20 the user agreed to a 30-minute-per-frequency runtime
acceptance budget for the selected 9,432-cell, epsilon_r=3, 28 MHz tree and
11-segment antenna at `rtol=1e-10`, with the enforced initial 512 MiB accounted
solver cap. The measured run took 20 min 25 s and peak RSS was 136.8 MiB.
This acceptance budget is not an enforced wall-time timeout. Verification runs
may take as long as needed, as authorized by the user.

M9 is complete. A multi-frequency N3 benchmark would add coverage but is not
an M9 exit requirement. The resource measurement and acceptance apply to this
28 MHz fixture; they do not establish performance or convergence over a
frequency range. M8's separate physical-validation gates remain unchanged.
