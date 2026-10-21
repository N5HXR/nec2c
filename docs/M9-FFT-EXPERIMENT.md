# M9 uniform-cube FFT experiment

This report records the initial isolated volume experiment on the 2 GB machine.
The FFT module and standalone benchmark remain separate targets; the later
production integration is described in
[M9-FFT-INTEGRATION.md](M9-FFT-INTEGRATION.md). `nec2c -s gmres` still uses
dense storage, while `nec2c -s fft` now selects the integrated backend.
Ground, loss, arbitrary cuboids and multiple active materials remain outside
this FFT scope.

## Operator and storage

For constant real relative permittivity, the existing current equation is

```
A J = epsilon_r J - (epsilon_r - 1) N J.
```

Equal cubes on a common integer lattice make each tensor block depend only
on the target-minus-source displacement. The operator scatters the three
current components onto the embedding grid, applies a zero-padded linear
convolution, and gathers only occupied cells. Empty sites remain vacuum.
Each padded axis is a power of two at least `2*naxis-1`; no periodic physical
boundary is introduced. Six symmetric tensor spectra and six vector grids
use twelve complex arrays. A dependency-free radix-two FFT uses precomputed
roots and strided transforms; it requires no external FFT workspace.

Self and near interactions use the existing certified `vie_cuboid_n`
integration with the dense solver's tolerances (relative `5e-7`, absolute
`5e-10`). Cubic permutation/reflection symmetry reduces kernel evaluations;
off-diagonal components vanish when either associated displacement is zero.

For well-separated cubes only, smooth quadrature independently integrates
`(k^2 I + grad grad) exp(-ikR)/(4*pi*R)` over the difference of source and test
coordinates. Along each axis the difference measure is the triangular weight
`h-|u|` on `[-h,h]`. Split at zero, tensor Gauss rules of orders 4, 6, 8, and 12
are compared, accepting component changes below `2e-12 + 2e-9*|Npq|`.
Otherwise the certified face integrator is used with far order 16. This is
volume integration, not a point-dipole approximation. The certified electrical
cell-size restriction `k*h <= 1` remains enforced.

The API requires explicit integer coordinates and rejects duplicate or
out-of-range cells. It does not silently snap an existing mesh onto a lattice.
All owned allocations are overflow-checked and compared with the supplied
byte budget before allocation. A plan is mutable and must not be shared
between simultaneous operator applications. Its Jacobi preconditioner uses
the complex self block. GMRES recomputes the unpreconditioned true residual.

## Reproduction and interpretation

```
make tests/vie_fft_test tools/m9_fft_benchmark
./tests/vie_fft_test
./tools/m9_fft_benchmark skeleton 1e-8 150
./tools/m9_fft_benchmark tree 1e-8 150
./tools/m9_fft_benchmark tree 1e-10 150
```

The disjoint boxes match `tools/m9_tree_sizing.py`, with pitch 0.25 m,
frequency 28 MHz, and epsilon_r=3. The full tree contains 9,432 cubes; the
216-cell skeleton omits the canopy. These are suspended free-space objects.
Excitation is an x-polarized plane wave travelling along +z, with exact
cell-averaged incident phase. The solved variable is `J/(i*omega*epsilon0)`;
reported integrated currents and transverse far-field sums omit that common
factor. Their complex phase origin is the embedding grid's lower corner.
They are not antenna feed impedances or realized gains.

The benchmark caps the plan plus caller arrays and GMRES workspace at
512 MiB. It reports wall-clock build, single matvec and solve times, actual
peak RSS on Linux, iteration count, operator applications and true residual.
Angular power integration uses 16 and 32 Gauss-Legendre polar nodes, each
with twice as many uniform azimuth points. It compares their convergence
and lossless extinction versus scattering. The RSS/timing checkpoint is
before angular postprocessing, which uses constant-size local storage.
Nonconvergence returns a nonzero exit code and an explicit status.

## Validation

The operator test independently assembles dense blocks through the certified
cuboid integrator (far order 16). It covers a one-cell grid, an asymmetric
embedding with holes, every Cartesian basis impulse, general complex input,
repeated application, static and dynamic fields, far offsets, `k*h=1`, and
the vacuum identity at epsilon_r=1.
It checks products against that dense matrix and recovers known solutions
through FFT-backed GMRES, requiring true residual at most `1e-11` and
solution error below `2e-8`. Low budgets, duplicate sites, coordinates outside
the grid, padded-size overflow and electrically oversized cells are rejected.

The FFT, kernel and GMRES test programs passed. The new FFT checks also
passed AddressSanitizer and UndefinedBehaviorSanitizer; leak detection was
disabled after LeakSanitizer failed in this sandbox. Strict C99 syntax/warning checks passed for the new sources.
An intentional one-iteration benchmark returned exit status 1, reported
`iteration limit`, and suppressed solved observables. The distribution build
also passed.

## Measurements on the current 2 GB machine

The following exploratory wall times were measured on 2026-10-19, with other
validation/development work running concurrently. They establish feasibility,
not a controlled performance comparison or a production runtime guarantee.
All cases converged using GMRES(30) and complex Jacobi preconditioning.

| Model | Target residual | Iterations | True residual | Build (s) | Solve (s) | Peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Skeleton, 216 cells | 1e-8 | 14 | 3.34e-9 | 7.91 | 11.96 | 14,884 |
| Skeleton, 216 cells | 1e-10 | 17 | 4.96e-11 | 8.86 | 18.89 | 14,856 |
| Tree, 9,432 cells | 1e-8 | 26 | 9.83e-9 | 40.48 | 238.22 | 126,384 |
| Tree, 9,432 cells | 1e-10 | 33 | 9.80e-11 | 39.78 | 242.96 | 129,840 |

The full-tree plan owns 100,742,992 bytes; GMRES workspace accounts for
28,538,944 bytes and caller arrays for 1,584,576 bytes. Their sum is
130,866,512 bytes (124.8 MiB). The measured tree matvec took 9.66 seconds.
This replaces the sizing estimate of 25.6 GB for just two dense volume
matrices with a measured approximately 123 MiB peak RSS for this experiment.
It does not establish the cost of adding conductors or refining the mesh.

For the skeleton, tightening from `1e-8` to `1e-10` changed the sampled complex
far fields by at most `3.51e-9`, or `8.30e-10` relative to the largest sampled
field; integrated normalized currents changed by `1.51e-9` on the same scaled
basis. At the tighter tolerance, 16-versus-32-node angular scattering power
changed by `4.13e-15` relative, and extinction/scattering differed by
`1.26e-14` relative. These checks establish iterative and integration
consistency; physical mesh convergence and antenna accuracy remain separate.

The full tree also converged at `1e-10`. The maximum sampled complex
far-field change between the loose and tight runs was `5.34e-07`, or `2.4e-09`
relative to the largest sample. Integrated normalized currents changed by
`1.02e-09` on that scaled basis. The tighter tree power calculation changed
by `7.32e-15` relative under angular refinement, and extinction/scattering
differed by `2.84e-11` relative. These are consistency checks across the
tolerance/fallback-order change described below, rather than a pure residual
sweep with a frozen operator. Its peak RSS was 129,840 KiB (126.8 MiB).

Raw measurements are retained in `tests/m9-fft-*.txt`. The first loose tree
and skeleton runs preceded angular-power postprocessing; the tighter runs
include it. The first loose runs used certified fallback far order 8, and
subsequent runs use order 16; the adaptive smooth far rule and near rules are
unchanged. Dense-reference product checks use far order 16 independently.

## Gates left at the isolated-experiment checkpoint

This experiment can establish whether regular-grid FFT is practical for the
selected tree. It does not close M9 by itself. Integration with the production
mesh requires explicit common-lattice validation. Antenna coupling must avoid
quadratic volume storage while retaining the validated NEC coupling blocks.
Dense-versus-accelerated wire, patch and mixed fixtures, residual/impedance/
pattern sweeps, a measured backend decision and an agreed runtime budget
remain required. The 2 ohm physical accuracy allowance is not a license for
FFT-versus-dense errors of that size; acceleration errors should be much
smaller than discretization errors. OpenEMS is optional for this operator
validation; the separate M8 physical-comparison gates remain unchanged.

The later integration note tracks implementation and validation of these gates.
