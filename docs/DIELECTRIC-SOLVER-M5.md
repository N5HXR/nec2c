# M5 dielectric-only dense JVIE solver

## Result

M5 is the first runnable dielectric result.  A free-space NEC deck containing
no `GW` or surface-patch geometry may load a `DF` resource, excite it with an
ordinary NEC plane-wave `EX` card, and request scattered fields with `RP`.
The run reports the polarization-current density in every active cell, a
test-only cell-averaged internal electric field, total electric dipole moment,
dense-system residual, extinction/scattering power balance, and scattered
far-field samples.

Wire, conducting-patch, and ground coupling remain deliberately rejected.
Those change the block system and belong to M6 and later milestones.  Merely
loading and meshing a `DF` resource beside conductor geometry without an
execution card remains valid, preserving the M2/M3 inspection workflow.

## Discrete system

For time dependence `exp(+j omega t)`, susceptibility
`chi = epsilon_r - 1`, and polarization current
`J = j omega epsilon_0 chi E`, the constant-cell system is

```text
[epsilon_r,t delta_ts I - chi_t N_ts] J_s
    = j omega epsilon_0 chi_t <E_inc>_t .
```

`N_ts` is the certified, volume-averaged cuboid Galerkin interaction from M4.
The plane-wave average is integrated exactly over each rectangular cell as a
center phase times three sinc factors.  Unlike a point-sampled RHS, this stays
consistent when a cell is electrically non-negligible.

The matrix is assembled densely and solved by complex LU factorization with
partial pivoting.  Only one kernel direction is evaluated for each cell pair;
the other follows from volume-scaled reciprocity,

```text
V_t N_ts = V_s transpose(N_st).
```

This is a reference implementation, not the eventual tree-scale algorithm.
Memory is quadratic and factorization work is cubic in the three-times-cell
count.  M3's cost report is therefore an important preflight check.

M4 retains its strict standalone interaction tolerances.  M5 uses a documented
whole-matrix budget of `5e-7` relative and `5e-10` absolute per interaction,
which remains well below observed volume-discretization errors.  See
`M5-REFINEMENT-REPORT.md` for the progressive study and the convergence case
that motivated this separation.

The far field uses the exact Fourier integral of each constant cuboid current:

```text
r E_s(rhat) = j omega mu_0 / (4 pi)
              [ (integral J exp(+j k rhat.r') dV') cross rhat ] cross rhat .
```

Extinction uses `0.5 Re integral(E_inc dot conj(J)) dV`; scattered power is
integrated from the far field on a deterministic 974-direction Fibonacci
sphere.  Their relative difference is printed for every lossless run.

## Deck interface

```text
CE
GE
FR 0 1 0 0 100.0
EX 1 1 1 0 90.0 0.0 90.0 0.0 0.0 0.0
DF "object.diel"
RP 0 91 1 0000 0.0 0.0 1.0 0.0 0.0 0.0
EN
```

M5 accepts NEC plane-wave types 1 (linear), 2 (right elliptic), and 3 (left
elliptic).  Frequency sweeps use the existing `FR` additive or multiplicative
semantics.  An `RP` card prints range-normalized complex `E_theta`, `E_phi`,
and bistatic RCS.  With no `RP`, `EN` still performs the solve and prints the
currents, dipole, residual, and power check.

## Automated validation

`tests/dielectric_solver.py` exercises the executable rather than calling
solver internals.  Its acceptance gates are:

- all-vacuum geometry produces zero unknowns, polarization, and scattering;
- the electrically small one-cell cube agrees with the analytical
  `3 chi/(epsilon_r+2)` depolarization result to 0.2%;
- the normalized residual is at most `1e-10`;
- the small-object dipole pattern has the analytical half-power sample and
  polarization-axis null;
- Cartesian polarization permutation preserves dipole magnitude to `1e-12`;
- a 32-cell staircase sphere agrees with Rayleigh polarizability to 8% and
  its volume-mean internal field is within 0.05 V/m of
  `3/(epsilon_r+2)`; these intentionally loose bounds include staircase-shape
  error and use the occupied-volume-equivalent radius;
- at requested-radius `k0*a = 0.8`, its total scattering agrees with an
  independent Lorenz-Mie calculator to 2%; and
- lossless extinction/scattering closes to `1e-4` for the cube, `2e-3` for
  the Rayleigh sphere, and 3% for the dynamic sphere.

The same 32-voxel dynamic sphere is compared to a genuinely independent DDA
solve from VoxScatter.  The pinned result in `tests/voxscatter-reference.csv`
was produced from upstream commit
`73f7bba2897ee6d8e6e34812adb40442440a419b` with
`tools/voxscatter_m5_reference.m`.  The total dipole magnitudes must agree
within 3%, a discretization-level comparison rather than an expectation that
DDA point interactions and cuboid Galerkin interactions are identical.

Additional fixed-boundary and staircase-boundary convergence results are in
`M5-REFINEMENT-REPORT.md`; the generator is `tools/m5_refinement.py`.

## Scope limits carried into M6

- free space only;
- real, scalar, positive relative permittivity;
- no embedded or nearby conductors;
- constant Cartesian basis/testing functions;
- M4's `k0` times longest-cell-edge certification limit;
- dense reference solution only; and
- internal-field output is diagnostic, not yet a stable machine interface.
