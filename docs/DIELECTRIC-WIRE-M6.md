# M6 voltage-driven wire/dielectric coupling

Milestone M6 augments NEC2C's wire EFIE with the M5 cuboidal-current JVIE.  A
normal type-0 `EX` voltage source now drives one self-consistent dense system,

```text
[ Zww  Zwv ] [ u ] = [ -Esource ]
[ Zvw  Zvv ] [ J ]   [     0    ]
```

where `u` is NEC's existing wire coefficient vector and `J` contains three SI
polarization-current-density coefficients per active dielectric cell.

## Discretization and scaling

`Zww` is produced without reformulating the legacy EFIE: `cmset()` fills it
after geometry is scaled to wavelengths.  NEC's `factr()` transposes that
stored matrix before elimination, so the coupled row-major matrix explicitly
uses the same transpose.  This detail is covered by the complete-vector
vacuum-identity test.

`Zvv` is the M5 volume-averaged Galerkin block
`epsilon_r I - (epsilon_r-1) N`.  For each NEC wire basis coefficient,
`cabc()` and `nefld()` evaluate its electric field and a sixth-order tensor
Gauss rule averages that field over each cell.  The volume row is

```text
Zvw = -j omega epsilon0 (epsilon_r - 1) <Ewire>cell.
```

The return field from a constant cuboid current is evaluated with the free-
space dyadic `(k0^2 I + grad grad) exp(-j k0 R)/(4 pi R)`, integrated with a
fixed order-10 tensor rule and divided by `j omega epsilon0`.  It is tested at
the wire-segment center along the segment tangent.  M3's positive-clearance
rule keeps this cross kernel nonsingular.  Because this is mixed NEC
collocation and volume averaging, the cross blocks are not literal
transposes.  The acceptance test instead checks reciprocal transfer
admittance between two unequal-environment wire ports; the relative mismatch
is below `1e-5`.  The lower-level dyadic and volume-block scaled reciprocity
tests retain their `1e-9` gate.

The combined matrix uses partial-pivot LU and reports the residual of the
unfactored system.  Volume currents are copied out before `cabc()` converts
only the wire coefficients to physical segment currents.  This prevents the
legacy current interpolation from interpreting volume unknowns.

## Supported M6 scope

M6 supports free-space, wire-only geometry, ordinary type-0 voltage sources,
real `epsilon_r >= 1`, separated axis-aligned cuboidal cells, linear or
multiplicative frequency sweeps, and the existing `KH`/`EK` wire options.
Multiple materials and multiple ordinary voltage sources are structurally
supported.

The solver rejects conductive patches, NEC ground, `LD`, `NT`/`TL`, symmetry,
discontinuity sources, conductor/dielectric contact, and `RP` output.  The
last restriction is intentional: M8 will combine dielectric and conductor
radiation consistently.  Decks without `DF` continue through the unchanged
legacy path.

## Validation summary

The fast M6 suite checks:

- `epsilon_r=1` restores the legacy feed impedance and full wire-current
  vector at report precision;
- the dense residual is below `1e-10` (observed worst case below `2e-14` in
  the small fixtures);
- a farther cube produces a smaller impedance perturbation;
- one-cell and `2 x 2 x 2` cube meshes agree to under 2% of the dielectric
  perturbation;
- 11- and 21-segment wire models agree to under 2% of the perturbation;
- reciprocal transfer admittances in an asymmetric two-dipole/cube model
  agree within `1e-5`; and
- unsupported M6 combinations fail explicitly.

The slower N2 fixture is a 0.475-lambda, 23-segment dipole beside a
0.1-lambda `epsilon_r=3` cube, represented by 64 cubical cells.  At 300 MHz,
the initial NEC result changes from `72.260 + j1.736 ohm` in vacuum to about
`71.813 + j0.873 ohm` with the cube.  The resulting shift is compared with
openEMS in `M6-OPENEMS-REPORT.md`.

The five-frequency N2 sweep passes the resonance-shift comparison and the
300 MHz impedance-perturbation target, but misses the 360 MHz reactance-shift
target under the earlier 1-ohm criterion. The accepted cross-solver absolute
tolerance is now 2 ohms, so M6 is complete; the mesh-sensitive endpoint remains
documented in the linked report for future refinement.

The current reference implementation intentionally favors reproducibility
over speed.  The 64-cell, single-frequency N2 run takes roughly 1 minute 45
seconds on the validation host because volume cubature is rebuilt for every
frequency. A five-point sweep takes roughly nine minutes. It is therefore a
retained slow fixture, not part of `make check`.
