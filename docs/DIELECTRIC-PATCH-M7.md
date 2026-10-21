# M7 conducting-patch/dielectric coupling

> Historical milestone record. The [revision-1 guide](DIELECTRIC-REVISION-1.md)
> states current supported scope; later milestones extend some restrictions below.


M7 extends the dense conductor/JVIE solve from wire-only geometry to NEC2
conducting surface patches and mixed wire/patch geometry. It remains a
free-space, separated-object feature: dielectric contact, NEC ground, loading,
networks, and final radiation-pattern contributions remain deferred.

## Coupled formulation

The conductor block is NEC2's existing `n + 2m` EFIE/MFIE matrix, transposed
into the row-major convention used by the dense coupled solver. For every NEC
conductor basis vector, `cabc` and `nefld` evaluate the electric field averaged
over each volume cell. This supplies both wire-to-volume and patch-to-volume
blocks without replacing NEC's established patch basis.

The reverse patch block uses the magnetic field of constant electric volume
current. With the project's `exp(+i omega t)` convention,

`H(r) = integral(grad(G(r,r')) x Jp(r')) dv'`,

where `G = exp(-i k R)/(4 pi R)`. A tenth-order tensor Gauss rule integrates
each separated cuboid. The two rows of each NEC patch use the same MFIE tests
as `cmss`: `-s t2.H` and `-s t1.H`, where `s` is NEC's patch orientation sign.
The existing clearance gate prevents singular conductor/volume configurations.

Plane-wave types 1 through 3 and ordinary wire voltage sources are supported.
The solver reports wire currents, feed impedance, Cartesian patch surface
currents, volume polarization currents, the full residual, and norms for each
of the four conductor/volume block directions.

## Fast validation

`tests/dielectric_patch.py` verifies:

- patch-only and mixed wire/patch vacuum identities;
- nonzero isolated wire-to-volume, volume-to-wire, patch-to-volume, and
  volume-to-patch blocks in a mixed model;
- patch-to-volume proportionality to susceptibility and material-independent
  volume-to-patch Green coupling;
- residuals no greater than `1e-10`;
- decreasing patch-current perturbation with separation;
- volume refinement from 8 to 27 cube cells; and
- surface refinement from a 2 by 2 to a 4 by 4 PEC plate mesh.

The area-integrated patch-current perturbation changes by about 4.23% under
the documented surface refinement. The 8-to-27-cell volume refinement changes
the patch-current perturbation by less than 5%.

The original plan called this a patch-block "reciprocity" test. A direct block
transpose comparison is not a valid reciprocity statement here: the JVIE row
tests electric field while NEC's patch MFIE row tests magnetic field. M7
therefore tests the analytically required material scaling, vacuum identity,
geometric symmetry through the plate fixtures, and the independent full-model
comparison instead of asserting a false EFIE/MFIE transpose identity.

## Independent comparison

The retained M7 model contains a 0.5 m center-fed dipole, a separated 0.4 m
square PEC sheet, and a separated 0.1 m `epsilon_r=3` cube at 300 MHz.
`tools/openems_m7_mixed.py` represents the same objects using openEMS 0.37.0,
a 10 mm lumped feed gap, PML boundaries, and identical vacuum/loaded meshes.
Both the 15 mm and 10 mm FDTD runs reached their `-50 dB` decay criterion.

| Solver | Vacuum Z (ohm) | Loaded Z (ohm) | Dielectric delta (ohm) |
| --- | ---: | ---: | ---: |
| NEC2C, 8 volume cells | 96.578 + j16.582 | 96.519 + j17.193 | -0.060 + j0.611 |
| NEC2C, 27 volume cells | 96.578 + j16.582 | 96.518 + j17.197 | -0.061 + j0.615 |
| openEMS, 15 mm cube mesh | 68.808 + j125.361 | 68.497 + j125.194 | -0.312 - j0.167 |
| openEMS, 10 mm cube mesh | 69.296 + j125.536 | 69.057 + j125.417 | -0.240 - j0.119 |

The absolute impedances differ substantially because the NEC delta-gap and
FDTD lumped-gap feeds, wire representations, and zero-thickness sheet
discretizations are not equivalent. The dielectric perturbations agree within
0.18 ohm in resistance and 0.74 ohm in reactance, inside the accepted 2-ohm
cross-solver tolerance. The reactance perturbations have opposite signs, so
this is retained as a refinement caveat rather than evidence of sub-ohm
accuracy. Raw values are in `tests/nec2c-m7-mixed.csv` and
`tests/openems-m7-mixed-{15mm,10mm}.csv`.
