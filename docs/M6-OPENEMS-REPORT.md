# M6 N2 openEMS comparison

This report records the independent Layer-2 check for the first coupled
wire/dielectric implementation. The comparator is openEMS 0.37.0 using its
EC-FDTD method; it shares no integral-equation code with NEC2C.

## Model

- free-space center-fed dipole, length `0.475 m`, radius `1 mm`, and 23 NEC
  segments;
- lossless cube, side `0.1 m`, center `(0.1, 0, 0) m`, `epsilon_r=3`;
- frequency samples 240, 270, 300, 330, and 360 MHz;
- openEMS dipole represented by two explicit 1 mm PEC cylinders and a 10 mm
  lumped-port gap;
- openEMS `PML_8` on all boundaries, domain approximately 1.2 m per axis,
  Gaussian excitation centered at 300 MHz, `1e-5` end criterion, and up to
  80000 time steps; and
- identical Yee meshes for each vacuum/loaded pair so the small impedance
  difference is not contaminated by a changed grid.

The retained data are `tests/openems-m6-n2-15mm.csv`,
`tests/openems-m6-n2-10mm.csv`, and `tests/nec2c-m6-n2.csv`;
`tools/openems_m6_n2.py` regenerates the openEMS pair. Both FDTD meshes reached
the configured energy-decay criterion.

## Results

At 300 MHz:

| Method | Vacuum Z (ohm) | Loaded Z (ohm) | Delta Z (ohm) |
| --- | ---: | ---: | ---: |
| openEMS, 15 mm cube mesh | 81.013 + j29.829 | 80.655 + j28.743 | -0.357 - j1.086 |
| openEMS, 10 mm cube mesh | 80.894 + j29.770 | 80.510 + j28.677 | -0.384 - j1.094 |
| NEC2C M6, 25 mm volume cells | 72.260 + j1.736 | 71.813 + j0.873 | -0.447 - j0.863 |

The 300 MHz dielectric perturbation differs by `0.064 ohm` in resistance and
`0.230 ohm` in reactance from the finer openEMS result. Linear interpolation
of the 270/300 MHz reactance samples gives:

| Method | Vacuum resonance | Loaded resonance | Dielectric shift |
| --- | ---: | ---: | ---: |
| openEMS, 15 mm | 289.921 MHz | 290.203 MHz | +0.282 MHz |
| openEMS, 10 mm | 289.928 MHz | 290.209 MHz | +0.281 MHz |
| NEC2C M6 | 299.426 MHz | 299.709 MHz | +0.283 MHz |

The absolute resonances differ because the FDTD lumped gap and NEC delta-gap
thin-wire source are not equivalent discretizations, but the independently
computed dielectric-induced resonance movement agrees to about `0.002 MHz`.

## Acceptance and retained caveat

The full dielectric-perturbation sweep is inside the accepted
component-wise `max(2%, 2 ohm)` target. At 360 MHz, NEC2C predicts roughly
`+0.062 - j3.319 ohm`, while the 10 mm openEMS mesh predicts
`+0.225 - j4.740 ohm`; the reactance-shift discrepancy is about `1.42 ohm`.
The openEMS 15-to-10 mm change is itself about `0.39 ohm` there, so the
high-frequency endpoint also needs more independent refinement.

The 1.42-ohm endpoint difference is inside the accepted 2-ohm cross-solver
tolerance, so M6 is complete. The implemented solver, fast Layer-1
suite, 300 MHz comparison, and resonance-shift result are useful evidence for
the block signs and scaling. Future validation should refine the NEC volume
mesh at 330/360 MHz,
refine the openEMS cube mesh again, and isolate source-model baseline effects.
