# M5 progressive-refinement report

Date: 2026-10-18

The reproducible driver is `tools/m5_refinement.py`.  It generates temporary
`DF` resources and NEC decks, invokes the public executable, and extracts the
reported dipole, scattered power, residual, and lossless power imbalance.

## Fixed cube, basis refinement

The same 4 cm cube with `epsilon_r=3` was represented by successively finer
equal cuboids at 100 MHz.  Thus this table changes the field basis without
changing the physical boundary.

| Cells per edge | Cells | `abs(p_y)` (C m) | Change from prior | Scattered power (W) | Residual | Power imbalance |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 6.803700005e-16 | — | 8.017130000e-12 | 1.559e-16 | 1.693e-6 |
| 2 | 8 | 6.868400005e-16 | +0.9510% | 8.170290000e-12 | 2.123e-16 | 1.686e-6 |
| 3 | 27 | 6.934700005e-16 | +0.9653% | 8.328710000e-12 | 5.154e-16 | 1.685e-6 |
| 4 | 64 | 6.968100005e-16 | +0.4816% | 8.409130000e-12 | 8.365e-16 | 1.685e-6 |

The one-cell value is therefore a useful smoke test, not the converged
physical-cube answer.  The last dipole increment halves, all systems solve
near machine precision, and power closure is stable rather than degrading.

## Rayleigh sphere, boundary refinement

Center-selected staircase spheres used requested `k0*a=0.08`.  Because simple
center inclusion alternately over- and under-fills a sphere, comparison uses
the radius of a sphere having the occupied cuboid volume.  This separates the
large, visible boundary-volume error from the volume-operator error.

| Voxels/diameter | Cells | Occupied/exact volume | Equivalent radius (m) | Rayleigh dipole error | Mean internal-field error | Residual | Power imbalance |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 8 | 1.90985932 | 2.481401964e-2 | +1.1556% | +1.1556% | 3.218e-16 | 1.670e-6 |
| 3 | 19 | 1.34397507 | 2.207125697e-2 | +2.1818% | +2.1818% | 5.650e-16 | 1.673e-6 |
| 4 | 32 | 0.95492966 | 1.969490044e-2 | +1.3027% | +1.3027% | 5.661e-16 | 1.677e-6 |
| 5 | 81 | 1.23758884 | 2.147280549e-2 | +1.2366% | +1.2366% | 9.309e-16 | 1.673e-6 |

The nonmonotone raw volume is the expected staircase-boundary oscillation.
After equivalent-volume normalization, all four representations are within
2.2% of the analytical Rayleigh dipole/internal-field result.  For production
spheres, Antenna Turtle should prefer boundary-aware cuboid generation rather
than center-selected voxels.

## Dynamic sphere, Lorenz-Mie refinement

The requested-radius size parameter is `k0*a=0.8`.  Exact Mie power is
evaluated at each occupied-volume-equivalent radius.

| Voxels/diameter | Cells | Occupied/exact volume | Equivalent `k0*a` | Scattered power (W) | Mie error | Residual | Power imbalance |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 19 | 1.34397507 | 0.88285028 | 5.614780000e-7 | -1.3372% | 4.750e-16 | 1.406e-7 |
| 4 | 32 | 0.95492966 | 0.78779602 | 2.909870000e-7 | +0.6660% | 8.285e-16 | 2.626e-7 |
| 5 | 81 | 1.23758884 | 0.85891222 | 4.896890000e-7 | +1.1912% | 1.084e-15 | 1.606e-7 |

All dynamic refinements are within 1.34% of exact Mie total scattering.  The
linear residual remains below `1.1e-15`, and extinction/scattering closure is
better than `2.7e-7`.

## Issue found and disposition

The first dynamic sweep exposed a face/edge-adjacent interaction whose
order-56/order-64 change was about `1e-7` relative.  It correctly failed the
standalone M4 kernel's `2e-8` certification threshold, even though that error
is orders of magnitude below the volume-discretization error in M5.

M4's strict default remains unchanged.  M5 now explicitly requests relative
and absolute interaction tolerances of `5e-7` and `5e-10`, respectively, and
reports the exact cell pair, relation, order, absolute change, and component
if convergence still fails.  The 32-cell case that already passed the strict
tolerance retained the same printed observables after this change.  This is a
deliberate separation between kernel certification and whole-model error
budget, not silent acceptance of an unconverged integral.

## Conclusion

These results are a credible starting point for M6: refinement changes are
small and trending downward for fixed geometry, analytical sphere errors are
at the expected low-single-percent staircase level, exact dynamic scattering
is within 1.34%, residuals are near machine precision, and lossless power
closure is substantially tighter than the physical discretization error.
They do not establish production-scale accuracy for arbitrary tree geometry;
that still requires per-model mesh refinement and later openEMS comparisons.
