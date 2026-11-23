# NEC4 Similarity M10D: Buried Radial-Wire Ground Screen

## Scope

M10D reproduces the buried radial-screen family in Figure 5 and the junction
qualification in the appendix of the 1984 Burke/King/Miller surface-wave
report.  The antenna is a `0.25 lambda` vertical monopole feeding `N` uniformly
spaced radial wires.  Monopole and radial radius are `1e-6 lambda`.  Each
radial leaves the common interface junction through a locally refined
transition cell and reaches the published `1e-5 lambda` burial depth; the
remaining screen is horizontal.  Screen radius, radial count, and the four
published grounds
`epsilon_r = 10-j1`, `10-j10`, `10-j100`, and `10-j1000` are varied.

The prescribed field and Relative Communication Efficiency normalization are
the same as M10C: `E_z` at the interface 1000 free-space wavelengths away,
divided by accepted input power and by a separate centre-fed `0.1 lambda`
reference dipole evaluated through the **same ground-field path** at the same
location.  The report's tabulated reference value is retained as an absolute
field diagnostic, but it is not substituted for the computed denominator.
This follows the report's definition and cancels the common surface-wave
propagation factor rather than mixing NECGS and nec2c field paths.

The primary model is the *buried* screen.  The report explicitly says that a
screen only `1e-5 lambda` above ground converged poorly unless its junction
segments were as short as the screen height.  That above-ground construction
is retained as a documented risk case, not silently mixed into the qualified
buried-screen result.

## Interface-junction contract

The common junction contains one upper-medium monopole ray and `N`
lower-medium radial rays.  It uses the multiwire form of the public NEC4
interface condition already qualified for a one-above/one-below crossing:

```text
upper endpoint factor = 1
lower endpoint factor = epsilon_ground
I'_upper / I'_lower = 1 / epsilon_ground
```

Current obeys Kirchhoff's law across all rays, while every lower-medium radial
receives the same complex-permittivity endpoint factor.  Because the
interface condition prescribes those ratios, the homogeneous local
junction-charge solve must not replace them or report a fallback.  Only the
one-upper/many-lower radial-screen topology is added here; arbitrary branched
interfaces remain unsupported.

## Numerical experiment

1. Unit tests cover the 1+1 interface group, frequency updates, endpoint
   factors, invalid unsplit wires, interface-plane wires, and regular-fan
   recognition.
2. A small explicit-radial model establishes source, junction, and result
   semantics before the published sweep.
3. Screen-radius curves sample sparse through dense screens and retain the
   published `N = 5, 10, 20, 50, 100` identities where computationally useful.
   The evidence file records every actually executed combination; it does not
   interpolate uncomputed NEC results.
4. Selected sparse, intermediate, and dense points are refined while the
   driven segment remains fixed.  Radial refinement changes the radial mesh,
   not the voltage-gap source model.
5. One- and multi-worker runs must preserve ordered numerical results.  Cold
   ground-table construction and warm antenna solves are timed separately.

## Acceptance gates

- the interface crossing has the documented permittivity-factor ratio and the
  buried regular fan is explicitly identified with no local-junction fallback;
- every retained model has finite fields, positive accepted input power, and
  finite feed impedance;
- the curves reproduce the published ordering: a denser screen never performs
  worse than the sparsest screen at the frozen review radii, moderate/poor
  ground shows the material screen benefit, and highly conducting ground shows
  a much smaller benefit;
- the sparse low-loss screen exhibits its published finite-radius optimum or
  oscillatory decline, while dense-screen curves approach a plateau;
- selected fixed-source-cell refinements stay within a tolerance frozen from
  the measured convergence sequence, with plot-reading uncertainty reported
  separately from numerical convergence;
- direct and cached M9 ground paths agree at selected points, and the full
  optimized regression suite remains green.

No empirical RCE scale, NEC4.2 grazing repair, or surface-impedance correction
is applied to make the curves match.  The Wait-Pope curves in Appendix A are
an independent trend check, not replacement physics for the wire solution.
