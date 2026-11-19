# NEC4 Similarity M10C Results

## Outcome

The published monopole/ground-stake behavior is reproduced across both
monopole heights, all four grounds, and 41 stake depths.  The test executes
328 complete antenna models plus focused refinement models.  The former
far-interface crash is fixed by the public Appendix J limit rather than an
observer-height offset.

## Published trends

| Ground | Predicted quarter-wave depth | 0.1λ peak | 0.25λ peak |
| --- | ---: | ---: | ---: |
| `10-j1` | 0.07896λ | 0.065λ, -2.539 dB | 0.055λ, 1.339 dB |
| `10-j10` | 0.07196λ | 0.070λ, -3.745 dB | 0.065λ, 1.826 dB |
| `10-j100` | 0.03363λ | 0.030λ, -4.537 dB | 0.030λ, -0.226 dB |
| `10-j1000` | 0.01112λ | 0.010λ, -8.075 dB | 0.010λ, -5.317 dB |

For `10-j1`, both curves place the sharp half-wave-in-ground null at
`0.155 lambda`, versus the analytical estimate `0.1579 lambda`.  The 0.1λ
curve reaches -20.085 dB there.  As loss increases, the peak and null flatten
in the same progression shown by the report.

The taller monopole exceeds the shorter one at every sampled depth.  The
minimum separation by ground is 3.846, 5.567, 4.311, and 2.758 dB in order of
increasing loss.

## Convergence and defects fixed

Refining all cells initially produced multi-dB changes because it also halved
the physical source segment.  That was rejected as an invalid convergence
comparison.  With the driven/interface cells fixed at `0.005 lambda` and only
the remaining cells refined to `0.0025 lambda`, the largest selected RCE
change is `0.0993 dB`, including the sharp low-loss null.

The first buried-stake run exposed a missing simultaneous
`z -> 0`, `Gamma_2 -> 0` limit in the production second-order evaluator.  The
generic formula put the saddle exactly on the `k_2` branch point and aborted.
The Appendix J `Lambda_2` form now evaluates all five field components at the
interface.  Its vertical field agrees with the independent 100-digit
small-height oracle to about `1.9e-8` relative before ordinary double-constant
roundoff.  A related source-on-interface branch strip is now owned by the same
limit, eliminating the failure seen when refinement placed a quadrature point
only `6.6e-6 lambda` above the interface.

Cold cache construction is materially more expensive than the warm antenna
sweep, especially for `10-j1000` (more than five minutes for the final cache
component on this system).  Once cached, the complete 328-model sweep takes
well under a minute with 32 independent one-thread solver processes.  These
are reported separately; cache construction is not presented as antenna solve
time.
