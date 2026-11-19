# NEC4 Similarity M10C: Published Monopole/Ground-Stake Curves

## Scope

M10C reproduces Figure 3 of the 1984 Burke/King/Miller surface-wave report:
0.1- and 0.25-free-space-wavelength vertical monopoles on a buried vertical
ground stake.  Stake depth is swept from zero through `0.2 lambda`; wire and
stake radius are the published `1e-4 lambda`.  The four published grounds are
`epsilon_r = 10-j1`, `10-j10`, `10-j100`, and `10-j1000`.

The prescribed field is `E_z` at the interface, 1000 wavelengths away.  RCE
is `10 log10((|E_z|^2/P_in)/R)`, where `R` is the corresponding Table 1 value
for the report's separate centre-fed 0.1-wavelength reference dipole.  The
base-fed zero-stake monopole is not substituted for that reference.

## Numerical contract

- The driven and interface-adjacent cells remain `0.005 lambda` long during
  convergence checks.  Only the remaining monopole and stake cells are
  refined from `0.005` to `0.0025 lambda`; changing the source cell is a
  different excitation model, not a mesh refinement.
- The complete original-NEC4 second-order interface field is retained.  No
  empirical RCE scale or NEC4.2 grazing correction is applied.
- Appendix J's `Lambda_2 = Gamma_2/z` limit is used when a buried source is
  observed exactly on the interface.  For a source on the interface, the same
  limit owns the numerically unresolved branch-point strip through
  `z/lambda <= 3e-4`; the outer expression is continuous within one percent
  at the qualified high-loss boundary.
- A cold material/frequency model-grid build and a warm antenna sweep are
  separate costs.  Exact-key persistent cache files may be reused.

## Acceptance gates

1. Every one of the 328 antenna points completes with a finite field and a
   positive accepted input power.
2. The measured RCE maximum lies within `0.03 lambda` of the report's stated
   approximate quarter wavelength in ground.
3. For `10-j1`, the sharp minimum lies within `0.01 lambda` of a half
   wavelength in ground.
4. At every common stake depth, the 0.25-wavelength monopole remains at least
   2 dB above the 0.1-wavelength monopole.
5. Fixed-source-cell refinement changes each selected peak and low-loss null
   by no more than 0.15 dB.
6. The equation-level interface limit and its source-interface transition
   remain covered by `tests/nec4_ground_asymptotic_test`.

Machine evidence is in
`tests/nec4-similarity/m10c-ground-stake.json`; the review plot is
`tests/nec4-similarity/m10c-ground-stake.png` and the executable gate is
`tests/nec4-similarity/m10c_ground_stake.py`.  A representative manual deck
at the low-loss 0.1-wavelength peak is
`tests/nec4-similarity/m10c-ground-stake.nec`.
