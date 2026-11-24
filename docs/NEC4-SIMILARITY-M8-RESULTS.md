# NEC4 similarity M8 results

Status: accepted by the user on 2026-11-17.

## Outcome

M8 adds a production-independent, arbitrary-precision electric-field oracle for
a point electric-current element near a planar interface. It evaluates all four
source/observer medium pairings (`12`, `21`, `11`, and `22`), both vertical and
horizontal source orientations, and assembles the complete Cartesian 3x3
electric dyadic. It is deliberately not connected to the production matrix
fill; that integration remains M9 work.

The implementation is based only on the public NEC-4 Theory Guide:

- Section 5.1.1 equations (5-1) through (5-7);
- Figure 5-2 Bessel and Figure 5-3 Hankel contours;
- Appendix F direct field formulas;
- Appendix G equations (G-1) through (G-20); and
- Appendix H equations (H-1) through (H-4) and the residual `1/R`
  subtraction.

## Numerical results

| Gate | Result |
|---|---:|
| derivative versus expanded equation identity | `3.78e-81` |
| Bessel contour versus homogeneous analytic field | `2.44e-15` |
| Bessel independent refinement | `1.16e-11` |
| Hankel pointwise Bessel identity | `3.94e-49` |
| compact Hankel versus Bessel overlap | `1.20e-10` |
| Appendix H literal versus cancellation-resistant density | `1.51e-90` |
| Appendix H large-`lambda` cancellation witness | `1.95e-69` |
| Appendix H integrated reconstruction near coincidence | `9.88e-12` |
| identical-medium Appendix G transformations | `2.14e-14` |
| direct Appendix F `21` versus Appendix G | `6.40e-65` |
| extrapolated tangential-E / normal-D interface conditions | `2.47e-10` |
| PEC image limit at `sigma=1e4 S/m` | `6.82e-5` |
| far-ray Fresnel residual at the farthest witness | `1.39e-2` |
| 65-point, five-medium qualification worst refinement | `3.90e-10` |

The PEC and Fresnel entries are limits rather than same-form identities. Their
residuals decrease monotonically as conductivity and ray distance increase,
respectively. The Fresnel comparison uses the public equation (5-35) dyadic
construction; its residual is expected to retain finite-distance near-field
terms.

The Appendix H test discovered and locked two cancellation-sensitive details:

1. `gamma1-gamma2` must be formed as
   `(k2^2-k1^2)/(gamma1+gamma2)`, not by subtracting the two squared roots.
2. The low-quality scan's final constant in the `S'phi12` expression is `+1`.
   That reading cancels the measured `1/R` coefficient and produces a bounded
   second remainder; the alternative `-1` reading diverges.

The axial `rho=0` path uses analytic Bessel limits. It does not evaluate
`J0'/rho` numerically and enforces the symmetry identities
`E_rho^V=E_z^H=0` and `E_rho^H=-E_phi^H`.

## Artifacts

- `tools/nec4_m8_equations.py`: equation, material, singularity, and analytic
  finite-medium layer;
- `tools/nec4_m8_direct.py`: Bessel/Hankel integration and all-pairing
  transformations;
- `tools/nec4_m8_oracle.py`: diagnostic point and line interface;
- `tools/nec4_m8_qualify.py`: convergence sweep and static-map generator;
- `tests/nec4-similarity/m8-*.json`: machine-readable gate evidence; and
- `docs/nec4-m8-convergence.png`: external-refinement error map.

## User acceptance commands

Run the compact automated gates:

```sh
python3 tests/nec4-similarity/m8_equations.py
python3 tests/nec4-similarity/m8_bessel.py
python3 tests/nec4-similarity/m8_hankel.py
python3 tests/nec4-similarity/m8_singularity.py
python3 tests/nec4-similarity/m8_pairings.py
python3 tests/nec4-similarity/m8_limits.py
python3 tests/nec4-similarity/m8_oracle.py
```

Run one direct diagnostic point in each pairing:

```sh
python3 tools/nec4_m8_oracle.py point --source-z -0.7 --observer-z 0.4 --output /tmp/m8-12.json
python3 tools/nec4_m8_oracle.py point --source-z 0.4 --observer-z -0.7 --output /tmp/m8-21.json
python3 tools/nec4_m8_oracle.py point --source-z -0.7 --observer-z -0.4 --output /tmp/m8-11.json
python3 tools/nec4_m8_oracle.py point --source-z 0.7 --observer-z 0.4 --output /tmp/m8-22.json
```

Run a short line crossing the interface. The exact `z=0` sample is omitted so
the tool never silently chooses a one-sided limit:

```sh
python3 tools/nec4_m8_oracle.py line --line-start -0.8 --line-stop 0.8 \
  --line-count 5 --output /tmp/m8-line.json
```

The deliberately expensive contour qualification is:

```sh
python3 tests/nec4-similarity/m8_hankel.py --full
MPLCONFIGDIR=/tmp/nec4-m8-matplotlib python3 tools/nec4_m8_qualify.py
```

M8 passes after review of the convergence map and a user verdict of `PASS`,
`FAIL`, or `PASS WITH CHANGES`.
