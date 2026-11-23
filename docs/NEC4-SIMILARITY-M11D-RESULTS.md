# NEC4 Similarity M11D: Direct Magnetic Field

## Outcome

M11D implements and qualifies an independently sampled magnetic field for the
public `adaptive` ground evaluator.  It does not compute H by differencing an
interpolated E table.  The direct H evaluator and its adaptive table are
complete; production routing remains M11E work.

The implementation follows the existing `exp(+j omega t)` passive-medium
convention and Maxwell equation

```text
curl(E) = -j omega mu H.
```

## Cylindrical kernels

The five Section 5 electric kernels form a cylindrical dyadic.  A horizontal
x-directed source has

```text
E_rho = E_rho_h cos(phi)
E_phi = E_phi_h sin(phi)
E_z   = E_z_h cos(phi).
```

Consequently H has four independent scalar kernels:

```text
H_phi(vertical source)
H_rho(horizontal source) sin(phi)
H_phi(horizontal source) cos(phi)
H_z(horizontal source) sin(phi).
```

For the positive-axis Bessel representation, let

```text
gamma_i = sqrt(lambda^2-k_i^2)
A       = exp(gamma_1 z' - gamma_2 z)
D       = k_1^2 gamma_2 + k_2^2 gamma_1
S       = gamma_1 + gamma_2
```

with outgoing sheets, `z'=-depth`, and `J`, `J'`, and `J''` denoting the
order-zero Bessel function and argument derivatives.  Direct analytical curl
of the E integrands gives

```text
H_phi_v = -A J' lambda^2 k_2^2 / (2 pi D)

H_rho_h = -A [(gamma_2-gamma_1) J' lambda^2/(rho D)
              + gamma_2 J lambda/S] / (2 pi)

H_phi_h = -A [(gamma_2-gamma_1) J'' lambda^3/D
              + gamma_2 J lambda/S] / (2 pi)

H_z_h   = -A J' lambda^2 / (2 pi S).
```

The `1/(2 pi)` factor is important: the published E normalization contains
ordinary frequency `f`, while Maxwell curl contains `omega=2 pi f`.  The
full-line Hankel densities carry half these coefficients.  Cylindrical-axis
limits are evaluated analytically instead of dividing by rho.

## Direct qualification

| Gate | Result |
|---|---:|
| Analytic H vs centered numerical curl of integrated E | `2.415e-5` |
| Identical-media closed-form Green function | `3.084e-5` |
| Independent Bessel/Hankel contour overlap | `2.283e-5` |
| Interface Richardson refinement | `1.374e-7` |
| Interface symmetric vs source-only path | `2.971e-9` |
| Fixed-depth transmission as lower medium approaches PEC | `1.466e-7` ratio |
| Cylindrical origin | exact analytic limits passed |
| 1-worker vs 4-worker sampling | bit-identical |

The homogeneous closed form is the reciprocal electric-current-to-magnetic-
field dyadic `H = (1+j k R) exp(-j k R) (I cross Rhat)/(4 pi R^2)`.  It fixes
the source/observer orientation and transverse signs independently of the
Sommerfeld derivation.  The interface curl check then qualifies the same
signs in the two-medium case.

## Adaptive H qualification

H uses its own anisotropic cubic table and direct-H withheld witnesses.  It
shares E's accepted F2 phase and inverse-square storage transform, but no E
samples or numerical E derivatives.

| Domain | Cells | Unresolved | Direct fallback witnesses | Worst independent error | Cold build |
|---|---:|---:|---:|---:|---:|
| grazing sheet | 2,671 | 0 | 0/128 | `2.024e-4` | 82.9 s |
| route-switch volume | 2,815 | 0 | 0/128 | `2.915e-4` | 234.3 s |

The frozen cell and external witness gate is `3e-4`.  The route-switch result
is deliberately retained near the gate as a regression witness.  No empirical
material correction was used.

## Scope boundary

- Direct and adaptive H are qualified numerical primitives in M11D.
- They do not change the current command-line result or production matrix path.
- M11E owns deterministic production routing, diagnostics, persistent cache
  treatment for anisotropic tables, and full-model timing.
- This is an independent NEC4-theory implementation, not a claim of NEC-4.2
  GN3 equivalence.

## Acceptance commands

```sh
make -j8 tests/nec4_ground_direct_test tools/nec4_m11d_qualify
./tests/nec4_ground_direct_test
python3 tests/nec4-similarity/m11d_magnetic.py
```

The retained evidence is `tests/nec4-similarity/m11d-magnetic.json`; the
summary plot is `tests/nec4-similarity/m11d-magnetic.png`.
