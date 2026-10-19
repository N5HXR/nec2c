# M4 cuboid Galerkin interaction kernel

## Scope and convention

This note defines the M4 interaction computed by `vie_cuboid_n()`. It is a
standalone numerical building block, not yet an enabled dielectric solve.

The time convention is `exp(+j omega t)` and the outgoing scalar Green
function is

```text
g(r,r') = exp(-j k0 |r-r'|) / (4 pi |r-r'|).
```

For a vector field `u`, define

```text
S[u](r) = integral_Vs g(r,r') u(r') dV'
N[u]    = curl curl S[u].
```

This agrees with the current-based JVIE convention used by Polimeridis et al.
and Georgakis et al. For electric polarization current `J` with units A/m^2,

```text
E_sca = (N - I) J / (j omega epsilon0).
```

M4 returns the volume-averaged Galerkin representation of `N`. Its entries are
dimensionless. Keeping `N` dimensionless avoids directly evaluating a
singular `1/omega` expression as frequency approaches zero. M5 will assemble
the purely dielectric current equation as

```text
[epsilon_r I - chi N] J = j omega epsilon0 chi E_inc,
chi = epsilon_r - 1.
```

## Piecewise-constant Galerkin reduction

Let `e_p` be Cartesian component `p`, constant on testing cuboid `Vt`, and
`e_q` be constant on source cuboid `Vs`. Appendix C of Georgakis' thesis
reduces the six-dimensional PWC interaction to the single closed-surface term

```text
N_pq = (1 / volume(Vt))
       integral_boundary(Vt) integral_boundary(Vs)
       (nt x e_p) dot (ns x e_q) g(r,r') dS' dS.
```

There is no center-point approximation. The formula includes the distributional
self term of `curl curl S`; consequently `(N-I)` must be used when recovering
the scattered electric field. For faces normal to axes `a` and `b`, with
normal signs `sa` and `sb`, the constant vector coefficient is

```text
delta(p,q) delta(a,b) sa sb - delta(a,q) delta(b,p) sa sb.
```

Only parallel face pairs contribute to diagonal tensor entries. Perpendicular
face pairs contribute to the corresponding off-diagonal entry. The raw
Galerkin matrix is reciprocal, so the volume-averaged result satisfies

```text
volume(Vt) N(Vt,Vs) = volume(Vs) transpose(N(Vs,Vt)).
```

The static self interaction of a cube is an analytical check:

```text
N(cube,cube,k0=0) = diag(2/3, 2/3, 2/3).
```

Its off-diagonal entries vanish by reflection symmetry. More generally, the
three static diagonal self entries of a rectangular cuboid sum to two.

## Reduced quadrature

Every face pair is reduced independently before numerical integration. Along
an axis where both rectangles vary, the two integrations become one difference
coordinate with the exact cross-correlation weight of the two intervals. Along
an axis where only one face varies, a shifted one-dimensional interval remains;
an axis where neither varies contributes a fixed squared separation. Parallel
faces therefore require at most a two-dimensional integral and perpendicular
faces at most a three-dimensional integral.

Intervals are split at every cross-correlation breakpoint and at zero. If the
origin is a corner of a reduced integration box, that box is partitioned into
radial Duffy sectors. The Jacobian cancels the weak `1/R` singularity for self,
face-, edge-, and corner-adjacent interactions. For a small positive fixed
gap, an additional power map clusters radial nodes near the resulting narrow
peak. This is required for continuity as a separated pair approaches contact.

The implementation classifies cell pairs as self, face-adjacent,
edge-adjacent, corner-adjacent, near, or far. Self, adjacent, and near
interactions compare increasing Gauss-Legendre orders from 4 through at most
64. The default stopping rule is componentwise equivalent to

```text
abs(new-old) <= 2e-11 + 2e-8 * max(abs(new),abs(old)).
```

Failure to satisfy it is an explicit kernel failure. Far interactions use a
fixed order of eight after oracle certification. The public options permit a
caller to request tighter tolerances or a different far order without changing
the formulation. The certified M4 domain requires `k0 * hmax <= 1` for both
cuboids, matching the M3 hard mesh limit; larger cells fail explicitly.
Distinct cuboids with overlapping interiors are also rejected, consistent
with the M3 geometry contract.

## Independent oracle and validation

`tests/generate_vie_reference.py` is independent of the C implementation. It
uses NumPy's separately implemented Gauss rule, vectorized reduced integrals,
and a separate Duffy transformation. It compares full tensors at orders 36 and
52 before emitting `tests/vie-reference.csv`. The checked reference set covers
dynamic self, face-, edge-, and corner-adjacent equal cuboids plus unequal near
and far cuboids. It varies dimensions, aspect ratios, separation, and
electrical size.

The singular surface evaluator was independently cross-checked against the
LGPL DIRECTFN quadrilateral implementation at `k0=2.1 rad/m`. The optional
`tools/directfn_surface_oracle.cpp` harness, built against an unmodified
DIRECTFN checkout at quadrature order 32, produced:

| Topology | DIRECTFN surface integral |
| --- | --- |
| same rectangle | `0.0033606454990071947 - j0.00059210721700740372` |
| common edge | `0.00094972271617248325 - j0.00055397084030742432` |
| common vertex | `0.00063384646893531579 - j0.00053738142780421466` |

The independent reduced oracle at order 128 differs by at most `5.2e-12` in
absolute value across those three cases. DIRECTFN is an oracle only; NEC2C
does not link it and gains no new runtime dependency.

The automated C test additionally checks the analytical static cube,
translation and Cartesian-permutation invariance, volume-scaled reciprocity,
electrical-scale invariance, continuity at face contact from gaps of `1e-3 m`
and `1e-5 m`, stable direct evaluation at `k0=0`, and deterministic rejection
of invalid options.

## Provenance

- A. G. Polimeridis, J. F. Villena, L. Daniel, and J. K. White, “Stable
  FFT-JVIE solvers for fast analysis of highly inhomogeneous dielectric
  objects,” *Journal of Computational Physics* 269 (2014), 280–296,
  <https://doi.org/10.1016/j.jcp.2014.03.026>.
- I. P. Georgakis, *Fast Integral Equation Methods for the Analysis of
  Electromagnetic Wave Interactions with Biological Tissue*, Appendix C,
  Skolkovo Institute of Science and Technology (2019),
  <https://back.skoltech.ru/storage/app/media/defenses/2019/ioannis-georgakis/revised/thesis31.pdf>.
- A. G. Polimeridis et al., DIRECTFN quadrilateral singular-integration
  library and algorithms, <https://github.com/thanospol/DIRECTFN> and
  <https://arxiv.org/abs/1703.08146>.
