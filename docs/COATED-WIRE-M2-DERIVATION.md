# M2 coated-wire parser and discrete-operator specification

Status: accepted 2026-11-05

Date: 2026-11-03

Scope: parsing, topology, metadata, and derivation only; no coated solve

## Delivered parser contract

M2 recognizes the frozen card:

```text
LD 7 TAG FIRST LAST ER OUTER_RADIUS TAN_DELTA SIGMA_DC
```

`LD 7` is stored in `coated_wire_model`; it never increments `zload.nload`,
enters `load()`, or modifies `zload.zarray`.  `LD -1` at the beginning of a new
loading sequence clears both legacy loads and coating designations, preserving
the traditional meaning of “clear loading.”  The optional loss fields default
to zero through the ordinary NEC card parser.  The two remaining floating
fields on the physical NEC record are reserved and must be zero, so a future
extension cannot silently reinterpret an old deck.

`--validate-only` resolves every designated tag to final, one-based global
segment numbers and emits those numbers in result schema 1 under the optional
`job.coated_wires` member.  An ordinary M2 execution is rejected with
`COATED_WIRE_OPERATOR_UNAVAILABLE`: accepting the card and silently solving the
bare structure would be a false result.  M3 will replace that gate with the
lossless operator.

The topology validator currently accepts disconnected coating domains and
contacts between domains whose four material/geometric card values are exactly
equal.  It rejects partial tags, duplicate tags, mixed conductor radii within a
tag, coated/bare contacts, unequal coated/coated contacts, surface-patch
terminations, and direct ground-plane terminations.  The last two constraints
are conservative M2 boundaries rather than claims that the physical structures
cannot be modeled later.

## Source model and phasor convention

Let a PEC wire of radius `a` carry axial current `I(s)` and have a concentric
coating of outer radius `b` and complex permittivity `epsilon_c`.  nec2c uses
`exp(+i*omega*t)`, so M4 will evaluate

```text
epsilon_c / epsilon_0 =
    ER - i*(ER*TAN_DELTA + SIGMA_DC/(omega*epsilon_0)).
```

Current continuity gives the free line charge

```text
lambda(s) = -1/(i*omega) * dI(s)/ds.                 (1)
```

Popovic, Djordjevic, and Kircanski derive the same relation first as wire
surface charge, assume a predominantly radial quasi-static coating field, and
obtain equal-and-opposite polarization-charge sheets at `a` and `b`.  Combining
the free charge with the inner bound charge leaves the following equivalent
line-charge weights:

```text
inner radius a:  (epsilon_0/epsilon_c) * lambda(s)
outer radius b:  (1 - epsilon_0/epsilon_c) * lambda(s). (2)
```

The weights sum to the original free charge.  Therefore the scalar-potential
kernel seen outside the coating is

```text
Kc = (epsilon_0/epsilon_c)*Ga
     + (1 - epsilon_0/epsilon_c)*Gb,                 (3)
```

where `Ga` and `Gb` are the same retarded free-space Green kernel evaluated at
the inner and outer equivalent radii.  Relative to the bare wire,

```text
Delta K = alpha*(Gb - Ga),
alpha = 1 - epsilon_0/epsilon_c.                     (4)
```

This is the two-radius specialization to be implemented first.  It has both
required null limits without an empirical adjustment: `alpha=0` when the
coating is vacuum, and `Gb-Ga=0` when `b=a`.  Richmond and Newman independently
derive the equivalent-polarization-current reaction formulation and show that
it changes the impedance matrix without adding unknowns.  The 1981 paper
reports best experimental agreement for coating thickness no larger than about
the conductor radius and relative permittivity below about 10; those are M3/M7
validation coordinates, not parser limits.

Mowete, Ogunsola, and Sandrolini later emphasize that a finite coating also has
an axial internal electric-field component.  Their result is why equations
(2)--(4) are a reduced model to validate, not a proof that arbitrary thick
sleeves are covered.

## Mapping to the NEC basis

On source segment `j`, NEC represents each connected basis contribution with
the existing coefficients returned by `trio()`/`sbf()`:

```text
I_n(s) = A_n + B_n*sin(k*s) + C_n*cos(k*s)            (5)
```

with the local orientation and coefficient signs already handled by those
routines.  The coating correction for source basis `n` is evaluated from its
charge basis

```text
q_n(s) = -1/(i*omega) * dI_n(s)/ds
        = -k/(i*omega) * (B_n*cos(k*s) - C_n*sin(k*s)). (6)
```

For an observation on test segment `i`, the additional electric field is

```text
Delta E_n(r_i) = -grad integral_support(n)
    [alpha_j * q_n(s') * (Gb(r_i,s') - Ga(r_i,s'))] ds'. (7)
```

The exact normalization and matrix sign will be expressed using nec2c's
dimensionless wavelength coordinates and `CONST1` in M3, then pinned by the two
null-limit tests rather than copied between incompatible phasor conventions.
Equation (7), not an `LD 2` series impedance, defines the required operator.

### Mutual terms

For distinct source and observation segments, apply (7) for every source
segment in the connected support exposed by `trio()`, project the resulting
field on the observation tangent exactly where `cmww()` projects `efld()`, and
accumulate with the corresponding `A/B/C` coefficients.  The correction is a
full wire-wire interaction: applying it only to diagonal matrix entries would
discard coated-wire mutual coupling and violate the source formulation.

### Self term

For `i=j`, retain the radius-separated difference `Gb-Ga` as one numerical
object.  Evaluating two large near-singular integrals and subtracting afterward
would lose precision as `b` approaches `a`.  M3 will use an analytic difference
where practical and paired quadrature otherwise.  The self test must cover
`b/a -> 1`, `ER -> 1`, segment refinement, and the matrix diagonal at the
existing thin-wire and extended-thin-wire settings.

### Endpoints

There are two equivalent implementation forms:

1. integrate the analytic derivative in (6) directly; or
2. integrate (7) by parts, producing a current-weighted kernel derivative plus
   the explicit boundary term

```text
[ I_n(s') * alpha_j * (Gb-Ga) ] at each support endpoint. (8)
```

The boundary term must not be dropped.  At a free wire end it represents the
terminal charge implicit in the current basis and is the reduced operator's
only first-order awareness of the coating termination.  It does **not** by
itself distinguish flush, extended-flat, or hemispherical dielectric caps;
M1's measured cap sensitivity remains part of the model uncertainty.  Direct
derivative quadrature is preferred initially because it retains the same term
without a cancellation-prone explicit endpoint subtraction.

### Junctions

At a wire junction, endpoint contributions from incident oriented supports
cancel only when the same `alpha`, `a`, and `b` apply and the existing NEC basis
enforces current continuity.  That is the mathematical reason for M2's equal-
coating junction rule.  A coated/bare or unequal-coating contact requires a
transition charge/basis treatment and is deferred.  A finite feed gap has no
geometric junction, so two independently coated arms may be excited across the
gap once M5 validates that source model.

### Vector-potential and radiation audit

For the straight, collinear M3 dipole, the azimuthally symmetric radial
polarization current has no axial vector-potential projection on the wire axis;
equation (7) is the promoted first operator.  This simplification is not valid
by assertion for bends, nonparallel mutual coupling, junction transitions, or
finite-end fields.  Before those geometries are enabled, compare the full
equivalent polarization-current reaction term against the scalar-kernel form.
Likewise, do not alter the radiation routines merely to improve a fit: first
measure whether direct radiation from the equivalent coating currents is below
the M3 error budget or must be included consistently.

## Primary evidence

1. J. H. Richmond and E. H. Newman, “Dielectric coated wire antennas,”
   *Radio Science* 11(1), 13–20 (1976),
   <https://doi.org/10.1029/RS011i001p00013>.
2. B. D. Popovic, A. R. Djordjevic, and N. M. Kircanski, “Simple method
   for analysis of dielectric-coated wire antennas,” *The Radio and Electronic
   Engineer* 51(3), 141–145 (1981),
   <https://www.worldradiohistory.com/UK/British-Institution-of-Radio-Engineers/1981/TREE-1981-03.pdf>.
3. A. I. Mowete, A. Ogunsola, and L. Sandrolini, “A quasi-static theory for
   dielectric-coated thin-wire antenna structures,” *PIER Letters* 20, 45–54
   (2011), <https://doi.org/10.2528/PIERL10121104>.

The first two papers define the implementation family.  The third is contrary
evidence used to define validation risks.  Secondary insulated-wire recipes are
not used as the derivation authority.

## M2 acceptance commands

```sh
python3 tests/coated-wire/m2_parser.py

./nec2c --validate-only \
  --result-json /tmp/coated-m2.json \
  -i tests/coated-wire/m2-valid.nec \
  -o /tmp/coated-m2.out
```

Expected: the test reports all parser/topology cases passed; the sample result
has `validation: "coated-wire-preflight"`, one resolved 11-segment domain, and
no execution rows.

## Automated verification

The completed implementation passes:

* four valid syntax/topology/metadata cases;
* seventeen stable invalid-card and invalid-topology diagnostics;
* the deliberate ordinary-execution gate;
* normal execution through every legacy loading type `LD 0` through `LD 5`;
* the immutable M0 bare-wire semantic and traditional-report hashes; and
* the full repository suite: 37 of 37 tests, with no skips or expected failures.
