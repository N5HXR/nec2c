# Coated-wire M7 practical bend envelope

Status: accepted; the qualified envelope is enforced by default.

## Question and decision rule

This addendum tests whether coated bends can be made practically reliable by
constraining the maximum direction change at a centerline vertex and the
minimum length of every bend chord relative to the outer coating radius.  The
initial candidate is `X=45 degrees, Y=4`: no vertex turns more than 45 degrees
and no chord is shorter than four coating radii.

Agreement is reported in two deliberately broader engineering bands:

| Band | Resonance | Resistance | Coating-shift error |
| --- | ---: | ---: | ---: |
| High confidence | <= 2% | <= 10% | <= 1 percentage point |
| Basically right | <= 5% | <= 25% | <= 2 percentage points |

These bands describe a converged comparison; they never substitute for
reference convergence.  Results outside 5% resonance agreement are qualitative
only.  This is intentionally less strict than the straight-dipole milestones,
because the goal is an explicitly bounded engineering approximation comparable
in spirit to NEC's documented ground-model limitations.

## Geometry battery

Every article is a center-fed, planar Z with a 474 mm centerline, 2 mm PEC
radius, and two nominal 90-degree bends.  The coated version has a 4 mm outer
radius, relative permittivity 3.2, and zero loss.  Keeping centerline length
constant prevents a fillet from appearing better merely because it shortened
the antenna.

The matrix contains a sharp control plus 1, 2, 3, 4, and 6 chords per corner at
bend radii of 1, 2, 4, 6, 8, and 12 coating radii: 31 geometries and 62 NEC
bare/coated sweeps.  A representative nine-case subset is independently solved
as a finite-radius PEC surface in MATLAB Antenna Toolbox.  OpenEMS uses true
cylinders, a Cartesian Yee mesh, a symmetric terminal-current contour, and the
same exact `1e-4` transient-energy criterion as the preceding milestones.

For a circular fillet with `N` chords and radius `R`, the chord-length ratio is

```text
L / b = 2 (R / b) sin(pi / (4 N)).
```

Consequently, the `X=45, Y=4` requirements are coupled.  The tested envelope
admits a one-chord chamfer at `R/b >= 4`, two chords at `R/b >= 6`, three chords
at `R/b >= 8`, and four chords at `R/b >= 12`.  Six chords cannot meet `Y=4`
inside this article's available half-span.

## NEC and MATLAB bare-wire result

The two independent MoM formulations agree across every representative case.
The maximum resonance difference is 0.432%; the maximum resistance difference
is 0.527 ohm.  This includes the sharp 90-degree control, for which MATLAB finds
317.685 MHz/49.497 ohms and NEC finds 319.057 MHz/49.074 ohms.

The NEC coating-induced resonance shift is also smooth across all ten cases
that satisfy `X=45, Y=4`, spanning only -6.596% to -6.521%.  This consistency is
encouraging, but is not independent coated-wire validation.

## OpenEMS convergence correction

The original M7 conclusion attributed a 10.85% sharp-bend discrepancy to NEC's
wire junction.  Mesh refinement and the third solver show that attribution was
incorrect:

| OpenEMS mesh | Resonance | Resistance |
| --- | ---: | ---: |
| 2 mm isotropic | 288.195 MHz | 42.824 ohm |
| 1 mm isotropic | 308.704 MHz | 46.631 ohm |
| 0.5 mm transverse, 1 mm axial | 313.373 MHz | 47.796 ohm |
| MATLAB surface MoM | 317.685 MHz | 49.497 ohm |
| NEC thin wire | 319.057 MHz | 49.074 ohm |

The finest practical FDTD result is 1.78% below NEC in resonance and 2.61%
below it in resistance.  It is 1.36% below MATLAB in resonance and 3.44% below
it in resistance.  Thus the bare sharp bend already fits the high-confidence
band; the 2 mm OpenEMS point was a coarse-grid artifact, not a valid reference.

## Rounded coated-wire result

The exact practical candidate, a one-chord 45-degree chamfer at `R/b=4`, was
then run as matched bare/coated OpenEMS pairs.  Its chord length is 5.657 outer
coating radii, so it satisfies both proposed constraints.

| Solver/mesh | Bare resonance | Bare resistance | Coated resonance | Coating shift |
| --- | ---: | ---: | ---: | ---: |
| OpenEMS, 2 mm | 282.354 MHz | 46.116 ohm | 269.664 MHz | -4.494% |
| OpenEMS, 1 mm | 302.270 MHz | 50.276 ohm | 283.401 MHz | -6.242% |
| OpenEMS, 0.5 mm | 306.377 MHz | 51.553 ohm | -- | -- |
| MATLAB surface MoM | 312.185 MHz | 53.528 ohm | -- | -- |
| NEC | 311.803 MHz | 53.001 ohm | 291.236 MHz | -6.596% |

The expensive 0.5 mm bare run used 61.85 million cells and took 96.5 minutes.
It is 1.74% below NEC in resonance and 2.73% below it in resistance, and 1.86%
and 3.69% below MATLAB respectively.  This establishes absolute geometry
convergence inside the high-confidence band without exceeding the two-hour
point budget.

At the finest *matched* grid, the coating-induced shift differs from NEC by
only 0.354 percentage point.  The OpenEMS resistance change is -10.283%, versus
-10.781% in NEC, a 0.498-point difference.  The shift also moves toward NEC
from the 2 mm to the 1 mm grid.  A coated 0.5 mm run would cost roughly another
two hours while adding little decision value, so it was not run.

## Recommendation

The evidence supports `X=45, Y=4` as a conservative initial engineering
envelope, with the strongest independent physical validation at the
one-chord, `R/b=4` boundary.  Bare NEC and MATLAB agree throughout the sampled
battery, the fine bare FDTD reference converges, and the matched coating effect
meets the high-confidence differential band.

This does not justify arbitrary coated bends.  The NEC operator remains a
thin-coating approximation, and only one member of the wider admissible set
received the costly coated FDTD comparison.  The recommended policy is:

* accept vertices of at most 45 degrees when every adjacent bend chord is at
  least four outer coating radii long;
* warn, but do not silently alter geometry, when a modeling front end offers
  automatic filleting;
* reject tighter coated bends by default while retaining the explicit
  research switch; and
* describe the envelope as an engineering approximation with 2% resonance,
  10% resistance, and 1-percentage-point coating-shift validation bands.

nec2c now applies that policy directly to the resolved wire centerline. NEC
subdivision does not make a geometric chord artificially short: the validator
walks the complete collinear run from each bend to the next bend or open end.
Qualified jobs are stamped `lossless-qualified-bend` or
`lossy-qualified-bend`, and their machine result records the observed maximum
turn, minimum chord/radius ratio, and validated limits. Tighter bends still
produce `COATED_WIRE_BEND_UNVALIDATED` unless the explicit experimental switch
is present. nec2c never changes the user's geometry; optional warning and
auto-fillet UX belongs in a modeling front end such as Antenna Turtle.

## Reproduction

The retained data, MATLAB configuration, comparison plot, and deterministic
checks are documented in `tests/coated-wire/README.md`.  Run the compact gate
with:

```sh
python3 tests/coated-wire/m7_bend_envelope.py
```
