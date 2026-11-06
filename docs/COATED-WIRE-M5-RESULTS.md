# Coated-wire M5 feed-model decision

Status: accepted by the user on 2026-11-06.

## Decision

The supported production form remains the ordinary NEC continuous wire with an
`EX 0` voltage source on its center segment.  The independent physical
reference is a pair of dielectric-coated conductor arms separated by a small
gap and driven by the axisymmetric OpenEMS ring-voltage source qualified in M1.

This is not an assertion that a real antenna has continuous metal through its
feed.  It says the NEC delta-gap source on a continuous current basis is the
appropriate reduced representation of the terminal discontinuity.  Adding a
second NEC geometry and differential-port syntax would add topology and source
unknowns without improving this validated article.

An electric-field source embedded inside uninterrupted PEC is not a useful
alternative OpenEMS reference: the PEC shorts it and the measured terminal
impedance is identically zero.  That failed construction is an OpenEMS modeling
constraint, not a reason to warn on the normal NEC `EX 0` form.

## Reference article and method

The article is the M1/M3 straight dipole: 475 mm total conductor length, 2 mm
conductor radius, 4 mm solid coating radius, relative permittivity 3.2, zero
loss, and flat coating ends flush with the metal.  Both OpenEMS cases used:

* cylindrical coordinates, eight azimuthal cells and the qualified `m=0`
  far-field reconstruction;
* 1 mm maximum axial pitch, 0.5 mm local radial pitch and four cells through
  the coating;
* 185--385 MHz in 0.5 MHz samples;
* the qualified `1e-4` (`-40 dB`) transient energy criterion; and
* nine current contours per arm pinned to explicit Yee-grid planes.

The voltage reference spans the centered physical gap.  OpenEMS measures its
lumped terminal current on the adjacent Yee-cell contour (`z=-1.5 mm` for the
1 mm gap and `z=-2 mm` for the 2 mm gap).  This contour displacement is not a
length of a defined transmission line, so applying a propagation correction
would manufacture a model not present in either solver.  The reported
de-embedding distance is therefore zero, and the 1--2 mm gap sequence directly
bounds the residual source and measurement-plane effect.

## Results

| Quantity | 1 mm physical gap | 2 mm physical gap | 1 mm minus 2 mm |
| --- | ---: | ---: | ---: |
| Resonance | 275.609354 MHz | 275.630067 MHz | -0.020713 MHz (-0.007515%) |
| Resistance at resonance | 64.286833 ohm | 64.328424 ohm | -0.041591 ohm |
| Maximum normalized-current difference away from feed | | | 0.002922 |
| Maximum pattern difference above -40 dB | | | 0.000691 dB |

The NEC continuous-source result is 276.642354 MHz and 64.632156 ohm.  Relative
to the 2 mm physical reference, that is `+0.3673%` in resonance and `+0.303732`
ohm in resistance.  The maximum normalized axial-current discrepancy away from
the feed is `0.005316`.  These are comfortably inside the M3 promotion limits
of 1% resonance, 2 ohms impedance and the observed current-shape tolerance.

The overlays make the result visually clearer than the scalar limits:

* [impedance and resonance](../tests/coated-wire/m5-impedance.png)
* [axial current](../tests/coated-wire/m5-current.png)
* [radiation pattern](../tests/coated-wire/m5-pattern.png)

The two physical-gap curves are effectively coincident.  The NEC current and
pattern also overlay the independent result; most of the remaining solver
difference is the already characterized reduced coating operator, not the feed
representation.

## Scope boundary

M5 certifies one source form only: a center-segment `EX 0` source on the
continuous, uniformly coated straight-wire domain admitted by M3/M4.  It does
not yet certify disconnected NEC arms, unequal hot/cold coatings, a source at a
coating transition, or a coated/bare junction.  Those are M6 topology tests.

No new run-time warning is appropriate for the certified form.  Unsupported
geometries continue to receive the existing hard diagnostic rather than a
silent bare-wire or extrapolated answer.

## Reproduction and acceptance

The dependency-free acceptance test is:

```sh
python3 tests/coated-wire/m5_feed.py
```

For manual NEC inspection:

```sh
./nec2c --backend builtin --jobs 1 --threads 1 \
  --result-json /tmp/m5-continuous.json \
  -i tests/coated-wire/m5-continuous-source.nec \
  -o /tmp/m5-continuous.out
python3 -m json.tool /tmp/m5-continuous.json | less
```

To regenerate the compact record and plots from retained OpenEMS solver
directories:

```sh
python3 tools/coated_wire_m5_acceptance.py \
  --gap1-dir /tmp/m5-openems-gap1 \
  --gap2-dir /tmp/m5-openems-gap2 \
  --output-dir tests/coated-wire
```

The full transient directories are intentionally not versioned.  The compact
JSON record preserves the impedance sweeps, current samples, patterns, mesh
controls and decision metadata required by the automated test.
