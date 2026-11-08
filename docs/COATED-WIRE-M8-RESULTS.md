# Coated-wire M8: thickness/permittivity validity envelope

## Outcome

M8 replaces the informal `t/a <= 1` expectation with a measured electrical-
radius envelope.  The useful predictor is

```
q = sqrt(ER) * k0 * b
```

where `b` is the dielectric outer radius.  Across the retained lossless,
flush-end grid, resonance error varies approximately as

```
error_percent = 8.1939*q + 0.0143     (R^2 = 0.9802)
```

There is no abrupt failure in resonance.  The strict boundary is instead the
point where resonance and feed resistance stop *both* meeting the 1% budget:

- largest tested passing `q`: 0.10908;
- first tested failing `q`: 0.12756;
- fitted 1% resonance crossing: 0.12030.

NEC2C therefore warns, but continues, when `q > 0.12`.  The warning is emitted
once per job while every coated-wire material evaluation is stamped with its
metric, threshold, and status.  This is a reference-envelope warning, not a
hard physical prohibition.

## Reference article and grid

The article is the M1/M3 split, center-fed 475 mm dipole with a 2 mm conductor
radius.  The primary references use the qualified axisymmetric OpenEMS
ring-voltage source, a 2 mm feed gap, flush dielectric ends, 1 mm axial and at
most 0.5 mm radial cells, eight angular cells, and a `1e-4` energy criterion.
The coating always has at least four radial cells.  NEC uses 95 segments and a
0.5 MHz sweep.

| ER | outer radius (mm) | `q` | resonance error | resistance error |
|---:|---:|---:|---:|---:|
| 3.2 | 3 | 0.0319 | 0.275% | 0.262% |
| 3.2 | 4 | 0.0413 | 0.378% | 0.856% |
| 3.2 | 6 | 0.0593 | 0.479% | 0.766% |
| 3.2 | 8 | 0.0765 | 0.606% | 0.570% |
| 10 | 6 | 0.1005 | 0.884% | 0.872% |
| 3.2 | 12 | 0.1091 | 0.825% | 0.476% |
| 10 | 8 | 0.1276 | 1.174% | 1.870% |
| 6 | 12 | 0.1414 | 1.152% | 3.120% |
| 10 | 12 | 0.1774 | 1.445% | 4.721% |

Resistance is the more sensitive quantity after the strict boundary.  Pattern
shape is much less sensitive: over the full grid, maximum directivity error is
0.062 dB and maximum normalized-pattern error is 0.83 dB for reference samples
above -40 dB.

At the largest strictly passing `ER=3.2`, `b=12 mm` article, normalized axial
current differs by at most 0.006 near the feed and 0.012 through the sampled
interior.  The probes one millimetre from the free metal ends differ by 0.067.
That last sample is also where the two observables are least interchangeable:
OpenEMS's contour includes dielectric displacement current around the ending
conductor, while NEC reports its wire-current basis approaching the open end.
It is retained as endpoint sensitivity, not treated as a 6.7% interior-current
failure.

The warning frequency follows directly from the metric:

```
f_warn_MHz = 0.12 * 299.8 / (2*pi*sqrt(ER)*b_metres)
```

For example, an `ER=3.2`, 12 mm-diameter housing (`b=6 mm`) reaches the warning
at about 533 MHz.  An `ER=10` housing of the same diameter reaches it near
302 MHz.

## Reference convergence

The worst primary point (`ER=10`, `b=12 mm`) was rerun with 0.5 mm axial and
0.25 mm radial resolution: 2.64 million cells and about 22 minutes.  Relative
to the primary grid, the refined OpenEMS result changed by only -0.0106% in
resonance and +0.00144% in resistance.  The retained 1.455%/4.720% NEC errors
are model drift, not reference-mesh uncertainty.

## Transfer to bends and a multiwire antenna

The straight-dipole result was challenged with a separate Cartesian OpenEMS
battery.  Each article was simulated bare and coated in both solvers.  The
primary metric is the difference between the NEC and OpenEMS *coating-induced
change*.  This pairing separates reduced-coating error from ordinary wire
staircasing, bend geometry, finite feed gaps, and multiwire-junction error.

| article | grid | `q` | resonance-shift error | resistance-change error | pattern-change error | current-change error |
|---|---:|---:|---:|---:|---:|---:|
| 45-degree dogleg, `b=6 mm` | 2 mm | 0.0496 | +0.119 point | -0.266 ohm | 0.226 dB | 0.00093 |
| 45-degree dogleg, `b=12 mm` | 2 mm | 0.0920 | -0.085 point | -0.855 ohm | 0.102 dB | 0.00125 |
| staged 90-degree route, `b=6 mm` | 1.5 mm | 0.0505 | +0.037 point | -0.498 ohm | 0.064 dB | 0.00174 |
| staged 90-degree route, 6 mm feed control | 1.5 mm | 0.0505 | -0.063 point | -0.627 ohm | 0.067 dB | 0.00102 |
| ground plane, four 45-degree down-sloped radials, `b=6 mm` | 1 mm | 0.0878 | -0.229 point | -0.733 ohm | 0.006 dB | 0.00003 |

The largest directivity-change error is 0.0112 dB.  OpenEMS's lossless NF2FF
power closure changes by at most 0.491 percentage point between a paired bare
and coated run; NEC reports unit efficiency for both.  This is treated as the
reference's numerical power uncertainty rather than dielectric absorption.

Coarse Cartesian scouting initially looked much worse.  The staged route's
resonance-shift error fell from 2.115 points on a 3 mm grid, to 1.386 points on
a 2 mm grid, to 0.037 point on the retained 1.5 mm grid.  Its pattern-change
error fell from 3.73 dB to 0.064 dB.  The 6 mm feed-gap control retained the
same conclusion, changing resonance-shift error by only 0.100 point.  The
initial apparent bend failure was therefore reference staircasing, not a
coated-wire operator failure.

The `b=12 mm` dogleg is the strongest thickness/geometry cross-test: `b/a=6`,
two qualified 45-degree bends, and `q=0.092`.  It remains comfortably inside
the straight-wire envelope.  The ground-plane article likewise supports the
coating response at `q=0.0878`, but its five-way
90-degree junction still requires `--coated-wire-junctions experimental`.
M8 does not promote that topology or relax the M7 bend constraints.

The ground-plane radials are 145 mm long and slope 45 degrees downward.  Its
Cartesian reference was repeated at 2, 1.5, and 1 mm.  The OpenEMS coating
resonance shifts were -11.554%, -12.115%, and -11.236%, respectively, versus
-11.465% in NEC.  This non-monotonic 0.879-point fine-grid movement is a
staircase-alignment uncertainty, not a converged Richardson sequence.  The
retained 1 mm comparison is the finest result, but the evidence assigns about
plus or minus 0.5 percentage point reference uncertainty to this sloped,
five-way junction.  The compact convergence record is
`tests/coated-wire/m8-groundplane45-convergence.json`.

The retained compact record and plot are
`tests/coated-wire/m8-geometry-analysis.json` and
`tests/coated-wire/m8-geometry-error.png`.  Regenerate the OpenEMS articles,
NEC comparisons, and combined record with the three
`tools/coated_wire_m8_geometry_*` scripts.

## Loss

The practical `ER=3.2`, `b=6 mm` article was tested at `tan(delta)=0.05`.
OpenEMS represents that material with constant conductivity, so NEC used the
identical center-frequency conductivity rather than a constant loss tangent.

- resonance error: 0.506%;
- resistance error: 0.581%;
- radiation-efficiency difference: -0.035 percentage point.

This intentionally large loss did not contract the observed envelope.

## Endpoint geometry is a separate uncertainty

For `ER=3.2`, `b=12 mm`, extending a flat dielectric cap one outer radius past
each metal end lowered the OpenEMS resonance from 242.447 MHz to 237.167 MHz,
a -2.178% physical shift.  The reduced operator is 0.825% high against the
flush reference but 3.069% high against the extended-flat article.

The supported interpretation is therefore a **flush metal/dielectric free-end
convention**.  Printed housings that continue beyond the conductor end can
produce more uncertainty than the reduced radial operator itself.  This is
reported separately rather than hidden by relaxing the operator budget to 5%.
M1 retains the smaller-radius flat, hemispherical, and recessed cap comparison.

## Reproduction and acceptance

The compact OpenEMS evidence and generated analysis are:

- `tests/coated-wire/m8-openems-thickness.json`;
- `tests/coated-wire/m8-variant-evidence.json`;
- `tests/coated-wire/m8-thickness-analysis.json`;
- `tests/coated-wire/m8-thickness-error.png`.

Run the retained acceptance check with:

```sh
python3 tests/coated-wire/m8_validity.py
```

The OpenEMS run directories are intentionally reduced to reproducible specs
and compact numerical results in the repository; the large transient field
files remain scratch artifacts.  `tools/coated_wire_openems.py` regenerates a
case, and the three `coated_wire_m8_*` tools regenerate the NEC maps, compact
comparison record, and figure.

M8 does not select the future VIE or layered cylindrical-kernel extension.
Those remain candidates for cases beyond the measured envelope.
