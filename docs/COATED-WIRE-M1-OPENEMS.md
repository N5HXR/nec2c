# Coated-wire M1 OpenEMS reference qualification

Status: M1 implementation complete; awaiting user acceptance.

Date: 2026-11-03  
Branch: `codex/coated-wire`

## Purpose

M1 establishes an independent FDTD reference for the first coated-wire article:
a straight, centre-fed, 475 mm dipole made from a 2 mm radius PEC conductor,
optionally in direct contact with a 4 mm outer-radius dielectric sleeve with
relative permittivity 3.2.  This milestone qualifies the reference calculation;
it does not yet add a coated-wire operator to NEC2.

## Reference construction

The primary model uses OpenEMS 0.37.0 and CSXCAD 0.7.0 in cylindrical
coordinates.  The radial coordinate starts at the axis, so there is no hollow
core or artificial inner boundary.  The PEC conductor occupies the core and a
finite-radius, z-directed ring voltage source spans the 2 mm centre gap.  This
avoids placing a lumped source on the cylindrical-coordinate singularity.

The source is uniform around azimuth.  Measured cylindrical near-field
components were axisymmetric to approximately `1e-7` relative RMS.  OpenEMS's
generic cylindrical NF2FF quadrature nevertheless showed an azimuth-dependent
far-field error when only eight physical azimuth cells were used.  The harness
therefore averages the measured cylindrical near field over azimuth, rejects a
non-axisymmetric result above `1e-5` relative RMS, and evaluates NF2FF with 512
azimuth samples.  Eight- and sixteen-cell physical simulations then agreed to
`2.5e-5 dB` in directivity and `1.6e-5 dB` in sampled pattern magnitude.

The selected Huygens surface is aligned with Yee-grid planes at 500 mm radius
and +/-600 mm axial extent.  Moving it from 250/400 mm to 500/600 mm changed
directivity by 0.036 dB and all reported pattern samples by less than 0.037 dB;
the earlier 100/300 mm surface was measurably too close.  The domain extends to
3200 mm radius and +/-3300 mm, terminated by eight-cell PML.

The dielectric thickness is resolved by at least four radial cells in every
screening calculation.  The selected endpoint contrast is repeated with eight
cells through the dielectric so that a claimed coated-wire validity limit is
not actually an FDTD through-thickness error.

## Geometry conventions under audit

The proposed reduced model's normative convention is a solid concentric sleeve
with flat dielectric faces flush with the conductor free ends and feed-gap
faces.  M1 compares four physical free-end geometries while keeping the feed
faces flush:

- `flush`: dielectric and conductor terminate on the same plane.
- `flat-overhang`: a solid flat dielectric cap extends one outer radius beyond
  the conductor.
- `hemisphere`: a solid hemispherical dielectric cap covers the free end.
- `recessed`: the conductor protrudes one outer radius beyond the dielectric.

These are OpenEMS reference geometries, not proposed NEC card variants.  If
their effect exceeds the reduced model's error budget, a later milestone must
decide whether endpoint corrections are needed.

## Convergence results

The initial cylindrical mesh sequence, before the final aligned-surface audit,
gave:

| Article | Coarse resonance | Fine resonance | Change | Coarse R at resonance | Fine R at resonance | Change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Bare | 295.659992 MHz | 295.815987 MHz | 0.05276% | 72.786263 ohm | 72.797845 ohm | 0.01158 ohm |
| Coated, flush | 275.979089 MHz | 276.061659 MHz | 0.02992% | 64.070819 ohm | 64.086058 ohm | 0.01524 ohm |

Both are within the M1 limits of 0.25% resonance change and 0.5 ohm change in
each impedance component.  Over the common-frequency samples for which
`|X| < 25 ohm`, the bare sequence changed by at most 0.117 ohm in resistance
and 0.409 ohm in reactance; the coated sequence changed by at most 0.075 ohm
and 0.235 ohm respectively.  The aligned-surface results, endpoint sweep,
feed/gap audit, decay audit, and eight-cell coating confirmation are recorded
in the versioned summaries described below.

The four-cell endpoint screening sweep produced:

| Free-end geometry | Resonance | R at resonance | Directivity | Resonance shift from flush | Max sampled-pattern change from flush |
| --- | ---: | ---: | ---: | ---: | ---: |
| Flush | 275.576650 MHz | 63.992099 ohm | 2.165759 dBi | reference | reference |
| Flat overhang | 273.523825 MHz | 63.947442 ohm | 2.167549 dBi | -2.052824 MHz (-0.745%) | 0.00270 dB |
| Hemisphere | 274.539915 MHz | 63.983400 ohm | 2.167295 dBi | -1.036735 MHz (-0.376%) | 0.00467 dB |
| Protruding conductor | 277.033327 MHz | 64.058827 ohm | 2.163385 dBi | +1.456678 MHz (+0.529%) | 0.00247 dB |

The end treatment is therefore physically material to resonance even though it
barely changes resistance or the normalized pattern for this article.  On the
eight-cell mesh, flush converged by 0.0084% in resonance, 0.092 ohm in
resistance, and 0.0024 dB in sampled pattern; flat overhang converged by 0.0008%,
0.075 ohm, and 0.0041 dB.  The refined overhang shift from refined flush is
-2.077938 MHz (-0.754%), while its normalized pattern differs by only 0.0011 dB
above -40 dB.  This promotes the flat-end contrast beyond the numerical error
budget.

The hemispherical result was also repeated at eight-cell resolution because
even the cylindrical grid staircases curvature in its `rho,z` meridian; it
cannot inherit the flat-boundary convergence result.  It converged by 0.0321%
in resonance, 0.088 ohm in resistance, and 0.00295 dB in sampled pattern.  The
refined hemisphere lies -0.9716 MHz (-0.353%) below refined flush.  Its measured
non-axisymmetric near-field energy ratio was `5.9e-6`, below the `1e-5` guard
(the flat cases were approximately `1e-7`).  This separation is
physically consistent with dielectric end loading: the cap perturbs the strong
open-end electric field and effective capacitance while leaving the article's
large-scale current distribution and radiation shape nearly unchanged.

### Feed, gap, and transient-decay audits

Reducing the physical centre gap from 2 mm to 1 mm changed resonance by
-0.006668 MHz (-0.0024%), resistance at resonance by -0.01065 ohm,
directivity by +0.00018 dB, and every sampled pattern value above -40 dB by at
most 0.00010 dB.  A 2 mm gap is therefore retained as the reference convention,
with gap sensitivity comfortably below the numerical budget.

Tightening OpenEMS's energy end criterion from `1e-4` (-40 dB) to `1e-5`
(-50 dB) changed resonance by +0.052695 MHz (+0.0191%), resistance at resonance
by -0.00047 ohm, directivity by -0.00179 dB, and sampled pattern by at most
0.01143 dB above -40 dB.  The less expensive `1e-4` criterion meets the M1
budgets; the tighter result quantifies its residual transient uncertainty.

The continuous-conductor (`embedded`) excitation experiment produced zero
terminal voltage and impedance: the PEC directly shorts the impressed axial
ring field, making the port unobservable.  Its far field is numerical residue,
not a valid antenna solution.  The harness now rejects this degenerate response
instead of reporting the first sweep frequency as a false resonance.  The
qualified reference consequently requires a physical split between the two
dipole arms; this is an excitation constraint, not evidence that coated and
uncoated conductors cannot meet in the eventual NEC model.

## Reproduction

Fast dependency-free harness checks:

```sh
python3 tests/coated-wire/m1_harness.py
```

A representative four-cell, flush-coated reference run is:

```sh
/home/josh/opt/openEMS/venv/bin/python tools/coated_wire_openems.py \
  --mode run --coordinates cylindrical --material coated --feed split \
  --cylindrical-source ring-voltage \
  --mesh-mm 1 --radial-mesh-mm 0.5 \
  --domain-radius-mm 3200 --domain-half-height-mm 3300 \
  --nf2ff-radius-mm 500 --nf2ff-half-height-mm 600 \
  --air-mesh-mm 40 --angular-cells 8 \
  --axisymmetric-nf2ff-samples 512 \
  --end-criteria 1e-4 --timesteps 720000 \
  --frequency-count 101 --threads 16 \
  --coating-end flush --output-dir /tmp/coated-wire-m1-manual
```

The repository retains compact JSON results and figures in
`tests/coated-wire/m1-openems-{caps,audits}*`.  Multi-gigabyte transient and
near-field dumps remain disposable `/tmp` artifacts.

## Result and acceptance summary

The selected coated reference is the split 2 mm gap, flat/flush sleeve,
eight-cell through-thickness mesh.  It gives 275.599687 MHz resonance,
64.083814 ohm resistance at resonance, and 2.164174 dBi directivity.  The bare
fine-mesh impedance reference gives 295.815987 MHz and 72.797845 ohm; the
aligned-surface bare pattern gives 2.151491 dBi.  The coating therefore lowers
this article's resonance by about 6.83%; deriving and validating the reduced
NEC operator that predicts that shift belongs to M2/M3.

| Qualification item | Worst relevant result | M1 limit | Result |
| --- | ---: | ---: | --- |
| Bare mesh resonance | 0.0528% | 0.25% | PASS |
| Coated flush mesh resonance | 0.0084% | 0.25% | PASS |
| Bare near-resonance R / X | 0.117 / 0.409 ohm | 0.5 / 0.5 ohm | PASS |
| Coated near-resonance R / X | 0.075 / 0.235 ohm | 0.5 / 0.5 ohm | PASS |
| Flush sampled pattern refinement | 0.0024 dB | 0.05 dB | PASS |
| Hemisphere sampled pattern refinement | 0.0030 dB | 0.05 dB | PASS |
| Huygens-surface sampled pattern | 0.037 dB | 0.05 dB | PASS |
| Physical 8/16-alpha-cell agreement after m=0 projection | 0.000016 dB | 0.05 dB | PASS |
| `1e-4` / `1e-5` decay criterion pattern | 0.0114 dB | 0.05 dB | PASS |
| Near-field non-axisymmetric energy | `5.9e-6` | `1e-5` | PASS |

The continuous-conductor source is explicitly rejected, not promoted.  Free-end
geometry is a physical model parameter outside the initial NEC card contract:
its approximately 0.35--0.75% resonance effects are much larger than the
reference uncertainty and must remain visible in M3/M7 comparisons.

## Promotion gate

The OpenEMS model is promoted as the M2/M3 validation reference only when:

- successive refinement changes resonance by less than 0.25%;
- each feed-impedance component near resonance changes by less than 0.5 ohm;
- sampled far-field magnitude changes by less than 0.05 dB;
- four- and eight-cell coating discretizations agree within those limits;
- feed construction, feed gap, decay criterion, endpoint convention, Huygens
  surface, and domain sensitivities are reported separately; and
- the user reviews the result package, runs the manual case, and records PASS.
