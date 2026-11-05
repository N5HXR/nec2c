# M3 lossless isolated-dipole operator

Status: accepted 2026-11-05

Date: 2026-11-05

## Delivered operator

M3 enables ordinary execution of `LD 7` for one lossless, isolated, straight,
collinear coating domain in free space.  It retains the original NEC current
unknowns and matrix dimensions.  Loss, ground, bends, surface patches,
multiple/disconnected coating domains, and coated-to-bare transitions remain
explicit errors.

For each source current piece, M3 adds only the scalar-potential correction
derived in M2:

```text
alpha = 1 - 1/ER
Delta K = alpha * (G_b - G_a)
```

The correction enters `cmww()` after the ordinary thin-wire field is projected
on the observation tangent.  The derivative of NEC's sine/cosine current basis
is integrated directly, so no empirical series impedance or velocity factor is
introduced.  The paired inner/outer kernel difference is evaluated at each
quadrature point before accumulation.  Adaptive 8/16-point Gauss quadrature
resolves self and adjacent terms without separately subtracting two completed
near-singular integrals.  The constant basis contributes no coating charge.

The normalization follows NEC's wavelength coordinates and existing
`CONST3 = -j*eta0/(4*pi)`.  It is the direct specialization of the two-potential
equation in Popovic, Djordjevic, and Kircanski (1981), with the vector-potential
term unchanged and the charge term split between radii `a` and `b`.

## Numerical invariants

`tests/coated-wire/m3_operator.py` verifies:

* exact execution-result identity with the bare wire for `ER=1`;
* exact identity when coating and conductor radii are equal;
* both null limits under NEC's ordinary and extended thin-wire kernels;
* stable first-order behavior as `b/a` approaches one;
* exact serial/parallel matrix-fill identity;
* built-in/LAPACK agreement within the repository's backend tolerances;
* exact reciprocal transfer current between mirrored test segments;
* convergence over 47, 71, and 95 segments; and
* stable diagnostics for loss, bends, and ground.

The immutable M0 bare-wire witness remains unchanged.  The uncoated 801-point
performance witness measured 0.131 seconds median in a seven-run local check;
inactive coating lookup occurs once per source segment, outside the observation
loop.

## OpenEMS comparison

The common article is a 475 mm dipole with 2 mm conductor radius.  The coated
case has 4 mm outer radius and `ER=3.2`.  NEC uses its continuous center-segment
source; OpenEMS uses the M1-qualified 2 mm split ring-voltage source.  M1 found
the 1/2 mm gap sensitivity negligible relative to this comparison.

| Quantity | NEC M3 | OpenEMS M1 | Difference |
| --- | ---: | ---: | ---: |
| Bare resonance | 295.952709 MHz | 295.815987 MHz | +0.0462% |
| Bare resistance at resonance | 72.260167 ohm | 72.797845 ohm | -0.537678 ohm |
| Coated resonance | 276.642354 MHz | 275.599687 MHz | +0.3783% |
| Coated resistance at resonance | 64.632156 ohm | 64.083814 ohm | +0.548342 ohm |
| Resonance shift | -6.5248% | -6.8341% | +0.3093 percentage point |
| Coated sampled directivity | 2.093163 dBi | 2.164174 dBi | -0.071011 dB |

The normalized coated pattern differs from the refined OpenEMS result by at
most 0.11345 dB above -40 dB.  Most of that maximum occurs at the first
five-degree sample near the axial null; the broadside region is substantially
closer.  These results meet M3's promotion target of the OpenEMS uncertainty
plus 1% in resonance and 2 ohms per impedance component.

The 47/71/95-segment coated resonance sequence is 276.124416, 276.424610, and
276.642354 MHz; resistance is 64.342710, 64.492366, and 64.632156 ohm.  Finer
segmentation is not promoted for this particular 4 mm-diameter conductor
because segment length then approaches or falls below conductor diameter,
outside the ordinary NEC thin-wire discretization regime.

## Acceptance evidence

The versioned machine record is `tests/coated-wire/m3-acceptance.json`.
The generated overlays are:

* `tests/coated-wire/m3-impedance.png`
* `tests/coated-wire/m3-current.png`
* `tests/coated-wire/m3-pattern.png`

Run the automated acceptance check:

```sh
python3 tests/coated-wire/m3_operator.py
```

For manual comparison, run the four supplied `m3-{bare,coated}-{dipole,sweep}`
decks using the commands in `tests/coated-wire/README.md`.  M3 is complete only
after explicit user review and a PASS verdict.
