# NEC4 Similarity M11A Results: Difficult-Domain Measurement Contract

## Outcome

M11A is complete and ready for user acceptance. It made no production
behavior change. It conditioned the direct-contour oracle and measured the
accepted original-NEC4-like evaluator over 80 deliberately difficult E-tensor
witnesses.

The proposed difficult-domain scaled-error target of `1e-4` is supportable:

- the production direct contour versus the tighter contour sequence has a
  worst scaled difference of `7.314e-5`;
- exchanging Bessel and Hankel contour ownership at the tighter settings has
  a worst scaled difference of `4.593e-6` where both are available;
- the largest conditioning difference is the lossless-water interface limit,
  so that limit is classified as difficult rather than ordinary.

The current asymptotic route is not supportable as an `adaptive` fallback over
the measured domain:

- all 56 asymptotic witnesses exceed `1e-4`;
- the median scaled E-tensor error is `11.80`;
- the worst scaled error is `1763.50`, for seawater near the 3-wavelength
  route boundary;
- all 36 warning-classified grazing witnesses fail, but so do all 20
  un-warned asymptotic witnesses;
- small-rho/large-z errors range from about `0.18` to `108` in the sampled
  cases.

These are tensor-level errors normalized by the largest fine-direct component,
not percentages in an antenna result. Values greater than one mean the largest
component error exceeds the magnitude of the largest reference component at
that point. Null-sensitive antenna observables may amplify or cancel these
errors, so M11E must still qualify full models.

## Why the result is credible

Each witness uses three comparisons:

1. the accepted `3 / 48 / 80` panel, quadrature, and cutoff sequence;
2. a tighter `5 / 80 / 120` sequence;
3. the tighter sequence with Bessel/Hankel contour ownership exchanged where
   both representations are available.

The reference grid spans average soil, low-loss water, seawater, and lossless
water. It covers the 3-wavelength route transition, 1 through 20 degree
grazing angles, 4 through 10 wavelength interface ranges, exact interface
limits, and small-rho/large-z cases.

The pattern matches the public GN3 report rather than being inferred from it.
Burke identifies lower-HF seawater interpolation, interfering interface rays,
the loss of the below-ground grazing contribution, and small-rho/large-z
asymptotics as weaknesses in the earlier GN2 path, and describes adaptive
sampling and separately represented interfering contributions as the later
remedies (Burke, *Ground Model Options in the NEC-4.2 Antenna Code*, pp. 1-2).
The report does not publish enough detail to reproduce GN3, so its numerical
curves remain corroboration rather than implementation data.

## Frozen M11B decisions

1. `adaptive` retains `1e-5` as the ordinary-domain E target and freezes
   `1e-4` as the difficult-domain E target.
2. The accepted direct contour is the fallback and reference in difficult
   regions.
3. The present asymptotic expression is never selected by `adaptive` merely
   because a geometric distance threshold was crossed. It may be selected
   only in a domain where withheld direct witnesses demonstrate overlap with
   the target.
4. Initial adaptive coverage must include both sides of the current
   3-wavelength boundary, exact-interface limits, and small-rho/large-z cases.
5. The existing `grid` mode and its original-NEC4 warning behavior remain
   unchanged.
6. Direct H remains independently deferrable. Any ambiguity is resolved in
   favor of the explicit NEC4 theory manual, not inferred NEC4.2 behavior.

## Retained evidence and acceptance command

- `tests/nec4-similarity/m11a-measurement.json` contains every complex
  fine-direct E component, route classification, conditioning error, and
  timing.
- `tests/nec4-similarity/m11a-measurement.png` plots the route discontinuity
  and grazing-angle error.
- `tools/nec4_m11_measure.c` is the independent measurement executable.

Run:

```sh
make -j16 tools/nec4_m11_measure
python3 tests/nec4-similarity/m11a_measurement.py
```

To regenerate the retained evidence:

```sh
python3 tests/nec4-similarity/m11a_measurement.py --collect
```

Expected verification line:

```text
M11A retained difficult-domain oracle and production map passed
```

## M11B entry gate

M11B may begin after user acceptance of the frozen `1e-4` difficult-domain
target, direct fallback policy, and rule that the current asymptotic evaluator
is not automatically inherited by `adaptive`.
