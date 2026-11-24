# NEC4 similarity M2: cancellation-resistant current representation

Status: accepted by the user on 2026-11-16.

Date: 2026-11-16
Branch: `codex/nec4-similarity`

## Result

M2 adds an explicitly selected `--wire-engine nec4ish` mode whose current on a
segment is stored and evaluated as

```text
I(s) = A + B sin(k(s-sj)) + C[cos(k(s-sj))-1].
```

This is the representation documented in section 3.2, equation 3-3 of the
public NEC-4 theory manual. The classic engine retains the NEC2/NEC3 equation
3-4 representation `A' + B sin(...) + C cos(...)` and its original operation
order. The two are algebraically identical for ordinary cases, with
`A = A' + C`, but the NEC4 form avoids repeatedly recovering a small center
current by subtracting nearly equal constant and cosine terms.

The new representation is carried through connected-basis construction,
matrix loads, voltage/current-discontinuity sources, networks, wire current
reports, near and far fields, radiation integrals, and the existing dielectric
wire/ground coupling consumers. `auto` deliberately continues to resolve to
the accepted NEC2 path during development. M2 does not add NEC4 junction,
end-cap, reduced-kernel, or ground physics.

## Mathematical and implementation provenance

The production change derives from NEC-4 theory manual section 3.2, equations
3-3 and 3-4 and the accompanying short-segment cancellation discussion. It
does not use a NEC4 executable, source distribution, or undocumented tuning
constant.

The coefficient algebra was transformed rather than computed in the old
basis and converted afterward. Small terms use cancellation-resistant
identities:

```text
cos(x)-1       = -2 sin^2(x/2)
1-cos(x)       =  2 sin^2(x/2)
cos(x)-cos(2x) =  2 sin(3x/2) sin(x/2).
```

The M1 independent `mpmath` oracle already freezes the first identity against
high precision through the extreme short-segment grid. Its maximum scaled
error is `2.647e-16`, comfortably inside the `1e-11` production target.

Representation-aware helpers are centralized in `wire_operator.c`; consumers
do not open-code assumptions about whether the third coefficient multiplies
`cos(x)` or `cos(x)-1`. The original NEC2 branches remain present and are
covered by the immutable M0 corpus.

## Numerical acceptance

The retained acceptance data are in `m2-acceptance.json`. Four end-to-end
models exercise distinct consumers:

| Model | Coverage | Maximum scaled NEC2/NEC4 delta |
| --- | --- | ---: |
| 101-segment center-fed dipole | currents, feed, radiation pattern | `1.967e-9` |
| 100-segment closed square loop | bends, closed topology, currents, pattern | `2.543e-11` |
| `NT` network fixture | networks and feed reporting | `3.647e-15` |
| 801-point real-ground sweep | loads, Sommerfeld fields, ordered sweep output | `2.848e-9` |

The comparison recursively covers every numeric value in the result document
except timing and backend identity. The frozen gate is `5e-9`. These tiny
differences are the expected rounding change from eliminating cancellation;
no later NEC4 physics is active.

For the conditioning test, a fixed 1 m, 21-segment dipole is driven at 3 MHz,
30 kHz, and 300 Hz. In the quasi-static region, ideal capacitive reactance has
approximately constant `-X f`. Taking the NEC4 3 MHz value as the reference:

| Engine | Drift at 300 Hz |
| --- | ---: |
| classic NEC2 representation | `1.360%` |
| NEC4 representation | `0.045%` |

That is a **30.4x reduction in representation-related drift**. At still lower
frequencies, errors elsewhere in the legacy field kernels eventually dominate;
those limiting branches belong to M4. M2 therefore improves conditioning
materially without claiming that it completes all low-frequency NEC4 work.

## Backward compatibility and performance

- `--wire-engine auto` and `--wire-engine nec2` still select the preserved
  classic engine.
- `--wire-engine nec4ish` is opt-in and advertised in capabilities.
- Existing NEC card syntax and result fields are unchanged.
- The immutable M0 free-space and real-ground hashes remain exact.

On the serial 801-frequency real-ground case, 15 repetitions pinned to one CPU
measured medians of `0.76 s` for NEC2 and `0.77 s` for NEC4. The observed 1.3%
difference is below the 10% regression gate and near the 0.01 s timer
resolution.

## Automated validation

The optimized build passed the complete suite: 49/49 tests, zero skips,
failures, or errors. An `-O1 -fsanitize=address,undefined` build passed both the
M2 acceptance check and immutable M0 baseline. Leak detection was disabled
because LeakSanitizer cannot operate under the app host's process tracing;
AddressSanitizer and UndefinedBehaviorSanitizer reported no defect.

No additional system software or library is required for M2.

## Manual acceptance

Run:

```sh
cd ~/homebrew/nec2c
python3 tests/nec4-similarity/m2_basis.py
```

Expected output includes:

```text
M2 ordinary maximum scaled delta 2.848e-09
M2 deep-short drift: nec2=1.360%, nec4=0.045%, improvement=30.4x
M2 representation equivalence and short-segment conditioning passed
```

Review the retained current overlay and short-segment conditioning figures:

- `tests/nec4-similarity/m2-current-equivalence.png`
- `tests/nec4-similarity/m2-short-conditioning.png`

The first should show coincident NEC2/NEC4 current magnitude on both the
straight dipole and closed loop. The second should show the NEC4 `-X f` trend
remaining close to unity as the classic representation drifts.

## Known limits

- The NEC4 mode remains development-only and explicit; `auto` does not select
  it yet.
- M2 changes representation only. Junction charge factors are M3, reduced
  wire kernels and limiting branches are M4, end caps are M5, and NEC4 ground
  physics begins in M7.
- Coated-wire and dielectric combinations are not newly claimed as NEC4-
  equivalent. M2 makes their coefficient consumers representation-safe, but
  their combined physics must be separately requalified before release.
- Extremely low-frequency failure after the demonstrated range reflects
  unmodified legacy field-kernel cancellation and is not hidden by this test.
