# NEC4 similarity M10A results

Status: implementation complete on 2026-11-19; pending user acceptance.

## What the audit found

M8 and M9 agreed with one another because they shared the same final matrix
adapter, but that adapter omitted the conversion from physical `E/I` in
ohm/metre to NEC's wavelength-normalized matrix coefficient in ohms.  The
adapter now applies exactly one `lambda0` conversion.

The buried-wire current basis also still used the free-space phase constant.
M10A now propagates the lower-medium complex wave number through coefficient
construction, current reconstruction, the reduced segment kernel, cap fields,
loads, sources, and output.  The corresponding source-medium intrinsic
impedance is used by the generalized reduced kernel.

Finally, section 4.1's interface-junction endpoint point charges are restored
with the documented `exp(+j omega t)` sign convention.  Their sign and Maxwell
similarity have independent unit tests.

## Acceptance evidence

The generated evidence is in `tests/nec4-similarity/m10a-normalization-gate.json`.

- SI/NEC Maxwell similarity: exact to the reported double-precision digits.
- Lossless `epsilon_r=4` medium-kernel similarity: below `2e-12` relative.
- Invisible vacuum/vacuum interface, wholly buried dipole: impedance
  `2.243e-11`, current vector `2.256e-11` relative difference from free space.
- Invisible interface, split wire crossing the plane: impedance `1.320e-11`,
  current vector `1.324e-11` relative difference from free space.
- 14 MHz seawater witness, accelerated versus direct ground evaluation:
  impedance `1.873e-3`, current vector `1.918e-3` relative difference.
- Seawater witness, one worker versus automatic parallelism: bit-identical
  impedance and ordered current vector.
- Fixed-feed seawater refinement (`upper + lower + two fixed segments`):
  `49+9+2 = 36.880886-j37.924199 ohm`,
  `98+18+2 = 38.733416-j39.240761 ohm`, and
  `147+27+2 = 38.815827-j38.968830 ohm`.  The two finest results differ by
  about 0.18% in resistance and 0.70% in reactance.
- Classic NEC2 regression: all six decks, eleven feature classes, and fifteen
  numerical samples unchanged.

The normalization change moves the original unqualified seawater result from
`14.585-j2554.222 ohm` to the converging, passive range above.  This is a
strong dimensional, limit, and convergence result; it is not being presented
as an NEC4 oracle comparison.  Published-oracle and realistic-antenna
qualification remains part of the broader M10 case matrix.

## User acceptance

Run from the repository root:

```sh
python3 tests/nec4-similarity/m10a_normalization.py
python3 tests/regression.py
```

The first command regenerates the evidence file and prints the direct/grid,
invisible-interface, parallel-identity, and segmentation results.  Acceptance
of M10A approves the absolute-normalization foundation; it does not yet accept
the complete M10 buried-antenna feature.
