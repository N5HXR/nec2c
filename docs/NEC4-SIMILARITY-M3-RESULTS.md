# NEC4 similarity M3: per-end charge factors and local junction solve

Status: accepted by the user on 2026-11-16.

Date: 2026-11-16
Branch: `codex/nec4-similarity`

## Result

M3 implements the documented NEC4 local charge solution at ordinary wire
junctions. The NEC4 engine now discovers each distinct connected endpoint
group and computes a separate complex, frequency-dependent charge factor for
every participating segment end. Radius steps and T, Y, and star junctions no
longer inherit the single NEC2 charge parameter of the segment being visited.

The solve is local: an M-member junction adds one dense `2M x 2M` system per
frequency, not new unknowns to the antenna-wide moment-method matrix. Ordinary
equal-radius two-member connections use the exact symmetry result. The classic
NEC2 engine, its real coefficient arrays, and its output remain selectable and
unchanged.

## Mathematical provenance

The implementation is independently derived from the public NEC-4 theory
manual:

- section 3.4, equation 3-42 supplies scalar-potential continuity at a
  multiwire junction;
- equations 3-43 and 3-44 define the local charge system and potential matrix;
- each ray uses the documented half-triangle charge basis and semi-infinite
  constant basis, matched at one-half and one segment length; and
- Appendix A, especially equation A-5, supplies the finite integral plus the
  outgoing half-line Hankel term.

The production solve retains the complete complex Green function and Hankel
term. It does not substitute the tempting electrostatic limit. Junction
factors are normalized by their complex mean because the common potential
scale cancels from the NEC connected-basis construction.

No NEC4 executable, source distribution, proprietary table, or undocumented
tuning constant was used.

## Implementation

`junction_charge.c` contains the local numerical operator. It uses adaptive
32-point Gauss-Legendre quadrature for the finite triangle and line integrals,
the analytic Hankel completion of the infinite line, and complex Gaussian
elimination with partial pivoting. It reports the number and maximum order of
junctions, scaled residual, pivot ratio, symmetry shortcuts, and fallbacks.
Strict model checking rejects a failed local solve; warning mode preserves a
stable legacy fallback and reports it.

The geometry connectivity is converted into endpoint groups with a union-find
topology pass. Segment direction is always oriented away from the junction,
including the sign convention encoded in NEC `ICON1`/`ICON2` links. Complex
per-end factors are carried through `sbf`, `tbf`, current reconstruction,
matrix consumers, sources, and networks. The wider coefficient storage is
needed because Appendix A produces small but real phase terms even for the
ordinary M3 fixtures.

## Independent local accuracy

`tests/nec4-similarity/m3_junctions.py` contains an independent 80-decimal
`mpmath` implementation. It evaluates its own quadrature and Hankel functions
and does not call the production local solver.

| Gate | Result | Target |
| --- | ---: | ---: |
| maximum production/oracle scaled error | `5.950e-15` | `<= 1e-11` |
| local scaled residual, four-ray star | `7.974e-17` | `<= 1e-11` |
| local permutation/reciprocity defect | `2.402e-16` | `<= 1e-11` |
| star pivot ratio | `4.236` | finite/diagnosed |

The unequal-radius straight step obtains factors
`0.3337871353 + 0.0013632223i` and
`1.6662128647 - 0.0013632223i`. The equal-radius two-member case reduces
exactly to the symmetric null result. Mirror members of the T junction agree
to rounding, and the unequal-radius four-ray star has the independently
predicted ordering.

## End-to-end behavior

The stepped-radius dipole was refined over 5, 9, 13, and 17 segments per arm.
The scaled complex-impedance errors relative to the 17-per-arm result are:

| Segments/arm | NEC2 | NEC4 M3 |
| ---: | ---: | ---: |
| 5 | `7.612e-2` | `4.954e-2` |
| 9 | `3.792e-2` | `2.611e-2` |
| 13 | `1.551e-2` | `1.142e-2` |
| 17 | reference | reference |

Both sequences are monotonic in the tested range, and the M3 sequence is
closer at every non-reference refinement. The test intentionally stops before
the current NEC2 reduced-kernel length/radius limit; completing that kernel and
its precision branches is M4.

The full-model reciprocal transfer-current check improves from `1.226e-3` in
NEC2 to `4.744e-4` in M3. It does not meet the local `1e-11` target because a
collocation voltage-gap excitation is not an exact reciprocal port observable
across unequal radii. This limitation is exposed rather than hidden; the
actual local operator's permutation reciprocity is `2.402e-16`.

The T, Y, and unequal-radius star fixtures all complete without fallback. The
historical circular center-junction crash now exercises 163 connected groups,
including 109 local solves and a maximum order of 18, with zero fallback, a
`2.153e-16` maximum scaled residual, and an `11.12` maximum pivot ratio.

## Compatibility, performance, and validation

- `--wire-engine auto` and `--wire-engine nec2` still use the preserved
  classic implementation.
- M3 remains opt-in through `--wire-engine nec4`.
- Existing cards and result schema are unchanged.
- The optimized complete suite passes: 51/51, with no skips or failures.
- An ASan/UBSan build passes M0, M2, M3, and the order-18 crash fixture; no
  sanitizer defect was reported.
- On the pinned serial 401-segment, nine-frequency real-ground fixture, 13-15
  measured runs rounded to the same `0.86 s` median for NEC2 and NEC4. A
  separate `/usr/bin/time -v` sample reported `9,728 KiB` peak RSS for both.
  M3 therefore shows no measurable ordinary-model regression at the timer's
  0.01 s resolution.

Retained evidence:

- `tests/nec4-similarity/m3-acceptance.json`
- `tests/nec4-similarity/m3-junction-convergence.png`
- `tests/nec4-similarity/m3-step-dipole.nec`
- `tests/nec4-similarity/m3-tee.nec`
- `tests/nec4-similarity/m3-y.nec`
- `tests/nec4-similarity/m3-star.nec`

No additional software or library is required.

## Manual acceptance

Run the complete M3 acceptance battery:

```sh
cd ~/homebrew/nec2c
make tests/junction_charge_test
python3 tests/nec4-similarity/m3_junctions.py
```

Expected final lines include:

```text
M3 local junction oracle maximum scaled error 5.950e-15
M3 local permutation defect 2.402e-16
M3 source-model reciprocity defect: nec2=1.226e-03, nec4=4.744e-04
M3 step convergence (NEC4 scaled to 17/arm): 5=4.954e-02, 9=2.611e-02, 13=1.142e-02, 17=0.000e+00
M3 junction solve, convergence, reciprocity, and crash gates passed
```

For direct manual inspection, run the stepped dipole and unequal-radius star:

```sh
./nec2c --wire-engine nec4 --ground-engine nec2 \
  -i tests/nec4-similarity/m3-step-dipole.nec -o /tmp/m3-step.out

./nec2c --wire-engine nec4 --ground-engine nec2 \
  -i tests/nec4-similarity/m3-star.nec -o /tmp/m3-star.out

grep -A1 "NEC4 JUNCTION CHARGE" /tmp/m3-step.out /tmp/m3-star.out
```

The step should report one non-symmetric local solve and no fallback. The star
should report 45 junction groups, one four-member local solve, 44 symmetry
shortcuts, and no fallback. Review
`tests/nec4-similarity/m3-junction-convergence.png`; the NEC4 curve should
remain below the NEC2 curve at each plotted refinement.

## Known limits

- The end-to-end convergence endpoint is the finest stable in-project result,
  not a licensed NEC4 oracle. Independent high precision covers the local
  operator; published or external complete-model witnesses remain future
  validation evidence.
- M3 deliberately stays inside the inherited wire-kernel validity envelope.
  M4 supplies the NEC4 reduced kernel and near-axis/short-segment branches.
- Free-end cap physics is unchanged until M5.
- NEC4 ground interaction is unchanged until the ground milestones.
- Junction diagnostics are in the human report and retained acceptance JSON;
  the final release milestone will consolidate validity diagnostics into the
  public machine-result contract.
