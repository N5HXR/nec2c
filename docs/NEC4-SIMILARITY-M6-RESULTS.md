# NEC4 Similarity M6 Results

## Outcome

M6 integrates the NEC4-similarity wire operator with NEC2C's existing
above-ground execution paths.  The free-end and source/load-gap flat-cap
fields now include the mirrored disk for perfect ground and use the same
Fresnel reflection dyadic as the existing reflection-coefficient
approximation.  The rigorous Sommerfeld path includes the direct and image
cap terms and clearly reports that its cap-specific Sommerfeld correction is
still provisional.  M8-M9 provide the public-theory interface Green-function
work needed to remove that qualification.

Unqualified combinations no longer mix operators silently.  Explicit NEC4
wire mode rejects LD 7 coated wires and dielectric VIE models.  Wire/patch
and loaded/network models remain enabled and are covered by M6 tests.

## Numerical gates

`tests/nec4-similarity/m6_combined.py` records the detailed evidence in
`tests/nec4-similarity/m6-acceptance.json`.

- PEC cap-image boundary: zero tangential electric-field defect at the five
  sampled ground-plane points.
- Reciprocity: scaled defects of `3.39e-20` (PEC), `1.52e-20` (RCA), and
  `4.29e-20` (Sommerfeld).
- Frequency concurrency: one-worker and two-worker result documents are
  exactly equal after timing fields are removed for all three ground modes.
- Ordinary thin-wire comparison: NEC2-to-NEC4 feed-impedance changes are
  `0.508%` (PEC), `0.141%` (RCA), and `0.509%` (Sommerfeld) for the frozen
  elevated-dipole witness.
- Wire/patch: 201 wire unknowns and four patch-current rows solve with finite
  feed impedance.
- Loaded/network: one-thread and eight-thread numerical delta is zero.
- All sweep executions produce ordered far fields and positive accepted
  power.
- The optimized regression suite passes 57/57, the M0-M6 milestone suite
  passes 9/9, and an out-of-tree ASan/UBSan build passes the complete M6
  corpus with LeakSanitizer disabled because it cannot run under the app host.

Adding the PEC cap image necessarily updates the M5 quarter-wave-monopole
artifact: that fixture uses `GN 1`, despite M5 otherwise being described as a
free-space preview.  The regenerated M5 evidence and plot are retained as the
current combined-operator result; the originally accepted pre-image figures
remain recorded in the M5 report and commit history.

## Performance gate

The nine-frequency, 401-segment real-ground fixture was run seven measured
times after two warmups with the builtin dense backend:

| Mode | Median wall time | Physics hash |
| --- | ---: | --- |
| `--jobs 1 --threads 1` | 1.884297 s | `34f754b7...31e1` |
| `--jobs auto --threads auto` | 0.346744 s | `34f754b7...31e1` |

Automatic concurrency is **5.43x faster** with an identical physics hash.
The accepted real-ground parallel speedup therefore remains effective under
the NEC4 wire operator.

## User acceptance

The supplied `tests/nec4-similarity/m6-user-above-ground.nec` is a 10.2 m,
101-segment horizontal dipole 10 m above `Er=13`, `sigma=0.005 S/m` ground.
It requests an 11-point 13.5-14.5 MHz impedance sweep and a 5-degree
hemisphere radiation grid at each frequency.

Run:

```sh
./nec2c --wire-engine nec4 --ground-engine nec2 \
  --backend auto --jobs auto --threads auto \
  --result-json /tmp/m6-user.json \
  -i tests/nec4-similarity/m6-user-above-ground.nec \
  -o /tmp/m6-user.out
```

Acceptance review:

1. Confirm all 11 frequencies complete and the report contains the explicit
   provisional Sommerfeld cap notice.
2. Review the feed resistance/reactance or compute VSWR for the intended feed
   impedance across 13.5-14.5 MHz.
3. Inspect the ordered 19-by-73 far-field grids and confirm the expected
   above-ground dipole pattern evolves smoothly across the sweep.
4. Optionally repeat with `--jobs 1 --threads 1` and require identical
   numerical result content apart from execution/timing metadata.

## Scope retained for later milestones

- The cap-specific rigorous Sommerfeld correction is not claimed complete.
- Radial-screen cap reflection uses the homogeneous-ground coefficient and
  emits a notice.
- Coated-wire/NEC4 and dielectric-VIE/NEC4 combinations require their own
  operator derivations and validation before enablement.
- Buried and interface-crossing wires remain M7-M10 work.
