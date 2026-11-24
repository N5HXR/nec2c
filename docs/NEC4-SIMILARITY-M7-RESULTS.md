# NEC4 similarity M7: homogeneous lower medium and interface topology

Status: accepted by the user on 2026-11-17.

## Scope and safety boundary

M7 adds the material and topology foundation needed by the buried-wire field
work in M8 and M9. It does **not** substitute the classic NEC2 ground field
for a buried or penetrating-wire answer. In this milestone,
`--ground-engine nec4ish` therefore requires all of the following:

- a finite Sommerfeld `GN 2` lower medium;
- `--validate-only --result-json FILE`;
- no segment that straddles the interface without ending at `z=0`;
- no complete wire segment in the interface plane; and
- no surface patches.

A production calculation requested with this backend fails explicitly with
`NEC4_M7_REFERENCE_ONLY`. Production direct interface fields begin in M8 and
the accelerated documented ground evaluator begins in M9.

## Public-theory implementation

The implementation follows the public theory guide's Chapter 3 formulation:

- the wire basis uses the wave number of the medium containing the wire
  (discussion following equation 3-3, PDF page 17);
- current is continuous across the planar interface;
- the derivative/current-charge parameter changes by the complex-permittivity
  ratio, with each side's basis factor set to the permittivity of its medium
  (equation 3-45 and surrounding text, PDF page 27); and
- crossings beyond 60 degrees from the interface normal are counted as audit
  metadata. They are not rejected or warned about: the guide reports that
  equation 3-45 had been used there with apparently good results.

For the lower medium, with the program's `exp(+j omega t)` convention,

```text
epsilon_r_complex = epsilon_r - j sigma/(omega epsilon_0)
n                 = sqrt(epsilon_r_complex)
k                 = (omega/c_0) n
eta               = eta_0/n
```

Negative NEC `GN` conductivity retains the classic meaning of directly
specifying the negative imaginary part of relative permittivity.

## Geometry contract

Classification uses a scale-aware interface tolerance:

```text
256 DBL_EPSILON max(1 m, maximum absolute model coordinate)
```

Endpoints within this tolerance of `z=0` are on the interface. All segments
must otherwise be wholly above or wholly below it. Coincident interface
endpoints are clustered with the same tolerance, making a crossing junction
only when the cluster contains at least one segment in each medium.

The human report and result JSON record segment counts, interface endpoints,
crossing junctions, highly oblique crossings, tolerance, complex permittivity,
refractive index, wave number, and intrinsic impedance.

## Automated results

Run:

```sh
make check TESTS='tests/nec4_medium_test tests/nec4-similarity/m7_medium.py'
```

The acceptance harness reports:

```text
M7 homogeneous medium: max analytic error 2.119e-16
M7 topology: buried 41 lower segments; stake 21 upper + 21 lower,
one crossing; unsplit/in-plane rejected; oblique audit passed
M7 interface condition and homogeneous reciprocity residuals are exact;
production field solve is explicitly deferred to M8/M9
```

The C unit test directly checks the crossing basis-factor ratio, scale-aware
endpoint clustering, unsplit and in-plane rejection, and both conductivity and
complex-permittivity loss inputs. The Python test independently reconstructs
the analytic material constants, vacuum null, lossless `epsilon_r=4` wavelength
and impedance scaling, homogeneous Green-function reciprocity, input topology,
machine diagnostics, and the production-solve safety gate.

The complete optimized repository suite passed `60/60`. A separate
`-O1 -fsanitize=address,undefined` build passed the M7 C unit test, the full M7
acceptance harness, and the immutable M0 baseline. Leak detection was disabled
because LeakSanitizer cannot operate under the app host's process tracing;
AddressSanitizer and UndefinedBehaviorSanitizer reported no defect. No new
system package or library is required.

## User acceptance

First build, then run the two supplied models through preflight:

```sh
./nec2c --wire-engine nec4ish --ground-engine nec4ish --validate-only \
  --result-json /tmp/m7-buried.json \
  -i tests/nec4-similarity/m7-buried-dipole.nec \
  -o /tmp/m7-buried.out

./nec2c --wire-engine nec4ish --ground-engine nec4ish --validate-only \
  --result-json /tmp/m7-stake.json \
  -i tests/nec4-similarity/m7-ground-stake.nec \
  -o /tmp/m7-stake.out
```

The first report should contain `0 UPPER, 41 LOWER`; the second should contain
`21 UPPER, 21 LOWER` and `1 CROSSING JUNCTIONS`.

The independent homogeneous reference can be inspected separately:

```sh
python3 tools/nec4_m7_reference.py \
  --frequency-mhz 7.1 --er 13 --sigma 0.005 \
  --distance 0.1 --distance 1 --distance 10 \
  --output /tmp/m7-reference.json
```

Expected for this deliberately lossy soil example: complex relative
permittivity `13 - j12.6585236511`, complex wave number
`0.5872136539 - j0.2386667641 1/m`, and intrinsic impedance
`81.93209786 + j33.30043254 ohm`.

User acceptance: **PASS** on 2026-11-17.
