# NEC4 similarity M1: independent oracles and operator boundary

Status: accepted by the user on 2026-11-16.

Date: 2026-11-16
Branch: `codex/nec4-similarity`

## Result

M1 establishes an implementation seam without changing wire physics. Both
`--wire-engine auto` and `--wire-engine nec2` resolve to the original NEC2
operator. Both ground selections available in M1 resolve to the original NEC2
ground implementation. Requests for the planned `nec4` or `adaptive` engines
fail explicitly; no incomplete physics can be selected accidentally.

The vtable covers basis evaluation, connected-basis assembly, electric field,
and magnetic field. Every external production call to `tbf`, `trio`, `efld`,
or `hsfld` is routed through it. Calls internal to the unchanged NEC2
implementation remain direct so that its numerical code is preserved.

The independent Python oracle imports no nec2c implementation. It evaluates:

- the cancellation-resistant `cos(x)-1` basis term with `mpmath`;
- the section 3.4 triangular-charge potential integral at 60 and 90 decimal
  digits;
- uniform flat-cap potential by direct multiprecision quadrature and its
  analytic axial limit;
- the exact, extended, and reduced constant-segment kernels by independent
  composite Gauss-Legendre surface quadrature; and
- a direct scalar Sommerfeld spectral integral against its closed Green-function
  identity.

The last check qualifies the independent infinite-interval spectral quadrature
machinery; the full two-medium electromagnetic tensor is deliberately M7 work.
The kernel grid excludes the singular surface `rho/a = 1`; M4 owns the public
singularity-extraction formulas and branch thresholds.

## Frozen numerical contract

The M1 measurements support the plan's initial production limits without
loosening them:

| Metric | Frozen limit |
| --- | ---: |
| Element/local-junction formula vs high precision | `1e-11` scaled |
| Direct Sommerfeld internal convergence | `1e-9` scaled |
| Well-conditioned linear-system residual | `1e-11` scaled |
| Reciprocity defect | `1e-11` scaled |
| Ordinary power-balance defect | `1e-6` scaled |
| Serial/parallel agreement | exact where ordered, otherwise `1e-12` scaled |
| Reference evaluator stability | `1e-13` scaled or better |

All scaled comparisons use `abs(a-b) / max(1, abs(a), abs(b))`, avoiding
meaningless relative blow-up at a true zero. Production comparisons may add a
problem-specific physical scale, but may not use a weaker limit silently.

Observed reference stability:

| Oracle | Maximum scaled delta/error |
| --- | ---: |
| Exact-kernel 48/80-point composite quadrature | `6.671e-16` |
| Junction 60/90-decimal evaluation | below JSON serialization floor |
| Cap direct/analytic and 60/90-decimal evaluation | below JSON serialization floor |
| Direct Sommerfeld/closed identity | `2.145e-30` |

The reference layer therefore has at least the required two decimal orders of
headroom relative to its corresponding production targets.

## Production consumer map

| Concern | Producers/consumers | M1 boundary consequence |
| --- | --- | --- |
| Basis topology and coefficients | `sbf`, `tbf`, `trio`, `cabc` | current reconstruction and source discontinuities dispatch through the selected operator |
| Matrix fill | `cmset`, `cmww`, `cmws`, `cmsw` | wire-wire and wire-patch field evaluations dispatch; process workers inherit selection |
| Sources and loads | `qdsrc`, matrix load terms, main excitation loop | charge-discontinuity sources dispatch; coefficient/load audit is required in M2-M5 |
| Networks | `network.c` calls `cabc` and reads A/B/C coefficients | representation migration must be explicit in M2 |
| Near fields | `nefld`, `nhfld` | per-wire E/H evaluations dispatch |
| Far fields and power | `radiation.c`, main power/current reports | these consume reconstructed A/B/C coefficients and are M2 representation consumers |
| Dielectric/coated extensions | `dielectric_wire.c`, `dielectric_ground.c`, `coated_wire.c` | they consume or temporarily replace A/B/C state; combinations stay disabled until requalified |
| Machine/human output | `result_json.c`, traditional report paths | classic output is frozen; requested/resolved engine metadata is deferred until schema-compatible integration |

## Public witness reconstruction

`m1-witnesses.json` records complete dimensionless geometry where the public
material supplies it and distinguishes equation-level references from values
that can only be digitized from plots. The strongest immediate witnesses are
the exact/extended/reduced kernel equations, the 0.5-lambda thick dipole, the
close PEC-image transmission line, the NEC4 junction-potential equations, and
the NEC4 ground-region error maps. Plot-derived points are tolerance bands, not
invented exact NEC4 values.

## Backward compatibility and performance

The M0 physics and normalized-report hashes pass under explicit
`--wire-engine nec2 --ground-engine nec2`. The CLI capabilities document now
advertises only `auto` and `nec2` as available, while separately identifying
`nec4` and `adaptive` as planned.

On `classic-frequency-sweep-801`, 15 serial repetitions measured:

| Executable | Median |
| --- | ---: |
| Archived M0 executable | 0.136136 s |
| M1 operator boundary | 0.138618 s |

The observed 1.82% difference is below the 10% regression gate and small enough
to be dominated by normal short-run host variation. Numerical results were
stable in both runs.

## Automated validation

The normal optimized build passed the complete suite: 48/48 tests, zero skips,
failures, or errors. An `-O1 -fsanitize=address,undefined` build also passed the
M1 acceptance check. LeakSanitizer itself is unavailable under the app host's
process tracing, so that run used `ASAN_OPTIONS=detect_leaks=0`; AddressSanitizer
and UndefinedBehaviorSanitizer reported no defect.

Run the milestone acceptance check:

```sh
MPLCONFIGDIR=/tmp/nec4-mpl \
  python3 tests/nec4-similarity/m1_operator.py
```

Expected final line:

```text
M1 operator boundary, independent oracles, and classic identity passed
```

Regenerate and inspect the retained oracle evidence:

```sh
MPLCONFIGDIR=/tmp/nec4-mpl \
  python3 tools/nec4_m1_oracles.py --self-test \
  --output-dir /tmp/nec4-m1
```

This writes `m1-oracles.json`, `m1-kernel-error-map.png`, and
`m1-sommerfeld-error-map.png`. No additional software or system library is
required on the current host.

## Known limits

- No NEC4 production physics exists in M1.
- The Sommerfeld reference currently qualifies scalar spectral quadrature, not
  the complete two-medium dyadic field.
- Singular on-surface kernel extraction and switching are reserved for M4.
- Junction constant-tail solution and cap field coupling are reserved for M3
  and M5; M1 validates independent integral building blocks.
- Existing coated-wire and dielectric extensions are not yet qualified with a
  future NEC4 wire representation.
