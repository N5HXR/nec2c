# NEC4ish Similarity M11E: Production Integration and Performance

## Outcome

M11E is implemented and ready for user acceptance.  It exposes the qualified
adaptive ground evaluator as an explicit production choice:

```text
--wire-engine nec4ish
--ground-engine nec4ish
--nec4ish-ground-model adaptive
```

`grid` remains the default pending acceptance.  `direct` remains the numerical
oracle.  The old `nec4` engine and option spellings are accepted only as input
compatibility aliases; capabilities, diagnostics, and reports always identify
the implementation as **NEC4ish**.  No licensed NEC source or executable was
used.

## Deterministic route hierarchy

The production evaluator never extrapolates an adaptive cell:

1. the accepted M9 near linear tables remain authoritative;
2. an accepted directional M11 cell is used inside its exact bounds;
3. the accepted ordinary M9 grid covers its qualified model volume;
4. non-grazing outer points use the published Appendix-J approximation;
5. unresolved or unqualified difficult points use the direct contour.

The retained directional regions are the M11C/M11D regions that were actually
qualified: the route-switch volume, the exact-interface line for E, and the
3-5 wavelength grazing sheet at source depth `0.02 Lg`.  A broader grazing
volume is deliberately not inferred from the sheet evidence.

Reports and machine diagnostics expose adaptive, direct, asymptotic, and
near-linear route counts.  The worker count used to construct a cold cache is
also reported; cached numerical results are independent of worker count.

## Persistent cache

The directional cache now records every k-d split axis.  Loading reconstructs
and validates the exact topology before accepting the checksum.  Cache
identity includes:

- algorithm and representation version;
- E versus H and tensor versus directional representation;
- frequency and both complex wave numbers;
- direct-contour panel, quadrature, and cutoff settings;
- bounds, error tolerance, depth limit, and cell budget.

Malformed topology, identity mismatch, checksum mismatch, and trailing data
are rejected.  Store is atomic.  E and H caches cannot be interchanged.

## Direct H production boundary

M11D's independently derived H kernels are connected for the pairing that was
actually proved: a lower-medium source and upper-medium observer.  Cylindrical
kernels are assembled into Cartesian H and integrated over the source segment
with the complex lower-medium current basis.

Same-medium and upper-to-lower H retain the legacy implementation.  M11E does
not invent unpublished transformations to make the scope appear broader.

## Retained full-model witness

The production deck contains two short horizontal wires separated by `3.1`
wavelengths over lossless water (`epsilon_r=81`).  The lower source depth is
exactly `0.02 Lg`; E and H near-field samples exercise the qualified grazing
sheet.  Adaptive and direct-oracle results agree as follows:

| Quantity | Difference from direct |
|---|---:|
| Input impedance | `2.622e-9` relative |
| All wire currents | `4.922e-5` L2 relative |
| Near E vector | `5.412e-9` L2 relative |
| Significant near H magnitudes | `0.00 dB` at report precision |
| Significant near H phases | `0.00 degree` at report precision |

The run used 288 adaptive matrix evaluations, 336 total adaptive E
evaluations, and 48 adaptive H evaluations with zero direct fallback.
One-worker and four-worker machine results were exactly equal.

## Performance

The hardest lossless-water cache is intentionally reported separately from
steady-state solving:

| Phase | Result |
|---|---:|
| Cold nine-table generation | `648.9 s` wall, `8384.6 s` CPU |
| Persistent payload | `79 MiB` |
| Warm cache startup | `0.096 s` |
| Warm adaptive end to end | `0.113 s` |
| Warm direct end to end | `0.238 s` |
| Warm end-to-end speedup | `2.10x` |
| Matrix-fill speedup | `1.48x` |
| Near-field speedup | `205.74x` |

The very large near-field speedup is specific to a witness deliberately placed
inside the qualified grazing sheet.  It must not be described as typical.

An accuracy-first control at 1,000 wavelengths on the interface lies outside
the qualified sheet.  Cache startup was only `0.101 s`, but direct E/H fallback
made the near-field phase take `230.62 s`.  This is expected and visible rather
than a silent use of the known-risk original asymptotic result.

## Acceptance

Fast retained gates:

```sh
make -j16 nec2c tests/nec4_ground_adaptive_test \
  tests/nec4_ground_adaptive_field_test tests/nec4_ground_matrix_test \
  tests/nec4_ground_runtime_test
./tests/nec4_ground_adaptive_test
./tests/nec4_ground_adaptive_field_test
./tests/nec4_ground_matrix_test
./tests/nec4_ground_runtime_test
python3 tests/nec4-similarity/m11e_production.py
python3 tests/result_json.py
```

Optional full cold regeneration takes about eleven minutes on the current
32-thread host and needs roughly 80 MiB of cache storage:

```sh
python3 tests/nec4-similarity/m11e_production.py --collect \
  --cache /tmp/nec4ish-m11e-cache
```

Evidence is retained in:

- `tests/nec4-similarity/m11e-production.nec`
- `tests/nec4-similarity/m11e-production.json`
- `tests/nec4-similarity/m11e_production.py`
