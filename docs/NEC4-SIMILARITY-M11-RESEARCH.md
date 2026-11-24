# NEC4 Similarity M11 Research Note: Adaptive Ground Evaluation

## Decision sought

M11 should build an independently specified ground evaluator named
`adaptive`.  It is informed by the public NEC-4.2 GN3 report, but it must not
be called GN3 or claim numerical equivalence to GN3.  The report describes
the architecture and validation cases but omits the adaptive tolerances,
region construction, branch-cut formulas, direct-H formulas, and higher-order
asymptotic formulas.

The proposed implementation obtains accuracy from public NEC4 equations,
independent Maxwell derivations, and the accepted M8 direct contour oracle.
It does not fit an offset to published plots or licensed-program output.

## Why M11 is justified

The accepted original-NEC4-like engine has a known boundary:

- its below-ground first-order grazing contribution loses two-ray
  interference;
- fixed interpolation regions can become inaccurate for seawater at lower HF;
- numerical curl magnifies small E-field errors in near H;
- the `10-j1000` reference-dipole witness is `5.596 dB` below the published
  absolute field, although this common factor cancels correctly in RCE.

The public GN3 report independently identifies these same classes of failure.
It says the later method uses adaptive cubic interpolation, separately
interpolated branch-cut terms where interference occurs, direct E and H, and
higher-order asymptotics.  Those statements are requirements and comparison
evidence, not a complete implementation specification.

## Proposed capability contract

The new CLI selection is:

```text
--nec4-ground-model adaptive
```

During M11 it remains opt-in.  `grid` retains accepted original-NEC4 behavior,
including its grazing warning.  M12 may select `adaptive` automatically only
inside the domain qualified here and only after explicit user approval.

`adaptive` means:

1. Direct-contour samples are generated until a declared interpolation error
   target is met; sample density is not fixed by material category.
2. Total-field interpolation is used only where its independently checked
   error is acceptable through the cell.
3. Interfering contributions are represented separately where total-field
   phase removal crosses or approaches a null.
4. The direct contour remains the fallback and oracle when no accelerated
   representation passes its error test.
5. An asymptotic path is used only where it overlaps the direct oracle within
   the frozen measurement contract.  Extending cached direct/adaptive tables
   is preferred to using an unqualified asymptotic expression.
6. E and H receive independent validation.  H is derived directly from the
   spectral Green function/Maxwell curl; it is not inferred by numerically
   differencing an interpolated E table.

## Milestone decomposition

### M11A - difficult-domain measurement contract

- Freeze a dimensionless domain spanning ordinary soil, low-loss dielectric,
  seawater/lower-HF, interface limits, two-ray nulls, small-rho/large-z, and
  the current direct/asymptotic switch.
- Establish converged direct-contour E values with two independent contour or
  precision sequences.
- Map current `grid` error and discontinuity without changing production
  routing.
- Reproduce the published 5 m seawater-whip geometry and short-dipole
  interference geometry as completely as the report permits.

### M11B - adaptive total-field E tables

- Implement hierarchical cubic cells with deterministic refinement.
- Use midpoint/face/volume witnesses excluded from the interpolation fit to
  estimate cell error.
- Retain exact interface and quasistatic subtractions already accepted in M8
  and M9.
- Add versioned persistent caching including tolerance, bounds, material,
  frequency, representation, and algorithm version.

### M11C - interference-safe representation

- Expose independently derived contour contributions before they are summed.
- Detect cells whose total field approaches a null or fails phase-smoothness
  tests.
- Compare total-field and separated representations against withheld direct
  samples; select only a passing representation.
- Extend adaptive/direct coverage through grazing ranges where the current
  asymptotic path does not overlap the oracle.

### M11D - direct magnetic field

- Derive the magnetic spectral kernels from Maxwell's equations using the
  same time and passive-medium conventions as the E evaluator.
- Qualify homogeneous, PEC, reciprocity/duality, interface, and numerical-curl
  cross-checks before production use.
- Add adaptive H tables only after direct H is independently accepted.

### M11E - production integration and performance

Implementation and retained results are documented in
`NEC4-SIMILARITY-M11E-RESULTS.md`.  User acceptance passed on 2026-11-23.

- Add deterministic routing, diagnostics, cache lifecycle, and worker-count
  identity.
- Compare uncached generation, cached startup, matrix fill, near-field, and
  full-model time separately.
- Preserve `grid` and direct-oracle choices as explicit fallbacks.

## Candidate accuracy contract

M11A measures what the accepted oracle can support before these numbers are
frozen.  The initial targets are deliberately accuracy-first:

| Quantity | Ordinary domain | Difficult/interference domain |
|---|---:|---:|
| Complex E tensor, scaled error | `1e-5` | `1e-4` |
| Route-boundary discontinuity, scaled | `1e-5` | `1e-4` |
| Direct H tensor, scaled error | `1e-5` | `3e-4` |
| Full-antenna field magnitude vs direct | `0.02 dB` | `0.05 dB` |
| Full-antenna field phase vs direct | `0.2 degree` | `0.5 degree` |
| Input impedance vs direct | `0.05%` | `0.1%` |

Scaled error uses a component-family scale from the direct field and analytic
singular terms, not pointwise relative error at a physical null.  Every table
cell also has an absolute floor tied to its dimensional field scale.

The targets may be relaxed only from measured oracle conditioning, with the
reason and affected region documented before implementation proceeds.

## Acceptance witnesses

- Public GN3-report 5 m whip: 5 MHz, 10-segment 5 m monopole, 0.2 m/6-segment
  ground stake, 5 mm wire radius, `epsilon_r=81`, `sigma=4 S/m`; the published
  GN2/GN3 impedances are corroboration because the complete deck is absent.
- The same geometry with `sigma=0.01 S/m`.
- Short vertical dipole from `0.002` to `0.012 lambda` over lossless water,
  with interface-range sweeps through and beyond the old 2/3-lambda switch.
- The M10B absolute interface field and M10D RCE controls.
- Small-rho/large-z, direct/asymptotic overlap, two-ray null, E/H curl, PEC,
  homogeneous-medium, reciprocity, and worker-identity grids.

## Stop and deferral rules

- If a direct-contour reference cannot be converged independently, that region
  is not accelerated and remains warning-backed direct-only or unsupported.
- If separated contour terms cannot be derived from public equations, M11 may
  extend direct/adaptive total-field coverage but will not imitate GN3's
  undocumented split.
- If direct H cannot meet Maxwell/PEC/homogeneous checks, E improvements may be
  accepted independently while H is deferred.
- No empirical material-specific correction is allowed.
- Failure to deliver M11 does not invalidate accepted M10; it produces a
  deferral report and M12 proceeds with the original-NEC4-like engine.

## User decisions before implementation

1. Approve the public name `adaptive` and the explicit non-GN3 claim.
2. Approve the accuracy-first direct fallback, even when a cold run becomes
   substantially slower.
3. Approve M11A as a measurement milestone that may revise numerical targets
   before M11B implementation.
4. Decide whether direct H is required for M11 completion or may be accepted as
   an independently deferrable sub-feature.

## Approved decisions

Approved by the user on 2026-11-23:

1. The public name is `adaptive`, with no GN3-equivalence claim.
2. Accuracy-first direct fallback is preferred even when cold execution is
   substantially slower.
3. M11A may revise the provisional targets from measured oracle conditioning
   before M11B begins.
4. Direct H is best-effort and independently deferrable.  When H behavior is
   ambiguous, the explicit 1992 NEC4 theory description is the governing
   baseline; inferred NEC4.2/GN3 behavior is not imported into this feature.
5. A licensed NEC4.2 numerical oracle is neither available nor required for
   this project.  Any later access to one belongs to a separate validation or
   refinement effort.
