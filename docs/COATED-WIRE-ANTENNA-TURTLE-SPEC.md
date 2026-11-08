# Candidate Antenna Turtle coated-wire integration

This is a handoff specification, not an implementation commitment in nec2c.
All Antenna Turtle code, tests, documentation, and release decisions belong in
the Antenna Turtle project.  The purpose of this note is to preserve the
validated nec2c contract and suggest an idiomatic starting point for that work.

## Candidate DSL

Reuse Antenna Turtle's dielectric material definition and add one scoped wire
property form:

```racket
(define petg
  (dielectric name petg er 3.2
              loss (narrow 0.015 435)))

(with-wire-coating material petg
                    diameter (* 12 millimeters)
  body ...)
```

`with-wire-coating` applies only to wire edges created dynamically by `body`.
It does not change `wire-diameter`: the existing wire diameter remains the
metal diameter, while `diameter` is the dielectric's outside diameter.  A
scoped form is preferred over persistent `wire-coating-on`/`off` state because
it makes accidental partial coverage less likely.

Suggested first example:

```racket
#lang antenna-turtle

(model title "PETG-encapsulated dipole"
       frequency 435
       max-length (/ (model-wavelength) 40)
       wire-diameter (* 2 millimeters))

(define petg
  (dielectric name petg er 3.2
              loss (narrow 0.015 435)))

(move up (/ (model-wavelength) 4))
(feed)
(with-wire-coating material petg
                    diameter (* 12 millimeters)
  (with-push
    (with-pen
      (spread forward (/ (model-wavelength) 2.06)))))
```

The exact surface spelling is an Antenna Turtle design decision.  The
important semantic distinction is between metal diameter, coating outside
diameter, and dielectric material.

## Required semantics

1. The coating covers every selected wire edge completely.  Partial coverage
   of an individual edge is not representable.  The DSL may split a requested
   path into separately tagged coated and uncoated edges before export.
2. Coating diameter must be at least the metal diameter and normally greater.
   Both values are physical diameters; nec2c's `LD 7` field receives half the
   coating diameter as `OUTER_RADIUS`.
3. The material must have finite `ER >= 1`, nonnegative loss tangent, and
   nonnegative bulk conductivity.
4. The coating property is captured when an edge is created.  Later changes to
   wire or coating state must not retroactively change existing edges.
5. The exporter must assign NEC tags so every `LD 7` designation covers a
   complete tag.  A tag must not mix coated and uncoated edges or two different
   coating specifications.
6. Strict junction policy is the default.  Wires meeting at a node must have
   identical coating parameters.  Coated/bare and unequal-coating junctions
   may be exposed only through an explicit experimental option that maps to
   nec2c's `--coated-wire-junctions experimental` behavior.
7. Ordinary models outside `with-wire-coating` must export exactly as before.

## NEC export contract

For each fully coated tag, emit:

```text
LD 7 TAG 0 0 ER OUTER_RADIUS TAN_DELTA SIGMA_DC
```

`OUTER_RADIUS` is in meters and is derived from the DSL's outside diameter.
A constant loss tangent maps directly to `TAN_DELTA`.  Constant bulk
conductivity in S/m maps to `SIGMA_DC`; both may be present.

For a frequency-dependent material, the exporter must not silently freeze the
material at one frequency.  A future implementation should either:

* evaluate the material and execute one single-frequency NEC deck per sample,
  then combine the ordered results; or
* reject the requested sweep with a clear unsupported-dispersion diagnostic.

The first behavior is preferable when Antenna Turtle already owns sweep
orchestration.  A narrow-band constant-tangent material is sufficient for the
initial integration.

## Diagnostics and results

nec2c remains authoritative for validity and topology diagnostics.  Antenna
Turtle should preserve or display these machine-result fields:

* coating inner and outer radii;
* `ER`, effective loss tangent, and conductivity;
* `sqrt(ER)*k0*b` validity metric;
* `within-reference-envelope` or
  `warning-outside-reference-envelope` status; and
* experimental-junction status when applicable.

An optional Antenna Turtle preflight may calculate the same metric for earlier
feedback, but it must not replace the nec2c result.  The current validated
warning boundary is `sqrt(ER)*k0*b > 0.12`.

## Suggested downstream acceptance cases

1. Bare model export is byte-for-byte unchanged.
2. A lossless coated dipole exports one complete-tag `LD 7` card and reproduces
   the accepted nec2c resonance shift.
3. A lossy coated dipole produces a VSWR sweep and nonzero dielectric loss.
4. An inside-envelope model completes quietly; an outside-envelope model
   completes with the retained warning.
5. A strict unequal-coating junction fails clearly; explicit experimental mode
   runs and remains visibly marked experimental.
6. A radiation grid preserves ordered samples and exposes the coating validity
   metadata alongside the pattern result.

The downstream project should use the versioned NEC decks and expected values
under `tests/coated-wire/` as integration fixtures rather than duplicating the
OpenEMS validation campaign.
