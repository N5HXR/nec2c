# NEC4 Similarity M10B: One-Sided Interface Endpoint and Surface-Wave Reference

## Scope

M10B qualifies the smallest published realistic antenna in the 1984
Burke/King/Miller surface-wave report before M10 proceeds to buried dipoles and
ground stakes.  The model is a centre-fed vertical dipole, 0.1 free-space
wavelength long, whose lower end lies on the air/ground interface without
electrical contact.  Its field is sampled on the interface 1000 wavelengths
away.

This sub-component adds no NEC card syntax.  It makes the already documented
free-wire-end condition legal in NEC4-similarity topology and repairs the
field-evaluation adapter so arbitrary NE/NH observation points are not treated
as matrix observation segments.

## Public reference and model uncertainty

Table 1 of `reports/surface-waves.pdf` publishes
`|E_reference|^2 / P_input` for relative ground permittivities `10-j1`,
`10-j10`, `10-j100`, `10-j1000`, and perfect ground.  The report does not state
the reference dipole radius.  It does state `1e-4 lambda` for the related
monopole/ground-stake models, so M10B uses that radius and records it as an
inference rather than a fact.  Thirty-one equal segments give a minimum
segment-length/radius ratio above 32.

## Frozen behavior

1. A unique interface point with exactly one incident upper- or lower-medium
   wire endpoint is a supported free end.  It is not an electrical crossing,
   gets no interface charge-ratio factor, and retains the ordinary free-end
   current/cap treatment.
2. Branches, unsplit crossings, and complete segments in the interface plane
   retain their existing errors.
3. Arbitrary field points infer their medium from physical `z`, using the
   interface tolerance and the upper-medium limiting value at `z=0`.
4. Matrix observations continue to use the classified segment medium.  The
   classic point-field path retains the historical `ij=1` basis/end selector.
   The NEC4 point-field path evaluates the complete interface segment field
   directly; it must not subtract one direct-field implementation and add a
   second one, since that destroys direct/image cancellation at long range.
5. Entering the accepted M9 asymptotic path at 15 degrees or less sets one
   deduplicated `NEC4_M9_GRAZING_INTERFERENCE` warning, even when the first
   risky interaction occurs during a deferred NE/NH request rather than matrix
   fill.
6. M10B does not silently replace the original-NEC4 grazing formulation with
   a later NEC4.2 correction.  A grazing warning does not relax numerical
   acceptance where a published original-NEC4 reference exists.

## Acceptance gates

- Unit topology: the one-sided endpoint is accepted and counted; a normal
  two-sided crossing remains unchanged.
- Published control: the classic exact-ground path reproduces each Table 1
  value within 16 percent under the documented radius inference.  The looser
  bound covers the report's omitted radius/discretization details; four of the
  five results must remain within 5 percent.
- NEC4 solution: reference input impedance differs from the classic solution
  by less than 0.1 percent.
- Surface-field accuracy: the NEC4 interface result is within 5 percent of the
  published value while producing exactly one conservative grazing warning.
- Interface limit: the complete-field result remains finite and continuous as
  the observation height approaches the interface from above.
- Parallel identity: the complete NEC4 result is exactly identical at one and
  automatic workers once the frequency/material tables are cached.
- Legacy regression: the full existing suite remains green.

The machine-readable evidence is
`tests/nec4-similarity/m10b-surface-reference.json`; the executable gate is
`tests/nec4-similarity/m10b_surface_reference.py`, and the `10-j1` manual
witness is `tests/nec4-similarity/m10b-surface-reference.nec`.
