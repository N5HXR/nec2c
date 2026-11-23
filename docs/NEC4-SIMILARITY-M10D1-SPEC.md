# NEC4 Similarity M10D.1: Radial-Screen Discrepancy Diagnosis

## Goal

Locate the source of the Burke Figure 5 radial-screen RCE discrepancy before
changing production mathematics.  This is a diagnostic extension of M10D,
not an empirical curve-fitting exercise.

## Questions

1. Does the solved radial-screen current and feed impedance agree with the
   NECGS appendix evidence (Figures A1--A3)?
2. Is the RCE error introduced by the radial-screen current solution, the
   ground-wave propagation path, or the final field normalization?
3. Does the error begin abruptly at a ground-kernel routing boundary, or grow
   continuously with range, conductivity, burial depth, or radial density?

## Work plan

1. Digitize Burke Figure 5 and appendix Figures A1--A3, retaining plot-reading
   uncertainty and source coordinates.
2. Record the raw screen and reference-dipole `|Ez|^2 / P` values as well as
   their RCE ratio; do not diagnose from dB values alone.
3. Compare `Zscreen - (38.6 + j22.2) ohm` with Figures A1--A3.
4. Sweep observation range across the direct/asymptotic transition using one
   unchanged solved antenna, and compare against a reference dipole evaluated
   at exactly the same ranges.
5. Audit radial-current rotational symmetry, junction KCL, and the prescribed
   interface derivative relation from the computed solution.
6. Force the direct-contour and asymptotic ground paths independently at
   selected points, then compare their overlap before proposing a correction.
7. Repeat the decisive cases with tighter radial segmentation until the
   diagnostic quantity changes by no more than 0.05 dB (or document why that
   target is impractical).

## Decision gates

- **Current-solution defect:** appendix impedance or junction residuals fail.
- **Propagation defect:** impedance passes, but RCE changes with propagation
  route or direct/asymptotic implementations disagree in their overlap.
- **Normalization defect:** raw screen and reference fields agree separately,
  but their reported RCE does not.
- **Source ambiguity:** results stay within digitization or undocumented-model
  uncertainty; no correction is made without a second independent oracle.

## Acceptance evidence

- A reproducible JSON diagnostic bundle with incremental checkpoints.
- Digitized published curves with uncertainty estimates and provenance.
- A written root-cause statement supported by at least two independent
  discriminators.
- A proposed correction and regression tests, reviewed before implementation.

## Diagnostic outcome

The Figure 5 discrepancy was localized before a forced-route range sweep was
needed.  The screen numerator used the NEC4 ground-field path while the
denominator was copied from the report's NECGS absolute-field table.  The
result therefore exposed a common propagation-path scale difference as though
it were radial-screen error.  Evaluating the report's reference dipole through
the same path changes the denominator by `-0.163`, `-0.635`, `-1.156`, and
`-5.596 dB` for the four grounds, respectively.  Those corrections account for
the observed offsets, including the full high-conductivity discrepancy.

A fresh high-conductivity screen solve independently reproduces the corrected
retained-data value within `0.00025 dB`.  The production field mathematics is
unchanged: the absolute grazing-field discrepancy remains warning-backed
evidence for later NEC4.2-style work.  The remaining direct/asymptotic overlap
experiments belong to that later absolute-field project, not to Figure 5 RCE
qualification.

The requested closure pass also completed the two remaining numerical checks.
Sensitive low-loss curves converge at `0.01 lambda` radial pitch, with a
further halving changing RCE by no more than `0.0183 dB`.  All nine digitized
Appendix A NECGS impedance anchors pass their scan-reading uncertainty bands.
Both the converged controls and the deliberately coarse failure witnesses are
retained as regression evidence.
