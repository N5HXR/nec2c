# Dielectric volume extension: revision 1 complete

Date: 2026-10-21. This is feature revision 1 of the NEC2C dielectric extension;
the upstream program version remains 1.3. M1–M11 and M13 are accepted. M12
(conductor/dielectric contact and penetration) is intentionally deferred.

## Announcement

Revision 1 of NEC2C's dielectric-volume extension is complete. NEC wire and
conducting-patch currents now couple self-consistently to dielectric volume
polarization, enabling antenna/dielectric models with passive scalar loss,
frequency-dependent material tables and touching dielectric materials.

The release includes dense and regular-grid FFT-backed solvers, free-space and
PEC operation, dielectric contact with PEC ground, and GN 2 real-ground support
for wires and positively separated volumes. Existing decks without dielectric
resources retain the legacy NEC path.

Analytical references, 24 regression groups, native refinement studies and
independently refined OpenEMS comparisons support the accepted scope. A user
guide, runnable canonical examples, material-assumption reporting and a register
of deferred features accompany the implementation.

Conductor/dielectric contact, anisotropy and broader real-ground modes remain
future work. Example tree materials are explicit assumptions; this release
does not claim calibrated vegetation properties or full-size lossy real-ground
tree performance.

## Validation record

* The revision-1 solver binary recorded by the final M13 checkpoints has SHA256
  `be00d866ec457a3aa81ecb4f3af8345455d89854c1dc67c3cdca2fdd640e0baf`.
  A rebuild may have a different binary hash; checkpoint hashes identify the
  actual builds used for their measurements.
* M13-01/02/03/04 retain independent OpenEMS decay, material/XML, mesh/domain and
  component-wise impedance audits. The agreed component tolerance is
  max(2 ohm, 2% of the corresponding reference perturbation component).
* M13-04's largest component discrepancy is 0.829761 ohm. M13-03's largest is
  1.97370 ohm; the refined hardest case gives 1.69712 ohm. These are benchmark
  observations, not universal model error bounds.
* M13-05 retains five assumed-material tree scenarios, separate nominal mesh
  refinement, material-provenance reporting and a remote-ground limit check.
* The canonical example runner verifies completion, residuals and power checks
  for four small reference models. A separately recorded distribution check
  verifies the source archive can build, test and install from scratch.

This release note records feature completion and provides announcement text.
It does not itself publish a remote release or change the upstream version.

## Navigation

[User guide](DIELECTRIC-REVISION-1.md) ·
[Examples](../examples/dielectric-r1/README.md) ·
[Deferred work and simplifying decisions](DIELECTRIC-DEFERRED-WORK.md) ·
[Material uncertainty guidance](M13-05-TREE-MATERIAL-GUIDANCE.md) ·
[Design and acceptance history](DIELECTRIC-VOLUME-DESIGN.md)

## Release verification (2026-10-21)

`make distcheck` passed from a fresh out-of-tree source archive: 24 tests passed,
zero failed, and build/install/uninstall/cleanup checks completed. All four
canonical examples also passed with the archive's rebuilt binary and matched
the retained numerical reference observations. The tested archive's implementation,
test/tool scripts and NEC/material inputs matched the staged sources. Subsequent
changes were documentation clarifications and CSV whitespace configuration only.
Local documentation links and staged whitespace checks pass. Raw release logs
are at `/tmp/dielectric-r1-distcheck.log` and
`/tmp/dielectric-r1-archive-examples.log` on the validation VM.
