# M13 loss and real-ground handoff

Status: M13 COMPLETE (2026-10-21) for passive scalar materials and the agreed
GN 2 wires-first real-ground scope. M13-01 through M13-05 are complete. All 24
regression groups pass. Production-tree material assumptions and sampled
sensitivity reporting are documented in
[M13-05-TREE-MATERIAL-GUIDANCE.md](M13-05-TREE-MATERIAL-GUIDANCE.md).
GN 0/real-ground patches/contact expansion, anisotropy and magnetic media remain
outside the completed scope; M12 is separately deferred.

## Milestone boundary (2026-10-20)

The user selected a lossless proof of concept through M11. M10 covers positive-gap
PEC ground (`GN 1`); M11 covers dielectric contact with that perfect ground.
Real ground belongs to M13 alongside dielectric loss and advanced material
behavior. Its physical motivation includes substrate absorption, even though the
ground is a background model rather than a dielectric-resource material.

M12, conductor/dielectric contact, penetration and embedding, is deferred to a
later revision. It retains its identifier and proposed criteria for tracking.
It is not a prerequisite for M13. The current revision keeps supported positive
clearance between wire/patch conductors and dielectric volumes. M11's contact
against an ideal PEC ground plane is a separate boundary case. The planned path is M10, M11, then M13.
Extending PEC-ground contact to real ground, or adding conductor/dielectric
contact later, requires explicit scope and validation rather than automatic
inheritance from earlier milestones.

## Work already available

Real-ground work developed during M10 is retained as preliminary M13 work:

- a passive nonmagnetic half-space reflected E/H operator for points, finite
  cuboids and mixed point/cuboid interactions;
- positive-gap GN 2 wire/volume coupling with dense and GMRES solvers;
- consistent excitation, near/far fields and substrate-power accounting;
- an optional reflected-wire-basis cache bounded at 16 MiB;
- analytical Fresnel, quasistatic, PEC-image, reciprocity, matched-background,
  dipole-work/flux and Maxwell-curl checks;
- legacy-vacuum compatibility and frequency-sweep state checks; and
- a completed 64-cell, 300 MHz NEC/OpenEMS perturbation checkpoint.

The implementation and its limitations remain recorded in
[M10-GROUND-FORMULATION.md](M10-GROUND-FORMULATION.md), now explicitly a combined
PEC record and historical real-ground handoff. Retain those derivations and
results rather than duplicating or removing tested code.

The checkpoint in [m10-nec-real-300.json](../tests/m10-nec-real-300.json) differs
from the completed coarse OpenEMS dielectric impedance perturbation by
0.27573 ohm, below the agreed 2-ohm allowance at that one checkpoint. The power
defect is 0.0876%. This does not establish full M13 acceptance: the fine loaded
OpenEMS result, broader refinement and source/domain convergence remain open.
The existing M10-prefixed filenames are retained for provenance.

The real-ground subset currently excludes GN 0, real-ground patches,
contact/buried geometry, radial screens, cliff ground and ground FFT. GN 2
wires first was an implementation increment, not acceptance of every NEC ground
abstraction or contact configuration.

## Completed boundary and later work

M13-01 validates constant scalar loss/conductivity and absorbed power;
M13-02 validates frequency tables and extrapolation refusal; M13-03 validates
touching materials; M13-04 validates combined lossy volumes and GN 2 real ground;
M13-05 reports tree-material assumptions, sampled sensitivity and numerical
refinement separately, and records the final exit-criterion review.

The accepted real-ground subset is GN 2 with wires, strictly positive ground
and conductor clearances, and dense or assembled GMRES. GN 0, real-ground
patches, contact/burial, screens/cliff and ground FFT are separate future work.
Anisotropic or magnetic media also require a later scoped extension.
The historical M10 checkpoints below are provenance; current physical
acceptance is in M13-04. Full-size lossy or real-ground tree resource and
empirical vegetation-material claims are not inferred from these fixtures.

## M13-01 completed increment (2026-10-21)

Completed constant isotropic complex permittivity/conductivity in free space
and over PEC. Anisotropy is excluded from this revision's current scope;
wood grain/moisture directionality requires a later tensor-material extension.
Experiment identifiers use the actual milestone and step: M13-01, M13-02, etc.
Stage independent OpenEMS jobs alongside development; returned evidence gates
acceptance of the step rather than blocking implementation within it. See
[M13-01-SCALAR-LOSS.md](M13-01-SCALAR-LOSS.md) for material syntax, tests and the
self-contained first laptop package.

M13-01 is COMPLETE: scalar loss/conductivity, separate absorption accounting,
zero-loss/sweep limits, dense/GMRES/FFT checks, six native benchmark pairs and
six independently refined OpenEMS verification runs pass. All 20 regression
groups pass. Retained evidence is in
[m13-01-acceptance.json](../tests/m13-01-acceptance.json) and
[m13-01-openems-review.json](../tests/m13-01-openems-review.json).
Frequency tables and touching materials are completed by M13-02/03 below.
Remaining real-ground coverage is a later increment; M13 is not complete.

## M13-02 completed increment (2026-10-21)

Frequency-dependent isotropic material tables are now implemented with resource
version 3, linear ER/EI-or-conductivity interpolation, whole-sweep range checks,
and a conservative table-wide mesh bound. Native analytical/refinement checks and all 21 regression groups pass;
all 12 independent OpenEMS verification runs pass decay, material, mesh/domain
and component-wise impedance checks. Maximum component error is 1.33655 ohm.
Evidence is retained in `tests/m13-02-openems-review.json`. See
[M13-02-MATERIAL-TABLES.md](M13-02-MATERIAL-TABLES.md).

## M13-03 completed increment (2026-10-21)

Touching passive scalar materials, cell-local FFT permittivity, independent
interface/coated-sphere checks, all 22 regression groups and split-pillar
OpenEMS acceptance pass. All 11 native benchmark/refinement comparisons satisfy
the component-wise impedance gate. Largest component error is 1.97370 ohm;
the joint refined hardest case gives 1.69712 ohm. The 12 returned verification
runs pass explicit decay, two-material/interface XML and mesh/domain audits.
See [M13-03-TOUCHING-MATERIALS.md](M13-03-TOUCHING-MATERIALS.md) and
[m13-03-acceptance.json](../tests/m13-03-acceptance.json).
M13 remains in progress for real-ground integration and acceptance.

## M13-04 completed increment (2026-10-21)

Validate combined scalar dielectric and substrate loss for the retained GN 2
wire/volume implementation with strictly positive ground and conductor gaps.
All 23 regression groups, six native antenna pairs and all 12 independent OpenEMS verification runs pass.
Worst component discrepancy is 0.829761 ohm; maximum mesh/domain changes are
0.280831/0.029949 ohm. Native analytical, material/sweep, power, reciprocity and
refinement checks cover the combined-loss model. GN 0, real-ground patches,
contact/burial, screens/cliff and ground FFT remain separate work. See
[M13-04-REAL-GROUND-WIRES.md](M13-04-REAL-GROUND-WIRES.md).

## M13-05 completed increment (2026-10-21)

The production-tree guidance and material-report tool attach provenance,
geometry/ground assumptions, exclusions and material scenarios to results.
Five synthetic PEC-tree material cases and a separate nominal refinement pass;
localized-source remote-ground E/H limits pass; all 24 regression groups pass.
This closes M13 for the agreed scope without another OpenEMS run. Existing
M13-01/02/03/04 independent refinement evidence is retained. See
[M13-05-TREE-MATERIAL-GUIDANCE.md](M13-05-TREE-MATERIAL-GUIDANCE.md) and
[m13-05-acceptance.json](../tests/m13-05-acceptance.json).
