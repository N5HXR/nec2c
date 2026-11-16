# NEC4-similarity project plan

Status: approved by the user on 2026-11-16; M0-M3 accepted and M4 in progress.

Date: 2026-11-16  
Branch: `codex/nec4-similarity`  
Starting commit: `9f7f108456a28b3f961403d9bad2769ad6c55fcf`

## Project objective

Implement an independently derived, numerically validated **documented NEC4 wire formulation** in this nec2c fork while preserving the complete classic NEC2 path. The project includes the documented NEC3/NEC4 homogeneous-half-space formulation for wires above ground, buried in ground, and crossing the planar air/ground interface.

The primary target is the 1992 public theory manual. Later public reports and the NEC5 validation manual are validation evidence and may justify independently derived improvements, but they do not silently redefine the NEC4 target.

The intended outcome is:

- materially better behavior for radius changes, nonsymmetric and unequal-radius junctions, electrically short segments, free ends, and source/load gaps;
- support for wires in either medium and wires penetrating a single planar interface;
- a reproducible original NEC4-style Sommerfeld/interpolation/asymptotic ground implementation;
- a conditional, separately named adaptive ground implementation inspired by the published NEC4.2 GN3 description if it can be specified and validated independently;
- unchanged NEC2 card meaning and result compatibility for existing models; and
- explicit model diagnostics where the thin-wire or ground approximation is outside its validated envelope.

This is not a NEC5 project. Mixed-potential wire bases, general open PEC surfaces, NEC5 wire/surface junctions, and NEC5-equivalent behavior are deferred.

## Governing decisions

### Compatibility and selection

The classic implementation remains a first-class engine, not a deleted reference. The provisional command-line contract is:

```text
--wire-engine auto|nec2|nec4
--ground-engine auto|nec2|nec4|adaptive
--model-check warn|strict|off
```

Exact names are frozen in M0 after user review.

- During development, explicit `nec4` selection is required.
- After final acceptance, `auto` selects the validated NEC4 engine for supported models and falls back to NEC2 only for a documented unsupported combination. Every fallback is reported in human and machine output.
- `--wire-engine nec2 --ground-engine nec2` retains the original implementation and is the opt-in route to old behavior.
- Existing NEC cards keep their meanings. In particular, existing `GN` decks are not reinterpreted merely because a new backend exists.
- A possible `GN 3` card is added only if its public syntax and semantics are verified. The command-line backend remains available regardless.
- Machine results record requested and resolved wire and ground engines, fallback reason, validity diagnostics, and ground-table provenance without removing or renaming existing fields.

### Numerical identity and accuracy

- Classic NEC2 tests and normalized outputs remain unchanged under explicit classic selection.
- No performance optimization may relax a numerical tolerance.
- Serial and parallel assembly of the same new operator must use the same per-entry arithmetic and produce bit-identical matrices where practical. Otherwise the acceptance metric is a scaled matrix difference of `1e-13` or less with identical reported results to the frozen output tolerances.
- Production double-precision formulas are validated against independent high-precision or direct-quadrature oracles. The oracle code must not call the production implementation.
- Published NEC4/5 plots are tolerance-band witnesses, not fabricated exact values. Each reconstructed deck records geometry, interpretation, digitization uncertainty, and source citation.

### Scope boundary

Included:

- free-space wires;
- PEC ground and the existing reflection-coefficient ground path;
- original documented NEC4 Sommerfeld/interpolation/asymptotic ground behavior;
- homogeneous lossy half-space with wire sources and observers above or below it;
- wires crossing that single planar interface;
- far and near fields, patterns, power accounting, loads, networks, and compatible classic patches insofar as the documented NEC4 wire formulation supports them;
- performance parity with the fork's accepted LAPACK, frequency, fill, field, and real-ground parallel work.

Excluded unless a later plan revision is approved:

- multilayer ground;
- terrain or nonplanar interfaces;
- NEC5 mixed-potential formulation;
- a new surface-EFIE solver or a claim of accurate wire/surface junction impedance;
- NEC4 insulated-wire equivalence; the existing `LD 7` coated-wire model remains independently named and must be separately requalified before combination;
- proprietary NEC4/5 source, binaries, internal tables, or licensed documentation;
- exact NEC4.2 GN3 compatibility when the public material does not specify an implementation detail.

## Provenance and clean implementation rules

1. Production equations are traceable to public documents by section/equation reference.
2. No licensed NEC4 or NEC5 binary or source is obtained or executed.
3. Publicly posted NEC4/5 numerical reports may be used as black-box witnesses when their provenance and complete input are recorded.
4. Exact, extended, and reduced thin-wire kernels receive distinct names. The production NEC4 mode uses the documented reduced kernel; exact or extended kernels may be test oracles and future engines.
5. GN3-inspired work is implemented from public mathematics and independently selected error targets. Missing NEC4.2 thresholds or tables are not guessed and labeled proprietary-equivalent.

## Acceptance protocol

Every milestone ends at a hard user gate. The next milestone does not begin until the user returns `PASS`, `FAIL`, or `PASS WITH CHANGES`.

Each acceptance bundle contains:

1. a design and implementation summary;
2. changed files and mathematical provenance;
3. automated regression, convergence, sanitizer, and applicable race-test results;
4. numerical errors against independent oracles and the previous accepted milestone;
5. phase timing and peak memory for affected workloads;
6. known limitations, warnings, and unresolved discrepancies;
7. retained JSON/CSV data and figures; and
8. exact commands and models for manual user validation.

User acceptance always includes at least one executable model. Plot-heavy milestones provide a compact summary plus machine-readable files so voluminous output need not be pasted into chat.

## Measurement contract

M1 freezes exact tolerances after the independent oracles and representative conditioning are measured. The initial targets are:

| Quantity | Initial target |
| --- | --- |
| Element or local-junction formula vs high precision | `<= 1e-11` relative away from true zeros; scaled absolute metric near zeros |
| Direct Sommerfeld integral internal convergence | `<= 1e-9` scaled relative |
| Linear-system scaled residual | `<= 1e-11` for well-conditioned acceptance cases |
| Reciprocity defect | `<= 1e-11` for reciprocal, well-conditioned cases |
| Power-balance defect | `<= 1e-6` ordinary cases; separately justified tolerance for extreme loss/conditioning |
| Serial/parallel result agreement | exact where arithmetic order is preserved; otherwise `<= 1e-12` scaled |
| Production ground field vs direct oracle | original NEC4 engine: characterized envelope; adaptive target `<= 0.5%` E and `<= 1%` H outside field nulls |
| Region-switch continuity | adaptive target `<= 0.5%` scaled field discontinuity |
| Ordinary-model performance | no unexplained regression above 10% versus the same operator's accepted baseline |

Convergence is judged against a stable limit or independent solution, not by agreement between two adjacent coarse meshes alone. Every frequency-domain comparison covers complex impedance/current/field, not only magnitudes.

## Milestones

### M0 - contract, branch, and immutable baseline

Status: accepted by the user on 2026-11-16.

Goal: freeze scope, interfaces, provenance, source state, and the legacy reference before numerical work.

Work:

- approve this plan and unresolved design choices;
- create `codex/nec4-similarity` from the approved main-branch commit;
- record compiler, dependencies, CPU, linked libraries, and default runtime policy;
- build and hash an archived classic executable and normalized output corpus;
- record the branch base and preserve explicit classic selection in every later milestone;
- freeze representative free-space, PEC, real-ground, network, patch, coated-wire, dielectric, VSWR-sweep, and radiation-grid regressions;
- inventory public reference models from the theory, GN3, accuracy, and NEC5 validation reports.

Acceptance metrics:

- clean build and complete existing suite pass;
- archived classic executable reproduces the frozen corpus;
- no production numerical source change;
- user runs one ordinary wire model and one current real-ground model with the archived executable.

### M1 - independent numerical oracles and operator boundary

Status: accepted by the user on 2026-11-16.

Goal: establish trusted mathematical references and a clean engine interface before replacing any physics.

Work:

- add independent high-precision/reference evaluators for basis transformations, junction potentials, reduced/exact/extended wire kernels, cap fields, and direct Sommerfeld integrals;
- map all production consumers of basis coefficients, fields, sources, loads, networks, currents, power, and output;
- introduce a wire-operator interface while routing all models through the unchanged NEC2 implementation;
- define dimensionless error metrics and freeze the numerical acceptance grids;
- reconstruct the strongest sufficiently specified published NEC4/5 witnesses.

Acceptance metrics:

- reference evaluators demonstrate stable convergence at least two decimal orders tighter than production targets;
- classic output remains frozen;
- oracle/reference code is structurally independent of production kernels;
- user reviews kernel and Sommerfeld error maps and runs the oracle self-test.

Possible administration request: high-precision C oracle work may require Ubuntu packages `libmpfr-dev` and `libgmp-dev`. Development will stop and ask before relying on them; a Python `mpmath` oracle is an acceptable slower alternative.

### M2 - cancellation-resistant basis and current representation

Status: accepted by the user on 2026-11-16.

Goal: implement the documented NEC4 short-segment basis representation without yet changing junction or field physics.

Work:

- replace cancellation-prone coefficient algebra in the NEC4 engine with the documented constant/sine/`cos(ks)-1` representation;
- carry the representation consistently through current reconstruction, loads, networks, near/far fields, and reports;
- retain NEC2 storage and behavior in the classic engine;
- add low-frequency and electrically short-segment conditioning tests.

Acceptance metrics:

- coefficient/current reconstruction agrees with high precision to the M1 limits;
- ordinary NEC4-engine cases agree with the algebraically equivalent classic solution before later physics changes;
- low-frequency sweeps remain finite and smooth materially farther than the classic representation;
- user runs a frequency-scaled dipole/loop case and reviews current plots.

### M3 - per-end charge factors and local junction solve

Status: accepted by the user on 2026-11-16.

Goal: implement the documented NEC4 junction treatment.

Work:

- build stable junction topology with distinct factors at each segment end;
- implement the Appendix A scalar-potential matrix and documented local solve;
- cover equal/unequal radii, T/Y/star junctions, radius steps, symmetry, and interface-independent junctions;
- publish junction conditioning and conservation diagnostics.

Acceptance metrics:

- equal-radius two-member null case reduces appropriately;
- scalar-potential continuity, current continuity, reciprocity, and local residual meet M1 tolerances;
- stepped-radius models converge toward the independent/reference limit rather than the known NEC2 limit;
- no regression in the circular-mesh center-junction crash fixture;
- user runs stepped dipole and multiwire-junction models.

### M4 - NEC4 reduced thin-wire kernel and precision branches

Goal: implement the documented surface-current/axis-match wire interaction and its cancellation-resistant numerical evaluation.

Work:

- implement the NEC4 reduced thin-wire operator behind the engine boundary;
- implement and validate Appendix B near-axis, short-segment, far-distance, and related limiting branches;
- derive double-precision switch thresholds from error maps rather than copying single-precision constants blindly;
- retain exact and extended kernels as reference-only implementations;
- emit validity diagnostics for segment-length/radius, close parallel wires, and overlapping bend geometry.

Acceptance metrics:

- interaction grid meets M1 oracle tolerance with no discontinuity spike at a branch switch;
- half-wave thick dipole, square loop, close two-wire line, and radius-step convergence studies reproduce the documented qualitative ordering;
- serial/parallel fill agreement meets contract;
- user runs a thick dipole refinement and reviews accuracy/validity plots.

### M5 - end caps, sources, loads, and complete free-space engine

Goal: close the NEC4 wire surface consistently at free ends and discontinuities and qualify the complete free-space engine.

Work:

- implement flat end-cap charge and field terms at free ends;
- implement corresponding source-gap and load-gap corrections;
- audit voltage sources, current-source behavior in scope, distributed/applied-field sources, loads, networks, and maximum-coupling port handling;
- diagnose source segments that are strongly mismatched to adjacent lengths;
- complete free-space power and field validation.

Acceptance metrics:

- cap and source fields meet independent oracle limits;
- source-length sweeps reproduce the published error trend and identify a validated range;
- resonance, impedance, currents, patterns, reciprocity, and power pass the frozen free-space corpus;
- ordinary thin uniform dipoles do not materially regress;
- user runs source-refinement, free-end, and loaded/network examples.

M5 is the first candidate for a useful **free-space NEC4-similarity preview**, but it is not the final project release.

### M6 - existing above-ground modes and combined-feature audit

Goal: integrate the accepted NEC4 wire engine with PEC ground, reflection-coefficient ground, and the fork's existing NEC2 Sommerfeld path before adding buried wires.

Work:

- qualify image/reflected interactions under the new basis and cap conventions;
- validate near/far fields, radiation patterns, power, and real-ground parallelism;
- test patches, networks, coated wires, and dielectric coupling and explicitly disable combinations not yet requalified;
- preserve current performance concurrency and JSON-safe worker-output behavior.

Acceptance metrics:

- above-ground canonical and Antenna Turtle models meet frozen numerical tolerances;
- real-ground serial/parallel results agree under the contract;
- no unexplained loss of accepted real-ground speedup;
- unsupported combinations fail or warn explicitly rather than mixing operators silently;
- user runs a VSWR sweep and radiation grid for an above-ground production model.

### M7 - homogeneous lower medium and interface topology

Goal: introduce the geometry and material foundations for buried and penetrating wires without yet relying on accelerated ground tables.

Work:

- allow wires wholly below the planar interface and wires split exactly at crossings;
- implement medium-dependent complex wavenumber, intrinsic impedance, loss, basis parameters, and output metadata;
- implement the documented interface charge/junction condition for crossing wires;
- define geometry tolerances for points on the interface and deterministic splitting requirements;
- start with homogeneous-medium and direct-field tests.

Acceptance metrics:

- vacuum-equals-vacuum and identical-media null limits reproduce the free-space engine;
- homogeneous lossless/lossy medium results satisfy wavelength and impedance scaling checks;
- crossing-interface charge/current conditions meet local residual and reciprocity targets;
- geometry at, above, and below the interface is classified deterministically;
- user runs a buried dipole and a ground-stake crossing model through geometry/preflight and direct reference mode.

### M8 - direct Sommerfeld field engine for all medium pairings

Goal: create a slow, accurate production-independent reference for source and observer above/above, below/below, and opposite sides of the interface.

Work:

- implement contour integration, singularity subtraction, branch/pole handling, reciprocity transforms, and boundary-condition transforms from the public theory;
- evaluate E directly and add a direct H reference where public derivation supports it;
- validate near-coincident, interface, high-conductivity, low-loss, grazing, small-rho, and cancellation regions;
- retain convergence and error estimates for every sampled reference value.

Acceptance metrics:

- direct integrals converge to the measurement-contract target under independent contour/quadrature refinements;
- Maxwell/interface boundary conditions and reciprocity hold across the qualification grid;
- agreement with analytic quasistatic, PEC, far-ray, and identical-media limits;
- user reviews field maps and runs a selected direct-mode point/short sweep.

The direct engine may be too slow for general matrix fill. It is retained as a correctness oracle and optional diagnostic path.

### M9 - documented NEC4 ground tables and asymptotics

Goal: implement the original documented NEC4 production ground evaluator for above, buried, and crossing interactions.

Work:

- implement two- and three-parameter table regions, singular-term extraction, transformations, linear/model-based interpolation, and documented table extents;
- implement the documented first/second-order and surface-wave asymptotic paths and switch regions;
- add exact-key table caching with material, frequency, extent, precision, and algorithm version in the key;
- parallelize table generation only after serial correctness is accepted.

Acceptance metrics:

- characterize error against M8 across the complete qualification grid, including field-null-aware metrics;
- no unexplained discontinuities at table/subregion/asymptotic switches;
- reproduce documented ground-report trends and sufficiently specified values within their uncertainty;
- cached and freshly generated tables give identical results and corrupt/incompatible caches are rejected;
- user runs the 5 m seawater whip/ground-stake witness and one near-interface field sweep.

### M10 - end-to-end buried and interface-crossing antennas

Goal: qualify the complete documented NEC4-similarity feature on realistic antennas.

Cases include:

- buried horizontal and vertical dipoles;
- a monopole with ground stake;
- a wire crossing the interface normally and at an angle;
- radial ground screens above and at small burial depth;
- a long wire supporting a surface-attached mode;
- near-ground horizontal wires;
- loss and power deposition in the lower medium;
- near E/H fields, far fields where physically defined, VSWR sweeps, and radiation grids.

Acceptance metrics:

- convergence, reciprocity, residual, passivity, and power meet frozen per-case limits;
- direct M8 and accelerated M9 results agree within the characterized envelope;
- published NEC4/5 witnesses and analytical limits are met within declared uncertainty;
- 1/2/4/8/16/32-worker tests retain numerical identity and useful scaling on large cases;
- user runs at least one buried/crossing model and one above-ground regression and reviews plots.

### M11 - conditional adaptive ground update

Goal: determine whether the public GN3 description and public field mathematics support a materially more accurate production evaluator without claiming exact NEC4.2 GN3 reproduction.

Proposed work:

- adaptive cubic table construction driven by explicit error estimates;
- separate branch-cut representations in interference regions;
- automatic choice between total-field and separated representations;
- direct H evaluation rather than numerical curl of interpolated E;
- improved large-distance treatment at grazing incidence and small-rho/large-z cases;
- persistent, versioned table cache.

Entry gate:

- M8-M10 must be accepted;
- every production formula must have public provenance or an independently derived reference;
- the user approves the final name and accuracy contract after reviewing an M11 research note.

Acceptance metrics:

- meet the adaptive E/H targets in the measurement contract over the frozen difficult-domain grid;
- eliminate the documented GN2 H-field glitches and preserve two-ray interference beyond the old asymptotic switch;
- outperform M9 accuracy materially for seawater/lower-HF and near-interface H cases;
- no material ordinary-case performance regression when cached; uncached table cost is reported separately;
- user runs seawater whip and short-dipole/interface-interference cases.

If these conditions cannot be met, M11 produces a research report and is closed as deferred. The accepted M10 NEC4 ground engine remains a valid project outcome.

### M12 - release candidate, default policy, and merge

Goal: produce the reviewable final feature without obscuring classic behavior or unresolved limits.

Work:

- run the entire legacy and new suite under sanitizers and applicable race tools;
- run end-to-end performance and memory characterization;
- freeze automatic engine-selection and fallback rules;
- publish validity guidance, examples, diagnostics, source provenance, and machine-result schema additions;
- prepare a compact manual acceptance bundle;
- only after user PASS, commit, push, and merge through the agreed workflow.

Acceptance metrics:

- explicit NEC2 mode reproduces the M0 corpus;
- validated NEC4 mode becomes `auto` for supported cases only if the user approves;
- every unsupported or out-of-envelope case has an actionable diagnostic;
- no unreviewed dependency or system-administration requirement;
- user performs the final free-space, above-ground, buried/crossing, VSWR, radiation-grid, and performance smoke tests.

## Planned manual acceptance portfolio

The final portfolio will be assembled incrementally and retained as ordinary NEC decks plus scripts:

- uniform and stepped-radius dipoles;
- tapered, T, Y, star, loop, and close-parallel-wire geometries;
- source-segment and load-gap sweeps;
- PEC and real-ground monopoles/dipoles;
- near-interface horizontal wire and long surface-wave wire;
- buried dipoles and ground stakes;
- normal and angled interface crossings;
- radial screens;
- VSWR sweeps and radiation grids from Antenna Turtle-exported models;
- deliberately invalid thin-wire, close-wire, and wire/surface-junction cases to verify diagnostics.

## Risks and stop conditions

| Risk | Response |
| --- | --- |
| Theory scan ambiguity or transcription error | Cross-check equations against dimensional limits, high precision, cited public papers, and independent derivation before production use. |
| No licensed NEC4 oracle | Use equation-level oracles, analytic limits, NEC5's independent published comparisons, openEMS where appropriate, and community reports with full decks. Claims remain capability-specific. |
| GN3 underspecification | Implement original documented NEC4 ground first. Treat adaptive ground as independently specified and separately named. |
| Ground interpolation error near cancellation/nulls | Use complex scaled absolute metrics and direct-integral witnesses rather than naive relative error. |
| Severe runtime of direct integrals | Retain direct mode as oracle; use cached accelerated tables for production. |
| Existing extensions assume NEC2 basis details | Disable combinations until explicitly requalified; never silently combine incompatible operators. |
| Scope growth into NEC5 or general surfaces | Require a new approved plan. |

## Approved design decisions

1. **Default policy:** after M12 acceptance, `auto` selects NEC4 for supported cases, while explicit `nec2` preserves old behavior. Until M12, new behavior is opt-in.
2. **GN3 naming:** call the independently derived later engine `adaptive`, not `gn3`, unless subsequent research establishes a complete public specification.
3. **Validity handling:** warn and continue by default, provide `--model-check strict` to turn validated-envelope warnings into errors, and `off` only to silence warnings. Numerically undefined geometry remains an error in every mode.
4. **Interface scope:** use one planar homogeneous half-space, matching documented NEC4. Multilayer soil and terrain remain future projects.
5. **Release granularity:** allow an accepted M5 free-space preview and M10 documented-ground preview on the feature branch, but merge only once M12 is accepted.

## Project completion rule

The project is complete when M12 receives an explicit user PASS and the accepted branch is merged. M11 may finish as either an accepted adaptive engine or an accepted deferral report; inability to claim exact GN3 equivalence does not invalidate the documented NEC4 result.
