# NEC2C-N5HXR integration enhancement specification

Status: **proposal for parallel solver development**

Companion client plan:
[NEC2C-DIELECTRIC-MODELING-PLAN.md](NEC2C-DIELECTRIC-MODELING-PLAN.md)

Target baseline: `nec2c-n5hxr 1.3.1`, dielectric feature revision 1.

## 1. Purpose

Define small, backwards-compatible solver enhancements that make
NEC2C-N5HXR safe and maintainable as a subprocess backend for Antenna Turtle
and other modeling tools. The numerical formulation and current human-readable
report remain authoritative for people. This specification adds a stable
machine contract, explicit capability discovery, and preflight diagnostics so
clients do not have to reverse-engineer report layout or duplicate solver
rules.

No item in this document changes the result of a traditional deck without a
`DF` card. New command-line options are opt-in. Existing `.nec`, `.diel`, and
human `.out` files remain accepted and produced as they are today.

## 2. Division of responsibility

NEC2C-N5HXR owns:

- dielectric resource grammar and electromagnetic interpretation;
- box meshing, active-cell allocation, and clearance rules;
- ground and conductor compatibility;
- solver selection, convergence, numerical residuals, and refusal;
- coupled power accounting;
- complex field/current results and their units; and
- stable machine-readable diagnostics for those facts.

The client owns:

- its modeling vocabulary and turtle transforms;
- generation of valid global axis-aligned boxes and material records;
- scenario selection and batch orchestration;
- conversion of raw complex fields to plots and comparisons;
- application-specific visualization and reports; and
- preservation of solver provenance alongside results.

The client may catch obvious errors early, but the solver remains the final
authority. A solver error must identify the offending resource object where
possible.

## 3. Enhancement priorities

### Priority 0: document and freeze existing producer contracts

Before adding interfaces, retain regression fixtures for:

- deck-relative `DF` path resolution;
- resource-format versions 1, 2, and 3;
- canonical material/region/box ordering;
- exact-touch versus interior-overlap behavior;
- mesh defaults and subdivision summaries;
- conductor, patch, and ground clearance diagnostics; and
- the current `nec2c-n5hxr 1.3.1` version string.

These already appear to work and should not be reimplemented in the client.

### Priority 1: machine-readable capability discovery

Add:

```text
nec2c --capabilities-json
```

It writes one JSON object to stdout and exits without reading a deck. The
object is intended for feature negotiation, not marketing metadata:

```json
{
  "program": "nec2c-n5hxr",
  "program_version": "1.3.1",
  "schema": "nec2c-capabilities",
  "schema_version": 1,
  "features": {
    "dielectric_feature_revision": 1,
    "dielectric_resource_versions": [1, 2, 3],
    "result_json_versions": [1],
    "validate_only": true,
    "solvers": ["dense", "gmres", "fft"],
    "grounds": ["free-space", "pec", "sommerfeld-norton"]
  },
  "limits": {
    "nec_card_characters": 132,
    "default_max_cells": 1000,
    "fft_accounted_storage_bytes": 536870912
  }
}
```

Rules:

- JSON goes to stdout; warnings/errors go to stderr.
- Success exits zero. Failure emits no partial JSON and exits nonzero.
- Unknown fields are permitted. Existing field meanings do not change within
  one schema version.
- Arrays communicate availability, not that every combination is compatible.
- `-v` remains unchanged for people and legacy probes.

### Priority 1: versioned machine result stream

Add:

```text
nec2c --result-json result.json -i model.nec -o model.out
```

Decision: schema version 1 is file-only. Stdout result streaming is deferred.
The file form cannot mix with process diagnostics, supports atomic completion,
and survives retained-artifact workflows. The human `.out` remains generated.

The top-level format is:

```json
{
  "schema": "nec2c-results",
  "schema_version": 1,
  "producer": { ... },
  "input": { ... },
  "jobs": [ ... ],
  "diagnostics": [ ... ],
  "completion": { "status": "ok" }
}
```

The solver writes to a temporary sibling and renames the completed document so
a killed process cannot leave an apparently complete result. Once option
processing has accepted a writable result destination, handled validation and
solver failures should publish a complete failure document and exit nonzero.
This is best-effort, not an absolute guarantee: a signal, crash, out-of-memory
condition, forced termination, or result-filesystem failure may leave no
document. A client always treats process status as authoritative and accepts a
result file only when it exists, parses completely, and contains a terminal
`completion` object. `completion.status` is `ok`, `validation-error`,
`solver-error`, or `internal-error`.

All numbers are JSON numbers only when finite. A non-finite numerical result is
a solver failure, not a quoted `NaN`. Complex values use objects:

```json
{ "re": 1.25, "im": -0.5 }
```

Every physical array declares units once at its containing section. Angles are
degrees, frequency is MHz, coordinates are meters, impedances are ohms, E is
V/m, H is A/m, range-normalized far E is V, wire current is A, and volume
polarization-current density is A/m^2.

### 3.1 Legacy no-DF coverage

Decision: result schema version 1 works for legacy no-`DF` jobs as well as
active dielectric jobs. It does not need to encode every historical
human-report table immediately. Its required legacy subset is the data already
consumed by Antenna Turtle: execution/frequency identity, feed impedance, wire
and patch currents, power, `RP` results, and rectangular `NE`/`NH` fields. This
gives the client one result path and provides a no-`DF` regression oracle
without turning schema version 1 into a transcription of every NEC table.

The absence of a dielectric resource is represented by
`"dielectric": null`; it is not an error and does not select a new numerical
path. Outputs outside the required schema subset may remain human-output-only
in version 1, but the result document must list any requested output omitted
from machine data so omission cannot be mistaken for an empty physical result.

## 4. Result schema version 1

### 4.1 Producer and input provenance

```json
"producer": {
  "program": "nec2c-n5hxr",
  "program_version": "1.3.1",
  "dielectric_feature_revision": 1,
  "build_id": "optional source revision or release id"
},
"input": {
  "deck": "model.nec",
  "dielectric_resource": "model.diel",
  "dielectric_resource_version": 2,
  "deck_sha256": "optional",
  "dielectric_sha256": "optional"
}
```

Paths may be basename-only to avoid leaking host directory names. Hashes are
strongly recommended for reproducibility but may be omitted in the first
implementation.

### 4.2 Jobs and repeated requests

One NEC deck can contain multiple frequencies, `NX` jobs, and repeated output
requests that trigger separate solves. Decision: each job contains one flat,
ordered execution array. Every execution repeats its `frequency_mhz`; there is
no unique-frequency wrapper. This preserves actual execution order, permits
interleaved or repeated frequencies, and lets the writer stream completed
executions without revisiting an earlier object:

```json
"jobs": [{
  "job_index": 0,
  "ground": { "kind": "free-space" },
  "dielectric": { ... },
  "executions": [{
    "execution_index": 0,
    "frequency_mhz": 1296.0,
    "excitation": { "kind": "driven-voltage" },
    "trigger": {
      "kind": "RP",
      "request_index": 0,
      "parameters": { ... }
    },
    "solver": { ... },
    "feeds": [ ... ],
    "wire_currents": [ ... ],
    "volume_currents": { ... },
    "power": { ... },
    "far_fields": [ ... ]
  }]
}]
```

`execution_index` is zero-based, unique within the job, and monotonically
increases in actual execution order. `request_index` identifies the source
control-card request within the job. Trigger parameters reproduce the
interpreted angular or spatial request sufficiently to validate sample counts
and axes. Do not make a client infer whether repeated current/power blocks
represent a new frequency or a new field request.

Each execution declares excitation semantics. Initial `kind` values are
`driven-voltage` and `plane-wave`. A plane-wave excitation additionally records
the interpreted incident direction, polarization or axial ratio, amplitude,
and phasor convention. Each job declares the interpreted ground model and its
parameters. These fields are required by the Antenna Turtle report and prevent
normalization or hemisphere rules from being inferred from absent quantities.

### 4.3 Geometry/resource summary

Each job reports the interpreted dielectric model:

```json
"dielectric": {
  "materials": [{
    "id": 1,
    "law": "constant-sigma",
    "er": 3.0,
    "sigma_s_per_m": 0.01
  }],
  "regions": [{ "id": 1, "material_id": 1 }],
  "boxes": [{
    "box_id": 1,
    "region_id": 1,
    "bounds_m": { "min": [0.1, -0.1, 0.2], "max": [0.3, 0.1, 0.5] },
    "divisions": [4, 4, 6],
    "cell_size_m": [0.05, 0.05, 0.05],
    "first_cell_id": 1,
    "cell_count": 96
  }],
  "mesh": {
    "max_kh": 0.5,
    "max_aspect": 2.0,
    "max_cells": 1000,
    "active_cells": 96,
    "volume_unknowns": 288,
    "total_unknowns": 299,
    "projected_dense_matrix_bytes": 1430416
  }
}
```

`box_id` is a stable canonical ordinal even though the resource grammar keys
boxes by region. Cell IDs are stable within one interpreted resource and fixed
mesh sweep. Table-valued materials additionally report the bracketing samples
and interpolated properties used at each frequency.

### 4.4 Solver convergence

```json
"solver": {
  "requested": "dense",
  "used": "dense",
  "relative_tolerance": 1e-12,
  "maximum_iterations": 300,
  "iterations": 0,
  "reported_residual": 2.1e-14,
  "true_residual": 2.1e-14,
  "converged": true,
  "accounted_storage_bytes": 1430416,
  "elapsed_seconds": 0.028
}
```

For direct solves, iteration fields may be zero. `used` must never imply a
fallback that did not occur. If fallback is ever introduced, it requires an
explicit command-line policy and a diagnostic.

### 4.5 Feed/input parameters

```json
"feeds": [{
  "source_index": 0,
  "tag": 1,
  "segment": 6,
  "voltage_v": { "re": 1.0, "im": 0.0 },
  "current_a": { "re": 0.01, "im": -0.002 },
  "impedance_ohm": { "re": 96.1538, "im": 19.2308 },
  "accepted_power_w": 0.005
}]
```

Multiple sources remain distinct. If accepted power is a system total rather
than independently attributable, omit the per-feed value and report only the
power total below.

### 4.6 Wire and volume currents

Volume-current output must scale to at least the documented FFT-sized models.
Avoid repeating property names and units for every cell. Use a columnar row
encoding (or an equivalently compact versioned representation):

```json
"wire_currents": [{
  "tag": 1,
  "segment": 1,
  "center_m": [0.0, 0.0, 0.1],
  "length_m": 0.02,
  "current_a": { "re": 0.003, "im": -0.001 }
}],
"volume_currents": {
  "encoding": "rows-v1",
  "columns": [
    "cell_id", "box_id", "region_id", "material_id",
    "center_x_m", "center_y_m", "center_z_m",
    "size_x_m", "size_y_m", "size_z_m",
    "jx_re_a_per_m2", "jx_im_a_per_m2",
    "jy_re_a_per_m2", "jy_im_a_per_m2",
    "jz_re_a_per_m2", "jz_im_a_per_m2"
  ],
  "rows": [
    [1, 1, 1, 1, 0.125, -0.075, 0.225, 0.05, 0.05, 0.05,
     0.1, 0.02, 0.0, 0.01, -0.02, 0.03]
  ]
}
```

Returning box/region/material IDs with each cell avoids reconstructing
ownership from printed ordering. Cell bounds or size are necessary for
correct visualization and integration. The result writer should stream rows
to its temporary file rather than retaining a second complete copy in memory.
This is required for first-release polarization-current rendering in Antenna
Turtle and for later FFT-sized environmental models.

### 4.7 Coupled power

```json
"power": {
  "basis": "driven-accepted-power",
  "accepted_w": 0.005,
  "air_radiated_w": 0.0038,
  "lower_halfspace_flux_w": 0.0,
  "dielectric_absorbed_w": 0.0011,
  "structure_loss_w": 0.0,
  "network_loss_w": 0.0,
  "balance_defect_w": 0.0001,
  "balance_defect_relative": 0.02,
  "angular_quadrature": { "theta_count": 0, "phi_count": 0 }
}
```

Use `null` or omit a quantity that is not defined for the excitation/ground;
never reuse zero to mean unavailable. Keep lower-half-space flux separate from
ground absorption and dielectric absorption.

For a plane-wave execution, accepted power is undefined. Decision: set
`accepted_w` to `null` (or omit it if the schema marks it optional) and use
`"basis": "plane-wave-scattered-balance"`. Record the incident-wave
normalization in the execution's `excitation` object. Report only physically
defined scattered radiation, lower-half-space flux, dielectric absorption, and
the solver's corresponding balance/defect quantities. The balance excludes the
unperturbed incident wave's substrate absorption, matching the documented
revision-1 semantics. Do not invent an aperture area, accepted power, gain, or
efficiency for a plane wave.

### 4.8 Complex far field

```json
"far_fields": [{
  "request_index": 0,
  "kind": "driven-total",
  "normalization": "range-times-electric-field",
  "units": "V",
  "samples": [{
    "theta_deg": 0.0,
    "phi_deg": 0.0,
    "r_e_theta": { "re": 0.1, "im": 0.02 },
    "r_e_phi": { "re": 0.0, "im": 0.0 }
  }]
}]
```

The solver should return the physical complex field it computed. It may also
provide gain derived from accepted power as a convenience, but raw complex
components and accepted power remain mandatory. State the phasor sign
convention (`exp(+i omega t)`) in schema documentation. For plane waves,
identify scattered versus total fields explicitly.

### 4.9 Near fields and validity

```json
"near_fields": [{
  "request_index": 0,
  "field": "electric",
  "kind": "driven-total",
  "units": "V/m",
  "samples": [{
    "position_m": [0.0, 0.0, 1.0],
    "value": {
      "x": { "re": 0.1, "im": 0.0 },
      "y": { "re": 0.0, "im": 0.0 },
      "z": { "re": 0.2, "im": 0.1 }
    },
    "valid": true
  }]
}]
```

Decision: driven coupled `NE`/`NH` results use `kind: "driven-total"`, because
they contain the complete field radiated by the driven conductor/dielectric
system and there is no separate incident plane wave. Plane-wave `NE`/`NH`
results use `kind: "plane-wave-scattered"`, matching revision 1's exterior
scattered-field semantics. The same naming applies to corresponding far-field
sections where applicable.

Current revision-1 rules reject observation points too close to cells. If a
whole request is rejected, issue an error. If future output permits a mixed
grid, every sample needs an explicit validity value and optional diagnostic;
invalid samples must not masquerade as physical zeros.

## 5. Structured diagnostics

Every warning/error exposed through the JSON contract has:

```json
{
  "severity": "error",
  "code": "DIELECTRIC_CONDUCTOR_CLEARANCE",
  "message": "wire tag 1 segment 4 is too close to dielectric box 2",
  "context": {
    "job_index": 0,
    "wire_tag": 1,
    "wire_segment": 4,
    "box_id": 2,
    "actual_gap_m": 0.012,
    "required_gap_m": 0.025
  }
}
```

Codes are stable tokens. Messages may improve without a schema change.
Priority codes include:

- `DIELECTRIC_RESOURCE_PARSE`
- `DIELECTRIC_RESOURCE_VERSION`
- `DIELECTRIC_BOX_OVERLAP`
- `DIELECTRIC_CONDUCTOR_CLEARANCE`
- `DIELECTRIC_GROUND_CLEARANCE`
- `DIELECTRIC_OBSERVATION_CLEARANCE`
- `DIELECTRIC_MATERIAL_RANGE`
- `DIELECTRIC_CELL_LIMIT`
- `DIELECTRIC_SOLVER_INCOMPATIBLE`
- `DIELECTRIC_FFT_LATTICE`
- `SOLVER_NOT_CONVERGED`
- `POWER_BALANCE_EXCEEDED`

Diagnostics should include resource line/column and material/region/box IDs
where applicable. This is substantially more useful than teaching every client
to parse English errors.

## 6. Validation-only operation

The current `DF` plus `EN` no-excitation deck performs useful mesh/cost
preflight. Expose that intent directly:

```text
nec2c --validate-only --result-json validation.json -i model.nec -o model.out
```

Validation performs complete deck/resource parsing, material evaluation over
the requested sweep, meshing, cell budget checks, conductor/patch/ground
clearance checks, and requested-solver compatibility checks. It does not
factor a matrix or execute an excitation.

The result contains resource/mesh summaries and diagnostics. Exit zero means
the interpreted model is admissible for the requested solver, not that its
physics or discretization is accurate.

If implementing the command-line option would duplicate existing control flow,
an equivalent explicit `--preflight` mode that internally constructs the
no-excitation request is acceptable.

## 7. Optional enhancements after the minimum contract

### 7.1 Exact requested resource pitch

`MAX_KH` is an accuracy bound, not a pitch request. For comparisons and FFT
lattice construction, a producer may need predictable common cells. A future
resource version could add an explicit pitch/lattice record, for example:

```text
MESH PITCH .025 ORIGIN 0 0 0 MAX_CELLS 10000
```

This should be added only if the solver can define boundary divisibility,
touching-region consistency, and refusal semantics precisely. Antenna Turtle
must not emulate it by perturbing physical box dimensions.

### 7.2 Rotated cuboids

Do not add rotated-box syntax merely for client convenience unless the volume
kernel, quadrature, mesher, clearance tests, ground images, and FFT policy can
support it as a genuine physical primitive. The initial client can express
right-angle axis permutations without new solver work. Arbitrary rotation is a
later numerical feature, not a serialization tweak.

### 7.3 Stable source labels

The current resource grammar uses numeric IDs. A future optional quoted label
on material/region/box records would improve diagnostics, provided labels have
no electromagnetic meaning and numeric IDs remain authoritative. Until then,
clients keep their own ID-to-name manifest.

### 7.4 Direct gain fields

The solver may report accepted-power-normalized component and total gain in the
machine result, particularly for ground cases where normalization subtleties
belong in the solver. Raw complex range-normalized fields must remain present.
Any gain property declares its normalization and reference power.

## 8. Testing and compatibility

Add golden-schema tests for:

1. `--capabilities-json` validity and stable required fields.
2. A legacy no-DF deck with machine output, verifying unchanged human output
   and electrical results.
3. The canonical driven free-space dielectric example, including impedance,
   volume currents, accepted/radiated/absorbed power, residual, and complex
   `RP` fields.
4. PEC contact and positive-gap cases.
5. A GN 2 case with lower-half-space flux kept separate.
6. Version-3 interpolation and out-of-range refusal.
7. Repeated `RP`/`NE`/`NH` executions at one frequency with unambiguous
   execution indexes.
8. Every structured diagnostic code used by Antenna Turtle preflight.
9. Interrupted execution, proving incomplete JSON is not accepted as complete.
10. Locale independence and paths containing spaces.

Schema fixtures should be checked with a strict JSON parser and a small
published JSON Schema if maintaining one does not become the primary
implementation burden. Numerical fixtures compare values with documented
tolerances rather than textual floating-point formatting.

## 9. Parallel delivery plan

### N5-D1: discovery and output skeleton

- Implement capability JSON.
- Implement atomic result JSON with producer/input/job hierarchy.
- Emit the required legacy no-DF subset: execution identity, feeds, wire/patch
  currents, convergence, power, `RP`, and `NE`/`NH` data.
- Preserve human output and all existing tests.

### N5-D2: dielectric result coverage

- Emit canonical material/region/box/cell metadata.
- Emit volume currents, dielectric power, ground flux, and repeated-request
  execution identity.
- Add structured diagnostics for resource, mesh, clearance, and solver gates.

### N5-D3: validation mode

- Add explicit validation-only operation.
- Report subdivisions, cost, solver compatibility, and all available
  diagnostics without solving.
- Supply Antenna Turtle with retained success/failure fixtures.

### N5-D4: optional numerical/API improvements

- Evaluate explicit lattice/pitch controls.
- Evaluate direct gain output, stable labels, and mixed near-field validity.
- Consider rotated primitives only as a separately validated solver feature.

## 10. Handoff contract with Antenna Turtle

The companion Antenna Turtle plan adds no hidden mandatory result fields.
Its client requirements are explicitly covered here: execution and request
identity, frequency, excitation and ground semantics, feed impedance, accepted
power for driven models, raw complex far and near fields, wire/patch and volume
currents, cell geometry and ownership, dielectric absorption and separated
ground flux, solver convergence, structured diagnostics, and input/producer
provenance. Antenna Turtle retains its own numeric-ID-to-user-name manifest.
Accepted-power gain, realized gain, directivity, polarization, and phase views
are client derivations from those raw quantities; they need not be duplicated
by schema version 1.

Parallel work can begin with checked-in fixtures before either side is fully
implemented:

1. Freeze capability JSON example and result-schema version 1 examples.
2. Retain one successful canonical dielectric result and one structured
   clearance failure from NEC2C-N5HXR.
3. Antenna Turtle builds its decoder and adapters against those fixtures.
4. Antenna Turtle independently implements box/resource generation against the
   existing resource-format specification.
5. Integration swaps the fixture producer for the live binary and compares
   decoded objects.
6. Text parsing remains a temporary compatibility path, not the permanent
   interchange contract.

The only hard dependency for Antenna Turtle's first export milestone is the
existing `.diel` contract. Machine output blocks reliable result integration,
not geometry/API development. Conversely, NEC2C-N5HXR can implement the result
contract using canonical examples without depending on the turtle API.

## 11. Acceptance criteria

This integration enhancement is complete when:

1. A client can determine dielectric resource versions, result schema, and
   supported solvers without parsing the version banner.
2. A successful dielectric run exposes all physical data needed for impedance,
   absolute pattern/polarization, convergence, power, and volume-current views
   without parsing the human report.
3. Repeated executions and frequencies are unambiguous.
4. Validation failures identify stable object IDs and measured versus required
   constraints.
5. Incomplete/failed execution cannot be mistaken for a complete result.
6. Traditional decks and human output remain backwards compatible.
7. Canonical JSON fixtures are exercised in NEC2C-N5HXR and consumed by an
   Antenna Turtle integration test.
