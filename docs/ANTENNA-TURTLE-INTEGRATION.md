# Antenna Turtle integration contract for NEC2C-N5HXR

Status: implemented and tested contract for `nec2c-results` schema version 1.

This document describes the interface emitted by this source tree. It is the
consumer contract for Antenna Turtle; it intentionally does not describe
proposed fields that are not emitted. The normative structural companion is
[`nec2c-results-v1.schema.json`](nec2c-results-v1.schema.json). Consumers must
validate required fields and types, ignore unknown fields, and reject an
unknown `schema_version`.

## 1. Process interface

Capability discovery:

```text
nec2c --capabilities-json
```

This writes one JSON object to stdout and exits without opening a deck. It
cannot be combined with deck, result, or validation options. The object has
`schema: "nec2c-capabilities"`, `schema_version: 1`, and advertises dielectric
resource versions, result versions, solver names, supported dielectric ground
families, and fixed limits. Unknown capability fields are forward-compatible.

Normal execution with machine results:

```text
nec2c --result-json result.json -i model.nec -o model.out
```

Validation without solving:

```text
nec2c --validate-only --result-json validation.json -i model.nec -o model.out
```

`--result-json` is file-only; `-` is not a valid destination. The traditional
human report is still written to `-o`. Neither the numerical path nor the
human report format changes when no new option is used. `-s`, `-r`, and `-n`
remain dielectric-solver controls and do not select a different classic NEC2C
solution path.

The producer writes a sibling temporary file and renames it only after writing
the terminal `completion` object. A signal abandons the temporary document.
Handled validation and solver failures publish a complete failure document
when the result filesystem remains usable. The client must require all three:

1. process exit status zero;
2. a completely parseable result file; and
3. `completion.status == "ok"`.

Any other combination is a failed run. Do not consume a stale result at the
destination from an earlier invocation; remove or uniquely name destinations
before starting a run.

## 2. Capability values

Version 1 currently reports:

- program `nec2c-n5hxr`, dielectric feature revision 1;
- dielectric resource versions 1, 2, and 3;
- result schema version 1;
- solvers `dense`, `gmres`, and `fft`;
- dielectric grounds `free-space`, `pec`, and `sommerfeld-norton`;
- `validate_only: true` and `legacy_result_output: true`;
- 132 NEC card characters, default 1000 cells, and a 512 MiB FFT accounted
  storage limit.

Availability arrays are not a cross-product guarantee. In particular,
Sommerfeld–Norton dielectric coupling is wire-only and does not support FFT;
PEC permits contact at `z = 0`; Sommerfeld–Norton requires strictly positive
dielectric and conductor clearance. Ground screens, GN 0 dielectric coupling,
and cliff-ground dielectric output are unsupported.

## 3. Document and job model

The top-level object contains exactly one ordered `jobs` array plus global
diagnostics and terminal completion:

```json
{
  "schema": "nec2c-results",
  "schema_version": 1,
  "producer": {
    "program": "nec2c-n5hxr",
    "program_version": "1.3.1",
    "dielectric_feature_revision": 1
  },
  "input": {"deck": "model.nec"},
  "jobs": [],
  "diagnostics": [],
  "completion": {"status": "ok"}
}
```

`input` is deck provenance only and uses a basename. Hashes and build IDs are
not emitted by this implementation. Dielectric resource name, format version,
materials, regions, and boxes are job-local because an `NX` deck can change
resources. A job without `DF` has `dielectric: null`.

Every job always contains:

- zero-based `job_index`;
- `dielectric` object or `null`;
- `mesh_states`, `ground_states`, `executions`, and `omitted_requests` arrays.

A validation job additionally has `validation: "dielectric-preflight"` or
`"classic-input-only"`. Validation executions are empty.

`NX` closes the current job and releases its states. Identifiers restart at
zero in the next job.

## 4. Identity and ordering

All indexes are zero-based except physical NEC segment, patch, box, and cell
identifiers, which are one-based as shown in their fields.

- `execution_index` increases in actual execution order within a job.
- `request_index` is one shared increasing sequence for result-producing
  `XQ`, `RP`, `NE`, and `NH` cards. A deferred solve caused by `EN` receives
  the next index. All frequencies and excitations produced by one card share
  its request index.
- Every realized dielectric mesh receives a new increasing `mesh_state_id`.
  A later `FR` invalidates the current mesh, and the next validation or solve
  realizes a new state even when subdivisions are identical.
- An execution with dielectric data references its exact mesh state. A classic
  execution uses `mesh_state_id: null`.
- Ground states are realized on execution. A later `GN` or `GD` forces a new
  increasing `ground_state_id` at the next execution even when its interpreted
  values match an earlier state.
- Dielectric cell identity is the pair `(mesh_state_id, cell_id)`.

Only states actually realized by validation or execution are serialized.

## 5. Execution contract

An execution always has identity, frequency, state references, excitation,
trigger, and solver fields. Successful executions also always contain all
result sections:

```text
material_evaluations  array
feeds                 array
wire_currents         array
patch_currents        rows-v1 object
volume_currents       rows-v1 object or null
power                 object or null
far_fields            array
near_fields           array
```

Empty means applicable but containing no values. `null` means not applicable
or unavailable. In particular, classic runs use `volume_currents: null`;
classic plane-wave runs use `power: null`; and an undefined scalar such as an
open-circuit feed impedance is `null`. Arrays are never replaced by `null`.

`status` is `ok` or `failed`. A failed execution retains identity, trigger,
excitation, state references, and any available solver metadata, but discards
partial physical arrays. Earlier successful executions remain present; later
unattempted executions are absent. The process exit and top-level completion
remain authoritative.

The excitation `kind` is `driven-voltage` or `plane-wave`. Every excitation
declares `phasor_convention: "exp(+i omega t)"`. Plane waves additionally
contain propagation direction, the complex electric field at the origin, EX
type, angles, and axial ratio.

Trigger parameters are interpreted values:

- `XQ`: mode;
- `RP`: mode, theta/phi counts, starts and steps in degrees, and range in m;
- `NE`/`NH`: coordinate mode, three counts, start vector, and step vector;
- deferred `EN`: an empty parameter object.

## 6. Solver metadata

Classic NEC2C executions report both `requested` and `used` as
`legacy-nec-direct`. Their tolerance, iteration count, residuals, and
convergence flag are `null`; execution `status: "ok"` records successful
completion.

Dielectric executions report `dense`, `gmres`, or `fft`. The direct dielectric
solver reports zero iterations and a computed true residual. Iterative solves
report the configured relative tolerance and limit, performed iterations,
reported/true relative residual, and convergence. `accounted_storage_bytes`
is the projected dense matrix size or the FFT solver's accounted storage.
`elapsed_seconds` is solver CPU time measured by this producer.

All JSON numeric tokens are finite. Detection of a non-finite required result
fails the execution and changes an otherwise successful completion to
`solver-error`; no `NaN` or infinity token is emitted.

## 7. Physical results and units

Complex values are objects with finite numeric `re` and `im` members.
Frequencies are MHz, Cartesian coordinates and lengths are m, patch area is
m², voltage is V, impedance is ohm, wire current is A, patch surface current
is A/m, volume polarization-current density is A/m², E is V/m, H is A/m, and
range-normalized far electric field `r_e_theta`/`r_e_phi` is V. Angles are
degrees.

Patch and volume currents use `encoding: "rows-v1"`. Consumers must map row
values through the accompanying ordered `columns` array, not hard-code an
unversioned positional layout. Patch rows include identity, center, area, both
tangent vectors, and Cartesian complex surface current. Volume rows include
cell/box/region/material identity, center, size, and Cartesian complex current
density.

`material_evaluations` is populated for tabulated materials on every
successful dielectric execution. It contains the zero-based lower and upper
sample indexes for that material, interpolation fraction, interpolated real
permittivity, loss kind/value, and complex relative permittivity. The complete
source table remains in the job dielectric object.

RP mode 0 is serialized. Its far-field block declares normalization and units.
Other classic RP modes remain in the human report and create a nonfatal
`omitted_requests` entry with code `UNSUPPORTED_RP_MODE`. They must not be
interpreted as an empty physical field.

Rectangular and spherical NE/NH requests serialize Cartesian observation
positions and complex Cartesian components. Each sample currently has
`valid: true`; a request failing clearance is rejected rather than emitting
mixed-validity samples.

Power `basis` is `driven-accepted-power` or
`plane-wave-scattered-balance`. Driven dielectric executions expose
accepted, air-radiated, lower-half-space flux when applicable, dielectric
absorption, and balance defect. Plane-wave accepted power is `null`. Classic
driven executions expose the existing NEC input/radiated/structure/network
budget; no new classic integration is introduced. Classic plane-wave power is
`null`.

When a computed relative balance defect exceeds `0.003`, the document includes
a nonfatal `POWER_BALANCE_EXCEEDED` warning. The run may still complete `ok`.

## 8. Validation and diagnostics

`--validate-only` fully parses the deck and dielectric resource, validates the
material sweep, realizes the mesh, checks cell limits, box overlap,
conductor/patch-to-box clearance, requested near-field clearance, planar
ground compatibility, and FFT lattice/storage compatibility. It does not
factor a matrix or execute an excitation. A no-DF deck is accepted as
`classic-input-only` and has no mesh states or executions.

Diagnostics have `severity`, stable `code`, mutable human `message`, and an
object `context`. Current dielectric diagnostics include job and resource line
when known. Clients must branch on `code`, never message text. Implemented
codes include dielectric parse/version, overlap, conductor, ground,
observation, material-range and cell-limit errors, solver incompatibility,
FFT lattice failure, nonconvergence, nonfinite result, and the power-balance
warning.

Completion status is one of `ok`, `validation-error`, `solver-error`, or
`internal-error`. A warning does not change completion status.

## 9. Required consumer behavior

An Antenna Turtle adapter should:

1. probe capabilities once per executable identity and require result version
   1 plus the needed dielectric resource version, solver, and ground family;
2. give every invocation a unique result path;
3. require zero process exit, complete JSON, matching schema/version, and
   completion `ok`;
4. preserve producer, deck provenance, job dielectric metadata, mesh state,
   ground state, solver metadata, and diagnostics with derived plots;
5. join cells only by `(mesh_state_id, cell_id)` and outputs by
   `(job_index, execution_index, request_index)`;
6. distinguish `null`, empty arrays, failed executions, omissions, and absent
   unattempted executions exactly as defined above; and
7. ignore unknown fields but never an unknown enum value or schema version in
   a field it consumes.

## 10. Verification performed in this tree

The automated `tests/result_json.py` contract test covers capability output,
strict JSON parsing, classic and dielectric success documents, classic
validate-only, classic direct-solver null semantics, dielectric currents and
RP samples, monotonic request/mesh/ground state IDs across repeated cards,
handled validation failure publication, and byte-equivalent classic human
reports with and without result observation (excluding the pre-existing run
time line). The full project test suite remains the regression authority for
classic NEC2C and dielectric numerical behavior.
