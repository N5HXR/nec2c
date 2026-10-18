# Broadband networks and Antenna Turtle integration

## Status

This is a forward-looking design for the N5HXR `nec2c` fork.  It is not an
implementation commitment and must not change the behavior of existing NEC-2
decks.  The immediate objective is to let one `nec2c` frequency sweep include
the measured, frequency-dependent behavior of a passive transformer, unun,
balun, choke, transmission-line fixture, or other reciprocal two-port.

The design also records a wider roadmap for making `nec2c` a stronger,
reproducible backend for the Antenna Turtle wire-antenna DSL.

## Motivation

Classic NEC-2 `NT` cards describe a non-radiating reciprocal two-port through
three complex short-circuit admittance terms: `Y11`, `Y12 = Y21`, and `Y22`.
They are excellent for a fixed-frequency network, but their six values are
constant throughout an `FR` sweep.  This is insufficient for a practical RF
transformer: magnetizing inductance, leakage inductance, winding resistance,
core loss, and parasitic capacitance all vary materially across an HF or VHF
band.

At present Antenna Turtle can run NEC2C sweeps and report feedpoint impedance,
VSWR, radiation patterns, and gain, but its NEC exporter emits wires, ground,
and voltage sources only.  It deliberately treats simultaneous feed sources
as active impedances rather than pretending they are general S-parameters.
The backend should gain a general two-port facility first; transformer-oriented
DSL forms can then be thin, convenient wrappers around it.

## Goals

1. Preserve the numerical behavior and deck compatibility of existing `NT`
   and `TL` cards.
2. Evaluate a network matrix at every requested frequency within one NEC2C
   process and one `FR` sweep.
3. Support practical device data through Touchstone 1.0 two-port files and a
   self-contained tabular form.
4. Return enough provenance to reproduce results: network identity, source
   data hash, interpolation/extrapolation policy, evaluated Y matrix, and
   warnings at every frequency.
5. Provide a clean Antenna Turtle mapping for ideal, equivalent-circuit, and
   measured networks.
6. Build an accuracy workflow around segmentation, kernel selection, and
   convergence rather than implying that a richer network model alone makes a
   simulation accurate.

## Non-goals

- Do not copy, translate, or reverse engineer licensed NEC-3/NEC-4 code.
- Do not claim that a two-port represents common-mode suppression or all
  physical effects of a balun.
- Do not silently extrapolate measured network data.
- Do not change classic NEC-2 card meanings or existing report text unless an
  explicit extension was requested.
- Do not model a transformer as radiating wire/toroid geometry in this phase.

## Existing solver insertion point

The current source already has the correct sweep lifecycle.

1. `main.c` reads `FR`, `NT`, and `TL` cards into `netcx` arrays.
2. For every `FR` point it computes `save.fmhz` and `data.wlam`, fills and
   factors the electromagnetic interaction matrix, and resets
   `netcx.ntsol = 0`.
3. `netwk()` in `network.c` forms and factors the small network equation
   matrix for that frequency.
4. Fixed `NT` values are copied into the small matrix.  `TL` values are
   already converted at run time using `data.wlam`.

Therefore this feature does **not** require a new outer sweep loop, a new
electromagnetic method-of-moments matrix, or one process per frequency.  It
requires a general network evaluator called while the per-frequency network
matrix is assembled.

## Internal architecture

Replace the implicit `ntyp` conditionals with an internal network descriptor.
The storage may remain array-based initially to minimize churn, but its API
should be explicit:

```c
typedef struct {
  complex double y11;
  complex double y12;
  complex double y22;
} nec2c_y_matrix;

bool network_evaluate(const network_descriptor *network,
                      double frequency_mhz,
                      nec2c_y_matrix *result,
                      network_diagnostic *diagnostic);
```

Network variants:

- `NETWORK_STATIC_Y`: legacy `NT`; returns the six stored scalar values.
- `NETWORK_TRANSMISSION_LINE`: legacy `TL`; preserves its present wavelength
  calculation and crossed-line behavior exactly.
- `NETWORK_TABLE_Y`: new embedded frequency table.
- `NETWORK_TOUCHSTONE`: new imported two-port data.

`netwk()` must receive evaluated Y values only.  It should not know whether
they originated in a legacy card, line formula, table, or file.  This isolates
the new behavior and makes static regression testing straightforward.

Evaluate every network before both network reporting and assembly.  The report
must print the actual Y matrix used at the current frequency, never merely the
raw source-file value.

## Deck extension

Reserve new NEC2C-extension mnemonics rather than overloading `NT`.  The
precise fixed-column grammar should be settled before implementation, but the
semantics are mandatory:

```text
NB  <tag1> <segment1> <tag2> <segment2> <network-id> <options...>
NF  <network-id> <deck-relative-touchstone-path> <interpolation> <extrapolation>
```

`NB` binds a two-port connection to a network ID. `NF` supplies the broadband
provider.  A deck-relative path is resolved against the directory of the input
deck, not the current working directory.  This makes retained Antenna Turtle
artifact directories portable and avoids dependence on the process launch
location.

Alternative syntax is acceptable if it retains these properties:

- unambiguous two-port connection endpoints;
- a named provider usable by diagnostics;
- deck-relative resource resolution;
- an explicit interpolation and extrapolation policy; and
- no reinterpretation of a classic `NT` card.

For self-contained decks, define an embedded table provider.  Each sample
contains frequency plus all three complex Y terms.  The source parser may use
continuation records or a keyword/value grammar because seven real values do
not fit naturally in the six floating-point fields of an original NEC card.

```text
NY <id> <MHz> <ReY11> <ImY11> <ReY12> <ImY12> <ReY22> <ImY22>
```

The initial implementation must sort samples by frequency, reject duplicate
frequencies, and require at least two points for interpolation.

## Touchstone support

Start with a deliberately narrow, well-tested reader:

- two-port (`.s2p`) only;
- Touchstone 1.0 option lines;
- frequency units Hz/kHz/MHz/GHz;
- `S`, `Y`, and `Z` parameters;
- `RI`, `MA`, and `DB` data representations;
- real positive common reference impedance (`R 50` etc.).

Read the complete file during deck parsing.  Fail early with a source line
number for malformed records.  Store normalized complex samples in SI units,
with their original format and file hash retained only for diagnostics.

For measured S parameters with equal real port reference impedance `Z0`,
convert at an evaluated frequency using:

```text
Y = (I - S) (I + S)^-1 / Z0
```

For a diagonal unequal real reference impedance matrix `Z0`, use the
reference-normalized form:

```text
Y = Z0^(-1/2) (I - S) (I + S)^-1 Z0^(-1/2)
```

Do not form an explicit inverse in production code; solve the 2x2 complex
system and reject a poorly conditioned `I + S` matrix.  Z-parameter input is
converted by solving `Y * Z = I`.  Native Y input avoids conversion.

Touchstone 2.0, complex or unequal port-reference conventions, more than two
ports, and mixed-mode files are useful later work, not prerequisites.

## Interpolation and validation

Default behavior must favor trustworthy data over a deceptively smooth graph.

- Default interpolation: linear in the real and imaginary parts of each
  complex matrix entry.
- Optional later interpolation: shape-preserving cubic interpolation of real
  and imaginary components, never magnitude in dB plus wrapped phase.
- Default extrapolation: error.
- Optional explicit policies: `clamp` and `linear`, each producing a warning
  and being recorded per result point.
- Reject a requested frequency outside the measured range unless an explicit
  extrapolation policy permits it.
- Check finite values, 2x2 conversion conditioning, and reciprocity.

NEC-2 `NT` stores one off-diagonal value, so a broadband matrix must be
reciprocal within configured absolute and relative tolerances.  If
`Y12 != Y21`, default to an error.  An explicit `symmetrize` option may use
`(Y12 + Y21) / 2` only with a diagnostic warning and the original mismatch
recorded in metadata.

Passivity should be reported, not silently repaired.  For each evaluated
matrix, test that the Hermitian conductance matrix is positive semidefinite to
within numerical tolerance.  A warning is appropriate for noisy measurements;
an error should be available for automated optimization workflows.

## Transformer semantics

### Ideal transformation

For a turns ratio `n = Ns / Np`, an ideal transformer maps the secondary
impedance to the primary as:

```text
Zprimary = Zsecondary / n^2
```

Thus a nominal 49:1 EFHW unun has `n = 7`; a 2450-ohm antenna-side impedance
maps to 50 ohms at the radio side. A nominal 4:1 impedance balun has a 2:1
voltage ratio; 200 ohms maps to 50 ohms.

An ideal transformer does not alter antenna currents, fields, or radiation
efficiency.  It can be implemented as Antenna Turtle reference-plane
post-processing before any solver extension exists.  A perfectly ideal
transformer is singular as a simple finite Y-parameter network, however; do
not approximate it with arbitrary huge admittances.  Use a dedicated ideal
reference-plane transform or a finite equivalent circuit.

### Equivalent circuits

An RLC transformer builder should create a reciprocal frequency-dependent
two-port from:

- turns ratio;
- magnetizing inductance and core-loss resistance;
- primary and secondary leakage inductance;
- primary and secondary winding resistance;
- optional interwinding/shunt capacitance; and
- optional fixture or termination parasitics.

Derive the two-port through a small modified-nodal-analysis circuit solver,
then hand its Y matrix to the common broadband-network evaluator.  Keep this
circuit solver independent of the NEC electromagnetic matrix.

### Measured networks

The preferred practical transformer model is a calibrated Touchstone file.
It captures losses and parasitics without relying on a guessed equivalent
circuit.  Document the measurement reference planes, fixture de-embedding,
port order, common-mode termination, calibration method, and valid frequency
range alongside the `.s2p` file.

## Baluns, ununs, and common mode

A two-port accurately represents differential impedance transformation and
insertion loss. It cannot fully represent all physical behavior of a balun.

- A 49:1 EFHW unun may rely on the feedline or counterpoise as part of the RF
  return system; that conductor can radiate.
- A current balun has differential and common-mode behavior and is normally a
  three- or four-terminal problem.
- A transformer measurement made in a 50-ohm fixture does not automatically
  characterize its common-mode impedance in an installed antenna system.

Near-term treatment: model the antenna differentially with the broadband
two-port, then include the exterior coax/counterpoise as explicit wire
geometry when it is electrically important.  Model a common-mode choke as a
frequency-dependent series load on that explicit path.

Future treatment: add a general external multiport network/mixed-mode solver
that couples antenna port impedances to differential and common-mode circuit
ports.  Do not present a two-port transformer result as a prediction of
common-mode suppression.

## Antenna Turtle API sketch

The DSL should use a general network primitive rather than introducing
transformer-only machinery first:

```racket
(broadband-network
  name "efhw-49-1"
  from "radio"
  to "antenna"
  touchstone "measurements/efhw-49-1.s2p"
  interpolation 'linear
  extrapolation 'error)

(transformer
  name "loop-balun"
  primary "radio"
  secondary "loop"
  impedance-ratio 4
  model ideal)
```

`transformer` is syntactic/API sugar over either:

1. an ideal reference-plane transformation;
2. an equivalent-circuit provider; or
3. a `broadband-network` provider backed by Touchstone data.

For a measured sweep, Antenna Turtle copies the Touchstone sidecar into the
retained analysis artifact directory, emits the network extension with a
deck-relative resource path, launches **one** NEC2C run, and returns for every
frequency:

- antenna-side and radio-side impedance;
- VSWR, S11, and return loss at the explicitly named reference plane;
- evaluated Y matrix and input-data provenance;
- extrapolation, reciprocity, passivity, and conditioning diagnostics; and
- available/accepted/dissipated power where the model permits it.

The existing Antenna Turtle multi-feed distinction remains important: active
impedance under simultaneous drives is not a Z/S parameter matrix.  Full
multiport measurements require one excitation/termination configuration per
port and are a separate extension.

## Broader Antenna Turtle backend roadmap

### Phase A: confidence in the existing NEC-2 solution

1. Add a segmentation audit at every relevant frequency. Flag segment lengths
   above 0.05--0.1 wavelength, very short segments, poor segment alignment on
   close parallel wires, interior crossings, overlaps, and risky junctions.
2. Add kernel policy `thin`, `extended`, and `auto`. Keep `thin` as the
   compatibility default. `auto` recommends or emits `EK` when the local
   segment-length/radius ratio calls for it, while warning about short segments
   at bends.
3. Add opt-in mesh convergence runs with increasingly fine segmentation and
   report deltas in feed impedance, gain, F/B, and lobe direction.
4. Segment broad sweeps from their highest requested frequency, while
   preserving a unique center segment for center feeds.

### Phase B: circuit/network capabilities

1. Implement the broadband two-port feature in this document.
2. Expose `LD`, `NT`, and `TL` through typed Antenna Turtle forms: lumped
   R/L/C loads, distributed conductor loss, transmission lines, traps, and
   matching networks.
3. Add source/reference-plane and power-budget metadata.
4. Offer stable machine-readable output (for example JSON sidecar output)
   while retaining classic text output unchanged by default.

### Phase C: ground and advanced models

1. Keep NEC-2 for conventional free-space and above-ground wire models.
2. Continue routing buried wires, ground stakes, complex dielectric objects,
   and conductor/insulator solids to Antenna Turtle's OpenEMS backend.
3. If warranted, create an independently derived and publicly documented
   two-medium Green-function backend. Validate it against published theory,
   measured antennas, and licensed tools where available; do not derive it
   from restricted NEC implementations.
4. Add explicit-feedline and mixed-mode network support for common-mode work.

## Compatibility requirements

- A deck without extension cards must produce the same numerical results as
  the baseline version, within normal platform floating-point tolerance.
- Legacy `NT` and `TL` output remains unchanged unless a new verbose network
  reporting option is enabled.
- Existing command-line arguments remain valid.
- The extension is opt-in and clearly identified in `-v` output and report
  metadata.
- Unknown extension cards fail clearly; they must never be treated as comments
  or silently ignored.

## Test and validation plan

### Unit tests

- Static-Y evaluator matches legacy `NT` values exactly.
- Transmission-line evaluator matches present `TL` results at several
  frequencies, including crossed-line behavior.
- Table sorting, duplicate rejection, exact-point selection, interpolation,
  and every extrapolation policy.
- Touchstone parsing across all supported frequency units and `RI`, `MA`, and
  `DB` formats.
- S-to-Y and Z-to-Y conversion against analytic two-port fixtures.
- Reciprocity, passivity, non-finite, and conditioning failures.

### Integration tests

- A legacy deck with `NT` produces a frozen impedance/pattern fixture before
  and after the refactor.
- A synthetic known RLC two-port swept through `FR` agrees with independent
  circuit calculations at every point.
- An embedded-table network and equivalent `.s2p` provider produce the same
  result.
- A transformed load verifies the ideal 49:1 and 4:1 reference-plane rules.
- Antenna Turtle retains the deck, Touchstone sidecar, stdout/stderr, parsed
  results, and provenance after a requested analysis.

### Acceptance criteria

- No process-per-frequency workaround is required for broadband two-port
  sweeps.
- Requested frequencies outside measured coverage fail by default.
- A user can identify from a result exactly which matrix was used at any
  frequency.
- Classic decks pass regression tests unchanged.
- Measured transformer results agree with independently calculated circuit or
  VNA fixture expectations within stated measurement uncertainty.

## References

- NEC-2 network and transmission-line control cards:
  <https://www.nec2.org/part_3/control.html>
- NEC-2 transmission-line card semantics:
  <https://www.nec2.org/part_3/cards/tl.html>
- NEC-2 loading card semantics:
  <https://www.nec2.org/part_3/cards/ld.html>
- NEC-2 structure and extended-thin-wire modeling guidance:
  <https://www.nec2.org/part_3/cards/secii.html>
- LLNL's NEC-4.2 catalog entry, documenting the distinction between licensed
  later NEC capabilities and this independently maintained NEC-2 fork:
  <https://softwarelicensing.llnl.gov/product/nec-v42>
