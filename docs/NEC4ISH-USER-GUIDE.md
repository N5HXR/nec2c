# NEC4ish user guide

## Identity and scope

NEC4ish is an independent implementation based on publicly documented theory.
It contains no NEC4 source code and was not developed against a licensed NEC4
executable.  The name is intentionally displayed as **NEC4ish** everywhere.
The old spelling `nec4` is accepted only as a command-line compatibility alias.

NEC4ish improves the wire basis, junction/end treatment, and homogeneous
half-space interaction used for wires at, crossing, or below a planar ground
interface.  It is not a claim of byte-for-byte or universal numerical identity
with any version of NEC4.

## Default and compatibility controls

Both selectors default to `nec2`:

```text
--wire-engine auto|nec2|nec4ish
--ground-engine auto|nec2|nec4ish
```

An invocation with no engine flags therefore retains classic NEC2 behavior.
NEC4ish is deliberately opt-in while the implementation accumulates broader
field experience.  Explicit `auto` uses the qualified NEC4ish wire operator
for ordinary wire models.  It
uses the NEC4ish ground operator only for finite Sommerfeld ground with a wire
endpoint below the interface.  Otherwise it retains the faster NEC2 ground
path.  The automatic decision is printed in the report and recorded in the
machine-result JSON.

The default and its fully explicit equivalent are:

```sh
nec2c --wire-engine nec2 --ground-engine nec2 -i model.nec -o model.out
```

Opt into the selective heuristic with:

```sh
nec2c --wire-engine auto --ground-engine auto -i model.nec -o model.out
```

Unsupported combinations do not silently mix unqualified physics:

| Model feature under explicit `auto` | Automatic wire | Automatic ground | Reason |
|---|---|---|---|
| Ordinary wires, free space/PEC/classic ground | NEC4ish | NEC2 | Classic ground is already qualified and faster |
| Ordinary wires wholly above finite Sommerfeld ground | NEC4ish | NEC2 | Above-ground optimized path |
| Wire endpoint below finite Sommerfeld interface | NEC4ish | NEC4ish | Documented interface formulation |
| Surface patches | NEC4ish | NEC2 | NEC4ish patch/ground coupling is not qualified |
| Coated wires | NEC2 | NEC2 | Combined coated-wire/NEC4ish basis is not qualified |
| Dielectric VIE | NEC2 | NEC2 | Combined VIE/NEC4ish basis is not qualified |

An explicit unsupported request fails with an actionable diagnostic.  It does
not fall back behind the user's back.

## Ground evaluators

When the NEC4ish ground engine is active:

```text
--nec4ish-ground-model grid|direct|adaptive
--nec4ish-ground-cache DIRECTORY|none
```

- `grid` is the production default.  It is the best balance for sweeps and
  ordinary antenna work.
- `direct` evaluates the documented contour integrals and is the numerical
  oracle.  It is intentionally expensive.
- `adaptive` is an independently qualified accuracy-first extension for the
  difficult grazing/interface domain.  A cold cache may take about eleven
  minutes and 79 MiB for the hardest retained material; warm runs are fast.

The original documented NEC4 grazing approximation can be badly wrong in its
risk region.  NEC4ish detects the known geometry and issues a warning.  The
`grid` path preserves the documented behavior rather than applying an
unpublished empirical correction.  Select `adaptive` or `direct` when the
warning matters to the analysis.

## Validity diagnostics

```text
--model-check warn|strict|off
```

- `warn` is the default: report validated-envelope concerns and continue.
- `strict` turns qualified validity violations into errors.
- `off` suppresses advisory warnings, but cannot make undefined geometry legal.

Guidance is labelled by evidence class:

- **theory-mandated**: required by the published formulation or a mathematical
  limit;
- **empirically qualified**: bounded by retained convergence/oracle tests;
- **practitioner-advisory**: useful modeling guidance without an independently
  established failure boundary.

## Machine-result additions

Each job contains an `engines` object:

```json
{
  "wire": {"requested": "auto", "resolved": "nec4ish"},
  "ground": {"requested": "auto", "resolved": "nec2"},
  "selection_reason": "classic-ground-mode"
}
```

Automatic requests also emit one `ENGINE_AUTO_SELECTION` informational
diagnostic with the same values.  Explicit engine requests retain the metadata
but do not add this informational diagnostic or an automatic-selection line to
the human report.

## Performance expectations

The NEC4ish basis is an accuracy feature, not a blanket speed feature.  On the
M12 host it cost about 2.1-2.8x wall time on representative 401/801-segment
wire cases and about 5-10% memory on the larger cases.  Existing frequency,
fill, field, FFT, LAPACK, and coupling parallelism remains active.  Small jobs
can show a large ratio while still taking only a few hundredths of a second.

For maximum legacy speed:

```sh
nec2c --wire-engine nec2 --ground-engine nec2 ...
```

For accuracy-first difficult interface work:

```sh
nec2c --wire-engine nec4ish --ground-engine nec4ish \
  --nec4ish-ground-model adaptive ...
```

## Known limits

- One planar homogeneous half-space; no terrain or multilayer soil.
- Coated-wire and dielectric-volume combinations use NEC2 automatically.
- Patch interaction with NEC4ish ground is not qualified.
- The `grid` evaluator intentionally retains the documented original grazing
  limitation and warns when detected.
- Published comparisons and independent equation-level oracles replace a
  licensed NEC4 executable; claims are capability-specific.
