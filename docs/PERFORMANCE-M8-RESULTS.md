# Milestone 8: reuse, batching, caching, and field outputs

Status: implementation complete; awaiting user acceptance.

## Delivered behavior

Ordinary classic RP samples are evaluated into fixed indexed buffers with
OpenMP and are then consumed by the original serial gain, normalization,
power-integration, report, plot, and JSON loops. This preserves source
summation within each angle and preserves all externally visible ordering.
Ground-wave `RP` mode remains serial because its legacy evaluator mutates
shared scratch state. Automatic classic field evaluation is conservatively
capped at four threads; an explicit `--threads N` request remains authoritative.

Volume-only dielectric far fields, scattered-power quadrature, and exterior
near-field grids use the same indexed-compute/serial-consume design. Coupled
conductor-volume far fields and power quadrature are also parallel where the
field evaluator is reentrant. Coupled near fields remain serial because NEC's
legacy `nefld` and `nhfld` routines use shared scratch storage. Small work is
kept serial, and every allocation or overflow failure falls back to the
preserved serial path.

The Sommerfeld interpolation grid now has a one-entry in-process exact-key
cache. Its key is the exact complex relative permittivity, which completely
defines the wavelength-normalized table. Repeating a frequency/material state
is a hit; any bitwise frequency/material change that changes the complex
permittivity is a miss. Timing JSON records hit and miss counts.

The classic matrix factorization and solved currents were already retained
through an excitation/output-card sequence. M8 completes the analogous network
lifecycle: NT/TL preparation, reduced factorization, pivots, equation maps,
and source state persist for the lifetime of the primary factorization.
Independent unit endpoint right-hand sides are submitted through the existing
multi-RHS solver in one batch, while `cabc` conversion and all output remain in
their original order. Rebuild is forced when `ntsol` is reset or cached
dimensions/topology counts differ. This also repairs the previous lifecycle in
which `ntsol` requested reuse after its local prepared arrays had been freed.

All caches are process-local and bounded to one current state. No persistent
on-disk cache was introduced; the roadmap still requires a separate design
review before any persistent cache is considered.

## Numerical and safety coverage

`tests/field_parallel.py` requires exact normalized report and result JSON for:

* a 101-segment classic model comparing serial with the automatic parallel
  field plan;
* a 216-cell volume-only model comparing one and four workers; and
* an identical-frequency real-ground sweep with exactly one Sommerfeld cache
  miss and one hit, plus a changed-frequency sweep with two misses.

`tests/network_cache.py` exercises an NT network over a 3 by 2 plane-wave
angle sweep. It requires six ordered network solutions, one preparation miss,
five hits, and exact output across one and four threads. A direct focused
built-in/LAPACK comparison differed by at most `6.94e-17`. Existing coupled
free-space, connected-patch, PEC, real-ground, FFT backend, legacy regression,
and backend-equivalence tests remain unchanged and pass.

The complete Automake suite passes 36 of 36 groups. A distribution build,
test, install, uninstall, and archive check passes. An independent build with
OpenMP and LAPACK disabled compiles without warnings and runs the preserved
serial field path. AddressSanitizer plus UndefinedBehaviorSanitizer passes the
focused field and network cache regressions (leak detection disabled because
it cannot initialize under the host tracing mode).

The clean-room AEG review found no ready implementation to copy. It confirmed
the architectural boundary used here: immutable geometry/current inputs,
point-local field outputs, and serial ownership of reports and reductions. No
AEG source or expression was used.

## Development performance evidence

Three measured repetitions after one warmup produced the following focused
results on the 32-logical-CPU acceptance host:

| Case | 1 thread | 16 threads | Phase speedup | Overall speedup |
|---|---:|---:|---:|---:|
| 801-segment, 16,471-point classic RP | 0.238311 s | 0.041265 s | 5.78x | 3.26x |
| Coupled angular power quadrature | 0.022610 s | 0.003546 s | 6.38x | 3.46x complete solve |

The original Antenna Turtle loop Yagi also improved in a focused five-repeat
development comparison, though its roughly 20--60 ms field phase is much more
sensitive to CPU frequency and desktop scheduling. The larger 801-segment
fixture is therefore the acceptance timing gate; the loop Yagi remains the
recommended application-level visual model.

## User acceptance

Run the deterministic regressions first:

```sh
python3 tests/field_parallel.py
python3 tests/network_cache.py
python3 tests/coupling_parallel.py
```

Expected final lines are:

```text
classic and volume-only far fields preserve exact ordered output under parallel evaluation; exact Sommerfeld setup reuse passed
prepared network state is exact across six ordered plane-wave excitations
free-space, connected patch, PEC, and real-ground coupling is bit-identical at 1, 2, and 4 workers; automatic FFT is exact
```

Then run the seven-repeat scaling study:

```sh
python3 tools/m8_field_benchmark.py \
  --threads 1,2,4,8,16,32 \
  --warmup 2 --repeat 7 \
  --output /tmp/m8-field.json
```

The benchmark fails if normalized numerical differences exceed `1e-10`, the
16-thread classic field-compute speedup is below 1.5x, the coupled power phase
is below 1.2x, or the coupled complete solve is below 1.2x. It records every
sample, phase and overall medians, peak RSS, numerical deltas, and validation
decision.

For a manual application-level pattern, run the Antenna Turtle loop Yagi at
one thread and with the default automatic policy:

```sh
python3 tools/performance_benchmark.py \
  --case antenna-turtle-loop-yagi --backend auto \
  --jobs 1 --threads 1 --repeat 7 --warmup 2 \
  --output /tmp/m8-yagi-serial.json

python3 tools/performance_benchmark.py \
  --case antenna-turtle-loop-yagi --backend auto \
  --jobs 1 --threads auto --repeat 7 --warmup 2 \
  --output /tmp/m8-yagi-auto.json
```

No additional software or library installation is required for this milestone.
