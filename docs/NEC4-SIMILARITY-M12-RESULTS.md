# NEC4ish similarity M12: release candidate

## Outcome

M12 was accepted with a user PASS on 2026-11-24.  Classic NEC2 remains the
default with no flags.  NEC4ish requires explicit engine selection, while
`--wire-engine auto --ground-engine auto` is an explicit opt-in to the
qualified selection heuristic.  The production NEC4ish ground default remains
`grid`; the very expensive cold `adaptive` cache is an explicit accuracy-first
choice.

## Release gates

- The preserved M0 NEC2 corpus is unchanged under default and explicit selection.
- Automatic and corresponding explicit engines produce equal projected
  machine results for free-space, buried real-ground, patch, coated-wire, and
  dielectric witnesses.
- All historical coated-wire and dielectric null-limit suites now request
  NEC2 explicitly, making their intended comparison stable.
- Automatic choice is exposed in both the human report and job-level JSON.
- Unsupported explicit combinations fail rather than silently falling back.

## Sanitizers and concurrency tools

An isolated GCC AddressSanitizer/UndefinedBehaviorSanitizer build passed the M0,
M4, M7, M9, M11E, M12, result-JSON, and frequency-parallel gates.  It found one
real defect: a patch-only model passed a null zero-length endpoint array to
`memset`.  The operation is now guarded and the same sanitizer portfolio passes.

GCC ThreadSanitizer cannot start on this host because its shadow-memory mapping
conflicts with the runtime environment.  Valgrind Helgrind reports libgomp
barrier internals as conflicts and is not a reliable OpenMP race oracle here.
Concurrency acceptance therefore uses exact ordered 1/2/4-worker identity plus
the sanitizer and existing parallel regression suites.  This tool limitation
is disclosed rather than represented as a clean race-tool result.

## End-to-end cost of the accuracy default

Sequential medians use nine measured runs after two warmups:

| Case | Opt-in auto NEC4ish | Explicit NEC2 | Time ratio | Peak RSS change |
|---|---:|---:|---:|---:|
| Tiny near-field witness | 0.0360 s | 0.00070 s | 51.4x | +1.8 MiB |
| Antenna Turtle loop Yagi | 0.0798 s | 0.0283 s | 2.82x | +31.2 MiB |
| 801-segment scaling case | 0.1280 s | 0.0605 s | 2.11x | +9.8% |
| 401-point real-ground sweep | 0.3436 s | 0.1518 s | 2.26x | +5.4% |

The tiny-case ratio is dominated by fixed NEC4ish work; its absolute cost is
36 ms.  The larger cases characterize the practical accuracy cost.  Existing
frequency parallelism remains active and previously demonstrated about 5.4x
on the retained sweep witness.

M11E separately measured warm adaptive versus direct at 2.10x end-to-end and
205.7x for the deliberately difficult near-field phase.  Cold adaptive cache
generation took 648.9 s and 79 MiB, which is why M12 leaves `grid` as the
automatic ground evaluator.

## Provenance

The implementation derives from the public NEC4 theory guide, public Burke
papers/reports, the public NEC5 validation manual as independent context, and
independent analytic/numerical oracles.  No NEC4 binary or source was acquired
or used.  The public name is always **NEC4ish**; `nec4` survives only as an
input compatibility alias.

## Manual acceptance

Build and run the compact release gate:

```sh
make -j16
python3 tests/nec4-similarity/m12_release.py
python3 tests/regression.py
```

Exercise the retained geometry/ground portfolio:

```sh
python3 tests/nec4-similarity/m6_combined.py
python3 tests/nec4-similarity/m9_production.py
python3 tests/nec4-similarity/m10c_ground_stake.py
python3 tests/nec4-similarity/m10d_radial_screen.py
python3 tests/nec4-similarity/m11e_production.py
```

Optional full suite:

```sh
make check
```

For manual NEC report inspection, compare automatic and preserved behavior:

```sh
./nec2c --wire-engine auto --ground-engine auto \
  -i tests/nec4-similarity/m7-buried-dipole.nec -o /tmp/m12-auto.out
./nec2c -i tests/nec4-similarity/m7-buried-dipole.nec -o /tmp/m12-default.out
./nec2c --wire-engine nec2 --ground-engine nec2 \
  -i tests/nec4-similarity/m7-buried-dipole.nec -o /tmp/m12-nec2.out
```

M12 received its required user PASS.  The feature branch may now be committed,
pushed, and merged.
