# Canonical dielectric revision-1 examples

These small models exercise the supported interfaces without OpenEMS. From the
repository root after building `nec2c`:

```sh
python3 examples/dielectric-r1/run.py --output-dir /tmp/nec2c-r1-examples
```

Use `--binary /absolute/path/to/nec2c` for another build. The script resolves
inputs beside itself, writes logs and raw outputs under the selected output
directory, and produces `summary.json` with the solver hash and key diagnostics.
It fails if a job fails, diagnostics are missing/nonfinite, residual exceeds
1e-9, absorption is negative, or power defect reaches 0.3%. These checks exercise
examples; physical acceptance remains in the numbered milestone evidence.

| Example | Purpose | Solver | Frequency |
| --- | --- | --- | --- |
| `free-space.nec` + `.diel` | Driven wire, lossless ER=3 cube, RP far field and exterior NE/NH fields | dense | 300 MHz |
| `pec-materials.nec` + `.diel` | Two touching materials on z=0 PEC, table interpolation and conductivity/EI loss | fft | 280, 300, 320 MHz |
| `real-ground.nec` + `.diel` | Lossy cube and driven wire with positive gaps above ER=13/sigma=.02 GN 2 ground | gmres | 300 MHz |
| `plane-wave.nec` + `.diel` | Dielectric-only lossy scattering and absorption | dense | 100 MHz |

The source is 1 V for driven examples. Geometry is in meters. Materials are
arbitrary demonstrations, not measured wood, foliage or soil datasets. The
PEC sweep interpolates ER=2/sigma=.01 at 280 MHz to ER=4/sigma=.03 at 320 MHz;
the second material is constant ER=2/EI=.3. Both material regions touch PEC and
each other, while every conductor retains positive clearance. The real-ground
cube is strictly above z=0; this does not demonstrate buried/contacting geometry.

To run just the real-ground case:

```sh
./nec2c -s gmres -r 1e-10 -i examples/dielectric-r1/real-ground.nec -o /tmp/r1-real.out
```

`DF` paths are relative to the deck, not the shell's current directory. Keep
each `.nec`/`.diel` pair together when copying the examples.

## Reference observations

`reference-summary.json` records the revision-1 run and solver hash. Rounded
feed impedances and dielectric absorption are:

| Model/frequency | Input impedance (ohm) | Dielectric absorption (W) |
| --- | --- | --- |
| Free space, 300 MHz | 74.13219 + j10.03804 | 0 |
| PEC materials, 280 MHz | 58.53213 - j49.22054 | 5.58566e-8 |
| PEC materials, 300 MHz | 73.16271 + j12.80129 | 7.92480e-8 |
| PEC materials, 320 MHz | 90.85383 + j73.26262 | 3.63868e-8 |
| GN 2, 300 MHz | 73.60612 + j12.16553 | 2.68617e-7 |
| Plane wave, 100 MHz | No driven feed | 3.18494e-9 |

Small floating-point differences between builds are expected; this is a
reference record, not a bit-for-bit portability requirement. The free-space
RP/NE/NH cards execute three requests, so its summary repeats the feed/current
power result three times at the same frequency. Near-field blocks exclude the
incident wave when present; far-field values are complex rE, not a conventional
NEC gain table. See [the user guide](../../docs/DIELECTRIC-REVISION-1.md).

## Conducting-patch reference

The validated mixed wire/patch/cube reference remains at
[`tests/dielectric/m7-mixed-cube.nec`](../../tests/dielectric/m7-mixed-cube.nec),
with its adjacent material file. It exercises the conducting-patch basis and
self-consistent coupling in free space:

```sh
./nec2c -s dense -i tests/dielectric/m7-mixed-cube.nec -o /tmp/r1-mixed-patch.out
```

Its comparison/refinement history is in
[the M7 report](../../docs/DIELECTRIC-PATCH-M7.md). The revision-1 patch and
FFT tests retain this reference; its M7 filename is provenance rather than a
restriction to old solver behavior. Real-ground patches remain excluded.

## Larger tree references

The small lossy-tree material-sensitivity witness is reproducible with:

```sh
python3 tools/m13_05_tree_experiment.py --workdir /tmp/r1-small-tree --output-dir /tmp/r1-small-tree-report
```

It solves five 156-cell assumed-material cases and one 1,248-cell nominal
refinement. The sampled material range is separate from mesh sensitivity and
is not an empirical confidence interval.

The historical 9,432-cell lossless tree plus wire fixtures remain in
`tests/dielectric/m9-tree.nec`/`.diel` (free space) and
`tests/dielectric/m11-pec-tree.nec` (PEC, sharing `m9-tree.diel`):

```sh
./nec2c -s fft -r 1e-10 -i tests/dielectric/m11-pec-tree.nec -o /tmp/r1-large-pec-tree.out
```

The large model can take tens of minutes per frequency. Its accepted historical
PEC runtime target is 40 minutes, with about 240 MiB measured peak RSS on the
constrained VM. This command is optional and is not part of the quick runner.
No comparable full-size lossy/GN 2 tree resource certification is claimed.

See [deferred features and simplifying decisions](../../docs/DIELECTRIC-DEFERRED-WORK.md)
and [tree-material guidance](../../docs/M13-05-TREE-MATERIAL-GUIDANCE.md).
