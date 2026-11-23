# NEC4 Similarity M10D Results: Buried Radial-Wire Screens

## Outcome

M10D adds the public NEC4 interface condition to a balanced radial-screen
junction: one upper-medium monopole may feed three or more nearly horizontal,
azimuthally balanced lower-medium radial wires.  Every buried ray receives the
complex-permittivity endpoint factor; the already-prescribed group bypasses
the unrelated homogeneous local junction-charge solve.  Generic branched
interfaces remain rejected.

The retained sweep reproduces the important published behavior without an
empirical correction:

- `epsilon = 10-j1` has the reported finite-radius optimum and strong
  standing-wave oscillation for sparse screens;
- `epsilon = 10-j10` closely reproduces the scanned Figure 5 levels and
  ordering: the 5-wire curve settles near 3 dB, while the 20-wire curve
  settles near 4.6 dB;
- `epsilon = 10-j100` has the expected monotone approach to a roughly
  `1.5--2.3 dB` plateau with the published density ordering;
- `epsilon = 10-j1000` is nearly flat at roughly `0.5--0.8 dB`, matching the
  published plot instead of the spurious `-5 dB` seen during initial review.

## M10D.1 discrepancy diagnosis

The initial result used a NEC4-path screen field as the numerator but Burke's
published NECGS reference-dipole field as the denominator.  That is not the
report's definition of RCE, which compares both antennas at the same location
over the same ground.  The NEC4-path reference differs from the published
absolute field by an increasingly large common propagation factor:

| Ground | NEC4-path reference minus published |
|---|---:|
| `10-j1` | `-0.163 dB` |
| `10-j10` | `-0.635 dB` |
| `10-j100` | `-1.156 dB` |
| `10-j1000` | `-5.596 dB` |

Using that same-path reference denominator removes the false screen offset.
A fresh `10-j1000`, `N=20`, `0.1 lambda` screen run gives `0.66156 dB` RCE;
the retained-data re-normalization gives `0.66181 dB`, differing by only the
retained sweep's printed precision.  No antenna current, impedance, ground
kernel, or fitted scale factor was changed.

The absolute-field discrepancy remains real and is retained as separate M10B
evidence.  Both antennas enter the original-NEC4 grazing-risk region and emit
the warning.  Because the common propagation factor cancels in correctly
formed RCE, it does not invalidate Figure 5.  Improving that absolute field is
future NEC4.2-style saddle/surface-wave work, consistent with the accepted
original-NEC4 behavior boundary.

## Numerical gates

The production sweep contains all four report grounds, complete `N=5` curves,
curves through `0.4 lambda` for `N=10` and `N=20`, and `0.1 lambda` density
anchors for `N=50` and `N=100`.  Selected radial meshes were halved from a
maximum `0.02 lambda` pitch to `0.01 lambda`:

| Ground | Radials | Radius | RCE change |
|---|---:|---:|---:|
| `10-j1` | 5 | `0.2 lambda` | `-0.447 dB` |
| `10-j10` | 20 | `0.2 lambda` | `+0.450 dB` |
| `10-j10` | 100 | `0.1 lambda` | `-0.024 dB` |

The frozen refinement gate is `0.50 dB`.  These differences are invariant
under the same-path normalization correction.  The two largest changes occur at
sharp sparse/intermediate-screen curve structure; the dense anchor is much
less mesh-sensitive.  Unit tests also require exact endpoint-factor ratios,
positive accepted power, finite impedance and fields, no interface-junction
fallback, and rejection of an unbalanced branch.

## Closure: segmentation and Appendix A

The initial `0.02 lambda` pitch was sufficient for broad curve topology but
was not converged in the low-loss standing-wave region.  Three sensitive
five-radial points and the moderate-loss refinement point were recomputed at
`0.02`, `0.01`, and `0.005 lambda` maximum pitch:

| Ground, N, radius | `0.02 lambda` | `0.01 lambda` | `0.005 lambda` | final change |
|---|---:|---:|---:|---:|
| `10-j1`, 5, `0.3 lambda` | `-1.130 dB` | `+1.068 dB` | `+1.070 dB` | `+0.0017 dB` |
| `10-j1`, 5, `0.4 lambda` | `-3.431 dB` | `-1.378 dB` | `-1.359 dB` | `+0.0182 dB` |
| `10-j1`, 5, `0.5 lambda` | `-0.426 dB` | `+0.847 dB` | `+0.848 dB` | `+0.0011 dB` |
| `10-j10`, 20, `0.2 lambda` | `+4.817 dB` | `+5.267 dB` | `+5.269 dB` | `+0.0017 dB` |

The final impedance changes are at most `0.434 ohm`.  This establishes
`0.01 lambda` as the practical precision threshold for these radial screens;
halving it again changes RCE by at most `0.0183 dB`.  A `0.02 lambda` pitch can
move a sparse-screen standing-wave lobe by `1.27--2.20 dB` and its impedance
by `57--150 ohms`.  Those coarse cases are now retained deliberately as
failure witnesses rather than accepted precision results.

Nine independent impedance anchors compare `Delta Z = Zscreen -
(38.6+j22.2) ohm` with the plotted NECGS symbols in Appendix A1--A3 at
`0.1 lambda` radius.  Every resistance and reactance value falls inside a
digitization interval that includes symbol thickness and scan-reading error.
The high-conductivity anchors are particularly close; for N=20, 50, and 100
the computed resistance changes are `8.493`, `6.428`, and `5.701 ohms`, while
the digitized symbol centers are `8.2`, `6.3`, and `5.6 ohms`.  The comparison
uses the NEC symbols, not the approximate Wait-Pope solid curves.

The retained closure evidence and its numerical gates are in
`tests/nec4-similarity/m10d-closure.json`; the convergence and appendix plot
is `tests/nec4-similarity/m10d-closure.png`.  The fast regression is:

```sh
python3 tests/nec4-similarity/m10d_closure.py
```

Recollecting all closure cases is intentionally an extended validation run:

```sh
python3 tests/nec4-similarity/m10d_closure.py --collect --threads 16
```

## Reproduction

Fast retained-evidence check:

```sh
python3 tests/nec4-similarity/m10d_radial_screen.py
```

Focused M10D.1 acceptance, including a fresh high-conductivity reference and
screen solve (normally about one to two minutes with a warm ground cache):

```sh
python3 tests/nec4-similarity/m10d_radial_screen.py \
  --acceptance --threads 16
```

The final witness must report a same-path RCE in the published `0.5--0.9 dB`
band and a fresh-minus-retained difference no greater than `0.01 dB`.

Long production recollection (up to about two hours on the development
machine):

```sh
python3 tests/nec4-similarity/m10d_radial_screen.py --collect --threads 16
```

Manual representative witness:

```sh
mkdir -p tmp/m10d/manual-cache
./nec2c --wire-engine nec4 --ground-engine nec4 --model-check warn \
  --threads 16 --nec4-ground-model grid \
  --nec4-ground-cache tmp/m10d/manual-cache \
  -i tests/nec4-similarity/m10d-radial-screen.nec \
  -o /tmp/m10d-radial-screen.out
```

The witness should report one `RADIAL-SCREEN JUNCTION`, one
`INTERFACE-PRESCRIBED` junction, no fallback, and the M9 grazing warning.
The retained curve plot is
`tests/nec4-similarity/m10d-radial-screen.png`.
