# M13-01: constant isotropic dielectric loss

Status: COMPLETE (2026-10-21). Native analytical, refinement, independent
OpenEMS and regression gates pass. M13 as a whole remains in progress.

This step adds passive, scalar, nonmagnetic material loss in free space and
above or touching PEC ground, with the existing conductor clearance rules.
Anisotropy, dispersion tables, additional real-ground modes and material
interface extensions are later work. An isotropic approximation cannot resolve
wood grain direction; selecting its effective scalar properties remains a
modeling assumption.

## Resource version 2

Version 1 remains supported unchanged. Version 2 accepts either:

```
NEC2C-DIELECTRIC 2
MATERIAL 1 ER 3 SIGMA 0.02
REGION 1 MATERIAL 1
BOX 1 .05 -.05 0 .15 .05 .65
MESH MAX_KH .5 MAX_ASPECT 2 MAX_CELLS 10000
END
```

or `MATERIAL 1 ER 3 EI 0.4`. `ER` is the real relative permittivity,
`EI` is the nonnegative magnitude of its negative imaginary part, and `SIGMA`
is conductivity in S/m. Either loss term may be omitted (zero loss). Specifying
both is rejected to avoid ambiguous input. ER remains finite and >=1;
negative or nonfinite loss is rejected. Version 1 rejects loss fields.

With NEC's exp(+i omega t) convention:

```
epsilon_relative = ER - i*(EI + SIGMA/(omega*epsilon0))
J = i*omega*epsilon0*(epsilon_relative - 1)*E
P_abs = sum_cells volume*omega*epsilon0*(-Im(epsilon_relative))*|E|^2/2
```

The cell-averaged E is recovered from J using the constitutive relation.
Absorption is included in extinction/accepted-power balances alongside
radiation and, where applicable, lower-half-space flux. Material absorption
is reported separately; it is not inferred from the global power defect.
Mesh electrical size uses `k0*sqrt(abs(epsilon_relative))*h`. ER=1 with
nonzero loss remains active material. Conductivity is reevaluated at every
frequency. FFT retains its uniform-material and cubic-lattice restrictions,
now allowing a uniform complex scalar permittivity.

## Development and acceptance

`tests/m13_loss.py` checks analytic electrically-small cube polarization,
independent Joule absorption, passive signs, EI/conductivity equivalence,
frequency-sweep state, zero-loss limits, conductive ER=1, invalid inputs,
dense/GMRES/FFT agreement and coupled power conservation over free space/PEC.
Existing regression groups must continue to pass.

`tools/m13_01_acceptance.py` produces cached native counterparts to the laptop
pillar, with exact binary/deck/material provenance. It supports 280/300/320 MHz,
25 mm and 12.5 mm volume meshes, and wire refinement. Physical comparisons use
each induced delta R/delta X component within max(2 ohm, 2% of the corresponding
OpenEMS reference component). Complex-vector discrepancy is diagnostic only.

## OpenEMS experiment M13-01

`tools/run-m13-01-loss.sh` runs the existing M11 PEC-contact pillar with ER=3
and conductivity=0.02 S/m. It stages a loaded pilot, coarse/fine vacuum/lossy
pairs and a fine pair with 150 mm extra side/top domain clearance. The shared
runner records conductivity in cache identities and supplies CSXCAD kappa in
S/m. Legacy zero-conductivity cache identities remain unchanged.

A self-contained download is staged at `/tmp/m13-01-openems.tar.gz` on rxtx.
Unpack it on the laptop and run `bash run-m13-01-loss.sh` in the existing
OpenEMS Python environment. It packages `m13-01-results.tar.gz` automatically,
including partial results on failure. Upload that archive to `/tmp` on rxtx.
The package README contains the same instructions; no session-only arguments
are necessary. Simulations run independently of native development. Returned
results gate acceptance of M13-01, not implementation work during this step.

### Patch discretization diagnostic

An exploratory 300 MHz free-space closed-patch cube at grazing plane-wave
incidence (`EX 1 ... 90 0 90`) gives power defects 0.6121%, 0.4352%, and
0.3266% for 2x2, 3x3, and 4x4 patches per face without dielectric loss. With
sigma=0.02 S/m in the small contact fixture, the corresponding defects are
0.6114%, 0.4346%, and 0.3261%. This inherited surface-discretization effect
converges with patch refinement and is not used to relax loss acceptance.
The oblique 37/23/19-degree regression gives 0.0248% free-space and 0.0139%
PEC power defects with loss. All step regression tolerances remain unchanged.

## Native pillar evidence (2026-10-21)

All six cached vacuum/lossy pairs pass true residual <=1e-10, positive
absorption and <0.3% total power defect. Retained values and binary hashes are
in [m13-01-acceptance.json](../tests/m13-01-acceptance.json).

| MHz | Cube pitch | Wire segments | Induced delta Z (ohm) | Absorbed power |
| ---: | ---: | ---: | ---: | ---: |
| 280 | 25 mm | 11 | 2.50837 + j2.29049 | 0.59833 mW |
| 300 | 25 mm | 11 | 2.54731 + j0.54561 | 0.43452 mW |
| 320 | 25 mm | 11 | 2.32949 - j1.87403 | 0.26651 mW |
| 300 | 12.5 mm | 11 | 2.61877 + j0.64390 | 0.43900 mW |
| 300 | 25 mm | 21 | 2.59301 + j0.41629 | 0.42913 mW |
| 320 | 25 mm | 21 | 2.29654 - j2.11294 | 0.26328 mW |

Volume refinement changes induced impedance by 0.12152 ohm; wire refinement
changes it by 0.13716 ohm. The 3328-cell case takes 558.02 s wall / 40.43 s
solver CPU and peaks at 28,948 KiB RSS. Refined-wire power defect is 0.01724% at 300 MHz.
These are native consistency results; the independent OpenEMS comparison
below supplies the physical acceptance check. The pinned benchmark predates only the lossless
zero-frequency guard and output-label clarification, which do not change the
positive-frequency equations or power calculations.

## Returned OpenEMS experiment M13-01

The returned archive `/tmp/m13-01-results.tar.gz` contains a pilot and six
verification runs. All six explicitly reach <=-60 dB decay without source
clipping or a step cap. Loaded XML confirms isotropic ER=3 and electric
conductivity kappa=0.02 S/m. Meshes retain >=2.5 mm minimum cells, <=2:1 grading,
and fewer than two million cells. Offline port DFT reproduces saved impedances.

| MHz | Fine, enlarged-domain OpenEMS delta Z (ohm) | Mesh refinement change | Domain enlargement change |
| ---: | ---: | ---: | ---: |
| 280 | 2.98857 + j1.96099 | 0.21851 ohm | 0.01378 ohm |
| 300 | 3.07830 - j0.16063 | 0.20595 ohm | 0.02801 ohm |
| 320 | 2.59591 - j3.41560 | 0.41371 ohm | 0.01838 ohm |

All six native cases pass component-wise max(2 ohm, 2% of reference component).
The largest component error is 1.54158 ohm for the 11-segment 320 MHz case;
21-segment refinement reduces its resistance/reactance errors to
0.29937/1.30267 ohm. No additional OpenEMS runs are needed for this step.
Archive hash, explicit decay, mesh audit, native comparisons and convergence
are retained in [m13-01-openems-review.json](../tests/m13-01-openems-review.json)
and [m13-01-openems-pairs.csv](../tests/m13-01-openems-pairs.csv).

## Step closure

All 20 regression groups pass. The initial run exposed a stale FFT test object
from a bootstrap dummy dependency file when the FFT API changed to complex
permittivity. Explicit header dependencies now rebuild both API consumers;
the rebuilt operator test also verifies complex-material products and solves
against dense references and rejects active permittivity. The final aggregate
summary is 20/20 PASS, with source distribution and whitespace checks passing.

M13-01 is complete. This accepts constant passive isotropic loss in free space
and over PEC, retaining existing conductor-clearance restrictions. It does not
complete dispersion tables, touching multi-material extensions or the remaining
real-ground acceptance work. Anisotropy remains excluded. No further OpenEMS
execution is required for M13-01.
