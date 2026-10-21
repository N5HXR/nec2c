# M13-02: frequency-dependent scalar material tables

Status: COMPLETE (2026-10-21).
Implementation, native analytical/refinement checks, all 21 regression groups
and independent OpenEMS witness acceptance pass. M13 overall remains in progress.

## Resource version 3

Versions 1 and 2 remain supported unchanged. Version 3 permits constant
materials with the existing ER/EI/SIGMA syntax, or sampled materials:

```
NEC2C-DIELECTRIC 3
MATERIAL 1 TABLE
SAMPLE 1 FREQ_MHZ 260 ER 4 SIGMA .01
SAMPLE 1 FREQ_MHZ 300 ER 3 SIGMA .02
SAMPLE 1 FREQ_MHZ 340 ER 2 SIGMA .03
REGION 1 MATERIAL 1
BOX 1 .05 -.05 0 .15 .05 .65
MESH MAX_KH .5 MAX_ASPECT 2 MAX_CELLS 10000
END
```

Each TABLE requires at least two distinct frequencies. Rows may be unordered
and may precede their material declaration; the loader sorts by material and
frequency. Orphan rows, rows attached to a constant material, duplicate
frequencies, and missing rows are rejected. All rows in one table must use
the same loss representation: EI or SIGMA. Frequencies are finite positive
MHz, ER is finite and >=1, and loss is finite and nonnegative.

A table sample must differ from vacuum (ER=1 and zero loss is rejected).
This preserves a nonzero constitutive susceptibility throughout interpolation;
use constant ER=1 for an omitted vacuum region. Lossless tables with ER>1 and
conductive tables with ER=1 and positive loss are supported. Anisotropic and
magnetic media remain outside this increment.

Between consecutive rows, ER and the chosen loss parameter are interpolated
linearly in ordinary frequency, not log frequency. With exp(+i omega t):

```
epsilon_relative(f) = ER(f) - i*EI(f)
# or
epsilon_relative(f) = ER(f) - i*SIGMA(f)/(2*pi*f*epsilon0)
```

Interpolating conductivity therefore differs from interpolating imaginary
permittivity. No extrapolation is permitted. Every declared table must cover
every requested FR frequency; the whole additive/geometric, ascending/descending
sweep is checked before any solution from that sweep is produced. Only floating-point endpoint
roundoff is clamped. Source samples are immutable, and each frequency evaluates
its own material value without retaining interpolated state. Later FR cards
invalidate the earlier mesh and rerun coverage checks; NX releases table data.

These tables specify steady-state frequency-domain values. They do not fit or
certify a causal broadband Debye/Lorentz constitutive law.

## Fixed sweep mesh

Meshing cannot inspect only the permittivity at the highest sweep frequency:
a table may have a larger interior value. The implementation bounds ER and
the imaginary-permittivity magnitude separately across all table rows, then
uses `sqrt(hypot(max_ER,max_imaginary_magnitude))` as the material index bound.
For EI tables each component is piecewise linear. For SIGMA tables, linear
SIGMA(f)/f is monotone within each positive-frequency interval, so its endpoint
values bound the loss component. Multiplying by the highest requested k0
bounds electrical cell size throughout the sweep.

This table-wide bound is conservative and can include rows outside the
requested sweep. The mesh remains fixed across frequencies; unused distant
high-valued rows may increase cell count. Existing mesh, conductor clearance,
FFT lattice/uniform-material restrictions and the 512 MiB solver cap remain.

## Native verification

`tests/m13_tables.py` covers knots and off-knot values, EI and SIGMA tables,
independent small-cube polarizability, absorption and zero-loss behavior,
dense/GMRES/FFT comparisons, ascending/descending/geometric sweeps, preflight
extrapolation rejection, interior high-permittivity mesh peaks, canonical
ordering/forward references, malformed tables, coupled wire/PEC power, and
NX table-data cleanup/state isolation.

`tools/m13_02_acceptance.py` runs cached native counterparts of the M13-02
pillar with the table above. Its mesh setting preserves the selected cubic
pitch using the conservative table bound. Binary/deck/material identities,
residual, power, impedance and resources are retained per configuration.

## OpenEMS experiment M13-02

A self-contained package is staged at `/tmp/m13-02-openems.tar.gz`. In the
existing laptop OpenEMS environment, unpack and run:

```
bash run-m13-02-tables.sh
```

Return `m13-02-results.tar.gz` to `/tmp` on rxtx. The script includes all
arguments, automatically packages partial results, and resumes completed cases.
It runs 280 MHz (ER=3.5, SIGMA=.015) and 320 MHz (ER=2.5, SIGMA=.025) witnesses,
each with a pilot, coarse/fine vacuum-loaded pairs and a larger-domain fine
pair. The 300 MHz knot reuses M13-01's ER=3, SIGMA=.02 evidence.

Each FDTD configuration uses the constant material matching its designated
witness frequency. Other reported frequency rows do not compare with the
interpolated table and must not be promoted to table acceptance. The reviewer
retains witness identities to keep physically different comparisons separate.
Accept only normal explicit decay, independently refined meshes/domains and
component-wise induced delta R/delta X errors within max(2 ohm,2% of the
corresponding reference component). Development continues while runs are
outstanding; returned evidence gates step acceptance.

Retained native evidence and completed acceptance checks are recorded in
[m13-02-acceptance.json](../tests/m13-02-acceptance.json). The 300 MHz table knot
exactly reproduces the retained M13-01 constant-material induced impedance.

## Native pillar checkpoint

The 416-cell 280/300/320 MHz table cases and 320 MHz volume/wire refinements
all pass residual <=1e-10, nonnegative absorption and <0.3% power defect.
At 320 MHz, 12.5 mm cubes give 3328 cells, delta Z=1.96120-j1.17894 ohm,
0.34066 mW absorption and 0.18437% power defect. The run takes 504.55 s wall
and 38.12 s solver CPU, peaking at 28,272 KiB RSS. Volume refinement changes
induced impedance by 0.11259 ohm. With 21 wire segments and 25 mm cubes,
delta Z=1.85381-j1.43644 ohm and power defect is 0.00903%; wire refinement
changes induced impedance by 0.18813 ohm. These native consistency results
do not replace the independent off-knot OpenEMS witnesses.

## Native development closure

All 21 aggregate regression groups pass, including the prior M13-01 loss and
M11 contact coverage. Source distribution and whitespace checks pass. The
independent M13-02 witness decay/material/mesh/domain audit and impedance
comparison also pass. No later M13 increment is started by this closure.

## Independent acceptance closure (2026-10-21)

All 12 returned verification runs explicitly reach -60 dB decay (final values
-65.05 to -60.04 dB), before the step cap, without source clipping. XML checks
confirm the isotropic witness permittivity/conductivity. All 14 returned meshes,
including the pilots, have minimum steps at least 2.5 mm, adjacent step ratios
at most 2 (within rounding), and fewer than two million cells.

Only each run's designated witness frequency is compared with the table:
280 MHz uses ER=3.5, sigma=0.015 S/m; 320 MHz uses ER=2.5, sigma=0.025 S/m.
Other spectral rows are retained for audit and are not table witnesses. The
300 MHz knot reuses the accepted M13-01 ER=3, sigma=0.02 S/m reference.
All five native benchmark/refinement cases satisfy the component-wise
max(2 ohm, 2% of the corresponding reference component) rule. The largest
component discrepancy is 1.33655 ohm, in refined-volume reactance at 320 MHz.
The complex-vector discrepancy is diagnostic only.

OpenEMS mesh refinement changes induced impedance by 0.20946 ohm at 280 MHz
and 0.33993 ohm at 320 MHz. Enlarging the fine-mesh domain changes it by
0.02750 and 0.02612 ohm, respectively. Retained results, archive hash,
material/decay/mesh audits and comparisons are in
[m13-02-openems-review.json](../tests/m13-02-openems-review.json), with spectral
pairs in [m13-02-openems-pairs.csv](../tests/m13-02-openems-pairs.csv).
These steady-state witnesses validate table evaluation at sampled frequencies;
they do not establish a causal broadband material fit.
