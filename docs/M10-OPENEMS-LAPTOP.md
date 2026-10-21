# M10 OpenEMS staging from rxtx to the laptop

The development tools run on the rxtx cloud VM. OpenEMS runs on the user's
laptop, through the user's existing OpenEMS Python environment. No simulation
has been run on the laptop by the development agent, and no laptop runtime
estimate for the revised mesh had been measured when staging began. Returned
pilot and verification timings are reviewed below.

## Current milestone scope (2026-10-20)

M10 PEC-only validation is complete; M11 continues the lossless proof of concept with dielectric
contact against PEC ground. Real-ground jobs and their pending verification
belong to [M13](M13-LOSS-AND-REAL-GROUND.md), not to M10 acceptance. Existing
real-ground results and runner modes remain available as preliminary M13 work.
M12's conductor/dielectric contact and penetration are deferred to a later
revision and do not block M13 on separated geometries. No additional real-ground
laptop job is needed solely to close M10 or M11. The historical staging commands
below are retained; real-ground commands are optional M13 verification work.

## Battery-limited pilot

The portable script defaults to `--mode plan`, which uses only Python's
standard library and runs no solver. Start with geometry preparation, then a
single bounded timing pilot. The wall budget covers worker startup, imports,
geometry creation and simulation; termination may add up to two seconds of
grace. The revised pilot uses one CPU thread and a 120,000-step ceiling so the
Gaussian pulse is not deliberately clipped at 5,000 steps. Its external wall
budget limits battery cost; it does not compute or certify impedance. Native solver descendants are terminated
with the worker's POSIX process group on timeout or Ctrl-C.

On the laptop (replace `rxtx` with your usual SSH destination if necessary):

```sh
scp rxtx:/tmp/m10-openems-bundle.tar.gz .
tar -xzf m10-openems-bundle.tar.gz
cd m10-openems
# Use/activate your existing Python environment with openEMS and CSXCAD.
python3 openems_m10_ground.py --mode plan
python3 openems_m10_ground.py --mode prepare --seconds 30 --workdir results
python3 openems_m10_ground.py --mode pilot --seconds 30 --threads 1 --workdir results
```

`prepare` writes geometry XML and `setup.json` with the actual mesh cell count,
minimum cell sizes, maximum adjacent-cell ratios, vacuum/material-adjusted
CFL estimates and a deliberately labeled memory heuristic. These are not
the actual Rennings_2 engine timestep; the engine log remains authoritative.
It performs no FDTD time stepping. Both preparation and the pilot have a
2,000,000-cell guard; exceeding it stops before running the solver. If the
pilot times out, that is a useful cost result, not a verification failure.
Return its files and defer longer work until power is available.

```sh
tar -czf m10-openems-results.tar.gz results
scp m10-openems-results.tar.gz rxtx:/tmp/
```

Return `last-invocation.json`, every `setup.json` and the `worker.log` files even
if setup/imports/solver execution failed. They show the installed OpenEMS
version, actual grid, speed and stopping behavior when available. An absent
OpenEMS Python module usually means the existing OpenEMS environment needs
activation; no installation on the cloud VM is required.

## Verification runs, after reviewing the pilot

These commands are supplied for later use; they are not required during the
flight. Coarse invocations below have a shared three-minute wall budget across vacuum
and loaded runs; fine pairs use four minutes, reflecting the returned timings. Identical completed configurations are cached so an
interrupted pair can resume. A changed mesh/mode gets its own directory.

```sh
python3 openems_m10_ground.py --mode verify --ground pec --mesh-mm 20 --seconds 180 --workdir results
python3 openems_m10_ground.py --mode verify --ground pec --mesh-mm 15 --wire-cell-mm 2.5 --seconds 240 --workdir results
python3 openems_m10_ground.py --mode verify --ground real --mesh-mm 20 --seconds 180 --workdir results
python3 openems_m10_ground.py --mode verify --ground real --mesh-mm 15 --wire-cell-mm 2.5 --seconds 240 --workdir results
```

A run can hit its time-step cap before transients decay. Every exported result
is therefore labeled **unreviewed**. The solver energy-decay logs, unused
primitive warnings, mesh convergence, domain/PML convergence and impedance
comparison must be assessed before accepting M10 evidence. Merely returning
zero or exporting a CSV does not close an exit gate. Increase budgets only
when appropriate for the laptop; there is no assumption that 180 seconds is
enough. Pilot data will let us estimate this rather than guess.

The revised model is a vertical 475 mm, **5 mm radius (10 mm diameter)** dipole
centered at z=600 mm,
with a 10 mm lumped feed gap and an er=3 cube spanning x=[50,150], y=[-50,50],
z=[550,650] mm. Each worker writes `reference.nec`, `reference.diel` and a no-DF
`reference-vacuum-legacy.nec` with the same 5 mm physical radius and material
parameters. The NEC delta-gap source and OpenEMS finite 10 mm feed still
require a source/segmentation sensitivity check in the eventual comparison. PEC uses the z=0 boundary;
real ground fills the lower domain to z=0, epsilon_r=13, conductivity
0.02 S/m. The dielectric and vacuum runs retain identical Yee meshes.
The substrate reaches the lateral and bottom PML. Domain depth and PML error
have not yet been measured. Three frequencies (280,300,320 MHz) come from one
broadband run, rather than three separate simulations. The Gaussian bandwidth
is now 240 MHz (previously 90 MHz); the nominal pulse is about 11.94 ns rather
than 31.83 ns, without changing the three analysis frequencies.

The finite-ground NEC deck path is still disabled while coupling, fields and
power are integrated consistently. These are staged independent comparison
candidates, not claimed passing NEC/OpenEMS comparisons. Final refinement
levels may change after reviewing their actual mesh and sensitivity.

The script uses the documented [OpenEMS Python Run options](https://docs.openems.de/en/latest/concepts/simulation.html)
for thread control. Its wall timer is external to OpenEMS; the solver's
simulated-time and step limits are distinct from elapsed wall time.

## Review of the first returned prepare/pilot results

Results were supplied under `/tmp/m10-results/` on rxtx. The original grid
had 219,177 geometric cells (63 x 49 x 71), minimum spacings 1/1/5 mm and
maximum neighboring-cell ratios 3/19/6 on x/y/z. The 1 mm lines were explicit
wire/feed anchors, not accidental floating-point slivers. The smoother did
not enforce the requested grading around all the fixed anchors.

The native operator used 230,400 engine nodes (64 x 50 x 72), an actual
Rennings_2 timestep of 1.22394 ps and 20.29 seconds for 5,000 steps. Thus the
old 120,000-step ceiling projects to 486.96 seconds per run, 973.92 seconds
per pair, assuming unchanged throughput. This is a cap-cost projection, not
an estimate of when energy convergence would occur. A three-minute shared
budget could not complete a pair at the full old step limit.

The old pilot simulated 6.12 ns, while the requested Gaussian pulse required
31.83 ns / 26,007 steps. OpenEMS clipped the excitation and reached the step
cap without decay. This is timing evidence only. Its log also reports zero
active cells in the dynamic RLC extension. This is not by itself evidence of
a missing 50-ohm resistor: the [OpenEMS implementation](https://raw.githubusercontent.com/thliebig/openEMS/master/FDTD/extensions/operator_ext_lumpedRLC.cpp)
leaves a pure parallel resistor to the regular engine rather than the dynamic
RLC extension. Feed/probe behavior still needs the usual convergence review. No impedance inference is made
from the truncated pilot.

Revision 2, including the user's authorized thicker wire:

- uses radius 5 mm consistently in both model descriptions;
- retains the finite port box, enlarged to the wire radius;
- builds a deterministic graded mesh, balancing intervals rather than
  leaving short end slivers and bounding adjacent growth to 2:1;
- keeps transverse refinement separate from physical radius with
  `--wire-cell-mm`; the default is 5 mm, with 2.5 mm for the next refinement;
- avoids a forced feed-center z line and uses a shorter broadband pulse;
- uses a full-source step ceiling for the pilot, with the external wall timer
  controlling its cost;
- includes mesh revision, wire radius, refinement spacing and pulse bandwidth
  in job identities, so old pilot results cannot be silently reused.

For the revised real-ground 30 mm/5 mm pilot, the generated grid is
51 x 50 x 70 = 178,500 geometric cells, minimum spacings 5/5/10 mm,
maximum growth 2:1 and vacuum CFL estimate 11.12 ps. Actual engine timestep,
throughput and decay remain to be measured on the laptop. The thicker antenna
is a changed verification fixture, not a prediction that its impedance equals
the old thin-wire result. No existing M6/M9 antenna fixture was altered.

Inspect any returned run without OpenEMS installed:

```sh
python3 openems_m10_ground.py --analyze /path/to/returned/results
```

## Second returned batch: completed pairs and one timeout

The user supplied revision-2 results in `/tmp/m10-results-2/`. The new real
pilot finished normally in 7.63 s at a 10.472 ps timestep, a substantial
improvement over the original 1.224 ps setup. The verification jobs show:

| Ground / bulk / wire-cell spacing | Vacuum engine time | Loaded engine time | Status |
| --- | ---: | ---: | --- |
| PEC / 20 / 5 mm | 14.79 s | 14.68 s | Both completed well below step cap |
| PEC / 15 / 2.5 mm | 64.05 s | 64.17 s | Both completed, last logged decay below -64 dB |
| Real / 20 / 5 mm | 20.60 s | 20.89 s | Both completed well below step cap |
| Real / 15 / 2.5 mm | 88.33 s | incomplete | Shared 180 s invocation budget expired |

Times exclude worker imports/setup/postprocessing. The fine real vacuum job
used 91.43 s wall time; the remaining loaded job was terminated after 88.59 s
wall time. Its last progress line reached 6,216 steps and -57.71 dB, and no
`result.json` or normal engine-completion line exists. This is an incomplete
run, even though it is close to the requested -60 dB criterion. Coarse jobs
also stopped far below the 120,000-step limit without clipping/cap warnings;
normal early termination is consistent with their configured decay criterion.

At 300 MHz, the dielectric-induced impedance changes are:

| Ground | Coarse delta Z | Fine delta Z |
| --- | --- | --- |
| PEC | -0.3041 - j1.4450 ohm | -0.2953 - j1.5666 ohm |
| Real | -0.3678 - j1.4583 ohm | -0.3415 - j1.5802 ohm (**provisional**) |

The largest coarse/fine PEC dielectric-perturbation difference is 0.2631 ohm
across 280/300/320 MHz. The real fine entries above were reconstructed from
its surviving raw port traces; they are **not accepted convergence evidence**.
The independent offline DFT reproduces all seven completed jobs' native
stored impedances within 1.6e-13 ohm, retaining the current-probe half-step
phase. Removing the last 10% of the incomplete loaded trace changes its
impedance by 0.069/0.107/0.179 ohm across the three frequencies. This is a
window-sensitivity diagnostic, not an upper bound on its remaining tail error.

Absolute feed impedance remains strongly mesh-sensitive: the PEC vacuum
change between these meshes is 2.895+j16.092, 5.415+j18.099 and
9.743+j20.709 ohm at the three frequencies. Consequently the absolute feed
solution is not demonstrated converged to 2 ohms. Small perturbation changes
do not establish convergence of the baseline or equivalence of the finite
OpenEMS feed to NEC's source. The change from PEC to real ground in vacuum
is more stable: at 300 MHz it is 0.6075-j1.1624 ohm coarse and
0.6654-j1.2996 ohm fine. Domain/PML sensitivity and completed real-ground
coupled NEC comparisons remain outstanding.

The cloud-only NEC extended-kernel wire check at 300 MHz gives PEC vacuum
Z=78.646+j21.780 ohm with 11 segments, versus 79.942+j21.497,
80.417+j20.791 and 80.621+j20.031 with 21/31/41 segments. The OpenEMS fine
vacuum value is 93.610+j51.837 ohm. Wire segmentation alone does not reconcile
this discrepancy; the wire/feed discretization and finite-source equivalence
need review before accepting absolute-impedance evidence. NEC's 11-segment
real-ground legacy vacuum result is 79.085+j21.081 ohm; enabling dielectric
GN 2 wire/volume integration is now implemented. Its independently computed
300 MHz, 64-cell cloud checkpoint differs from the completed coarse OpenEMS
dielectric perturbation by 0.27573 ohm. Its power defect is 0.0876%; see
[the formulation record](M10-GROUND-FORMULATION.md) and
[`tests/m10-nec-real-300.json`](../tests/m10-nec-real-300.json). This single
checkpoint does not establish full real-ground physical acceptance.

A cloud-only, 64-cell NEC PEC solve of the returned 300 MHz benchmark
converged with GMRES. Its vacuum Z is 78.646047+j21.779558 ohm, loaded Z
78.328756+j20.610552 ohm and dielectric delta -0.317292-j1.169006 ohm.
The difference from the fine OpenEMS dielectric delta is 0.3982 ohm in
complex magnitude (0.0220 ohm resistance, 0.3976 ohm reactance). This is below
the existing 2-ohm perturbation allowance at this one checkpoint. It does not
establish absolute-impedance convergence, VIE volume refinement, domain
convergence, all frequencies or real-ground coupled solver agreement. Run
metadata and values are retained in `tests/m10-nec-pec-300.json`.

Normalized review data are versioned as `tests/m10-openems-ground-pairs.csv`
and `tests/m10-openems-ground-review.json`, with incomplete estimates labeled
explicitly. Reproduce the port/log review with:

```sh
python3 tools/m10_review_results.py /tmp/m10-results-2 --output-dir /tmp/m10-review
```

No further battery-powered refinements are requested now. The missing real-
ground run is now M13 evidence and does not block M10 or M11. If completing it
is convenient later, rerun **on the laptop in the original
working folder**, using the same existing `results` directory:

```sh
python3 openems_m10_ground.py --mode verify --ground real --mesh-mm 15 --wire-cell-mm 2.5 --seconds 120 --workdir results
```

The exact completed vacuum job will be reused. Only the loaded job restarts
from the beginning (approximately 90 s on this laptop), not from its last
FDTD step. Retain the archived partial traces already on rxtx. This rerun
closes the missing M13 decay result, not the remaining source/domain/solver
validation gates. It is not an M10 prerequisite.

### PEC feed and domain sensitivity (2026-10-20)

Two additional vacuum/loaded pairs isolate remaining source and boundary choices.
Keep the baseline 20 mm bulk grid, 5 mm near-wire spacing and 5 mm wire radius:

```sh
python3 openems_m10_ground.py --mode verify --ground pec --mesh-mm 20 --wire-cell-mm 5 --domain-margin-mm 200 --seconds 120 --workdir pec-domain-check
python3 openems_m10_ground.py --mode verify --ground pec --mesh-mm 20 --wire-cell-mm 5 --feed-gap-mm 5 --seconds 120 --workdir pec-feed-check
```

The domain variant extends the side boundaries from ±650 to ±850 mm and
raises the top from 1400 to 1600 mm; the PEC plane stays at z=0. The feed
variant shortens the finite source from 10 to 5 mm while preserving outer
wire endpoints, radius, and dielectric geometry. Both sample 280/300/320 MHz.
The budget covers both members of a pair; rerun an interrupted command to
reuse its completed member. Normal decay must still be checked before
acceptance. Default configuration hashes remain unchanged. Variant physics
are included in job identities and offline pairing to prevent cross-pairing.

The returned sensitivity archive was accepted on 2026-10-20. All four runs
reached -60 dB decay normally. Maximum changes in dielectric-induced impedance
were 0.00574 ohm for the larger domain and 0.00630 ohm for the shorter feed.
The complete audit is `tests/m10-pec-sensitivity-review.json`; paired results
are in `tests/m10-pec-sensitivity-pairs.csv`. These close the remaining
PEC-only M10 acceptance gates. No additional laptop runs are required for M10.
