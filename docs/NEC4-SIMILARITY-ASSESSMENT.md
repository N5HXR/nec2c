# NEC4/5 Similarity: Technical Assessment

Status: research draft for objective-setting; no implementation plan or code change is authorized by this document.

## Source and interpretation

The primary source is Gerald J. Burke, *Numerical Electromagnetics Code - NEC-4, Method of Moments, Part II: Program Description - Theory*, UCRL-MA-109338 Pt. II, January 1992 (the local `nec4-theory.pdf`). It is a 150-page scan containing the theory volume, including the equation-level wire, junction, end-cap, and ground formulations and appendices for their numerical evaluation. Three later primary sources refine the assessment:

- Burke, *Accuracy of Reduced and Extended Thin-Wire Kernels*, LLNL-PROC-409033, November 25, 2008 (`reports/accuracy.pdf`), defines and compares exact, extended and reduced kernels and identifies NEC4 as a point-matched EFIE code using the reduced thin-wire kernel.
- Burke, *Ground Model Options in the NEC-4.2 Antenna Code*, LLNL-TR-490316, July 6, 2011 (`reports/gn3.pdf`), describes the later GN3 ground option and briefly notes EX6 current sources and dynamic array allocation.
- Burke, *NEC-5 Validation Manual*, LLNL-SM-791163, September 24, 2019 (`reports/NEC5 Validation Manual 092419.pdf`), compares NEC4 with the independently formulated mixed-potential NEC5 and documents both agreement and known NEC4 failure modes.

The manual describes the original NEC-4 theory. It does not by itself establish exact behavioral equivalence with every later NEC 4.1 or 4.2 implementation detail. LLNL still distributes NEC 4.2 under a separate end-user license. Consequently, this project should use capability-specific language such as **NEC4-similarity wire engine** or **documented NEC4 wire formulation**, and should not claim to be NEC 4.2 without independent validation and a separate decision about naming.

The manual is sufficient to implement a substantial, independently derived numerical method. It is not sufficient by itself to reproduce every NEC4-family feature or user-interface behavior:

- the core wire basis, local junction solve, free-space field kernel, low-frequency approximations, and flat end-cap model are specified mathematically;
- the NEC3-derived buried/cross-interface half-space formulation is also described, but is a much larger subsystem;
- the manual says insulated wire is included, but does not provide a self-contained insulated-wire derivation comparable to the other major algorithms;
- input-card details and some convenience features belong in the separate User's Guide;
- later NEC 4.2 maintenance changes are not necessarily represented in this 1992 document.

### Are undocumented tuning parameters a blocker?

Not for the recommended wire-accuracy project. The manual exposes most of the choices that might otherwise be mistaken for secret tuning:

- Sec. 3.4 fixes the printed junction expansion at `n = 2`, uses the junction-segment length for the charge basis and match points, and explains which junctions require the local solve.
- Appendix B provides the complete branch logic for the cancellation-resistant field formulas and prints the five switch constants used by the single-precision program: `C1=10`, `C2=C3=0.01(2*pi)`, `C4=0.01`, and `C5=0.04`. It explicitly suggests increasing `C1` and decreasing `C2` through `C5` by two orders of magnitude as a first guess for double precision.
- Appendices A, C, and E specify the junction-potential, segment-integral, and end-cap evaluations, including series or quadrature alternatives and their truncation behavior.

There are implementation choices, but they can be selected against high-precision arithmetic and convergence rather than fitted to proprietary NEC4 output. For example, the Appendix B branch boundaries can be generated automatically by comparing the fast double-precision formulas with MPFR or long-double reference evaluations over a parameter grid, then tested for continuity across every boundary.

The larger buried/cross-interface ground project has more empirical machinery. Section 5 prints the interpolation regions, fitting functions, sample locations and counts, and asymptotic strategy, but also says those functions were selected by studying interpolation errors, reports roughly four-percent error in one illustrated fit, and acknowledges cases with larger errors. Reproducing a particular NEC4 build's exact ground results would therefore be harder than implementing a controlled method from the documented equations. This is another reason to keep that work separate from the initial wire engine.

The principal uncertainty in the wire engine is version drift, not missing mathematics: later NEC 4.2 behavior may incorporate choices or corrections made after the 1992 manual. Community-published NEC4/5 results can reveal such differences, but the implementation should remain traceable to the manual unless an independently justified correction is adopted and documented.

The 2011 GN3 report is different: it is a requirements-and-results description, not an equation-level specification. It says GN3 uses adaptive-density cubic interpolation over more regions, sometimes interpolates the two branch-cut integrals separately, tests total versus branch-cut interpolation accuracy to select a representation, evaluates both E and H directly, and uses improved higher-order asymptotics. It does not print the adaptive error targets, region construction, switching thresholds, branch-cut formulas, direct-H formulas, or higher-order asymptotic formulas. An exact independent reproduction of NEC 4.2 GN3 is therefore not supported by this report alone.

## What changed between NEC2 and documented NEC4

The differences fall into three layers. NEC4 retains the NEC3 capabilities, so a comparison directly against NEC2 must distinguish NEC4-specific numerical changes from NEC3 functionality inherited by NEC4.

### Layer 1 - NEC4 wire-accuracy formulation

| Change | Manual basis | Practical effect | Initial risk |
|---|---|---|---|
| Cancellation-resistant current expansion | Sec. 3.2, Eq. 3-3 replaces `A + B sin(ks) + C cos(ks)` with `A + B sin(ks) + C[cos(ks)-1]` | Avoids loss of precision on electrically short segments and very low-frequency models | Medium; algebra touches basis construction, current reconstruction, fields, loads, sources, and output |
| Per-end charge factors | Secs. 3.2 and 3.4 | Allows the two ends of a segment to participate in different junction charge conditions | Medium-high; the classic code assumes one radius-derived factor in several intertwined routines |
| Local junction charge solve | Sec. 3.4, Eqs. 3-42 through 3-44; Appendix A | Replaces NEC2's Wu-King logarithmic rule at unequal-radius and nonsymmetric junctions with a small scalar-potential MoM solve | Medium-high; localized mathematics, but topology and basis coefficients must be exact |
| NEC4 thin-wire kernel | Sec. 2.1 and Sec. 4.1 | Places tubular current on the wire surface while matching the boundary condition on the axis; drops artificial end-point-charge fields | High; changes every wire-wire matrix element and several near-field paths |
| Precision-safe segment fields | Sec. 4.1 and Appendix B | Uses limiting expansions in small radial distance, small electrical length, and large distance to prevent cancellation | Medium-high; well specified, but branch boundaries and continuity require careful tests |
| Closed flat wire ends | Sec. 4.4, Eqs. 4-35 through 4-40; Appendix E | Adds a constant-charge flat end cap to free ends and cap fields at voltage-source gaps; analogous corrections apply to loads | High; affects basis end conditions, matrix fill, and excitation/load terms |
| No NEC2 extended thin-wire kernel mode in the new formulation | Secs. 2.1 and 4.1 | NEC4 uses the extended-boundary-condition thin-wire model plus end caps instead of NEC2's optional EK approximation | Medium as a compatibility/UI issue; old EK behavior should remain available in the classic engine |
| Correct network-port handling in maximum coupling | Sec. 7.5 | Correctly shorts network/transmission-line ports while extracting the two-port admittance | Low-medium and separable from the wire engine |

These changes are coupled. In particular, a NEC4-like kernel without the end-cap model does not represent the closed surface required by its extended boundary condition. Likewise, the local junction solve produces per-end charge factors that must be consumed by the revised basis algebra.

#### Kernel terminology, accuracy and validity

LLNL-PROC-409033 removes an important ambiguity in the phrase "extended boundary condition." The documented NEC4 wire operator is still a **reduced thin-wire kernel (RTWK)**: current is represented on the cylindrical surface and the point-matched field is evaluated on the axis. It is not the later **extended thin-wire kernel (ETWK)** and is not an exact cylindrical-surface kernel.

The 2008 paper's ETWK should also not be assumed identical to the legacy NEC2 `EK` card. The paper defines ETWK by exact cylindrical integration of the singular `1/R` term plus a reduced approximation for the remainder; the classic NEC2 extended kernel is described in the 1992 manual as a low-order series approximation to the surface kernel. They need separate names and validation in this fork.

The report gives three directly implementable reference operators:

- exact kernel: integrate `exp(-jkR)/R` over axial position and circumference;
- ETWK: integrate the singular `1/R` part accurately over the cylinder and use the reduced-filament approximation for the phase-bearing remainder;
- RTWK: use the reduced single axial integral with the circumferential location collapsed into `R0`.

It also derives a stable ETWK singularity extraction. The axial `1/R` integral is written in logarithmic form with sign selected to reduce cancellation, and the circumferential singular term integrates analytically to `4*pi*log(max(a,rho))`; only the smooth remainder needs low-order Gauss-Legendre quadrature.

The accuracy consequences should become explicit model limits and tests:

- RTWK error depends strongly on segment-length/radius ratio and becomes poor when `Delta/a` is of order one or less.
- ETWK error is nearly independent of `Delta/a`, is limited mainly by phase approximation as `a/lambda` grows, and remains stable at very small `Delta/a` in the reported dipole tests.
- A half-wave dipole with `a=0.01 lambda` becomes unstable in the NEC5 RTWK experiment by 80-320 segments (`Delta/a` from `0.625` to `0.156`), while ETWK remains stable through 640 segments (`Delta/a=0.078`). NEC4's sinusoidal basis converges faster than NEC5's triangular basis, so those exact instability points are not direct NEC4 thresholds.
- Close parallel wires expose a different thin-wire limitation: neither RTWK nor ETWK represents circumferential current concentration. A patch model agrees more closely with the ideal transmission line as spacing approaches diameter.
- Tight bends can place short segment cylinders inside one another; the report observes similar RTWK/ETWK stability differences for a square loop.

For strict documented NEC4 similarity, RTWK remains the required kernel. ETWK or the exact kernel is a valuable independent oracle and a possible later modern engine, but substituting it silently would produce a different method rather than a NEC4-similarity implementation.

#### What the NEC5 validation manual adds

The NEC5 validation manual is not another NEC4 implementation specification. Its value is that NEC5 reaches many of the same physical results through a substantially independent formulation: mixed potentials, piecewise-linear wire bases, triangular/rooftop surface bases, and different boundary-condition enforcement. Agreement between NEC4 and NEC5 is therefore stronger evidence than agreement between two closely related NEC wire codes.

It adds practical acceptance limits that should be part of this feature:

- It states a simple RTWK modeling rule: accuracy degrades for segment lengths below about five wire radii, and the current solution can blow up when segment length is much shorter than radius. NEC4's approximate end caps extend the usable range but do not remove the limit.
- Radius-step convergence is a discriminating test. The NEC2 current-on-axis/surface-match convention converges to the wrong solution, while the NEC4 surface-current/axis-match convention converges toward the correct radius-step behavior.
- Source discretization is part of the numerical method, not merely input syntax. The NEC4 voltage is distributed approximately over its source segment; making that segment shorter than its neighbors increases input-admittance and power-balance error. Acceptance models should keep the source segment comparable to adjacent segments and should separately test deliberately nonuniform cases.
- The apparent convergence of a narrowing source gap can include changing physical gap capacitance. A distributed source occupying a fixed physical width is a useful independent convergence experiment.
- Close parallel wires remain outside the thin-wire approximation when circumferential current concentration matters. The manual shows deterioration as spacing approaches wire diameter; a NEC4-like kernel cannot fix that by segmentation alone.
- Electrically tiny passive loops expose a distinct NEC4 point-matching failure: sampled scalar-potential gradients do not telescope exactly around a closed loop. The manual reports that multi-point segment integration postpones failure, while a mixed-potential formulation remains stable much farther. This is not solved by merely porting the documented NEC4 wire kernel.
- NEC4 wire-to-surface junction impedance is explicitly described as unreliable and sensitive to wire-segment length relative to patch size. Pattern results can be reasonable while connected-antenna impedance is not. This should be an exclusion or a prominently diagnosed experimental mode, not a release claim for the initial wire engine.
- For ordinary above-ground wires, NEC4 and NEC5 agree well over demanding height ranges. The manual also supplies horizontal-wire, long-wire surface-wave, small-loop, and radial-ground-screen families that can seed the real-ground acceptance suite. Its NEC5 potential-table ground model is described as similar in accuracy to NEC4.2 GN3, but it does not fill in GN3's missing implementation equations.

The resulting objective is more precise: reproduce documented NEC4 behavior in its valid domain, improve the cases NEC4 was designed to fix, and diagnose known-invalid geometries. “NEC4 similarity” should not be presented as inheriting NEC5's open-surface, junction, close-wire, or electrically tiny-loop robustness.

#### What the 2017 practitioner comparison adds

Steve Stearns' *Antenna Modeling for Radio Amateurs* presentation
(`reports/nec-tips.pdf`) is practitioner guidance rather than an implementation
specification or a numerical oracle. It nevertheless sharpens the validation
method in several useful ways:

- Slide 55's NEC2 rules (`Delta < lambda/20`, `Delta > 2d` preferred, roughly
  equal lengths at junctions, aligned segmentation on close parallel wires,
  and clearance above ground) are useful stress-grid coordinates. They are not
  automatically NEC4 hard limits; the NEC4 feature exists partly to extend
  several of those boundaries.
- Slide 63's `lambda/200` minimum height for horizontal wires with the NEC2
  Sommerfeld-Norton model is a useful above-ground crossover witness, not a
  restriction on the documented NEC4 buried/crossing formulation.
- Slide 73 recommends geometry checks, average-gain checks, segmentation and
  source perturbations, PEC-versus-real-ground loss comparison, and an
  independent method. These become acceptance techniques rather than a single
  universal percentage tolerance.
- Slides 87--93 show that multiple solvers can have similar conductance while
  their susceptance differs by an almost linear offset attributed to different
  delta-gap shunt capacitances. Cross-solver validation must align the feed
  model and compare admittance components before treating an impedance
  difference as an antenna-operator error.
- Slides 64--66 caution that soil conductivity and relative permittivity vary
  with frequency and are not independent. A broad sweep with constant material
  parameters is a declared numerical experiment, not a universal physical soil
  prediction.
- Slide 43 identifies the public multiradius bridge-current method as a useful
  future comparator for very close wires, acute junctions, and wire grids. It
  is not silently substituted for the documented NEC4 formulation.

The presentation supplies no raw result tables or generally applicable error
bars. It therefore supports perturbation tests and model-validity diagnostics,
but cannot justify an arbitrary global `0.5%`, `1%`, or `5%` safe-error cutoff.

#### Expanded feature framing: NEC4/5 similarity

NEC5 techniques are useful enough to justify an umbrella **NEC4/5 similarity** project, provided it contains distinguishable numerical modes rather than a single ambiguous claim:

1. A **documented NEC4 mode** implements the 1992 sinusoidal-basis wire formulation, local junction solve, reduced thin-wire kernel, numerical limiting forms, and flat end caps. This is the traceable compatibility and convergence target.
2. A **modern mixed-potential mode** independently applies public mixed-potential EFIE, piecewise-linear wire bases, extended or exact thin-wire kernels, and telescoping potential differences. This is the path to improved electrically-small-loop behavior and a cleaner foundation for later surface coupling.
3. A later **wire/surface mode** can use public triangular/rooftop surface bases and singular wire-surface junction bases. This is a substantially larger solver project and should not be implied by the initial wire-engine release.

The NEC5 validation manual describes these methods and cites the foundational public literature, but it is a validation document rather than a complete equation-by-equation NEC5 implementation specification. Therefore the second and third modes should initially be called **NEC5-inspired** or **mixed-potential**, with their own published mathematical provenance and validation. “NEC5-equivalent” should be reserved for a later finding that public specifications and independent tests support that claim.

This split also avoids forcing one compromise basis on every use case. The NEC4 sinusoidal basis converges efficiently for ordinary wire antennas, while the mixed-potential formulation offers more robust discrete conservation properties and a natural route to surfaces. Both can share geometry, excitations, solvers, parallel assembly, output, and much of the validation infrastructure while retaining separate operator implementations.

### Layer 2 - NEC3 capabilities inherited by NEC4

| Change from NEC2 | Manual basis | Practical effect | Initial risk |
|---|---|---|---|
| Wires below ground and wires crossing an air/ground interface | Introduction; Secs. 3.4, 5, 7.4.2; Appendices F-J | Extends the Sommerfeld solution from above-ground sources/observers to either medium and supplies interface charge conditions | Very high; geometry, basis wavenumbers, source/observer transforms, interpolation tables, asymptotics, fields, power, and patterns all change |
| Wires and MFIE surfaces in an infinite dielectric or lossy medium | Introduction and Sec. 4 | Permits a homogeneous non-free-space background | High; complex/material wavenumbers propagate throughout the classic wire/surface engine |
| NEC3/4 Sommerfeld table and asymptotic machinery | Sec. 5 and Appendices F-J | Evaluates cross-interface and near-interface fields over a much larger parameter space | Very high; effectively a replacement/extension of the classic ground subsystem |

#### NEC 4.2 GN3 update

LLNL-TR-490316 identifies the original documented Sommerfeld/interpolation method as `GN2` and the later NEC 4.2 method as `GN3`.

| GN3 change | Motivation and effect | Specification completeness |
|---|---|---|
| Cubic interpolation with additional regions and adaptively selected sample density | Improves difficult ground-parameter ranges, especially seawater at lower HF | Algorithmic intent only; tolerances and region-generation rules are not given |
| Separate interpolation of two branch-cut integrals where needed | Preserves two-ray interference and avoids division through field zeros | Selection principle is described; equations and decision thresholds are not given |
| Accuracy test choosing total-field or branch-cut interpolation | Adapts the representation to interference regions | Error metric and acceptance threshold are not given |
| Direct evaluation of both E and H | Avoids magnification of E-field interpolation error when GN2 obtains H from a numerical curl | Direct-H formulation is not given |
| Higher-order asymptotics along the interface | Preserves interference near grazing incidence that GN2's first-order below-ground ray loses | Higher-order formulas are not given |
| New small-rho/large-z asymptotic treatment | Replaces a GN2 interpolation workaround outside the usual large-argument Bessel regime | Formulation and switch logic are not given |
| Table caching | Amortizes table construction that can take tens of seconds | Straightforward behavior; file format is not described |

The report supplies useful external witnesses: a 5 m, 10-segment whip on a 0.2 m, 6-segment seawater ground stake at 5 MHz and 5 mm wire radius gives `58.76 - j765.42` ohms with GN2 and `58.88 - j765.56` ohms with GN3 for `epsilon_r=81`, `sigma=4 S/m`. It also documents qualitative GN2 E/H glitches and loss of interface interference when GN2 switches to asymptotics beyond two wavelengths. The prose does not fully specify the excitation and deck, so the impedance is corroboration rather than a reproducible acceptance oracle until the original model is found.

This is not a prerequisite for the free-space/above-ground NEC4 wire-accuracy layer. It should be treated as a separate objective rather than allowed to make the initial feature unbounded.

### Layer 3 - capabilities mentioned but not fully specified here

| Capability | Assessment |
|---|---|
| Insulated wire | The preface and introduction say it is part of NEC4, and Secs. 3.2/4.1 allow a distinct current-expansion wavenumber, but this volume does not present a complete insulation operator. The fork's existing `LD 7` coated-wire model is independently derived and validated; it should not be relabeled NEC4 insulation. Additional primary literature or the relevant code/user documentation is needed. |
| Generalized combinations of incident waves and voltage sources | Public NEC4 development literature identifies this as a new feature, but the theory volume does not define its input semantics. It is mostly excitation orchestration after the underlying linear system is available. |
| Monostatic RCS angle sweep | Primarily a command/output workflow around existing scattering calculations; it needs the User's Guide for syntax and exact behavior. |
| Exact NEC 4.2 behavior | Not guaranteed by a 1992 NEC4 theory volume. It needs a licensed executable as an oracle, published update notes, or explicitly accepted independent behavior. |

## Mapping to the current fork

The current repository still contains the recognizable NEC2 wire core, now surrounded by performance, dielectric-volume, and coated-wire extensions.

| NEC4 concern | Current implementation area | Observed gap or overlap |
|---|---|---|
| Junction topology | `geometry.c:conect()` and `data.icon1/icon2` | Existing topology is adequate for enumerating junction members, although it will need a stable junction object or precomputed per-end factors |
| NEC2 junction/basis rule | `calculations.c:sbf()`, `tbf()`, `trio()` | Uses the radius/log-derived `aj`, `ap` factors characteristic of NEC2; this is the central replacement point for per-end NEC4 factors |
| Current reconstruction | `calculations.c:cabc()` and current/output consumers | Stores separate constant/sine/cosine coefficients and reconstructs the center current as `A+C`; a cancellation-resistant representation must be carried consistently |
| Wire field kernel | `fields.c:efld()`, `eksc()`, `ekscx()`, `gx()`, `gxx()` | Implements the NEC2 thin-wire and optional EK paths; it does not implement the complete NEC4 surface-current/axis-match kernel, precision branches, or cap fields |
| Matrix fill and fields | `matrix.c`, `fields.c`, `calculations.c`, parallel fill paths | Every classic matrix-fill implementation must share one reviewed kernel to prevent serial/parallel numerical divergence |
| Free ends, sources, and loads | `calculations.c:tbf()`, excitation assembly in `main.c`, load handling | Some NEC2 free-end coefficient logic exists, but not the documented NEC4 cap charge fields at free ends/source gaps/load gaps |
| Real ground | `ground.c`, `somnec.c`, `fields.c` | Implements the NEC2 above-ground model; classic geometry explicitly rejects below-ground segments |
| Volume dielectric | `dielectric_*`, `vie_*` | Valuable validation/coupling infrastructure but not a substitute for NEC3/4 buried thin wires |
| Coated wire | `coated_wire.c` and `LD 7` integration | Independent thin-coating operator with validation limits. It overlaps the user problem addressed by NEC4 insulation, but not the documented algorithm |
| Dense and iterative solvers | `dense_backend.c`, `dielectric_fft.c`, `vie_iterative.c` | Existing LAPACK/parallel infrastructure is reusable; the principal work is operator correctness rather than linear algebra |

The safest architecture is a new wire-physics mode behind a single operator interface, while preserving the current NEC2 engine unchanged as a compatibility and regression oracle. Trying to mutate `sbf/tbf/efld` piecemeal without an engine boundary would make it difficult to know whether a discrepancy comes from basis algebra, kernel evaluation, end caps, ground, or the existing performance paths.

## Feasible project objectives

### Objective A - documented NEC4 wire-accuracy subset (recommended first project)

Implement the 1992 documented formulation for free-space, PEC ground, and existing above-ground real-ground geometries:

1. precision-safe basis representation and reconstruction;
2. per-end junction factors and the local junction charge solve;
3. NEC4 thin-wire kernel and Appendix B numerical branches;
4. flat end caps for free ends, voltage sources, and loads;
5. corrected maximum-coupling network-port behavior;
6. unchanged NEC2 syntax and an explicit CLI selector retaining the classic engine.

This objective targets the largest practical NEC2 accuracy limitations without first rebuilding the entire half-space engine. It is substantial but bounded and can be validated incrementally.

### Objective B - NEC3/4 buried and cross-interface wires

After Objective A is stable, extend the classic half-space formulation to homogeneous-medium, buried-wire, and crossing-wire cases using Secs. 3.4/5/7.4.2 and Appendices F-J. This is a separate major project with its own geometry, operator, power, pattern, and convergence program.

### Objective C - broader NEC4-family workflow similarity

Add user-facing excitation/RCS conveniences and investigate insulated-wire behavior. This requires the NEC4 User's Guide and additional primary sources. Exact NEC 4.2 compatibility would additionally require a licensed validation oracle and a policy decision about behavioral targets.

## Recommended definition of success

For Objective A, success should mean **independent implementation of the documented 1992 NEC4 wire formulation**, not source compatibility with proprietary NEC4 and not an unsupported claim of exact NEC 4.2 equivalence.

The new engine should:

- preserve every existing NEC2 card and result schema;
- be selectable by command line, with the classic NEC2 path preserved bit-for-bit where practical;
- initially remain opt-in until its acceptance suite is complete;
- converge more reliably than classic NEC2 for electrically short segments, radius steps, multiwire junctions, small segment-length/radius ratios, and source-adjacent refinement;
- preserve reciprocity and power balance to explicit tolerances;
- produce identical results across serial and parallel matrix assembly;
- have no material regression for ordinary uniform-radius thin-wire models;
- publish capability boundaries instead of silently applying unsupported combinations.

After acceptance, making the new engine the default can be a separate user decision, with a CLI flag retaining classic behavior.

## Validation strategy available to this repository

1. **Equation-level unit oracles:** implement the manual's Appendix A, B, and E equations independently from production code and test junction matrix elements, limiting field forms, and cap fields.
2. **Null limits:** equal-radius two-member junctions, long-enough segments, vanishing radius, and ordinary uniform dipoles should approach the established solution appropriately.
3. **Convergence studies:** radius-step dipoles, tapered wires, T/Y/star junctions, close-wire junctions, electrically small loops/structures, thick-wire monopoles, source refinement, and loaded segments.
4. **Structural invariants:** reciprocity, Kirchhoff current continuity, scalar-potential continuity at solved junctions, passivity, power balance, symmetry, rotation/translation invariance, and deterministic serial/parallel identity.
5. **Independent implementations:** the MIT-licensed `momwire` project can serve as a secondary implementation and test-generator, but claims or later-version adjustments must not replace the printed equations without an explicit behavioral decision.
6. **Published NEC4/5 witnesses:** reconstruct sufficiently specified models from the NEC5 validation manual and other public reports. Treat digitized plots as tolerance bands rather than invented exact values, and preserve citations and reconstruction assumptions with each deck.
7. **Physical and full-wave corroboration:** analytical canonical cases and openEMS remain useful, particularly where two MoM implementations may share assumptions.

## Decisions needed before milestone planning

1. Is the target **documented NEC4 (1992) wire physics**, or should later NEC 4.2 black-box behavior be treated as normative when it differs from the manual?
2. Is Objective A the first bounded project, leaving buried/cross-interface wires for a later project?
3. Should the new engine begin as `--wire-engine nec4-similar` (name provisional), with `--wire-engine nec2` preserving the original path?
4. Which published NEC4/5 reports should be elevated from corroborating plots to reconstructed acceptance decks? The project will not acquire or run the licensed NEC4 binary; public equations, published numerical witnesses, analytic cases, `momwire`, and openEMS form the oracle set.
5. Should the first release explicitly exclude coated wires from the new engine until the independent `LD 7` operator is revalidated on top of the revised basis/kernel?

## Provenance and licensing guardrail

The implementation should be written from the public mathematical description, with equation and section references in design notes and tests. The proprietary NEC 4.2 source must not be used. LLNL's current licensing page describes NEC 4.2 as separately licensed software and lists buried and insulated wires among its capabilities: <https://softwarelicensing.llnl.gov/product/nec-v42>.

This is an engineering provenance recommendation, not legal advice. Before advertising compatibility under the NEC4 name, the project should separately review license and trademark implications and use conservative capability-specific wording.
