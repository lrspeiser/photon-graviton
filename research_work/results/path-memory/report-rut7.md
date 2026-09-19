# RUT-1 stage 7: the cold annulus is collectively unstable within ten periods, the warm ones are quiet for thirty, and the cold ring's spectrum is certified to the axis

Protocol `protocol-rut7.md` (e59ef6f) and `protocol-rut7-amendment-1.md` (cd8bbd8), both declared and pushed
before the campaign. Campaign driver `rut7.py`, its units of work `rut7_tasks.py`, modules `spectrum.py`,
`population.py`, `pairrun.py`; archive `rut7-results.json` and `rut7-series/` (19 files, 9.6 MB); suite job
`rut7_checks.py`. The campaign was 411 tasks, 61,000 CPU-seconds, 47 minutes on 22 workers.

This stage is the execution of items 2 to 5 of the owner's handoff in `research_work/annulus_sampling`. Their
corrected node sampler is the sampler of record at every source-sampling call site. Stage 6 is untouched:
`equilibrium.py` and `rut6.py` are byte-identical, its archive still reproduces, and **H6 stays failed**.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | **passes.** `rut7_checks.py` re-hashes the code, the owner's package, both protocol files and all 19 series files; recomputes 32 items from scratch (E1, E2, E4, the E3 control, two spectra, all five populations, V1 and V3 on all five, V2, V4, V5, V7, the declared part R run); replays 2 T0 prefixes of one stationarity run and one live and frozen pair; and recomputes every gate and reading that is derived from the series. Zero differences, 48 seconds. |
| **numerical verification** | **passes: 13 of 13 gates, and every negative control is rejected.** No deviation from the protocol as amended. |
| **scientific outcome** | **a legitimate negative result for the cold annuli and a bounded positive one for the warm.** With the physics held fixed, the annulus at dE = 0.010 roughly doubles its radial spread and its residual dispersion in 30 periods at every body count; at dE = 0.020 and 0.030 every declared quantity is quiet at 1,024 bodies. **No population is claimed stable**: thirty periods cannot exclude a slow mode, and the linear mode calculation that could has not been run. |

Two of my declared part X readings turned out to be badly designed, and they are reported as they fell with
the defect named beside them, not repaired after the fact. They are in "Where two declared readings were
wrong-headed" below.

## Part E. The cold-ring spectrum, certified down to the axis

The strip is Re s from −0.008 to 1.5 and |Im s| ≤ 6, for m = 0 to 16, for four rungs at the stage 5
configuration (32 writers, w/R = 0.2, the 10% label). Its left edge lies to the left of the imaginary
axis, so marginal and slowly growing roots are inside the contour. Every sub-rectangle gives two
independent answers from one set of contour evaluations: Beyn's moments locate, and the argument principle
counts.

| check over all 1,768 sub-rectangles (4 rungs × 17 m × 2 partitions × 13) | worst case | declared |
|---|---|---|
| Beyn's count with multiplicity equals the winding number | every one | every one |
| moment capacity saturated | never | never |
| smallest gap in the singular spectrum at the chosen rank | 4.8×10⁵ | ≥ 10³ |
| polished residual | 8.4×10⁻¹⁴ | < 10⁻¹⁰ |
| distance the polish moved a root | 6.7×10⁻⁸ | < 10⁻⁴ |
| largest phase step in the winding integral | 0.386 | ≤ π/8 = 0.393 |
| winding number's distance from an integer | 3.8×10⁻¹⁵ | — |
| the two partitions, which share no cut, return the same root set | to 8.4×10⁻¹⁵ | to 10⁻⁸ |

The declared edge rule fired once, in both partitions of the one-stage ring at m = 15: a *stable* root at
−0.00702 − 0.00695i lay 0.00098 from the left edge, the edge moved to −0.010 and the strip was swept again.
Both sweeps are archived.

**The gates and their negative controls.**

| gate | result | the control it had to reject |
|---|---|---|
| E1 known spectra | exp(−s) = ½: three roots, **six with multiplicity**, winding number six, unpolished error 1.7×10⁻¹⁴. Zero-lag quartic: 46 located of 46 expected over m = 0–16, unpolished error 2.8×10⁻⁹. Kepler control: all four roots located at every m tried and none called unstable | the stage 6 rank rule on the first problem returns **five** of six: rejected |
| E2 raw robustness | the *unpolished* located roots move 3.2×10⁻¹⁴ at most under doubled quadrature and under a shifted contour (one-stage m = 2, 5; two-stage m = 1, 2, 3) | stage 6's quadrature moves the two-stage m = 2 root by **1.2×10⁻⁵** against its own doubling: rejected |
| E3 completeness | the table above | a synthetic T(s) holding a root whose residue is 10⁻⁹ of its neighbour's: Beyn's count 1, winding number 2 — **flagged**, not dropped |
| E4 edge handling | a root 10⁻⁶ inside each of the four edges of a sub-rectangle is found at all four, with the winding integral still resolved to π/8 | the stage 6 solver, which filters and then polishes, **loses it at the left edge**: rejected |

One thing here is worth knowing beyond the gate. With the root only 10⁻⁶ from the contour the located value
is still accurate to 10⁻¹⁵ on three edges. Beyn's method is a rational filter: a quadrature error multiplies
a pole's residue by the same factor in every moment, so the Hankel structure survives and the eigenvalue
comes out exact. What stage 6 lost at an edge, it lost to its order of operations, not to the contour.

### What the strip holds

| rung | roots with multiplicity | unstable | marginal (\|Re s\| ≤ 10⁻⁷) | fastest |
|---|---|---|---|---|
| no attraction (Kepler control) | 68 | **0** | 68 | — |
| zero-delay field | 46 | 22 | 24 | m = 4, 0.155 T0 |
| one-stage | 55 | 34 | 2 | m = 5, 4.45 T0 |
| two-stage | 71 | **24** | 22 | m = 2, 8.58 T0 |

**The two-stage ring is unstable at every m from 1 to 15, on four different branches, and only one of them
is fast.** Stage 6's "exactly one unstable mode", withdrawn in 70aa1c0, is replaced by this list:

| branch | where | growth | what it is |
|---|---|---|---|
| the known mode | m = 2, s = 0.018552 + 2.096535i | **8.58 T0** | the inertially stationary two-lobed pattern of stages 4 to 6 (pattern speed 0.0005 Ω) |
| memory pole | m = 1, s = 5.58×10⁻⁴ + 1.048779i | 285 T0 | a **static lopsided** mode, Im s = Ω to 2×10⁻⁶: a fixed dent in the written field. Continuation cannot reach this branch |
| co-rotating | every m from 1 to 15, Im s from 0.0065 down to 0.0008 | 732 T0 at m = 1, 1,533 at m = 2, 2,432 at m = 3, rising to 1.6×10⁵ at m = 15 | a pattern that drifts slowly against the bodies (pattern speed 0.994 to 0.9999 Ω), each with a stable mirror root. This is the root stage 5's first solver found and reported as "neutral, 732 periods" |
| epicyclic | m = 3 to 9, Im s ≈ 2.4686 | 852 T0 at m = 3, 1.4×10⁴ at m = 4, to 1.2×10⁶ at m = 9 | the radial oscillation, very weakly driven |

At m = 0 and at m = 16 — the Nyquist mode of 32 writers — nothing is unstable. The one-stage ring has its
own static lopsided m = 1 root at 7.3×10⁻⁴ (217 T0), which no earlier stage had seen.

**Twenty-two roots of the two-stage ring, counted with multiplicity, sit inside the resolution floor** and are
called marginal, not stable: the m = 0 double root at the origin (the ring's neutral displacement, located
at 1.7×10⁻⁸ because a double root is known only to the square root of the precision), and twenty roots at
m = 6 to 16 — the epicyclic roots, and at m = 16 the slow pair too — whose real parts are 10⁻⁸ to 10⁻¹⁰ of
either sign. At this precision their sign is not known, and it is not claimed.

**Continuation, kept as a cross-check only.** Stage 5's continuation returns a duplicated root at every m
from 1 to 15 for the two-stage ring. Deflation recovers the lost partner in every case, and it is the
unstable co-rotating root each time. Even so continuation still lacks five certified two-stage roots — the
static m = 1 mode among them — and two one-stage ones. It is not a gate and nothing depends on it.

## Part S. The population, defined, and drawn faithfully

    f(E, L) = A exp[−(L−L0)²/2dL²] exp[−(E−E_circ(L))²/2dE²]   for |L−L0| ≤ 4dL and 0 ≤ E−E_circ(L) ≤ 4dE,

zero otherwise, with E_circ(L) the global minimum over radius of Φ + L²/2r², A fixed by a total mass of 1,
L0 = 1, dL = 0.06, w = 0.2, τ_keep = 10 T0, τ_form = 3 T0. The truncation is part of the definition.

| population | family | dE | α | α × mass | mass-weighted support | its spread over the population | mean radius | radial width | σ_r |
|---|---|---|---|---|---|---|---|---|---|
| A_cold | A, matched support | 0.010 | 0.03434181 | 0.034342 | 10.000% | ±16.6% | 0.92314 | 0.04780 | 0.0940 |
| A_warm | A, matched support | 0.030 | 0.02943086 | 0.029431 | 10.000% | ±20.3% | 0.95571 | 0.09097 | 0.1659 |
| B_cold | B, fixed physics | 0.010 | 0.03156101 | 0.031561 | 9.289% | ±15.9% | 0.92940 | 0.05041 | 0.0940 |
| B_mid | B, calibrates α | 0.020 | 0.03156101 | 0.031561 | 10.000% | ±19.4% | 0.93838 | 0.06850 | 0.1348 |
| B_warm | B, fixed physics | 0.030 | 0.03156101 | 0.031561 | 10.607% | ±21.4% | 0.94931 | 0.08647 | 0.1661 |

In family A the coupling is re-solved for each temperature, so the two members differ in total writing by
17%: that family asks only whether supported populations exist across temperatures. **In family B the mass,
the coupling, the footprint and the response times are fixed, and the support that results is reported, not
forced.** A "10% support" is a population mean, with a mass-weighted spread of ±16 to ±21% of Newtonian: across
A_cold the memory force runs from 68% of Newtonian inward on one flank to 27% outward on the other.

| gate | result | the control it had to reject |
|---|---|---|
| V1 the measure | the owner's `verify_distribution` on all five: radial marginal to 6.5×10⁻¹⁶, mass and first two radius moments to 10⁻¹⁵, source unmutated | the stage 6 weights are rejected on all five, with a mean-radius shift equal to Var(R)/Mean(R) to 10⁻¹⁵: 0.0025 (A_cold), 0.0087 (A_warm), 0.0027 (B_cold), 0.0050 (B_mid), 0.0079 (B_warm) |
| V2 drawn moments | 200,000 bodies from **each of two samplers** — the owner's node sampler and a continuous rejection sampler that places no body on a node: worst \|z\| of any mean 2.3 against quadrature and 2.3 between the samplers; every variance within 1% | a draw re-weighted by R is off by 7.1σ (coldest) to 19.8σ in mean radius: rejected on all five |
| V3 the support fits | every populated orbit is bound, inside the solver's domain and inside 0.9 of the box; pericentre ≥ 0.455, apocentre ≤ 2.797; E_circ's radius increasing and continuous in L | the stage 6 five-width support reaches the edge of the domain for both warm annuli and for B_mid: rejected |
| V4 discretization | five refinements of each population, separately and together: α, mean support, mean radius, width and peak field move ≤ 3.2×10⁻⁵ (declared 10⁻⁴); profiles ≤ 1.0×10⁻⁵ of their peak (declared 2×10⁻³) | stage 6's label — support at the argmax node, unbracketed secant — moves its α by **23%** and **19%** under the same radial refinement: rejected |
| V5 the drawn source writes the field | rms error 0.63% → 0.052% of the peak field from 256 to 16,384 bodies (coldest), fitted exponents −0.44 to −0.58 | the R-re-weighted sample through the same gate: exponents −0.03 to −0.10, outside [−0.7, −0.3]: rejected on all five |
| V6 the draw is stationary | 16,384 bodies in the frozen field for 30 T0, four realizations of each population: worst single \|z\| 2.8, worst ensemble \|z\| 2.1, mean-radius fluctuation ≤ 1.2 standard errors | kicked outward by 0.15: ensemble \|z\| 90 to 183. Velocities scaled by 1.05: 104 to 141. R-re-weighted: 7.1 to 18.1. All rejected, 60 runs |
| V7 restart | bodies, E and C saved at 2 T0 and resumed: positions, velocities, both fields, every record and the step count identical **bit for bit** at 4 T0 | a restart from the bodies alone differs by 1.1×10⁻²: rejected |
| V8 the ensemble | 284 runs, every declared one completed; all 100 pairs start from byte-identical bodies | — |

**What the amendment changed, and whether it mattered.** Both original control wordings are archived beside
the amended ones. As first worded, V5's control needed an error floor above 0.5% of the peak field: the cold
annuli's floors are 0.31% and 0.33%. V6's needed the re-weighted sample to fail a per-realization test: for
the cold annuli it fails in one realization of four and two of four. So as worded both gates would have been
recorded as failed for every cold population, for a reason that had nothing to do with the sampler: the
stage 6 defect shifts a cold annulus by only Var(R)/Mean(R) = 0.0025. Put through the gate it imitates, the
biased sample is caught every time — its error stops falling with N (V5), and its relaxation is small but
always in the same direction, 7σ on the ensemble (V6).

## Part R. The refinement stage 5 quoted, now archived

The unseeded two-stage cold ring, 30 T0, the m = 2 mode emerging alone from the numerical floor, fitted
over the declared 10–30 T0. Prediction: 0.018552 and 2.096535.

| grid spacing, step | growth rate | against prediction | frequency | against prediction |
|---|---|---|---|---|
| w/5, 0.01 | 0.018839 | +1.55% | 2.0965772 | 2.0×10⁻⁵ |
| w/10, 0.01 | 0.018857 | +1.64% | 2.0965764 | 2.0×10⁻⁵ |
| w/5, 0.005 | 0.018871 | +1.72% | 2.0965744 | 1.9×10⁻⁵ |
| control: no primed field | 0.045117 | +143% | 1.7681588 | 16% — **rejected** |

R1 passes: all three within the declared 2% and 0.5%. Refining the grid or the step moves the rate by 0.1%
to 0.2%, so the simulator is converged at the declared resolution and the common +1.6% is not
discretization. It is the fit window: over 20–30 T0 (exploratory, from the archived series) the three give
0.018541, 0.018564 and 0.018572 — **−0.06%, +0.06%, +0.11%**. Stage 5's unarchived table showed ±0.8%
"scatter"; archived and declared, the scatter is under 0.2%.

## Part X. Identical bodies under a live field and a frozen one

Five populations × (8, 8, 4) realizations at 64, 256 and 1,024 bodies × two runs: 200 runs of 30 T0, all
completed. Each pair starts from byte-identical bodies. The frozen run holds the constructed field fixed and
is the stationary reference; the live run lets the same bodies write. **D_Q is the difference of the two
runs' changes** between the first 2 T0 and the last 5 T0, over a scale that belongs to the initial
population, with its standard error over realizations. "Dispersion" is what is left after a mean flow, a
breathing term and the m = 1 and m = 2 streaming are fitted out.

### What happened, at 1,024 bodies (D in percent; mean ± standard error over four pairs)

| quantity | A_cold | B_cold | B_mid | A_warm | B_warm |
|---|---|---|---|---|---|
| mean radius | −0.99 ± 0.12 | −0.63 ± 0.04 | +0.06 ± 0.03 | +0.03 ± 0.12 | +0.03 ± 0.08 |
| radial spread | **+119.7 ± 1.8** | **+73.2 ± 3.4** | +0.87 ± 0.21 | −0.28 ± 0.46 | −0.33 ± 0.46 |
| residual radial dispersion | **+119.6 ± 2.0** | **+72.4 ± 4.4** | +0.19 ± 0.10 | +0.35 ± 0.66 | +0.60 ± 0.55 |
| residual azimuthal dispersion | +5.5 ± 1.5 | +7.6 ± 0.7 | +0.13 ± 0.10 | −0.08 ± 0.11 | −0.11 ± 0.17 |
| rms v_r (reported under that name) | +135.2 ± 2.4 | +88.9 ± 1.6 | +0.26 ± 0.10 | +0.40 ± 0.62 | +0.63 ± 0.57 |
| mean angular momentum | **−3.52 ± 0.06** | **−2.01 ± 0.05** | −0.021 ± 0.005 | −0.0018 ± 0.0003 | −0.0049 ± 0.0016 |
| spread of angular momentum | −15.9 ± 0.7 | −10.0 ± 1.0 | −0.07 ± 0.06 | −0.01 ± 0.01 | −0.04 ± 0.02 |
| support felt by the bodies | +11.0 ± 1.1 | +9.0 ± 0.9 | −0.84 ± 0.16 | −0.13 ± 1.38 | −0.52 ± 1.06 |
| m = 2 streaming (over initial σ_r) | +110.9 ± 13.1 | +100.5 ± 10.1 | +1.2 ± 1.3 | +0.1 ± 1.4 | +0.4 ± 0.4 |

The frozen run's own change is under 0.8% in every row above the streaming one, and under 2.1% of σ_r in
that one: the reference is stationary, so D is the live field's doing.

**The cold annuli are destroyed as cold annuli, at every body count.** Radial spread D for A_cold is
+116 ± 6, +113 ± 4 and +120 ± 2 percent at 64, 256 and 1,024 bodies; for B_cold +87 ± 5, +70 ± 4 and
+73 ± 3. This time it is heating and not a stream miscounted as heat, which was stage 6's mistake: the
*residual* dispersion, with the m = 2 stream fitted out, rises by the same factor.

**The warm annuli do nothing resolvable.** Every declared quantity is consistent with zero at 1,024 bodies.

### The declared readings, as declared

| population | mean radius (1%) | radial spread (10%) | residual radial dispersion (10%) | residual azimuthal dispersion (10%) | support (10%) | discreteness exponent |
|---|---|---|---|---|---|---|
| A_cold | unresolved | **exceeds the bound** | **exceeds the bound** | quiet | unresolved | collective, p = −0.04 ± 0.02 |
| B_cold | quiet | **exceeds the bound** | **exceeds the bound** | quiet | unresolved | collective, p = −0.11 ± 0.03 |
| B_mid | quiet | quiet | quiet | quiet | quiet | unresolved |
| A_warm | quiet | quiet | quiet | quiet | quiet | unresolved |
| B_warm | quiet | quiet | quiet | quiet | quiet | unresolved |

"Quiet" means |mean D| plus two standard errors fits inside the bound at 1,024 bodies; it is a bound on a
drift over thirty periods and not a stability result. B_mid's spread, +0.87 ± 0.21%, is resolved from zero
at four standard errors and is small.

**The declared collective-mode reading, m = 2, growth of |C₂| over 15–30 T0:** A_cold **absent**
(−0.0113 ± 0.0022 at 1,024 bodies), B_cold **absent** (−0.0156 ± 0.0063), B_warm absent, A_warm and B_mid
unresolved. Of twenty population-and-mode cells, none reads "present".

### Where two declared readings were wrong-headed

They are reported above exactly as they fell. Neither is repaired here; both need a new declaration.

**1. The mode window was set after the event.** I put the fit window at 15–30 T0 so that the shot-noise
structure of a freshly live field would have saturated. The cold annuli grow *and saturate* before 15 T0,
so the window sees the decay that follows, and the rule reads "absent" for a mode that doubled the
population's dispersion. "Absent" is the rule's output and is false as a description. I had the evidence to
know better: replaying stage 6's archived cold run had given an e-folding of 5.3 to 5.7 T0, and I sized the
window on the cold *ring's* 8.58 T0 instead.

**2. The discreteness exponent cannot separate the two things it was declared to separate.** I declared
that heating falling as N^p with p below −½ is discreteness and p above it is collective. But a collective
mode seeded by shot noise starts at an amplitude ∝ N^−½, so until it saturates the heating it causes falls
as 1/N too — exactly like two-body relaxation. The reading is sound only after saturation, which is where
the cold annuli are (p ≈ 0: collective, correctly). For B_mid the spread D is 10.9 ± 3.0, 5.6 ± 1.2 and
0.87 ± 0.21 percent at 64, 256 and 1,024 bodies — steeply falling — and that pattern is what *either*
hypothesis predicts. The declared rule says "unresolved" because the 1,024-body change is not two standard
errors from zero; the honest statement is that this reading could not have resolved it.

### Exploratory: the growth phase the declared window missed

**Exploratory, and labelled so in the archive.** The rule was chosen after the series were seen, it is not a
declared reading, and it is not a gate. For each population and body count: the geometric mean over
realizations of |C_m|/C₀; a fit from 6 T0, after the imprint of shot noise on the freshly live field, to the
first time that mean reaches half its peak; one fit per realization. If that window is shorter than 2 T0
there is no growth phase. The same three numbers are applied to every cell; none is tuned to a population.

| population | N | fit window (T0) | m = 2 growth rate | e-folding | pattern speed | \|C₂\|/C₀ at 6 T0 | peak | at |
|---|---|---|---|---|---|---|---|---|
| A_cold | 64 | 6–9.5 | 0.038 ± 0.008 | 4.2 T0 | +0.028 ± 0.004 | 1.2×10⁻³ | 5.3×10⁻³ | 13 T0 |
| A_cold | 256 | 6–11.0 | 0.045 ± 0.002 | 3.5 T0 | +0.031 ± 0.002 | 6.7×10⁻⁴ | 5.2×10⁻³ | 15 T0 |
| A_cold | 1,024 | 6–14.25 | 0.048 ± 0.007 | 3.3 T0 | +0.026 ± 0.001 | 2.0×10⁻⁴ | 4.8×10⁻³ | 18 T0 |
| B_cold | 64 | 6–9.0 | 0.053 ± 0.012 | 3.0 T0 | +0.029 ± 0.008 | 4.9×10⁻⁴ | 2.4×10⁻³ | 14 T0 |
| B_cold | 256 | 6–11.5 | 0.038 ± 0.013 | 4.1 T0 | +0.041 ± 0.007 | 3.9×10⁻⁴ | 2.7×10⁻³ | 16 T0 |
| B_cold | 1,024 | 6–14.5 | 0.037 ± 0.006 | 4.3 T0 | +0.045 ± 0.006 | 2.4×10⁻⁴ | 2.9×10⁻³ | 19 T0 |
| B_mid, A_warm, B_warm | all | none | — | — | — | 0.8 to 5.5×10⁻⁴ | — | the amplitude peaks during the imprint and does not grow afterwards |

This is what a mode of the underlying distribution looks like, and it is how the protocol said one would
be told from noise: **a growth rate that does not depend on the body count** (0.038, 0.045, 0.048), from a
seed that falls with N (exponent −0.65 for A_cold), **saturating at the same amplitude** (5×10⁻³ of the
axisymmetric field) but later the more bodies there are — 13, 15 and 18 T0. Its pattern speed is 0.03 to
0.04 against an orbital frequency of 1.17: **nearly stationary in the inertial frame**, like the cold
ring's mode and like the modes found by replaying stage 6. It is about twice as fast as the cold ring's
8.58 T0, and I do not know why; an annulus with a ±16% spread of support across it is not a ring.

The history at 1,024 bodies shows the sequence — the field pattern grows, the stream grows with it, and
the stream then decays into dispersion:

| t (T0) | A_cold: spread, residual σ_r, \|C₂\|/C₀ | B_cold | B_mid | A_warm |
|---|---|---|---|---|
| 0 | 0.0480, 0.0891, 0 | 0.0499, 0.0903, 0 | 0.0683, 0.1277, 0 | 0.0913, 0.1535, 0 |
| 6 | 0.0488, 0.0887, 2.0×10⁻⁴ | 0.0513, 0.0872, 2.4×10⁻⁴ | 0.0686, 0.1268, 2.4×10⁻⁴ | 0.0910, 0.1545, 1.4×10⁻⁴ |
| 10 | 0.0495, 0.0902, 7.4×10⁻⁴ | 0.0516, 0.0893, 3.0×10⁻⁴ | 0.0689, 0.1275, 2.0×10⁻⁴ | 0.0907, 0.1525, 6.6×10⁻⁵ |
| 14 | 0.0639, 0.1087, 2.4×10⁻³ | 0.0572, 0.0956, 1.3×10⁻³ | 0.0692, 0.1268, 1.3×10⁻⁴ | 0.0878, 0.1551, 4.4×10⁻⁵ |
| 18 | 0.0908, 0.1574, 4.8×10⁻³ | 0.0755, 0.1262, 2.7×10⁻³ | 0.0686, 0.1267, 1.3×10⁻⁴ | 0.0916, 0.1517, 4.2×10⁻⁵ |
| 25 | 0.1035, 0.1867, 2.1×10⁻³ | 0.0883, 0.1491, 1.8×10⁻³ | 0.0697, 0.1263, 1.0×10⁻⁴ | 0.0911, 0.1523, 2.7×10⁻⁵ |
| 30 | 0.1057, 0.1995, 2.2×10⁻³ | 0.0873, 0.1567, 8.8×10⁻⁴ | 0.0681, 0.1273, 1.8×10⁻⁴ | 0.0901, 0.1536, 2.2×10⁻⁵ |

**The cold annuli end where the warm one starts.** A_cold finishes at a spread of 0.106 and a residual
dispersion of 0.20, against A_warm's initial 0.091 and 0.154: the instability heats the population until
it is warmer than the warm annulus, and stops. In the warm annuli |C₂| *decays* after the imprint, by a
factor of seven in A_warm.

### Temperature, with the physics fixed

This is the comparison stage 6 claimed and had to withdraw: there the warm population carried 44% more
total writing. In family B the mass, α, the footprint and both response times are identical and only dE
changes.

| dE | σ_r | support that results | radial spread D at 1,024 | residual dispersion D | m = 2 growth phase |
|---|---|---|---|---|---|
| 0.010 | 0.094 | 9.29% | **+73.2 ± 3.4%** | **+72.4 ± 4.4%** | e-folding 3 to 4 T0 at every N |
| 0.020 | 0.135 | 10.00% | +0.87 ± 0.21% | +0.19 ± 0.10% | none in 30 T0 |
| 0.030 | 0.166 | 10.61% | −0.33 ± 0.46% | +0.60 ± 0.55% | none in 30 T0 |

**Velocity dispersion suppresses the collective mode, with everything else held fixed**, and the transition
lies between σ_r = 0.094 and 0.135. The colder population is the less supported one here (9.3% against
10.6%), so the stronger coupling is not what destabilises it. The declared summary of this was a straight
line through three points — the slope of the spread D against dE is −3.0 ± 0.5 per unit dE, six standard
errors — which is the right sign and the wrong shape: it is a threshold, and three points do not locate it.

**The bodies lose angular momentum, and it goes nowhere.** Mean L falls 3.5% and 2.0% in the cold annuli
while the mode is active, independent of N. In B_mid it falls 0.45%, 0.21% and 0.021% at 64, 256 and 1,024
bodies — shot noise exerting a torque. The field in this model carries no angular momentum of its own, so
this is a loss with no ledger. It is the reciprocal accounting the owner's item 6 asks for, and it has not
been done for any of these runs.

## What I got wrong on the way

* **I read a peak at a grid node, again.** V4 failed on three populations as first coded, because the "peak
  field" was the largest node value, which moves 1.3 to 2.0×10⁻⁴ under radial refinement while the profile
  agrees to 10⁻⁵. It is the mistake this stage records against stage 6's label. With the parabola's vertex
  the peak moves 1.4×10⁻⁷. Fixed in the code and written into the amendment.
* **I sized two negative controls on the warm prototype.** They could not be met for a cold annulus as
  worded. Found by running the driver's pieces on the real populations, stated in the amendment with the
  numbers, and corrected there before the campaign with no threshold changed.
* **I sized the mode window on the wrong system**, and declared a discreteness reading that cannot
  discriminate before saturation. Both are above, reported as they fell.
* **A file's hash depended on which worker finished first.** Two assemblies of the same finished tasks gave
  different bytes for two series files, because a dictionary kept the order the tasks completed in. Keys are
  now sorted; two assemblies are byte-identical in all 19 series files and in the results file.

## What this stage does not establish

* **Thirty periods.** B_mid and the warm annuli are quiet for thirty periods at 1,024 bodies. The cold
  ring's spectrum above shows what thirty periods cannot see: unstable roots with growth times of 285 to
  10⁵ periods at every m. Only the linear mode calculation around these populations — the transport of a
  perturbation along the unperturbed orbits, coupled to the memory response — can exclude a slow mode in a
  warm annulus. **It is declared from stage 6 and it has not been run.**
* **Where the threshold is**, or what sets it. Three temperatures, one width, one footprint, one pair of
  response times, one radius.
* Part E certifies the cold *ring's* strip. It says nothing about an annulus's spectrum, and the measured
  cold-annulus rate (3 to 4 T0) is not the ring's (8.58 T0).
* The exploratory growth rates are exploratory. A declared measurement needs a window sized on the system
  it measures, fixed beforehand.
* Writing proportional to mass is still a label with a declared mass; the dissipation reservoir is
  unidentified; the kernel is instantaneous in space; no energy or angular-momentum budget is carried.
* No observational comparison is made. The owner's handoff reserves the local databases for explicitly
  declared comparisons, and none is declared.

## Next, in the owner's order

Item 6 of the handoff: the linear mode calculation around these populations, whose first test is now
concrete — it must return an m = 2 root near 0.04 to 0.05 with a pattern speed near 0.03 for the cold
annuli, and nothing faster than about 0.005 for B_mid and the warm ones; then targeted full-state
perturbations, formation from an empty field, and the reciprocal accounting through growth and saturation.
Each gets its own protocol file, with its mode window sized on the annulus and not on the ring.
