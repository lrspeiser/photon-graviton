# RUT-1 stage 9: the corrected response calculation

Protocol `protocol-rut9.md` (d4a5646), pushed with the code before any stage 9 run. Library `warm_response.py`,
units of work `rut9_tasks.py`, driver `rut9.py`, archive `rut9-results.json`, series `rut9-series/`, suite job
`rut9_checks.py`. Stage 8's library, stage 7's population builder and the owner's sampler are used by import,
unchanged. **Stage 8's archive, report and checks job are not touched by this stage: L3 and L3b are failed
there, and `rut8_checks.py` still exits non-zero.** This stage's status is its own.

This is the owner's ruling on 3a80fec carried out as far as it can be without a new simulation of the live
model: the third ring-limit verification, as a numerical repair with each annulus resolved on its own; the same
checks on B13, B24 and the near-threshold branch; stage 8's forecasts recalculated beside the originals; and,
for B13's 22% shortfall, the boundary-derivative audit, the eigenvalue sensitivity and a time-domain check on
B13's own eigenfunction. The fourth thing asked for, the seeded full-state experiment, is stage 10.

The methods are established ones and are not claimed: the matrix method for the linear response of a
distribution function in action-angle variables, and Beyn's contour-integral eigensolver with an
argument-principle count ([P2024] and [Beyn2012] in the owner's prior-art audit). What is this project's is the
interaction they are applied to.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | **passes.** `rut9_checks.py` re-hashes the code, the protocol, stage 8's predictions file and every series file; confirms the archive was made from the committed protocol; recomputes from scratch gate B1 on one population, the widest ring-limit annulus with its exact-ring target and its root under the declared representation and under stage 8's node rule, the slowest resolved member of the family scan, the sensitivity of B13's root and the declared replay batch of the time-domain check, body for body; and re-derives all eight gates from the archived task results. 15 of 15 comparisons with zero differences. |
| **numerical verification** | **passes: 8 of 8 gates, every one with its negative control rejected.** This is this stage's own status. Stage 8's L3 and L3b are failed in stage 8's archive and stay so. |
| **scientific outcome** | archived as it fell. The omitted boundary terms move B13's growth rate by +0.006%; the repaired representation moves no population's root by more than 2.6×10⁻⁵ of |s|; with both, the ring limit is approached as +23.45%, +4.79%, +1.13%, +0.28%; B13's own eigenfunction, imposed at its predicted complex rate on 8 batches of 200,000 bodies of the simulated population, implies a shift of the growth rate of −0.30% ± 0.24%, where a uniform change of the response of −23% would be needed to account for the shortfall. **The 22% is not in the linear response of the smooth equilibrium, to the accuracy of these checks. What it is, this stage cannot say.** |

## The correction: what the sharp bounds contribute

The population that is built, drawn and simulated is exactly zero outside |L − L0| ≤ 4dL and E − E_c(L) ≤ 4dE.
Stage 8's `warm_modes.py` differentiates the Gaussian inside those bounds and omits what the bounds themselves
contribute — an approximation I made knowingly in its prototype, judged negligible without calculating, and
did not disclose; the owner found it. With f0 = g Θ(E_max(L) − E) Θ(L − L_lo) Θ(L_hi − L) and
dE_max/dL = Ω_c(L), the response gains three families of boundary rows, each a one-dimensional integral over
orbits on the bound:

| family | F_l on it | measure |
|---|---|---|
| energy edge, E = E_c(L) + 4dE | (m Ω_c − ν_l) g | (2π)² dL / Ω_r |
| lower L edge | + m g | (2π)² dE / Ω_r |
| upper L edge | − m g | (2π)² dE / Ω_r |

The lower energy bound is J_r = 0, which nothing crosses, and carries no term. The cutoff is not replaced by a
taper anywhere in the calculation.

**Gate B1 — without the rows, the new library is stage 8's.** On stage 8's quadrature and nodes:

| population | orbit arrays | max \|T − T_stage8\| | archived root recovered to | with the rows, max \|ΔT\| (the control) | the rows move the root by | growth rate |
|---|---|---|---|---|---|---|
| A_cold | identical | 1.2×10⁻¹⁶ | 6.9×10⁻¹⁸ | 6.3×10⁻⁴ | +2.15e-06 -3.21e-05i | **+0.004%** |
| B_cold | identical | 4.9×10⁻¹⁷ | 6.9×10⁻¹⁸ | 3.4×10⁻⁴ | +2.16e-06 -1.89e-05i | **+0.005%** |
| B13 | identical | 3.0×10⁻¹⁷ | 6.9×10⁻¹⁸ | 2.1×10⁻⁴ | +2.00e-06 -1.29e-05i | **+0.006%** |
| B16 | identical | 3.9×10⁻¹⁷ | 1.4×10⁻¹⁷ | 1.3×10⁻⁴ | +2.94e-06 -8.97e-06i | **+0.015%** |

The rows' norm is 2.8×10⁻⁴ to 1.0×10⁻³ of the interior's for the energy edge and 2.4×10⁻⁵ to 4.0×10⁻⁵ for the two angular-momentum edges; almost all of the movement is in the frequency and almost all of it comes from the energy edge. **For the populations the omitted terms were negligible in the growth rate.** That is a result of calculating them; it was not known before.

**Gate B2 — the rows are the limit of a resolved taper.** Every bound replaced by a quintic smoothstep, on dedicated panels, with ∂f/∂E and ∂f/∂L taken as central differences of the function so that none of the algebra above enters. The tapered response minus the sharp interior, as a fraction of the rows' norm:

| population | bounds tapered | width 0.08 | width 0.04 | width 0.02 | Richardson limit, against the rows | the same limit against the interior alone (the control) |
|---|---|---|---|---|---|---|
| A_cold | energy | 1.84% | 0.87% | 0.43% | **1.6×10⁻⁵** | 100.0% |
| A_cold | energy, lower L, upper L | 1.91% | 0.91% | 0.44% | **1.9×10⁻⁵** | 100.0% |
| A_cold | lower L | 1.67% | 0.77% | 0.37% | **6.1×10⁻⁵** | 100.0% |
| A_cold | upper L | 6.46% | 3.06% | 1.49% | **5.6×10⁻⁵** | 100.0% |
| B13 | energy | 0.93% | 0.45% | 0.22% | **1.2×10⁻⁵** | 100.0% |
| B13 | lower L | 7.41% | 4.13% | 2.20% | **9.5×10⁻⁴** | 100.0% |
| B13 | upper L | 5.10% | 2.41% | 1.17% | **1.8×10⁻⁴** | 100.0% |
| B13 | energy, lower L, upper L | 1.12% | 0.56% | 0.28% | **7.0×10⁻⁵** | 100.0% |

## Part N. The numerical representation

**Gate N1 — every annulus resolved on its own (passes).** Relative movement of the m = 2 root, |Δs|/|s|, under each of the nine declared refinements; the bound is 10⁻⁴ each and 5.0×10⁻⁴ for their sum, a tenth of the ring-limit tolerance. The control is stage 8's node rule (spacing 0.1, margin 0.1) put through the same three kernel refinements; it must be caught on annuli 2, 3, 4 of 4, where it has three or four nodes.

| | dL 0.04, width 0.0279 | dL 0.02, width 0.0133 | dL 0.01, width 0.0066 | dL 0.005, width 0.0033 |
|---|---|---|---|---|
| kernel nodes, basis size | 18, 13 | 15, 12 | 14, 11 | 13, 11 |
| m = 2 root | 0.04467014 − 0.01235641i | 0.03452302 − 0.00519212i | 0.03258567 − 0.00435633i | 0.03213527 − 0.00417770i |
| kernel spacing 0.025 | 3.1×10⁻¹⁰ | 4.1×10⁻¹¹ | 3.0×10⁻⁹ | 3.4×10⁻⁹ |
| kernel margin 0.5 | 3.4×10⁻¹⁰ | 1.9×10⁻¹⁰ | 5.4×10⁻¹⁰ | 3.3×10⁻¹⁰ |
| centres shifted by half a spacing | 5.7×10⁻¹⁰ | 1.6×10⁻¹⁰ | 2.6×10⁻⁹ | 1.6×10⁻¹⁰ |
| rank cutoff 10⁻¹² | 3.6×10⁻¹⁰ | 2.2×10⁻¹⁰ | 3.1×10⁻⁹ | 3.2×10⁻¹⁰ |
| rank cutoff 10⁻⁸ | 9.1×10⁻⁹ | 8.4×10⁻⁸ | 6.7×10⁻⁹ | 3.0×10⁻⁷ |
| orbit quadrature 288 × 96 | 8.5×10⁻¹² | 1.8×10⁻¹¹ | 1.0×10⁻¹¹ | 2.1×10⁻¹¹ |
| 256 angle nodes | 7.5×10⁻¹³ | 1.5×10⁻¹² | 1.8×10⁻¹¹ | 1.6×10⁻¹¹ |
| harmonics to ±12 | 0 | 0 | 0 | 0 |
| the population's own grid refined | 3.4×10⁻⁶ | 1.1×10⁻⁶ | 1.8×10⁻⁶ | 1.9×10⁻⁶ |
| **sum of the nine** | **3.4×10⁻⁶** | **1.1×10⁻⁶** | **1.8×10⁻⁶** | **2.2×10⁻⁶** |
| the mode equation between the nodes | 2.6×10⁻¹⁰ | 9.6×10⁻¹¹ | 1.0×10⁻⁹ | 1.4×10⁻¹⁰ |
| *control*: nodes | 5 | 4 | 3 | 3 |
| *control*: largest movement under its own kernel refinements | 1.6×10⁻⁴ | 6.0×10⁻³ | 4.1×10⁻³ | 7.5×10⁻³ |
| *control*: the mode equation between its nodes | 7.3×10⁻⁶ | 2.2×10⁻⁴ | 2.7×10⁻³ | 3.3×10⁻³ |

The harmonics move the annuli's roots by exactly nothing because on orbits this nearly circular every row beyond l = ±8 falls below the 10⁻¹⁵ pruning threshold, and the matrices are identical: that variant is a test on the populations (below), not here.

**Gate N2 — the ring limit (passes).** Against the exact-history ring of stage 5 at the annulus's own mean radius and total writing rate:

| dL, dE | radial width | annulus, m = 2 | exact ring | growth rate | frequency | ring at R = 1 instead (the control) |
|---|---|---|---|---|---|---|
| 0.04, 0.004 | 0.02793 | 0.04467014 − 0.01235641i | 0.03618503 − 0.00604944i | **+23.449%** | -6.31e-03 | +112% |
| 0.02, 0.001 | 0.01335 | 0.03452302 − 0.00519212i | 0.03294411 − 0.00451168i | **+4.793%** | -6.80e-04 | +79% |
| 0.01, 0.00025 | 0.00660 | 0.03258567 − 0.00435633i | 0.03222170 − 0.00421366i | **+1.130%** | -1.43e-04 | +73% |
| 0.005, 6.25e-05 | 0.00329 | 0.03213527 − 0.00417770i | 0.03204621 − 0.00414354i | **+0.278%** | -3.42e-05 | +71% |

The differences fall at every step and the last is 0.278% against a bound of 0.5%. The observed orders between successive annuli, a diagnostic with no bound, are 2.15, 2.05, 2.02.

**Why stage 8's sequence bent.** Two defects, both repaired here. The under-resolved basis: N1's control shows stage 8's rule misplacing the root by 2.6×10⁻⁴, 4.1×10⁻³ and 7.4×10⁻³ of |s| on the three narrow annuli. And the omitted boundary rows: in the scratch prototype — not archived, and stated as such — the same four annuli without the rows read +23.37%, +4.70%, +1.03% and +0.18%, heading for −0.10% and not for zero, because the rows' share of the response does not shrink as the annulus narrows. With both repaired the observed order settles towards 2, which is what a finite-width correction even in the widths would give. I state that as an observation. It is not a bound, and the gate does not use it.

**Gate N3 — every population root resolved (passes).** The same table for the populations; the control is a deliberately coarse rule (spacing 0.3, margin 0.1).

| | `A_cold` | `B_cold` | B13 | B16 | dE 0.017 | dE 0.018 |
|---|---|---|---|---|---|---|
| kernel nodes, basis size | 24, 16 | 25, 17 | 27, 18 | 33, 21 | 35, 22 | 37, 23 |
| m = 2 root | 0.05087809 − 0.05203586i | 0.04131699 − 0.07746983i | 0.03216180 − 0.08786641i | 0.02022862 − 0.09795266i | 0.01558459 − 0.10132103i | 0.01058554 − 0.10475641i |
| kernel spacing 0.025 | 1.5×10⁻¹¹ | 2.3×10⁻¹² | 1.3×10⁻¹² | 2.0×10⁻¹² | 2.1×10⁻¹² | 2.2×10⁻¹² |
| kernel margin 0.5 | 1.5×10⁻¹¹ | 2.8×10⁻¹² | 1.2×10⁻¹² | 6.5×10⁻¹³ | 6.6×10⁻¹³ | 6.7×10⁻¹³ |
| centres shifted by half a spacing | 6.4×10⁻¹¹ | 4.4×10⁻¹¹ | 2.1×10⁻¹² | 5.0×10⁻¹² | 5.7×10⁻¹² | 6.4×10⁻¹² |
| rank cutoff 10⁻¹² | 1.6×10⁻¹¹ | 5.8×10⁻¹² | 2.9×10⁻¹² | 1.2×10⁻¹² | 1.1×10⁻¹² | 1.0×10⁻¹² |
| rank cutoff 10⁻⁸ | 2.0×10⁻⁹ | 1.6×10⁻⁹ | 1.0×10⁻⁹ | 6.4×10⁻¹⁰ | 5.8×10⁻¹⁰ | 5.4×10⁻¹⁰ |
| orbit quadrature 288 × 96 | 4.0×10⁻¹⁰ | 7.2×10⁻¹⁰ | 1.9×10⁻⁹ | 3.5×10⁻¹⁰ | 3.2×10⁻⁹ | 1.9×10⁻⁷ |
| 256 angle nodes | 9.9×10⁻¹² | 2.5×10⁻¹¹ | 2.2×10⁻¹¹ | 1.6×10⁻¹¹ | 2.9×10⁻¹¹ | 1.4×10⁻¹⁰ |
| harmonics to ±12 | 3.4×10⁻¹⁶ | 5.0×10⁻¹⁶ | 1.4×10⁻¹⁴ | 1.8×10⁻¹¹ | 5.5×10⁻¹¹ | 1.3×10⁻¹⁰ |
| the population's own grid refined | 2.6×10⁻⁵ | 1.9×10⁻⁵ | 2.1×10⁻⁵ | 2.3×10⁻⁵ | 2.4×10⁻⁵ | 2.5×10⁻⁵ |
| **sum of the nine** | **2.6×10⁻⁵** | **1.9×10⁻⁵** | **2.1×10⁻⁵** | **2.3×10⁻⁵** | **2.4×10⁻⁵** | **2.5×10⁻⁵** |
| the mode equation between the nodes | 5.7×10⁻¹¹ | 2.2×10⁻¹¹ | 1.7×10⁻¹¹ | 1.5×10⁻¹¹ | 1.7×10⁻¹¹ | 1.9×10⁻¹¹ |
| *control*: nodes | 3 | 4 | 4 | 5 | 5 | 6 |
| *control*: largest movement under its own kernel refinements | 1.8×10⁻¹ | 1.8×10⁻¹ | 1.3×10⁻¹ | 9.4×10⁻² | the root is lost | the root is lost |
| *control*: the mode equation between its nodes | 1.4×10⁻¹ | 4.6×10⁻² | 1.0×10⁻¹ | 1.5×10⁻¹ | 1.9×10⁻¹ | — |

**Gate N4 — the root-free conclusions survive refinement (passes).** 100 strips, 1,320 rectangles, for dE 0.019, B_mid, B24, A_warm, B_warm: the base representation in both partitions and, in the first, the basis, quadrature and harmonic refinements. Winding numbers summed over all of them: 0; roots located: 0; rectangles in which the count and the location disagree: 0 strips. The control — the m = 2 strip of dE 0.018, the slowest root stage 8 resolved — finds 0.010586 − 0.104756i (base); 0.010586 − 0.104756i (basis); 0.010586 − 0.104756i (quadrature); 0.010586 − 0.104756i (harmonics), agreeing to 1.9×10⁻⁷. **This is a statement about the declared domain — m = 1 to 4, growth rates of at least 0.006, the declared frequencies — and not a stability result.**

## Part F. Stage 8's forecasts, recalculated beside the originals

**The forecasts of record are stage 8's, pushed in aaf0bfe before anything was simulated.** What follows tests
their robustness and replaces nothing.

| dE | forecast of record (aaf0bfe) | recalculated | without the boundary rows | growth rate against the forecast |
|---|---|---|---|---|
| 0.010 | 0.041315 − 0.077451i | 0.04131699 − 0.07746983i | 0.04131483 − 0.07745090i | +0.0052% |
| 0.011 | 0.038562 − 0.080961i | 0.03856402 − 0.08097761i | 0.03856219 − 0.08096102i | +0.0048% |
| 0.012 | 0.035511 − 0.084429i | 0.03551322 − 0.08444356i | 0.03551140 − 0.08442896i | +0.0052% |
| 0.013 | 0.032160 − 0.087854i | 0.03216180 − 0.08786641i | 0.03215980 − 0.08785351i | +0.0062% |
| 0.014 | 0.028502 − 0.091239i | 0.02850394 − 0.09125019i | 0.02850166 − 0.09123874i | +0.0080% |
| 0.015 | 0.024528 − 0.094596i | 0.02453033 − 0.09460582i | 0.02452772 − 0.09459567i | +0.0106% |
| 0.016 | 0.020226 − 0.097944i | 0.02022862 − 0.09795266i | 0.02022568 − 0.09794369i | +0.0146% |
| 0.017 | 0.015581 − 0.101313i | 0.01558459 − 0.10132103i | 0.01558132 − 0.10131317i | +0.0210% |
| 0.018 | 0.010582 − 0.104750i | 0.01058554 − 0.10475641i | 0.01058199 − 0.10474961i | +0.0336% |
| 0.019 | — | — | — | — |
| 0.020 | — | — | — | — |
| 0.021 | — | — | — | — |
| 0.022 | — | — | — | — |
| 0.023 | — | — | — | — |
| 0.024 | — | — | — | — |

The zero-growth crossing, **extrapolated** from the resolved roots with Re s < 0.03 — nothing is resolved below a growth rate of 0.006, and no stability boundary is measured by this:

| polynomial | through the forecasts of record | through the recalculated roots |
|---|---|---|
| straight line | 0.02044 | 0.02044 |
| quadratic | 0.01993 | 0.01993 |
| cubic | 0.01991 | 0.01991 |

## Part S. What the root is sensitive to

| population | root | +1% of the whole response moves it by | an actual re-solve | first-order error | with H(s) held fixed in T′ (the control) | growth rate per 1% |
|---|---|---|---|---|---|---|
| A_cold | 0.050878 − 0.052036i | +2.705e-04 -1.586e-04i | +2.696e-04 -1.581e-04i | 0.34% | off by 175% | +0.53% |
| B_cold | 0.041317 − 0.077470i | +2.839e-04 -1.472e-04i | +2.829e-04 -1.467e-04i | 0.35% | off by 121% | +0.69% |
| B13 | 0.032162 − 0.087866i | +3.065e-04 -1.459e-04i | +3.054e-04 -1.454e-04i | 0.37% | off by 111% | +0.95% |
| B16 | 0.020229 − 0.097953i | +3.372e-04 -1.435e-04i | +3.359e-04 -1.431e-04i | 0.39% | off by 100% | +1.67% |

B13's mode is carried by the l = −1 radial harmonic: +1% of it moves the root by 3.41e-04, the next largest (l = 0) by 4.3e-06. It is the combination 2Ω_θ − Ω_r, about 1.1 for these orbits against a mode frequency of 0.09: nowhere near a resonance, which is why the orbit quadrature converges so easily and why a floor set by resonant denominators does not bear on this root. +1% in τ_keep, in H(s) alone, moves the growth rate by +0.13% and +1% in τ_form by −0.52%. **A uniform change of the response of −23% would be needed to lower B13's rate by 22%.**

## Part E. B13's own eigenfunction, in the time domain

Gate E1 passes. Bodies of B13 drawn with the owner's sampler move in the frozen potential plus ε e^{γt} Re[δC_m(r) e^{i(mθ+βt)}], δC_m the predicted eigenfunction and γ + iβ the predicted root, ε = 0.0005, from 8 e-foldings before t = 0, with +ε and −ε on identical bodies. The kernel-weighted second harmonic of the difference, projected on the left eigenvector, is the root the sampled, sharply truncated population would have, to first order.

| batch | bodies | implied shift of the growth rate | of the frequency | worst node | with the prediction scaled by 0.77 (the control) |
|---|---|---|---|---|---|
| 0 | 200,000 | +0.48% | -1.7e-04 | 0.95% | +22.4% |
| 1 | 200,000 | −1.20% | -8.3e-06 | 1.03% | +20.7% |
| 2 | 200,000 | −0.72% | +1.5e-04 | 0.71% | +21.2% |
| 3 | 200,000 | +0.37% | +2.8e-04 | 1.15% | +22.3% |
| 4 | 200,000 | −0.51% | -1.2e-04 | 1.51% | +21.4% |
| 5 | 200,000 | −0.50% | +8.8e-05 | 1.23% | +21.4% |
| 6 | 200,000 | −0.90% | -5.6e-05 | 0.97% | +21.0% |
| 7 | 200,000 | +0.55% | +1.2e-04 | 1.04% | +22.5% |
| **mean** | | **−0.30% ± 0.24%** | +3.6e-05 ± 5.4e-05 | 1.07% | +21.6% |

On ONE draw of 50,000 bodies, so that the differences carry no sampling noise: as declared +0.144%; with half the step +0.132% (a change of 0.0120% of the rate); with twice the amplitude +0.071% (a change of 0.073%, against a bound of 0.2%). The regenerated field over the imposed one, averaged over the nodes that carry the mode, is 0.9957 -0.0003i.

## What this says about B13's shortfall, and what it does not

B13 grew 22% more slowly in stage 8's simulations than its forecast, at 5.4 standard errors. This stage tested, on the forecast side, everything the owner named that needs no new simulation:

| candidate | what it does to B13's growth rate |
|---|---|
| the boundary terms of the sharp truncation, which stage 8 omitted | +0.006% |
| the kernel basis, its coverage, placement and rank; the orbit quadrature; the angle nodes; the harmonics | under 1.9×10⁻⁹ of \|s\| each |
| the population's own discretization | 2.1×10⁻⁵ of \|s\| |
| the linear response of the simulated, sampled, sharply truncated population to the mode itself | −0.30% ± 0.24% |
| *for scale:* a uniform change of the whole response that would produce −22% | −23% |

**The largest of them, taken with two standard errors, is 0.78% — a factor of 28 short of the shortfall.** The forecast of record, 0.03216, is what the linear theory of the smooth equilibrium gives, to 2.1×10⁻⁵ of |s| numerically and to the accuracy of the time-domain check physically. What this does **not** say is that the live model grows at that rate: a simulation differs from the linear theory of a smooth equilibrium by its finite number of bodies, by the noise its mode has to be read against, by its background's own slow evolution and by nonlinearity, and none of those is tested here. That is stage 10, and no outcome of it is assumed.

## What I got wrong on the way

* **The omitted boundary terms**, above: an approximation of stage 8 that I made knowingly and did not disclose,
  found by the owner. They turn out not to matter to the populations' growth rates and to matter at 0.1% to
  the ring limit — which I could not have known without calculating them, and did not calculate.
* **My first time-domain test measured its own amplitude.** Six prototype batches at ε = 0.002 implied a shift of
  −0.63 ± 0.19%, every batch negative, and I first wrote that down as a small offset of unknown origin. On one
  draw of bodies the implied shift is exactly quadratic in ε: it was the finite amplitude of the test. The
  declared ε is a quarter of it. I found this before declaring and not after, because the refinement runs were
  made on identical bodies, where a difference carries no sampling noise.
* **A longer switch-on is not a refinement of that kind**, and I nearly declared it as one: the same bodies
  evolved for longer reach t = 0 in other orbital phases, which is another sample of the response, and the
  result moves by a batch's noise.
* **My first negative control for the sensitivity formula was not a control.** I meant to show that the right
  eigenvector in place of the left one gives the wrong answer. T(s) is complex symmetric, so the left vector is
  the conjugate of the right, the right vector is nearly real, and the "wrong" formula is right to 0.4%. The
  declared control holds the memory transfer function fixed in dT/ds instead, which is off by more than 100%.
* **A control that was barely rejected**: a coarse node rule 0.2 apart moved B13's root by 1.2×10⁻⁴ against a
  bound of 10⁻⁴. A control that close to its bound tests nothing; the declared one is 0.3 apart and moves it by 2%.
* **The guard on my push did not guard.** The push of the protocol was rejected because the owner had pushed
  a commit of their own in the meantime, and my command piped the push through `tail`, whose exit status is
  what the rest of the chain saw: the local branch moved although the remote had not. Nothing was forced and
  nothing was lost — the one unpushed commit was rebased onto the owner's and pushed as a fast-forward — but the
  pattern I have used since stage 3 had this hole in it all along.

## What this stage does not establish

* **Whether the live model's B13 grows at the predicted rate.** This stage shows what the shortfall is not. The
  seeded, full-state experiment is stage 10.
* **Stability of anything.** Part N4 says that five populations hold no root with a growth rate of at least 0.006
  for m = 1 to 4 in the declared frequencies, and that the statement survives refinement of the basis, the
  quadrature and the harmonics. An e-folding of 26.5 periods is the floor; slower modes, other angular numbers
  and every damped mode are outside it.
* **A measured stability boundary.** The zero-growth crossing is an extrapolation of the resolved branch, and
  its value depends on the polynomial at the 3% level.
* **Anything about stage 8's status.** L3 and L3b failed and stay failed in stage 8's archive. This stage's
  ring-limit gate is a new gate on a corrected calculation with its own status.
* Writing proportional to mass is still a label with a declared mass; the dissipation reservoir is
  unidentified; the kernel is instantaneous in space; the bodies' lost angular momentum has no ledger; and no
  observational comparison is made.

## Next, in the owner's order

Stage 10, the seeded full-state experiment, under its own protocol: a four-fold symmetric sample, which the
simulator's grid keeps free of m = 2 noise, seeded with the predicted eigenmode in bodies, excitation and field,
at more than one amplitude, in several sampling realizations and with the step and the grid refined. Then the
quiet region's scope — the resonant denominators treated properly below the present floor, more angular
numbers, controlled disturbances — and then formation from an empty field with the reciprocal energy and
angular-momentum budget carried through. The owner's prior-art programme (`research_plan/prior-art/`) asks,
after these, for the published response-matrix and contour-solver benchmarks to be reproduced with this code;
that is the right external check of parts B and N and it has not been made.
