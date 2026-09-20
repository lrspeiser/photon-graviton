# RUT-1 stage 10: the seeded, full-state warm-mode experiment

Protocol `protocol-rut10.md` (7762282), pushed with the code before any stage 10 run. Units of work
`rut10_tasks.py`, driver `rut10.py`, archive `rut10-results.json`, series `rut10-series/`, suite job
`rut10_checks.py`. Stage 9's corrected response calculation, stage 7's simulator and population builder and the
owner's sampler are used by import, unchanged.

This is the fourth of the four things the owner asked for on B13's 22% shortfall, the central discrepancy of
stage 8. Stage 9 ran the other three and found that the shortfall is not an error of the linear calculation.
The owner's question for this one: **does the measured rate approach the predicted rate as the seed amplitude
and the numerical error decrease?**

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | **passes.** `rut10_checks.py` re-hashes the code, both protocols, stage 8's predictions, stage 9's archive and every series file; **replays the exact preparation of one seeded run** — the quiet sample drawn again, moved again by the predicted mode's field, the field seeded again — and finds the bodies and the complete field identical to the saved ones, bit for bit; evolves a declared prefix of that run and of the unseeded quiet run and compares them with the series record for record; and re-derives every fit, gate, diagnostic and reading. 14 of 14 comparisons with zero differences. |
| **numerical verification** | **passes: 6 of 6 gates.** |
| **scientific outcome** | archived as it fell. **B13, seeded with its predicted eigenmode in a quiet sample, grows at the forecast rate: 0.03180 ± 0.00109 at 4,096 bodies against 0.03216 (−1.1% ± 3.4%), with pattern speed +0.04388 ± 0.00060 against +0.04393**; 0.03126 ± 0.00159 at 1,024. By the declared rule that is *reproduced quantitatively* (the archive's word for it is `confirmed`; see below). Stage 8's unseeded reading was 0.0252 ± 0.0013. The rate does not depend on the seed's amplitude or sign, on the step or on the grid; what limits it is the spread between sampling realizations, 14% at 1,024 bodies and 7% at 4,096. B16: 0.02155 ± 0.00085 against 0.02023, +6.6% ± 4.2% — inside the rule through its three-standard-error clause and outside its 5% clause, which is compatibility and not a precise test. **The 22% shortfall does not appear when the mode is seeded in a quiet sample.** That places it in the unseeded measurement — a mode read out of shot noise of its own size — or in what shot noise does to the population; this stage does not separate the two. |

## The construction

**A quiet sample.** N/4 bodies from the owner's sampler, the sampler of record, each joined by its images under
rotation by 90°, 180° and 270°, which are exact in floating point. The population is axisymmetric, so this is a
sample of it; its m = 1, 2 and 3 harmonics vanish identically; and the simulator's Cartesian field is symmetric
under the same rotations, so the live system keeps the symmetry and an m = 2 pattern can grow only from what is
put there. It is a quieter sampling of the azimuth and of nothing else: the sample still holds N/4 independent
orbits, and that is what limits the measurement.

**The seed is the predicted eigenmode in the full state.** The bodies carry δf because they have been moved by
the mode's own field, imposed at its predicted complex rate from eight e-foldings before t = 0 — the first-order
response of this sample of this sharply truncated population, which is the operation stage 9's E1 verified.
The force-producing field carries δC with its analytic gradient. The excitation carries
δE = (s + 1/τ_keep) δC, which is what a mode growing as e^{st} requires. Every run's initial bodies and the
complete seeded fields are in the series, with the mode's coefficients, nodes and root.

**Q1 the sample (passes).** 16 quiet samples: the largest departure of the mean radius, the mean angular momentum or the radial velocity variance from the owner's node distribution is 2.41 standard errors of the N/4 independent bodies (bound 4); the largest |S_m|/S₀ for m = 1, 2, 3 is 1.8×10⁻¹⁷ (bound 10⁻¹²). The control, an ordinary draw of the same size, carries an m = 2 source of at least 4.7×10⁻³.

**Q2 the symmetry survives the live run (passes).** Unseeded, quiet, 40 T0: the largest |C_m|/C₀ at any record is 5.6×10⁻¹⁶ (m = 1), 9.0×10⁻¹⁵ (m = 2), 2.5×10⁻¹⁶ (m = 3) against a bound of 10⁻¹⁰. The control, an unseeded ordinary draw, reaches 5.1×10⁻⁴ within 8 T0.

**P1 the seeded field (passes).** Target 10⁻⁶: |C₂|/C₀ at t = 0 is 1.000000e-06; E₂/C₂ = 0.048077 -0.087867i against s + 1/τ_keep = 0.048077 -0.087866i. The control, the equilibrium relation δE = δC/τ_keep, gives 0.015915 -0.000000i.

## The growth of the seeded mode

**The window, declared before any run:** from the first record at which |C₂|/C₀ exceeds three times its value
at t = 0 to the last record before it first exceeds 2×10⁻⁴.

**G1 linear in the seed (passes).** While the +10⁻⁵ run is under the ceiling (57 records) the ratio of its amplitude to the +10⁻⁶ run's departs from its initial value by at most 0.09% (bound 1%); the −10⁻⁶ run has the +10⁻⁶ run's amplitude to 6.7×10⁻¹¹ and a phase π away to 4.2×10⁻¹¹. The control: followed above the ceiling to 2.5×10⁻³, the same ratio falls by 23%.

**G2 the step and the grid (passes).** Halving the step changes the first realization's rate by 0.004% and a grid of w/8 by 0.017% (bound 0.5%). The control, a grid of w/1.25, changes it by 7.6%.

| population | bodies | realization | window (T0) | e-foldings | growth rate | pattern speed | rms residual of log amplitude |
|---|---|---|---|---|---|---|---|
| B13 | 1,024 | 0 | 5.50 – 23.75 | 4.11 | 0.03563 | +0.04409 | 0.025 |
| B13 | 1,024 | 1 | 5.25 – 25.75 | 4.16 | 0.03183 | +0.04529 | 0.030 |
| B13 | 1,024 | 2 | 6.00 – 28.75 | 4.14 | 0.03019 | +0.04414 | 0.056 |
| B13 | 1,024 | 3 | 6.50 – 31.50 | 4.16 | 0.02565 | +0.04420 | 0.046 |
| B13 | 1,024 | 4 | 5.00 – 23.75 | 4.12 | 0.03499 | +0.04148 | 0.030 |
| B13 | 1,024 | 5 | 5.50 – 25.00 | 4.15 | 0.03402 | +0.04204 | 0.024 |
| B13 | 1,024 | 6 | 5.75 – 24.75 | 4.15 | 0.03420 | +0.04455 | 0.033 |
| B13 | 1,024 | 7 | 6.25 – 33.00 | 4.16 | 0.02357 | +0.04738 | 0.089 |
| B13 | 4,096 | 0 | 5.75 – 28.25 | 4.15 | 0.02900 | +0.04321 | 0.016 |
| B13 | 4,096 | 1 | 5.00 – 26.25 | 4.19 | 0.03130 | +0.04539 | 0.010 |
| B13 | 4,096 | 2 | 5.25 – 24.50 | 4.12 | 0.03409 | +0.04422 | 0.004 |
| B13 | 4,096 | 3 | 5.50 – 25.50 | 4.12 | 0.03280 | +0.04268 | 0.004 |
| B16 | 4,096 | 0 | 9.25 – 37.75 | 4.17 | 0.02318 | +0.04902 | 0.042 |
| B16 | 4,096 | 1 | 9.00 – 38.00 | 4.14 | 0.02272 | +0.04663 | 0.070 |
| B16 | 4,096 | 2 | 8.75 – 41.00 | 4.15 | 0.02075 | +0.04949 | 0.098 |
| B16 | 4,096 | 3 | 10.25 – 44.25 | 4.17 | 0.01955 | +0.05440 | 0.067 |

## The declared reading

Against the forecast of record — stage 8's, pushed in aaf0bfe before anything was simulated.

| population | bodies | realizations | growth rate | forecast of record | difference | pattern speed | forecast | difference | reading |
|---|---|---|---|---|---|---|---|---|---|
| B13 | 1,024 | 8 | **0.03126 ± 0.00159** | 0.03216 | −2.80% ± 4.94% | **+0.04415 ± 0.00065** | +0.04393 | +0.50% ± 1.47% |  |
| B13 | 4,096 | 4 | **0.03180 ± 0.00109** | 0.03216 | −1.13% ± 3.40% | **+0.04388 ± 0.00060** | +0.04393 | −0.12% ± 1.36% | **confirmed** |
| B16 | 4,096 | 4 | **0.02155 ± 0.00085** | 0.02023 | +6.55% ± 4.20% | **+0.04988 ± 0.00163** | +0.04897 | +1.86% ± 3.33% | **confirmed** |

**The archive says `confirmed` where the protocol says *reproduced quantitatively*.** The rule applied is the protocol's, word for word; the label in the driver is a leftover of stage 8's vocabulary that I did not change, and the archive is not edited to fix it. B13 meets the rule through its 5% clause in both the rate and the pattern speed, with a standard error of 3.4% of the forecast, inside the 5% power requirement (at 1,024 bodies it is 4.9%: barely). **B16 meets it only through the three-standard-error clause**: its mean is +6.6% from the forecast, 1.6 standard errors, with a standard error of 4.2%, just inside the power requirement. That is the same kind of pass as stage 8's B16 and I describe it the same way: compatible, not a precise confirmation. Nearer the end of the branch the rate is a steeper function of everything, and four realizations at 4,096 bodies are not enough.

## The owner's question

| run (first realization, 1,024 bodies) | growth rate | pattern speed | window (T0) |
|---|---|---|---|
| seed +10⁻⁶ | 0.03563 | +0.04409 | 5.50 – 23.75 |
| seed +10⁻⁵ | 0.03423 | +0.04327 | 5.50 – 14.00 |
| seed −10⁻⁶ | 0.03563 | +0.04409 | 5.50 – 23.75 |
| step halved | 0.03564 | +0.04409 | 5.50 – 23.75 |
| grid w/8 | 0.03564 | +0.04409 | 5.50 – 23.75 |
| grid w/1.25 (G2's control) | 0.03292 | +0.04552 | 5.75 – 25.75 |
| WRONG: field seeded, bodies not prepared | 0.03467 | +0.04372 | 7.75 – 26.50 |
| WRONG: δE = δC/τ_keep | 0.03610 | +0.04434 | 6.50 – 24.75 |

**Does the measured rate approach the predicted rate as the seed amplitude and the numerical error decrease?** It does not move with either, and it sits on the prediction within the spread of realizations. A seed ten times larger follows the same records to 0.09% while it is under the ceiling (its own window is shorter, which is why its fitted rate differs). The opposite sign is a rotation of the same run by 90° and agrees to rounding: it tests that the simulator is rotation-covariant, which the quiet start relies on, and is not an independent measurement. Halving the step and refining the grid change nothing; a grid coarse enough to matter is caught. Four times the bodies moves the mean from 0.03126 ± 0.00159 to 0.03180 ± 0.00109: towards the forecast 0.03216, and not significantly. **The two wrong preparations barely move the rate over this window** (−2.7% for bodies that were not prepared, whose window starts two periods later, and +1.3% for the equilibrium relation in the excitation): an instability finds its own mode, so the growth rate is not a sensitive test of how the seed was prepared. What tests the preparation is stage 9's E1; what the full-state seed buys here is a clean exponential from the first record.

## Four declared diagnostics, with no verdict

**Where the window starts.**

| population, bodies | window from | realizations | growth rate | pattern speed |
|---|---|---|---|---|
| B13, 1,024 | 3 × the seed | 8 | 0.03126 ± 0.00159 | +0.04415 ± 0.00065 |
| B13, 1,024 | 10 × the seed | 8 | 0.03170 ± 0.00152 | +0.04442 ± 0.00067 |
| B13, 1,024 | 30 × the seed | 7 | 0.03216 ± 0.00141 | +0.04477 ± 0.00086 |
| B13, 4,096 | 3 × the seed | 4 | 0.03180 ± 0.00109 | +0.04388 ± 0.00060 |
| B13, 4,096 | 10 × the seed | 4 | 0.03169 ± 0.00120 | +0.04394 ± 0.00065 |
| B13, 4,096 | 30 × the seed | 4 | 0.03178 ± 0.00117 | +0.04423 ± 0.00061 |
| B16, 4,096 | 3 × the seed | 4 | 0.02155 ± 0.00085 | +0.04988 ± 0.00163 |
| B16, 4,096 | 10 × the seed | 4 | 0.02204 ± 0.00123 | +0.04969 ± 0.00187 |
| B16, 4,096 | 30 × the seed | 4 | 0.02202 ± 0.00172 | +0.04974 ± 0.00195 |

**Stage 8's own rule on the seeded runs.** Stage 8's own growth-phase rule applied to the 8 seeded quiet runs of B13 at 1,024 bodies: window 6 – 34.25 T0, rate 0.03163 ± 0.00161, pattern speed +0.04419 ± 0.00065.

**Stage 8's unseeded runs read again with B24's shot noise as a floor.** B24 holds no unstable root, so its
m = 2 pattern is shot noise; its ensemble-mean power is a *proxy* for the noise under the others, which a
population nearer instability amplifies more. Rough, and declared so that it is made once:

| population, bodies | stage 8's window (T0) | stage 8's declared rate | rate of the total m = 2 power | with B24's power subtracted | B24's power over the total at the window's start |
|---|---|---|---|---|---|
| B13, 256 | 6 – 9.75 | 0.0365 ± 0.0068 | 0.0294 | **0.0342** | 24% |
| B13, 1,024 | 6 – 14.75 | 0.0252 ± 0.0013 | 0.0259 | **0.0284** | 46% |
| B16, 256 | 6 – 8.25 | 0.0346 ± 0.0171 | 0.0210 | **0.0296** | 36% |
| B16, 1,024 | 6 – 11.75 | 0.0270 ± 0.0071 | 0.0339 | **0.0400** | 57% |

**The live background.** The predicted root moves by +0.2152 +0.3683i per unit relative change of a uniformly scaled C₀: +0.1% in the equilibrium field is +0.67% in the rate and −0.42% in the pattern speed.

| run | live C₀ on the ring over the window, against its primed value | implied shift of the rate | of the pattern speed |
|---|---|---|---|
| B13:q256:k1:t+1e-06:base | +0.068% | +0.00015 | -0.00012 |
| B13:q256:k4:t+1e-06:base | −0.083% | -0.00018 | +0.00015 |
| B13:q256:k0:t+1e-06:base | −0.116% | -0.00025 | +0.00021 |
| B13:q256:k5:t+1e-06:base | −0.168% | -0.00036 | +0.00031 |
| B13:q256:k0:t+1e-05:base | −0.075% | -0.00016 | +0.00014 |
| B13:q256:k0:t-1e-06:base | −0.116% | -0.00025 | +0.00021 |
| B13:q256:k2:t+1e-06:base | −0.001% | -0.00000 | +0.00000 |
| B13:q256:k6:t+1e-06:base | −0.268% | -0.00058 | +0.00049 |
| B13:q256:k3:t+1e-06:base | −0.008% | -0.00002 | +0.00001 |
| B13:q256:k7:t+1e-06:base | −0.060% | -0.00013 | +0.00011 |
| B13:q1024:k1:t+1e-06:base | +0.055% | +0.00012 | -0.00010 |
| B13:q1024:k0:t+1e-06:base | −0.127% | -0.00027 | +0.00023 |
| B13:q1024:k3:t+1e-06:base | −0.089% | -0.00019 | +0.00016 |
| B13:q1024:k2:t+1e-06:base | +0.003% | +0.00001 | -0.00001 |
| B16:q1024:k2:t+1e-06:base | −0.131% | — | — |
| B16:q1024:k0:t+1e-06:base | +0.033% | — | — |
| B16:q1024:k1:t+1e-06:base | +0.152% | — | — |
| B16:q1024:k3:t+1e-06:base | +0.046% | — | — |

**The unseeded quiet population.**

| quantity (unseeded quiet run) | first 5 T0 | last 5 T0 | change |
|---|---|---|---|
| radial spread | 0.05589 | 0.05624 | +0.63% |
| residual radial dispersion | 0.10014 | 0.10069 | +0.55% |
| residual azimuthal dispersion | 0.05443 | 0.05441 | −0.03% |
| L mean | 1.00447 | 1.00446 | −0.00% |
| support on bodies | 0.10268 | 0.09844 | −4.14% |

*Where the window starts:* at 4,096 bodies nothing moves (0.03180, 0.03169, 0.03178; pattern speeds +0.04388, +0.04394, +0.04423). *Stage 8's rule* read the seeded quiet runs at 0.03163 ± 0.00161: the rule is not biased when there is no shot noise under it. *With B24's shot noise as a floor*, stage 8's unseeded B13 at 1,024 bodies moves from 0.0252 to 0.0284 — B24's power is 46% of the total where that window starts — which is part of the way to the seeded 0.0313; but the same subtraction takes B16 at 1,024 bodies to 0.0400, twice its forecast, so **the proxy is too crude to carry a conclusion**, and I draw none from it beyond this: the unseeded windows began where the mode and the noise were the same size. *The live background* drifted by −0.27% to +0.07% on the ring over the windows, worth at most 1.8% of the rate by the proxy of a uniform scaling: it does not explain the spread between realizations (B16's rows carry no implied shift because I declared the sensitivity for the gated population only). *The unseeded quiet population*, 40 T0 at 1,024 bodies: its radial spread changes by +0.63% and its residual radial dispersion by +0.55% between the first and the last five periods, where stage 8's ordinary draw of the same population gained 44% in radial spread against its frozen twin. **With m = 1, 2 and 3 removed by symmetry and m = 0 and 4 left in, B13 does not heat**: what damaged it in stage 8 was the m = 2 mode, as the linear calculation says. The mean support felt by the bodies falls by 4.1% over the same run; in two prototype runs it moved in opposite directions, it is a mean over 256 independent orbits of a quantity that varies strongly along each, and I read nothing into it.

## What I got wrong on the way

* **I reported a pattern-speed offset and a drift that were three realizations.** The protocol says, from a prototype of three runs at 1,024 bodies and one at 4,096, that the pattern speed was "3.6% high in every realization — and it drifts", and that I did not know why. With eight and four realizations the pattern speed is +0.5% ± 1.5% and −0.1% ± 1.4% from the forecast and does not drift. There was nothing to know: four runs on the same side is a one-in-eight event. I had stated it as an observation with its sample size and drawn no conclusion, which is the only reason it did no harm. The same prototype suggested a spread of 3% between realizations; it is 14%.
* **The archive's label for the reading is not the protocol's word.** `confirmed` in the archive is the protocol's *reproduced quantitatively*; the rule is the protocol's. A leftover name in the driver, caught only when the campaign printed it.
* **The live-background diagnostic has no implied shifts for B16**, because I declared the sensitivity task for the gated population only.
* **My first step refinement refined nothing**, and the prototype caught it before the declaration: the simulator's step is min(h_max, η√(r³/GM)), η is the binding limit, and halving h_max alone reproduced the base run to four figures for the wrong reason. The declared refinement halves both.
* **The first control I tried for the step-and-grid gate could not be caught**: steps four and eight times coarser change the rate by 0.06%, and a grid twice as coarse by 0.5%. The declared control is a grid of w/1.25, which changes it by 7.6%.

## What this stage does not establish

* **Why the unseeded reading of stage 8 was 22% low.** A mode read out of shot noise of its own size, or a population that shot noise has changed: not separated here. The rough noise-floor diagnostic accounts for part of it for B13 and fails for B16.
* **B16's rate to better than 4%.** Its reading passes on its error bar.
* **The limit of many bodies.** Two body counts, consistent with each other and with the forecast; the spread between realizations is the limit, and a quieter sampling of the actions, which this stage did not attempt, is what would reduce it.
* **Anything about an ordinary sample's noise**, its saturation, modes other than m = 2, or growth slower than the horizon can show. A four-fold symmetric sample is a device for measuring one mode; it is not the model's generic state.
* **Stability of anything.**
* Writing proportional to mass is still a label with a declared mass; the dissipation reservoir is unidentified; the kernel is instantaneous in space; the bodies' lost angular momentum has no ledger; and no observational comparison is made.

## Next, in the owner's order

The owner's milestone has three clauses: one corrected response calculation that passes its limiting and population-specific checks — stage 9; that reproduces a deliberately seeded warm-population mode quantitatively — by the declared rule this stage does, for B13, and only weakly for B16, and whether that meets the clause is the owner's judgement and not mine; and that identifies a supported quiet region with a stated timescale and mode-domain bound — not done. That is next: the resonant denominators treated properly below the present floor of 0.006, more angular numbers, and controlled disturbances of the populations that hold no root, for which the quiet start and the full-state seed built here are the tools. Then formation from an empty field with the reciprocal energy and angular-momentum budget carried through, and the owner's prior-art programme (`research_plan/prior-art/next-comparisons.md`), beginning with the published response-matrix and contour-solver benchmarks.
