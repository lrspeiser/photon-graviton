# RUT-1 stage 4: the two-stage response delays the failure and changes its character — the ring is unstable to a large-scale mode

[protocol-rut1.md](protocol-rut1.md), stage 4 declared in 8e2d681 before any of these runs. Twenty runs to 100 reference periods, the equation frozen, every threshold and every refinement fixed in advance. **All four gates pass.**

**The answer to the question stage 4 was built to ask.** The radial motion left in the two-stage runs is neither initial adjustment nor bounded breathing. It is a **large-scale collective instability** — angular modes m = 2 and m = 3–4 — that grows from a cold, mature, supported ring, saturates into a hot state, and carries essentially all of the angular-momentum loss. The fine structure at the writer spacing, which maturation was designed to suppress, stays suppressed by five to eight orders of magnitude throughout. Maturation solved the fast-wake problem and exposed a slow one, which is the character your transfer-function analysis predicted.

## Over 100 periods the one-stage model fails everywhere

| configuration | one-stage | two-stage |
|---|---|---|
| 16 writers, w/R 0.1 | unresolved centre at **10.0 T0** | completed 100 |
| 32 writers, w/R 0.2, seed 1 | unresolved centre at **27.2 T0** | completed 100 |
| 32 writers, w/R 0.2, seed 2 | unresolved centre at **38.1 T0** | completed 100 |
| 16/0.1 primed with a mature field | unresolved centre at **12.8 T0** | completed 100 |
| 32/0.2 primed with a mature field | unresolved centre at **35.0 T0** | completed 100 |

Both 32-writer one-stage rings had *completed* 20 periods in rut3; neither survives 40. A 20-period horizon was not long enough to see the outcome, which is why this stage exists. Every one-stage run's late window is therefore "not reached", and the no-memory controls are flat — mean radius constant to 2×10⁻⁴, dispersion 0.021, angular momentum 1.0000 — for all 100 periods.

## The instability, measured

The primed runs are the clean experiment: a cold ring launched in equilibrium with an already-mature field. The m = 2 band of the force-producing field grows exponentially out of the numerical floor:

| t (T0) | 10 | 30 | 45 | 55 | 65 | 100 |
|---|---|---|---|---|---|---|
| m = 2 band power, r = 0.9 | 2.2×10⁻¹⁰ | 2.1×10⁻⁸ | 7.2×10⁻⁷ | 6.5×10⁻⁶ | 1.0×10⁻⁵ | 4.8×10⁻⁶ |
| radial dispersion | 0.0069 | 0.0108 | 0.026 | 0.087 | **0.32** | 0.32 |

Power rises about a hundredfold per twenty periods — an **amplitude e-folding time of roughly 8.7 periods** — then saturates near 60 T0, at which point angular momentum drops from 0.993 to 0.852 and the ring settles into a hot, contracted, supported state. The 16-writer primed ring does the same on a faster clock, heating between 25 and 35 periods, and then **keeps losing angular momentum**: 0.844 → 0.831 → 0.813 → 0.797 → 0.775 through the late window.

**It is not a grid artefact.** Two independent checks. The declared grid-rotation control — the same initial condition turned 0.3 rad on the Cartesian grid — agrees with the unrotated run on every late-window quantity, all within tolerance:

| | unrotated | rotated | difference | tolerance |
|---|---|---|---|---|
| mean radius | 0.8752 | 0.8751 | 1.7×10⁻⁴ | 7.7×10⁻³ |
| radial dispersion | 0.2650 | 0.2787 | 1.4×10⁻² | 6.8×10⁻² |
| support | 0.1314 | 0.1341 | 2.7×10⁻³ | 7.1×10⁻³ |
| angular momentum | 0.9181 | 0.9075 | 1.1×10⁻² | 4.3×10⁻² |

And an **exploratory** cross-check (run after declaration, so not a gate): the stage 4R spectral implementation — a periodic Fourier box with exactly reciprocal deposition and readout, no Cartesian grid at all — grows the same mode from the same primed state, with m = 2 power rising ×10.6 then ×10.1 per ten periods against the grid model's ×10.5 and ×9.1, saturating at 59 T0 against 60, and losing the same angular momentum.

## Where the angular momentum goes, by angular mode

Net late-window torque summed over bodies, by band:

| run | total | m = 1 | m = 2 | m = 3–4 | m = 5–8 | m = 9–16 | m = 17–32 |
|---|---|---|---|---|---|---|---|
| two-stage 32/0.2 s1 | −0.576 | +0.009 | −0.060 | **−0.554** | +0.029 | −2×10⁻⁴ | +6×10⁻⁷ |
| primed two-stage 32/0.2 | −5.05 | +0.021 | **−4.51** | −0.555 | +0.007 | −0.011 | +3×10⁻⁴ |
| two-stage 16/0.1 s1 | −0.402 | +0.022 | +0.151 | **−0.346** | −0.230 | +0.001 | −0.001 |
| primed two-stage 16/0.1 | −1.298 | +0.120 | −0.015 | −0.396 | −0.420 | **−0.587** | +1×10⁻⁵ |

During the m = 2 saturation, that single band carries **89%** of the torque. Away from it, m = 3–4 dominates. The writer-spacing scale is invisible for the 32-writer rings (m = 17–32 at 10⁻⁴–10⁻⁷). One caveat: for the 16-writer primed run the largest single contribution sits in the m = 9–16 band, which *contains* that ring's own spacing mode m = 16 — these bands cannot separate m = 16 from m = 9–15, and finer bands would be needed to say which.

## The declared classification, and where it is conservative

| run | verdict | resolved flags |
|---|---|---|
| two-stage 16/0.1 s1 | formation adjustment, then stationary | positive support |
| two-stage 32/0.2 s1 | **candidate heating or instability** | secular migration, excursion growth, positive support |
| two-stage 32/0.2 s2 | formation adjustment, then stationary | positive support |
| primed two-stage 16/0.1 | **candidate heating or instability** | migration, angular-momentum drift, dispersion growth, excursion growth, positive support |
| primed two-stage 32/0.2 | **candidate heating or instability** | migration, excursion growth, positive support |
| no-memory ×3 | stationary | — |

Read this carefully, because the rule is deliberately strict. In the empty-field runs the late-window trends are real but often **unresolved**: the 32/0.2 seed-1 ring loses 0.0139 ± 0.0014 of its angular momentum over the late window, which is 10σ statistically, yet the numerical uncertainty measured from its own half-step and half-spacing refinements is 0.0057, and 3× that is 0.0171 — larger than the change. So it is not claimed. The refinements diverge from the base run because the saturated state is chaotic, not because the integration is poor, and that chaos is folded into the uncertainty honestly rather than argued away. **The primed runs carry statistical resolution only**, since numerical uncertainty was measured for the three empty-field configurations that were refined, and this is stated in their rows rather than hidden.

Support is positive and resolved in every two-stage run: **+0.091, +0.131, +0.132** late-window mean for the empty-field runs and +0.100, +0.136 for the primed ones — consistent with the 9.6–14.1% of your exploratory runs. The cost is equally clear: those same rings sit at mean radius 0.84–0.88, radial dispersion 0.27–0.29 against the no-memory control's 0.021, and angular momentum 0.86–0.92.

## Gates

| gate | requirement | measured |
|---|---|---|
| ε-consistency residual | < 10⁻³ relative, halving the step reduces it ≥2× | **6.0×10⁻⁴**; ratios 2.79, 2.14, 6.79 |
| torque identity L(t) − L(0) = ∫∂C/∂θ dt | < 10⁻¹⁰ | **8.8×10⁻¹⁴** |
| band attribution at the last checkpoint | torque 2%, work 5% | **0.12%**, **1.97%** |
| grid rotation | within 3× numerical uncertainty or the seed spread | all four quantities within |

The ε-residual is the check that the moving body, the sampled force, the sampled potential and the field update are mutually consistent; it holds to 6×10⁻⁴ even in the hot saturated states and converges under refinement. The torque identity holds to machine precision.

## What this does and does not establish

**Shown.** Over a horizon four times longer than rut3's, the one-stage model fails in every configuration, including two that had "completed" at 20 periods — a 20-period result was not a result. The two-stage response survives all 100 periods in every configuration and delivers 9–14% measured inward support, but it is **unstable to a large-scale m = 2 / m = 3–4 mode** that grows with an amplitude e-folding time near 9 periods even from a cold mature ring, saturates into a hot contracted state, and carries up to 89% of the angular-momentum loss in a single band. Fine structure at the writer spacing is genuinely suppressed. The instability survives a grid-rotation control and is reproduced by an independent spectral discretization.

**Not shown.** The mode's **pattern speed was not measured** — only band powers were recorded, not their phases — so the hypothesis that this is the slow phase-lagged response your transfer function predicts (quadrature ratio above one for variation slower than 34.4 T0) is **untested**. That is the single most informative next measurement, and it decides whether a longer formation time would help or hurt. Nothing here is a robustness map: one planar ring per configuration, three configurations, two seeds, a fixed maturation time, an instantaneous kernel, no field energy budget for the grid model, and a writing rate set from a label. A verdict classifies this frozen equation under these conditions.

## Reproduce

```bash
python research_work/results/path-memory/rut4.py --canonical   # tens of minutes on many cores
python research_work/results/path-memory/rut4_checks.py        # the suite job
```
