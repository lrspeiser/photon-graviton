# Extra event stretching does not close the brightness gap

The preceding goal turn tested deposited-particle support. This run returns to the first objective and links the existing brightness fit to an independently constrained event-stretch exponent. It makes no new parameter fit and exposes why changing time alone is not an adequate repair in the current Euclidean-beam branch.

## Known flux accounting connects the observables

Let S_E=1+z, S_t=(1+z)^b, and D=ln(1+z)/alpha for the retained propagation law. For conserved photon number, fixed source luminosity and Euclidean beam area,

    F=L/[4pi D^2 S_E S_t].

Our earlier optional photon-removal fit used S_t=S_E and survival P=(1+z)^(-epsilon), giving F=L/[4pi D^2(1+z)^(2+epsilon)]. These fluxes are identical if b=1+epsilon. This is an algebraic restatement of known energy/rate accounting, not a new physical mechanism. The equivalence includes the earlier calibrator attenuation adjustment when applied consistently.

The earlier 466-row brightness training fit gives epsilon=0.481189 with formal conditional error 0.102815. Reinterpreting all its extra dimming as timing therefore gives b=1.481189. This is not a measured timing exponent or a claim that all brightness observations require that value. It depends on the fixed alpha, distance convention, photon-count and beam assumptions, and on using standardized magnitudes whose reduction was not redone for this time law.

## Timing and frozen transfer

[The DES timing analysis](https://arxiv.org/abs/2406.05050) reports b=1.003 with statistical scale 0.005 and systematic estimate 0.010, conditional on its source-duration and light-curve methodology. We retain those distinct uncertainties rather than manufacturing a Gaussian joint significance.

| Timing prescription | Exponent b | Duration multiplier at z=1 | Archived 35-row spectral-aging diagonal sum |
|---|---:|---:|---:|
| Baseline | 1.000 | 2.000 | 26.949 |
| DES central value | 1.003 | 2.004 | 26.967 |
| Brightness-only reinterpretation | 1.481 | 2.792 | 47.058 |

For the high-redshift spectral-aging subset, the same score rises from 3.654 to 24.984 under the brightness-derived exponent. Those predictions use aging rate=(1+z)^(-b) with no refitting to the 35 observations. All rows and original errors are retained in time-only-results.json. The data were already exposed, and the scores are descriptive diagonal comparisons, not a combined independent likelihood or theory-exclusion probability.

Moreover, the equivalent brightness revision already worsened its own 494-row farther-distance transfer score from 424.482 to 432.554 on the unchanged covariance. Renaming epsilon as extra time does not erase that failure or make the earlier brightness fit new evidence.

## Required geometry if timing is retained

More generally use beam distance D_G, whose square sets illuminated area: F=L/[4pi D_G^2 S_E S_t]. To produce the same extra dimming while retaining b_DES, a conditional target is

    D_G/D=(1+z)^[(b_brightness-b_DES)/2].

At z=1 this is 1.18025 in beam distance, or 1.39299 in beam area. The expression is a rearranged requirement, not a derived metric. It must also be applied consistently to calibrators, angular sizes and lensing. The previous inhomogeneous-clock beam calculation offers a framework for computing such an area, but has not supplied a universal field that produces this target. Since the original brightness revision fails farther transfer, this target itself is not an accepted final function.

The consequence is specific: stronger event stretching alone does not reconcile this fixed-distance, fixed-source-calibration branch. A different propagation geometry, actual photon survival or justified source model would need to be specified and evaluated jointly with timing. This does not prove an expanding universe or exclude every hypothetical time mechanism; it identifies an internally linked observational constraint that the next model must respect. All six goals remain open.

Reproduce with `python research_work/results/brightness-distance-consistency/time-only.py`. The result verifies the spectral-data hash and reproduces the archived b=1 score before evaluating the brightness-derived exponent.
