# Calibrators cover the bright nearby sample much more strongly

Applying the original input-quality requirements to all 10,737 staged galaxies leaves 7,375 eligible measurements and 3,362 documented exclusions. All 73 original calibrators are reproduced exactly. We now have explicit per-galaxy coverage flags for deciding where the conditional distance-error calibration is being transported into sparsely checked measurement conditions.

| Corrected apparent i magnitude | Eligible catalog galaxies | Indicator calibrators |
|---|---:|---:|
| Below 12 | 518 | 53 |
| 12 to below 14 | 3,108 | 12 |
| 14 to below 16 | 3,181 | 7 |
| 16 to below 18 | 564 | 1 |
| 18 or fainter | 4 | 0 |

Larger magnitude means fainter apparent light. The calibrator median is 10.90, versus 14.02 for the eligible parent. About 73% of calibrators but 7% of the parent lie below magnitude 12. These counts are not known inclusion probabilities: they describe crossmatches to an existing convenience calibration, not a documented random sampling process.

## Width coverage is broader than joint measurement coverage

Only 23 parent rows lie outside the calibrator range of corrected log linewidth, 323 outside its apparent-magnitude range, and none outside its inclination range. Marginal minimum/maximum ranges therefore give a misleadingly permissive picture of joint coverage.

Of 7,375 eligible rows, 588 (8.0%) lie inside the calibrator convex hull in log linewidth and apparent magnitude, and 452 (6.1%) inside the hull after adding inclination. These are geometric interpolation diagnostics, not validated physical or selection boundaries. All calibrators lie in both hulls, as required.

Crucially, these low percentages do NOT establish that the physical luminosity–rotation relation fails for the other galaxies. A more distant galaxy naturally appears fainter at the same rotation speed even if that relation is universal. Apparent-magnitude hull membership must not become a hard distance-estimator acceptance rule. What is weakly tested is transport of measurement errors, dust corrections and the selected calibration's discrepancy distribution into the much fainter sample. The hull also does not encode sky coverage, morphology, environment, surface brightness or observing strategy.

## Why a larger catalog alone cannot identify voids

The known sampling identity is lambda_observed(r,L,...)=S(r,L,...)*lambda_true(r,L,...), where S is the probability of entering the catalog. This is standard selection/thinning mathematics, not our new physics. With unknown S, a shortage of observed galaxies can mean lower true density, lower detection/selection probability, or both. Many choices of S and true density produce the same catalog.

Observed brightness and linewidth distributions alone do not identify the spatial selection function. Nor can the ratio of 73 calibrators to 7,375 parent rows supply inverse-probability weights without a sampling model. Our earlier core/broad residual mixture likewise does not identify the missing-galaxy distribution. Assigning wider independent distance errors cannot resolve unknown systematic incompleteness.

For stage 2, retain directly estimated distances where available and treat the broader TF catalog as conditional tracer information. Before calling a catalog gap a void, obtain the source survey footprints and measurement/target-selection rules, model observability of an explicitly chosen tracer population, and test that selection model against counts and measurement distributions. Calibration and targeting histories must be represented separately. A restricted region with a documented selection function may support a defensible test even if the entire combined catalog does not. The present hull is not such a region and must not substitute for that work.

## Reproduction and verification

Run `python research_work/results/tf-calibration-support/run.py`. `support.csv` records all 7,375 eligible identities and flags; `excluded.csv` records the remaining 3,362 identities and reasons. Every one of 10,737 source identities has exactly one disposition. All 73 calibration values of x and corrected magnitude agree with the prior calculation to 1e-12; both hulls contain every calibrator. All eligible rows are accounted for by the declared magnitude bins. Aggregate results and input hashes are in results.json.

No redshift, predicted distance or observed residual was used in assigning support. No new fit, weight, exclusion rule or physical void coefficient was adopted. The previous turn answered the user's arrival-time question and did not advance this selection calculation; this turn supplies the missing parent-support evidence. The redshift calculator and conditional messenger-delay conclusions remain unchanged. All four goal stages remain incomplete.
