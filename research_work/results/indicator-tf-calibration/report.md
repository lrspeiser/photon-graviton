# A working distance estimator, with inadequate precision for small voids

The staged brightness/rotation measurements now have an actual indicator-anchored calibration and group-excluded distance predictions. This moves the distance-versus-environment test beyond source inventory. It does not yet produce a reliable physical void map or validate the photon-companion mechanism.

Of 109 matched galaxies, 73 pass the declared input requirements, spanning 3.12–36.81 Mpc in 57 catalog groups. There are 46 TRGB, 17 Cepheid, 9 SBF and 1 maser calibrators. Thirty-five lack required measurements and one fails the declared inclination/width requirements. All exclusions are recorded, and no residual-based rejection was used. No recession velocity or redshift outcome is used in the calibration.

## Equations and provenance

The known empirical Tully–Fisher relation is M_i=a*x+b, where x=log10[Wmx/sin(i)]-2.5. This is not a new formula or a derivation of gravity. The source convention gives corrected apparent magnitude m_i=icmag-Ai. The [source table documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJ/902/145?format=html&tex=true) states that icmag already includes foreground, K and aperture corrections, but excludes internal inclination-dependent dust attenuation. We subtract Ai once.

The added conditional photometry is m_i=M_i+5log10(D)+25+5*p*alpha*D/ln(10). The previously fitted alpha is fixed at 0.0002488993286382367/Mpc. p=0 is an inverse-square diagnostic; p=1 assumes the proposed shared spectral/event stretching, photon conservation and static geometric dilution. The extra term is a consequence of those postulates, not an established new interaction. The source's existing rest-band corrections are retained provisionally. Published indicator distances are stipulated here; their own brightness response and common calibration are not jointly refitted.

## Results

Each eligible galaxy is predicted using a fit that excludes every calibrator with its CF4 group identifier. These are exploratory cross-validation predictions, not genuinely withheld observations.

| Quantity | Inverse-square control p=0 | Shared stretch p=1 |
|---|---:|---:|
| Fitted slope a | -7.5679 | -7.5736 |
| Fitted intercept b | -20.7665 | -20.7744 |
| Training RMS, mag | 0.9472 | 0.9486 |
| Group-excluded RMS, mag | 0.9943 | 0.9958 |
| Median absolute fractional distance discrepancy | 24.89% | 24.95% |
| Distance prediction/reference ratio, 16th percentile | 0.7402 | 0.7409 |
| Ratio, median | 1.0051 | 1.0055 |
| Ratio, 84th percentile | 1.3538 | 1.3563 |

These percentiles describe this sample, not universal confidence intervals. A roughly 25% discrepancy corresponds to 5 Mpc at 20 Mpc or 10 Mpc at 40 Mpc as a scale illustration; 40 Mpc is already outside this calibration range. Such uncertainty cannot be ignored when assigning 1–5 Mpc environment boundaries. A probabilistic broad environment estimate might still be possible, but hard void labels would overstate the information.

The largest p=1 error is NGC4424: a -5.265-mag residual and predicted/reference distance ratio 0.0888. The staged row has Wmx=55 km/s, inclination 65 degrees, icmag=10.92 and Ai=0.21. It is retained. This is a failure of the present estimator for that object, not grounds to discard its measurement after inspecting the residual. UGC7599 and KK149 have ratios 2.57 and 2.43. A single RMS or Gaussian scatter would hide these tails. Physical suitability, inclination/systematic linewidth errors and source photometry quality need independent assessment before a revised sample rule is adopted.

For p=1, group-excluded RMS is 0.635 mag for Cepheids, 0.819 for SBF, and 1.137 for TRGB. The single maser point is not a precision test. Differences may reflect population, selection, kinematics or calibration; this calculation does not distinguish them. Predictor-error propagation alone gives median magnitude uncertainty 0.186 mag and is insufficient to characterize the observed spread.

## Verification and limits

All 73 galaxies receive exactly one group-excluded prediction per scenario; every regression has full rank. Independent bracketed roots and Lambert W distance inversions agree to 1.50e-14 relative error. Both scenarios and all 146 predictions are archived in cross-validation.csv; source hashes and results are in results.json. Reproduce with `python research_work/results/indicator-tf-calibration/run.py`.

The ordinary least-squares fit is deliberately transparent and is not a full errors-in-variables or selection likelihood. Slope attenuation, common indicator zero points, dust errors, source selection, inherited redshift-dependent corrections and correlated group membership remain. The CF4 group identifiers may themselves incorporate velocities; only those identifiers are used to keep groups out of the training fold. They never supply distances or fitted predictors. Shared zero-point correlations still cross folds. Source SDSS/WISE quality flags exist but were not included in the previously staged feature subset.

The immediate next step is an input-quality and physical-suitability audit, followed by a declared uncertainty/selection model and renewed exploratory validation. Any revised cut must be based on independent measurement quality or physical applicability, with the current failures preserved. Only then should broader distance probabilities be propagated into environment inference. Do not retune alpha or select a photometry scenario from this tiny difference. The root redshift calculator is unchanged, the physical cause remains unfinished, and all four goal stages remain incomplete.
