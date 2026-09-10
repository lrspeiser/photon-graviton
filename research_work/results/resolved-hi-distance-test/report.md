# Resolved rotation does not yet deliver a substantially better distance ruler

We tested whether spatially resolved rotation measurements improve the indicator-anchored distance estimator, using identical galaxies and fitting rules for each velocity statistic. The result is a modest difference in median discrepancy and essentially unchanged overall scatter. This is useful negative evidence against immediately replacing the existing linewidth estimator and declaring the void-map problem solved.

## Matched observations

The [Ponomareva et al. 2016 atlas](https://arxiv.org/abs/1609.00378) supplies 32 nearby galaxies with resolved HI measurements. Fifteen match the previous 73 eligible indicator calibrators by normalized NGC/IC identity. NGC4605 lacks a measured flat rotation speed, leaving 14 objects in 14 CF4 groups for the paired comparison, spanning 3.12–26.55 Mpc. All 32 atlas identities and match decisions are archived; the 59 previous calibrators outside the common subset remain outside this comparison, not rejected as bad galaxies.

Only kinematic speeds, quoted errors, inclination and profile classifications were extracted for the test. Atlas systemic velocities, distances and halo quantities do not enter the calculation. Although the source article discusses dark-matter interpretations, this comparison uses empirical rotation measurements without adopting a halo model. Kinematic deprojection and disk assumptions remain measurement-model dependencies.

## Same-sample comparison

All three fits use the same corrected i-band magnitudes, stipulated indicator distances, fixed prior alpha and p=1 conditional photometry. The known empirical formula is M_i=a*x+b. Predictors are x=log10(Wmx/sin i)-2.5, log10(2Vmax)-2.5, or log10(2Vflat)-2.5. Neither the TF relation nor the logarithmic substitutions are new physics. Each prediction excludes its galaxy's entire group. There is no clipping, no outcome-selected photometry scenario and no alpha refit.

| Predictor | Group-excluded magnitude RMS | Median absolute distance discrepancy |
|---|---:|---:|
| Integrated linewidth | 0.4754 mag | 15.16% |
| Maximum resolved rotation speed | 0.4880 mag | 14.30% |
| Outer flat rotation speed | 0.4779 mag | 13.92% |

The outer-speed fit's predicted/reference distance ratios at the 16th, 50th and 84th percentiles are 0.830, 0.974 and 1.258. With only 14 exposed objects these are sample summaries, not calibrated population confidence limits. Median discrepancy favors Vflat slightly while RMS does not; this is not evidence of a statistically established improvement.

The prior fit trained on the larger 73-object sample already had only 0.4182-mag RMS on this overlap under its original group exclusions, compared with 0.9958 mag across all 73. Thus the lower scatter of this overlap predates replacement of the velocity statistic. Coverage, sample composition and training size differ; no causal contribution can be assigned from that comparison alone. Reporting a reduction from 25% to 14% as the benefit of resolved rotation would be misleading.

## What this changes for the void test

Resolved observations remain useful for diagnosing disturbed kinematics, but the available overlap is small and contains only 19.2% of our 73 calibration objects. The new measurements do not certify undisturbed or circular motion; asymmetric profiles remain included. They also do not create a dense environment map along the longer galaxy sightlines.

We should carry this paired result into a joint uncertainty/selection model instead of choosing the best-looking statistic. The next useful comparison needs more matched resolved measurements or a physically specified mixture of suitable and disturbed systems anchored by independently estimated distances. Distance-error dependence on environment must remain explicit. Neither a quality-image cut nor substituting Vflat has demonstrated that it removes that confounding.

The 14% discrepancy is still a few Mpc at a 20-Mpc distance as a scale illustration, before common calibration and selection effects. Hard 1-Mpc void boundaries would continue to overstate the precision. Coarser probabilistic environments may be testable, but need a selection model and cannot treat catalog gaps as empty space.

## Verification and limitations

The source PDF's Table 5 on page 10 was visually inspected against the HTML extraction, including NGC2841 (Vmax=325, Vflat=290 km/s), NGC4258 (242, 200), NGC4414 (237, 185), and NGC4605's absent Vflat. Extraction retains all 32 unique identities and missing values; all 42 paired predictions have exactly one group exclusion and full-rank fits. Bracketed distance roots and Lambert W inversion agree to 3.59e-16 relative error. The raw-source and calibration-input hashes are recorded.

Reproduce with `python research_work/results/resolved-hi-distance-test/run.py`. Inputs are the cached source HTML and previously committed calibration predictions. Coverage, individual predictions and results are committed beside the script.

Quoted velocity errors are archived but not used in this initial unweighted regression. Errors in predictors, distance and dust, inclination systematics, common zero points, sample selection, and source brightness-correction dependencies remain. The same exposed data cannot be relabeled as a fresh holdout. No complete physical candidate, environment law or full-model freeze was established. The root redshift equation remains unchanged and all four goal stages remain incomplete.
