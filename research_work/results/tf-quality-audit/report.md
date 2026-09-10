# Image quality does not resolve the distance-estimator failures

The 73 previously calibrated galaxies were matched to the source's image-quality ratings, dust uncertainties and published inclination-corrected linewidths. All prior predictions remain unchanged. This is a diagnostic on exposed data, not a new holdout or a refitted model.

| Optical image quality | Galaxies | Existing prediction RMS (mag) | Median absolute distance discrepancy |
|---|---:|---:|---:|
| 3 | 1 | 0.341 | 17.0% |
| 4 | 16 | 1.038 | 30.5% |
| 5, highest | 56 | 0.991 | 23.5% |
| At least 4 | 72 | 1.002 | 25.1% |

The quality scores are the catalog's image assessments, not ratings invented from our residuals. Even the highest-quality stratum retains nearly one magnitude of RMS discrepancy. These conditional stratum summaries use the original group-excluded predictions; they are not performance estimates for newly trained quality-cut models.

The median diagonal measurement-only error scale is 0.236 mag after including indicator distance uncertainty, the source's conservative 0.05-mag photometry error, dust uncertainty and linewidth error. Twenty of 73 residuals exceed three of their individual scales. This is not a significance test: intrinsic scatter, fitted-parameter uncertainty, common zero points, covariance and selection are omitted. It establishes that quoting the listed measurement errors alone would leave substantial dispersion unexplained.

Our reconstructed log linewidth differs from the published value by at most 0.00277 dex, equivalent to about 0.021 mag at the existing slope. This compares fixed predictions rather than a complete refit; it is much smaller than the largest residual. Median and maximum differences in reconstructed log-linewidth uncertainty are 0.000243 and 0.00198 dex. The discrepancy is consistent in scale with rounded input reconstruction; its origin is not independently proved.

## A physically informative failure

NGC4424 has image-quality scores 5 in both optical and WISE, yet a -5.265-mag prediction residual. Its reconstructed log linewidth is 1.783087 versus the published 1.783, so this discrepancy is not explained by our inclination-conversion arithmetic.

Independent optical and integral-field observations describe a disturbed stellar disk, noncircular molecular-gas motions, low stellar rotation and substantial random-motion support in NGC4424. This supplies an observational reason that a simple circular-rotation luminosity ruler may be inappropriate. We do not use the paper's modeled distance to replace the indicator measurement. See [Cortes, Kenney and Hardy, 2006](https://arxiv.org/abs/astro-ph/0511081).

A subsequent study explicitly found that disturbed cluster galaxies can have gas linewidths reduced by environmental effects, undermining gas-based Tully-Fisher distances. Its distances depend on a different dynamical/TF construction and are not imported as independent facts. See [Cortes, Kenney and Hardy, 2008](https://arxiv.org/abs/0803.3638).

This literature check was targeted after viewing the largest residual. It neither certifies the remaining galaxies nor licenses deleting this one and calling the new fit validated. NGC4424 and every other failure remain in the archived calibration.

## Consequence for the void test

Our inference is that distance errors may be related to environment, rather than merely random. If surroundings disturb a galaxy's gas, an erroneously inferred distance could correlate with surroundings and imitate the very redshift effect we are looking for. This is a potential confound, not a measured void signal or proof that TF distances are universally unusable.

For clarity, let Dhat=D+deltaD and y=ln(1+z). Under a true distance-only law y=alpha*D, the residual formed with the estimated distance is exactly y-alpha*Dhat=-alpha*deltaD. This is ordinary algebra, not a new physics formula. An environment-dependent mean deltaD therefore produces an environment-dependent redshift residual even with no new void physics. If the same distances locate environment tracers, the environment labels can be wrong as well. Random-error propagation alone cannot fix either bias.

Before fitting a void coefficient, use spatially resolved kinematic/suitability information across the sample, or a distance-estimator mixture that explicitly allows disturbed/non-rotating systems. Such a mixture must retain uncertainty and selection, and its environmental bias needs calibration from indicator anchors rather than target redshift residuals. A known independent starting resource is [Ponomareva et al. 2016's resolved-HI calibrator sample](https://arxiv.org/abs/1609.00378), which contains 32 galaxies with Cepheid/TRGB distances and resolved kinematic information. It is a possible input for the next audit, not already a validated replacement sample. The reported source calibration and distance assumptions must still be checked.

## Reproduction and status

Run `python research_work/results/tf-quality-audit/run.py`. The query retrieves only PGC, QSflag, QWflag, e_Ai, logWmxi and e_logWmxi from J/ApJ/902/145/table1. All 10,737 source identities are unique; all 73 calibration matches and finite measurement scales are verified. Joined values, source hashes and outputs are archived. The initial run exposed a NumPy-integer JSON serialization issue; conversion to a native Boolean count fixed it and the complete rerun succeeded.

Source query: https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJ/902/145/table1&-out=PGC,QSflag,QWflag,e_Ai,logWmxi,e_logWmxi&-out.max=unlimited . Source field conventions and the conservative photometry uncertainty are documented in the previously cached J/ApJ/902/145 ReadMe, notes 2, 4 and 10.

No redshift parameter, transfer mechanism, distance, or holdout status changed. The quality-only shortcut is unsupported by this audit. The next stage-2 action is a sample-wide physical-suitability/error model, with environmental distance bias treated explicitly; the four-stage goal remains incomplete.
