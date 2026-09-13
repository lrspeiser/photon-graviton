# Source starting-phase sensitivity of the SN 2006mk light-curve pilot

This fit tests whether fixing the first spectral phase at -3.4 days artificially produced the previous preference for stretched light curves. It keeps the same 25 unanchored measurements, empirical source, legacy filters and one measured brightness anchor per band. It fits one shared starting phase to both bands for each fixed event-time factor.

## Results

| Time map | Starting phase before fit | Fitted starting phase | Fitted chi-square | Fitted chi-square + log determinant |
| --- | ---: | ---: | ---: | ---: |
| spectral_pair | -3.40 | -2.890 | 36.769 | 41.587 |
| event_stretch_1_plus_z | -3.40 | -3.258 | 36.649 | 41.386 |
| unchanged_arrival_intervals | -3.40 | -4.556 | 200.054 | 204.358 |

The spectral-derived stretch and 1+z stretch become almost indistinguishable by this conditional score (difference about 0.20). Unchanged timing remains substantially worse under the same fixed empirical source. Freeing the origin of the source-time axis does not repair its faster fading. These are not calibrated significance statements or confidence intervals.

The fitted phase bounds are -10 to +5 days. All solutions lie inside them. Every previously included photometric point stays within the template support throughout this parameter box, so no model-dependent data selection changes the score. This is exploratory fitting on exposed observations; the previously unadjusted source predictions remain separately recorded.

## Likelihood and formula provenance

Use phase=q+(t_observer-t_first_spectrum)/A, where q is now fitted and A remains fixed in each case. This is a known affine time parameterization, not a companion mechanism. No independent prior on q is asserted; its deviation from a published or independently estimated spectral age is not assigned a significance here.

The anchor covariance depends on the source ratios. Accordingly the minimized Gaussian residual score is chi-square + log(det C), with the common Gaussian constant omitted. The earlier report listed chi-square only as a descriptive diagnostic. Including the determinant avoids interpreting a model with larger propagated anchor variance as equally predictive without accounting for that variance. The original independent-photometric-error assumption remains conditional.

One observed anchor in each band means this is still a temporal-shape comparison, not absolute brightness or a prediction of cross-band luminosity normalization. No event width, dust law, source-population variation or propagation coefficient is fitted here. In particular, a broader intrinsic supernova history can imitate part of an event stretch; this starting-phase check does not remove that ambiguity.

## Numerical verification

Integrating the source once at each tabulated phase and then interpolating is equivalent to integrating its linearly interpolated spectrum. An off-grid check agrees to relative tolerance 1e-12. At fixed q=-3.4, all three summed chi-square scores reproduce the previous direct-integration calculation to 1e-8. The rank-one covariance formula is used unchanged from that checked calculation.

Profiles on 301 and 601 starting-phase nodes, each followed by a local bounded refinement near its best node, agree in best score to below 1e-6. All 601 profile values per hypothesis are retained, including inferior values and boundaries. A repeated complete run reproduces the result file exactly. These verify this numerical comparison; they do not calibrate the omitted observational or source systematics.

## Research consequence

The current observations remain compatible with a stretched empirical source history and are poorly described by the tested unchanged-time version of that same source. Neither result derives the stretching from photon-companion conversion. A joint physical theory still needs a source-variation model, calibrated passbands, uncertainty shared with spectral ages, independent distance/redshift tests and absolute energy accounting. All six original objectives remain open.

[Previous full in-support curve comparison](../sn2006mk-lightcurve-pilot/report.md) | [Source assumptions](../hsiao-source-audit/report.md)
