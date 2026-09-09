# Focused temporal-history test

9 September 2026. Exploratory analysis of previously inspected Pantheon+ observations; not independent validation.

The most promising one-parameter change is an evolving propagation history in flat, fixed material geometry. With the original nearby coefficient and Cepheid calibration retained, the best fit is p = 0.263906. It improves the covariance-weighted brightness fit more than either negative curvature alone or a constant magnitude offset, while substantially reducing the trend in residuals with redshift. This supports prioritizing the history branch for further testing; it does not demonstrate atomic protection or a sustainable physical field.

## Model and data

Assume material atomic frequencies are fixed, propagation is homogeneous and nondispersive, photons are conserved, and physical space is static and flat. The proposed rule is dn/dt = gamma n^p, with n(today) = 1 and kappa = gamma/c = 0.000077315 per million light-years. Measured redshift and neighboring signal arrival intervals both scale by n(observer)/n(emitter). Integration gives R(z) = [(1+z)^p - 1]/(kappa p), with the p=0 limit R=ln(1+z)/kappa. Luminosity distance is (1+z)R and angular diameter distance is R under these assumptions. Atomic protection and shared photon/tensor propagation remain assumptions, not newly derived interactions.

The fit uses 960 noncalibrator supernovae with zHD >= 0.1, the released full STAT+SYS covariance, and the previous Cepheid-host calibration M=-19.24888454. Calibration uncertainty and cross-covariance are propagated as before. In the implemented observer-frame convention, DL=(1+zHEL)R(zHD). This convention and the published standardization retain processing assumptions from the original data analysis. Kappa is fixed; uncertainty in that fitted coefficient is not propagated here.

## Results

| Candidate | Added fitted parameters | Chi-square | Change in AIC relative to original |
|---|---:|---:|---:|
| Original exponential, flat | 0 | 898.475 | 0 |
| Constant magnitude offset only | 1 | 871.552 | -24.922 |
| Negative curvature only | 1 | 873.041 | -23.433 |
| History only, p=0.263906 | 1 | 850.201 | -46.274 |
| History plus offset, diagnostic | 2 | 840.647 | -53.828 |

Lower chi-square and AIC favor a fit within the stated common likelihood. AIC values are descriptive comparisons, not probabilities of a physical theory. The models with one added parameter can be compared without a difference in parameter-count penalty. The previously evaluated fixed expanding-universe reference has chi-square 838.483, but was not optimized; the new history branch has not outperformed that reference.

The conditional one-parameter profile interval for history-only p is 0.226425–0.301224 at Delta chi-square=1, and 0.190289–0.336896 at Delta chi-square=3.84146. These intervals condition on fixed kappa, model assumptions and the released covariance; they are not complete physical uncertainty budgets.

The four correlated redshift-bin mean residuals change from [0.1383, 0.2004, 0.2182, 0.3631] mag to [0.0862, 0.1081, 0.0628, 0.1144] mag. The last-minus-first contrast falls from 0.2249 to 0.0283 mag. The constant-bin-residual diagnostic falls from 36.743 to 6.465. No naive post-fit significance is assigned to that diagnostic. A constant offset alone leaves the trend unchanged; curvature-only gives a diagnostic of 12.523.

The remaining positive residual means the history-only model still predicts sources systematically brighter than their calibrated observations. Allowing a separate constant offset yields p=0.224240 and an offset of +0.101119 mag, with chi-square 840.647. This two-parameter diagnostic indicates that the history preference survives a free zero point. It does not establish that the calibration is wrong, nor identify the origin of the offset. No redshift-dependent luminosity evolution was modeled.

## Behavior and sensitivity

At redshift 1, the history-only fit gives R=3.01622 Gpc and DL=6.03244 Gpc in the ideal common rest frame. This is 9.73% farther than the original flat exponential relation. At redshift 2, R=5.05397 Gpc and DL=15.16190 Gpc, a 16.01% increase in distance. The arrival-duration factor remains exactly 1+z for neighboring events; changing p does not introduce a separate duration-stretch parameter.

Separate fits below z=0.6 (831 objects) and at or above z=0.6 (129 objects) give p=0.33211 and p=0.33000 respectively. These subsets share calibration and systematic correlations; they are not independent holdouts. The joint optimum can differ from both subset optima because it includes cross-subset covariance. Excluding z>=1 gives p=0.27128. Thus the preference for positive p is not supplied only by the 25 highest-redshift objects. The subset checks are descriptive and were performed after inspecting the full sample.

The analytic distance integral was checked against direct numerical quadrature for four p values and three redshifts; the maximum absolute discrepancy in the dimensionless integral was 4.44e-16. A 401-point scan over the specified parameter ranges checked the scalar optimizer. These verify the calculation, not the model's physical validity.

## Recommendation and limits

Prioritize the flat evolving-history model as the smallest promising next candidate. Retain the history-plus-offset result as a calibration sensitivity diagnostic, and keep curvature-only as a documented alternative. There is no present reason to add CMB thermalization, rotation or gravity stabilizers to this brightness comparison.

The next observational milestone is an independent brightness validation with the selected formula and fitting protocol frozen, followed by a carefully sourced angular-distance check. The relation DL=(1+z)DA remains fixed in this model for every p; changing the history cannot rescue a failure of that relation. Published angular distances often inherit cosmological or source-physics assumptions, which must be assessed before treating them as independent measurements.

Stable atomic ratios, a microscopic interaction, and a field trajectory that persists with energy conservation and acceptable stability remain unresolved. This test identifies a better phenomenological history, not a finished nonexpanding cosmology. No new angular-distance, duration, gravitational-wave or CMB dataset was fitted.

## Reproduction

Extract focused_history_test.zip and run `python focused_history_test/run_test.py` with NumPy, SciPy and pandas installed. Included data and calibration code reproduce this analysis offline. Importing the original brightness module regenerates its original baseline result files in shared_interaction_test. The input supernova data are the Pantheon+ release used previously, repository commit c447f0fea703fcd0fff57de5000947b5ca81286b: https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b . The archive manifest records hashes of included files.
