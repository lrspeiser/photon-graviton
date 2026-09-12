# Joint detector-flux likelihood with shared calibration uncertainty

The forward calculator now has a complementary normalized likelihood for comparing predicted native FLUXCAL to measurement arrays. The likelihood retains correlations induced by a common survey/filter zero point. It is a checked statistical component, not an observational fit or a completed joint-light theory.

## Statistical assumptions and known formulas

Let y be measured native flux, m the model prediction, D the diagonal measurement-error covariance, and delta the vector of residual band zero-point offsets. For small offsets, the native-system response is linearized as

\[
y=m+J\delta+\epsilon,\qquad
J_{ik}=-0.4\ln(10)m_i1_{\operatorname{band}(i)=k}.
\]

Assume delta is Gaussian with covariance V and independent Gaussian measurement noise epsilon with covariance D. These are explicit statistical assumptions, not new physical laws. Marginalizing delta gives the known Gaussian relation

\[
C=D+JVJ^T,\qquad
\log L=-\tfrac12\left[N\log(2\pi)+\log\det C+(y-m)^TC^{-1}(y-m)\right].
\]

The model determines J, so C changes with the prediction. The determinant term must remain in the likelihood. Dropping it would compare probability densities with inconsistent normalization.

The implementation uses the known matrix determinant identity and Woodbury inverse with a small Cholesky solve over the calibration dimensions. This avoids forming an observation-by-observation dense covariance for normal use. The independent verification deliberately forms that full covariance instead.

## How to use this component

Concatenate observations from all events that share the calibration before calling the likelihood. Multiplying separate event likelihoods with separately marginalized copies of the same offset would incorrectly allow different calibration errors for each event. For multi-survey samples, retain the relevant cross-survey covariance and expand band indices accordingly.

The code accepts the covariance explicitly; it does not silently assume independent bands. The tests load the measured DES5YR zero-point subblock. Predicted fluxes come from the combined transport calculator with declared artificial sources; noisy test values are newly generated artificial observations. Original observed fluxes remain untouched.

## Independent checks

Nine cases cover 4, 12 and 40 observations and zero, nominal and four-times-nominal zero-point covariance. The small-matrix implementation agrees with an independent dense multivariate-normal calculation to maximum absolute log-likelihood difference 5.2e-13.

A separate three-measurement, one-band calculation integrates over the shared offset numerically instead of using a covariance formula. It agrees with the analytic marginal likelihood within numerical precision. This verifies the linearized Gaussian model, not the exact nonlinear magnitude-offset model.

A constructed example with independent 2% measurement errors and the measured g-band zero-point component gives:

| Number of measurements | Correct uncertainty on average | If calibration is incorrectly independent |
| ---: | ---: | ---: |
| 1 | 2.0762% | 2.0762% |
| 10 | 0.8429% | 0.6565% |
| 100 | 0.5921% | 0.2076% |

The shared contribution remains when averaging. These figures illustrate uncertainty propagation; they are not measured precision for the DES sample.

## Remaining scope

The likelihood currently assumes diagonal measurement noise before adding calibration covariance. Additional observational covariance, source-population variation, uncertain distances, source-time offsets, dust, passband shifts and selection require explicit treatment. The approximation linearizes zero-point offsets; its accuracy relative to full nonlinear offset marginalization remains unverified.

This component can score the full predicted flux-time array, which is necessary for joint timing and brightness inference. It does not independently determine intrinsic luminosity or duration evolution, and cannot remove the exact source/transport ambiguity identified earlier. It also does not contain a redshift measurement likelihood or a physical receiver-sector model. All six scientific demonstrations remain open.

```powershell
python research_work/results/joint-light-likelihood/verify.py
```
