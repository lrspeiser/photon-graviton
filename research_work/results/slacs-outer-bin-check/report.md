# Predicting outer stellar motions from inner bins

## Outcome

We withheld the outermost radial bin of each of the same six training galaxies from the parameter fit. The extra-force branch has a smaller absolute covariance-adjusted outer residual in five of six systems, although the difference in J1621+3931 is negligible. It performs worse in J1204+0358. Several outer-bin mismatches remain substantial even after the extra force is included.

These are already-exposed observations. This is an exploratory leave-one-bin-out robustness diagnostic, not the untouched prediction required by goal 6, and not a completed joint gravity test.

| Galaxy | Observed outer Vrms | Baseline raw prediction | Extra-force raw prediction | Baseline conditional residual | Extra-force conditional residual |
|---|---:|---:|---:|---:|---:|
| J0037-0942 | 242.91 | 231.25 | 239.34 | +3.91 | +1.16 |
| J1112+0826 | 255.30 | 232.53 | 237.68 | +5.62 | +4.53 |
| J1204+0358 | 219.37 | 218.96 | 222.75 | +0.20 | -1.43 |
| J1402+6321 | 255.57 | 234.85 | 242.00 | +5.63 | +3.96 |
| J1621+3931 | 225.72 | 219.71 | 226.47 | +0.67 | -0.66 |
| J1630+4520 | 240.80 | 222.77 | 229.52 | +4.83 | +3.06 |

The first three numerical columns are km/s. The last two are residuals divided by the **conditional measurement-error** standard deviation. They are not calibrated significance levels: fitted-parameter uncertainty and model/PSF/light-profile uncertainty are omitted. A positive residual means the measured outer value exceeds its conditional predicted mean. The raw physical predictions are shown separately because correlated-error conditioning shifts that mean.

## What was fitted

Only the inner N-1 bins constrain each model's total mass and constant beta, with the same mass bounds, widened beta range [-2,0.45], fixed published stellar profiles, force coefficients and physical cutoff as the preceding run. Starting mass is the common 10^11 solar masses, rather than an earlier fitted mass. Both models use the same starts. All fitted betas are interior to the declared range.

The held-out outer measurement does not enter the objective. Its covariance with the retained bins does enter the subsequent predictive calculation, as required for correlated Gaussian measurement noise. No lensing measurement selects parameters. The same fit also produces an archived lens-angle prediction, but this report focuses on the outer-bin diagnostic.

## Known Gaussian conditioning, not a new physical formula

Partition the measurements into retained inner bins t and omitted outer bin h. At a fitted physical prediction mu, the conditional measurement model is

    E[y_h | y_t, parameters] = mu_h + C_ht C_tt^(-1) (y_t-mu_t)
    Var[y_h | y_t, parameters] = C_hh - C_ht C_tt^(-1) C_th.

The diagnostic residual divides y_h minus that conditional mean by the square root of that variance. A Cholesky solve fits the inner data with the full retained covariance. The exact fixed-parameter Gaussian identity

    chi2_full = chi2_inner + residual_conditional^2

is verified to below 1e-7 for every fit. This checks that the correlated errors are handled consistently; it does not establish predictive uncertainty after estimating parameters. No arbitrary diagonalization, extra error inflation or independent treatment of correlated bins is used.

## Numerical and scientific limits

The radial/angle/deprojection integration is doubled at every fitted solution. Fractional predicted Vrms changes remain below 4.18e-5, much smaller than the reported measurement differences. All optimizer starts, predictions, conditional means/errors, factorization checks and refinement changes are retained.

The predictor uses estimated mass and beta as fixed values. A proper predictive distribution must average over their uncertainty from the inner data, with declared priors or a calibrated frequentist procedure. Therefore counting these standardized residuals as rejection probabilities would be premature. The narrow constant-beta spherical source model, fixed seeing and image profiles remain additional assumptions.

The six systems are a small, selected and previously examined training overlap. Cross-object measurement covariance is not supplied here. The empirical extra-force relation still lacks a source/capture derivation and independent stellar-population mass calibration. These limitations apply even where the extra-force prediction is closer.

## Next step

Propagate the inner-fit mass/orbit uncertainty into the outer prediction before assigning coverage or probabilities. Keep the same held-out-bin rule and retain every system; do not remove difficult galaxies or tune an outer-only force. Independent validation and the full six-goal program remain open.

Run `python research_work/results/slacs-outer-bin-check/run.py` to reproduce. Previous flagged and missing-profile exclusions remain explicitly listed.
