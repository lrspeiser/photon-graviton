# Quick residual check: a radial-shape problem, not a uniform gravity shortage

13 September 2026. No fits or downloads. All samples previously exposed.

## Main result

The exact-third model is almost unbiased when averaged across entire galaxies
(-0.51 km/s), but that average hides underprediction inside and overprediction
outside. Increasing the shared amplitude cannot fix both signs. This diagnoses
the stipulated profile plus fixed ordinary-matter baseline; it does not prove
that the one-third exponent or companion physics is the cause.

Negative bias means predicted speed is too low. Each galaxy has equal weight
within a bin; the same galaxy may appear in several radial bins.

| Radius / disk scale length | Galaxies / points | Mean error km/s | Descriptive 95% interval | Speed RMSE km/s |
|---|---:|---:|---:|---:|
| Below 1 | 137 / 744 | -6.73 | -10.72 to -2.96 | 25.78 |
| 1 to below 3 | 149 / 1122 | -7.49 | -11.06 to -3.88 | 25.92 |
| At least 3 | 134 / 1284 | +10.41 | +5.39 to +15.43 | 33.50 |

The intervals resample galaxies, not their full measurement/calibration errors.
They are conditional descriptive intervals with no multiple-comparison correction.
Distances, inclinations and stellar mass-to-light ratios remain fixed. The input
distances remain scenario facts, including those originally inferred with an
expansion-based calibration; this is not an independent distance-law test.

The mean sign pattern persists in train (-8.30, -7.63, +8.88 km/s), validation
(-3.50, -6.56, +11.88) and test (-5.44, -7.94, +13.89). These are the same old
partitions, not independent new validation. Among 122 galaxies with both inner
and outer coverage, the average outer-minus-inner fractional residual is +12.79
percentage points; 80/122 have a positive difference. Paired comparison reduces
the concern that different galaxies alone explain the radial pattern, but does
not establish a universal trend or causal mechanism.

## Matched controls by radius

Entries are mean signed error / speed RMSE, in km/s, on the same radial samples.
NFW is the existing restricted shared scaling, not a general dark-matter test.

| Model | Inner | Middle | Outer |
|---|---:|---:|---:|
| Ordinary matter | -12.55 / 28.76 | -36.11 / 45.16 | -57.81 / 66.56 |
| Exact-third companion | -6.73 / 25.78 | -7.49 / 25.92 | +10.41 / 33.50 |
| Fitted simple MOND | +2.57 / 23.49 | -4.83 / 21.66 | -5.45 / 20.47 |
| Restricted shared NFW | +1.75 / 23.37 | -3.15 / 30.91 | +8.88 / 47.34 |

MOND has lower speed RMSE in all three aggregate radial bins. These descriptive
scores do not replace the original equal-galaxy logarithmic fitting objective.

## Surface brightness, gas proxy and concentration of errors

| Subsample | Galaxies | Companion mean error km/s | Speed RMSE km/s |
|---|---:|---:|---:|
| Central disk brightness below 100 | 34 | +2.96 | 23.28 |
| Brightness 100 to below 500 | 44 | -0.09 | 21.74 |
| Brightness at least 500 | 71 | -2.43 | 34.32 |
| Atomic-gas proxy below 0.2 | 43 | -5.42 | 35.80 |
| Gas proxy 0.2 to below 0.5 | 42 | +5.80 | 30.52 |
| Gas proxy at least 0.5 | 64 | -1.35 | 21.15 |

Brightness units are solar luminosities/pc squared. The gas proxy is
1.33 MHI/(1.33 MHI+0.5 L3.6), not a complete gas fraction: it omits molecular gas
and a distinct bulge mass-to-light ratio. Higher km/s errors in bright systems
partly reflect speed scale; log errors and intervals are retained in the JSON.
Brightness-bin signed-error intervals all cross zero. Gas-proxy biases are not
monotonic. These marginal summaries therefore do not justify changing the
one-third intensity exponent on their own, or separate environmental corrections.
The variables are correlated and no multivariate causal fit was performed.

The worst five galaxies contribute 24.38% of the equal-galaxy squared-speed loss;
the worst ten contribute 40.78%. Using the actual logarithmic objective, the
corresponding shares are 23.14% and 35.93%, with different worst objects. The top
squared-speed contributors include NGC5985, NGC0801 and NGC0289; the top logarithmic
contributors are NGC3741, UGC06667 and NGC2915. Errors are concentrated but not
explained by a single bad target, and both over- and underpredictions matter.
No object was removed.

## Six lens galaxies: the benefit remains concentrated

Saved full-covariance all-motion scores, Chabrier proxy, with the same four fitted
stellar freedoms and lens-calibrated mass. These are fitted observations, not
predictions. Lower chi-squared is better.

| Galaxy | Transferred companion | Matched stellar-only | Difference |
|---|---:|---:|---:|
| J0037-0942 | 10.44 | 32.42 | -21.97 |
| J1112+0826 | 18.67 | 18.03 | +0.64 |
| J1204+0358 | 14.33 | 13.95 | +0.39 |
| J1402+6321 | 64.97 | 105.96 | -40.98 |
| J1621+3931 | 3.41 | 3.09 | +0.31 |
| J1630+4520 | 1.74 | 1.65 | +0.09 |

The Salpeter proxy gives the same improvement/worsening classification and is
recorded separately in the JSON. Only two systems improve; J1402 still contributes
about 57% of the companion chi-squared. Saved signed velocity residuals are
negative on average in all six, but a low average stellar speed is not directly
a missing-mass estimate: lens normalization and orbital anisotropy are coupled.
The original full covariance remains represented by its saved likelihood score;
unweighted signed means are supplementary diagnostics, not replacement scores.
No radial rebinning of lens data or new covariance reconstruction was performed.

## What to try next, and what not to infer

The lowest-cost useful follow-up is a single shared size/retention modification
that can increase inner acceleration while decreasing it at the sampled outer
radii. Merely increasing the amplitude is inconsistent with the average residual
signs. Merely moving all mass inward while preserving total mass is also
insufficient in the spherical deposited-component approximation: inward transfer
can only increase or preserve enclosed mass, so it cannot lower outer circular
speed. This is ordinary enclosed-mass gravity, not new companion physics.

One diagnostic candidate is the one-parameter homology

    rho_s(r) = s^(-2) rho_reference(r/s),  0 < s <= 1
    M_s(<r) = s M_reference(<r/s)

For a finite-density center, enclosed mass grows approximately as r cubed, so
inner acceleration increases by approximately 1/s squared; beyond most stored
mass it decreases by approximately s. Intermediate radii, including hollow
profiles, require actual evaluation. The total stored mass becomes s times the
reference total, so this is NOT mass-conserving inward migration: the missing
stored energy must remain traveling, escape or be transferred explicitly. These
identities are standard change-of-variable mathematics; the tied size/retention
rule is an untested project hypothesis, not a derived capture law or an established
novel physical law. At s=1 it returns the unchanged exact-third reference.

This candidate has not been fitted or adopted. If explored next, fit only its
single shared parameter on the existing training galaxies, retain the original
controls and all comparison targets, and separately check lens effects. Do not
claim success from restoring the average bias alone or rename reduced inventory
as energy conservation without its receiving sector. The current quick check
supports investigating radial shape; it does not establish this particular fix.

## Reproduction and checks

Run `python research_work/results/companion-extensions/shared-residual-map.py`.
The protocol fixes bins before this diagnostic calculation. The executable
reproduces all archived split RMSE and logarithmic scores to 1e-12 relative
tolerance, checks 149 galaxies/3150 radii and matches the lens-control reference
scores and observed vectors. It records source hashes, all bins/splits, per-galaxy
losses and both population-proxy lens comparisons. Bootstrap seed is 13092026.
The calculation uses saved predictions and introduces zero fitted parameters.
