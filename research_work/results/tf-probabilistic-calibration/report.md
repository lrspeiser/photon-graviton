# Probabilistic TF calibration: better tail score, unstable generalization

We now have an executable conditional error model for the 73 indicator-anchored galaxies, including reported measurement uncertainty and an allowance for unusually large discrepancies. It retains every galaxy. It does not yet give validated distance probabilities for the environment map.

## Statistical assumptions, not new physics

The known empirical relation remains y=a*x+b+epsilon, with x=log10(Wmx/sin i)-2.5. We use the published corrected log linewidth rather than reconstructing it from rounded values. y is the previously derived corrected absolute i magnitude under the p=1 shared-stretch photometry, with the redshift coefficient alpha fixed. The model neither estimates alpha nor derives a photon interaction.

The Gaussian baseline has residual variance v+s^2. The alternative uses the known Gaussian-mixture density (1-f)*N(0,v+s^2)+f*N(0,v+25s^2), where v=sigma_mu^2+0.05^2+e_Ai^2+a^2*sigma_x^2. The broad/core intrinsic-width ratio five and f<=0.5 are declared analyst choices, not measured properties of disturbed galaxies. A broad-component probability is not a physical classification.

The measurement-error convolution assumes independent Gaussian errors and a locally flat prior for true x. It therefore does not solve the population selection or full errors-in-variables problem. Common zero points, parameter uncertainty, environment-dependent bias and photometry/distance calibration dependencies remain absent.

## Group-excluded predictions

For each of 57 groups, all its calibrators are omitted from fitting. Both models then predict every omitted galaxy's noisy absolute-magnitude measurement. All 73 galaxies receive one prediction per model. These are exposed exploratory validation cases, not fresh holdouts.

| Diagnostic | Gaussian | Core plus broad component |
|---|---:|---:|
| Full-sample intrinsic core scale | 0.860 mag | 0.621 mag |
| Full-sample broad fraction | 0 | 0.0338 |
| Sum of group-excluded log densities | -112.069 | -102.243 |
| Measurements inside nominal central 68% interval | 57/73 | 47/73 |
| Measurements inside nominal central 95% interval | 69/73 | 67/73 |
| Median 68% residual half-width | 0.885 mag | 0.681 mag |
| Group-excluded magnitude RMS | 0.994 mag | 0.985 mag |
| Median central-inversion distance discrepancy | 24.32% | 22.03% |

Higher log density means the model assigns more probability to what was observed. The mixture's 9.827-point improvement is descriptive, not a Bayes factor or evidence for a physical redshift model. A post hoc accounting check shows 9.539 points come from NGC4424 alone, with only 0.288 points net improvement across the other 72 galaxies. No observations were removed or refitted for that check. The gain cannot support a claim of broad improvement across environments.

The mixture's nominal 95% interval covers 91.8% of these measurements. Sampling fluctuations, correlated folds and omitted parameter uncertainty prevent a formal calibration verdict from these counts alone. They do not prove that the intervals are suitable for a new population.

## Tail instability remains visible

When NGC4424's group 41220 is excluded, the mixture optimum reaches the declared f=0.5 boundary, with core s about 0.176 mag. This is substantially different from the full-sample 3.38% broad fraction and 0.621-mag core. Its broad component is now about 0.88 mag wide. Different starts also reach a distinct interior local solution with worse likelihood; all start objectives and selected optima are archived. The boundary should not be hidden or interpreted as a precise measured fraction of disturbed galaxies.

NGC4424 still has a -5.311-mag residual under the mixture prediction and remains outside the predicted central 95% interval. Its predicted/reference central distance ratio is 0.0870. Giving this extreme observation a less tiny probability improves log score without making the distance estimate correct. Its broad responsibility is numerically one, but the model has not learned the physical reason for the failure.

## Why this is not yet a distance posterior

The residual intervals include uncertainty in the held-out indicator measurement. A new galaxy without an indicator distance does not have that measurement. Inverting these intervals unchanged would misstate its uncertainty. A distance posterior also needs the population distribution, selection, unknown linewidth distribution, correlated calibration and fitted-parameter uncertainty. Environmental bias must be represented explicitly; a common zero-mean wide tail can absorb dispersion while leaving a false environment correlation intact.

The next useful step is to constrain that population/selection model from the parent measurement sample and the indicator-anchor selection, with environmental disturbance labels or independently justified bias bounds. The current core/broad fit is an intermediate likelihood component, not a complete joint likelihood or a license to assign hard void boundaries. Both candidate error models remain archived; neither is adopted as a validated map generator.

## Verification and goal status

All 348 optimization starts across 116 fits report success, and one selected fit reaches a parameter boundary. This does not prove global optimality. All 146 predictions retain the original groups; NGC4424 is explicitly present in both. Mixture interval roots reproduce their target CDF values to 1.58e-13. The initial fit and the later score-concentration diagnostic were both run successfully, with source hashes recorded. Reproduce using `python research_work/results/tf-probabilistic-calibration/run.py`.

No redshift outcome was fitted or newly queried, no galaxy was rejected, and alpha remains unchanged. This advances stage 2 by supplying and testing an explicit conditional uncertainty model. It does not close selection/environment bias, complete a physical candidate, or establish a model freeze and genuinely withheld test. All four goal stages remain incomplete.
