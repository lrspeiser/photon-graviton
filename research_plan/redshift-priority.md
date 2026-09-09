# Current priority: predict observed redshifts accurately

User direction: 9 September 2026. Make agreement with observed redshifts the immediate research priority. Preserve the broader companion-energy program, but defer new capture/halo research unless needed to make or check redshift predictions.

## Latest constraint

The user requires a physical approach independent of cosmic expansion, not novelty claimed for an existing exponential relation. Follow [the mechanism-first requirements](mechanism-first-redshift.md). Empirical tuning is a diagnostic; the central task is to derive alpha from a shared interaction and measured or independently specified inputs.

## Objective

Develop the simplest shared environmental-time propagation law that predicts catalog redshifts at the user-stipulated published distances, quantifies unresolved scatter, and improves prediction on data excluded from fitting. Do not equate a close fitted curve with confirmation of the mechanism. No dark-matter, expansion or Big-Bang explanation is adopted.

Root: z_time = exp(integral alpha ds) - 1, with alpha = partial_t n / c0 in the current ray candidate. Constant-alpha reference: 0.0002488993286382367 per Mpc. All conversion energy lost goes to companions in the declared common energy standard.

## Work queue, in order

1. Audit the 164-group residuals against distance, sky direction and independent available environment measurements. Record catalog origin, reference frame, selection and group membership. Explain whether deviations have repeatable structure before changing the formula. Existing labels are exposed and cannot be reused as blind evidence.
2. Freeze a small candidate set before fitting it. Compare the constant-rate baseline to a smooth distance-dependent rate as an empirical diagnostic. Add an environmental term only when its line-of-sight input is independently available, not reconstructed from the redshift residual. A candidate form is alpha(s,t) = a0 + a_env f_env(s,t), with a specified dimensionless, independently measured f_env and positive total alpha over its declared domain. This is optional, not an adopted physical law. Do not assign an adjustable rate to each object.
3. Run grouped cross-validation on the exposed sample with all preprocessing and parameter fitting inside each training fold. If selecting models or regularization, use inner folds and reserve outer folds for evaluation. Group correlated objects and sky regions; report that this remains exploratory because research has already seen the sample. Compare to the linear mathematical control as well as constant-alpha. Retain all rows and report failed candidates.
4. Audit additional nearby and more distant catalogs before using their outcomes. Keep the user-adopted distance values, but flag any distances derived directly from redshift: those rows cannot independently validate a redshift-distance relation. Quarantine a genuinely unexposed sample with fixed selection and overlap exclusions before fitting the final model. Do not claim coverage of all distances from the current 10.2-93.2 Mpc sample.
5. Specify an observation model: spectral shift includes conversion plus motion, endpoint and measurement effects. Combine redshift factors consistently with declared frame standards. Do not tune an independent velocity for each object to erase residuals. Derive or independently constrain nuisance effects and audit any alternative-theory assumptions in supplied corrections.
6. Freeze the chosen formula, parameters, intended distance range and uncertainty model before opening the fresh evaluation sample. Report predictions, signed residuals, RMS, median/mean absolute error, bias, interval coverage and residual trends by distance and environment. Compare paired errors with uncertainty that respects group/sky correlations. Set precision targets from the catalog and independently justified nuisance model before fresh evaluation, not from whichever result appears favorable.
7. For any successful empirical candidate, solve for a realizable time field producing the required alpha. Check clock standards, different wavelengths, line widths, image fidelity, whole-event timing and brightness using the same law and energy account. A fitted environmental or distance dependence does not supply this derivation automatically.
8. Update the paper, full object-by-object comparison table, plots and GitHub main with both improvements and limitations.

## Completion criteria

A reproducible, frozen common model with demonstrated predictive improvement on fresh data relative to the existing baselines; documented performance across the actual sampled distance ranges; credible uncertainty coverage and no hidden per-object tuning. Report the remaining precision floor and physical gaps even if empirical prediction improves. Full fundamental-theory confirmation and galaxy gravity remain separate, unfinished requirements. If candidate extensions do not improve predictive performance, record that outcome and investigate missing inputs rather than escalating fit flexibility indefinitely.

## Baseline

The 164 previously exposed groups have inherited train/validation/test sizes 104/35/25. Their exponential RMS residuals are 457.6/437.1/415.4 km/s in c times redshift units. The inherited test linear-control RMS is 414.7 km/s. These are descriptive starting values, not success thresholds or new independent validation. See the version 0.2 manuscript and paper analysis outputs.
