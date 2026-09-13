# Uniform background: retention versus galaxy motions

13 September 2026. Controlled comparison using the existing 149-galaxy sample.

## Result

A uniform continuing background can preserve deposits in the chosen threshold kinetics, but cannot simply be added while preserving the fitted rotation curves. The background required to retain 90% of every initially equilibrated deposit increases the frozen validation/test velocity RMSE from 32.49/23.59 to 36.88/30.67 km/s. Refitting the shared amplitude on training galaxies reduces the damage, but the 90% background case still scores worse than an identically refitted no-background control: 31.54/22.52 versus 30.82/21.41 km/s.

This constrains the uniform additive input interpretation with the existing spatial profile. It does not exclude an environment-dependent background, different transport geometry, protected storage, or the companion concept. Nor does the comparison measure any actual external background.

## Hypothesis and derivation

Keep the exact-third reference eta(X)=X^(1/3)/(1+X^(1/3)) and interpret the existing X=L/R_d^2 proxy as a local contribution. Add one common nonnegative dimensionless input B. Initially the source is X+B; after the local component is removed, it is B. Assume linear addition of these rates and the same threshold capture/release kinetics. Both the source-to-rate interpretation and the universal B are new hypotheses.

For initially equilibrated stores, the perpetual retained fraction is

\[
\mathcal R(X,B)=\frac{\eta(B)}{\eta(X+B)}.
\]

The individual threshold populations decrease monotonically after this source reduction, so the late-time ratio is the minimum. This follows from the known linear rate solution, not a new fundamental law. Use the exact reference response here, rather than the finite-cutoff version in the previous floor calculation. Also unlike that experiment, the background contributes to the initial equilibrium as well as the final one.

For each chosen retention target, solve the smallest individual B and select the largest required value as the common background. B is determined by this requirement, not fitted to velocities. All sample input values, including the exposed validation/test inputs, enter the requirement; this is not an untouched-system prediction.

| Required retention in every input | Minimum common B | B/local input range |
|---|---:|---:|
| 50% | 0.210789 | 0.00836-12.10 |
| 90% | 11.256188 | 0.44668-646.20 |
| 99% | 123.890283 | 4.91637-7112.33 |

UGC11914 sets the largest requirement in all three cases. The low-intensity systems become background dominated. This removes much of the sensitivity of retention to local luminosity density: d ln eta(X+B)/d ln X = X/[3(X+B)(1+(X+B)^(1/3))]. That derivative is ordinary differentiation of the reference curve, not a microscopic explanation.

## Rotation calculation

Hold the published project baryonic inputs and reference deposition geometry fixed. For each measured radius, substitute only the retention amplitude:

\[
v_B^2=v_b^2+\frac{\eta(X+B)}{\eta(X)}(v_0^2-v_b^2).
\]

This is a controlled hypothesis about how background loading enters the existing gravity model. Actual external illumination could alter the deposited spatial profile; that is not derived here. The amplitude substitution is an algebraic consequence of the model's linear extra-gravity normalization. No baseline parameter is changed in the frozen comparison.

| Retention target | Frozen train RMSE | Frozen validation RMSE | Frozen test RMSE |
|---|---:|---:|---:|
| No background control | 29.03 | 32.49 | 23.59 |
| 50% | 29.68 | 32.79 | 24.20 |
| 90% | 36.33 | 36.88 | 30.67 |
| 99% | 42.09 | 41.89 | 36.52 |

All entries are km/s, using equal-galaxy mean squared velocity residuals before taking the square root. These are descriptive errors, not likelihoods with propagated observational uncertainties. The 90% case changes an individual predicted speed by as much as 62.43 km/s. It multiplies the extra-gravity component, not total gravity, by up to 3.36.

## Can a shared amplitude compensate?

As a separate diagnostic, multiply the extra-gravity term by a nonnegative a and minimize the equal-galaxy velocity MSE using only the 89 training galaxies. Keep the other two reference constants fixed. Apply the same procedure to the B=0 control. All solutions lie inside the stated search interval 0<a<4.

| Retention target | Fitted C/original C | Train RMSE | Validation RMSE | Test RMSE |
|---|---:|---:|---:|---:|
| No background control | 0.747109 | 26.87 | 30.82 | 21.41 |
| 50% | 0.724837 | 27.06 | 30.82 | 21.46 |
| 90% | 0.570148 | 28.59 | 31.54 | 22.52 |
| 99% | 0.482222 | 28.92 | 31.75 | 22.76 |

The 50% variant has a tiny validation improvement (~0.002 km/s) with worse training and test scores; it is not selected as a winner. The 90% and 99% variants worsen all three scores against the matched adjusted control. Adjusting the amplitude improves the original control too. Consequently, comparing an adjusted background model only against the original frozen control would give a misleading impression of improvement.

This diagnostic uses velocity MSE rather than silently claiming to reproduce the original fitting objective. There is one re-estimated amplitude per scenario, two frozen reference constants, and a background value fixed by a declared retention target and the sample input range. The targets are exploratory choices. This does not amount to a full reoptimization, uncertainty analysis, or fresh observational test. None of these amplitude adjustments is adopted into the reference theory.

## Energy and alternatives

Under the equilibrium energy-capacity convention, final turnover is E_cap lambda B[1-eta(B)]. No physical watts follow without a capacity and rate. Any proposed external flux must support this continuing capture and release. A uniform physical radiation field need not produce the same dimensionless B in every galaxy: geometry, capture cross sections and spectra could change the mapping. This calculation isolates that simplest common-rate option.

The next physical alternatives are a background tied to a measured environment, capture that depends on local properties, or a storage law with less release. Each needs an independently specified rule. Allowing a separate freely chosen B for every galaxy would obscure the issue rather than derive external supply. Refitting spatial parameters could change the scores but cannot itself establish a source or capture mechanism.

## Reproducibility

Run `python research_work/results/companion-extensions/background-retention.py`. The executable verifies original catalogue/archive hashes, checks all 3150 accepted radii and measured velocities, reproduces frozen baseline scores, solves the three retention requirements across 149 inputs, evaluates 596 galaxy/scenario curves, and separately performs four training-only amplitude fits. Results record both raw-data and model-input hashes. No observed data are modified. The known Hill mathematics and rate solutions are distinguished above from proposed physical interpretations.

Related: [constant-floor calculation](threshold-floor-report.md), [finite-threshold comparison](threshold-galaxies-report.md). Full results: `background-retention-results.json`.
