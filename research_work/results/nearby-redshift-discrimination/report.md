# What can the nearby redshift comparison distinguish?

**The eight-galaxy sample has little expected ability to distinguish our frozen exponential rule from its frozen linear control under the tested uncertainty assumptions.** Their prediction differences are only 6.4-14.4 percent of the per-object distance-plus-spectroscopy uncertainty scale. This is a post-evaluation sensitivity calculation, not a new validation result or a reason to erase the observed residuals.

## Known formulas and assumptions

Established first-order uncertainty propagation gives the exponential rule's distance-error contribution:

`sigma_distance_in_cz = c * alpha * exp(alpha*D) * sigma_D`.

The spectral-center error is multiplied by the same exact frame factor used for its measured redshift. We first combine these in quadrature and assume independent Gaussian errors, fixed model coefficients and no galaxy motion. Finite symmetric distance perturbations check the linear approximation. Adopted distances and earlier predictions remain unchanged; quoted errors are used only for this sensitivity model. Stipulating those distances as exact fictional facts is a different convention and does not make unmodeled galaxy motions vanish.

For two fixed Gaussian predictions with the same assumed covariance C, define

`d_squared = (prediction_1 - prediction_2)^T inverse(C) (prediction_1 - prediction_2)`.

If one fixed rule generates the data, the expected log likelihood ratio favoring it is d_squared/2. This is established simple-hypothesis Gaussian algebra, not a new physics formula. We verify the identity with 100,000 synthetic noise realizations. It is not the observed likelihood ratio of the galaxies and is not a posterior model probability.

## Results and practical meaning

With independent quoted distance/spectral errors and zero motion scatter, d_squared is approximately **0.10265** and the expected log likelihood ratio is only **0.05133** natural-log units. The square root, 0.3204, is a separation measure under this model; it must not be advertised as an observed significance or a cosmological exclusion.

For perspective, adding illustrative independent motion scatter changes the expected discrimination as follows. These values are scenarios, not measured or fitted galaxy-motion dispersions.

| Assumed extra motion scatter (km/s) | Squared rule separation | Expected log likelihood ratio |
|---:|---:|---:|
| 0 | 0.10265 | 0.05133 |
| 100 | 0.01064 | 0.00532 |
| 300 | 0.00140 | 0.00070 |

We also test illustrative correlations of 0.5 and 0.9 among projected distance errors, with every original marginal distance uncertainty preserved. Correlation effects are not necessarily monotonic: common errors can obscure a common shift while leaving some differences better constrained. None of the nine scenarios supplies substantial separation here. No covariance scenario is claimed to reconstruct the actual TRGB/SBF calibration or the galaxies' correlated motions.

The current two rules also have separately calibrated slopes. Their difference is therefore not purely exponential curvature. The saved per-object table separately records how far the exponential lies above its own linear tangent. At these distances that curvature is smaller still. Agreement with a nearly linear distance trend cannot establish the proposed microscopic meaning of alpha.

## Consequence for the larger goal

We should not treat additional small differences between nearby exponential and linear scores as progress toward proving photon conversion. The measured discrepancies remain real comparison results; the inference about their cause remains underdetermined. Seven of the eight galaxies were also below the original calibration range.

The remaining nearby candidates can still provide consistency checks once their source identities and selection are understood. They should not be consumed as decisive theory-validation holdouts without a prediction that actually separates the candidate from alternatives at their distances and uncertainties. No new outcomes were opened in this calculation.

The central research priority is a common physical mechanism that jointly predicts spectral shifts, event timing, companion energy retention and gravity, with detector standards and source/environment parameters specified. Current small-scale propagation constructions have documented failures, and favorable stellar gravity fits do not close those gaps. A more distant sample alone would also require non-redshift distance provenance, calibration and motion/selection treatment; extrapolating our fitted alpha is not a first-principles derivation.

No distances, parameters or prior scores are changed. The total photon-supply budget remains deferred, and the final Cepheid test remains unopened. This calculation changes the interpretation and prioritization of future tests; it does not create a strong case for the theory.

## Reproduction

Run `run.py`. `results.json` contains all eight uncertainty projections, the nine declared covariance scenarios, the synthetic identity check and the hashes of the unchanged evaluated inputs. The source computation checks positive-definite covariance matrices, finite perturbation agreement and the analytic expected-likelihood identity.

- [Frozen eight-galaxy comparison](../nearby-redshift-crosscheck/report.md): original predictions, residuals, source conventions and exposure record.
- [Shared-speed obstruction](../shared-speed-obstruction/report.md): the physical issue that cannot be repaired by parameter-only tuning within that Hamiltonian family.
