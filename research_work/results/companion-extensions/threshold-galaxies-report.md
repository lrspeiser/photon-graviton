# Finite threshold mixtures on the frozen galaxy sample

13 September 2026. Existing exposed SPARC partitions; no refitting or new blind sample.

## Outcome

A finite threshold distribution can approximate the reference closely on the existing 149 galaxies and 3150 rotation measurements. A 12-decade threshold range changes predicted speeds by at most 0.640 km/s; an 18-decade range reduces that maximum to 0.0634 km/s. A six-decade range changes some predictions by 7.03 km/s and worsens aggregate validation and test velocity RMSE. This does not resolve the reference model's existing discrepancies.

The one-third reference remains the primary empirical prescription. The mixture is an inverse-designed interpretation using known fractional-response mathematics, not an independently derived capture law.

## Frozen comparison

The executable verifies hashes of the original input files, reads the existing catalogue and rotation decomposition, and reproduces every stored radius and observed speed. It retains the original ordinary-matter prescription, distance choices, companion shape, scale, opacity and normalization. It does not run an optimizer.

Let eta_0 be the reference retention and eta_W the normalized mixture truncated to thresholds [10^-W,10^W]. In the frozen model, eta multiplies the companion density while the shape is unchanged. Therefore linearity of its gravitational acceleration gives

\[
v_W^2(R)=v_b^2(R)+\frac{\eta_W}{\eta_0}\left[v_0^2(R)-v_b^2(R)\right].
\]

This rescaling is a derived consequence of the existing amplitude-only retention prescription, not a new gravity law. The source proxy X=(L_3.6/10^9)/(R_d/kpc)^2 spans approximately 0.01742 to 25.1995 in this sample. The threshold mixture was defined before this comparison; the three cutoff choices are sensitivity cases, not fitted parameters or a selected winner.

The labels train, validation and test preserve the historical 89/29/31 split. All partitions have already been examined in this research, so these results are not fresh out-of-sample evidence.

## Results

RMSE is the square root of the mean of per-galaxy mean squared speed residuals, matching the original equal-galaxy weighting.

| Threshold range | Train RMSE (km/s) | Validation RMSE (km/s) | Test RMSE (km/s) | Maximum speed change (km/s) |
|---|---:|---:|---:|---:|
| Exact reference | 29.0251 | 32.4949 | 23.5908 | 0 |
| 10^-3 to 10^3 | 29.1489 | 33.0540 | 23.9520 | 7.0273 |
| 10^-6 to 10^6 | 29.0325 | 32.5406 | 23.6186 | 0.6398 |
| 10^-9 to 10^9 | 29.0258 | 32.4993 | 23.5935 | 0.0634 |

The test logarithmic residual metric improves slightly for some truncations even while velocity RMSE worsens; both metrics are retained in the JSON. No general improvement is claimed. Broadening the distribution approaches the reference by construction, which is not evidence that astrophysical thresholds actually span that range.

## Implications

The broad finite mixture is numerically compatible with the existing reference at a level small compared with its overall residuals. It supplies neither improved rotations nor a source-energy budget. The required threshold population still needs an independent environmental explanation. Its equilibrium occupied capacity is not automatically the fraction of incoming companion energy retained permanently.

The next physical task is to derive a distribution or a release history and test a consequence beyond reproducing the reference curve. A future galaxy test must freeze those rules before inspecting genuinely new systems. Lensing, vertical motions, propagation and energy supply remain separate requirements; this amplitude-only substitution does not solve them.

## Verification

All original input hashes match. The same 149 galaxies and 3150 accepted measurements are retained. The unmodified curve reproduces the stored finer-grid split RMSE values within 1e-8 km/s, and its pointwise speeds within 1e-10 km/s. All 596 galaxy/model records include predicted curves and per-galaxy metrics. No raw observations or prior fitted parameters are changed.

Files: `threshold-galaxies.py`, `threshold-galaxies-results.json`. Related: [threshold-mixture derivation](threshold-mixture-report.md).
