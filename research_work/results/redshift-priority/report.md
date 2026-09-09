# Redshift priority: first predictive comparison

The frozen three-candidate comparison retained all 164 groups and left out each of 65 existing sky tiles in turn. The smooth distance-dependent rate did not demonstrate improvement. No new preferred formula is adopted.

| Candidate | Out-of-fold RMS (km/s) | MAE (km/s) | Bias (km/s) |
|---|---:|---:|---:|
| linear | 448.09 | 326.78 | -20.02 |
| constant | 448.38 | 326.68 | -23.55 |
| smooth | 452.52 | 330.91 | -14.06 |

These are c times redshift residuals, not measured velocities. The smooth-minus-constant RMS difference is +4.14 km/s; the descriptive paired sky-tile bootstrap interval is [-0.72, 9.26] km/s. The linear control is also indistinguishable at this precision. This comparison differs from the paper's inherited train/test split and must not be treated as an improvement over that split's RMS.

## Residual/input audit

The original frozen constant fit has a distance correlation of 0.139 and Cartesian sky-direction correlations of 0.256, 0.093 and 0.164. These are descriptive correlations on exposed data, not evidence for a void law. Distance-quartile RMS grows from 366.6 to 506.0 km/s, but this does not identify a cause: selection, calibration, motions and the propagation law remain possible contributors. The direction axes are equatorial Cartesian components, not independently selected physical preferred directions. Adjacent tiles and common calibrations can retain correlations; existing small tiles are not a guarantee of independent errors.

The source is Cosmicflows-4 CDS J/ApJ/944/94 table 2, as documented in the recovered protocol. It selected method-specific surface-brightness-fluctuation distances with 10-150 Mpc range and modulus uncertainty greater than zero and at most 0.30 mag. All 68 earlier pilot groups were excluded; 213 measurements reduced to 164 representatives by minimum SBF error per dominant-PGC group, ties by PGC. No target-redshift cut was used. Actual retained distances span 10.2-93.2 Mpc. The current fixed-distance contract supersedes the historical luminosity-distance remapping: that remapping and historical expansion predictions are not used here.

Redshifts use the published CMB frame convention. The audit has not independently reconstructed the frame correction or underlying SBF calibration. Group definitions also retain catalog assumptions. No independent line-of-sight density, void fraction, gravitational potential or time-field measurement exists in these input columns. An environmental coefficient cannot yet be fitted honestly from this sample alone.

## Observation model and next work

In consistently declared standards, use (1+z_observed) = (1+z_time)(1+z_motion)(1+z_endpoint), with measurement error applied to the measured quantity. This is a factorization convention requiring care about which frame effects were already corrected, not a derivation of the unknown factors. Do not fit one motion per group or count CMB-frame correction twice. Independently constrained motion/calibration uncertainty is needed before declaring a physical precision floor or credible prediction intervals.

Next audit catalog metadata and independent environmental/flow inputs, and identify unexposed distance samples with no target-redshift-derived distance. Freeze fresh selection and exclude known groups before inspecting outcomes. No fresh sample has yet been quarantined, no final uncertainty model selected, and no first-principles field or all-distance accuracy claim is established. The broader goal remains active.

## Reproduce

Run `python research_work/results/redshift-priority/run.py` from repository root. The runner verifies the catalog hash, trains every fold without its held-out tile, saves all 164 out-of-fold predictions and every fold's parameters, and records the pre-run protocol hash. All candidates were reported without data-driven model selection; no nested selection was performed or required in this fixed comparison. Bootstrap intervals are conditional resampling of out-of-fold residuals, not a full refitting bootstrap or independent observational validation.
