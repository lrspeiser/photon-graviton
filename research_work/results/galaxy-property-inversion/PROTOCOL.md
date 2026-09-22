# JR-2: fixed-formula inverse galaxy-property analysis

21 September 2026 (America/Los_Angeles). Baseline main df6bb0a20062aa3bbcdfecccbecba1cc3a80ae7f. Freeze the exact R10 formula, every universal parameter, and original JR-1 evidence. This is post-data exploratory inverse inference, not independent validation or a new measurement. This protocol was written locally before execution; its GitHub publication occurs while the locally declared calculation is running.

## Theory first and user goal

The companion hypothesis seeks one connected explanation of motion and lensing. Here the user requests a different direction of inference: hold the formula fixed, find the galaxy properties needed to match the observations, and identify shared changes. This does not promote R10 to the full microscopic theory. Preserve original observations, profiles, errors and archived predictions; write altered properties only into separately labelled counterfactual results.

## Scope

All 149 archived SPARC galaxies and six SLACS galaxies. The supplied JR-1 forward model does not implement clusters; no cluster result is inferred. Run all galaxies and separately summarize the original greater-than-20-percent fractional-RMS subgroup, so the sample is not selected only after seeing improvements.

## Disk-property inversions

Individually vary distance by factors 0.2 to 5, inclination from 5 to 90 degrees, stellar mass-to-light normalization by factors 0.1 to 10, and gas normalization by factors 0.25 to 4. These wide ranges diagnose what the formula demands, not credible observational errors.

Also fit all four jointly within two quoted marginal catalog errors for distance/inclination, 0.30 dex for stellar mass-to-light and 0.20 dex for gas. Penalize catalog distance/inclination departures with their quoted errors, and stellar/gas log offsets with working scales of 0.15/0.10 dex. The latter are explicit sensitivity assumptions, NOT per-object measurements. A wide joint fit without these penalties diagnoses remaining shape disagreement. Retain all bounds, optimizer status and unresolved cases.

Distance factor d changes radii and half-light radii by d, photometric stellar and gas masses by d^2, and tabulated baryonic squared speeds by d. Population/gas multipliers propagate into both ordinary gravity and the unchanged source-linked companion law. Changing inclination preserves the reconstructed projected velocity v_obs*sin(i_catalog) and its error, comparing with v_model*sin(i_new). This is an approximate major-axis projected proxy, not re-extraction of raw 2D velocity maps or re-deprojection of the image geometry. No independent radial corrections or companion amplitudes are fitted.

Use the exact R10 disk/spheroid relation. Record one-at-a-time and joint changes, physical fractional errors, original and fitted properties, standardized catalog shifts, and a common four-parameter property-calibration control. A common shift is an in-sample diagnostic, not a measured recalibration. Correlations with inclination, surface brightness and data quality are descriptive, not corrected for all selection/systematic effects.

## Lens-property inversions

Freeze angular light profiles, observed stellar V_rms and covariance, and Einstein-angle image-model summaries. Start with R10's inferred stellar mass and orbital beta per object. At fixed conditional distance geometry, compare: mass alone fitted to the ring; mass alone fitted to stars; mass and beta jointly; and mass/beta fitted to stars followed by the additional external focusing needed for the ring. Allow stellar mass offsets within 0.30 dex and beta in [-1,0.35]. Published mass errors are reference comparisons, not a complete posterior under the alternate propagation convention. Positive Jeans moments do not prove a valid distribution function.

Every changed mass enters ordinary gravity, the companion source, motion and light bending together. The joint objective retains the previous 5-percent working deflection scale, not a measured lens error. No new gravity coefficient or per-object light-response multiplier is permitted.

Required effective external convergence is kappa = 1 - (Dls/Ds)*alpha_internal(theta_observed)/theta_observed. This algebraically closes one ring condition, not full lens images or a physical environment model. Report the alternative required Dls/Ds ratio and its physical range. Do not claim a foreground/background object exists because its required focusing was calculated. The environment approximation assumes negligible effect on the local stellar dynamics; a real spatial model must justify it.

## Interpretation

A residual-driven property shift is not evidence that the catalog is wrong. Distinguish a small independently allowed shift, a shift allowed only by stipulated sensitivity bounds, and an extreme diagnostic adjustment. Opposite corrections for different objects are not one common systematic. Even a perfect inverse fit would not establish microscopic conversion, adequate source energy, Solar-system behavior or cluster physics. Verify original data and formula hashes remain unchanged.
