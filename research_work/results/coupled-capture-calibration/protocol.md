# Fit the coupled growth model, not its no-feedback approximation

Declare before computation. Keep p=2 or3, Wstar=40000(km/s)^2, radius100kpc, the refined axisymmetric ordinary field, isotropic thin irradiation and permanent in-place retention. Integrate d rho/ds=C*[W_total/(W_total+Wstar)]^p from rho=0 to s=1 using the preceding quasistatic Poisson feedback model. C alone is adjustable, with units Msun/kpc^3 and fixed search range10^4..10^8.5.

Use only previously exposed542 training-star radii and12 Jeans-proxy bins. Fit the unweighted RMS of averaged squared circular-speed predictions. No reserved stars, vertical inferences or lensing observations enter the objective. All observational assumptions and energy/support limitations from the previous report persist.

Use400 radial nodes,64 angular nodes, even multipoles through24 and128 midpoint steps for calibration. Scan13 equally spaced log10 C points on[4,8.5], retain the entire scan, and locally refine each interior scan minimum over its two neighbors with bounded scalar optimization. Compare scan edges as well; flag an edge optimum instead of extending bounds silently. Require optimizer success and report all p variants.

At the chosen C, rerun the same spatial grid with256 steps, and an800-node/128-angle/multipole48 grid with256 steps. Require each prediction difference from the fitted run below0.5km/s; retain failures and refine if necessary. Evaluate C times0.99 and1.01 on the calibration grid to expose sensitivity; do not use that sensitivity to alter the chosen best fit. Record mass, solar deposit potential, saturation range and all12 residuals. Numerical precision is not a confidence interval.

These are training calibrations within a finite-cutoff, thin-source, quasistatic model. An improved RMS is not a proof of source-energy supply, stability under changing illumination/cutoff, physical support, persistent redshift or unseen-data performance.
