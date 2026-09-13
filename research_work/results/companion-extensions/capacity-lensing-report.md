# Capacity branches transferred to stellar motions and lensing

13 September 2026. Matched transfer to the six previously used lens systems.

## Outcome

Neither local filling nor full recycling resolves the joint motion/lensing problem. With their amplitudes frozen from SPARC training, local filling has lens-angle fractional RMS 13.29%-13.46%. Full recycling gives 13.05%-13.16% after refinement, only slightly below the matched adjusted reference at 13.07%-13.19%. Recycling improves the absolute lens residual in only three of six systems under each population assumption. Its known SPARC test error is substantially worse than the reference, and the outer stellar-bin discrepancies remain large.

This is a conditional comparison on already-exposed targets, not a blind validation. No change is adopted into the reference. The remaining ~13% lens-angle residual is not an observational uncertainty or a calibrated significance measure.

## Consistent profiles and parameters

Use A=2 C0, with the corrected normalization from the original fit. The three profiles are

\[
\rho_{\rm ref}=A\eta(X)Jg,\qquad
\rho_{\rm local}=A\eta(XJ)g,\qquad
\rho_{\rm recycle}=A\eta(X)g,
\]

where g=[1+(r/a)^2]^-2 and eta is the exact-third reference Hill response. The local branch uses the prescribed full-opacity field, while recycling is the stationary isotropic same-channel return limit. They remain different physical branches; this comparison does not merge them into one mechanism.

The C0 multipliers are 0.747108933 for the adjusted reference, 0.352688398 for local filling, and 0.277799531 for full recycling. They are read directly from earlier galaxy training-only amplitude fits. The original unadjusted reference, multiplier one, is also rerun as a verification control. Other shared capture parameters are held fixed. No lens target sets a capture amplitude, exponent, profile scale or optical parameter.

For each system and each population proxy, only stellar mass and constant anisotropy are fitted to the inner stellar Vrms bins using their released covariance. The outermost bin is held out of that nuisance fit and evaluated with the same covariance-conditioned residual as the historical test. The same total mass profile supplies spherical stellar accelerations and the weak-field deflection integral. This is a consistent effective gravitational response within the conditional model, not a relativistic derivation of companions.

## Aggregate comparison

Each row contains six systems. Inner chi-squared and the sum of squared covariance-conditioned outer residuals are shown alongside lens-angle fractional RMS. They measure different residuals and are not pooled into one undocumented score. Neither stellar population is selected as the preferred result.

| Profile | Population proxy | Inner chi-squared sum | Outer standardized residual square sum | Lens-angle RMS (%) |
|---|---|---:|---:|---:|
| Original reference | Chabrier | 42.16 | 44.32 | 13.923 |
| Original reference | Salpeter | 41.55 | 42.70 | 14.094 |
| Adjusted reference, refined | Chabrier | 45.81 | 55.02 | 13.068 |
| Adjusted reference, refined | Salpeter | 45.27 | 53.42 | 13.187 |
| Local filling | Chabrier | 41.63 | 50.41 | 13.287 |
| Local filling | Salpeter | 40.86 | 48.69 | 13.461 |
| Full recycling, refined | Chabrier | 42.90 | 54.93 | 13.045 |
| Full recycling, refined | Salpeter | 42.34 | 53.59 | 13.161 |

Reducing the reference amplitude itself improves lensing compared with the original reference, while worsening the outer stellar residuals. Therefore crediting that entire lens improvement to a capacity mechanism would be misleading. Compared with the matched adjusted reference, recycling improves lens RMS by only 0.0231/0.0260 percentage points. Its outer residual sum is marginally lower for Chabrier and marginally higher for Salpeter. Local filling improves the inner and outer stellar summaries relative to that adjusted reference but worsens lens RMS. No branch simultaneously removes these discrepancies.

The corresponding previously measured SPARC test RMSE values are 21.41 km/s for the adjusted reference, 25.19 for local filling and 30.88 for recycling. These cross-system tradeoffs matter: a tiny lens-only improvement does not establish superiority when the galaxy test comparison worsens.

## Per-system predictions

The machine-readable result records every lens angle, stellar mass, anisotropy, inner chi-squared and conditional outer residual for each branch and population. The following angles use refined reference/recycling runs and the standard local run; each model cell lists Chabrier / Salpeter, in arcseconds. The catalogue angle is the published SIE model summary, not a direct model-free observable.

| System | Catalogue | Adjusted reference | Local filling | Full recycling |
|---|---:|---:|---:|---:|
| J0037-0942 | 1.530 | 1.79882 / 1.80314 | 1.82746 / 1.83555 | 1.82015 / 1.82553 |
| J1112+0826 | 1.490 | 1.45716 / 1.45962 | 1.45385 / 1.45669 | 1.44564 / 1.44756 |
| J1204+0358 | 1.310 | 1.51696 / 1.51861 | 1.50598 / 1.50737 | 1.50052 / 1.50151 |
| J1402+6321 | 1.350 | 1.62664 / 1.62835 | 1.62599 / 1.62797 | 1.62158 / 1.62297 |
| J1621+3931 | 1.290 | 1.36251 / 1.36405 | 1.35460 / 1.35487 | 1.35208 / 1.35258 |
| J1630+4520 | 1.780 | 1.72369 / 1.72729 | 1.72783 / 1.73264 | 1.71540 / 1.71857 |

## Data and geometry limitations

Both population cases use a stellar-population-mass-to-3.6-micron-luminosity conversion; they do not supply directly measured rest-frame 3.6-micron luminosities for these targets. X remains a proxy. Angular distances and lens-source distance ratios use the project's existing nonexpanding optical prescription, which is itself conditional and has not been independently established by this test. Published lens angles summarize an image model; no raw-image lens likelihood or source reconstruction is fitted here.

The spherical deprojection, fixed light components, common mass-to-light ratio and constant anisotropy assumptions are inherited. Source histories, external fluxes, three-dimensional shape, stellar-population gradients and the microscopic release/support laws are not derived. A full comparison to appropriate halo and MOND models requires matched inputs and information budgets beyond this capacity-branch transfer. Previous halo comparisons remain separate records.

## Numerical and provenance checks

The original-capacity command exactly reproduces the historical summary and per-system results to the checks recorded by the verifier. All new runs retain identical measured stellar arrays, population mappings, catalogue angles and conditional geometry. Each nuisance fit has successful optimizer starts and no anisotropy boundary flag.

Because the recycling/reference lens difference is small, those two branches are rerun with 8001 rather than 4001 radial points, 256 rather than 128 projection angles, 512 rather than 256 deprojection points and 192 rather than 96 capture angles. The maximum lens-angle change is 7.73e-5 arcseconds; the maximum change in conditional outer prediction is below 0.00076 km/s. Fractional lens RMS changes by at most 3.33e-5 (0.00333 percentage points), smaller than the small branch difference. This supports numerical resolution of the difference, not its observational significance or robustness to the much larger model uncertainties. The local branch is not independently refined here.

The radial mass grid extends beyond 21919 capture scales in each system, retaining the historical exterior continuation. No new physical halo cutoff is inferred from this numerical boundary. Formation, dynamical stability, radiation forces, photon timing and lensing-geometry validity are not tested by these convergence checks.

## Reproduction

Run `python research_work/results/isotropic-galaxy-transfer/lensing.py --capacity-branch=BRANCH` for original, reference, local and recycling. Run reference and recycling again with `--capacity-refine`. Then run `python research_work/results/companion-extensions/capacity-lensing.py` to verify and summarize. Outputs use new `capacity-...-optics-results.json` names and preserve the historical files. The summary records hashes for raw/derived inputs, model code, source calibrations and all six output runs. See the [pre-execution protocol](capacity-lensing-protocol.md).

The result is an executed cross-scale assessment with explicit failures and small tradeoffs. It does not complete the broader theory or identify a preferred physical capacity mechanism.
