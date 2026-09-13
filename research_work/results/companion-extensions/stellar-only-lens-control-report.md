# Matched stellar-only control for the lens-calibrated motion fits

13 September 2026. Same four stellar freedoms and orbital constraints, with the companion density removed. All observations are exposed and fitted.

## Outcome

The current companion profile improves the aggregate motion fit relative to the matched stellar-only control. Total chi-squared is 175.09 for stellar-only, versus 113.56/111.27 for companions under Chabrier/Salpeter proxy assumptions: reductions of 35.1% and 36.4%. The improvement is concentrated in J0037-0942 and J1402+6321. The other four systems fit slightly better without the companion profile.

This distinguishes added-density effects from the improvement supplied by flexible stellar modeling. It does not establish a radiation-derived origin, the one-third exponent, or a joint physical solution. J1402+6321 remains poorly described even with companions. The result is conditional on the adopted spherical geometry, stellar gradient/orbit forms, population proxies, light profiles and nonexpanding optical geometry.

## What is matched

Both versions use the same six galaxies, forty motion bins, covariance, PSF, light profiles and catalogue lens angles. Both fit four stellar nuisance parameters: mass-to-light gradient h, central and outer anisotropy, and orbital transition radius. The same parameter bounds and necessary central/finite-radius orbital conditions are imposed. All motion bins are optimized, including the formerly reserved outer bin.

At each gradient, the lens angle calibrates total stellar mass. Removing companion bending therefore requires recalculating the stellar normalization rather than leaving the previous stellar mass fixed. The same recalculated mass enters the stellar accelerations and lensing. No unrelated lensing multiplier is fitted.

The control contains only the modeled stellar mass as a gravitational source. It does not contain separate gas or central black-hole components, which are also absent as separate components from the companion comparison. Thus this is a matched stellar-only ordinary-matter control, not an exhaustive baryonic model for every galaxy.

The companion version has a specified extra density with shared parameters trained in prior galaxy work; those parameters are not optimized on these lens systems. The stellar-only control removes that component. Equal stellar freedoms do not imply equal total theoretical parameter counts: the companion branch still inherits its global parameters and physical assumptions. No likelihood-based significance or model-comparison probability is assigned here.

## Formula and provenance

Introduce a control switch s, equal to one for the recorded companion profile and zero for stellar-only:

    g_total(r) = g_stars(r; M_star, h) + s*g_companion(r).
    alpha_required = alpha_stars(M_star, h) + s*alpha_companion.

These are standard superposition operations in the existing spherical effective-gravity/weak-field lens calculation, not new photon-conversion equations. Solving the lens equation determines M_star for each h and s. The change is a diagnostic ablation of the proposed component; it is not a revised energy-conversion mechanism. Four stellar parameters are then optimized against the full motion covariance under either switch value.

The fitted companion density remains the fixed exact-third reference, including its corrected A=2 C0 normalization. Turning its effective gravity off for the control does not silently convert its stored energy into stars or alter the companion branch's energy ledger. It defines a separate comparison model.

## Results

| Population proxy | Stellar-only total chi-squared | Companion total chi-squared | Reduction with companions | Galaxies improved by companions |
|---|---:|---:|---:|---:|
| Chabrier | 175.088 | 113.565 | 35.14% | 2 / 6 |
| Salpeter | 175.088 | 111.274 | 36.45% | 2 / 6 |

The two proxy labels produce the same stellar-only physical problem: their differing luminosity proxy only supplied the companion normalization in this comparison, and the stellar mass is lens-calibrated. The duplicated runs therefore check the zero-companion limit; they are not independent observations.

| Galaxy | Population | Stellar-only chi-squared | Companion chi-squared | Stellar-only minus companion |
|---|---|---:|---:|---:|
| J0037-0942 | Chabrier | 32.416 | 10.444 | 21.972 |
| J0037-0942 | Salpeter | 32.416 | 9.739 | 22.678 |
| J1112+0826 | Chabrier | 18.025 | 18.668 | -0.643 |
| J1112+0826 | Salpeter | 18.025 | 18.700 | -0.675 |
| J1204+0358 | Chabrier | 13.945 | 14.333 | -0.388 |
| J1204+0358 | Salpeter | 13.945 | 14.424 | -0.479 |
| J1402+6321 | Chabrier | 105.958 | 64.974 | 40.984 |
| J1402+6321 | Salpeter | 105.958 | 63.243 | 42.715 |
| J1621+3931 | Chabrier | 3.094 | 3.409 | -0.315 |
| J1621+3931 | Salpeter | 3.094 | 3.426 | -0.332 |
| J1630+4520 | Chabrier | 1.650 | 1.737 | -0.087 |
| J1630+4520 | Salpeter | 1.650 | 1.742 | -0.093 |

Positive differences favor the current companion profile descriptively; negative differences favor the control. The aggregate gain is driven by the two systems with large positive differences, not uniform improvement across the sample.

| Galaxy | Stellar-only h | beta0 | beta_infinity | ra/Re | Conditional outer residual |
|---|---:|---:|---:|---:|---:|
| J0037-0942 | -0.7667 | 0.3750 | 0.9500 | 0.5536 | 3.617 |
| J1112+0826 | -0.2198 | 0.3750 | -2.0000 | 0.4658 | 3.203 |
| J1204+0358 | 5.1264 | -2.0000 | 0.9500 | 1.4526 | 2.906 |
| J1402+6321 | -0.6874 | 0.3750 | 0.9500 | 0.8630 | 6.954 |
| J1621+3931 | -0.6334 | 0.2252 | 0.4206 | 0.7659 | -0.442 |
| J1630+4520 | 0.0865 | 0.2919 | -2.0000 | 0.5916 | 0.421 |

The stellar-only outer residual-square sum is 80.51 and its inner contribution is 94.58. With companions these are 51.87/50.84 and 61.70/60.44. These separate contributions are correlated-model diagnostics, not independent evidence categories. Conditional outer residuals are based on the released covariance; they are not model-independent exclusion levels.

Both versions retain boundary dependence. The control reaches an outer anisotropy endpoint in five systems and the beta0=-2 bound in J1204+0358. It does not reach the stellar-gradient bounds. These differences in stellar parameters reflect competing fits, not observed changes in stellar populations or orbits.

## Execution and verification

The control uses the previous 27 starts plus the companion variable-radius optimum evaluated in the stellar-only potential, for 28 starts per case. Twenty-seven or twenty-eight converge; all successful objectives and failure messages are retained. Every chosen result is no worse than that feasible reference configuration after recalculating its lens-required stellar mass. Local minimization does not establish a global optimum over all orbit families.

The code explicitly reevaluates the recorded reference parameters with companion gravity restored. It reproduces the variable-radius reference objectives exactly at recorded precision and retains the earlier nested-model checks. The full-covariance decomposition into inner and conditional-outer contributions agrees to within 4.27e-14. The stellar-only lens equation closes below 1e-10 fractional error, and stellar masses and motion moments are positive.

All control fits pass the imposed central and refined finite-radius slope checks; the minimum finite margin is 0.0777244. Positivity of the full distribution function and dynamical stability remain unproved. The two zero-companion population runs differ in predicted motions by at most 1.40e-6 km/s and in stellar mass by 5.01e-8 fractionally, consistent with optimizer tolerances for the identical control problem.

## What we learn

The added profile is doing useful gravitational work in this conditional comparison, particularly for J0037-0942 and J1402+6321. The recent improvements were not entirely due to extra stellar freedom. However, two improved galaxies do not establish the density's proposed radiation source, its retention law or universality. Four systems do not prefer the added component in this fit, and the dominant residual remains unresolved.

A matched dark-halo or alternative-gravity comparison would need the same stellar freedoms, inputs and treatment of lensing calibration, together with explicit global and local parameter budgets. This control does not provide that ranking. Independently constrained stellar populations and physically constructed stellar distribution functions would reduce ambiguity before interpreting any residual difference as new physics. Photon propagation, absolute companion supply and persistence remain separate unsolved requirements.

## Reproduction

Run `python research_work/results/companion-extensions/gradient-orbits.py --slope-constrained --all-motion-bins --free-orbit-radius --stellar-only`. The command writes stellar-only-lens-control-results.json and preserves previous results. The file contains the zero-density switch, all observed and modeled motion bins, matched companion objectives, control objectives evaluated at reference stellar parameters, orbital conditions, optimizer outcomes and input hashes. Legacy comparison fields retain their identified companion-model context. Pre-execution choices are in stellar-only-lens-control-protocol.md.
