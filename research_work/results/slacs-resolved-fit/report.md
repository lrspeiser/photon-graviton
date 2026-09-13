# Resolved stellar-motion fit followed by lensing prediction

## Result: a partial empirical improvement, not a solution

The fixed extra-force relation reduces the radial kinematic discrepancy in all seven release-approved training galaxies. The sum of per-galaxy chi-square falls from 205.74 to 111.22 for the same 47 measured radial bins and 14 fitted mass/orbit parameters. Several individual fits nevertheless remain poor. The comparison does not provide a successful joint motion-and-lensing explanation.

The predicted-to-catalog lens-angle ratio has median 1.075 for the extra-force branch, versus 1.034 for the mass-following-light baseline. Their RMS fractional angle discrepancies are 13.68% and 13.24%, respectively. Lensing was not used in the fit. Thus improved radial kinematics do not translate into improved overall lensing in this pilot.

| Training galaxy | Baseline radial chi-square | Extra-force radial chi-square | Extra-force fitted beta | Extra-force lens-angle ratio |
|---|---:|---:|---:|---:|
| J0037-0942 | 27.41 | 10.49 | 0.158 | 1.082 |
| J1112+0826 | 56.13 | 37.84 | 0.040 | 0.866 |
| J1204+0358 | 7.61 | 5.59 | 0.073 | 1.122 |
| J1402+6321 | 54.13 | 28.76 | 0.201 | 1.075 |
| J1538+5817 | 30.62 | 15.36 | 0.029 | 1.252 |
| J1621+3931 | 6.31 | 3.45 | 0.258 | 0.934 |
| J1630+4520 | 23.53 | 9.73 | 0.294 | 0.865 |

These are descriptive training results. The chi-square sum assumes block-diagonal per-galaxy covariance and is not assigned a global significance: cross-galaxy calibration covariance, source-model uncertainty and PSF/light-profile errors are not included. The nominal bin-minus-parameter count is 33, not a calibrated goodness-of-fit reference distribution. The force models are not selected by a new likelihood-ratio significance claim. No validation or test score is opened.

## Data and assumptions

We use the audited KCWI radial V_rms profiles and full within-galaxy covariance matrices. [The pinned public release](https://github.com/TDCOSMO/TDCOSMO2025_public/tree/d7f38db341f68be1df0d9ac1fc528c45113f94cf/ExternalLenses/SLACS) marks seven of our training objects usable. J0330-0020 remains explicitly excluded because its use flag is zero, with the underlying reason unresolved; it is not counted as a successful fit or discarded based on our residuals.

We retain the prior pilot's spherical Hernquist light sizes, constant mass-to-light ratio, static conditional distances, equal metric potentials, fixed empirical extra-force coefficients and cutoff 20a. Here we use the release's 0.8-arcsecond Gaussian FWHM and actual angular annuli, rather than the earlier SDSS aperture. Constant orbital beta is fitted in [-0.5,0.45], and mass in [10^7,10^14] solar masses. All reported beta optima are interior to that chosen range. These are working bounds, not measured physical exclusions; formal beta intervals have not been calibrated.

The seven galaxies are a selected overlap sample, not a representative survey. A spherical nonrotating second-moment model remains approximate for galaxies with rotation or flattening. Neither the optimized baseline mass nor the optimized extra-branch stellar mass is independently established from stellar populations. The baseline can absorb extra gravity in its free normalization.

## Formula status and fitting procedure

All Jeans, probability, projection and Gaussian-likelihood equations are known mathematics. The only extra force is the previously fitted empirical relation g_c=A a_star (g_b/a_star)^p, continued beyond its declared cutoff with inverse-square acceleration. It is not a photon-supply or capture-derived deposit profile.

We compute the constant-anisotropy pressure

    nu sigma_r^2 = r^(-2 beta) integral_r^infinity nu(s) g(s) s^(2 beta) ds.

For each annulus, the seeing probability equals the probability of reaching its outer circular aperture minus that of reaching its inner aperture. Both the luminosity weight W and anisotropic projection weight T are differenced **before** forming the second-moment ratio:

    V_rms,annulus^2 = integral r^2 nu sigma_r^2 (W-beta T) dr
                     / integral r^2 nu W dr.

Subtracting cumulative dispersions themselves would be incorrect. The observed radial quantity is V_rms, so the likelihood residual is observed V_rms minus the square root of this second moment, not a residual in squared velocity.

For residual vector r and released covariance C, we use

    log L = -0.5 [r^T C^(-1) r + log det C + N log(2 pi)].

The calculation uses a Cholesky factor, without explicitly inverting C or adding extra scatter. Three optimizer starts vary beta; the best successful solution is retained and every start is recorded. Mass and beta are fitted **only to radial kinematics**. The unchanged weak-field lensing solver then predicts the Einstein angle from that same total force. No independent lensing multiplier, per-object force exponent or lensing-driven orbit choice is introduced.

## Numerical verification

Flux-weighted annular second moments reconstruct the previously checked cumulative aperture at beta=-0.3,0,+0.3 to relative difference below 2.3e-16. For every final fit, radial and angular quadrature are doubled from 4001/128 to 8001/256 nodes at the fitted parameters. The largest velocity change is 3.18e-5 fractionally, and the largest covariance-whitened change has norm 0.00374. These are tiny compared with the observed discrepancies. This is a fixed-parameter numerical refinement, not a complete optimizer or astrophysical uncertainty audit.

## What remains

This provides a reason to keep investigating the empirical force shape: its radial residuals decrease consistently in this small training overlap. It does not identify photon conversion as the cause, demonstrate the model's physical source, or make the remaining radial and lensing discrepancies disappear.

The next physical/modeling work is to replace the crude light profile with documented photometry, constrain stellar population masses and investigate flattening, rotation, orbital flexibility and PSF uncertainty. These changes must be justified by observations rather than chosen simply to force lensing agreement. A capture-derived deposit source, its metric response and energy funding remain missing from this pilot. The broader redshift/timing, Milky Way, clusters and distinguishing-prediction objectives remain open.

## Reproduction

Run `python research_work/results/slacs-resolved-fit/verify.py` and `python research_work/results/slacs-resolved-fit/run.py`. Protocol choices, fitted parameters, per-bin predictions, optimizer starts, refinement results, exclusions and input/code hashes are archived. No external release notebook or serialized posterior is executed.
