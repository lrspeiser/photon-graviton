# Refitting resolved motions with published stellar-light components

## Outcome

For the same six training galaxies, replacing the older Hernquist approximation with the published image components lowers the extra-force radial chi-square from 95.86 to 85.87. Its lensing RMS fractional discrepancy changes from 10.60% to 11.60%, so the better light inputs do not solve the joint motion/lensing problem. The previous seven-object totals are not used for this comparison: J1538+5817 has no supplied detailed components and remains explicitly missing.

| Same six systems | Radial chi-square sum | Lens-angle RMS fractional discrepancy | Orbit-boundary fits |
|---|---:|---:|---:|
| Earlier light profile, free-mass baseline | 175.12 | 11.56% | 0 |
| Earlier light profile, extra force | 95.86 | 10.60% | 0 |
| Published components, free-mass baseline | 170.60 | 13.93% | 4 |
| Published components, extra force | 85.87 | 11.60% | 0 |

There are 40 radial bins and 12 fitted parameters in each branch. These descriptive sums use released within-object covariance, without cross-object covariance or image-profile uncertainty. They are not assigned a formal global significance. Four baseline solutions reach beta=-0.5, limiting interpretation of the score difference between force models. A wider tangential-orbit investigation remains necessary; the boundary is retained rather than called an interior best fit.

| Galaxy | Baseline radial chi-square | Extra-force radial chi-square | Extra-force beta | Extra-force predicted / catalog lens angle |
|---|---:|---:|---:|---:|
| J0037-0942 | 28.57 | 7.30 | -0.075 | 1.028 |
| J1112+0826 | 48.80 | 33.42 | -0.125 | 0.834 |
| J1204+0358 | 16.76 | 7.17 | -0.129 | 1.032 |
| J1402+6321 | 46.22 | 24.15 | 0.100 | 1.044 |
| J1621+3931 | 6.79 | 3.42 | 0.143 | 0.866 |
| J1630+4520 | 23.45 | 10.40 | 0.244 | 0.823 |

## What is fixed and what is fitted

The six image profiles, their component radii, relative amplitudes and Sersic indices are fixed to the [pinned published metadata](https://github.com/TDCOSMO/TDCOSMO2025_public/blob/d7f38db341f68be1df0d9ac1fc528c45113f94cf/ExternalLenses/SLACS/slacs_all_params.csv). They are not adjusted to radial or lensing residuals. Their interpretation was checked in the preceding light-profile audit.

The calculation makes a new explicit spherical deprojection of the equal-area radial profiles, discarding their angular flattening for this conditional comparison. Each galaxy has a single constant mass-to-light ratio shared by its components; no component is assigned a separately fitted gravity normalization. Both tracer weighting and ordinary-matter acceleration use that same deprojected profile.

The total mass and constant orbital beta are fitted to radial V_rms only, with the previous bounds and starts. The extra-force coefficients, metric response, conditional static distances, released annular covariance and 0.8-arcsecond PSF are unchanged. The extra-source cutoff remains **20 times the old Hernquist a in physical kpc**, isolating the light-profile change instead of silently moving the capture boundary. This boundary is still a trial source assumption, not a measured capture radius.

J0330-0020 remains excluded by the release's use flag, with its underlying reason unresolved. J1538+5817 is separately excluded because detailed component parameters are unavailable. Neither is counted as a success. No validation or test score is calculated.

## Known mathematics and conditional physics

For each fixed projected Sersic profile I(R), spherical Abel inversion gives

    nu(r) = -1/pi integral_r^infinity I'(R)/sqrt(R^2-r^2) dR
          = -1/pi integral_0^infinity I'(r cosh u) du.

This is a known geometric inversion, not a photon/graviton equation. We compute the cumulative ordinary-matter fraction from 4 pi integral nu r^2 dr and obtain g_b=GM(<r)/r^2. A tiny finite-quadrature discrepancy in total luminosity is normalized consistently and recorded. The source integration truncates at an exponentially faint tail; its numerical limits are tested, not treated as a physical stellar edge.

The previous constant-beta Jeans equation and seeing-convolved annular second-moment projection then predict V_rms. The empirical g_c=A a_star(g_b/a_star)^p remains an unproven photon interpretation. Under the same equal-potential weak-field assumption, the total force predicts a spherical Einstein angle. We use the full released Gaussian covariance in V_rms; lensing is excluded from mass/orbit fitting.

A free mass scale derived from stellar motions is not an independent stellar-population mass. The baseline can absorb additional gravity; neither branch here establishes how much ordinary matter is really present. A positive Jeans second moment is also not proof of a globally nonnegative, stable stellar distribution function.

## Numerical evidence

The Abel inverse reproduces an analytic Gaussian density to 4.45e-16 relative error. Doubling deprojection order on a two-component profile changes meaningful densities by less than 2.3e-14. Its integrated spherical luminosity agrees with the analytic projected luminosity to 5.53e-6. A separate incomplete-gamma projected-mass calculation verifies the baryonic Sersic Einstein root to better than 1e-4.

For every fitted galaxy and branch, radial nodes, angular nodes and deprojection order are doubled at the fitted parameters. Maximum fractional V_rms change is below 3.23e-5. These checks rule out a large integration-resolution explanation for the remaining residuals; they do not validate sphericity, isotropy/constant anisotropy, seeing, source geometry or the companion mechanism.

## Interpretation and next work

The empirical extra-force shape continues to reduce radial discrepancies within this restricted model family. Remaining radial residuals, orbit-boundary fits and mixed lensing predictions prevent claiming a unified solution. Next work should quantify tangential-orbit freedom and image-shape uncertainty, constrain stellar population masses and replace the empirical extra source with calculated capture. Do not tune the published light distribution or assign a separate lensing gain to force agreement.

All six scientific goals remain open, including the source-energy, redshift/timing, Milky Way and cluster requirements outside this pilot. Code, protocol, per-bin predictions, optimizer starts, exclusions and source hashes are archived. Reproduce with verify.py and run.py in research_work/results/slacs-component-refit.
