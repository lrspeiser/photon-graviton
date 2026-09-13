# Stars plus NFW: same lenses, two geometries

Changing the geometry materially improves the joint stellar-motion/lens-angle fit for NFW in this sample. NFW also struggles in the problematic systems under our retained optical branch, so those mismatches cannot be attributed solely to the fixed companion deposit profile. This is not proof that standard geometry is uniquely correct or that NFW completely fits all six galaxies.

## Matched comparison

All rows below impose the catalog lens angle exactly, and score all stellar bins with the same released covariance, light-profile shape, PSF/apertures, stellar mass bounds and orbital beta bounds. Lower is better. There is no held-out lens prediction or calibrated goodness-of-fit probability.

| Galaxy | Flexible companion, our geometry (Chabrier) | NFW, our geometry | NFW, standard benchmark |
|---|---:|---:|---:|
| J0037-0942 | 177.04 | 88.59 | 20.20 |
| J1112+0826 | 30.08 | 33.92 | 11.63 |
| J1204+0358 | 131.27 | 54.38 | 4.05 |
| J1402+6321 | 251.73 | 195.85 | 43.07 |
| J1621+3931 | 3.83 | 5.24 | 4.76 |
| J1630+4520 | 7.90 | 11.17 | 1.57 |

NFW has target-fitted halo amplitude and scale, whereas the companion mixture preserves its one-third inventory but fits target-specific redistribution weights. These are compatibility diagnostics with different freedoms, not an equal-parameter forecast contest. The NFW stellar mass is not tied to a Chabrier/Salpeter mass prior, so its fit is population-independent under this broad allowed range.

The standard benchmark improves all six NFW scores. J0037, J1204 and J1402 show especially large reductions, while J1402 still has appreciable residuals. Replacing companion deposits with NFW under our geometry does not eliminate the original pattern of tension. Thus our specific optical branch needs further testing; this is not a rejection of every nonexpanding model.

## Why distances change the comparison

The known spherical lens relation is theta_E=(D_ls/D_s)*alpha_hat(D_l*theta_E). Stellar dynamics and the image aperture also depend on the physical length D_l*theta. The same apparent image and velocities therefore constrain a different mass distribution when the optical distance ratios change. We rebuild every physical light radius, aperture and PSF scale under each geometry; we do not just rescale the final angle.

| Galaxy | Our D_ls/D_s | Standard D_ls/D_s | Standard NFW largest motion residual / plotted error |
|---|---:|---:|---:|
| J0037-0942 | 0.7747 | 0.6532 | 3.54 |
| J1112+0826 | 0.6472 | 0.5234 | 1.40 |
| J1204+0358 | 0.8203 | 0.7055 | 0.72 |
| J1402+6321 | 0.6320 | 0.5433 | 6.10 |
| J1621+3931 | 0.6710 | 0.5531 | 0.97 |
| J1630+4520 | 0.8005 | 0.6388 | 1.09 |

The standard alternative is flat FLRW with H0=70 km/s/Mpc, Omega_m=0.3 and Omega_Lambda=0.7, with radiation neglected. Its distances use the known comoving integral and D_ls/D_s=1-chi_l/chi_s. These are comparison assumptions only, not inserted into our fictional universe. Our branch retains its prior fitted alpha and regular optical response. The old conditional-geometry.json is static Euclidean and was not mislabeled as standard cosmology.

## Parameter limitations

| Galaxy | Geometry | Halo fraction of deflection | r_s / R_e | Beta | Bound flags |
|---|---|---:|---:|---:|---|
| J0037-0942 | companion_regular | 0.5190 | 0.0277 | -2.0000 | beta |
| J0037-0942 | standard_flat_FLRW | 0.1141 | 1.3683 | -0.1990 | none |
| J1112+0826 | companion_regular | 0.2409 | 3.2058 | -0.0033 | none |
| J1112+0826 | standard_flat_FLRW | 0.5516 | 9.4444 | 0.4500 | beta |
| J1204+0358 | companion_regular | 1.0000 | 0.0678 | -2.0000 | beta, fraction |
| J1204+0358 | standard_flat_FLRW | 0.2446 | 0.0578 | -2.0000 | beta |
| J1402+6321 | companion_regular | 0.9543 | 0.0512 | -2.0000 | beta |
| J1402+6321 | standard_flat_FLRW | 0.2127 | 0.0122 | -2.0000 | beta |
| J1621+3931 | companion_regular | 0.3312 | 1.5047 | 0.4500 | beta |
| J1621+3931 | standard_flat_FLRW | 0.2762 | 100.0000 | 0.0084 | scale |
| J1630+4520 | companion_regular | 0.2371 | 2.9526 | 0.2586 | none |
| J1630+4520 | standard_flat_FLRW | 0.4577 | 37.3172 | 0.4500 | beta |

Several fits require compact halos or orbital/scale/fraction boundaries. In particular the standard J1621 fit reaches r_s/R_e=100, so the scale is not a precise inference. Very small NFW scale radii and almost halo-only solutions should not be presented as plausible cold-dark-matter halos without additional stellar-population and halo-population priors. No such priors were imposed here, to make the flexible mass-profile comparison explicit. The halo fraction is a projected bending fraction, not a three-dimensional total mass fraction.

NFW is the known density profile rho proportional to 1/[x(1+x)^2]. No new dark-matter formula is proposed. It is applied without a cosmological critical-density normalization; the amplitude is fitted directly. Positive mass and shared gravitational response are assumed for both lensing and dynamics. An untruncated NFW is the conventional local fitting profile here, not a finite total halo inventory comparable to the conserved companion reservoir.

## What this resolves and what remains

The same-data comparison has now been executed. It shows that optical geometry is a material part of our current mismatch, not that adding dark matter automatically cures it. It does not yet isolate the stellar/orbit approximation or propagate source-image, mass-to-light, distance and anisotropy uncertainty. A companion run under the standard geometry has not been done, so this is not a complete two-by-two test and does not prove its one-third profile would succeed with those distances.

The next test should constrain a shared revision of our optical geometry against lensing together with the existing brightness/redshift requirements, while preserving the one-third capture law. It must not choose a separate distance multiplier for each lens. A further nonspherical/radially varying stellar-orbit model may be necessary; the NFW standard residuals show that geometry alone is not a demonstrated complete resolution. All six goals remain open.

## Verification and sources

Exact-angle closure is checked to fractional 1e-10. Analytic projected NFW mass agrees with independent force quadrature to maximum relative error 1.04e-12. Fits use direct Jeans coefficients, multiple starts and positive halo/stellar mass constraints. Full predictions, fitted parameters and boundary flags are in nfw-geometry-results.json.

Known equations: [Wright and Brainerd, NFW lensing](https://arxiv.org/abs/astro-ph/9908213); [Hogg, cosmological distance measures](https://arxiv.org/abs/astro-ph/9905116). Protocol: nfw-geometry-protocol.md. Run nfw-geometry.py then nfw-geometry-report.py.
