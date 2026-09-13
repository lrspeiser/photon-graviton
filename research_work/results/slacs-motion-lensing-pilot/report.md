# SLACS stellar-motion to lensing transfer pilot

## Outcome

We predicted lensing for 33 training-role early-type galaxies after fixing each mass normalization from its spectroscopic stellar dispersion. No Einstein angle was used to fit those masses or to change the previously frozen extra-gravity coefficients. The empirical extra-force template does not improve this descriptive comparison under the stated assumptions.

| Force model | Outer extra-source cutoff | Median predicted / catalog angle | Median absolute fractional discrepancy | RMS fractional discrepancy |
|---|---:|---:|---:|---:|
| Mass follows light, free mass normalization | Not applicable | 1.017 | 12.55% | 16.41% |
| Frozen empirical extra force | 20a | 1.092 | 13.70% | 19.51% |
| Frozen empirical extra force | 100a | 1.099 | 13.40% | 19.89% |

Every selected system has a numerical Einstein root. The baseline is duplicated across the two cutoff settings in the machine-readable output as a consistency check; those duplicates are not independent observations. This is an exploratory training result, not a formal significance ranking or a held-out success. Neither validation nor test systems are scored.

**Crucial interpretation:** the baseline mass normalization is fitted from gravity-sensitive motions, not measured independently from a stellar population. It can absorb missing gravity. Its lensing agreement therefore does not demonstrate that ordinary matter alone suffices. Likewise, the extra-force branch uses the empirical rotation relation, not calculated photon supply and capture. The comparison tests transfer of a radial force shape, not the cause of the extra gravity.

## Observational inputs

[Bolton et al. (2008), SLACS V](https://arxiv.org/html/0805.1931v1), Table 4/5 and section VI.3, provide the archived image-model sizes, dispersions and SIE lens summaries. We use its 3-arcsecond-diameter spectroscopic aperture and approximate the quoted median 1.4-arcsecond seeing by a circular Gaussian. Actual per-object seeing is not supplied here. Table 4 radii use the intermediate axis. We retain all training early-type systems with measured positive dispersion, including catalog dispersion-quality flags; we do not select systems using lensing residuals.

The existing static-Euclidean conditional distances are retained, including their extrapolation of the fitted logarithmic distance-redshift rule and omission of peculiar/endpoint corrections. They are not independent geometric distances. SIE Einstein angles are image-model summaries with ellipsoidal assumptions, not raw image positions or exact spherical critical radii. Their full covariance is unavailable in this table. Consequently no chi-square probability or formal model preference is assigned.

## Formula provenance and assumptions

**Known stellar model:** use the spherical [Hernquist (1990) profile](https://doi.org/10.1086/168845), with tracer density nu proportional to a/[r(r+a)^3], mass-following-light acceleration g_b=GM/(r+a)^2 and a=R_e/1.8153. It is an approximation to the catalog de Vaucouleurs light model. Sphericity, isotropic equilibrium orbits, constant mass-to-light ratio, no separate gas or black hole, and no environmental field are explicit assumptions.

**Existing empirical project fit, not a derived deposition law:**

    g_c = A a_star (g_b/a_star)^p

The exact archived A, p and a_star are read from joint-galaxy-audit/results.json and copied into this result. Beyond r_t=20a or 100a, g_c(r)=g_c(r_t)(r_t/r)^2. These two cutoff choices are declared sensitivity cases, not measured boundaries or fitted values. Stellar motions and light use the same total g_b+g_c and the same cutoff.

**Known isotropic Jeans relation:**

    nu(r) sigma_r^2(r) = integral_r^infinity nu(s) g(s) ds.

The aperture prediction weights that pressure and tracer density over the line of sight and the seeing-convolved circular fiber. For an isotropic system, the local line-of-sight variance equals sigma_r^2. For projected radius R, the Gaussian probability of entering the fiber is the noncentral chi-square CDF with two degrees of freedom, threshold (R_ap/s_psf)^2 and noncentrality (R/s_psf)^2. Averaging this probability over each spherical shell yields W(r). Then

    sigma_ap^2 = integral r^2 nu(r) sigma_r^2(r) W(r) dr
                 / integral r^2 nu(r) W(r) dr.

These are established projection and probability identities, not new gravity formulas. For fixed size and cutoff, sigma_ap^2=C_b m+C_c m^p, with m=M/(10^11 solar masses). We solve this monotonically for mass using only the observed dispersion. The catalog fitted spectral Gaussian width is approximated by this model second moment.

**Known weak-field lensing, conditional equal-potential response:**

    alpha_hat(b) = 4/c^2 integral_0^infinity g(sqrt(b^2+z^2)) b/sqrt(b^2+z^2) dz
    theta_E = (D_ls/D_s) alpha_hat(D_l theta_E).

The equality of temporal and spatial potentials is a postulate here, not a companion-field derivation. There is no independent lensing gain. Each row contains the predicted angle, fitted mass and the changes obtained from measured dispersion plus/minus its quoted error. Those endpoints propagate **only** dispersion uncertainty and are not full prediction intervals.

## Numerical checks and reproducibility

The global isotropic Hernquist aperture reproduces the scalar-virial value GM/(18a) to relative error 7.53e-5. Doubling radial and angular resolution for a representative finite aperture changes its two force coefficients by at most 3.14e-5. A separate analytic projected Hernquist mass verifies the baryonic Einstein root to better than 1e-7. These checks concern numerical implementation, not adequacy of spherical geometry for flattened real lenses or convergence of every nuisance choice.

Input and code hashes are recorded in results.json; protocol.json records the pre-calculation choices. Run verify.py and run.py in this directory from the repository root using Python. No dark-halo masses, cosmology-corrected luminosities or observed lens angles enter mass normalization.

## Next scientific requirement

Constrain stellar masses independently and replace or marginalize the spherical/isotropic, shared-seeing and cutoff assumptions. Add the image likelihood and capture-derived spatial source. Only then can the same-source motion/lensing demonstration be claimed. The current residuals should guide training work; reserve validation and test scores until those choices are frozen. All six research objectives remain open.
