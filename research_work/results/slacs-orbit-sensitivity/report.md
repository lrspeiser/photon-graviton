# Stellar-orbit sensitivity of the SLACS lensing pilot

## Finding

The predicted median extra-force lensing excess changes from 15.5% to 3.7% across the declared common orbital assumptions, compared with 9.2% for isotropy. Orbital uncertainty is therefore material to the previous result. The isotropic discrepancy is not a robust rejection of the force law, while selecting the most favorable orbit assumption after seeing the result would not establish a companion explanation.

| Model | Constant beta | Median predicted / catalog angle | Median absolute fractional discrepancy | RMS fractional discrepancy |
|---|---:|---:|---:|---:|
| baryons | -0.3 | 1.0609 | 14.28% | 17.84% |
| baryons | +0.0 | 1.0166 | 12.55% | 16.41% |
| baryons | +0.3 | 0.9809 | 12.08% | 15.90% |
| empirical_extra | -0.3 | 1.1548 | 15.89% | 22.37% |
| empirical_extra | +0.0 | 1.0922 | 13.70% | 19.51% |
| empirical_extra | +0.3 | 1.0368 | 11.53% | 17.24% |

All six branches retain the same 33 training systems and have numerical Einstein roots. Under each matched orbit assumption the extra-force branch still has larger RMS discrepancy than the free-mass baseline. Its beta=+0.3 median absolute discrepancy is slightly smaller, illustrating why a single favorable summary is not sufficient to select a model. No formal statistical preference is assigned and no beta is chosen for subsequent validation from this experiment.

## What changed

Only the constant orbital anisotropy changes: beta=-0.3, 0, +0.3. We refit the mass normalization from the same measured dispersion in each case, then predict lensing. Negative beta favors tangential motion, positive beta favors radial motion, and zero is isotropic. These are declared sensitivity values, not measured orbital properties or an observational confidence interval. They are applied identically to every training galaxy, without selecting an orbit pattern from its individual lens angle.

We retain the original Hernquist approximation, typical 1.4-arcsecond Gaussian seeing, 1.5-arcsecond fiber radius, conditional static distances, equal metric potentials and frozen empirical force coefficients. The extra-source cutoff is 20a for every branch. The prior 100a cutoff comparison remains separate. Validation and test systems are not scored.

## Known formulas, not new physics

The constant-anisotropy Jeans and projection relations are established stellar-dynamics mathematics; see [Mamon and Boue (2010)](https://academic.oup.com/mnras/article/401/4/2433/1127116). Only those local dynamics relations are used, not a dark-halo interpretation or an expanding-universe distance rule.

Define beta = 1 - (sigma_theta^2+sigma_phi^2)/(2 sigma_r^2). For tracer density nu and inward acceleration magnitude g:

    d(nu sigma_r^2)/dr + (2 beta/r) nu sigma_r^2 = -nu g

    nu sigma_r^2 = r^(-2 beta) integral_r^infinity nu(s) g(s) s^(2 beta) ds.

At polar cosine u along the observer's direction, the projected local variance is sigma_r^2 [1-beta(1-u^2)]. The seeing/aperture kernel W from the previous pilot must therefore be supplemented by

    T(r) = integral_0^1 aperture_probability(r sqrt(1-u^2)) (1-u^2) du.

The numerator of the aperture second moment weights pressure by W-beta T, while its luminosity denominator still weights nu by W. This is important: merely changing radial pressure while keeping the isotropic projection would give the wrong answer. The code evaluates both changes.

The photon interpretation of the extra-force formula remains unproven. These orbital relations neither create a capture law nor establish a nonnegative full phase-space distribution for the combined potential. Solving Jeans moments is a necessary equilibrium calculation, not a proof of a stable distribution function.

## Numerical evidence

The beta=0 predictions reproduce the preceding pilot to better than 1e-10 relative difference. For all three beta choices, the full-aperture baryonic calculation recovers the same scalar-virial mean line-of-sight variance GM/(18a), with relative error below 1.5e-4. That invariant checks both the pressure integration and anisotropic projection. Original finite-aperture convergence checks remain in the pilot; this run does not establish convergence for every geometry or all anisotropy laws.

All 198 predictions, dynamically inferred masses, source hashes and numerical checks are in results.json. Run `python research_work/results/slacs-orbit-sensitivity/run.py` from the repository root. Protocol choices were recorded before executing the comparison.

## Consequence for the research program

The same central spectroscopic dispersion permits different mass normalizations and lens predictions when the orbit distribution is unknown. We should constrain this freedom with spatially resolved stellar kinematics and independently inferred stellar masses before interpreting the lens residual as a gravity-law failure or success. Further work also needs realistic shape, source capture and image likelihoods.

The mass-following-light baseline is normalized from dynamics and can absorb missing gravity. Its performance does not prove independently measured baryons suffice. The extra-force branch is still an empirical transfer test, not a calculated companion-deposition model. All six scientific demonstrations remain open.
