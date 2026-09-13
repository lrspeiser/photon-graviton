# KCWI radial input and covariance audit

## Result

The eight existing training-role lenses contain 54 radial V_rms measurements. All eight covariance matrices have the expected dimensions, finite entries, numerical symmetry, strictly positive eigenvalues and successful Cholesky factors. Radial bins start at zero and are contiguous with strictly positive widths. Every measured V_rms is positive. This establishes usable numerical inputs, not a validated gravity model.

The largest inter-bin correlation is 0.8433 and the largest matrix condition number is 36.95. Using diagonal errors alone would discard material covariance. The metadata specify Gaussian PSF FWHM=0.8 arcseconds for all eight profiles, rather than the 1.4-arcsecond SDSS seeing approximation in our earlier integrated-aperture pilot. We retain the released common value without calling it a separately measured PSF for each galaxy.

## Data roles and quality flag

The source release marks J0330-0020 with flag_ifu_kcwi=0. The other seven training systems have flag 1. We retain J0330-0020 in this audit and do not silently discard it or include it as an approved fit input. The flag's scientific reason must be understood before deciding its role in the resolved fit. Positivity of a covariance matrix does not override an observational quality flag.

No validation or test profile has been fitted or numerically scored. Metadata files contain the wider sample, but numerical array checks and the analysis-ready JSON here are restricted to the eight pre-existing training systems.

## Verified conventions

The pinned [TDCOSMO preprocessing notebook](https://github.com/TDCOSMO/TDCOSMO2025_public/blob/d7f38db341f68be1df0d9ac1fc528c45113f94cf/ExternalLenses/SLACS/kinematic_sample_slacs_preprocessing.ipynb) directly loads binned_Vrms and its covariance. It constructs centered IFU_shells apertures from the angular bin edges and passes psf_ifu_kcwi as a Gaussian FWHM. It does not square the released velocities or covariance before its likelihood interface. Our arrays preserve those conventions: arcseconds, km/s, and (km/s)^2.

The covariance-file header says the radial covariance was calculated from 200 draws of the Voronoi-bin covariance. Its finite estimation uncertainty is not assessed by the positive-definiteness check. We do not add a Hartlap factor or other correction without establishing the correct sampling assumptions. We also do not add an arbitrary extra diagonal error or count the documented systematic contribution twice.

The inspected notebook loads the released radial products; it does not regenerate the original Voronoi-to-radial averaging. Thus the dimensionally inconsistent averaging expression in the README remains a documentation issue. We have not used that expression to alter the data. Reconstructing the original aggregation requires additional source information.

## What was acquired and excluded

We additionally archived and verified the pinned preprocessing notebook and slacs_all_params.csv against source Git blob hashes and byte counts. Their SHA-256 hashes are in results.json. Only observational seeing and use-flag fields are imported from the metadata into the training arrays. The notebook was read as JSON text, never executed. Its fitted lens slopes, cosmology, external-convergence distributions and posterior objects are not adopted as facts in our model.

The original 29-file acquisition hashes are rechecked before the training arrays are read. results.json retains all bin edges, velocities, covariance matrices, PSF entries, flags and numerical diagnostics. The original input files remain unchanged.

## Next calculation

For the eligible training profiles, predict seeing-convolved annular second moments with the same light profile, orbit assumptions and total force used for lensing. Use the full multivariate Gaussian residual likelihood with its covariance; do not reuse cumulative aperture dispersions as if they were annular values. Infer mass and orbital behavior from kinematics, then compare the lens angle without a separate lensing gain.

Before making a substantive conclusion, test radial numerical convergence, investigate the flagged object and propagate light-profile and seeing uncertainties. A spherical Jeans solution remains a conditional approximation, and the empirical extra-force formula remains short of a capture-derived source. All six scientific objectives remain open.

Reproduce with `python research_work/results/slacs-resolved-input-audit/run.py`. This audit downloads only the two pinned supporting files, checks local training arrays and writes its own result.
