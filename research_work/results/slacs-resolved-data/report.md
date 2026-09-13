# Resolved stellar kinematics acquired for the lensing comparison

## What is now local

We archived 28 KCWI radial data products (one velocity profile and one covariance matrix for each of 14 SLACS lenses), plus the release documentation. Every downloaded byte count and Git blob hash matches the pinned source tree; SHA-256 hashes and cache locations are recorded in manifest.json.

Primary release: [TDCOSMO2025_public, SLACS data](https://github.com/TDCOSMO/TDCOSMO2025_public/tree/d7f38db341f68be1df0d9ac1fc528c45113f94cf/ExternalLenses/SLACS), commit d7f38db341f68be1df0d9ac1fc528c45113f94cf. The observational study is [Knabel et al., Spatially Resolved Kinematics of SLACS Lens Galaxies I](https://arxiv.org/abs/2409.10631). The public products used here are azimuthally rebinned radial profiles, not full two-dimensional maps or detector cubes.

## Existing roles preserved

| Existing role | Systems in this release |
|---|---:|
| Training | 8 |
| Validation | 3 |
| Test | 2 |
| Not assigned in our original modeled-lens table | 1 |

The eight training systems are J0037-0942, J0330-0020, J1112+0826, J1204+0358, J1402+6321, J1538+5817, J1621+3931 and J1630+4520. These provide a subset of the existing 33-system training pilot for the next resolved comparison. The unassigned object stays unassigned; it is not silently placed into training or claimed as an untouched prediction.

No profile fit or new lensing score was computed in this acquisition. The source metadata and example J0037-0942 profile have been inspected, and a published integrated-kinematics table was read during discovery. Therefore we make no new claim of pristine observational blindness. Reserved role assignments and uncomputed evaluation scores remain preserved.

## What the products measure

The release identifies the radial columns as inner/outer angular bin edges in arcseconds and binned V_rms in km/s. This combines mean streaming motion and random velocity dispersion. The accompanying covariance includes correlated uncertainties according to the documentation; it must not be replaced by independent per-bin errors.

The README's displayed shell-averaging formula appears dimensionally inconsistent as written. We will verify the source preprocessing before reproducing its binning. The intended standard local second moment is V_rms^2=V^2+sigma^2, with luminosity-weighted second moments combined before taking a square root. That identity is known kinematics, not a new companion postulate. Published released values are preserved without modification here.

The radial profile and a dispersion extracted from one integrated spectrum are different observables. Neither should be substituted for the other. For rotating galaxies, our existing spherical nonrotating Jeans model is an explicit approximation; adopting V_rms does not automatically validate spherical symmetry or the orbital assumptions.

## Why this advances the joint test

The earlier pilot used one aperture dispersion to set one mass scale. Several radial measurements can constrain how the potential and orbital projection change with radius, reducing the freedom that moved the median lensing discrepancy from 3.7% to 15.5%. They do not automatically remove mass-orbit degeneracy, especially after azimuthal averaging.

Before fitting, we must check covariance dimensions and numerical validity, angular bins and PSF conventions, spatial light weights and centering, and how the released second moments were formed. Then predict those same seeing-convolved annular measurements under our common force rule. Infer mass/orbit parameters using training kinematics and predict lensing without selecting a separate lensing gain. A full capture-derived source and independent stellar-mass constraints remain necessary for the physical interpretation.

No cosmological posterior, dark-halo mass fit, external-convergence distribution or precomputed dynamical scaling is imported as a physical fact. The paper's cosmography motivation does not prevent using the measured stellar spectra and angular binning in our hypothetical universe, provided preprocessing assumptions are audited.

## Reproduction and scope

Run `python research_work/results/slacs-resolved-data/acquire.py`. It downloads only the explicit documentation and radial CSV products into the ignored generated-data directory, verifies immutable source hashes, and writes the tracked manifest. It neither executes release code nor unpickles release likelihood objects.

This is completed data acquisition and identifier matching, not covariance validation, a resolved fit, an anisotropy measurement or completion of any of the six scientific demonstrations.
