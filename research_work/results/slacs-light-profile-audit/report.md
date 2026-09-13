# Published stellar-light profiles: input correction before gravity refitting

## Finding

Six of the seven resolved-fit training galaxies have single- or double-Sersic image-component parameters in the pinned release. Their reconstructed half-light radii reproduce the tabulated modern radii to within 0.029%. Relative to the older sizes used in our Hernquist pilot, the largest change is +52.23% for J0037-0942. This is a material source-profile difference, not a new force parameter fitted to lensing.

| Galaxy | Previous half-light radius (arcsec) | Reconstructed published profile (arcsec) | Change |
|---|---:|---:|---:|
| J0037-0942 | 2.190 | 3.334 | +52.23% |
| J1112+0826 | 1.500 | 1.623 | +8.22% |
| J1204+0358 | 1.470 | 1.544 | +5.06% |
| J1402+6321 | 2.700 | 2.721 | +0.77% |
| J1538+5817 | 1.580 | Missing | Not calculated |
| J1621+3931 | 2.140 | 2.402 | +12.24% |
| J1630+4520 | 1.960 | 2.012 | +2.64% |

J1538+5817 lacks component parameters in this release. It stays missing; its older model can remain explicitly labeled for comparison but must not be represented as a newly recovered detailed image profile. We do not substitute a fallback radius or invent a second component.

## Provenance and convention check

Inputs come from the pinned [TDCOSMO2025 SLACS metadata](https://github.com/TDCOSMO/TDCOSMO2025_public/blob/d7f38db341f68be1df0d9ac1fc528c45113f94cf/ExternalLenses/SLACS/slacs_all_params.csv), verified against the prior acquisition SHA-256. Only the existing seven training systems are evaluated. No lensing residual enters component selection or radius calculation.

The release preprocessing passes these amplitudes, Sersic radii and indices into its light model. We interpret amplitude as surface brightness at the component radius and use the documented lenstronomy product-average/equal-area radius convention. [The official profile documentation](https://lenstronomy.readthedocs.io/en/latest/_modules/lenstronomy/LightModel/Profiles/sersic.html) and [Sersic utilities](https://lenstronomy.readthedocs.io/en/latest/_modules/lenstronomy/LensModel/Profiles/sersic_utils.html) specify these conventions and the approximation b_n=1.9992n-0.3271. The historical fitting package version is not verified; the agreement with the independently tabulated combined radii is a practical consistency check, not full reproduction of the original imaging fit.

All double components have the same ellipticity within each system. Under the product-average convention, ellipse area inside radius R is pi R^2. The component light fractions can therefore be combined using the spherical radial expression without an extra axis-ratio factor. That equal-area accounting does not make the physical galaxies spherical.

## Known formulas used

These are established Sersic and incomplete-gamma identities, not new companion/time laws:

    I_j(R) = amp_j exp[-b_j ((R/R_j)^(1/n_j)-1)]
    L_j = 2 pi amp_j R_j^2 n_j exp(b_j) Gamma(2n_j) / b_j^(2n_j)
    L_j(<R) / L_j = P(2n_j, b_j (R/R_j)^(1/n_j)).

P is the regularized lower incomplete gamma function. The combined half-light radius solves sum L_j(<R) = 0.5 sum L_j. Relative amplitudes are retained; arbitrary luminosity units cancel in the normalized fractions. No absolute stellar mass follows from that normalization. All component parameters, ellipticities, light fractions and reconstruction differences are stored in results.json.

The b_n approximation means the component parameter is only approximately an exact mathematical half-light radius, particularly for n=1. We keep the documented approximation for this interpretation audit rather than silently replacing it with an exact gamma-inverse constant.

## Consequence for the last fit

The previous chi-square improvement used an approximate Hernquist profile and older radii. It cannot be interpreted independently of this light-profile mismatch. A larger or more extended stellar distribution changes both the tracer weighting and the baryonic acceleration. Merely rescaling the reported mass is insufficient; the radial Jeans prediction must be recomputed using the same updated profile that supplies the gravitational mass distribution.

The next comparison should deproject the published components under an explicit geometry assumption, refit only the radial kinematics, and predict lensing again. Preserve the prior extra-source cutoff in physical units when isolating the light-profile change, so changing a size does not silently change a capture boundary. Missing J1538+5817, component covariance, stellar population mass-to-light ratios, dust/bandpass differences and three-dimensional shape remain uncertainties. Better fits are not guaranteed and must not be obtained by adjusting the image profiles to the desired gravity result.

This completes an observational input audit, not a new dynamical fit or validation of the companion mechanism. All six objectives remain active.

Reproduce with `python research_work/results/slacs-light-profile-audit/run.py`.
