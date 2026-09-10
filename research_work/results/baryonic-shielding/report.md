# Can baryonic shielding create the required uneven deposit pattern?

**The tested shielding rule improves the match to the target source shape but does not reproduce it.** Its exponentiated RMS logarithmic residual falls from 6.01 for constant loading to 4.49 with shielding. This is an equal-weight comparison of model values at 240 grid points, not a measured uncertainty, stellar likelihood or independent observational success.

This follows the constant-loading audit by giving the incoming companions a specified spatial transport law. It retains one common absorption coefficient and one common normalization; it does not assign a separate loading to each location.

## Hypothesis and formula provenance

Assume a stationary, isotropic companion intensity at a spherical boundary, straight propagation, no internal emission or scattering, and absorption coefficient per unit length `beta(x)=kappa*rho_b(x)`. Deposits remain attached to their local ordinary matter, and all locations share the same exposure duration. Ordinary matter is held fixed. These are optional physical assumptions for this test, not consequences of the user's general deposit concept.

Known absorption transport gives:

`Sigma(x,n) = integral_from_x_to_boundary rho_b(x-s*n) ds`
`I(x,n) = I_boundary * exp[-kappa*Sigma(x,n)]`
`q_model(x) = C * <exp[-kappa*Sigma(x,n)]>_directions`.

Here Sigma is ordinary-matter column density, kappa is absorption cross section per ordinary mass, and q is the extra equivalent source per ordinary mass. The angular brackets average the incident directions. The exponential and angular transport are known mathematics; identifying the absorbing sector as companions and specifying beta are hypothetical. No novel fundamental equation is claimed.

C contains the accumulated exposure, capture normalization and assumed conversion from retained energy to an ordinary Newtonian source. Fitting it independently is a source-shape diagnostic; it does not measure the available photon energy or break the source-to-gravity normalization degeneracy. The constant-loading comparison is the earlier shape model, not physical capture at exactly zero absorption with a finite fixed exposure.

## Numerical setup

The ordinary density uses the same analytic bar, nuclear, stellar-disk, gas-disk and softened-central components as the preceding audit. Their sum is checked against all 240 archived ordinary-density values. Rays end at a 30-kpc sphere. A 60-kpc boundary is an explicit sensitivity test; these radii are modeling choices, not measured halo boundaries.

Angular quadrature uses 128 directions for the coarse run and 512 for the finer runs. Radial quadrature uses 64 and 128 nodes on each segment, splitting the path at the disk midplane before integration. The zero-optical-depth angular normalization is verified. The search scans kappa from 1e-12 to 1e-5 kpc squared per solar mass, then refines around the lowest scanned interval. Reported optima lie inside that scan; this is not a proof of a global optimum over every possible capture law.

For each kappa the common log-normalization is solved analytically by the mean log target minus mean log attenuation. The loss is the mean squared log residual across the selected grid. It does not use observational error weights, volume weights, a selection function or held-out stars.

| Calculation | Best kappa (kpc²/Msun) | RMS log residual | exp(RMS log residual) |
|---|---:|---:|---:|
| Constant loading | — | 1.7939 | 6.0129 |
| coarse | 2.31959e-09 | 1.5044 | 4.5016 |
| fine | 2.31998e-09 | 1.5016 | 4.4887 |
| outer60 | 2.31941e-09 | 1.5017 | 4.4892 |

Only 23.3% of sampled loadings fall within a factor of two of the target. The largest mismatch is a factor of 102.8. Those summaries depend on the declared grid and loss function; they are not rejection probabilities.

## Resolution and boundary checks

Holding the coarse fitted parameters fixed, finer quadrature changes predicted log loadings by at most 0.01830, roughly 1.85% multiplicatively. Refitting barely changes kappa or the overall residual. Moving the boundary to 60 kpc leaves the fitted residual factor at approximately 4.489. Thus the tested resolution and boundary changes are much smaller than the remaining factor-level mismatch. More refinement would be necessary for precision capture-rate inference; no such precision claim is made.

The target effective source itself has up to 5.4% fine/finer density changes in the previous audit. Ordinary-matter uncertainties, more general outer illumination and ray bending are not covered by these numerical checks.

## What remains physically unresolved

Shielding can reduce the flux per unit ordinary matter in the interior, giving the desired direction of uneven loading. In this specified geometry, a single absorption coefficient does not supply the full required spatial pattern. Merely adjusting the global exposure cannot fix that shape mismatch.

A next physical candidate would need an independently specified change, such as capture dependence on the gravitational environment, internally generated companions with a modeled stellar-emission history, a supported reservoir separate from ordinary matter, or a different gravitational response. Giving each region an arbitrary cross section would make agreement a definition rather than a prediction.

This stationary snapshot also omits deposition-induced motion, recoil, drag and changing baryonic mass. It cannot serve as a self-consistent history when accumulated loading is substantial. It does not derive the boundary bath from photon redshift, test total photon supply, or show that ordinary clocks and event timing work. The broader theory remains unvalidated; no holdout has been opened.

## Reproduction

Run `run.py`, then `report.py`. Input hashes, both quadrature levels, the boundary sensitivity, fitted parameters and frozen-parameter comparisons are retained. `predictions.json` records the required and predicted loading at every probe. Large column-density arrays remain in the ignored cache. The existing target potential, source and stellar catalogs are unchanged.
