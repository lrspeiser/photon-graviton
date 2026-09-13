# Shared supernova optical geometry does not yet reconcile lensing

The previous turn improved supernova brightness with an observer-regular beam-area postulate. This run uses that same optical law to change lens distances, image/aperture physical scales and lens-source propagation consistently. The supernova and SPARC capture parameters are frozen; ordinary stellar masses and stellar anisotropies alone are refitted to inner stellar bins.

## Lens-to-source geometry from the same beam equation

Let y(lambda)=D_A be the observer's angular diameter distance from the regular-area law, lambda=f/alpha and f=z/(1+z). A conditional isotropic optical background has scalar Jacobi equation y''=-R_opt y. A solution with a vertex at the lens follows by known reduction of order:

    D_ls=(1+z_l)y_l y_s integral_lambda_l^lambda_s d lambda/y(lambda)^2,
    D_ls/D_s=(1+z_l)y_l integral_lambda_l^lambda_s d lambda/y(lambda)^2.

The factor 1+z_l normalizes the beam's opening angle using photon energy at the lens. Omitting it would use the observer's affine normalization as if it were local to the lens. D_ls is not generally D_s-D_l.

This is established second-order ODE mathematics applied to a postulated optical background, not a new field equation or proof of a global nonexpanding spacetime. Background shear is omitted and thin localized lens perturbations use the same ordinary-plus-deposit weak-field bending as before. No law producing the background, companion transfer and lens perturbations together has yet been derived.

The actual lens/source redshifts determine conditional path lengths via the unchanged alpha. These are distinguished from the revised angular diameter distances. All angular light profiles, PSF/apertures and the equivalent-exponential morphology mapping use the revised angular distance consistently. Independently stipulated catalog distances are not changed. The six lens geometries were previously redshift-inferred, not independently measured.

## Predicted lens angles

Interception branch, with no lensing data in the inner-stellar fit:

| Galaxy | Catalog SIE angle | Earlier geometry prediction | Shared optical-law prediction |
|---|---:|---:|---:|
| J0037-0942 | 1.530 | 1.552 | 1.797 |
| J1112+0826 | 1.490 | 1.210 | 1.453 |
| J1204+0358 | 1.310 | 1.333 | 1.514 |
| J1402+6321 | 1.350 | 1.408 | 1.630 |
| J1621+3931 | 1.290 | 1.150 | 1.368 |
| J1630+4520 | 1.780 | 1.415 | 1.717 |

Angles are in arcseconds. Several earlier deficits improve, but others become excesses. Catalog SIE radii remain model summaries rather than exact spherical rings; no full lens likelihood or uncertainty propagation is supplied.

| Diagnostic | Earlier interception geometry | Shared optical-law geometry |
|---|---:|---:|
| Inner stellar chi-square sum | 47.614 | 48.917 |
| Outer conditional residual-square sum | 57.978 | 63.039 |
| Lens fractional RMS | 12.352% | 13.112% |
| Anisotropy-boundary fits | 0 | 0 |

The transparent comparison under the revised geometry has inner sum 41.723, outer sum 48.150 and lens RMS 13.810%. All fits use the original bounds, starts and covariance. Scores are descriptive plug-in results on exposed systems, not independent significance statements.

## Assessment

The successful supernova brightness correction cannot be treated as unrelated to lensing. Once the distance factors are linked through the same conditional beam equation, the combined capture model still fails to deliver convincing joint agreement. This result does not prove that the optical law or every capture model is impossible; it exposes the specific consequences of combining these two current candidates.

The equivalent-disk mapping for early-type galaxies, spherical source, ordinary stellar mass assumptions and unproved particle support remain material limitations. Stellar supply, a common global spacetime/interaction and independent source distances also remain missing. No parameter was tuned to make the improved individual lens predictions look like a global success. All six goals remain open.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/lensing.py --regular-optics`, then `python research_work/results/isotropic-galaxy-transfer/regular-optics-check.py`. The latter solves the lens-vertex Jacobi initial-value problem independently; its distance ratios agree with the integral solution within 8.4e-12 relative. This verifies the chosen optical completion, not its astronomical truth. Old lensing results are preserved; revised geometries, all stellar predictions, nuisance fits and input hashes are in regular-optics-results.json.
