# Directional moments of catalog-fed deposition

## Result

Keeping source directions in the unchanged catalog pilot produces a noncentral deposit-power distribution. THEMIS gives a three-dimensional centroid offset 74.25 kpc from NGC4486 and an orthographic sky-plane offset 38.91 kpc; DL14 gives 73.70 and 38.65 kpc. The projected second-moment minor/major RMS ratio is 0.9041 and 0.9046 respectively. These describe where power is delivered in this prescribed scenario, not measured cluster structure or a lensing peak.

## Assumptions and provenance

All 814 usable source positions, distances and SED luminosities, the 1 Mpc receiver radius, archived empirical alpha and hypothetical kappa R=10 are inherited unchanged from cluster-catalog-pilot. No capture coefficients or density profiles are refitted. The catalog's luminosities remain model-dependent and its completeness and source histories remain unknown. Point emitters, constant opacity, straight rays, no intervening capture and a stationary receiving sphere remain conditional assumptions, not first-principles physics. Full energy-momentum/gravitational backreaction and support are absent.

Known geometry and moment identities are applied to the conditional transfer solution. The deposited power element on a ray is proportional to kappa C(s) ds. Relative to the source-center axis u, let z be the longitudinal coordinate and b its perpendicular distance. Azimuthal symmetry for each individual isotropic source permits exact azimuth integration:

    first moment = integral z dP * u
    second moment = integral (b^2/2) dP * Identity
                    + integral (z^2-b^2/2) dP * outer(u,u).

Summing sources gives centroid X=integral x dP / P and covariance Cov=integral xx^T dP/P-XX^T. These are established moment formulas, not uniquely proposed laws. Unlike radial averaging, they retain actual source directions. They do not uniquely reconstruct the full spatial distribution.

The sky basis is east/north perpendicular to the line from Earth to the catalog center. Projection is orthographic in physical coordinates, not an exact finite-distance angular lens calculation. The eigenvalues of the projected covariance yield the quoted RMS axis ratio; this is not an isophotal ellipticity or fitted mass-profile axis ratio.

## Verification

Angular and along-ray integration use 64 and 128 Gaussian nodes each. Total power agrees with the independent completed-ray catalog ledger within 4.7e-14 relative. Centroid refinement changes are below 1.7e-12 Mpc and covariance changes below 1.3e-12 Mpc^2; these are numerical sensitivity estimates, not observational precision. Covariances have nonnegative eigenvalues. A central isotropic source has zero first moment and equal longitudinal/per-axis transverse second moments, explicitly checked. Both catalog model versions and hashes are recorded in results.json.

## Why this matters and what it does not establish

The current all-direction supply is not spatially uniform because the contributing sources are not distributed uniformly. Thus a spherical capture boundary does not imply a centered spherical deposit field. The offset and shape emerge without fitting a lens map. Their astronomical accuracy is untested and may change with distance uncertainties, missing sources, luminosity corrections and historical evolution.

A centroid is not the position of maximum density, maximum shear or a reconstructed lensing peak. Calling this a prediction of a 39 kpc lensing offset would be incorrect. With stationary geometry and constant illumination, deposits with identical retention histories would inherit these power moments; evolving systems need explicit history and transport. The next step toward a lens test is a spatially resolved source, a specified response of both metric potentials, and forward prediction of image distortions with uncertainties. No final observational holdouts were used and all six physical goals remain open.

Run `python research_work/results/cluster-deposit-directionality/run.py`.
