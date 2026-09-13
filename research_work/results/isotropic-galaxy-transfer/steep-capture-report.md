# Faster-declining capture worsens frozen galaxy predictions

**Decision: do not adopt the exponent-3 capture profile.** Fitting its shared parameters to 89 galaxies worsens training and both subsequent comparison groups relative to the retained exponent-2 formula. This is evidence against this particular modification, not against all concentrated deposits or all companion mechanisms. All six goals remain open.

## Formula

The tested postulate is

    kappa(r)=k0/[1+(r/a)^2]^3, a=s*R_disk
    rho_deposit(r)=C*J(r)/[1+(r/a)^2]^3
    C=(k0/c) integral u_external dt.

The retained formula has exponent 2. Both use known all-direction radiative transfer for J and known Newtonian enclosed-mass gravity for rotation. The profile is an assumed capture susceptibility, not new mathematical integration or a first-principles interaction. The steeper tail is the proposed physical change; no extra fitting parameter is introduced.

Fit only the 89 training galaxies and freeze:

| Shared parameter | Fitted value |
|---|---:|
| C | 24,191,166.45 Msun/kpc^3 |
| k0 | 0.1092896327 /kpc |
| s=a/R_disk | 5.810229602 |

The fit increases the scale relative to the retained s=4.03298. Therefore the steeper tail does not guarantee a smaller deposit density at every measured radius. The exposure C/k0 stays linked to attenuation as required by the assumed transport energy bookkeeping; no stellar source history or cosmic age was supplied or inferred.

## Predictions on other galaxies

Errors below give each galaxy equal weight. The objective is mean squared logarithmic speed residual, not an observational chi-square. The ordinary matter baseline and published inputs are unchanged.

| Sample | Galaxies | Retained speed RMS (km/s) | New speed RMS (km/s) | Retained log RMS | New log RMS |
|---|---:|---:|---:|---:|---:|
| Training | 89 | 30.70873 | 31.45867 | 0.141004 | 0.142989 |
| Validation | 29 | 33.18419 | 33.86901 | 0.119980 | 0.121644 |
| Test | 31 | 23.99377 | 24.70532 | 0.098545 | 0.101195 |

Only 6/29 validation galaxies and 4/31 test galaxies have lower km/s errors. On logarithmic error the counts are 6/29 and 3/31. This degradation is not just an aggregate change driven by one object. These are previously exposed partitions; no claim of untouched validation is made.

The exponent-3 transparent control also worsens: validation/test RMS are 36.34679/43.79298 km/s, compared with the retained transparent values 35.13687/42.98709. Interception still helps relative to that control, but changing the tail does not improve the retained intercepted model.

Every radial observed/predicted speed is saved in `steep-capture-predictions.json`. Per-galaxy errors and comparisons are in `steep-capture-comparison.json`. These comparisons inherit the existing distance, inclination and ordinary-matter mass-to-light assumptions. They are not uncertainty-marginalized likelihoods or identification of the companion mechanism.

## Verification and consequence

Three optimizer starts converge to the same training loss; all succeed without reaching parameter bounds. The inherited finer integration changes group RMS by less than 0.001 km/s. Five independent upstream integral checks agree with the analytic recurrence within 3.5e-11 relative error. Input hashes and equality of observed radii/speeds with the previous predictions are checked.

A faster-declining capture tail alone is not a supported solution to the current motion/lensing mismatch. Keep the exponent-2 reference and preserve this failed candidate. No lens calculation is presented for this branch, because it already fails to improve its declared rotation comparison. Deposit formation, retention, energy supply and a common physical optics/gravity field remain unresolved.

Reproduce:

    python research_work/results/isotropic-galaxy-transfer/run.py --steep-capture
    python research_work/results/isotropic-galaxy-transfer/steep-capture-compare.py
