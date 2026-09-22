# JR-8: resolved-map source response and exact spatial susceptibility

Declared locally before the calculations below; published to GitHub with the completed results. Baseline main 20255778863150392baef9702b6ebb5df35d1342. Inputs are the exact three FITS products and 12-cluster extract in the attached JR-7 archive. This is exploratory physical construction with observed map shapes, not a calibrated gravity fit or a fresh blind sample. No archival inputs or previous results are overwritten.

## Deliverables

1. Infer smooth positive continuum/Halpha source-shape maps from their measured spaxels, respecting masks, photometric ellipse and resolution; obtain stellar/gas velocity harmonics on actual stellar spectral bins, not repeated pixels. Report radial variation and the distinction between projected gradients and physical pattern speed.
2. Feed the observed spatial shapes into the existing positive, energy-accounted 3D transport operator. Apply the SAME four JR-7 rate laws. Compute density, common force/lensing readouts, and paired phase controls. No per-object gravity normalization is fitted.
3. Derive an exact adjoint sensitivity of the stationary companion content to local gas changes. Verify it by independent finite differences. Locate positive/negative regions and whether a scalar gas-content statistic loses sign information.
4. Use the actual cluster electron-density profiles to calculate where the local rate balance changes sign as a UNIVERSAL reference density is swept. Do not infer cluster gravity or fitted density constants from this exercise.

## Observations and reconstruction

Select continuum SPX_MFLUX >0, SPX_SNR>=5 and Halpha Gaussian flux positive with SNR>=5 and zero flux mask. Read named channels from FITS headers. Fixed named products must reproduce JR-7 archived arrays exactly. FITS input WCS ellipse is an assumed thin-disk geometry, not measured physical depth. Coordinates use R/Re and intrinsic-ellipse azimuth. Aperture is min(2Re,90th percentile of valid radius), not chosen from gravitational residuals. Positive log-flux reconstruction uses six radial Gaussian basis functions and angular harmonics m<=2; fixed regularization and a common valid footprint. Full original maps, rejected pixels, mask coverage and reconstructions are retained. Gaussian vertical lifts conserve integrated source inventory: source scale height .20Re; gas scale height .20Re primary and .45Re sensitivity. Halpha is emissivity, NOT total gas mass. Gas proxy proportional to emissivity is primary, its square root is a separate sensitivity; neither is relabeled a measurement of gas density. Each source is normalized to unit source power and unit proxy gas inventory in model units. Absolute photon supply, gas-density calibration and gravitational normalization are not supplied by these shape-only runs.

## Spatial solves

Rates are exactly k+=.5 exp(a h), k-=.5 exp(b h^2), h=rho_proxy/(rho_proxy+.15), with (a,b)=(0,0),(2,0),(0,4),(2,4). Radiation-like transport D=.5, loss=.5; companion D=.05, loss=.05. Reflecting radial/vertical boundaries and periodic angle. Currents are prescribed 0,+.6,-.6 in model units; they are NOT measured stellar pattern speeds. The source and gas pattern are frozen snapshots; their long-term persistence is not derived. Compare gas phase shifts 0,90,180,270 degrees at fixed radial inventory; these are counterfactual controls, not observed changes. Primary mesh 10x12x24; selected mixed phase pairs 0/90 and each proxy/height use finer 14x16x48 and, where needed, a further refinement. The whole source and gas map rotate together for covariance checks. Forces and deflections use the same softened 3D source, G=c=1, with softening .15Re. Changes in physical units are not inferred.

## Adjoint

For A y = j, E_C=w^T C, solve A^T psi=(0,w). At fixed source and transport,
 dE_C/d ln rho_j = (psi_C-psi_P)_j [a k+ P - 2b h k- C]_j h_j(1-h_j).
This expression is exact for this finite-volume model. Verify cell perturbations and fixed-total-inventory paired perturbations at both signs. Integrated finite perturbations must use refitted field solutions, not merely a linearized prediction.

## Numerical interpretation

Require positivity to roundoff, relative linear residual<1e-8 and channel energy imbalance<1e-8. Adjoint finite-difference relative discrepancy below 1e-5 for resolved nonzero derivatives. Preserve every refinement difference; call a paired sign resolved only if its magnitude exceeds 3 times its finest paired-grid change. No acceptance probability is assigned to map rotations or per-pixel regression. Reconstructed light/gas covariance is separated from their shared radial profile. No new lens, Solar-system or cluster gravitational success is claimed without corresponding calibrated observables.

## Executed scope clarification

The final finer-grid campaign covers the primary mixed family, primary gas proxy and height, and the 0/90 phase pair on all three objects at 14x16x48 and 18x24x64. Other proxy/height sensitivities remain on the primary mesh and are not described as equally refined. Publication does not retroactively enlarge the executed scope.
