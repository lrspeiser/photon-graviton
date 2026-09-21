# JR-1 amendment 1: source geometry, not a separate lensing gain

21 September 2026. Exploratory, declared after R0-R6 training/validation development and before the added fits. No SPARC comparison-role or two-lens out-of-fit results have been scored in this run.

R0's fixed smooth envelope gives SPARC train/validation mean galaxy RMSE 16.95/19.99 km/s but underpredicts the fitted lens stellar motions (mean RMSE 107.26 km/s). Jointly fitting the initial global mass/size law improves the fitted lens block while worsening SPARC. R6 reaches 29.04/38.51 km/s on SPARC and 21.50 km/s on fitted lens motions, with 20.38 percent lens-angle RMS. These are model-development residuals, not an observational acceptance.

This motivates a geometry-dependent source proxy. The revised law allows different efficiency for spheroidal and disk source components, applying ONE coefficient to all spheroids (including SPARC bulges), never fitting an independent amplitude for an individual lens or a photon-specific multiplier:

    S = u_star * [Mdisk + w_sph * Msph] + Mgas
    f_sph = Msph/(Mdisk+Msph)
    rc = c * Re * C^dc * 10^(-core_sph*f_sph).

C remains based on actual input baryonic mass, NOT efficiency-weighted S. SPARC Msph is the archived spherical bulge mass inferred from its bulge contribution; disk mass is integrated photometric surface density. For spherical SLACS deprojections the stellar component is treated as spheroidal. This source decomposition is a modeling assumption, not a measured microscopic production rate or proof that geometry alone drives the observed discrepancy. Source-population systematics can masquerade as efficiency changes.

Permit universal log10(w_sph) in [-1,2] and core_sph in [-2,2]. Preserve all earlier results. Fit both a power-law occupation version and the steady production/release version with radial orbital anisotropy. Constants and remaining bounds are inherited. No mass/light datum or uncertainty is changed. No transfer score is used to choose this revision.

The resulting sphere/disc distinction must eventually be derived from spatial interactions and phase-space source structure, not elevated from a successful regression to fundamental physics. A field predicted around an entire galaxy is not yet an algorithm for assigning independent clouds to each star; do not claim Solar-system consistency from this stage.
