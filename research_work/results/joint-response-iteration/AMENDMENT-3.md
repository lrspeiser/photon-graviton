# JR-1 amendment 3: post-transfer population follow-up

21 September 2026. R9 has now been frozen and its two out-of-fit lens predictions exposed: stellar fractional RMS averages 3.208 percent, Einstein-angle RMS is 13.570 percent, with both angles too small. Those results are retained unchanged. The 31 SPARC comparison galaxies have also been exposed. No subsequent formulation can claim those same results are new independent confirmation.

R9's four fitted lenses all use positive individual stellar-population mass offsets (0.0757 to 0.1624 dex), while their shared stellar scale is 0.8956. This suggests that imposing one stellar-population normalization on both disks and spheroids can push the per-lens mass nuisance parameters in the same direction. It is a modeling hypothesis, not a discovery of a different physical IMF or proof of the companion mechanism.

R10 permits ONE extra shared spheroid stellar-population normalization u_sph:

    Mstar_true = u_star * (Mdisk + u_sph*Msph)
    S = u_star * (Mdisk + w_sph*u_sph*Msph) + Mgas.

The spheroid fraction used for the geometry-dependent core/shape becomes u_sph*Msph/(Mdisk+u_sph*Msph). Thus changing the population normalization changes source gravity, production proxy and shape consistently. For spherical SLACS models all stellar mass is spheroidal. The log10(u_sph) bound is [-0.3,0.4], with an explicit zero-centered width-0.15 dex working population prior. This prior represents a declared nuisance range, not a fresh stellar-population measurement. The original four individual log-mass priors and orbital priors remain unchanged.

Fit only the original 89 SPARC training galaxies and four fit lenses, retaining all observational velocities, covariances, photometric shapes and lens angles. Score all samples afterward. The follow-up is residual-driven exploratory development, not a replacement of R9's independent-of-this-fit transfer record. Do not reinterpret any reduced error as blind validation or select a new per-object companion normalization. Keep the same equal-potential matter/light coupling and perform frozen-parameter resolution checks on the follow-up.
