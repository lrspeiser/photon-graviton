# Matched rotation comparison: companion, MOND, NFW

Declared before fitting. Use the same 149 SPARC galaxies, 3150 accepted radial points, published distances/inclinations and stellar mass-to-light choices 0.5 disk/0.7 bulge. Reuse the 89/29/31 train/validation/test split. The splits are exposed; no blind evidence claim. Fit with the inherited equal-galaxy mean squared log10 speed residual. Report km/s RMS, log RMS, and a separate mean per-galaxy diagonal measurement-error score using published velocity errors. The latter omits shared distance/inclination/mass-model covariance and intrinsic scatter, so it is not a calibrated goodness-of-fit probability.

Predictions to compare:

1. Ordinary matter only, no fitted parameter.
2. The existing exact-one-third companion candidate, three shared fitted parameters, frozen from its training-only run. No rerun/refit to comparison targets.
3. MOND simple algebraic interpolation: g=0.5*g_b+sqrt(0.25*g_b^2+a0*g_b), v=sqrt(r*g). Fit one shared a0 on training, in [1e-12,1e-8] m/s^2. Also show fixed a0=1.2e-10 as a conventional benchmark. This is an established phenomenological approximation, not a full nonspherical AQUAL/QUMOND field solve or an external-field model. Reference: https://academic.oup.com/mnras/article/455/1/449/985618 for interpolation context.
4. NFW halo: rho=rho_s/[x*(1+x)^2], M=4*pi*rho_s*rs^3*[ln(1+x)-x/(1+x)], x=r/rs. This is a known halo profile (https://arxiv.org/abs/astro-ph/9611107). To make galaxy-to-galaxy predictions without fitting test galaxies, explicitly POSTULATE rs=s*R_disk and rho_s=rho0*X^p, X=(L3.6/10^9 L_sun)/(R_disk/kpc)^2. Fit three shared parameters on training: log10 rho0 in [-2,12] Msun/kpc^3, log10 s in [-1,2], p in [-2,2]. This is our restricted NFW scaling benchmark, not a comprehensive LCDM prediction or abundance-matching model. No expansion parameters or cosmological density normalization are needed.

Additionally show an ordinary two-parameter NFW flexibility diagnostic: fit rho_s,rs to inner 60% of radial points in each validation/test galaxy, then predict its outer 40%. Compare those outer points with frozen MOND/companion/global NFW. This uses target-specific inner velocities, unlike the frozen models, and must be labeled separately rather than counted as equal-information superiority. Bounds log10 rho_s [-2,12], log10 rs/kpc [-2,3], three starts. No outer target is used in a per-galaxy fit.

These benchmarks compare selected rotation prescriptions only. They do not rank the full theories on cosmology, lensing, Solar System tests or graviton microphysics. Parameter counts, boundary fits and limitations must accompany results. Do not identify the previously fitted empirical force baseline as MOND. Preserve all outputs and original candidate constants.
