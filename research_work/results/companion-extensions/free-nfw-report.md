# NFW scale and strength released in the lens/motion comparison

13 September 2026. Bounded six-parameter-per-galaxy descriptive fit. This materially qualifies the earlier restricted NFW comparison.

## Outcome

Allowing each NFW halo's scale and strength to vary reduces the total all-motion chi-squared to 49.71, compared with 111.27-113.56 for the fixed companion profile with four stellar nuisance parameters. The restricted NFW shape comparisons had totals 162.62-194.57. Thus their worse scores did not establish a general advantage of the companion radial shape over fitted NFW halos.

The largest change is J1402+6321: its motion chi-squared becomes 4.49, versus 63.24/64.97 for the companion cases. This improvement requires a substantially different halo contribution and stellar configuration. The halo supplies about 76.8% of the required bending at the catalogue radius; the stellar gradient reaches h=9 and the orbit transition reaches its lower bound of 0.1 Re. These are fitted, boundary-dependent values, not independently measured properties.

The result does not establish that this is the correct dark halo or a stable physical solution. It adds two local fitted parameters per galaxy, giving 36 local parameters for forty motion measurements, with lens angles separately consumed to calibrate stellar mass. All stellar outer-anisotropy endpoints reach bounds; several halo scales do too. Raw chi-squared alone cannot rank explanatory power or establish independent predictive success.

## Equations and provenance

Retain the established NFW form and analytic enclosed mass:

    rho(r) = rho_s/[x(1+x)^2], x=r/rs,
    M(<r) = B [ln(1+x)-x/(1+x)], B=4*pi*rho_s*rs^3.

This is known halo phenomenology, not new photon-companion physics; see [Navarro, Frenk and White](https://arxiv.org/abs/astro-ph/9611107). Its cosmological formation scenario is not imposed as a fact of the hypothetical universe. The ideal profile is untruncated here: it has a logarithmically divergent total mass, so the fit is not a finite energy inventory or a formation calculation.

Let alpha_required be the physical deflection required by the catalogue angle and adopted distance ratio. Let d_star(h) be the stellar deflection per unit total stellar mass for a gradient h, and d_NFW(rs) the deflection per unit B. Fit a halo bending fraction f and scale rs, and set

    B = f*alpha_required/d_NFW(rs),
    M_star = (1-f)*alpha_required/d_star(h).

These are algebraic lens-calibration relations within the existing spherical weak-field calculation. They ensure positive stellar mass for f<1 and use the same stellar and halo sources for the Jeans acceleration and light bending. The chosen range is 0<=f<=0.95; it is a search bound, not a measured baryonic fraction. F is a fraction of bending, not a total halo mass fraction or capture efficiency.

The four stellar variables remain h, beta0, beta_infinity and log(ra/Re). Add log(rs/Re) in the range 0.1<=rs/Re<=100 and f in [0,0.95]. All six variables fit the full motion covariance. The companion normalization, capture scale and luminosity proxy do not constrain the new halo; earlier fits supply starting configurations only. The stellar mass-to-light gradient still transitions at Re, and all previous stellar bounds and necessary orbital conditions remain enforced.

## Data and parameter accounting

Use the same six exposed galaxies, forty motion bins, PSF, light profiles and conditional nonexpanding geometry. Both historical population-proxy cases have identical geometry and light inputs, verified before fitting. With freely determined halo strength there is one physical NFW fitting problem per galaxy, not two independent problems. Its six results are compared against both historical companion proxy assumptions.

The companion comparison has four local stellar nuisance parameters per galaxy plus inherited shared companion parameters calibrated in prior galaxy work. This NFW comparison has six fitted local parameters per galaxy. Catalogue lens angles calibrate M_star in both. No cosmic concentration-mass relation, independent halo-mass prior or independent stellar-population gradient constraint is supplied. No p-value or information-criterion winner is inferred from the nominal parameter counts, boundaries and raw residual totals.

## Scores and fitted values

| Galaxy | Companion Chabrier chi-squared | Companion Salpeter chi-squared | Free NFW chi-squared | NFW conditional outer residual | Halo bending fraction |
|---|---:|---:|---:|---:|---:|
| J0037-0942 | 10.444 | 9.739 | 9.265 | 0.319 | 0.4901 |
| J1112+0826 | 18.668 | 18.700 | 17.641 | 3.070 | 0.6431 |
| J1204+0358 | 14.333 | 14.424 | 13.945 | 2.906 | 0.0000 |
| J1402+6321 | 64.974 | 63.243 | 4.494 | 1.251 | 0.7683 |
| J1621+3931 | 3.409 | 3.426 | 2.820 | -0.221 | 0.3694 |
| J1630+4520 | 1.737 | 1.742 | 1.548 | 0.333 | 0.7728 |

The new inner chi-squared contribution is 30.02 and the conditional outer residual-square sum is 19.69, totaling 49.71. The respective companion contributions are 61.70/60.44 and 51.87/50.84. The outer bins are fitted in every comparison here; none is a withheld prediction.

| Galaxy | h | beta0 | beta_infinity | ra/Re | rs/Re |
|---|---:|---:|---:|---:|---:|
| J0037-0942 | 9.0000 | 0.3750 | 0.9500 | 0.2395 | 8.0785 |
| J1112+0826 | -0.8000 | 0.3750 | -2.0000 | 0.4114 | 0.1000 |
| J1204+0358 | 5.1264 | -2.0000 | 0.9500 | 1.4526 | Unidentified: halo absent |
| J1402+6321 | 9.0000 | -0.1242 | 0.9500 | 0.1000 | 8.8354 |
| J1621+3931 | 7.6071 | 0.3750 | 0.9500 | 1.0020 | 100.0000 |
| J1630+4520 | -0.8000 | 0.2969 | -2.0000 | 0.5088 | 0.1000 |

J1204+0358 chooses f=0, recovering the stellar-only control. Its halo radius is unidentifiable in that limit; the optimizer's stored numerical value is not a physical inference. J1112+0826 and J1630+4520 reach rs/Re=0.1; J1621+3931 reaches rs/Re=100, about 877 kpc in the adopted geometry. These boundary values expose unresolved scale dependence, not measured halo radii. No bending-fraction upper bound is reached.

All six outer-anisotropy endpoints are at either -2 or 0.95. Four stellar gradients reach a bound, and J1402+6321 reaches ra/Re=0.1. Even though total fit quality improves substantially, these stellar choices still need independent evidence and physical orbital completion. The remaining conditional outer residuals for J1112+0826 and J1204+0358 are about 3.07 and 2.91, respectively; these are model/covariance-dependent diagnostics, not universal significance measures.

## Numerical method and verification

The unit NFW deflection is tabulated at 257 logarithmic rs/Re values and interpolated in log deflection for the initial search. Enclosed mass and force use the analytic expression directly, with a small-x series for numerical stability. The three best successful candidates for each galaxy are polished using direct deflection quadrature, and all final reported scores use direct quadrature.

Use 25 starting configurations per galaxy: six previous NFW shape/population solutions, the stellar-only fit, and eighteen fixed grid starts. Twenty-four or twenty-five converge in each galaxy, and all three direct polishing runs succeed per galaxy. Multiple local minima and all successful objectives are recorded. This is a bounded search and does not prove a global optimum.

At the old NFW parameters, direct evaluation reproduces stored objectives to within 2.92e-13; the zero-halo control is reproduced exactly at recorded precision. Every final objective is no worse than its best tested feasible starting point. At final parameters the initial interpolation differs from direct deflection by at most 4.25e-8 fractionally; final results do not rely on that approximation. Lens closure is within 2.23e-16 fractionally, and total/inner/conditional-outer covariance decomposition agrees within 1.78e-15.

Masses and predicted second moments remain positive. All final profiles satisfy the imposed central and refined 8193-point density-slope/anisotropy condition. That is still only necessary for the specified separable augmented-density class: no nonnegative distribution function or stable equilibrium has been constructed. The model's extrapolated central light profiles and spherical approximation remain material assumptions.

## What this changes

The earlier NFW comparison was restricted in a consequential way. Releasing the inherited halo normalization and prescribed scale reverses the raw residual comparison. The records of that restricted test remain useful as shape diagnostics, but must not be cited as evidence that companions outperform dark matter generally.

The main unresolved companion comparison is now more specific: a fixed galaxy-transferred companion profile plus four stellar freedoms is being compared with a halo whose amplitude and size are fitted separately for every object. A next controlled test can allow the same amplitude and radial-scale freedoms to the companion profile, clearly labeling it as an inverse diagnostic rather than a successful energy-supply prediction. Whether the required changes can be generated by transport, capture and retention with shared constants is a separate physical requirement.

These results do not validate an expanding universe or prove a dark-matter particle exists. They show that a more freely adjusted additional mass distribution can fit these motion and lens-calibration constraints substantially better, at the cost of flexibility and boundary dependence. Photon conversion, absolute energy supply, event timing and stable retention remain unresolved by either this stellar fitting exercise or the restricted companion comparisons.

## Reproduction

Run `python research_work/results/companion-extensions/free-nfw.py`. The independent diagnostic script writes free-nfw-results.json and preserves all earlier fits. It records six unique systems, every motion bin, halo and stellar parameters, optimizer results, constraints, numerical checks and input hashes. Pre-execution choices are in free-nfw-protocol.md.
