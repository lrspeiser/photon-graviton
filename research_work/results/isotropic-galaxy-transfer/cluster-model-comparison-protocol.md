# Cluster comparison scope and equations

Declared before executing the new NFW/MOND fits. Reuse Coma's six vector-reconstructed Kubo Figure 2 shear points, retaining negative values and published plotted errors. Fit bins 1–3; freeze and predict bins 4–6. These outer bins have already been explored, so this is not blind validation. Use radius ratios, positive amplitude and a scale in [0.01,100] times the first radius for each two-parameter shape. Minimize the inner diagonal squared residual sum. The point-mass control has only one amplitude. No expansion-derived physical radius is adopted.

## Formulas and ownership

* Companion shapes reuse the project's proposed deposition profile rho proportional to J(r/a)/(1+(r/a)^2)^2, with optical depth T=0 or 100 fixed in the earlier Coma pilot. Angular transport and projection are known mathematics; the companion interpretation is our hypothesis. Neither T was derived from cluster observations. Both amplitude and a are cluster fits; this is NOT the frozen galaxy parameter test.
* Known NFW: rho=rho_s/[x(1+x)^2]; dimensionless enclosed mass m=ln(1+x)-x/(1+x), and acceleration g=m/x^2. It is fitted here as a total shear shape. An actual baryons-plus-dark-halo test requires measured baryons separately.
* Known simple MOND interpolation: g=g_b/2+sqrt(g_b^2/4+a0*g_b). For this explicitly simplified compact-baryon comparison, g_b=GM_b/r^2, r_M=sqrt(GM_b/a0), x=r/r_M, giving g/a0=1/(2x^2)+sqrt(1/(4x^4)+1/x^2). A point mass is not Coma's measured gas profile. Allowing both normalization and r_M to fit relaxes the relationship imposed by known geometry and a fixed a0, so this is a shape diagnostic, not a test of the a0 fitted to galaxies.
* Known spherical weak-lensing projection, with equal lensing potentials stipulated for MOND: deflection proportional to 2 integral_0^infinity g(r)b/r dz; tangential shear proportional to deflection/b minus its b derivative. Using z=b tan(theta), the dimensionless shape is 2 integral_0^(pi/2) cos(theta)[g(r)-r g'(r)] dtheta. Geometry and overall factors are absorbed into the fitted amplitude. MOND's nonrelativistic acceleration equation alone does not specify light bending. The equal-potential completion is an assumption here, not a new derivation.

Check this projection against the Newtonian point-mass result 4/b^2 and double its quadrature nodes. Use weak shear rather than reduced shear, bin-center sampling rather than source-weighted bin averages, and diagonal errors rather than unavailable covariance. These approximations and the reconstructed input prevent precision likelihood claims.

## Why this cannot test one-third retention yet

Our proposed retention use of the known Hill function is eta=X^(1/3)/(1+X^(1/3)). For a fixed cluster X, eta only multiplies the density normalization. Fitting a free amplitude B absorbs eta exactly: B_effective=B*eta. Every positive eta has the same best shape fit. Thus this experiment cannot identify the exponent, verify the external energy supply, or show that the fitted mass is physically deposited companion energy.

A physical comparison needs the cluster's gas and stellar distributions (which set the ordinary gravity), a declared cluster replacement for the disk luminosity/scale proxy X, source geometry (which converts mass to shear), and the supply/capture/retention normalization frozen without fitting the cluster lensing amplitude. Without these, a good curve can conceal the wrong amount or location of gravity. The same restrictions must apply to the MOND and halo alternatives. Use that model to fit a designated cluster and predict another with no target-specific rescue adjustments.

The current repository has Coma shear, plus Virgo/M87 supply pilots without matched lensing actuals. Do not count the latter as a second successful cluster test, or convert inferred dark-matter masses into independent observations. No cosmic age, size or formation time is assumed here. All six goals remain open.

Sources: [Coma measurements, Kubo et al.](https://arxiv.org/abs/0709.0506); [NFW profile](https://arxiv.org/abs/astro-ph/9611107); [MOND review and relativistic qualifications](https://arxiv.org/abs/1112.3960). The local input provenance is in ../cluster-observation-readiness/kubo-figure-data.json.
