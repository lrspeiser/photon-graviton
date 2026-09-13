# Frozen outer-gravity and lensing predictions from the fitted companion reservoirs

13 September 2026. No new fit and no new outer observational data. These are conditional extrapolations to test next.

## Outcome

The large J1621+3931 reservoir produces a concrete, unusual outer signature. The fitted model has essentially no additional enclosed mass in the inner few hundred kiloparsecs, while its projected companion convergence there is approximately 0.0847. At large radii its implied circular speed rises to about 4246 km/s near 9.17 Mpc. It predicts negative tangential shear at the sampled radii 1 and 3 Mpc and positive shear at 10 Mpc.

This is not an observed speed or shear signal. It is the consequence of preserving the freely fitted density, scale and geometry rather than retuning them outside the inner data. It supplies specific independent tests of the extended-reservoir branch. A good inner fit does not establish that such a reservoir exists or can be supplied with energy.

The other five profiles have much smaller implied circular-speed peaks, approximately 430-636 km/s over the scanned 0.01-10000 kpc range. None of the six profiles shows additional tangential or radial critical curves beyond the ordinary inner pair in that finite scan. Absence of another critical ring therefore does not by itself reject the diffuse large reservoir.

![Frozen outer predictions; markers are calculated radii and connecting lines guide the eye. These are not measured outer data.](outer-companion-predictions.png)

## Equations and provenance

For the frozen spherical effective mass profile,

    v_c^2(r)=G[M_stars(<r)+M_companion(<r)]/r.

This is standard circular-orbit dynamics, not a prediction for the speed of each pressure-supported star or satellite. Actual velocity distributions also need an orbit/population model and environmental treatment.

Project the same stellar and companion densities to surface density Sigma and cylindrical mass M_2d. Under the same weak-field thin-lens response and source geometry as the fit,

    Sigma_crit=c^2/(4*pi*G*Dl*(Dls/Ds)),
    kappa=Sigma/Sigma_crit,
    mean_kappa=M_2d/(pi*R^2*Sigma_crit),
    gamma_t=mean_kappa-kappa,
    g_t=gamma_t/(1-kappa),
    lambda_t=1-mean_kappa,
    lambda_r=1-2*kappa+mean_kappa.

These are established circular-lens relations; see [Bartelmann and Schneider, Weak Gravitational Lensing](https://arxiv.org/abs/astro-ph/9912508). Their application to the effective companion density is conditional on the model's gravitational response. The expanding-universe calibration in other applications of these equations is not adopted as a premise here; the recorded conditional distances and source-distance ratio are retained.

Zeros of lambda_t and lambda_r locate tangential and radial critical curves in this spherical lens model. They are not a reconstruction of real image shapes. Shear/convergence depend on the source distance; the tabulated values apply to the source geometry used by the original fit. Observations of other background populations require their appropriate geometry.

The mass projection uses M_2d(R)=integral_0^(pi/2) M_3d(R/cos(theta))*cos(theta) dtheta. Surface density is independently integrated along the line of sight. The identity dM_2d/dR=2*pi*R*Sigma supplies a numerical consistency check. None of these equations is a newly derived photon-conversion law.

## Outer signatures for the six profiles

Peak speeds are maxima on the stated radial scan, not observed maxima or uncertainty estimates. Angular scales and full projections at every sampled radius are in the result JSON.

| Galaxy | Peak scanned circular speed (km/s) | Peak radius (kpc) | v_c at 100 kpc | v_c at 1000 kpc |
|---|---:|---:|---:|---:|
| J0037-0942 | 429.6 | 51.6 | 389.8 | 147.7 |
| J1112+0826 | 569.0 | 44.2 | 491.0 | 181.1 |
| J1204+0358 | 506.9 | 0.776 | 116.1 | 36.8 |
| J1402+6321 | 635.8 | 51.6 | 572.9 | 217.0 |
| J1621+3931 | 4245.7 | 9.17e+03 | 200.7 | 63.8 |
| J1630+4520 | 509.0 | 1.58 | 196.8 | 63.6 |

## Why the very large reservoir can hide from inner motions

In the fitted J1621+3931 profile, strong attenuation suppresses deposited density near the center; most stored mass lies in a very extended envelope. Spherical exterior shells exert no net Newtonian acceleration on interior stars, but light from a background source still traverses their projected gravitational field. The calculation therefore permits almost stellar-only inner accelerations together with a nearly constant added surface-density contribution to lensing.

That is a geometrical degeneracy of this fitted mass distribution, not proof that companions form a real shell. The total mass is poorly determined by the inner fit because changes far outside it can preserve a similar central projected contribution. This explains the earlier enormous mass inventory without assuming that stars inside it must already orbit at thousands of km/s.

| Radius (Mpc) | Angular radius (arcmin) | Implied v_c (km/s) | Companion convergence | Total reduced tangential shear |
|---|---:|---:|---:|---:|
| 0.01 | 0.046 | 352.75 | 0.084744 | 0.297231 |
| 0.03 | 0.137 | 316.27 | 0.084746 | 0.074613 |
| 0.1 | 0.456 | 200.72 | 0.084769 | 0.009849 |
| 0.3 | 1.369 | 116.39 | 0.084974 | 0.000984 |
| 1 | 4.562 | 63.76 | 0.087505 | -0.001454 |
| 3 | 13.686 | 1520.51 | 0.116744 | -0.015756 |
| 10 | 45.620 | 4236.46 | 0.015023 | 0.034307 |

The negative shear at 1 and 3 Mpc denotes radial image stretching in this sign convention, not negative mass or repulsive gravity. It occurs where local projected density exceeds its interior average. By 10 Mpc the shear is tangential again. The discrete values do not locate the exact zero crossings.

Near 3 Mpc, corresponding to approximately 13.7 arcmin for the adopted lens distance, the predicted reduced shear is about -0.0158. At 10 Mpc, roughly 45.6 arcmin, it is about +0.0343. This broad signed pattern is a prospective test with appropriately selected background sources. Foreground structure, source-redshift distributions, nonspherical geometry and environmental masses must be modeled before comparing to actual data. The negative sign is not unique evidence for companion particles; other projected mass arrangements can produce radial shear.

A nearly uniform projected component can have little shear even when convergence is nonzero, since gamma_t is the difference between local and mean convergence. Shape measurements alone do not directly measure an arbitrary constant convergence. This is why both the radial extent and departures from uniformity matter for an independent test.

## Critical curves and weak-field diagnostics

Every model has one tangential and one radial critical curve found in the finite scan. The tangential radius reproduces the consumed catalogue lens calibration; it is not a fresh prediction. No additional outer critical curve is found between 0.01 and 10000 kpc. A finite grid and this spherical approximation cannot prove that no other image structures exist under different source conditions.

For J1621+3931 the maximum scanned compactness 2GM/(r*c^2) is about 4.01e-4. The computed near-central potential magnitude is about 4.79e-4 c^2, with an omitted exterior-companion potential bound below 6.07e-8 c^2. Thus the very large mass is sufficiently diffuse that these diagnostics remain small relative to unity. They do not prove a relativistic completion, stationary support or dynamical stability. The adopted isolated-system and thin-lens approximations also become important when extrapolating over megaparsecs; this report does not establish their accuracy for the actual environment.

## Verification and numerical correction

The companion mass integral uses 8001 radial points and 192 incoming-direction angles, with the numerical outer radius inherited from the fit. Integrated masses change by less than 4.32e-6 fractionally from the fitted-grid inventories. The earlier analytic exterior-mass bound still describes the unintegrated tail; the numerical outer radius is not a measured halo boundary.

A first unsplit angular mass projection failed the refined central critical-curve check for the most extended profile. The final calculation splits the integration at capture-scale crossings, resolving the distant-envelope contribution. The final recorded outputs supersede that unsuccessful preliminary calculation.

At reported radii, increasing projection quadrature from 192 to 384 changes kappa by at most 2.26e-11 and mean kappa by at most 1.88e-10. Refining critical-curve evaluation to order 768 gives eigenvalue residuals below 1.22e-8. The independent projected-mass derivative check agrees in kappa to within 3.07e-7, and the catalogue mean-kappa normalization is reproduced within 1.55e-6. These are numerical checks of a frozen model, not observational error bars or evidence that the model is correct.

The plotted markers are the seven reported radii; lines merely connect them. The peak-speed estimates use the separate 801-point scan. Surface-density and mass contributions are recorded separately for stars and companions where relevant, allowing the extended component to be distinguished from the stellar signal.

## Consequence and next observational comparison

The extreme free fit now makes concrete outer predictions instead of merely demanding a large unspecified energy budget. A comparison with background-source shear and with independently modeled outer tracers can test whether this reservoir is allowed. No new outer observations have been obtained or tested in this calculation; neither rejection nor confirmation is claimed.

The reference exact-third law remains unchanged. These predictions concern the freely fitted branch and cannot be attributed to a derived universal capture mechanism. Supply, redistribution, stopping and stability must still account for the density profile if it survives independent gravitational constraints.

## Reproduction

Run outer-companion-predictions.py and plot-outer-companions.py in research_work/results/companion-extensions/. The result JSON records all seven radii per system, source geometry, convergence/shear, enclosed mass, circular speeds, critical radii, potential diagnostics, numerical checks and input hashes. Choices are in outer-companion-predictions-protocol.md.
