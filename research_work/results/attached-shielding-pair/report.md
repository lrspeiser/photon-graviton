# Can ordinary-matter shielding reverse the loading order?

**Yes, for the two critical locations in the fixed model, with a fitted absorption coefficient.** This is pairwise feasibility, not an independent prediction, a whole-galaxy fit, or an observed confirmation.

The existing target requires loading 1.19439 per unit ordinary mass at R=8 kpc, z=1 kpc, phi=0, versus 0.198146 at R=0.5 kpc, z=0.1 kpc, phi=pi/4. Their ratio is 6.02784. For a nondecreasing local capture efficiency with depth and common duration, the shallower location would need at least this exposure contrast. That is a conditional necessary requirement; changing histories or redistribution changes the argument.

## Postulates and known transfer mathematics

Optional project postulates, not claimed novel: companions enter isotropically at a spherical boundary, travel on straight rays, and are absorbed with opacity per length beta=kappa*rho_b. The coefficient kappa is common and independent of depth in this diagnostic. Captured energy stays attached locally. We omit self-gravity feedback, focusing, scattering, receiver saturation and moving matter.

Known absorption law: I(x,n)/I_boundary=exp[-kappa*Sigma(x,n)], where Sigma is the integral of ordinary mass density along the incoming ray. The angle-averaged exposure fraction is T(x)=(1/4pi) integral exp[-kappa*Sigma] dOmega. With common exposure time and constant kappa, loading per ordinary mass is proportional to T. We solve T_shallow/T_deep=6.02784 for kappa. Thus matching this ratio is calibration, not validation.

The density includes the same stellar bar, nuclear components, stellar and gas disks, and softened center used in the attached-loading audit. Densities at both endpoints reproduce their earlier values to relative 1e-12. No dark-halo density is introduced.

## Results and resolution

| Boundary | Angular nodes (mu x phi) | Fitted kappa (kpc^2/Msun) | Shallow exposure fraction | Deep exposure fraction |
|---|---|---:|---:|---:|
| 30 kpc | 32 x 64 | 9.68440e-10 | 0.913243 | 0.151504 |
| 30 kpc | 64 x 128 | 9.67569e-10 | 0.913691 | 0.151579 |
| 60 kpc | 32 x 64 | 9.68468e-10 | 0.913201 | 0.151497 |
| 60 kpc | 64 x 128 | 9.67656e-10 | 0.913560 | 0.151557 |

The coarsest 16 x 32 calculations are also retained. Doubling the last angular resolution and segment quadrature changes fitted opacity by about 0.09%; changing the outer boundary from 30 to 60 kpc changes the finest estimate by about 0.009%. These are descriptive convergence comparisons, not predeclared validation gates. Angular mean columns converge less tightly because dense, strongly extinguished sightlines contribute little to transmitted intensity. Neither boundary is a universe-size assumption.

At the fitted coefficient the outer location loses about 9% of the incident companion intensity and the inner one loses about 85%. Absorption removes companion energy from rays into the postulated receiver; it does not create energy. No photon supply, total stored gravitational mass, capture momentum balance or stable storage configuration is calculated here.

## What follows

The earlier monotone-loading failure depended on equal exposure. This calculation demonstrates that the declared ordinary-matter geometry can violate that assumption in the needed direction for one critical pair. It does not establish that the same coefficient reproduces all heights, bar angles, galaxies or lensing. A capture coefficient that itself grows with depth would need a new self-consistent opacity calculation; the present constant coefficient is not that model.

Next, hold the pair-calibrated coefficient fixed and calculate exposure throughout the existing 240-point development grid. Compare the resulting relative loading with the target using one common normalization, retaining all failures. This will test whether pairwise shielding survives the full geometry rather than fitting every point separately. Independent observations and physical supply remain subsequent requirements. All six research objectives remain open.
