# Hollow reservoirs require more than isotropic collisionless particle support

13 September 2026. Necessary-condition calculation for an optional particle interpretation. No refitting, new observational comparison or change to the exact-third reference.

## Outcome

Four of six freely fitted reservoirs, five of six transferred Chabrier reference reservoirs, and all three extended J1621 fixed-scale branches fail a necessary condition for spherical, stationary, bound collisionless particles with isotropic velocities. Their deposited density rises outward somewhere. The compact free J1204 and J1630 profiles, transferred J1204 profile and compact J1621 scale-0.3 profile have no sampled violation; that is not proof of a valid distribution function.

This is a constraint on the deposited particles themselves. Earlier stellar orbital slope tests concerned a different component, so passing those tests does not establish companion support. Traveling companions are also not the bound particle population assessed here.

## Why a hollow profile fails this completion

Assume a Newtonian spherical attractive potential, nonnegative total gravitating mass, relative potential Psi=Phi(infinity)-Phi(r), and bound isotropic deposited-particle distribution f(epsilon)>=0 depending only on relative energy epsilon=Psi-v^2/2. The standard velocity integral gives

    rho_d(Psi)=4*pi*sqrt(2)*integral_0^Psi f(epsilon)*sqrt(Psi-epsilon) d epsilon,
    d rho_d/d Psi=2*pi*sqrt(2)*integral_0^Psi f(epsilon)/sqrt(Psi-epsilon) d epsilon >= 0.

Since dPsi/dr=-G*M_total(<r)/r^2<0 where enclosed mass is positive, d rho_d/dr<=0. Positive outward density slope contradicts a nonnegative f(epsilon). This is a necessary condition, not a full Eddington inversion. The relation is established spherical distribution-function mathematics, not newly invented companion physics; see [spherical DF formalism](https://galaxiesbook.org/chapters/I-04.-Equilibria-of-Collisionless-Stellar-Systems_6-Spherical-distribution-functions.html) and [official galpy documentation](https://docs.galpy.org/en/v1.12.0/tutorials/distribution_functions/spherical_dfs.html).

The result assumes a bound collisionless population with isotropic velocities and no additional non-gravitational forces. It does not reject anisotropic angular-momentum-dependent distributions, coherent wave stresses, collisional matter, ongoing flux or modified force laws. Those alternatives require explicit equations; their availability is not evidence that they work.

## When the prescribed capture profile becomes hollow

Let x=r/ac and T=k0*ac. For the existing isotropic-illumination optical-depth expression, expansion about x=0 gives

    tau(x,mu)=T[pi/4 + mu*x - (3*pi/8)*(1-mu^2)*x^2 + O(x^3)],
    J(x)=exp(-pi*T/4)[1+(pi*T/4+T^2/6)*x^2+O(x^4)],
    rho_d(x)/rho_d(0)=1+B*x^2+O(x^4),
    B=pi*T/4+T^2/6-2.

A hollow center occurs for B>0, equivalently

    T > [-3*pi/2 + sqrt((3*pi/2)^2+48)]/2.

For the frozen k0, this is T>1.833275, or ac>8.990920 kpc. The series uses angular averages of mu and mu^2 and the known Taylor expansion of attenuation. This is a new explicit derivation within the project's prescribed transport profile, not a first-principles derivation of photon conversion or a claim that the mathematics is novel. The factor D and one-third retention amplitude cancel, so adjusting normalization cannot remove the hollow shape. Lowering ac or k0 changes both attenuation and gravity and requires renewed comparison with observations.

For every examined nonzero profile, the central-series classification agrees with the numerical outward-slope result. The threshold is sufficient to identify a central failure; being below it alone is not a universal guarantee of a positive f or monotonicity at every radius.

## Executed results

| Recorded profile | Isotropic necessary-condition test | Density-peak radius (kpc) |
|---|---|---:|
| J0037-0942 free | Fails | 15.99 |
| J1112+0826 free | Fails | 10.76 |
| J1204+0358 free | No violation found | 4.137e-05 |
| J1402+6321 free | Fails | 13.13 |
| J1621+3931 free | Fails | 3526 |
| J1630+4520 free | No violation found | 7.416e-05 |
| J1621 scale 0.3 | No violation found | 0.0002632 |
| J1621 scale 10 | Fails | 149.8 |
| J1621 scale 30 | Fails | 688.8 |
| J1621 scale 100 | Fails | 3526 |
| J0037-0942 transferred | Fails | 11.68 |
| J1112+0826 transferred | Fails | 3.517 |
| J1204+0358 transferred | No violation found | 0.000683 |
| J1402+6321 transferred | Fails | 8.304 |
| J1621+3931 transferred | Fails | 8.514 |
| J1630+4520 transferred | Fails | 5.682 |

For monotone cases the listed peak is the innermost sampled radius, not a measured central core radius. The extended J1621 branches peak near 150, 689 and 3526 kpc. Refined log-density slopes are strongly positive below those peaks, so this is not merely roundoff near a flat center. A density peak need not enclose half the mass; the full results record sampled mass fractions without interpreting them as observed inventories.

## Numerical checks

The coarse grid uses 4097 logarithmic radii over 1e-4<=r/ac<=1000 and 192 incoming angles. Refinement doubles both to 8193 and 384. A flagged interval requires dln(rho)/dln(r)>1e-3 and rho greater than 1e-12 of its peak. All classifications agree after refinement. Log-sum-exp prevents extreme attenuation from causing numerical underflow. The analytic central criterion independently agrees with every classification. These are necessary-condition violations, not complete distribution-function inversions. No-violation cases remain unproved.

## Independent-observation search and its consequence

[Gavazzi et al. (2007), SLACS IV](https://arxiv.org/pdf/astro-ph/0701589), Table 1, lists 22 lenses and does not include J1621+3931. Its weak-lensing result is an ensemble measurement. The published source weighting and cosmological geometry also differ from this project's conditional single-source predictions. It therefore cannot be assigned to J1621 as an individual measured outer profile. J1402+6321 and J1630+4520 overlap our six-system sample, but a stacked result is not either object's own shear curve.

[Treu et al. (2009), SLACS VIII](https://arxiv.org/pdf/0806.1056), studies neighboring-galaxy density and environment. Such counts do not directly measure the proposed companion convergence or its mass. This targeted search did not yield an individual outer-shear likelihood for J1621; it does not establish that none exists. No companion branch was accepted or rejected using those papers' ensemble/environment summaries. The support calculation above supplies a separate mathematical constraint while direct outer comparison remains open.

## Concrete alternatives

1. Angular-momentum-supported particles: use f(epsilon,L) or an orbit superposition that avoids excessive central crossing. An ideal collection of randomly oriented circular orbits can maintain a shell-like density kinematically; a formation mechanism and collective stability must still be demonstrated.
2. Interacting or wave-supported deposits: derive radial and tangential stresses from actual state variables, rather than inserting whichever pressure reproduces the chosen density. Include energy cost and escaping radiation.
3. Time-dependent reservoirs: evolve capture, inward movement and release. A snapshot can be hollow without equilibrium, but its lifetime and observable motions must be calculated.
4. Revise opacity/capture shape: the derived hollow-center threshold identifies which parameter combination causes the problem. Any revision must preserve or re-evaluate the successful motion/lensing limits and shared-parameter transfer.

The most direct next particle calculation is whether a tangentially supported construction meets the required angular momentum and stability conditions without an unaccounted energy source. No isotropic support is adopted by default. The existing fits and one-third reference are preserved as phenomenological profiles, with this newly established limitation attached.

## Reproduction

Run reservoir-isotropic-support.py. The protocol, input hashes, coarse/refined outputs and analytic threshold are retained in this directory. No data or fit values were altered. This report extends the paper supplement beyond PDF v1.3.
