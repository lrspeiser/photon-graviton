# Tangential circular-orbit support of the fitted reservoirs

13 September 2026. Explicit optional deposited-particle completion; no new fits or observations.

## Outcome

A spherical ensemble of circular orbits, with random orbital planes/phases and paired directions, can kinematically maintain the tested hollow densities. All ten nonzero recorded profiles admit positive circular speeds and pass the individual radial stability criterion in their frozen stellar-plus-deposit Newtonian potentials. This avoids the isotropic-particle monotonicity failure without changing the density profile.

It is a singular, maximally tangential construction, not a smooth isotropic velocity distribution. It does not show that particles can form from incident companion waves, arrive on these orbits, or remain collectively stable when the whole reservoir responds. The calculated kinetic energy is small relative to the reference mass-energy inventory but not zero; for the largest J1621 branch it is about 1.05e60 J, or 6.97e-5 of M_d c^2.

## Established construction and derivation

The circular-orbit ensemble is a known Einstein-cluster type of construction, used here only in its Newtonian leading-order form. It is not new companion physics. See [Boehmer and Harko, On Einstein clusters as galactic dark matter halos](https://arxiv.org/abs/0705.1756) and [Modelling Einstein cluster using Einasto profile](https://arxiv.org/abs/2311.18622). Stability results for other specified clusters are not transferred to our fitted profiles.

At radius r, put the deposited mass element dM_d on circular orbits with

    v_c^2(r)=r Phi'(r)=G M_total(<r)/r,
    <v_r^2>=0,  <v_theta^2>=<v_phi^2>=v_c^2/2,
    p_r=0,  p_theta=p_phi=rho_d v_c^2/2.

The spherical stationary stress equation is

    dp_r/dr + (2p_r-p_theta-p_phi)/r = -rho_d Phi'(r).

Substitution closes this equation exactly. Each orbit stays at a fixed radius; uniform phases and orbital-plane directions preserve the prescribed spherical density. These are standard mechanics and moment equations, applied to the hypothetical density. Closing the moment equation alone is not proof of formation or collective stability; the ideal orbit ensemble supplies a stationary measure rather than a regular finite-width distribution function.

For a small radial displacement of an individual particle at fixed specific angular momentum j, the effective potential is Phi+j^2/(2r^2). At its circular orbit,

    j_c^2=G M_total(<r) r,
    omega_r^2=Phi''+3Phi'/r=G[M_total+r M_total']/r^3,
    omega_r^2/Omega^2=1+dln(M_total)/dln(r) > 0.

Positive density makes this positive in the adopted potential. A particle then oscillates after a small radial displacement rather than immediately running away. This is frozen-potential radial stability only. It does not include reservoir self-response, nonspherical modes, counterstream interactions or collective relaxation. The frequency identity is established epicyclic mechanics, not a new force law.

## Stored motion and angular momentum

    K_d = (1/2) integral v_c^2 dM_d,
    L_absolute = integral r v_c dM_d,
    L_vector = 0 for paired orbital directions.

Zero net rotation does not mean individual angular momentum vanishes. The result JSON records the absolute inventory and its mass-weighted specific value. Orbital motion need not consume continuous power in a collisionless stationary potential; feeding is required for losses or further growth, not simply because particles keep moving.

| Profile | Maximum sampled circular speed (km/s) | Kinetic energy inside grid (J) | K_d/(M_d c^2) |
|---|---:|---:|---:|
| J0037-0942 free | 429.6 | 5.97e+53 | 7.61e-07 |
| J1112+0826 free | 569.0 | 1.77e+54 | 1.3e-06 |
| J1204+0358 free | 506.9 | 1.79e+52 | 1.07e-06 |
| J1402+6321 free | 635.9 | 3.18e+54 | 1.61e-06 |
| J1621+3931 free | 4245.8 | 1.05e+60 | 6.97e-05 |
| J1630+4520 free | 509.0 | 2.3e+52 | 1.18e-06 |
| J1621 scale 0.3 | 392.4 | 6.06e+51 | 7.39e-07 |
| J1621 scale 10 | 885.4 | 9.4e+55 | 3.03e-06 |
| J1621 scale 30 | 1873.6 | 8.06e+57 | 1.35e-05 |
| J1621 scale 100 | 4245.8 | 1.05e+60 | 6.97e-05 |

The extreme J1621 branch has an orbital period near its density peak of about 9.81 billion years; the scale-10 and scale-30 branches give about 1.91 and 4.34 billion years. These are computed orbital periods, not imposed universe ages, formation times or evidence that an equilibrium has been reached. Monotonic compact profiles peak at the innermost numerical point, so their recorded peak periods have no measured-core interpretation.

These calculations identify the existing effective density with inertial mass only to leading nonrelativistic order. If that mass denotes rest mass, kinetic and binding energy must enter the total gravitational accounting; if it already denotes total energy-equivalent mass, the rest-mass assignment must be adjusted. The kinetic ratios quantify one correction rather than silently adding it to the fitted lens source. Kinetic energy can arise through settling and binding-energy release; K is not automatically an additional external photon requirement on top of the full binding ledger. No relativistic lensing correction is claimed to have been calculated here.

## Capture remains a separate problem

Under ordinary local energy-momentum kinematics, one isolated null packet cannot become a single massive particle in vacuum: its invariant mass is zero. Capture therefore requires interacting packets, a recoiling reservoir, another outgoing channel or changed microscopic laws with an explicit conservation rule. Oppositely directed incident packets can have nonzero combined invariant mass, but they still need an interaction and the correct orbital momentum distribution.

For an illustrative unattenuated isotropic flux crossing a sphere, let mu be the inward radial direction cosine. The flux-weighted density is p(mu)=2mu, with 0<=mu<=1. A null packet assigned E approximately m c^2 has specific angular momentum j_in approximately r c sqrt(1-mu^2). The fraction with j_in<=j_c is then (v_c/c)^2. At the density peak this is approximately 2.58e-6, 1.06e-5 and 5.42e-5 for the extended scale-10, 30 and 100 branches.

This fraction is not capture efficiency and is not the actual attenuated angular field. Even a packet below that threshold does not automatically circularize; exactly matched angular momentum is a sharper condition. The estimate demonstrates why arbitrary incoming rays cannot be assigned unchanged to slow circular orbits. Momentum must be exchanged among packets, matter, outgoing waves or the reservoir; cancellation of the vector total does not describe that process or account for its energy.

## Finite-domain bounds and verification

The refinement uses 16001 radial mass points and 384 incoming angles versus 8001/192. Maximum relative changes in stored mass and K are below 1.44e-6 and 1.28e-6. A finite-difference derivative of j_c^2 independently checks the radial-frequency formula; the maximum refined relative difference is 1.64e-4 and all frequencies remain positive. The inherited stellar grid is unchanged.

All integrated inventories refer to the stated finite numerical domain. With rho_d<=D ac^4/r^4 and an upper total mass M_u=M_stars+M_d(<R)+4*pi*D*ac^4/R, exterior bounds are

    K_d(>R) <= pi G M_u D ac^4/R^2,
    L_absolute(>R) <= 8*pi D ac^4 sqrt(G M_u/R).

These follow by elementary integration and bound only the tail, not model validity. The largest kinetic tail is at most 0.0441% of the interior K; the absolute-angular-momentum tail can still reach 17.3% of its interior value. Small mass-tail error therefore does not imply equally small angular-momentum-tail error. The kinetic inventory remains well bounded while the angular inventory requires the stated qualification.

## Consequence and next distinguishing calculation

Tangential support is a mathematically viable stationary alternative to the excluded isotropic completion at this idealized level. It preserves the existing density and leading-order lens/motion predictions, but does not derive them. A physical model must next show how conversion and interactions produce the necessary distribution of angular momentum and binding energy, and test perturbations allowing the deposited density itself to move. A broadened orbit distribution or explicit wave stress may be more realistic than exact circular shells; either must preserve the hollow profile for a quantitatively adequate lifetime.

The one-third law remains the empirical reference. No density parameters were fitted, no outer observations were compared and no claim of full stability or successful capture is made. The paper supplement records this result after the v1.3 PDF snapshot.

## Reproduction

Run circular-reservoir-support.py and repeat with --refine. Both JSON outputs preserve the ten profiles, domains, energy and angular-momentum inventories, individual-orbit checks and input provenance. Choices are recorded in circular-reservoir-support-protocol.md.
