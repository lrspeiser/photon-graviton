# Distant stellar supply with competing capture

## Result

A homogeneous comparison calculation shows how an unlimited spatial population of emitters can deliver finite power to each receiver without imposing a universe size. Other receivers capture companions en route. This preserves the user's external-source argument while preventing the same energy from being allocated to every cluster.

For positive conversion and capture coefficients, the long-time radiation bath can approach a finite state even though permanent deposits continue growing. The calculation does not determine actual luminosity density, capture abundance, cosmic source history or the gravity of the resulting reservoir. It establishes neither observed cluster sufficiency nor shortage.

## Postulates and known transport mathematics

Use Euclidean static geometry, a homogeneous luminosity density J (power/volume), photon and companion speed c, fractional conversion alpha per length, and effective companion removal beta per length. Companions lose no energy during free propagation. Removal transfers their energy to deposits elsewhere; it is not destruction. Assume independent, dilute randomly distributed receivers, stationary coefficients and no photon absorption or re-emission. Homogeneity neglects correlations between stars and wells and does not describe the detailed internal cluster kernel.

For a unit emitted photon-energy pulse, the established linear transfer equations are

    dP/dD = -alpha P
    dC/dD = alpha P - beta C
    dE_deposited/dD = beta C.

Starting with P=1, C=0, their solution is

    P(D)=exp(-alpha D)
    C(D)=alpha/[beta-alpha] * [exp(-alpha D)-exp(-beta D)]
    E_deposited(D)=1-P(D)-C(D).

At beta=alpha, C(D)=alpha D exp(-alpha D). At beta=0, C(D)=1-exp(-alpha D). The surviving companion term is NOT generally [1-exp(-alpha D)] exp(-beta D): companions can be generated late along the path, after other companions have been captured. Conversion, propagation and capture must be solved together.

These equations are known coupled transport mathematics applied to optional companion postulates, not a unique microscopic derivation.

## Integrating the stellar shells

Let a particular receiver have effective absorption cross-section sigma. Distant shells emit 4 pi D^2 J dD and the receiver intercepts sigma/(4 pi D^2). Thus

    deposited power from shell = sigma J C(D) dD.

For a source population illuminated for a finite duration T, with H=cT in this uniform simultaneous-start comparison,

    P_receiver(T)=sigma J integral_0^H C(D) dD.

For alpha,beta>0 and unlimited established illumination,

    integral_0^infinity C(D) dD = 1/beta
    P_receiver(infinity)=sigma J/beta.

This is the infinite homogeneous medium limit, using coarse-grained transport down to small distances; a real excluded receiver volume, correlated local emitters and near-field geometry need matching to the detailed cluster calculation. It is not an exact full-catalog cluster prediction. Independence of alpha in the steady limit does not mean conversion is irrelevant: the approach requires both alpha c T and beta c T to be sufficiently large. Zero alpha gives zero supply and does not commute with an infinite-time limit.

## Capture abundance fixes the competition

For receiver classes i with number densities n_i and cross-sections sigma_i, the dilute random-encounter approximation gives

    beta = sum_i n_i sigma_i
    P_i = sigma_i J/beta
    sum_i n_i P_i = J.

For identical receivers, P=J/n. More cross-section increases capture but also shortens the distance companions travel before another receiver captures them. Cross-section and free-travel distance therefore cannot be adjusted independently in this homogeneous model. For a mixed population the allocated fraction is n_i sigma_i/beta.

For a uniform spherical receiver, the known straight-chord absorption average gives

    sigma=pi R^2 A(kappa R)
    A(tau)=integral_0^1 2y [1-exp(-2 tau y)] dy
          =1-[1-(1+2 tau)exp(-2 tau)]/(2 tau^2).

At tau=0.1,1,10, A=0.12385,0.70300,0.99500. These are effective companion absorption cross-sections, not measurements of real cluster opacity. An illustrative two-class allocation is included only to check that total allocated power is J; its densities and radii are not observational estimates.

## Finite histories and energy conservation

For an initially empty homogeneous medium with constant J switched on at T=0,

    du_gamma/dT=J-alpha c u_gamma
    du_companion/dT=alpha c u_gamma-beta c u_companion
    du_deposit/dT=beta c u_companion.

Consequently u_gamma+u_companion+u_deposit=JT. At late times, u_gamma approaches J/(alpha c), u_companion approaches J/(beta c), and deposited energy grows at rate J. Finite radiation energy is not a stationary total universe. Stellar fuel must supply the emitted energy; this model does not authorize an eternal constant stellar luminosity or generate energy by recycling it while keeping deposits.

The run retains 25 finite-history cases with alpha c T=0.1,1,10,100,1000 and beta/alpha=0,0.01,0.1,1,10. These dimensionless histories are not chosen cosmic ages. Direct shell integration agrees with companion bath energy, and an independent augmented matrix exponential verifies photon/companion/deposit evolution. Maximum normalized matrix difference is 5.69e-16. Energy ledgers, nonnegative energies and three independent cross-section integrals pass. These numerical tests do not establish physical capture microphysics or observational agreement.

## Next constraints

Measure or bound J with a consistent brightness law; estimate receiver abundances and cross-sections from the same capture model; include real source histories, spatial correlations and intervening capture; then compare predicted photon backgrounds and deposited fields with observations. The stationary photon energy J/(alpha c) provides an additional consistency check, not permission to hide an arbitrarily bright stellar source population. Its spectrum and event timing remain separate unresolved requirements.

The original finite catalog is still a partial inventory, not an upper bound on supply. This result supplies a conservation-consistent framework for extending it. All six objectives remain open and no final holdouts were opened.

Run `python research_work/results/cluster-competing-capture/run.py`.
