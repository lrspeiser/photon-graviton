# The clock candidate needs a physical companion source

**Result:** the existing prescribed lapse can redshift and stretch an event with no source of a separate companion stress tensor. Its optical success therefore does not calculate companion production. The previous two-channel report has been corrected where it described the time-associated endpoint energy difference as necessarily delivered work.

This is a derivation within an existing candidate, not another data-resolution or numerical-calibration exercise. It identifies exactly which extra equation must connect the redshift calculation to the cluster capture calculation. All six objectives remain open.

## Known identities and candidate assumptions

Use c=1 and signature (-,+,+,+). Four-velocities satisfy u.u=-1. In the existing hypothetical clock prescription,

    ds_metric^2 = -dt^2/n(x,t)^2 + dx^2 + dy^2 + dz^2
    u = n partial_t
    d tau = dt/n.

This is established lapse geometry; interpreting a particular n as a radiation-generated void field is our unproved physical proposal. Standard covariant stress conservation and observer projections are reviewed in [Gourgoulhon's 3+1 notes](https://arxiv.org/abs/gr-qc/0703035) and [Tong's treatment of matter and energy conservation](https://davidtong.org/teaching/general-relativity/grhtml/S4). The following product-rule derivation does not assume expansion, a cosmic age, dark matter, or Einstein's equations for the still-unspecified source of n.

Define a physical transfer four-vector Q_gamma by

    nabla_mu T_gamma^(mu nu) = Q_gamma^nu
    J_gamma^mu = -T_gamma^(mu nu) u_nu.

The exact product rule gives

    nabla_mu J_gamma^mu
      = -u_nu Q_gamma^nu - T_gamma^(mu nu) nabla_mu u_nu.

The first term is local energy supplied to photons by another sector in this observer frame. The second changes the energy current because the observers used to measure it vary through spacetime. It need not be a new particle source. A nonzero second term does not imply a nonzero Q_gamma.

For these fixed-position observers the spatial metric is time independent: expansion and shear vanish, even if n changes. Their acceleration is a_i=-partial_i ln(n). Decompose the radiation stress into local energy density, flux F and spatial stress. Contracting with nabla u yields

    T_gamma^(mu nu) nabla_mu u_nu = F dot a
    nabla_mu J_gamma^mu = -u dot Q_gamma - F dot a.

This agrees with the previous directional result: a geodesic photon measured locally obeys d ln(E)/d ell = -a dot e = e dot grad ln(n), where e is its unit spatial direction and ell is path length. Thus E can change while Q_gamma=0. For locally isotropic radiation F=0, this observer term vanishes at that event. This does not prove zero integrated redshift on different evolving ray histories, nor exclude all gravitational energy exchange descriptions.

## Apply the result to the retained optical candidate

The prescribed test profile n=1+beta*t*sin^2(pi*x), with endpoints x=0 and x=1 and c=1, already gives the exact event and wavelength factor S=exp(beta/2). Its endpoint clocks both have n=1. Treating its photons as freely propagating geometric-optics radiation gives Q_gamma=0 throughout: the Hamiltonian contains no conversion collision term or separate emitted species.

For beta=0.2, that known solution has

    1+z = S = exp(0.1) = 1.105170...
    E_observed/E_emitted = exp(-0.1) = 0.904837...
    separate companion production from the stated equations = 0.

These values follow from the already solved profile; no new simulation or observational pass is claimed. The geometry remains prescribed and lacks a demonstrated self-consistent source. This example establishes a limited but decisive point: even equal endpoint clock rates and persistent measured redshift do not determine a separate companion-production current.

In curved spacetime, covariant matter conservation is not a general global scalar-energy conservation law. A complete model must state its field dynamics, boundary fluxes and conserved energy where such an energy exists. We are not claiming that gravity cannot exchange energy, or that the measured redshift is unreal. We are rejecting the shortcut from a measured energy ratio to a local population of lossless companions.

## What an explicit receiver would change

**Proposed sector structure, using known interacting-stress mathematics:** introduce companion production C^nu and capture K^nu. A simple three-sector convention is

    nabla T_gamma = -C
    nabla T_companion = C - K
    nabla T_deposit = K.

These signs make the displayed total stress divergence zero. If the independent time field exchanges energy or momentum too, it requires its own tensor and exchanges; it cannot be omitted from that sum. For future-directed C and the local observer u, the production power is P_companion=-u.C >= 0. The deposited source is P_capture=-u.K, not alpha times the observed photon energy by definition.

For example, a collinear kinetic drag may be postulated as Dp^mu/dlambda = -eta E p^mu, where E=-u.p and eta>=0 has inverse-length units. Its photon energy loss has the local form

    d ln(E)/d ell = -a dot e - eta.

The geodesic path remains unchanged under this parallel drag; only its momentum scale changes. This is a kinetic ansatz, not a coherent Maxwell interaction or a microscopic graviton amplitude. It adds genuine transfer but does not stretch event spacing beyond the existing ray/clock map. Consequently, if A_metric=ln(S) and A_conversion=integral eta d ell, it gives

    ln(1+z_total) = A_metric + A_conversion
    ln(S) = A_metric.

This reproduces the restricted channel-budget relation b=A_metric/[A_metric+A_conversion]. The prior observational comparison favors b near one under its stated source assumptions. Adding large drag to a metric already fitted to the full redshift therefore overproduces color change relative to timing. The same metric redshift cannot be counted once as metric stretching and again as this extra drag.

## Concrete implication for the research program

The simple additive repair is insufficient: a near-zero conversion share preserves timing but provides no demonstrated large companion source; a large conversion share changes the observed redshift/timing relation. This is not an energy-supply exclusion for the whole theory, because unrestricted histories, coupling strength and other mechanisms remain open.

A viable completion must generate the phase/arrival map and the nonzero receiving current together from one interaction, with actual clock standards included. The earlier coherent-receiver model establishes that energy-conserving coherent frequency conversion can be constructed, but does not yet provide spatial event transport. The lapse model provides event transport but no separate companion current. Neither can borrow the other's successful result without deriving their combined dynamics. The next model must therefore change the interaction itself, rather than append a deposit ledger to unchanged geodesics.
