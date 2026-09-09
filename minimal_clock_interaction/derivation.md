# Minimal atom–wave–cavity interaction test

This is an analytic exploration with numerical consistency checks, not a new observational fit. The previous best history p=0.2639061166 and gamma=7.7315e-11/year remain fixed. No complete microscopic or gravitational theory is established.

## A constraint independent of atomic details

Assume homogeneous, nondispersive, adiabatic propagation with frequency omega=ck/n(t), conserved spatial wave number k, and a cavity with fixed longitudinal mode number j and instantaneous length L(t). Let nu_A(t) denote a chosen atomic transition; it need not be constant. A wave emitted on that transition has observed spectroscopic redshift S=(n_o/n_e)(nu_A,o/nu_A,e). The cavity frequency is nu_C=j c/(2 n L). If nu_C/nu_A is exactly constant throughout the trajectory, n L nu_A is constant and therefore S=L_e/L_o.

Consequently, if material lengths also remain fixed, S=1. If S>1 while the cavity/atom ratio is protected, those lengths must decrease. For fixed geometric source separation R, the dimensionless separation R/L then increases by S. This is expansion measured in those material rod units, even if an auxiliary geometric metric is static. The identity does not prove the full theory is equivalent to standard cosmology. It does show that exact cavity protection cannot be combined with all of the fixed-material-length assumptions used in the previous brightness derivation.

Exact protection over the entire history is stronger than an observationally small present-day drift. A coupling with a stationary response today could evade a present-day derivative constraint without satisfying this identity across past epochs; its higher derivatives and past atomic spectra would then need testing. No such coupling is claimed here.

## Explicit leading-order interaction family

Use a preferred-frame electromagnetic Lagrangian density L_EM = epsilon_0 n^a E^2/2 - B^2/(2 mu_0 n^(2-a)), plus minimally coupled nonrelativistic charged particles with kinetic masses m_e=m_e0 n^b and m_N=m_N0 n^b. Hold charges, hbar, mass ratios and Pauli g factors fixed. The spin magnetic moments are g q hbar/(4m). This specifies a Coulomb/Pauli toy interaction; it is not a relativistic, radiatively protected theory. n(t) is a prescribed background in this test.

Permittivity times permeability is proportional to n^2, so waves have instantaneous speed c/n. The Coulomb interaction scales as e^2/(epsilon r). Leading hydrogenic scalings give:

- Atomic length a_B proportional to epsilon/m_e, hence n^(a-b).
- Optical frequency proportional to m_e/epsilon^2, hence n^(b-2a).
- Hyperfine frequency proportional to mu times the two magnetic moments times a_B^(-3), hence n^(2+b-4a).
- Cavity frequency, assuming its spacer tracks a_B, proportional to n^(b-a-1).

Thus the cavity/optical ratio scales as n^(a-1), while the hyperfine/optical ratio scales as n^(2-2a). Both are constant at this leading order only for a=1, independent of b. The observed redshift then scales as S=(n_o/n_e)^(b-1), and the rod length scales as n^(1-b), reproducing S=L_e/L_o.

| Choice | Atomic frequency | Material length | Observed redshift |
|---|---|---|---|
| a=1,b=1 | n^(-1) | Constant | Cancels |
| a=1,b=2 | Constant | n^(-1) | n_o/n_e |

The second row is a constructive leading-order compensation, but pays for it with changing rods. It does not preserve the fixed-material behavior assumed by our original signal model. Relativistic corrections, nuclear structure, radiative shifts, and actual solid-state spacer response have not been derived. With an independently fixed matter limiting speed, even the optical relativistic corrections can introduce new variation. Therefore neither row proves exact protection of real optical/hyperfine clock ratios.

## Observable tests and scope

For fixed atomic frequencies and fixed cavity length, the universal index model predicts d ln(nu_C/nu_A)/dt = -dot n/n. Today this is -7.7315e-11/year regardless of p. The history fit does not reduce this local prediction. This is a forecast, not a measured drift or an exclusion significance.

Kennedy et al., arXiv:2008.08773, compare a silicon cavity, an optical atomic clock and a hydrogen maser, demonstrating the relevance of these three different responses. Their oscillatory-field analysis is not a direct secular-drift likelihood for our model; the earlier paper review noted removal of overall linear drift. We have not fitted their raw frequency data here. Source: https://arxiv.org/abs/2008.08773 . Their physical sensitivity formulas also include corrections absent from this leading toy.

For stationary source and observer in the fixed-atom history branch, observer-time redshift drift today is dz/dt=gamma[(1+z)^(1-p)-(1+z)]. At z=1 it is -2.58495e-11/year. This remains a conditional forecast, not a new detection.

## What to prioritize

The history formula is an improved empirical candidate. This derivation improves the theoretical specification by exposing which assumptions compete. The next narrow task is to establish what present-day cavity data can actually constrain about the secular drift while testing atomic responses from the same interaction. Exact zero cavity drift should not be substituted for an actual likelihood.

If fixed local matter and protection are required, an interaction beyond a universal homogeneous index must specify why free and cavity waves respond differently. An extra propagation mode or a boundary-dependent response is a possible research direction, not yet a solution. Simply declaring that bound photons are exempt is not an interaction, and ordinary frequency selectivity alone cannot distinguish a cavity photon from a traveling photon at the same frequency. Screening is not reintroduced here.

If changing rods is accepted as an alternate-universe candidate, its brightness, angular-distance and matter dynamics must be rederived together; it cannot inherit the previous fit unchanged. This branch is documented, not adopted.

No general impossibility theorem has been proved for dispersive, nonlocal, environment-dependent, multimode or nonadiabatic theories. No curvature, rotation, CMB or new gravity sector was introduced in this calculation.

## Reproduction

Run check.py using Python 3; no extra packages are required. checks.json contains a nine-case exponent scan, exact cavity/redshift identity checks and redshift-drift forecasts. These are synthetic checks of algebra, not real-data validation. The maximum numerical error in the identity check was zero for the tested values.
