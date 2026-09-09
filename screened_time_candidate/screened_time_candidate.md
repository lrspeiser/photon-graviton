# A screened temporal field: a concrete candidate and its first tests

This is a speculative construction for the fictional, nonexpanding universe. Its purpose is a physical law governing observable timing and redshift, with a possible later connection to the CMB and gravity. The new numerical work below tests a mathematical subsystem, not astronomical observations. No new observational validation is claimed.

## 1. What we mean by a new law for time

Introduce a dimensionless field χ and a positive timing factor

\[n(x,t)=e^{\chi(x,t)}.\]

In a local region where ordinary gravitational effects can be neglected, let light obey

\[\omega(x,k,t)=c|k|e^{-\chi(x,t)},\qquad d\ell/dt=c/n.\]

One optical-metric representation is

\[ds_\gamma^2=-c^2e^{-2\chi}dt^2+d\mathbf{x}^2.\]

The corresponding optical lapse is dτγ=e^{-χ}dt. This labels a proposed propagation law; photons do not carry proper-time clocks. Material clocks retain their ordinary gravitational behavior at this stage, with dτmatter≈(1+Φ/c²)dt for stationary weak-field clocks in conventional coordinates. A full covariant model must specify both sectors and their coupling.

The observable hypothesis is therefore a difference between the timing of light propagation and material clocks. If we change every clock and propagation law identically through a homogeneous lapse, the change can be a coordinate relabeling. Calling that “new time” would not produce a measured redshift. Here χ changes a measurable relation. Whether this deserves the name a universal theory of time remains a question for the completed theory.

Distance does not itself act locally on time. Light accumulates the effect of the field along its trajectory. Hamilton's equations give, for this nondispersive local model,

\[\frac{d\ln\omega}{dt}=-\partial_t\chi,\qquad
\ln S=\int_{\rm ray}\partial_t\chi\,dt
=\frac1c\int_{\rm ray}\partial_t n\,d\ell.\]

S is the propagation frequency ratio ωem/ωobs. Identifying it with measured 1+z additionally requires identical source and observer atomic standards, and treatment of motion and ordinary gravity. Infinitesimal wave-packet intervals stretch by the same S under these assumptions. Finite-duration events require that the field change negligibly over the event.

A static spatial profile alone does not accumulate this propagation redshift. It can refract rays and alter travel times. Ordinary gravitational endpoint shifts must be included separately.

## 2. A conservative subsystem that permits rolling and suppression near matter

Our previous single, pinned, undriven field oscillated instead of maintaining a monotonic drift. Try a second field ψ that carries finite kinetic energy. For a fixed homogeneous matter density ρ, take

\[\mathcal L_{\rm fields}=\frac K2\dot\chi^2+\frac{K}{2\epsilon}\dot\psi^2
-\frac K2\omega_0^2(\chi-\psi)^2
-\frac K2\omega_\rho^2\chi^2.\]

Here K>0, ε>0, and ωρ² increases with density; for example ωρ²=αρ with α>0 in consistent units. Spatial gradient terms with positive energy are needed in the global version. This is an illustrative effective model, not the established chameleon model. Density-dependent suppression has precedents in scalar-field research [1], but those papers do not validate this construction or a static cosmology.

The local equations are

\[\ddot\chi+\omega_0^2(\chi-\psi)+\omega_\rho^2\chi=0,\]
\[\ddot\psi+\epsilon\omega_0^2(\psi-\chi)=0.\]

In an ideal empty patch, ωρ=0, the exact solution χ=ψ=γt has a freely rolling common mode. Its kinetic energy is finite per unit volume and constant. A continuous external drive is unnecessary for this isolated field solution. It still requires an initial velocity and does not explain why γ has its fitted value or the same sign across the universe.

In a dense patch, the slowly varying equilibrium is

\[\chi\simeq s(\rho)\psi,\qquad
s(\rho)=\frac{\omega_0^2}{\omega_0^2+\omega_\rho^2}.\]

This gives physical suppression in the equations rather than an imposed mask. ψ experiences a restoring force there: it is not an inexhaustible source of constant drift. Small ε gives it greater inertia and slows the response.

### Numerical result

Use time u=γt, ω0/γ=10, ε=10⁻⁴. Compare an empty patch to a dense patch with ωρ/ω0=1,000. Initialize χ=ψ=0 and dimensionless velocities χ′=s, ψ′=1. Solve by diagonalizing the symmetric, mass-weighted oscillator matrix. These parameter choices illustrate a mechanism; they are not measured screening parameters.

| Quantity at u=1 | Empty patch | Dense patch |
|---|---:|---:|
| χ | 1.000000 | 0.000000998333 |
| ψ | 1.000000 | 0.998334 |
| χ/ψ | 1 | 0.00000100000 |
| ψ′ | 1 | 0.995004 |
| Maximum relative field-energy error over the interval | below reported precision | 5.6×10⁻¹⁶ |

With the original local redshift slope, u=1 represents 12.934 billion years. The dense homogeneous patch suppresses χ by approximately a million relative to ψ while the latter's velocity changes by about 0.5%. This is mathematical feasibility within the subsystem, not evidence of actual screening on Earth.

All local mode frequencies squared are nonnegative. This establishes linear stability of these fixed-density homogeneous field oscillators. It does not establish stability across boundaries, in gravity, or after coupling to radiation and moving matter. In particular, joining rolling voids to screened galaxies introduces gradients that can remove or modify the common rolling mode. Initial conditions were chosen to favor the slowly varying solution; their natural production has not been demonstrated.

## 3. The redshift formula must come from the dynamics

In the homogeneous rolling limit, set χ=γ(t−tobs), so nobs=1. Then

\[R=\int_{t_{em}}^{t_{obs}}ce^{-\gamma(t-t_{obs})}dt
=\frac c\gamma(S-1),\]

and therefore

\[\boxed{1+z=1+\frac\gamma c R.}\]

This is exact for the homogeneous propagation model with fixed endpoint clock standards. It is not yet the global law for rays leaving and entering screened regions. It shares the first-order slope of our original relation

\[1+z=e^{\kappa R},\qquad \gamma=c\kappa,\]

but differs beyond first order. Using κ=0.000077314965955 per million light-years:

| Geometric distance, million light-years | Original exponential z | Homogeneous rolling-field z |
|---|---:|---:|
| 1 | 0.000077318 | 0.000077315 |
| 100 | 0.00776146 | 0.00773150 |
| 1,000 | 0.0803823 | 0.0773150 |
| 10,000 | 1.16658 | 0.773150 |

These are model calculations at the same slope, not fitted comparisons to distant observations. The earlier nearby-group exercise already examined a related linear-in-geometric-distance law and could not discriminate decisively. Its brightness conversion must not be imported automatically into this screened optical model; focusing, flux and observer clocks need to be derived again.

If we insist on the exact original exponential law in a homogeneous medium, an affine optical factor n=1+γ(t−tobs) supplies it. With n=e^χ, this requires

\[\chi=\ln[1+\gamma(t-t_{obs})].\]

For a single canonical homogeneous scalar with no radiation exchange, reconstructing a potential to support this trajectory gives

\[V(\chi)=V_0-\frac{K\gamma^2}{2}e^{-2\chi}.\]

It has V″/K=−2γ²e⁻²χ: negative curvature and a potential unbounded below as χ decreases. The affine n also reaches zero at a finite past coordinate time, 1/γ≈12.934 billion years before the reference epoch. This is a pathology of that simple realization, not evidence of a Big Bang and not a general impossibility proof for every completion. It argues against forcing the exact exponential fit into this particular canonical field model.

## 4. Energy and the CMB

“Time uses energy” becomes testable only after specifying what carries energy and exchanges it. In this construction χ and ψ carry energy. Time itself has not been assigned an independently measured energy density.

The two-field energy conserved in the calculation is

\[E_{fields}=\frac K2\dot\chi^2+\frac K{2\epsilon}\dot\psi^2
+\frac K2\omega_0^2(\chi-\psi)^2+\frac K2\omega_\rho^2\chi^2.\]

Once light is included, its lost energy must enter the fields, matter, or boundary fluxes. The equations above omit that backreaction. A full action or Hamiltonian must supply it; writing Q with opposite signs in two continuity equations alone is not a microscopic derivation. K is unconstrained here, so we have not shown that the reservoir has enough capacity or that its gravitational effect is acceptable.

There is a useful conditional CMB result. In a homogeneous medium, conserved photon wavevectors and conserved mode occupation turn an initially Planckian distribution into another Planckian distribution when every mode frequency is scaled by the same 1/n:

\[T(t)=T_i\frac{n_i}{n(t)}.\]

No expansion or Big Bang is required for this mathematical cooling operation. However, it preserves an existing thermal distribution; it does not explain its origin or select 2.72548 K. An initial temperature and field history could simply tune that number.

Also, the density of photon states changes in this optical model. With ω=ck/n in a fixed physical volume, the mode-counting result is uγ∝n³T⁴, hence uγ∝1/n for this homogeneous adiabatic evolution. It would be inconsistent to import the usual vacuum u∝T⁴ with a fixed coefficient while n changes. Locally screened instruments and propagation through boundaries require an explicit radiative-transfer calculation before comparison to the observed CMB intensity. Work on nonuniversal radiation couplings likewise warns that distribution functions and measured radiation properties can change [2].

Our earlier FIRAS residual-table test remains a useful target: equal-weight mixtures of received thermal temperatures separated symmetrically by ±1% fitted much worse than one temperature (approximate diagonal-error χ² 235.84 versus 45.02). That restricted test is not a confidence bound for this field theory. It says that a global implementation cannot freely mix substantially different received temperatures and still assume a blackbody.

The CMB anisotropies, angular peaks, polarization and proposed preferred-axis alignment remain unpredicted. A sphere or other topology could supply boundary conditions, but a preferred axis must emerge from specified conditions and then predict additional observations. Fitting an axis after inspecting the sky would not establish the mechanism.

## 5. Next milestone and success criteria

| Goal | What would count as progress | Present status |
|---|---|---|
| Physical screening | Solve the field across realistic dense objects and voids; predict suppression without setting χ=0 by hand | Homogeneous mechanism works; spatial problem remains |
| Observable time law | Derive clocks, frequency shifts and pulse stretching from one matter–light–field theory | Optical timing law specified; atomic completion missing |
| Sustained evolution | Rolling solution survives boundaries and backreaction over the required epoch | Stable isolated patches only |
| Energy conservation | Derive photon-to-field exchange and determine its effect on drift and gravity | Field subsystem conserves energy; full system not closed |
| Redshift prediction | Derive the distance and brightness laws, freeze parameters, test unused data | Two homogeneous laws now distinguishable in principle |
| CMB | Reproduce spectrum and intensity at screened detectors, then anisotropy and polarization | Thermal-shape preservation only in homogeneous limit |
| Galaxy dynamics | Derive stellar forces and lensing from the same fields | No new rotation-curve prediction |

The most useful next calculation is a spatially resolved void bounded by dense regions, with both fields evolving and photons exchanging energy with them. It should determine whether physical screening and accumulated stretching actually coexist in a connected system. Only then should we fit a global redshift law or claim a CMB prediction.

## Sources and reproducibility

1. Khoury & Weltman, *Chameleon Fields: Awaiting Surprises for Tests of Gravity in Space*, https://arxiv.org/abs/astro-ph/0309300. Precedent for density-dependent scalar behavior, not a source for our proposed two-field model.
2. van de Bruck, Morrice & Vu, *Constraints on Disformal Couplings from the Properties of the Cosmic Microwave Background Radiation*, https://arxiv.org/abs/1303.1773. Relevant caution about radiation-sector coupling and observed distributions.
3. Prior calculations in `time_first_principles/goals_and_results.md` and its accompanying reproducibility package supply the previously reported spectral-aging and FIRAS comparisons. No observations were reprocessed in this step.

Run `python test_candidate.py` with NumPy and SciPy. `results.json` contains all new numerical results. The normal-mode calculation uses 20,001 sampled times over 0≤γt≤1. The field model and analytic redshift derivations in this note are new speculative constructions, not conclusions of the cited literature.
