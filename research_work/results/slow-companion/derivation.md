# A joint energy-angle interaction with a slower scalar companion

## Why test this candidate?

The candidate review called for an interaction that links energy transfer and direction in one model. A photon can emit a slower wave while conserving energy and momentum in a preferred rest frame. Unlike the earlier point-target filters, the emitted wave's dispersion relation fixes its energy and momentum together. We test whether this can produce small energy transfers without scattering photons out of sharp images or shifting different colors at different rates.

This is an optional scalar-field comparison, not an ordinary spin-two graviton or an adopted cosmic medium. Subluminal scalar emission is a known type of Cherenkov process; [Dalang, Fleury and Lombriser](https://arxiv.org/abs/2109.10812) study a different gravitational coupling and emphasize the importance of the effective-theory cutoff. Their rates are not imported here. The following rate is derived from the electric-response action specified below, without expansion, an initial dark-matter population or a Big-Bang history.

## 1. Explicit action, state and validity range

Use units hbar=c_photon=1, in a specified preferred frame:

    L = (1+g phi) E_field^2/2 - B^2/2
        + (partial_t phi)^2/2 - v^2 |grad phi|^2/2,
    0 < v < 1.

Fields are Maxwell A_mu and one real scalar phi. The interaction is (g/2) phi E_field^2. It preserves electromagnetic gauge invariance, time/space translation and spatial rotation in this frame. Boost invariance is not assumed. The scalar has one polarization and dispersion omega=v|Q| about phi=0. Its massless choice is an explicit unprotected parameter choice; radiative corrections and a full material/gravity completion are missing.

The coupling g has inverse-energy dimension. Use a finite EFT cutoff Lambda with all participating energies and momenta well below Lambda, perturbative g and a decay rate small compared with photon energy. The numerical table reports rates divided by g^2 E^3, not physical rates at an adopted coupling. The perturbative domain must be rechecked for small v, where phase space enhances the event rate. No cutoff-crossing extrapolation is justified.

With D=(1+g phi) E_field and 1+g phi>0, the classical Hamiltonian density is

    H = D^2/[2(1+g phi)] + B^2/2
        + (partial_t phi)^2/2 + v^2 |grad phi|^2/2.

It is positive in this local domain. This is not a global nonlinear stability or quantum vacuum proof. A sustained radiation field sources phi; an exactly fixed phi=0 illuminated background has not been established. The tree-level calculation uses asymptotic fluctuations about that background, with initially empty companion modes.

The local energy equations are

    partial_t u_EM + div(E_field cross B) = -(g/2) partial_t phi E_field^2,
    partial_t u_phi + div(-v^2 partial_t phi grad phi) = +(g/2) partial_t phi E_field^2.

Thus the interaction exchanges energy with the scalar field, rather than losing it. The scalar equation is phi_tt-v^2 Laplacian(phi)=(g/2) E_field^2. Time-independent free propagation conserves the energy of a companion mode in this frame. Capture, permanent deposition and gravitational response do not follow from that fact.

## 2. Kinematics derives the energy-angle relation

For gamma(k) -> gamma(k') + phi(Q), let incident photon energy be E, surviving energy E'=xE, photon deflection cosine mu, and companion energy omega=(1-x)E. Translation symmetry gives Q=k-k'. The scalar dispersion yields

    (1-x)^2 = v^2(1+x^2-2x mu),
    mu(x) = [1+x^2-(1-x)^2/v^2]/(2x),
    (1-v)/(1+v) <= x <= 1.

Each event retains a photon and deposits its lost energy and momentum into the scalar mode. No infinitely heavy target absorbs unrecorded recoil. The emitted mode is subluminal in a preferred frame; its energy-momentum is not the null four-momentum of an ordinary graviton. Conservation here does not establish its stress tensor's gravitational coupling.

## 3. Matrix element and phase space

With canonically normalized modes, the vertex between transverse photons is

    M = g E E' (epsilon dot epsilon').

The general electric polarization factor is E epsilon_vector-k epsilon_0. Replacing a polarization by its four-momentum makes this factor vanish, providing the photon Ward check. Averaging over initial polarizations and summing over final polarizations gives

    mean |M|^2 = g^2 E^2 E'^2 (1+mu^2)/2.

Using the scalar measure d^3Q/[(2 pi)^3 2 omega] and the energy delta function, the angular Jacobian is

    |partial omega/partial mu| = v^2 E E'/omega.

The differential spontaneous rate per unit preferred-frame time is therefore

    dGamma/dx = g^2 E^3 x^2(1+mu(x)^2)/(32 pi v^2).

The outgoing particles are distinct, so no final-state identical-particle factor is inserted. This kernel includes all directions, not just photons admitted by a telescope aperture. Its polarization dependence must be retained in any full polarized transport calculation.

Define K(v)=integral x^2(1+mu^2)(1-x) dx/(32 pi v^2). The instantaneous fractional energy-loss rate for a photon of energy E is

    alpha_time(E) = K(v) g^2 E^3.

With photon speed one, it is also the fractional loss per path length between events. Doubling the photon energy increases this fractional loss rate by eight at fixed g and v. This exact scaling within the stated scale-free EFT is already unsuitable as an achromatic redshift law. Altering v or the overall coupling does not remove it. A frequency-dependent medium/cutoff or another operator would be a different calculation, not a successful fit of this kernel.

For a broad ensemble, d< E >/dt=-K g^2 < E^4 >, not -K g^2 < E >^4. The deterministic mean-field solution is therefore not an exact finite-distance prediction for the large-jump branch.

## 4. Small packets and sharp images pull in different directions

Use y=(1-x)(1+v)/(2v) on [0,1]. In the slow-wave limit, mu approaches 1-2y^2 while the unnormalized distribution approaches 1-2y^2+2y^4. This gives

    <1-x> -> (10/11) v,
    <1-mu> -> 46/77,
    <(1-x)^2>/<1-x> -> (46/35) v,
    alpha_time/(g^2 E^3) -> 1/(12 pi).

Small energy loss does not mean a small deflection: the slow wave can carry appreciable momentum with little energy. In the opposite limit v approaches 1 from below, the surviving-energy distribution approaches 3x^2, so mean loss approaches 1/4 even as directions become increasingly forward. The exact v=1 collinear endpoint is not used as a separate finite-rate derivation.

| Companion speed / photon speed | Mean fractional loss per event | RMS photon deflection per event |
| ---: | ---: | ---: |
| 0.000001 | 0.000000909 | 75.30 degrees |
| 0.01 | 0.008923 | 74.52 degrees |
| 0.1 | 0.07570 | 67.46 degrees |
| 0.5 | 0.20238 | 38.09 degrees |
| 0.99 | 0.24913 | 4.05 degrees |
| 0.999999 | 0.25000 | 0.04051 degrees |

These are single-event, all-direction moments. They are not observed image sizes, line widths, or a finite-distance diffusion simulation. An aperture-selected sample would have a different distribution and survival probability. The color dependence excludes this simple kernel as the requested common fractional drift even without claiming a universal numerical image-blur bound.

## 5. Inverse processes, timing and missing completion

The same Hermitian interaction permits gamma(k')+phi(Q)->gamma(k). For occupation numbers f_k, f_k' and n_Q, the forward-minus-reverse occupation factor for this channel is

    f_k (1+f_k') (1+n_Q) - f_k' n_Q (1+f_k).

It vanishes at a common Bose-Einstein temperature with zero chemical potentials and E=E'+omega; the script checks this identity. A populated companion bath can therefore stimulate emission and return energy to photons. The spontaneous kernel does not justify ignoring those effects throughout a long source history. No thermal companion bath is assumed here.

The preferred-frame stationary interaction also does not by itself stretch the separation between arbitrary emission events. Its delays, angular redistribution and detector response must be calculated; mean energy loss alone supplies no arrival-time dilation law. The phi-dependent electromagnetic response affects material standards, so atoms cannot be exempted without a common matter action.

## Decision and next work

This explicit joint kernel changes the earlier point-target setup and allows a nonzero surviving-photon conversion channel. It nevertheless fails the desired achromatic fractional drift in its specified domain. Keep the action, rate and tradeoff as a scoped negative result; do not adopt it or normalize it to galaxy redshifts.

Any extension must change the scale-free coupling/dispersion assumptions and derive the resulting medium response, fluctuations and energy destinations together. A microscopic model with a physically specified dispersion scale or symmetry is a meaningful next comparison; an arbitrary energy-dependent coupling chosen to cancel E^3 is not a derivation. Capture and gravity remain separate uncompleted requirements.

## Verification

The script checks random energy/momentum configurations, polarization sums, the electric Ward factor, the delta-function Jacobian, independent integration coordinates, limiting moments and bosonic detailed balance. These verify the analytic calculation, not observational success, nonlinear stability, radiative naturalness or a complete underlying medium.
