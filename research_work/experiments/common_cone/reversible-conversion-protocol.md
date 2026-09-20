# RC-1: reversible conversion, environmental response and a closed budget

Declared 20 September 2026 in response to the user's photon/radiation-to-graviton-and-back proposal, before these calculations. Continue SR-2 unchanged. This is a bounded consistency screen and derivation, not a galaxy/cluster fit or an ordinary-graviton identification.

## What conversion can and cannot change

In SR-1 with beta=0, the photon Hamiltonian is exp((1+s)U)|p|. Its derivative with respect to U is (1+s) times its energy. RB-1 derives the same leading source coefficient (1+s) for weak high-frequency massless field waves. Consequently merely exchanging equal energy between these two sectors leaves the leading scalar gravitational source unchanged. Repeated cycles do not multiply simultaneous stored energy. This does not exclude additional interaction-energy or nonlinear constitutive contributions; those must be derived, not inserted into a lens profile.

For diagnostic active-source coefficients chi_gamma and chi_g and converted energy fraction f, M_active=[(1-f)chi_gamma+f chi_g] E_stored/c^2, excluding additional interaction contributions. At chi_gamma=chi_g, conversion has no active-mass gain. These coefficients are diagnostic derivatives of an energy functional, not new free fit parameters in SR-1. Given a declared matter-energy allowance epsilon M c^2, any equivalent-energy explanation must satisfy chi_eff epsilon M >= M_required. For RB-1's finite-front example M_required=v^2(R-r0)/G. Evaluate epsilon={1,1e-3,1e-6}, f={0,.1,.5,1} and s={0,1}, with M=2e41 kg, v=200 km/s, R=6e20 m, r0=R/100. Show required chi_eff; do not count a chosen large chi as observational success.

## Explicit reciprocal local conversion Hamiltonian

Use two complex mode amplitudes a,b, environmental coordinate Q with momentum P, and dimensionless action units:

H = P^2/2 + Omega^2 Q^2/2 + (1+delta(Q)) |a|^2 + (1-delta(Q)) |b|^2 + 2 kappa(Q) Re(a* b),

delta(Q)=delta0 tanh Q, kappa(Q)=kappa0[1+A tanh Q], Omega=.2.

i adot=(1+delta)a+kappa b; i bdot=kappa a+(1-delta)b;

Qdot=P; Pdot=-Omega^2 Q-delta'(Q)(|a|^2-|b|^2)-2 kappa'(Q) Re(a* b).

Both conversion directions are mandatory. The real symmetric mode matrix gives exact norm N=|a|^2+|b|^2 conservation, and including both derivative forces gives autonomous total H conservation. Positive mode energy requires sqrt(delta^2+kappa^2)<1; the entire declared parameter box satisfies the stronger global bound sqrt(delta0^2+[kappa0(1+|A|)]^2)<1. This proves bounded local quadratic mode energy, not full spatial causality or gravitational stability. Q is a toy environmental degree of freedom, not a completed spin-carrying matter source. Initial radiation is prepared and counted; matter emission has not been derived.

If delta0=0 and the state starts purely in a, Re(a* b) stays zero. Changing conversion rate with Q alone then gives no conversion force on Q in this symmetric setup. At constant kappa and delta=0, |b|^2=sin^2(kappa t) for initial a=1,b=0: conversion can reach unity and reverse, without increasing total energy. Detuning/environmental force can change trajectories, but source recoil and interaction energy must remain counted. No desired orbit, lens map or target observation enters the local law.

## Preregistered numerical tests

- 24 reciprocal runs: kappa0={.05,.2,.4}, A={-.5,.5}, delta0={0,.3}, initial P={0,.2}; initial Q=-1,a=1,b=0,T=100. Use DOP853 rtol=1e-10,atol=1e-12,1001 output times. Require relative total-H and absolute N drift <1e-8; backward integration state error <1e-7.
- Four no-conversion controls: kappa0=0,A=0,delta0={0,.3},initial P={0,.2}. Receiving population must stay below1e-12.
- Three constant-conversion controls: kappa0={.05,.2,.4},A=delta0=0,P=0. Compare population to sin^2(kappa0 t), absolute error <1e-8.
- For the twelve delta0=.3 fixtures, deliberately omit the reciprocal environmental force. The full energy ledger must expose an error >1e-5. This wrong model is not a candidate.
- Verify source forces by 100 centered finite differences in Q (seed20260922), absolute error <2e-7. Check the analytic global eigenvalue bound.
- Record transferred energy, mode energies, signed interaction energy, peak/final receiving fraction and changes versus matched no-conversion Q trajectories. Test the symmetric delta0=0 no-force identity separately. Preserve first failures. Independently reconstruct archived H/N and analytic controls.

## Attribution, identification and next decision

Hermitian two-mode mixing, coherent conversion and Hamiltonian recoil are established mathematics. Photon/graviton conversion in an external magnetic field has established Gertsenshtein precedent: Palessandro and Rothman, https://arxiv.org/abs/2301.02072 . This effective Q-dependent oscillator model neither implements that mechanism nor inherits its conversion probability, spin or polarization properties. Existing project CWC-1 already tested reciprocal radiation-to-scalar energy exchange; RC-1 specifically examines reversible two-mode conversion and the no-gain equal-response limit, not a rediscovery of CWC-1.

Strongly unequal gravitational response cannot simply be assigned to ordinary massless spin-2 gravitons while keeping all standard assumptions: Weinberg's universal-coupling result has stated Lorentz-invariance and S-matrix assumptions, https://doi.org/10.1103/PhysRev.135.B1049 . A fictional extension must specify its different fields/assumptions and test them. We claim no historical uniqueness.

Conserved total energy includes matter, radiation, receiving field and interaction terms; baryonic rest mass alone need not be conserved when energy is emitted. Compatibility with conservation does not prove enough force, stable orbits, transparent images, correct clocks or a joint observational fit. The next full model would need spatial propagation, spin/angular momentum, emission/capture/release, residence time, common local light/field cone and one universal law for galaxies and clusters. Conversion rates cannot be adjusted separately to target each orbit.
