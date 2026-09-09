# A possible starting point without special time stretching

This is a conditional three-mode toy, not a derived photon-graviton theory. Its purpose is to show exactly what a candidate interaction must add to the fitted energy-loss law. No spin, gauge-invariant gravitational vertex, momentum-allowed rate or deposited state is established by it. Ordinary dynamical evolution in time is retained; no altered flow of time or void clock gradient is assumed.

## A time-independent inelastic interaction

Let a describe a photon mode of frequency omega_0, b a lower-frequency photon mode omega_1, and c a companion mode Omega. Stipulate a bosonic Hamiltonian with real coupling g:

```
H_free = hbar (omega_0 a†a + omega_1 b†b + Omega c†c)
H_int  = hbar g (a b† c† + a† b c)
omega_0 = omega_1 + Omega        [resonance condition]
```

The first interaction term destroys a high-frequency photon and creates a lower-frequency photon and a companion. The Hermitian-conjugate term permits the reverse process. Photon number a†a+b†b is conserved in this toy. At resonance H_free commutes with H_int: one event conserves the sum of the mode energies. A time-independent total H conserves total energy, including interaction energy, also away from resonance; away from resonance one may not claim the free energies alone are separately conserved.

Starting with one high-frequency photon and no other quanta, the dynamics stays within two states:

```
|A> = |1,0,0>
|B> = |0,1,1>
H / hbar = [[omega_0, g], [g, omega_1+Omega]]
P_B(t) = sin²(g t)              [on resonance]
mean photon energy = hbar omega_0 - P_B hbar Omega
mean companion energy = P_B hbar Omega
```

The numerical runner verifies these probabilities and total energy by matrix exponentiation for 12 scenarios. This gives a concrete interaction algebra that exchanges photon frequency and companion energy without changing clock standards. It does **not** show that photons and real gravitons possess this coupling in empty space. Momentum, polarization, gauge constraints and the source of g have not been supplied. An environment that supplies recoil or phase matching must appear in the ledger and dynamics.

The finite two-state calculation also does not establish that a full interacting field theory has energy bounded below, a stable vacuum or causal propagation. Those consistency tests must be made for a completed action rather than inferred from this restricted mode sector.

## What prevents this toy from being the complete explanation

1. Conversion is reversible and oscillatory. It is not exponential monotonic energy loss with travel distance.
2. It gives a mixture of two photon frequencies, not a single smoothly shifted sharp line.
3. A continuum of companion modes, escape, decoherence or other open-system approximation might produce an effectively one-way rate, but must be derived from a larger conserved system. The companion sector cannot be discarded from the total bookkeeping.
4. The value and energy/environment dependence of g are stipulated. Angular-momentum and gauge constraints can suppress or forbid a proposed physical vertex even when this matrix model works.
5. There is no capture or stored state, and no gravitational field equation. Calling c a graviton does not supply those properties.
6. Stationary conversion/delay kernels do not automatically stretch complete astronomical events. The timing derivation in report.md remains relevant.

## Connecting a derived event rate to the measured coefficient

Suppose a future calculation supplies W(E,epsilon), the differential event probability per path length and fractional energy transfer epsilon. For a forward channel that preserves one photon:

```
alpha(E) = integral_0^1 epsilon W(E,epsilon) d epsilon
drift(E) = -E alpha(E)
energy_diffusion(E) = E² integral_0^1 epsilon² W(E,epsilon) d epsilon
```

These are the first and second jump moments, before adding inverse processes and other channels. The drift is the conditional mean change per unit path; a broad ensemble obeys d< E >/ds=-<E alpha(E)> rather than necessarily -alpha(<E>)<E>. A differential diffusion approximation requires sufficiently small jumps and controlled higher moments.

For independent identical fractional transfers, W reduces to an event rate Gamma per length at epsilon_0, giving alpha=Gamma epsilon_0. The mean energy then follows the fitted exponential law exactly, while the existing Poisson formulas predict added spectral width. Many different pairs (Gamma,epsilon_0) give the same redshift, so line widths and other observations are needed to separate them. Reverse processes depend on mode occupations and couplings and cannot simply be set to zero without a physical limit.

The fitted alpha from this pass is therefore a target for a first-principles rate calculation. A candidate should predict W, inverse rates, outgoing directions, polarization, arrival times and companion propagation together before its capture and gravity predictions are trusted.
