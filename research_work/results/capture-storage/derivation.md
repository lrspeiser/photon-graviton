# Catching companion energy, retaining it and running out of storage

## Purpose in plain language

A companion passing through a galaxy is not yet a deposit. We need a receiving system that can catch its energy and hold it. Otherwise, integrating a capture rate forever could count energy that has already escaped or assume an unlimited supply of empty storage spaces.

This pass tests an optional two-step receiver. An easily excited state catches energy; a second, longer-lived state stores part of it. A physical example might involve internal states of existing matter or a collective state, but **no such material has been identified for our companions**. We do not introduce an independent dark-matter population or claim that ordinary gravitons already have these couplings.

The model answers three questions: how absorption relates to escape, where energy goes during the transfer, and how storage capacity limits accumulation. It extends the earlier [recoil and inverse-decay kinematics](../microphysics/action-specification.md), which did not calculate capture bandwidth or storage dynamics.

## 1. Why simply making the absorbing state permanent is not enough

The [Particle Data Group resonance formula, section 51.1](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-cross-section-formulae.pdf) relates a resonant cross section to the entrance and exit branching fractions and the resonance width. We use its isolated narrow-resonance form as a conditional scattering model. The following bandwidth/lifetime comparisons and storage calculations are our derivations, not observational results from that source.

Let an incident companion excite a receiver from g to b. Write Gamma_c for the energy width associated with returning a companion, Gamma_s for a transition to a storage channel, and Gamma for all widths combined. With a nearly constant momentum/spin prefactor sigma_unit,

    sigma_inclusive(E) / sigma_unit
      = Gamma_c Gamma / [4((E-E0)^2+Gamma^2/4)]
    sigma_shelving(E) / sigma_unit
      = Gamma_c Gamma_s / [4((E-E0)^2+Gamma^2/4)].

Here sigma_unit=4 pi zeta/k^2 in hbar=c=1; zeta accounts for the actual spin/polarization multiplicities. It has not been assigned to a specific companion/material. Inclusive resonance formation followed by an exit is not permanent capture. The shelving cross section describes the specified storage exit, which must itself be long-lived.

Integrating the narrow line while treating its prefactor as constant gives

    integral sigma_inclusive dE = sigma_unit pi Gamma_c/2
    integral sigma_shelving dE = sigma_unit pi Gamma_c B_s/2
    B_s = Gamma_s/Gamma.

The usual extension of the narrow line to the whole real energy-offset axis is used; this is not a threshold model or a broad-resonance result. A weak entrance coupling narrows the integrated acceptance even if the exact resonance peak remains large. For two exit channels the shelving peak is sigma_unit B_c B_s, at most sigma_unit/4, while its integrated area at fixed Gamma_c approaches sigma_unit pi Gamma_c/2 as Gamma_s increases. A peak and a band-integrated capture rate are different quantities.

For a single directly absorbing state, its exponential lifetime is tau=hbar/Gamma. Requiring survival S over time T implies Gamma <= -hbar ln(S)/T and Gamma_c <= Gamma. For a flat incident spectrum over a band B much wider than Gamma and centered on the resonance,

    average sigma_inclusive / sigma_unit <= pi[-hbar ln(S)]/(2 T B).

For the illustrative values S=0.9, T=10 billion years and B=10^-6 eV, the ratio is at most 3.45e-28. The corresponding width is at most 2.20e-34 eV. These are a bandwidth comparison and a finite-time survival requirement, not an adopted cosmic age or an absolute capture cross section. A precisely matched narrow coherent input, another receiver structure, a distribution of resonances, or a different interaction needs a separate analysis. Lifetime broadening, line overlap and threshold effects cannot be ignored outside this approximation.

Consequently, declaring the directly absorbing state permanent while retaining a finite ordinary reverse transition is inconsistent in this model. A protected second state is one way to separate easy entry from slow escape; it is not yet a proof that the desired receiver exists.

## 2. A reversible three-state Hamiltonian

Take states g,b,d with internal energies 0,E_b,E_d, where 0<E_d<E_b. Let c_k be companion modes, a_k modes that receive the shelving energy E_b-E_d, and r_k modes that receive eventual deposit leakage E_d. A rotating-wave comparison Hamiltonian is

    H0 = E_b |b><b| + E_d |d><d| + sum bath mode energies
    Hint = sum_k [g_k |b><g| c_k
                  + h_k |b><d| a_k
                  + l_k |d><g| r_k + Hermitian conjugates].

The mode couplings include the necessary dimensions and normalizations. Resonant transitions conserve the corresponding level-plus-mode energy. Hermiticity includes each reverse process. In the weak-coupling continuum/Markov limit, bosonic matrix elements sqrt(n) and sqrt(n+1) give rates

    g -> b: kappa_c n_c           b -> g: kappa_c(n_c+1)
    d -> b: kappa_s n_s           b -> d: kappa_s(n_s+1)
    g -> d: kappa_d n_r           d -> g: kappa_d(n_r+1).

Each kappa is a golden-rule expression proportional to the squared corresponding matrix element times the mode density, evaluated at its gap. Gamma=hbar kappa for the associated lifetime width. The actual matrix elements, mode densities and occupations have not been derived for a material. This is a reversible effective candidate, not a fundamental graviton interaction or a calculation of their numerical values.

The rate approximation assumes incoherent reservoirs, resolved levels, weak coupling and negligible memory. Dimensionless rate ratios used below can all be rescaled to rates much slower than E_b/hbar; rate unity does not mean a width equal to the energy gap. Finite recoil and target motion are omitted in this heavy-receiver internal-state calculation. Their energy and momentum must be restored in a physical model; an anchor or prescribed bath is not a free source or sink.

## 3. Energy accounting and an explicit example

Let p_g+p_b+p_d=1. With empty shelving/leak reservoirs (n_s=n_r=0), W=kappa_c n_c gives

    dp_g/dt = -W p_g + (W+kappa_c)p_b + kappa_d p_d
    dp_b/dt =  W p_g - (W+kappa_c+kappa_s)p_b
    dp_d/dt =  kappa_s p_b - kappa_d p_d.

The stored internal energy is U=E_b p_b+E_d p_d. Its balance is

    dU/dt = E_b W p_g                          incoming companions
            - E_b(W+kappa_c)p_b               returned companions
            - (E_b-E_d)kappa_s p_b            shelving radiation/receiver energy
            - E_d kappa_d p_d.                deposit leakage

The implementation includes occupations in all three reservoirs and records incoming and outgoing energy separately for each. At a common thermal temperature it recovers a stationary Gibbs distribution, checking that reverse transitions have not been removed by hand. The reservoirs are prescribed boundary conditions; recorded energy flow is not a solved closed-universe source history.

For an initially caught packet in b, choose E_b=1, E_d=0.9, kappa_s=9 kappa_c and kappa_d=0, with empty reservoirs. The competing exits give:

| Destination | Fraction of initially caught energy |
|---|---:|
| Long-lived internal energy | 0.81 |
| Reemitted companions | 0.10 |
| Shelving radiation or receiving system | 0.09 |
| Total | 1.00 |

These numbers follow from chosen parameters; they are not measured efficiencies. In particular, 90% probability of reaching the storage state is only 81% retained energy in this example. The 9% released during shelving must have a destination and a spectrum. Calling it heat or radiation does not dispose of its observational consequences. Choosing a smaller level separation can retain more energy per successful event, but it also changes the transition rate and thermal occupation; our free rate parameters do not prove that efficiency and rapid capture can be independently selected in a real system.

## 4. Permanent storage eventually fills

With kappa_d=0, n_s=0 and sustained nonzero companion illumination, every receiver in this three-state model eventually ends in d. Then p_d -> 1 and its further absorption stops. The maximum long-lived energy is E_d per receiver, or

    u_d,max = n_receiver E_d

per physical volume. Permanent storage does not mean infinite storage. More levels, field amplitudes, new receivers or a phase transition may change capacity, but must be explicitly specified and energetically supported.

For finite kappa_d, the steady state under constant W satisfies

    p_b = 1 / [1 + (W+kappa_c+kappa_s)/W + kappa_s/kappa_d]
    p_g = (W+kappa_c+kappa_s)p_b/W
    p_d = kappa_s p_b/kappa_d.

The numerical evolution agrees with these limits and with independent matrix exponentials. Even before full saturation, a fixed capture coefficient multiplied only by companion energy density misses the declining number of available receivers. A completed reservoir law needs the bright and occupied populations, not only u_companion and u_deposit.

## 5. A warm receiving bath can reopen the storage state

The Hamiltonian that allows b -> d plus a shelving quantum also allows d plus that quantum -> b. If n_s>0, an otherwise protected d state can therefore be reexcited and release a companion through b -> g, even when the direct leakage kappa_d is zero.

After switching off companion input, with kappa_d=0, set

    B=kappa_s n_s; C=kappa_s(n_s+1); A=kappa_c+C.

The coupled b,d populations have a slow decay rate

    lambda_slow = 2 B kappa_c /
                  [A+B+sqrt((A+B)^2-4 B kappa_c)].

It is positive whenever both reverse excitation and companion escape are allowed. The code checks this against the exact generator eigenvalues and tracks energy drawn from the receiving bath during reexcitation. Thus setting the direct leakage to zero alone does not establish permanence. A physical calculation must constrain the bath spectrum, temperature, collisions and other reopening channels.

## 6. Capacity and location constrain the proposed extra gravity

If deposits are internal energy in existing matter, that energy increases the receivers' total mass-energy; it must not also be counted as a second independent gravitating substance. Under the cold, ordinary-gravity approximation, the added rest-mass density is n_receiver E_d p_d/c^2. Pressure, recoil, binding stresses and lensing still need the complete stress tensor.

An illustrative single 1 eV slot per proton would add at most about 1.07e-9 of the original proton rest mass. Adding as much rest mass as the original protons would require about 9.38e8 eV of stored energy per proton. This is a capacity benchmark, not an observed dark-matter mass or a claim that every possible store has an eV capacity.

Moreover, with the same illumination and rates everywhere, p_d is the same everywhere and the deposits follow the receiver number density. Absorption by stars alone does not automatically create an extended galactic reservoir. The spatial locations and support of the actual receiving system must be included before claiming the desired rotation or lensing profile. Modified gravitational response remains a separate candidate requiring its own derivation; it is not inserted into this capacity comparison.

## What this pass establishes and what remains

We now have a reversible capture/storage candidate with a testable energy ledger, finite bandwidth, explicit leakage channels and finite capacity. It gives concrete requirements for replacing the earlier phenomenological deposition rate. It does **not** identify a material, calculate an ordinary-graviton capture cross section, derive a new population, or complete the photon-to-extra-gravity chain.

The next unified capture model must specify receiver identity and abundance, its energy spectrum and transition matrix elements, available empty capacity, its spatial distribution, all release/reverse channels, and the motion/lensing response. Recovered astronomical observations and all 32 original response areas remain in scope. These conditional mathematical checks are not astronomical validation.

Run `python -X utf8 research_work/results/capture-storage/check_capture_storage.py` or the full diagnostic runner. Saved numerical evidence is in `capture-storage-results.json`; fresh runs default to an ignored generated directory.
