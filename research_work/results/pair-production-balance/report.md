# A finite illuminated population, with creation and loss from one interaction

The optional scalar-pair model now has a common rate calculation for creation, inverse photon scattering and two-photon annihilation. Under the stated steady illumination it approaches a finite average population, rather than accumulating energy indefinitely. This result concerns one specified pair of bound modes, not the storage capacity of an entire galaxy.

## What changes physically

The same interaction lets a photon lose energy to create a pair or gain energy by removing one. In the point-electric-dipole approximation, the inverse cross section at the same incident photon energy is larger. Bosonic occupation factors then produce a stationary distribution. Spontaneous two-photon loss provides an additional drain, whose rate scales as the seventh power of the pair energy gap **when the transition polarizability is held fixed**.

That qualification matters: a tiny gap can give slow spontaneous loss, but the real bound-mode overlap and coupling are not known. No astronomical lifetime has been established, and no independent production and decay efficiencies have been inserted.

## The formulas

In units hbar=c=epsilon0=1, for photon energy E, pair gap delta and transition polarizability Lambda:

```
sigma_create = Lambda^2 E (E-delta)^3/(6pi), for E>delta.
sigma_inverse = Lambda^2 E (E+delta)^3/(6pi).
C = single-pair spontaneous decay rate = Lambda^2 delta^7/(1680 pi^3).
```

For photon flux F, let A=F sigma_create and B=F sigma_inverse. The rates at pair occupation n are A(n+1)^2 for creation and (B+C)n^2 for removal. They give

```
mean(n)_stationary = A/(B+C-A).
stored pair energy = delta mean(n).
```

These are leading heavy-store, point-response rates. The [derivation](derivation.md) gives field normalization, polarization factors, phase-space integrals, occupation assumptions, finite source accounting and the limits of the approximation.

## Illustrative results

The values below are stored energy divided by the incident energy of one photon. Each row describes a different stipulated gap; neither the coupling nor mode structure is fitted to observations.

| Pair gap / photon energy | Weak spontaneous loss, C/A=1e-6 | C/A=1 |
| ---: | ---: | ---: |
| 0.01 | 0.161708 | 0.009418 |
| 0.1 | 0.121096 | 0.054771 |
| 0.2 | 0.084210 | 0.059259 |
| 0.5 | 0.019231 | 0.018519 |

As the positive gap becomes very small and C/A becomes negligible, this model's equilibrium energy approaches E/6 per pair of modes. The number and overlap of actual modes remain unknown; multiplying by an invented number would not establish a halo energy supply.

At equilibrium, any residual energy taken from the incident light pays for escaping two-photon radiation while stored energy stays fixed. If spontaneous loss is negligible, creation and inverse scattering nearly cancel in the energy account. After illumination stops, a single pair has survival probability exp(-C t); a populated mode has the larger n^2 removal factor. The calculations include that dark evolution rather than imposing a permanent deposit.

## Numerical evidence

Twelve rate combinations were tested at two occupation cutoffs and integration tolerances: 24 finite illumination histories followed by 24 dark histories. The largest change in predicted mean occupation under refinement is 1.11e-7, below the 1e-6 gate. Final mean populations agree with the analytic stationary result within 2.48e-11. The maximum probability-normalization error is 1.04e-12 and the leading-model event-energy residual is 8.57e-13. Integrated omitted creation at the numerical boundary is at most 6.72e-8 events, below the 1e-6 cutoff gate.

Independent angular and frequency integrals agree with the rate coefficients within 2.23e-16 relative error. A separate single-pair decay control agrees with the exponential within 5.79e-11. An isotropic finite-seed fixture bounds omitted expected recoil heating, including all subsequent dark decays, below 2.24e-10 in photon-energy units. That is a controlled approximation check, not exact finite-mass momentum evolution.

The finite illumination interval has a declared mean input-energy budget; probabilities, event counts and energy accounts are ensemble averages, not a closed finite-number photon-source simulation; the initial companion seed is accounted separately and its origin remains unexplained. [Protocol](protocol.json), [code](check.py) and [saved results](pair-production-balance-results.json) make the diagnostic reproducible.

## Consequence for the time concept

A steady mean occupation supplies a steady mean density in this incoherent model. Its mere existence therefore does not maintain the changing propagation factor needed by our current time-stretch formula. Previously acquired redshift is not undone. Coherent evolution, spatially distinct production and storage, time-dependent sources and other interactions remain separate possibilities.

Next derive spatial modes and the coherent forward response from the same parent interaction, including energy exchange and source/detector clocks. The direct-event color problem, seed supply, weak-gravity coupling, mode count, supported halo structure and joint rotation/lensing remain unresolved. All 20 tasks and 32 observational areas remain in scope; these tests do not establish a completed graviton theory.
