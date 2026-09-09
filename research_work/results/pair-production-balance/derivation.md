# Creation and storage must use the same interaction

## Question and scope

The [bound-pair production calculation](../bound-pair-production/derivation.md) identified an optional charge-neutral scalar interaction that can create pairs in specified bound modes. It also contains inverse scattering and pair annihilation. This step derives their rates together and asks whether steady light keeps increasing the stored energy.

The calculation concerns one particle mode and one antiparticle mode in a seeded companion store, with a fixed positive pair gap delta. It is a weak-coupling, point-electric-dipole, heavy-store approximation with incoherent dilute incident light and freely escaping outgoing radiation. Bound mode profiles, the coupling's actual value, a galaxy's mode count, seed formation and ordinary-graviton identity are not supplied. The result must not be applied to a galaxy-sized target as though its spatial form factor were already known.

## 1. Specify the normalization before comparing rates

Set hbar=c=epsilon0=1. After projecting the parent operator onto the two bound modes, write its transition part as

```
H_int = (Lambda/2) (b^dagger d^dagger + b d) :E(X)^2:.
```

Lambda includes the underlying coupling and bound-mode overlap. It has units of energy^(-3) in this normalization. Holding it fixed while changing delta is a diagnostic comparison, not a derived constituent-mass scaling. Density terms in the parent operator are reserved for the forward/coherent response calculation; they cannot be counted as a second independent energy loss.

In a photon normalization box of volume V,

```
E(X) = i sum_(k,pol) sqrt(omega/(2V)) epsilon
       [a exp(i k.X) - a^dagger exp(-i k.X)].
```

For a zero-to-one-pair transition the relevant two-photon matrix-element magnitude is Lambda sqrt(omega1 omega2)/(2V), times the polarization overlap. The same magnitude occurs for one incoming/one outgoing photon and for two outgoing photons. Identical outgoing photons require a phase-space factor 1/2!, not a second independent coupling.

We use the transition rule 2pi times squared matrix element times final-state density; see [Tong's derivation of Fermi's golden rule, section 3.6.1](https://davidtong.org/pdfs/teaching/quantum-field-theory/qft3.pdf). The rates below are derived here for the proposed operator, not borrowed as an observed interaction.

## 2. Three rates from that operator

Let E be incident photon energy. At leading order in inverse receiving mass, pair creation leaves energy E-delta, while pair removal leaves E+delta. The cross sections for unit occupation matrix element are

```
sigma_create(E) = Lambda^2 E (E-delta)^3/(6pi), E>delta; zero below.
sigma_inverse(E) = Lambda^2 E (E+delta)^3/(6pi).
```

For fixed incident polarization, the final-polarization angular integral is 8pi/3. Combining it with the final photon density of states produces these coefficients. The inverse cross section is larger at the same incident energy because the higher-energy outgoing photon has more available phase space. This is a property of this point operator and fixed transition, not a universal preference for blue shifts.

The two-photon annihilation channel has total outgoing energy delta. Its polarization sum is 1+cos^2(theta12), whose double angular integral is 64pi^2/3. The remaining frequency integral is

```
integral_0^delta omega^3 (delta-omega)^3 d omega = delta^7/140.
C = Gamma_annihilate_one_pair = Lambda^2 delta^7/(1680 pi^3).
```

These phase-space factors include the heavy store as the momentum receiver. They are not the two-body decay phase space of an isolated freely moving pair with no store. Recoil corrections were kept exactly in the preceding kinematic calculation and are neglected at this leading order. The numerical protocol bounds their expected size for a finite illustrative seed.

At fixed Lambda, a small gap can make spontaneous loss very slow. It is therefore incorrect to infer a short lifetime merely from the existence of the annihilation channel. Conversely, an absolute astronomical lifetime cannot be quoted until delta, Lambda and actual mode overlaps are derived.

## 3. A population under illumination

Take isotropic monochromatic incident photons of energy E>delta and dilute number flux density F, with negligible final photon occupation. Define

```
A = F sigma_create(E),
B = F sigma_inverse(E),
C = Lambda^2 delta^7/(1680 pi^3).
```

For n particles and n antiparticles occupying the two modes, bosonic creation gives a squared matrix element (n+1)^2, while removal gives n^2. The pair-count transition rates are therefore

```
birth_n = A (n+1)^2,
death_n = (B+C) n^2.
```

The B process gives energy back to an incident photon; the C process emits two photons. They are distinct destinations within the total radiation account. The number difference of the two companion species remains zero. The general use of particle/antiparticle creation operators and their occupation factors follows [complex-scalar quantization](https://davidtong.org/pdfs/teaching/quantum-field-theory/qft2.pdf); the specific kinetic reduction is an assumption here.

This Markov population calculation is different from the preceding coherent two-frequency model. It assumes that outgoing photons escape and do not build up a stimulated reservoir. Incident photons all have energy above delta, so two of those photons cannot resonantly supply only the fixed gap in the heavy-store approximation. Incoming low-frequency photons whose energies sum to delta, stimulated final modes, optical thickness and coherent correlations would require additional terms. No such terms may be dropped when applying the model in an environment where they matter.

## 4. Steady light gives a finite population in this limit

The probability equation is

```
dp_n/dt = birth_(n-1) p_(n-1) + death_(n+1) p_(n+1)
          - (birth_n+death_n) p_n.
```

Detailed balance between adjacent states gives

```
r = A/(B+C) < 1,
p_n(stationary) = (1-r) r^n,
mean(n) = r/(1-r) = A/(B+C-A).
```

B>A for positive delta and this incident spectrum, so the geometric distribution is normalizable. Increasing the common coupling accelerates creation and loss together; it is not an independent knob that increases equilibrium occupation. Changing flux can change C/A, but at negligible spontaneous loss the limiting occupation depends on delta/E rather than the flux.

Write q=delta/E. When C/A tends to zero,

```
U_store/E = q / [((1+q)/(1-q))^3 - 1].
lim_(q->0) U_store/E = 1/6.
```

This is a conditional energy per specified pair of modes. It is not a bound on an entire galaxy or on a continuum with an unknown number of available modes. The zero-gap limit itself is singular: occupation diverges, so it cannot be adopted as a finite zero-energy particle store. The displayed finite-energy limit approaches it through positive gaps.

## 5. Energy flow at and away from equilibrium

Count creation events D, inverse-scattering events U and two-photon annihilations L. Each event changes pair energy by delta, so

```
delta [mean(n(t))-mean(n(0))] = delta [D(t)-U(t)-L(t)].
net Raman photon energy loss = delta [D-U].
two-photon outgoing energy = delta L.
```

The solver evolves all three event counters alongside the probabilities and checks this identity independently of normalization. Binding is already included in delta; it is not credited twice as constituent rest energy plus an additional release.

At stationarity, A mean((n+1)^2) = (B+C) mean(n^2). The residual Raman photon loss pays for the two-photon radiation C delta mean(n^2), while stored energy stays fixed. If C is negligible, forward and inverse Raman transfers nearly cancel; many individual scattering events do not imply net continuing storage. These events can still alter spectra and directions, which are not observationally tested by the population calculation.

The illumination lasts for a finite declared time. Its effective source area fixes a finite mean incident energy F A_source E T. An equivalent mean source-fuel account decreases by the mean incident power during that interval; a closed finite-photon-number source is not simulated. The area exceeds the summed scattering cross sections throughout the numerical occupation range. Mean outgoing energy equals that mean incident energy minus delta D plus delta U plus delta L. After switch-off no new incoming energy is supplied. Thus the maintained flux during illumination has an explicit finite mean input account. Probabilities, event counters and source energies in this kinetic calculation are ensemble averages, not individual photon histories.

The preexisting companion seed's rest energy is separate and unexplained by this calculation. The fixture uses isotropic independent incidence and a very heavy seed; the expected omitted recoil heating is bounded by

```
[(2E+delta)^2 (D+U) + delta^2 (L+mean(n_at_switch_off))]/(2 M_seed).
```

The final term includes all possible subsequent dark annihilations. Isotropy removes a mean directed impulse; this is a bound on expected kinetic heating, not a pathwise zero-recoil claim. The small bound checks this fixture's leading approximation. It does not replace a finite-mass, moving-store calculation for a galaxy. Energy accounting is exact in the leading kinetic model, with this separately reported recoil error for the finite illustrative seed.

## 6. Darkness does not mean an imposed permanent deposit

Once the incoming flux vanishes, A=B=0 and

```
d mean(n)/dt = -C mean(n^2).
```

For a single pair the survival probability is exp(-C t). For larger populations, the two-mode annihilation factor is n^2, so a simple exponential with a population-independent total loss rate is incorrect. The general bound mean(n(t)) <= mean(n(0)) exp(-C t) is tested. No-loss travel is a separate question; this is a bound-mode loss channel, not a propagation redshift assigned to free companions.

## 7. Numerical and theoretical limits

Twelve combinations of gap and C/A are run at maximum occupations 512 and 768 with tighter integration tolerances. Twenty-four finite illumination histories and twenty-four dark histories check probability, event-energy accounts, cutoff leakage and refinement. The omitted birth rate at the cutoff is integrated explicitly; it is not hidden by renormalizing probabilities. Analytic stationary distributions, a single-pair exponential control and independent angular/frequency quadratures provide further comparisons. Rates and conservation gates are specified in [protocol.json](protocol.json).

In this incoherent model a stationary pair population has a stationary mean density. If that density supplies an optical factor, its mere presence gives no continuing partial_t n. This connects to the earlier static-field limitation without undoing redshift acquired while the field was changing. It does not exclude a coherent evolving field, multiple spatial modes, separate production and protected storage regions, changing source histories, or a different interaction.

Next derive spatial receiving modes and the coherent forward response from the same parent operator, including the real clock/propagation law and its energy account. The direct-event color dependence remains unresolved. Mode count, seed supply, environmental coupling, long-lived halo support and simultaneous rotation/lensing predictions remain within the original 20 tasks and 32 requirements; none is completed by a stationary single-mode calculation.
