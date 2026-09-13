# Reciprocal scattering: flat forward redshift does not guarantee growing storage

13 September 2026. Conditional stationary-population result, not a galactic capacity or lifetime prediction.

## Result

Even an ideally achromatic forward transfer law has a reverse channel if the microscopic interaction is time-reversal invariant under the stated channel assumptions. In the inherited two-mode bosonic population model, steady monochromatic illumination gives finite occupation. With no spontaneous decay, forward and inverse energy transfers cancel at stationarity; adding a decay channel makes residual photon energy loss feed that channel rather than grow stored energy indefinitely.

The new result links the improved color target to the older storage problem. It does not assume the resonance construction is a completed microscopic interaction, and it does not exclude fresh modes, transport into protected states, evolving illumination, or a different mode statistics.

## Channel relation and provenance

**Known scattering consequence of time-reversal invariance**, equal channel degeneracies, consistently integrated reversed directions/polarizations, and the heavy-store limit: at the same total energy, photon momenta weight the forward/reverse cross sections. Writing forward photon energy E+Delta and reverse incoming photon energy E,

\[
(E+\Delta)^2\sigma_+(E+\Delta)=E^2\sigma_-(E).
\]

Here the store gains Delta in the forward transition and loses it in reverse. The assumption includes the participating reservoir state, not just the photon. Occupation factors are applied separately below. Finite store recoil, degeneracies, polarization restrictions and nonreciprocal backgrounds require their corresponding factors. Reciprocity has precise conditions and is not interchangeable with every intuitive notion of reversing a wave; see [Deak and Fulop](https://arxiv.org/abs/1108.5743). Standard scattering kinematics and phase-space factors are summarized in the [PDG kinematics review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf).

**Ideal forward target, not a new measured law:** let alpha_0=n_store*sigma_+(E)*Delta/E be constant. Thus sigma_+(E)=K*E over the energies used. Reciprocity then requires

\[
\sigma_-(E)=K\frac{(E+\Delta)^3}{E^2},\qquad
\frac{\sigma_-(E)}{\sigma_+(E)}=(1+\Delta/E)^3>1.
\]

This is not a fit of the inverse process. The same ideal forward rule fixes it at the related energy. The relation need only apply above the transition threshold where the flat target is stipulated; no unphysical extension below threshold is used.

For a forward coefficient that is not exactly flat, the ratio gains the factor alpha_+(E+Delta)/alpha_+(E). Therefore the prior finite-band resonance residuals and actual channel response would need to be inserted before claiming a quantitative result for that full ensemble. Here the exactly flat case is an explicit limiting test.

## Populated pair of modes

Reuse the conditional bosonic pair ladder from [production/loss balance](../pair-production-balance/derivation.md): an occupied state n has forward rate A(n+1)^2, reverse rate B*n^2 and independent spontaneous removal C*n^2. A and B use the same monochromatic photon flux, so B/A=(1+Delta/E)^3. C>=0 remains a diagnostic extra channel, not a calculated damping width. The C=0 case already establishes finite stationary occupation in this model.

**Known birth-death detailed balance:** define R=A/(B+C)<1. Then

\[
P_n=(1-R)R^n,\quad
\langle n\rangle=\frac{A}{B+C-A},\quad
U_{\rm pair}=\Delta\langle n\rangle.
\]

The stationary energy identity is

\[
\Delta\{A\langle(n+1)^2\rangle-B\langle n^2\rangle\}
=\Delta C\langle n^2\rangle.
\]

The left side is net energy removed from the incident radiation by forward and inverse scattering; the right side is energy leaving through spontaneous pair removal. Stored energy has zero time derivative at this stationary distribution. For C=0 both sides vanish, although individual scattering events still occur and can affect spectra.

## The one-third limit, correctly identified

Writing x=Delta/E and taking C=0 gives

\[
\frac{U_{\rm pair}}E=\frac{x}{(1+x)^3-1}\longrightarrow\frac13
\quad\text{as }x\to0^+.
\]

This is an algebraic limit of the stated reciprocal rates and bosonic occupancy model. It is **not** a derivation of our retention-law exponent, a fraction of all stellar energy captured, a statement about one photon per store, or a galaxy-wide energy budget. Illumination is a maintained bath, and the energy ratio is per specified pair of modes relative to the bath photon energy. Mode count, spatial overlap, initial seed energy and source duration remain unknown. The zero-gap endpoint has divergent occupation and is not a finite physical particle population.

The earlier point-response calculation approached E/6 per pair instead, because its inverse/forward ratio had a different small-gap slope. Changing the response changes this numerical limit without deriving the astrophysical retention function.

| Gap / bath photon energy | Inverse / forward coefficient | Stationary stored energy / bath photon energy, C=0 |
|---|---:|---:|
| 0.001 | 1.003003 | 0.333000 |
| 0.01 | 1.030301 | 0.330022 |
| 0.1 | 1.331000 | 0.302115 |
| 0.5 | 3.375000 | 0.210526 |

The larger-gap rows diagnose the ideal flat-rate hypothesis; they are not assertions that the small-gap Gaussian/resonance asymptotics remain valid there.

## Checks and scope

Sixteen cases combine four gaps with C/A=0,1e-6,0.1,1. Explicit sums of the geometric distribution verify normalization, occupation and stationary energy-flow identities against the analytic formulas. The occupation cutoff is chosen to leave a probability tail below 1e-15; no renormalization hides that tail. Relative energy-flow residuals are below 1e-10. These are equilibrium calculations, not time-dependent formation simulations or new astronomical fits.

This tests one pair of modes with inherited occupation factors, not the mode distribution that would supply a galaxy halo. Long storage in darkness cannot be inferred by setting C=0 as a diagnostic; an actual spontaneous rate must come from the common interaction. Likewise, breaking time-reversal symmetry or imposing a driven environment changes the assumptions but needs an explicit source and energy account.

## Decision

Preserve the resonance ensemble as an incomplete color-response candidate, and attach this storage restriction to it. A receiver cannot be counted as permanently accumulating solely because its empty-state forward coefficient has the desired color dependence. The next physical distinction is whether energy migrates into other, weakly illuminated states before inverse transfer, with that migration and its reverse channel derived from the same model. The previous spatial settling calculations provide conditional endpoints, not this microscopic protection mechanism.

The theory still needs source supply, mode creation, the time/event-duration mechanism, stable gravitational deposits and joint lensing/motion predictions. The exact-third reference remains an empirical hypothesis.

Reproduce with `python research_work/results/companion-extensions/reciprocal-storage.py`. [Source](reciprocal-storage.py), [complete results](reciprocal-storage-results.json).
