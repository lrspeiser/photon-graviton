# Resonant redshift: photon survival and energy destinations

13 September 2026. Conditional branching requirement, not measured extinction or a derived companion coupling.

## Result

The near-achromatic resonance response needs very selective channels if redshift comes from many tiny transfers. For the previous 2 eV optical example with 1e-8 eV per transfer and z=0.00766048, retaining 90% of photons requires the removal-to-useful-transfer rate ratio below 6.93e-8. This is roughly one removal per 14.4 million useful events. At z=1 the corresponding ratio is below 1.05e-9.

The 90% transmission condition is an illustrative allowance, not an astronomical measurement. These numbers constrain a possible implementation; the resonance model has not supplied its branching ratios. The calculation also does not identify every damping process with absorption: elastic radiation, lower-energy outgoing photons and dark loss must be separated.

## Channel definition and provenance

Use three destinations: elastic photon scattering, an outgoing photon with small energy loss Delta into companions, and removal of the photon from the tracked radiation population. The last can mean absorption or another excluded channel; its energy must still be recorded. Directional aperture losses require a separate geometric calculation.

**Known resonance-channel bookkeeping:** for an isolated resonance with common incoming channel, cross sections into final channels contain the corresponding branching fractions. Their ratio can therefore be written sigma_remove/sigma_inelastic=Gamma_remove/Gamma_inelastic when the same resonance denominator and channel normalization apply. See the [PDG resonance-formation formula](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-cross-section-formulae.pdf). This is not a derivation of a unitary completion for our broad oscillator ensemble.

Define r=kappa_remove/kappa_inelastic as the effective ratio after summing actual channels. **Hypothesis for this diagnostic:** r and the small inelastic gap Delta are constant along the path. A real ensemble generally requires energy-dependent partial widths; the constant-r case is a stated requirement to test, not a fitted microscopic law. The oscillator damping parameter cannot be identified with removal alone without specifying its outgoing channels.

## Drift equations and solution

**Known rate/energy accounting, imposed on the target fractional redshift:**

\[
\frac{dE}{ds}=-\alpha E,\qquad
\kappa_{\rm inelastic}=\frac{\alpha E}{\Delta},\qquad
\frac{dN}{ds}=-r\frac{\alpha E}{\Delta}N.
\]

Here N starts at one and denotes the surviving fraction, while E is the deterministic energy of a surviving photon in the small-gap drift approximation. Using E_f=E_0/(1+z), integration gives

\[
T=N_f=\exp\!\left[-\frac{rE_0}{\Delta}\frac{z}{1+z}\right],\qquad
r_{\max}=\frac{-\ln T_{\min}\;\Delta(1+z)}{E_0z}.
\]

These are analytic consequences of the stated equations, not a new fundamental law. They hold for nonuniform alpha as well when the integrated energy ratio and constant r assumptions hold. They do not include stochastic energy spread, survivor-selection effects, inverse transfer, stimulated channels, or source luminosity evolution. Small gaps justify testing the drift limit but do not remove the need for those effects in a full model.

## Explicit energy ledger

Per initial photon, define C as all energy transferred through the small-gap channel before either escape or removal. Define H as energy sent to the removed-photon sector. Then

\[
\frac{dC}{ds}=N\alpha E,\qquad
\frac{dH}{ds}=N r\frac{\alpha E^2}{\Delta},\qquad
E_0=N_fE_f+C+H.
\]

For r>0, C=(Delta/r)(1-T); its r->0 limit is E_0-E_f. H is the remaining energy. This includes useful transfers made by photons that are later removed. Simply counting the redshift of surviving photons would omit that contribution.

For 2 eV initial photons, gap 1e-8 eV, z=0.00766048 and T=0.9, the per-initial-photon energies are:

| Destination | Energy (eV) |
|---|---:|
| Surviving radiation | 1.786316 |
| Small-gap companion transfer | 0.014431 |
| Removed-photon sector | 0.199253 |
| Total | 2.000000 |

The removed sector receives much more energy than the intended redshift channel in this example. It cannot disappear from the account. It might feed another reservoir or be reradiated, but those possibilities need their own gravity, brightness and spectral predictions. Absorption energy is not evidence of additional gravity unless its retention and field response are specified.

## Numerical coverage

Three photon energies (0.01, 2, 100 eV), three gaps (1e-8, 1e-6, 1e-4 eV), and three redshifts (0.00766048, 0.1, 1) give 27 base cases. Each has four removal ratios: zero and 0.1, 1 or 10 times the 90%-transmission limit. All 108 drift ledgers close to numerical roundoff. Nine independent ODE integrations verify transmission and sector energies to better than 1e-9 normalized error.

| Redshift | Maximum removal/useful rate ratio, 2 eV and gap 1e-8 eV, T>=0.9 |
|---|---:|
| 0.00766048 | 6.93e-8 |
| 0.1 | 5.79e-9 |
| 1 | 1.05e-9 |

These are brightness-channel requirements only. Larger gaps relax the survival ratio but increase discrete-transfer spectral effects; none of the larger-gap cases is declared to pass the earlier linewidth or color tests. The z=1 calculation is a hypothetical integrated transfer target, not a newly fitted distant source or assumed expansion distance.

## Implication for the resonance candidate

The resonance population's nearly constant normalized rate remains a mathematical result. It is insufficient until the same couplings generate a suitably small photon-removal branch and acceptable elastic scattering, along with the desired inelastic rate. In a radiatively damped system, damping can populate photon-preserving channels, so positive damping alone does not prove excessive absorption. Conversely, declaring all damping harmless would omit its branching and angular effects.

The next microscopic step is a common-channel scattering model with widths fixed by couplings, rather than choosing damping, useful conversion and removal independently. Its optical response must be recalculated after enforcing those relations; the prior color cancellation is not guaranteed to survive. Matter clocks, event-duration stretching, photon supply, capture, permanent storage and the exact-third reference remain separate unresolved requirements.

Reproduce with `python research_work/results/companion-extensions/resonance-loss.py`. [Source](resonance-loss.py) and [complete results](resonance-loss-results.json).
