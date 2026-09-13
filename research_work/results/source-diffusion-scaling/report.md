# Physical source broadening: a conditional diffusion calculation

A minimal radioactive-heating diffusion model can approximate, but does not exactly reproduce, a uniformly time-stretched reference curve while retaining fixed nuclear lifetimes. This is a synthetic source calculation, not a bolometric, spectral or photometric fit to SN 2006mk.

## Model and provenance

Use the established class of one-zone, radiation-dominated diffusion approximations with freely expanding supernova ejecta, constant effective opacity and fully trapped radioactive heating. Ejecta expansion here is the local explosion, not expansion of the universe. A source-model example is discussed in [the Arnett-model appendix of the SN 2019odp study](https://doi.org/10.1051/0004-6361/202346313).

The energy equation is dE/dt=Q-L-E/t, with diffusive luminosity L=2tE/t_d^2. E/t accounts for work done as the ejecta expand. Eliminating E gives dL/dt=(2t/t_d^2)(Q-L). These are known approximate source-physics relations; they are not new companion equations or an exact radiative-transfer solution.

Adopt illustrative heating Q(t)=0.8 exp(-t/8.8 days)+0.2 exp(-t/111.3 days). The mean-life benchmarks represent nickel/cobalt-like heating; the 0.2 coefficient is a declared simplification, not an inferred isotope mixture or precision decay-chain treatment. The reference diffusion time is 13.5 days, also illustrative. Full trapping, constant opacity, no extra power sources and the one-zone approximation restrict the result.

## Scaling result

The desired stretching factor is 1.528877, carried over from the conditional SN 2006mk template-width fit solely as a target ratio. The synthetic target is L_reference(t/1.528877), not the observed supernova light curve.

- Multiplying diffusion time by this factor gives 20.63984 days. The peak-normalized curves differ with RMS 2.93% over 5-120 days after explosion.
- Fitting diffusion time and one common amplitude while leaving the decay times fixed gives 21.15867 days and amplitude 1.35029. The residual RMS is 2.30% of the target peak over that illustrative window.
- Rescaling both diffusion time and heating lifetimes by 1.528877 reproduces the uniformly stretched curve to a maximum numerical error 3.2e-11. Rescaling nuclear lifetimes is a mathematical control, not an allowed change to measured constants.

Therefore fixed nuclear clocks break exact template-time scaling in this model, but the resulting shape difference can be small. Its size relative to peak luminosity is not an observational rejection threshold. Source spectra, band fractions, trapping changes and intrinsic diversity may change the comparison.

For the usual conditional scaling t_d^2 proportional to kappa*M/v, multiplying diffusion time by 1.528877 requires kappa*M/v to increase by 2.33747. This is not a measured mass requirement: opacity and velocity can also change, and fitted template width need not equal diffusion time. No limit on physically allowed broadening is inferred here.

## Conservation and numerical checks

For both the reference and diffusion-broadened examples, independently integrating emitted radiation, ejecta expansion work and final stored radiation recovers the injected heating to relative error below 8e-11 through 120 days. This verifies energy bookkeeping in the adopted ODE. It does not validate full ejecta hydrodynamics, a nuclear-energy budget in physical units, or the actual supernova.

## Consequence for the hypothesis test

We cannot use the previous fixed-template score gap to exclude every unchanged-propagation source model. A slower physical source remains a possibility requiring a quantitative fit. Conversely, the present toy calculation does not demonstrate that SN 2006mk had such a source. A full test needs a physically constrained source family, spectra/velocity information, calibrated multi-band or bolometric flux and realistic covariance.

A uniform source-time rescaling is not automatically a first-principles solution. The energy and diffusion equations show which physical quantities must change and which independent clocks resist exact rescaling. This complements the existing late radioactive-tail audit and does not repeat its optical-tail fitting assumptions.

All six research objectives remain open. No companion-transfer coefficient, gravity deposition law, universe age or universe size was changed.

[Observed width ambiguity](../sn2006mk-intrinsic-width/report.md) | [Existing late-tail clock analysis](../radioactive-clock-feasibility/report.md)
