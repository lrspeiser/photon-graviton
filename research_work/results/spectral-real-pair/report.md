# First real SN 2006mk spectral-pair estimate

**Exploratory outcome:** the declared primary matcher estimates phases -3.4 and +16.7 days, a 20.1-day difference over 32.073696 observer days. The six continuum/host variants give 19.8-20.3 days. This is closer to the conditional 21.738983-day source interval under an event stretch of 1+z than to 32.073696 days under unchanged arrival intervals. It is not a calibrated significance test or a photon-companion validation.

## What was fitted

The two actual ESO flux arrays are mapped to the approximate rest wavelength grid using measured host z=0.4754, retaining observer dates unchanged. Weighted polynomial removal and cosine correlation against processed SNID library features select the best matching phase from each reference object; the median of the five best distinct objects is the estimator. SN 2006mk is absent from the reference identities. All matches and correlations are retained in results.json. The independently written matcher is not a replication of the SNID executable.

These data and the published aging result were already exposed during development. This is not blind confirmation or an independent replication of the entire published sample. Nearby template labels retain historical calibration assumptions; the library contains later observations than the original publication.

## Preprocessing sensitivity

| Polynomial degree | Include nuisance host basis | Estimated age span (days) |
| --- | --- | ---: |
| 2 | False | 20.1 |
| 2 | True | 19.8 |
| 3 | False | 20.3 |
| 3 | True | 20.0 |
| 5 | False | 20.2 |
| 5 | True | 20.2 |

The narrow range across these six choices is not a total uncertainty interval. The host basis is itself noisy, and shared extraction errors, wavelength/frame details, telluric residuals, intrinsic spectral diversity and template resolution remain incompletely modeled.

## Recovery at measured residual power

For each synthetic template and observation, set amplitude so its projected whitened signal power equals the measured quadratic-continuum residual power minus the expected diagonal-noise power. This resolves the previous mismatch between injected pre-projection RMS and observed post-projection RMS. It assumes that the positive excess structure is signal appropriate to the injected template; artifacts may violate that assumption.

The scaling uses the known relation A=sqrt(S_observed / ||P C^(-1/2) t||^2), where S_observed=||P C^(-1/2) f||^2-(n-k). This is a conditional statistical calibration rule, not a new physical equation. No observing dates or expected timing span set A.

Fourteen nearby pairs with labeled spans at least 20 days are reused over ten seeds. Entire target objects are excluded from their matches. Independent-pixel and illustrative AR(1) noise with adjacent correlation 0.5 retain the supplied marginal errors. The diagonal matcher is deliberately unchanged under correlated draws to expose sensitivity; rho=0.5 is not a measurement of the real covariance.

| Noise correlation | Mean span error (days) | RMS error (days) | Empirical error 5th/50th/95th percentiles (days) |
| --- | ---: | ---: | --- |
| 0 | 0.399 | 5.304 | -6.005, -0.650, 7.430 |
| 0.5 | 0.188 | 5.804 | -6.840, -0.700, 12.010 |

The 140 outcomes in each row reuse only 14 objects; they are not 140 independent supernovae. These distributions are not confidence intervals for SN 2006mk. In particular, do not divide its discrepancy from unchanged timing by the synthetic RMS and report a sigma rejection. Recovery uses a selected phase-span range and does not exhaust object-dependent bias.

## Consequence for the physical model

This restricted real-data calculation gives no reason to remove event stretching from the requirements. A stationary energy-loss rule with arrival time t_arrival=t_emission+T for a fixed path has derivative dt_arrival/dt_emission=1: losing photon energy alone does not produce the estimated slower spectral aging. A proposed companion/time mechanism must supply a causal arrival mapping as well as the wavelength change. This is a known kinematic identity, not a new theory.

A duration factor of 1+z is a comparison hypothesis here, not an input to the spectral phase estimate. The measured redshift only aligns wavelengths. The calculation neither establishes expansion nor identifies companion conversion as the cause. Joint brightness, full covariance and multiple independently calibrated source clocks remain required. All six research objectives remain open.

## Provenance and verification

FITS hashes are checked against the acquisition manifest; the library loader checks the archived template hash. Calibration excludes all same-object templates. A repeated full run must reproduce results.json byte for byte. Primary selection and all sensitivity settings are recorded in protocol.json, but were not registered against unseen observations.

[Input acquisition](../essence-spectrum-acquisition/report.md) | [Feature-strength audit](../spectral-feature-strength/report.md) | [Earlier conditional recovery](../spectral-observed-noise/report.md)
