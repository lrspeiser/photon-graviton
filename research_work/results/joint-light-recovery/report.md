# End-to-end joint redshift, timing and brightness recovery

The combined forward calculator and correlated-flux likelihood recover known parameters from two artificial controls. The data generator uses an independent analytic fixed-band expression, while the fitter uses numerical real-filter integration. This tests that the assembled pipeline can infer shared parameters; it is not evidence that the physical companion hypothesis fits astronomical data.

## Frozen control design

Four artificial sources sit at exactly known distances of 30, 300, 1000 and 3000 Mpc. Every source has the same emitted flat-frequency spectrum, absolute AB magnitude -19 and 30-day Gaussian FWHM. The spectrum has finite support from 1000 to 20000 angstrom, covering all relevant emitted passbands. There is no intrinsic source-population variation, dust or lensing.

Each source has observations in g/r/i/z at 16 observer times from -60 through 90 days: 256 flux values per control, plus four redshifts. Peak flux-to-error ratio is 30. Artificial zero-point errors are shared across all sources, drawn from the measured DES5YR covariance and applied in the same declared linear approximation as the likelihood. Redshift noise is Gaussian with standard deviation 0.0005. Distances have no uncertainty in this control.

Two paired cases use event exponents b=0 and b=1 with the same random draws. The conversion rate is 0.0002488993286382367 per Mpc. This value is borrowed from the earlier fit only as injected truth; the control does not validate its extrapolation to these large distances.

## Independent injection and common-parameter fitting

For the artificial flat-frequency spectrum, known radiometry gives a closed-form native-band peak flux proportional to D^-2 S^-b times the fixed calibration factor. The temporal width is 30 S^b days. We use that analytic expression to generate fluxes, without calling the forward integrator.

The fitting side instead evaluates the emitted spectrum through the actual filter curves and historical calibration conversion. At the injected parameters, those numerical predictions agree with the analytic source to relative error below 2.3e-15.

The fit has four common parameters: conversion rate, event exponent, luminosity normalization and intrinsic duration. There is no independently adjusted luminosity or duration for each source. The objective combines the normalized joint flux likelihood, including shared zero-point covariance, with a normalized Gaussian redshift likelihood. Three optimizer starts are retained.

These are explicit effective-model and statistical assumptions using known exponential-loss, radiometric and likelihood mathematics. They are not new first-principles derivations of conversion or event stretching.

## Results

| Injected b | Recovered b | Recovered alpha (Mpc^-1) | Common luminosity / injected | Common width (days) |
| ---: | ---: | ---: | ---: | ---: |
| 0 | -0.02101 | 0.00024891713 | 0.99738 | 30.2569 |
| 1 | 0.98208 | 0.00024891713 | 0.99891 | 30.2098 |

Both controls pass the pre-recorded implementation screening limits: rate error below 3%, exponent error below 0.1, and luminosity/width errors below 5%. All three starts converge to consistent interior fits. Those limits are engineering checks for this artificial exercise, not criteria for agreement with real observations or uncertainty coverage.

The negative fitted value near zero is allowed by the diagnostic optimizer bounds so the null control is not forced against a boundary. It is a noise-level deviation in the control, not evidence for physical time compression.

## What this establishes and what it does not

The code can predict and jointly infer color shift, duration scaling and detector brightness under one specified transport model with common source parameters. It no longer needs a separate fitted answer for every event in this controlled case.

This success relies on exactly known distances and identical sources with no evolution. Allowing arbitrary intrinsic duration and luminosity trends restores the previously derived ambiguity. Real supernova source populations, distances, calibration nonlinearities, passband changes, extinction, selection and peculiar motions remain to model or constrain. One paired seed cannot establish bias or interval coverage.

The ongoing 160-case timing calibration probes a different, broader set of source/noise realizations and already contains boundary failures. This successful integration control does not supersede those failures. No real flux fit has been performed; all six scientific demonstrations remain open.

```powershell
python research_work/results/joint-light-recovery/run.py
```
