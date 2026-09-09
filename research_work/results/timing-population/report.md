# Integrated timing estimator: first numerical pilot

## What was built and why

The revised estimator integrates possible explosion shapes, peak times and durations instead of reducing every event to a single fitted width. Uncertain positive amplitude and background are analytically integrated using the previously checked likelihood component. Every one of the original 98 cadence slots is retained, including incomplete light curves. All fluxes are artificial; actual DES flux values are not read.

This is needed because a missing portion of an explosion should leave uncertainty about its duration. The first estimator's duration errors could create redshift-dependent bias, and cuts on fitted coverage could discard events selectively. A population likelihood can carry duration uncertainty forward, but its numerical integrals must be reliable before its results mean anything.

## Model and calculation

At each observer-time log-width node x, average the positive-amplitude flux likelihood over four independent uniform nuisance priors: peak between -10 and 10 days, rise fraction between 0.15 and 0.65, and separate rise/fall powers between 1.2 and 3.2. These are declared empirical assumptions, not derived explosion physics. Scrambled Sobol integration approximates that average.

For population parameters a, b and sigma, each event contributes

    L_i(a,b,sigma) = integral L_i(x) p(x | a+b log(1+z_i), sigma) dx.

The implementation uses normalized trapezoidal weights on log-width from log(2) to log(256). This is a finite-grid approximation to a truncated normal population. Probability outside those limits is separately reported; grid and prior-boundary sensitivity remain necessary. A constant event likelihood remains constant after integration, so an uninformative event does not by itself create a preferred width or slope.

Maximum marginal likelihood fits a, b and log(sigma) from three starts, with recorded numerical bounds. No event time is divided by an assumed (1+z), no fitted duration is used as an exact datum, and no calibrated confidence interval is claimed.

## Frozen pilot and findings

One newly seeded realization (601) uses the original split-Gaussian-plus-shoulder shape, injected b=1, intrinsic log-width scatter 0.1, and nominal SNR 20 at all 98 original cadence patterns. The shoulder is not an exact member of the estimator's unimodal shape family. This is development on already exposed cadence data, not independent validation.

| Integration | Nuisance samples | Width nodes | Recovered b | Recovered intrinsic scatter |
| --- | ---: | ---: | ---: | ---: |
| Base | 256 | 161 | 1.09507 | 0.10664 |
| Refined | 512 | 321 | 1.10019 | 0.10546 |
| Independent scramble | 512 | 321 | 1.10559 | 0.10966 |

The slope of the true injected widths in this particular finite sample is 1.04684. This distinguishes sampling variation from fitting error; it is not a truth-based correction available for real events.

All fits converged to interior solutions, and continuous population probability beyond the width boundaries is below 4e-35. Changes in b are about 0.0051 and 0.0054, within the frozen 0.05 numerical tolerance. **The numerical pilot nevertheless fails:** posterior-weighted changes in the centered individual-event log-likelihood curves are 0.177 and 0.218, above the frozen 0.1 tolerance. An apparently stable aggregate slope does not establish accurate event probabilities or trustworthy uncertainty intervals.

The base/refined comparison changes both grids; the independent scramble changes only the nuisance integration realization at the same resolution. The latter establishes that shape/peak integration error still matters. It does not isolate which nuisance dimension causes it. Refinement or adaptive integration must be recorded and checked before repeated-trial calibration; the threshold must not be relaxed after this failure.

## Verification and limits

The vectorized amplitude/background likelihood agrees with 12 scalar comparisons using the independently verified component. The constant-likelihood population check passes, invalid errors are rejected, every input cache hash is checked, and all 98 events remain present. The saved result includes optimizer starts, source hashes, convergence comparisons and hashes for each compressed event-likelihood array. Arrays contain synthetic-model results and IDs, not actual observed fluxes.

The code ran the three integrations in roughly 1.6, 5.9 and 5.9 seconds on this machine, excluding common input setup. This is a cost measurement for this pilot, not a general performance guarantee. No default physics job is added.

No claim of calibrated bias, interval coverage, observed time dilation or an established physical redshift mechanism follows. Expanded SNR, independent shapes, source evolution, filter throughput, selection and real-flux calibration remain required. Intrinsic duration evolution can still mimic propagation stretching.

```sh
python research_work/results/timing-population/run.py
```

Requires the verified timing-foundation cache, Python, NumPy, SciPy and Astropy. Fresh outputs go under `research_work/generated/timing-population`.
