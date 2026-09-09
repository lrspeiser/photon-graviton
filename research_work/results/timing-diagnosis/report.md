# Why the first timing estimator failed

## Purpose and result

The first injection test rejected too many stretched, shoulder-shaped artificial events. This follow-up uses their known true durations to separate sampling variation, acceptance selection and fitted-width error. The goal is to identify the repair needed before a real-flux timing inference, not to reinterpret the failed test as a pass.

The diagnosis shows that width estimation itself can contribute more to the recovered timing slope than the rejected-event selection. Improving coverage handling alone would therefore be insufficient. No new estimator, threshold or observational exponent is fitted here.

## Exact attribution

For each previously recorded run, compute three ordinary least-squares slopes against log(1+z), each with an intercept:

    b_true_all      from all injected true log widths;
    b_true_accepted from true log widths of accepted events;
    b_fit_accepted  from their fitted log widths.

Then the following identity holds:

    b_fit_accepted - b_injected
       = (b_true_all - b_injected)
       + (b_true_accepted - b_true_all)
       + (b_fit_accepted - b_true_accepted).

The terms respectively describe finite-sample intrinsic variation, selection of accepted events, and fitted-width error on that same accepted sample. This is an algebraic attribution, not a general causal decomposition of every survey systematic. Interactions between the fitter and acceptance are included in the way the accepted sample is defined.

The artificial widths were reconstructed from the original shape functions and random seeds, including the same consumption of random noise draws. Twenty-four individual fits were replayed to verify the reconstruction, and all twelve original recovered slopes were recovered from the saved widths. The sum above closes to numerical precision.

## What the decomposition reveals

Selected stretched-event trials, all with injected b=1:

| Shape and seed | Intrinsic sample contribution | Acceptance contribution | Fitted-width contribution | Recovered b |
| --- | ---: | ---: | ---: | ---: |
| Gaussian, 501 | +0.0421 | +0.0551 | +0.3529 | 1.4501 |
| Gaussian, 502 | -0.0940 | -0.0558 | +0.1384 | 0.9886 |
| Shoulder, 501 | +0.0421 | +0.0501 | +0.4312 | 1.5235 |
| Shoulder, 502 | -0.0940 | -0.0911 | +0.1961 | 1.0110 |

The apparently accurate shoulder result near 1.011 is a cancellation of opposing errors. It does not show that its width estimates are unbiased. In all six stretched trials, fitted-width errors add positively to the slope, by approximately 0.138 to 0.431. The unstretched trials have a different error pattern. These statements describe this limited injected set, not a measured bias of all possible data or fitters.

Even a good median width is insufficient. In shoulder trial 501 with b=1, the median fitted/true width ratio is about 1.018, yet width errors change the timing exponent by about 0.431. Errors at particular redshifts have different leverage in a regression. The result file lists the three largest absolute fitting contributions per trial for diagnosis; it does not authorize removing those events from a later real-data analysis.

## Coverage and fitted coverage are not the same

For the three stretched-shoulder trials, the true half-maximum endpoints lie inside the observed time range for 81, 85 and 81 of the 98 events. The original estimator accepted only 71, 77 and 70. Thus 10, 8 and 11 events were rejected despite geometric coverage of the true endpoints.

This does not prove those widths can be accurately recovered: sparse points between the endpoints, noise and shape mismatch still matter. It does show that the earlier acceptance failure cannot be attributed solely to the survey having no temporal coverage. Conversely, a few other trials accepted an estimated width even when a true endpoint lay outside the observed range. Cuts based on fitted crossings are imperfect proxies for actual coverage.

The true injected widths are available only in simulation. Applying these truth-based checks directly to real supernovae would be impossible; a revised estimator must propagate uncertainty and test its selection behavior without access to hidden truth.

## Consequences for the next estimator

The next protocol must address both uncertain widths and incomplete coverage. Possible methods include fitting observed fluxes jointly with a flexible shape distribution, or carrying width likelihoods and censored intervals into the population fit. A new method must be tested against shapes other than its own fitting function.

Retain the original frozen gate and failed trials. In a revised test, examine recovery on the full injected population and on the retained population separately, including width-error dependence on redshift. Also check whether apparently accurate aggregate slopes result from opposing biases. More seeds, varying signal-to-noise, realistic wavelength matching and selection sensitivity are still required before a scientific timing claim.

A method must not be tuned by dropping the known influential events, using true widths as real-data corrections, or changing the original pass threshold after observing the failure. The diagnosis guides a new recorded method; it is not an independent validation sample.

## Reproduction and scope

With the verified DES cache and timing-injection dependencies installed:

```sh
python research_work/results/timing-diagnosis/diagnose.py
```

The script reads historical injection results and produces a new diagnosis under `research_work/generated/timing-diagnosis`. It does not modify those results or use the real photometric flux column. It records source hashes, all twelve attributions and twenty-four replay checks. No default physics job is added; the offline suite remains at 31 jobs. The estimator feasibility gate and the overall theory remain incomplete.
