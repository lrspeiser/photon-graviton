# Can the timing estimator distinguish stretched from unstretched events?

This is a small, frozen synthetic control experiment toward the joint redshift, timing and brightness test. It uses 98 previously exposed DES cadence/error patterns, with newly seeded artificial fluxes. No measured supernova flux is fitted. This tests the measuring method, not the physical companion hypothesis.

## Why this test is needed

A method can give stable numbers and still favor the wrong answer. We therefore create events whose timing rule is known: b=0 gives no duration stretch; b=1 multiplies duration by 1+z. The estimator is not told which answer was injected. It fits observer-time light curves without dividing times by 1+z.

The law W proportional to (1+z)^b is a known phenomenological statistical parameterization, not a new physical derivation. An intrinsic source-duration trend can mimic b; these controls deliberately have no such trend.

## Frozen design

Three seeds (701, 702, 703) each supply paired b=0 and b=1 cases. The paired cases share intrinsic-width draws, peak offsets and noise draws, isolating the change of propagation exponent. The source has an asymmetric Gaussian shape plus a shoulder, outside the estimator family. Intrinsic log-width scatter is 0.1 and nominal peak-to-median-error ratio is 20. No cadence slots are removed.

Each case uses 2048 Sobol shape samples and 321 duration-grid nodes. Before running, we specified refinement of the case with the largest absolute error from its injected exponent to 8192 samples with an independent scramble. This diagnostic case selection is not a validation split; it cannot guarantee all cases have converged.

The small-control thresholds are absolute median bias at most 0.15 for each truth, median b separation at least 0.7, successful interior fits with negligible width-boundary probability, and refinement changes at most 0.05 in b and 0.1 in the mean centered event log-likelihood. These are engineering acceptance criteria, not a statistical significance claim or thresholds for observational agreement.

| Injected b | Seed | True finite-sample slope | Recovered b | Difference from injected b |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 701 | 0.04817 | 0.11248 | +0.11248 |
| 0 | 702 | -0.12803 | -0.11262 | -0.11262 |
| 0 | 703 | -0.02669 | -0.00608 | -0.00608 |
| 1 | 701 | 1.04817 | 1.10971 | +0.10971 |
| 1 | 702 | 0.87197 | 0.90989 | -0.09011 |
| 1 | 703 | 0.97331 | 1.02213 | +0.02213 |

The true finite-sample slope differs from the injected population exponent because randomly drawn source widths happen to correlate with redshift. It is diagnostic information available only for artificial data, not a correction applied to real observations.

Median recovered b for no stretch: -0.006076578253726324. For injected stretch: 1.0221310524496083. Median separation: 1.0282076307033348.

## Numerical check and outcome

Refined case: b0-seed702. Change in b: 0.001343024059719311. Mean centered event log-likelihood change: 0.08211356826413518. Both estimates are retained in the result file.

**Small-control gate: PASS.** Median-bias gate: True; refinement gate: True; optimization/boundary gate: True.

## What this does not establish

Three realizations per truth cannot establish calibrated uncertainty coverage or broad robustness. The samples use one shoulder shape, one signal level, no intrinsic evolution, and simplified band matching. Survey selection, actual filter transmission, source population variation, amplitude-prior sensitivity and additional noise levels remain necessary. Artificial amplitude normalization removes the real brightness-selection problem.

Because each event amplitude is marginalized freely, this test supplies no evidence for the physical brightness law. Joint redshift/timing/brightness still needs calibrated flux predictions, source luminosity/distance information and a common transport mechanism. None of the six scientific demonstrations is marked complete.

## Reproduce

```powershell
python research_work/results/timing-null-stretch-controls/run.py
python research_work/results/timing-null-stretch-controls/summarize.py
```

## Residual estimator shift

Mean recovered-minus-true-finite-sample slope for b=0: +0.03344.
Mean recovered-minus-true-finite-sample slope for b=1: +0.04943.

Every base-fit recovered slope exceeds its corresponding true finite-sample slope in these six controls. The small-control median criterion can pass despite a residual positive shift. Three paired samples cannot establish its general size or origin; subsequent source-shape, prior and selection sensitivity tests must address it rather than treating this pass as an unbiasedness guarantee.
