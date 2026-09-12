# Supernova timing: shape-integration refinement

This calculation repairs a numerical prerequisite for the joint redshift, timing and brightness demonstration. It uses artificial light curves at the same 98 previously exposed DES cadence patterns. It does not measure real supernova timing, test brightness, or validate the companion mechanism.

## What changed and why

The previous pilot failed its numerical probability-accuracy threshold. We kept its synthetic seed, source shape, noise, priors and thresholds fixed, and increased shape/peak integration from 512 samples to 2048 and 8192. Both comparisons now use the same 321-node observer-time width grid, isolating shape integration. An independent 8192-point scramble checks sensitivity to the quadrature realization. The protocol was written before this run. The original failed result is retained.

The likelihood formulas are known statistical marginalization and numerical quadrature, not new physical equations. Times are not divided by redshift. All events remain included.

| Calculation | Recovered timing exponent b | Intrinsic log-width scatter |
| --- | ---: | ---: |
| base | 1.104584 | 0.108320 |
| refined | 1.103855 | 0.107937 |
| independent_scramble | 1.104743 | 0.107868 |

Injected exponent: 1. The finite-sample true-width slope is 1.0468366922915948. Agreement in one realization cannot establish population bias or uncertainty coverage.

| Comparison | Change in b (limit 0.05) | Mean event log-likelihood change (limit 0.1) |
| --- | ---: | ---: |
| base to refined | 0.000728 | 0.053124 |
| refined to independent_scramble | 0.000888 | 0.035619 |

**Original numerical pilot gate: PASS.** No threshold was relaxed. The gate also requires successful interior optimization and negligible population probability outside the width bounds.

Event-level diagnostics are saved separately. The gate concerns a mean over events; passing does not assert every individual event meets 0.1.

## What remains before interpreting observations

Separate width-grid and prior-boundary sensitivity, repeated injections of both b=0 and b=1 with different source shapes and signal levels, calibrated uncertainty coverage, actual filter transmission, source evolution and survey selection remain unverified. The estimator integrates an arbitrary brightness amplitude for each event, so this timing calculation cannot itself test the modelÃ¢â‚¬â„¢s brightness prediction. A joint test must retain physical flux calibration and independently constrain source luminosities and distances.

If this numerical gate fails, local adaptive integration or further convergence work is needed before observational inference. If it passes, proceed to the remaining convergence and injection checks; do not describe the model as matching real supernovae.

## Reproduce

From the repository root:

```powershell
python research_work/results/timing-population/run.py --protocol research_work/results/timing-quadrature-refinement/protocol.json --output research_work/results/timing-quadrature-refinement
python research_work/results/timing-quadrature-refinement/check_width_grid.py
python research_work/results/timing-quadrature-refinement/summarize.py
```

## Individual-event limits

| Comparison | Events above 0.1 | Largest change |
| --- | ---: | ---: |
| base to refined | 14/98 | 0.274472 |
| refined to independent_scramble | 5/98 | 0.155684 |

The mean-based gate passes even though some individual curves remain less stable. These events are retained and should receive attention in subsequent sensitivity checks; no individual-event accuracy guarantee is claimed.

## Separate duration-grid check

The nested 161-to-321-node check uses identical refined shape samples. Change in b: 1.8003915e-07; mean event log-likelihood change: 0.034599. Original limits: 0.05 and 0.1. Result: PASS. This does not check a still finer 641-node grid.

Reproduce with `python research_work/results/timing-quadrature-refinement/check_width_grid.py` before running this summary script.
