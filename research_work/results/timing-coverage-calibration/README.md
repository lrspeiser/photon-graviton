# Expanded synthetic timing calibration

This batch advances the joint-light test by checking estimator bias and nominal interval coverage across 160 artificial samples: 20 seeds in each of eight source-shape, signal-level and timing-law combinations. It uses the 98 previously exposed DES cadence/error patterns. No observed flux is fitted, and passing cannot validate the companion mechanism or its brightness law.

The protocol is fixed before results. Cases use plain asymmetric Gaussian or shoulder-bearing light curves, nominal signal-to-error ratios 5 and 20, and duration exponents 0 and 1. The same seeds are paired across cells; coverage is summarized per cell, not pooled as 160 independent trials. The worst absolute exponent error in each cell is subsequently refined with 8192 shape samples and an independent scramble. All base results remain in frequency estimates.

Coverage at the injected exponent can be measured by profiling the likelihood at that exponent and comparing to the nominal threshold. Full interval endpoints are not needed for this membership calculation. These are known statistical constructions, not new physical postulates. Twenty trials per cell give broad binomial uncertainty; they can reveal large calibration failures but cannot establish precise 95% coverage. Invalid/boundary fits count as failures to cover and are separately counted.

## Reproduce and resume

From the repository root:

```powershell
python research_work/results/timing-coverage-calibration/run.py
python research_work/results/timing-coverage-calibration/summarize.py
```

Do not launch a second copy while the first process is live. The execution tool's process/session is authoritative; a checkpoint alone does not establish that work is still running. After a confirmed stop, the same command resumes hash-verified completed cases. It refuses to combine changed code/protocol/input versions. A case interrupted before checkpoint completion is recomputed.

Arrays and the atomic checkpoint are stored in ignored `research_work/generated/timing-coverage-calibration`. On successful batch completion the runner writes `results.json` here; the summary command then generates `report.md`. The checkpoint is resumability data, not an observational result. The protocol, source and retained input hashes allow reproduction if the generated cache is lost.

This expanded run does not yet incorporate real bandpasses, intrinsic duration evolution, survey selection, or physically calibrated luminosities and distances. All six scientific demonstrations remain incomplete.
