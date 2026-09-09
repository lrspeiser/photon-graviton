# T01 baseline reproduction

All three latest scripts completed successfully in an isolated copy. The 271 available original snapshot files match their recorded hashes; the five missing files belong to historical paper construction. The original scientific inputs and outputs were unchanged after the run.

Python/package versions are in environment-lock.txt. python-docx is absent from this interpreter, but none of these three numerical scripts needs it. No installation or source-code change was necessary for the scientific run. A setup-only metadata check initially stopped before any calculation; the runner was corrected to record optional missing packages.

The run compared 30252 numerical values plus labels and structure across eight saved artifacts. The comparison thresholds were recorded before execution (relative 1e-5, absolute 1e-8). There were 170 exceedances: seven summary/parameter values in results.json and 163 predicted velocities in cv_predictions.csv. The other six artifacts pass the original thresholds. Thresholds were not relaxed to make the discrepancies disappear.

The substantive exceedances trace to the mixed-retention model's fold-3 fitted parameters, its CV predictions and derived bootstrap summaries. The maximum change across all saved CV predictions is 0.019801501 km/s. Mixed-retention overall RMSE changes from 28.1594467 to 28.1602225 km/s; log RMSE changes from 0.111256356 to 0.111257726. The model ranking under the declared log metric is unchanged, and the paired interval endpoints retain their negative signs. Small other mixed-model numerical differences remain below the declared comparison policy.

This isolates the numerical disagreement; its exact platform/optimizer cause has not been proven. It is consistent with optimizer sensitivity and is not evidence of new physics. Reproduced energy budgets, controls, cluster transfer and circulation do not change the previously reported conclusions. Do not claim bit-identical reproduction of the entire first-fit output.

Observed run times: run.py 2.20s; followup.py 9.76s; circulation.py 0.19s. These are timings on this machine, not universal guarantees; the runner uses a 900-second per-script timeout. The dependency graph documents the shared loader and its protocol-writing side effect.

Verdict: T01 completed with a documented, localized numerical reproducibility difference. Retain this baseline and its tolerance report. No physical path is closed by this result.

Reproduction harness: run_baseline_001.py. It deliberately refuses to overwrite an existing run-status record; for a future independent rerun choose a new run/output directory. The first run is terminal and complete, not an active background job.
