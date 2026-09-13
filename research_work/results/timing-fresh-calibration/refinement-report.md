# Eight-case timing integration refinement completed

The previous turn verified a live process and audited six completed cases. This turn confirmed terminal exit code zero, all eight required cases, and the complete saved-result audit. All six physical research objectives remain open.

The inherited selection rule picked the largest absolute fitted-exponent error in each shape/SNR/truth cell from the completed 160-case fresh-noise calibration. The runner and protocol were committed before outcomes at 6e79203. Each selected case retains all 98 exposed cadence slots and its original synthetic noise realization. Shape/peak integration increased from 2048 to 8192 Sobol samples, with a different scramble seed. The 321-node width grid was unchanged.

| Shape | SNR | True b | Seed | Absolute b change | Weighted event-curve change | Gate |
|---|---:|---:|---:|---:|---:|---|
| Split Gaussian | 5 | 0 | 1913 | 0.001974 | 0.004121 | pass |
| Split Gaussian | 5 | 1 | 1902 | 0.000545 | 0.003660 | pass |
| Split Gaussian | 20 | 0 | 1918 | 0.000646 | 0.092857 | pass |
| Split Gaussian | 20 | 1 | 1918 | 0.004531 | 0.063520 | pass |
| Shoulder | 5 | 0 | 1913 | 0.002444 | 0.005565 | pass |
| Shoulder | 5 | 1 | 1908 | 0.002574 | 0.005150 | pass |
| Shoulder | 20 | 0 | 1918 | 0.003022 | 0.086835 | pass |
| Shoulder | 20 | 1 | 1918 | 0.000491 | 0.064483 | pass |

The fixed limits were absolute b change <=0.05 and posterior-weighted mean centered-event-likelihood change <=0.1. All pass; the high-SNR curve metrics have less margin. Increasing sample count and changing scramble together tests sensitivity to that combined change, not an isolated estimate of either error contribution.

## Audit evidence

refinement-results.json contains fits, thresholds, case selection, source hashes and saved-array hashes. refinement-audit.json checks all eight cases against one atomic checkpoint snapshot. Input/code hashes, original and refined array hashes, expected selection, dimensions, optimizer status, parameter bounds and the probability outside width support pass. Every best objective recomputes with zero discrepancy, and both numerical metrics and gate decisions reproduce. Original synthetic generation was checked by recomputing each case's first event at the original integration settings before refinement.

The audit shares the population integration implementation and does not independently prove it correct or prove global optimization. Existing first-case progress and audit files are historical snapshots; they do not imply the process is still running. The execution handle returned terminal success for this batch.

## Meaning and next tests

This provides evidence that the fitted exponent is insensitive to the tested shape-integration refinement in these selected difficult realizations. It is not width-grid convergence, precision coverage certification, realistic source-population validation, or a physical explanation of time dilation. The known relation log W=a+b log(1+z)+scatter remains a statistical model, not a newly derived propagation law.

Original cases and their coverage estimates remain unchanged. These outcome-selected refinements are not an independent coverage sample. Fixed-truth likelihood-ratio acceptance has not been recomputed at the refined settings in this runner and remains a follow-up diagnostic, alongside separate width-grid checks. Actual source evolution, filter response, brightness selection and a jointly predictive redshift-duration-brightness mechanism still require work. No observed flux or final observational holdout was used.
