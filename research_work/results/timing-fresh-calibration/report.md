# Fresh-seed timing calibration completed

The preceding goal turn ruled out one proposed selection-code explanation on the available Milky Way parent subset. That unresolved survey issue remains open. Work now also advances objective 1 through independent-noise calibration of the existing timing estimator.

The protocol and runner were committed as 4811d90 before any new outcomes were examined. The batch uses seeds 1901–1920, two existing injection shapes, SNR 5 and 20, and true duration exponents b=0 and b=1: 160 cases of 98 events each. Seeds are independent within each cell; shared seeds across cells mean their outcomes are not independent of each other. The cadence slots and shape families are already exposed, so this is fresh-noise calibration, not an untouched observational test.

The unchanged event likelihood and boundary-capable population integrator are used, with the original estimator providing initialization. The population relation `log W = a + b log(1+z) + scatter` is a known statistical parametrization. It is not a derived photon–companion propagation law. Zero intrinsic scatter is permitted in fitting, as in the repaired estimator.

The frozen coarse screens require absolute mean exponent bias at most 0.1, at least 85% nominal-95% truth acceptance and no invalid fits in each 20-case cell. Exact binomial intervals will accompany the acceptance rates; twenty trials cannot demonstrate precise 95% coverage. No realization may be replaced because of an unfavorable result. Source hashes prevent resuming a batch with changed code or inputs.

The original process completed all 160 cases with exit code zero. The complete-result audit verifies every expected case, all frozen source and saved-array hashes, and recomputes objective values with maximum discrepancy zero. No observed flux values were used; no physical validation is claimed.

Subsequent numerical refinement, source-evolution/selection tests and a jointly derived redshift–duration–brightness law remain necessary. All six goals remain open.

## Running-result integrity audit

`python research_work/results/timing-fresh-calibration/audit.py` verifies a single atomic checkpoint snapshot, all frozen source hashes, each saved event-array hash, expected seed/cell membership, parameter bounds and probability outside the width support. It recalculates unrestricted and successful fixed-truth likelihood values from the arrays and checks the stored likelihood-ratio decisions. This checks recorded-result consistency; it does not independently prove a global optimum or correctness of the shared likelihood implementation.

The first audited snapshot contains two of 160 cases, with zero discrepancy in the recomputed best objective values. Both have valid fits. The bounded snapshot and its hash are recorded in audit-snapshot.json; later process progress may exceed it. A checkpoint file alone does not prove the process remains live. The auditor writes final cell summaries only when every expected realization is present. No incomplete calibration pass is reported.


## Complete batch outcome

All eight frozen coarse screens pass, with zero invalid fits. Nominal 95% truth acceptance is 18/20 for split-Gaussian SNR 5, b=1; 19/20 for split-Gaussian SNR 5, b=0 and shoulder SNR 5, b=1; and 20/20 in the other five cells. All four misses are retained. Cell mean exponent biases range from -0.04817 to -0.00047, below the preset absolute 0.1 limit.

The exact 95% binomial intervals are wide: 18/20 gives [0.6830, 0.9877], 19/20 gives [0.7513, 0.9987], and 20/20 gives [0.8316, 1]. Shared seeds correlate cells; do not treat the combined 160 outcomes as independent Bernoulli trials. Full per-cell biases, Monte Carlo standard errors and intervals are in summary.json. The initial two-case audit-snapshot.json remains a historical snapshot; the completed audit is recorded in summary.json.

This is a successful fresh-noise screening check of the repaired estimator, conditional on the existing injection families, cadence slots and likelihood assumptions. It is not precision coverage certification, a test of source-population evolution, or observational evidence for photon conversion. Numerical refinement on this batch and realistic source/selection modeling remain required before interpreting astronomical timing results. All six physical research objectives remain open.
