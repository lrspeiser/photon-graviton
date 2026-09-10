# Frozen diagnostic: distinguish histogram resolution from smooth population drift

Use the eight existing common-passing ordinary/candidate 978-Myr trajectories, with no new stellar outcomes or integrations. Preserve all previous numerical failures and selection caveats. Before reading the new diagnostic values, fix the following choices.

1. Compare the first and second 1000 noninitial time samples. As a sampling-resolution control, compare alternating samples across the same whole trajectory. Alternating samples are strongly dependent and do not supply independent evidence of equilibrium or a statistical null distribution.
2. Compute the previous six-dimensional histogram total variation for both splits, plus a Gaussian-kernel mean-discrepancy measure in Cartesian position and momentum. Use illustrative scales (0.5 kpc, 50 km/s) and (2 kpc, 100 km/s); these are diagnostic resolutions, not measurement uncertainties or fitted bandwidths.
3. Retain equal orbit weights and the previously fitted capped longest-window position/velocity weights for each model, unchanged. Do not optimize the new smooth metric.
4. Compute exact finite-sample kernel sums in blocks, including self-pairs; label these discrepancies of empirical distributions, not unbiased population estimators or significance levels. Report each individual orbit as well as both mixtures.
5. Check symmetry, positive-semidefiniteness within floating error, equality of zero-displacement samples, a known two-point analytic case and unchanged source hashes. The calculations may soften or reinforce the earlier histogram finding, but cannot establish stationarity, observed likelihood or model preference.

Purpose: if a fine histogram overstates drift relative to smooth measurements, improve the convergence diagnostic before treating its failure as a need for a larger orbit library. If smooth drift also survives, merely changing bins cannot repair the population. No target is tuned and no numerical discrepancy is translated into an observational rejection threshold.
