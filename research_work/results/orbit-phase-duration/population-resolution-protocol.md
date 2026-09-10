# Supplement after the strict long-path failure

Both four-unit integrations failed the predeclared pointwise accuracy gate. Do not change that result. Test whether the solver disagreement also materially changes empirical orbit occupations. This supplement is declared after that failure and is not an originally reserved test.

Reproduce the 2e-11 run for each failed model because the original driver retained only the final 2e-13 trajectory. Save this otherwise redundant run so that distributions at the two accuracies can be compared. Reproduce the recorded maximum full-path differences to relative tolerance 1e-7 (absolute 1e-9) before using the new output.

For intervals (0,1], (0,2], (0,4] and (3,4], retain 2000 uniformly spaced states per interval at each accuracy. Report pointwise position/velocity differences, old-bin total variation between the two solver occupations, and exact Gaussian-kernel MMD at the two earlier resolutions. Also report early/late MMD at each solver accuracy. No new optimizer, numerical MMD pass threshold, physical fit or equilibrium declaration. Differences between two solvers are not a certified bound on error relative to the exact orbit.

Purpose: decide whether strict individual-path failure also spoils the population statistic. If population scores are robust here, develop an appropriate population convergence requirement; do not silently replace the failed primary gate with a newly convenient passing threshold.
