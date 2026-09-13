# Frozen follow-up numerical refinement

The inherited calibration protocol already requires refining the largest absolute exponent error in each of eight shape/SNR/truth cells. Apply that rule to the completed fresh-seed batch, with deterministic label ordering for ties. Retain all eight results, including failures, and do not replace original cases in coverage estimates.

refine.py uses 8192 shape/peak integration samples and independent scramble seed 1402, versus the base 2048 samples and seed 1401. The 321-point width grid remains unchanged; this is not a width-grid convergence test. Synthetic fluxes are regenerated from identical seeds and cadence slots. The first event in each case is recomputed under the old settings and checked against its saved likelihood. Input, code and array hashes are checked, and resume requires unchanged hashes.

Use the inherited gates: absolute fitted-b change <=0.05 and mean posterior-weighted centered-event-curve difference <=0.1. The latter integrates the absolute difference against the refined event likelihood and fitted width population, using the boundary-capable integrator. No observed flux values or final holdouts are used. This is selected numerical sensitivity, not a new unbiased coverage sample, a global optimization proof, or evidence for a physical propagation cause.

The runner checkpoints completed cases under research_work/generated/timing-fresh-refinement and writes refinement-results.json only after all eight finish. Run `python research_work/results/timing-fresh-calibration/refine.py`. All six physical objectives remain open.
