# SE-R final report: direct-companion convergence trend fails

20 September 2026. All seven new runs are complete. The independent archive
audit passes50 checks, but the source-accuracy campaign fails the declared
decreasing-error requirement for the Y-channel spatial comparison. The small
absolute differences do not justify retroactively relaxing that requirement.

| Comparison | Wave-source relative change | Total-source relative change | Result |
|---|---:|---:|---|
|Mixed, n48/n56|0.758519%|0.00263718%|Pass|
|Y, n48/n56|0.0879441%|0.00007237%|Fail: trend|
|Mixed timestep|0.00000335%|0.00000000915%|Pass|
|Y timestep|0.00000204%|0.00000004041%|Pass|
|Mixed rotation|0.134938%|0.00077649%|Pass|
|Y rotation|0.153970%|0.00347369%|Pass|

The mixed wave-source difference falls from1.174305% at n40/n48 to0.758519%
at n48/n56. The Y difference instead rises from0.0584820% to0.0879441%.
Both latest differences are below the5% magnitude limit, but the second
does not establish the declared monotone convergence trend. There is no
approved source-accuracy pass for the whole campaign.

These source values are Hamiltonian derivatives under a uniform scalar
variation, with matter and field contributions included for the total. They
are not gravitational masses, measured accelerations, lensing multipliers or
observational fits. The campaign concerns the earlier six-field emitting
model, not the newly implemented joint local-transfer candidate.

All first-run arrays and traces remain in source-refinement-v1. The audit
independently reconstructs endpoint energies/source diagnostics, checks saved
traces, hashes and the declared comparisons; its result is in
source-refinement-audit.json. No dark matter or expansion was introduced,
and agreement with old gravity laws was not a success criterion. All twelve
broader goal requirements remain active.
