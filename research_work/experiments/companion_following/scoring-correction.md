# GF-1 scoring correction: complete declared lag search

A final code review found that best_correlation limited lag to half the stored
trajectory as well as 20 time units. A time-32 chain therefore searched only
through 16 rather than all delays through 20. This was an implementation
error, not a change in the declared scientific threshold.

Keep screen-v1, screen-completion-v1, refine-v1 and diagnostics-v1 unchanged.
Rescore saved chain trajectories through the full declared 20-time-unit window,
writing scoring-v2. Do not reintegrate those trajectories. Recompute selection
using the original selection rule, then execute refine-v2 and diagnostics-v2
for the corrected selections. Ring scores and equations are unchanged.
Keep the first 26 passing controls and add an independent synthetic lag-18
control that the old scorer fails. Run the controls again in controls-v2.

Initial counts of 218 phenomenological and 4 conservative chain passes were
reported before this correction and are superseded, not silently discarded.
Historical selections and their outcomes remain reproducible from their
recorded source commits. The new outcomes will be reported separately.

Correlation compares overlapping windows and is a heuristic of this exposed
fixture. Searching more lags increases the chance of a match; no p-value,
holdout success, finite signal speed, or astrophysical validity is implied.
