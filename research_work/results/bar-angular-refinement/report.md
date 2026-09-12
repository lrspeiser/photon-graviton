# Source-direction refinement

The source law, radial launch speed, fixed ordinary bar and integration horizon are unchanged. This compares deterministic angular quadratures, not different physical theories or fresh observational fits. Each launch sphere is normalized separately.

## Population fractions at approximately 244 Myr

| Launch radius (kpc) | Angular grid | Cohort ever entered (%) | Continuous population ever entered (%) | Continuous population inside now (%) |
|---:|---|---:|---:|---:|
| 1 | 6x12 | 86.209 | 63.615 | 0.5141 |
| 3 | 6x12 | 19.191 | 9.681 | 0.0466 |
| 1 | 8x16 | 80.348 | 54.367 | 0.4672 |
| 3 | 8x16 | 20.741 | 12.241 | 0.0389 |
| 1 | 12x24 | 67.197 | 48.609 | 0.3907 |
| 3 | 12x24 | 25.910 | 13.376 | 0.0583 |

All 208 new orbit status checks pass: **True**. Every new orbit is checked for Jacobi drift. Tighter-tolerance pairs cover every eighth direction (and any invariant failure), not every new orbit. The two grids add64 and144 trajectories across the two spheres. L40 is used throughout; prior L40/L64 agreement on the6x12 grid is not a new field-order test for every trajectory.

10 of 36 adjacent-grid statistical comparisons fail the predeclared development gates. The complete differences, including relative residence changes, are in summary.json. These are not confidence intervals.

| Comparison | R (kpc) | T | Failed statistic | Absolute change | Limit |
|---|---:|---:|---|---:|---:|
| 6 to 8 | 1 | 0.05 | cohort_ever_entered | 0.135897 | 0.05 |
| 6 to 8 | 1 | 0.1 | cohort_ever_entered | 0.233656 | 0.05 |
| 6 to 8 | 1 | 0.1 | continuous_ever_entered | 0.069427 | 0.05 |
| 6 to 8 | 1 | 0.25 | cohort_ever_entered | 0.058606 | 0.05 |
| 6 to 8 | 1 | 0.25 | continuous_ever_entered | 0.092472 | 0.05 |
| 8 to 12 | 1 | 0.05 | cohort_ever_entered | 0.053457 | 0.05 |
| 8 to 12 | 1 | 0.1 | cohort_ever_entered | 0.062156 | 0.05 |
| 8 to 12 | 1 | 0.25 | cohort_ever_entered | 0.131506 | 0.05 |
| 8 to 12 | 1 | 0.25 | continuous_ever_entered | 0.057585 | 0.05 |
| 8 to 12 | 3 | 0.25 | cohort_ever_entered | 0.051682 | 0.05 |

## What this means

A central encounter flag jumps at a grazing orbit, whereas time inside tends to zero. Consequently the encounter statistic can be harder to integrate over directions. Both matter for describing trajectories, but instantaneous gravity follows the actual spatial distribution, not the fraction with a past visit. A pass at the loose0.005 absolute residence threshold does not establish a small relative uncertainty.

At the final epoch,8-to12 residence changes are19.6 percent for the1kpc source and33.2 percent for the3kpc source relative to the finer values. In contrast, their total angular source normalizations change only0.0000904 percent and0.277 percent. Thus source normalization can be well integrated while a trajectory-dependent statistic remains poorly resolved. Do not repair this by tuning the injection amplitude.

The [population-to-gravity derivation](population-to-gravity.md) uses that spatial distribution directly. It applies known transport and Newtonian formulas under the declared cold-particle postulate; none of those standard formulas is claimed unique. Source amplitude and physical capture funding remain unspecified. The companion/time redshift mechanism is not tested here.

Next assess the three-dimensional force coefficients under direction refinement. Even agreement of enclosed central mass does not establish agreement of local forces in a nonspherical distribution. A self-gravitating formation calculation, full ordinary-matter field and joint observational tests remain required. All nine goals remain active.

Reproduce with run.py8 and run.py12 (pass the number as a separate argument), followed by export.py and document.py. Raw data retain primary and tighter-tolerance trajectories. No discarded failing direction, tuned source parameter or opened holdout is involved.
