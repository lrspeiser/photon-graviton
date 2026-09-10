# Can orbit weights produce a stable population?

Optimizing the mixture cannot make this eight-orbit library stable in joint position and velocity at the tested resolution. Over the longest available run alone, the smallest early/late total-variation distances are 0.448 for ordinary matter and 0.491 for the empirical extra-potential candidate. Those are finite-bin numerical bounds for these paths, not discrepancies with observed stars or a ranking of the gravity models.

Coarse radius/height distributions can nevertheless pass the illustrative 0.05 stability target. This explains why a stable-looking spatial summary is insufficient for the proposed bulge/disk comparison: velocity and angular-position information must be checked together.

## Formula provenance and calculation

All mathematics in this test is known: nonnegative mixtures, normalized histograms, total variation, and linear programming. No new force or time formula is introduced. The candidate potential is the existing empirical conservative extra response; the photon-deposition derivation is still missing.

For orbit k, let A_jk and B_jk be its normalized early/late occupation histograms for duration j. The mixture has one common weight vector w, with w_k >= 0 and sum(w_k) = 1. The diagnostic is:

    TV_j(w) = 0.5 * sum_bins |[(A_j - B_j) w]_bin|
    primary objective: minimize over w the maximum TV_j(w)

The primary objective uses the first 244, 489 and 978 Myr together. A secondary check, declared after the primary result, uses only the early/late halves of the 978-Myr run. It distinguishes failure to agree with short approximations from failure within the longest interval. The secondary check does not replace the original result.

The eight starts are the common numerically passing pairs [0,8,12,16,20,24,28,32]. The previous failed pair 4 stays excluded from this numerical use and its failure remains archived. These starts were selected during candidate training-library development; they are neither representative weights for all Galactic stars nor a neutral model-selection sample.

## Best achievable stability for these paths

Smaller TV means closer early/late binned fractions. The 0.05 target is an explicit numerical diagnostic choice, not an observational error bar. A TV of 0.49 does not mean 49% of real stars disagree with the theory.

| Model | Information retained | Maximum individual weight | Best worst TV across 3 durations | Best TV, longest run only | Longest-only TV, half sampling |
|---|---|---|---:|---:|---:|
| Ordinary matter | Radius and signed height | Unrestricted | 0.1140 | 0.0209 | 0.0241 |
| Ordinary matter | Radius and signed height | 0.25 | 0.2539 | 0.0359 | 0.0373 |
| Ordinary matter | Radius, height, bar angle | Unrestricted | 0.4560 | 0.1786 | 0.1843 |
| Ordinary matter | Radius, height, bar angle | 0.25 | 0.7625 | 0.2385 | 0.2385 |
| Ordinary matter | Position and 3 velocities | Unrestricted | 0.9512 | 0.4480 | 0.4500 |
| Ordinary matter | Position and 3 velocities | 0.25 | 0.9690 | 0.6125 | 0.6140 |
| Empirical extra potential | Radius and signed height | Unrestricted | 0.0537 | 0.0335 | 0.0406 |
| Empirical extra potential | Radius and signed height | 0.25 | 0.0545 | 0.0415 | 0.0448 |
| Empirical extra potential | Radius, height, bar angle | Unrestricted | 0.1991 | 0.1550 | 0.1579 |
| Empirical extra potential | Radius, height, bar angle | 0.25 | 0.4509 | 0.1998 | 0.2041 |
| Empirical extra potential | Position and 3 velocities | Unrestricted | 0.5174 | 0.4910 | 0.4840 |
| Empirical extra potential | Position and 3 velocities | 0.25 | 0.6882 | 0.5608 | 0.5575 |

With unrestricted weights, the optimum for joint position/velocity in the longest interval concentrates on a single orbit in each model. Even this permissive choice cannot get near 0.05. The optional 0.25 weight cap forces at least four nonzero weights and makes the mismatch larger: 0.6125 ordinary and 0.5608 candidate. That cap is a diagnostic control, not a physical prior.

In contrast, longest-interval radius/height TV is 0.0209 ordinary and 0.0335 candidate with unrestricted weights, and stays below 0.05 with the cap. Marginalizing over bar angle and velocity hides structure that remains time-dependent in the more complete distribution. These are separately optimized partitions; their fitted weights need not be the same.

## Numerical evidence for the bounds

The 48 optimization runs all pass primal and dual feasibility checks and independent reconstruction of the achieved total variations. The largest primal/dual objective gap is 6.33e-15; the largest constraint/stationarity/sign residual is 1.68e-14. These support numerical optimality for the stated finite matrices, not an interval-arithmetic theorem or a universal bound over all possible orbits.

Six analytic controls check identical distributions, common non-cancelling drift, and two opposing drifts that cancel in a mixture, each with and without the weight cap. Every-second-time sampling preserves the qualitative failure in joint position/velocity. Both refitted and frozen-weight half-sampling results are retained, so reoptimization is not hiding sampling sensitivity.

The inherited spatial bins cover radius, signed height and eight bar-angle sectors. The joint test adds the declared coarse radial, rotational and vertical velocity bins. Normalized histogram columns are checked before optimization. Sources and trajectory hashes are verified against the paired-duration audit. All weights, effective orbit counts, objective bounds and window-specific distances are saved in results.json.

## What follows, and what does not

Reweighting these eight paths cannot repair the full finite-window instability at this resolution. A later fit to observed stars would impose additional constraints, so it cannot evade this lower bound while keeping exactly these trajectories, bins and time comparisons. This does not prove that larger libraries, different justified resolution, longer integrations, or a different population construction will fail. Finite sample noise and bin choices remain part of the limitation.

This step implements and tests common population weights for numerical stability; it does not implement the survey-selected stellar likelihood. No observed number-density target, velocity distribution, mass parameter or gravity parameter was fitted, and no validation/final-test outcomes were read. These numerically optimized weights must not be reused as an inferred Galactic population.

The next stellar step needs a population representation whose joint spatial/velocity predictions converge with duration and phase sampling, using equivalent flexibility for ordinary and candidate gravity. Survey selection, distance-prior recycling, measurement covariance and ordinary-matter uncertainty must then enter the shared likelihood. A radius/height-only comparison is insufficient for the bulge question.

At the broader theory level, a stable stellar fit would still not identify a photon origin for extra gravity. The redshift/event-timing/clock relation, physical energy-momentum completion, capture and supported deposits, lensing, and the deferred total photon-supply calculation remain necessary. This result changes the population-modeling next step; it does not resolve those physical requirements or complete the research goal.
