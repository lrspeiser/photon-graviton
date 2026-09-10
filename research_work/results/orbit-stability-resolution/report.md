# Orbit-population drift survives a smooth-resolution check

The first and second halves of the existing orbit populations differ under a smooth position-and-velocity comparison as well as the previous histogram comparison. Thus the earlier finding is not solely an artifact of bin edges. However, the large histogram discrepancies are not percentages of gravitational error or evidence against either potential. Neither metric is calibrated here to observational uncertainty.

There is a second actionable finding: weights optimized for the previous histogram sometimes worsen the smooth comparison. Improving one numerical score alone is insufficient to declare the orbit population ready for a stellar-motion likelihood.

## Why this check matters

An orbit library represents a population by spending time along each trajectory. If its predicted distribution depends strongly on which part of the integration we use, a later comparison with measured stars can confuse incomplete orbital sampling with an incorrect gravity law. Conversely, bins with sharp boundaries can exaggerate differences between nearby states. We tested that possibility before spending more time on another library expansion.

The previous goal turn was progress: it completed and published the receiving-energy release audit. This turn returns to the independent observational requirement. No physical coefficient or star-by-star prediction was adjusted to improve agreement.

## Fixed inputs and methods

The protocol was written before calculating these new values. Inputs are the eight common passing ordinary-matter and candidate-potential trajectories with launch indices 0, 8, 12, 16, 20, 24, 28 and 32. Each contains 2000 noninitial samples spanning about 978 Myr. The previously failed index 4 remains excluded in both models; this audit does not repair it or claim the passing subset is representative. Initial starts were selected through candidate-based training coverage, so this is not a neutral comparison of the two gravity models.

Two partitions have exactly 1000 samples per orbit on each side:

- **Chronological:** the first half against the second half of the integration.
- **Interleaved control:** alternating samples across the entire integration. This checks sampling sensitivity; neighboring samples are strongly correlated, so a small value is not an independent equilibrium test or a valid statistical null distribution.

Weights are either equal across all eight trajectories, or the previously fitted longest-window histogram weights capped at 0.25 per orbit. Those weights are frozen, not reoptimized for this calculation. No new observed stellar outcomes, validation data or final holdouts are used.

The known total-variation histogram statistic is reused without changing the six-dimensional bins. The smooth comparison uses the known Gaussian-kernel maximum mean discrepancy (MMD), described in [Gretton et al., A Kernel Two-Sample Test (2012)](https://www.jmlr.org/papers/v13/gretton12a.html). We use its empirical-distribution distance, not that paper's independent-sample significance tests.

For rotating Cartesian coordinates x and the corresponding momentum-per-unit-mass components p, the known Gaussian kernel is

    k(a,b) = exp[-0.5 (|x_a-x_b|^2 / s_x^2 + |p_a-p_b|^2 / s_p^2)].

For the signed weighted empirical distribution d = P_early - P_late, the known quadratic identity is

    MMD^2 = sum_i sum_j d_i d_j k(i,j).

The implementation includes self-pairs and computes exact finite-sample kernel sums in blocks. The result is nonnegative within floating-point precision and is not an unbiased population estimator. MMD is dimensionless; it is not a percentage and is not numerically comparable to total variation. The illustrative resolutions are (0.5 kpc, 50 km/s) and (2 kpc, 100 km/s). They are neither observational error bars nor fitted parameters. These formulas are established statistics, not a new gravity law.

## Results

| Potential | Weights | Position / velocity scale | Chronological MMD | Interleaved MMD |
|---|---|---|---:|---:|
| Ordinary matter | Equal | 0.5 kpc / 50 km/s | 0.055722 | 0.001134 |
| Ordinary matter | Frozen histogram fit | 0.5 kpc / 50 km/s | 0.088495 | 0.000351 |
| Ordinary matter | Equal | 2 kpc / 100 km/s | 0.069545 | 0.000272 |
| Ordinary matter | Frozen histogram fit | 2 kpc / 100 km/s | 0.052446 | 0.000307 |
| Candidate extra potential | Equal | 0.5 kpc / 50 km/s | 0.072586 | 0.000250 |
| Candidate extra potential | Frozen histogram fit | 0.5 kpc / 50 km/s | 0.115737 | 0.000349 |
| Candidate extra potential | Equal | 2 kpc / 100 km/s | 0.074074 | 0.000270 |
| Candidate extra potential | Frozen histogram fit | 2 kpc / 100 km/s | 0.132830 | 0.000360 |

All chronological values exceed their corresponding interleaved controls. These ratios are descriptive, not confidence levels: the interleaved control is expected to look similar because it samples almost the same neighboring portions of each path.

The previous histogram values reproduce exactly at the displayed precision:

| Potential | Weights | Chronological TV | Interleaved TV |
|---|---|---:|---:|
| Ordinary matter | Equal | 0.780625 | 0.043500 |
| Ordinary matter | Frozen histogram fit | 0.612500 | 0.062250 |
| Candidate extra potential | Equal | 0.667250 | 0.054875 |
| Candidate extra potential | Frozen histogram fit | 0.560844 | 0.038203 |

For example, the candidate's fitted weights reduce chronological histogram TV from 0.667250 to 0.560844, while increasing fine-resolution MMD from 0.072586 to 0.115737 and coarse-resolution MMD from 0.074074 to 0.132830. This is a concrete warning against equating optimization of the old histogram with general population convergence. The ordinary-matter fit also worsens the fine smooth comparison, although it improves the coarse one.

For the candidate, index 12 has the largest individual chronological MMD at both resolutions (0.398119 and 0.555384). This is a diagnostic priority for studying phase coverage and integration duration, not grounds for deleting that orbit after seeing an unfavorable result.

## Verification

Run `python research_work/results/orbit-stability-resolution/run.py` from the repository root. Checks cover upstream hashes, previous capped longest-window histogram values, Gram symmetry and positive-semidefiniteness, identical distributions, a known two-point Gaussian-kernel distance, and a direct small two-orbit matrix calculation independent of the blocked reduction. Source hashes are checked again after computation. Full values, fixed weights and Gram matrices are in `results.json`; the protocol and script hashes are recorded there.

No integrations were rerun, new stars selected, physical parameters fitted or likelihoods evaluated. The apparent numerical precision describes reproducibility of these cached trajectories, not astrophysical certainty.

## What changes next

The finite-time population remains insufficiently demonstrated for a strong galaxy test. Merely making bins larger is not a defensible repair: a smooth metric still sees early/late changes, and shrinking all differences through extreme smoothing would conceal useful structure.

Retain both histogram and smooth diagnostics during library development. Next investigate phase/duration coverage for the most drifting trajectories, preserving difficult cases and using the same procedure for both potentials. A ready observational pipeline must then predict selected populations with distance and motion errors, ordinary-matter uncertainty and survey selection. Its relevant convergence tolerance must be connected to changes in those predicted observables, not chosen from MMD or TV after viewing the outcomes.

This report does not rank the gravity models, prove physical nonequilibrium, establish a stationary distribution or derive companion deposition. The source/receiver mechanism, source and detector clocks, capture/support, lensing and deferred total energy supply remain independent open requirements. The full research goal remains active; there is not yet a strong case for the unified theory.
