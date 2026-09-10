# Longer phase sampling helps the candidate orbit, but does not certify a population

Extending the selected trajectory from about 0.98 to 3.91 Gyr reduces the candidate potential's smooth early/late occupation discrepancy. Its overall occupation is also much less sensitive to solver tolerance than its late instantaneous position. This supports investigating longer phase averaging, while retaining explicit numerical uncertainty. It is not a successful galaxy fit, an equilibrium proof or a comparison favoring one gravity theory.

The ordinary-matter path develops substantially larger late-time numerical sensitivity. Both models fail the originally declared strict four-unit pointwise trajectory gate. That failure is preserved. A separate, clearly labeled follow-up measures population sensitivity instead of silently replacing the original gate with a convenient new pass criterion.

In plain language, two calculations can disagree about exactly where a star is at a particular distant time while agreeing much more closely about where it spends most of its time. The latter matters when predicting a whole population. We need to check both rather than demand exact long-term positions without examining the prediction we actually intend to use.

## Fixed experiment

The previous goal turn was progress: it published the smooth-resolution audit, which identified launch index 12 as the largest individual candidate drift case at both fixed resolutions. This is the reason for selecting it now. The choice is targeted, not random, and must not be used to infer typical Galactic behavior.

Both models start from the same previously used training-star state, source ID 4064133698926246784. The ordinary and candidate three-dimensional potentials, the bar speed and all physical coefficients remain unchanged. Each integration starts again at the original initial state and ends at 4 kpc/(km/s), approximately 3.91 Gyr, with 8001 saved samples.

The bar and potential are frozen over this mathematical integration. We are testing an orbit representation of a stationary model, not reconstructing the Milky Way's actual history over four billion years. No stellar velocity outcomes or holdouts are newly examined.

The known rotating-frame Hamiltonian equations remain

    J = |p|^2/2 + Phi(x) - Omega L_z,
    dx/dt = p - Omega cross x,
    dp/dt = -grad(Phi) - Omega cross p.

These are standard mechanics, not original equations of companion physics. The extra potential is the previously frozen empirical candidate, not a newly derived deposition process. J is the Jacobi integral, the conserved Hamiltonian in this rotating frame, not inertial-frame energy by itself.

## The strict path check fails

The predeclared gates compare consecutive DOP853 runs at relative tolerances 2e-9, 2e-11 and 2e-13. Required full-path differences are below 1e-4 kpc and 0.01 km/s, alongside a Jacobi-drift and field-resolution check. Final consecutive-run results are:

| Check | Ordinary matter | Candidate extra potential |
|---|---:|---:|
| Maximum position difference over full saved path | 1.945625 kpc | 0.040807 kpc |
| Maximum velocity difference over full saved path | 161.087451 km/s | 1.887244 km/s |
| Jacobi drift divided by 220^2 | 1.003e-9 | 1.224e-10 |
| Difference from old path over first unit: position | 5.493e-5 kpc | 1.649e-6 kpc |
| Difference from old path over first unit: velocity | 0.001297 km/s | 0.0000455 km/s |
| Strict full-path gate | Failed | Failed |

The candidate extra-force refinement difference along the sampled path is at most 0.11391%, below the unchanged 1% force criterion. The old first-unit trajectory checks pass in both models. Small Jacobi drift does not establish the full trajectory: a conserved quantity constrains motion without uniquely fixing a path. No diagnosis of physical chaos, resonance or a faulty force law is established by these tolerance comparisons alone.

The main driver therefore correctly records no physically interpreted occupation windows. Raw trajectory descriptions remain explicitly ungated diagnostics.

## Separate population-sensitivity follow-up

After observing the path failure, a supplementary protocol was written. The 2e-11 runs were reproduced because the original driver retained only the tightest trajectories. Their maximum differences from the 2e-13 runs reproduce the primary results within the declared numerical tolerance. Both outputs are now saved.

For each duration 1, 2 and 4, use exactly 2000 noninitial equally spaced samples, with 1000 in each chronological half. Also compare the final interval (3,4] at the original sampling cadence. The fixed Gaussian-kernel MMD and old histogram bins are inherited from the preceding audit. No weights or bandwidths are fitted.

MMD is established distribution-comparison mathematics, not a gravity formula; see [Gretton et al. (2012)](https://www.jmlr.org/papers/v13/gretton12a.html). Here it is a dimensionless distance between weighted empirical distributions, including self-pairs, not a percentage, significance level or certified error relative to an exact solution. Nearby time samples are correlated. Two tolerances cannot establish exact-solution convergence by themselves.

At the illustrative 0.5-kpc / 50-km/s resolution:

| Potential | Interval, time units | Early/late MMD, tight solver | MMD between solver occupations |
|---|---|---:|---:|
| Ordinary matter | (0,1] | 0.191240 | 0.00000560 |
| Ordinary matter | (0,2] | 0.141343 | 0.00004493 |
| Ordinary matter | (0,4] | 0.140679 | 0.027077 |
| Ordinary matter | (3,4] | 0.211307 | 0.083526 |
| Candidate extra potential | (0,1] | 0.398119 | 0.0000000785 |
| Candidate extra potential | (0,2] | 0.211252 | 0.0000001647 |
| Candidate extra potential | (0,4] | 0.124649 | 0.000832 |
| Candidate extra potential | (3,4] | 0.375803 | 0.003278 |

At the coarser 2-kpc / 100-km/s resolution, the candidate's early/late MMD likewise decreases from 0.555384 to 0.276314 to 0.138076 as the averaging duration increases. The four-unit discrepancy between its two solver occupations is 0.000355. Thus the improvement with duration is not confined to one kernel scale or plausibly explained by the measured difference between these two solver runs alone.

However, the candidate's final one-unit interval still has early/late MMD 0.561855 at the coarser resolution, with solver-occupation discrepancy 0.001396. A longer averaging interval helps; merely waiting until a later short interval does not remove the finite-window difference. This is consistent with a need to sample more of the orbit's phases, but it does not identify a resonance or prove an eventual limiting distribution.

The ordinary-matter late population is less numerically robust in this experiment. Between solver occupations, four-unit histogram TV is 0.0895, increasing to 0.233 in the last unit; the corresponding candidate values are 0.001 and 0.003. No observational error model is attached to these numbers, so none becomes a rejection threshold or evidence favoring companion gravity.

## Verification and artifacts

- `run.py` executes the primary paired extension and writes `results.json`, preserving both failed full-path gates.
- `population_resolution.py` executes the post-failure follow-up and writes `population-resolution.json`.
- Both protocols distinguish predeclared choices from the supplementary experiment. Both calculations completed without changing physical coefficients or the original pass criteria.
- Existing input hashes, old-path checks, Jacobi drift, full-path refinement differences, candidate force refinement and source integrity are recorded. The follow-up reproduces the original full-path differences before comparing occupations. Exact finite-sample kernel matrices use the previously verified blocked implementation and its symmetry/positive-semidefiniteness checks.
- Cached trajectories, including the reproduced coarser runs, remain under `research_work/data-cache/orbit-phase-duration/`; hashes are retained in the tracked result files.

## Consequence for the full goal

There is now evidence that longer averaging improves this difficult candidate orbit's population representation. There is also a demonstrated distinction between strict individual-path accuracy and accuracy of the quantities needed for population modeling. Neither finding supplies a completed stellar likelihood or a strong case for the theory.

The next observational step should validate population predictions against solver and force resolution, integration duration and phase coverage across a broader paired set. The relevant tolerances ultimately need to be tied to the selected stellar-motion predictions and their errors. We must retain challenging trajectories and use equivalent procedures for both potentials; one favorable orbit cannot replace the bulge/disk comparison.

The new supplement does not authorize silently promoting the failed primary trajectories into a certified stationary library. No new observational fit or holdout score is reported. Photon production/transport, source and detector timing, capture/support, lensing and the deferred total energy supply remain open requirements of the unified theory. The original goal remains active.
