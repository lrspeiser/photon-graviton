# Extending the environmental map: coverage, uncertainty and identifiability

2026-09-10. This advances stage 2 of the joint observational goal. It is a measured-position geometry and uncertainty analysis, not a redshift fit, certified void map or withheld validation. The preceding goal turn was progress: it established a nearby catalog's actual coverage limit. Here we extend beyond that limit using a different source already on disk.

**Main result:** 774 selected distance-indicator tracers reach 131.8 Mpc, but only 16 lie beyond 50 Mpc. More distant apparent gaps therefore cannot be interpreted as physical voids without modeling survey selection. Distance uncertainty also changes the environmental path feature, which must enter the eventual comparison jointly with distance uncertainty.

## Which coordinates were used

The archived Cosmicflows-4 individual table contains 55,877 galaxies. We extract the separate method columns, not its combined distance, velocity, inferred peculiar velocity or reconstructed density. The main prototype selects one available positive-error measurement per galaxy, preferring maser, TRGB, Cepheid then SBF, and removes the 164 exact target IDs from the tracer sample. Other group members remain, so neither sample nor calibrations are independent in the statistical sense.

The source registers these method columns to a common scale. Calling them independent distance indicators means they are not obtained by directly inverting the target redshift; it does not mean independent calibrations, an unmodified raw measurement, or a purely geometric sample. Published calibration assumptions remain stipulated pending a more detailed source audit. [Cosmicflows-4](https://arxiv.org/abs/2209.11238)

| Selected method | Tracers |
|---|---:|
| TRGB | 445 |
| SBF | 276 |
| Cepheid | 47 |
| Maser | 6 |
| Total | 774 |

Other individual-method counts are retained in results.json: TF 12,222, FP 42,223, SNIa 1,004 and SNII 94, alongside the above methods before per-galaxy preference and target removal. Counts overlap where galaxies have multiple methods. The dense TF/FP catalogs are not declared forbidden; they require a separate audit of observable corrections, calibration and selection before supplying our coordinates.

## Radial sampling and positional uncertainty

| Distance shell | Selected tracers |
|---|---:|
| 0–11 Mpc | 432 |
| 11–25 Mpc | 256 |
| 25–50 Mpc | 70 |
| 50–100 Mpc | 13 |
| 100–200 Mpc | 3 |

This is observational sampling density, not physical galaxy or matter density. Sparse sampling at larger distances can produce a rising apparent no-neighbor fraction even if the physical environment has no such trend. Removing exact target IDs prevents every target counting as its own environmental tracer but also changes selection and is not a substitute for a completeness model.

**Known distance-modulus mathematics applied to stipulated calibrations:**

\[
D=10^{(\mu-25)/5}\ {\rm Mpc},\qquad
\sigma_{D,\rm linear}=\frac{\ln10}{5}D\sigma_\mu.
\]

The selected tracer's median radial error scale is 0.445 Mpc; its 16th–84th population percentiles are 0.086–1.562 Mpc. Among 86 selected tracers beyond 25 Mpc, **84 have quoted linear radial errors above 1 Mpc**. These are individual quoted errors, not a full covariance or survey-wide systematic uncertainty. Exact asymmetric bounds from mu +/- sigma accompany every selected tracer.

## What the path prototype can distinguish

**Illustrative catalog-proximity feature, not a physical postulate:** V is path length with no selected tracer within radius r. Exact ray/sphere intersections and interval unions compute it at r=1,3,5 Mpc. There is no dark-mass weighting, reconstructed-velocity correction or fit to redshift. All radii remain visible rather than choosing whichever later fits best.

We propagate 64 fixed-seed realizations of independent Gaussian modulus errors in both tracers and targets. Real shared calibration and association covariance is not supplied; these are conditional sensitivity distributions, not posterior credible intervals. A separate common +0.1-mag shift is an illustrative calibration scenario.

| Radius | Median no-neighbor fraction | Median sampled V interval width, 16th–84th | Median projected interval width | Small/large normalized design singular value |
|---|---:|---:|---:|---:|
| 1 Mpc | 85.97% | 4.31 Mpc | 0.79 Mpc | 0.0471 |
| 3 Mpc | 54.74% | 4.74 Mpc | 2.25 Mpc | 0.1297 |
| 5 Mpc | 29.78% | 4.86 Mpc | 3.04 Mpc | 0.1977 |

The common calibration scenario changes median absolute V by 1.48, 1.27 and 1.10 Mpc respectively. These are feature sensitivities, not measured shifts in the universe.

**Known linear algebra:** for the empirical ansatz ln(1+z_transfer)=aD+bV, form columns D and V normalized separately to unit Euclidean norm. Equal row weights are used for this geometry diagnostic. A small singular-value ratio means the two columns are close to redundant. At r=1 Mpc, only 9.40% of the norm of V remains after projecting out its component proportional to D; at r=3 and 5 Mpc the fractions are 25.52% and 38.06%. A larger radius changes the proxy definition and cannot be selected merely to make this number better.

The projected uncertainty column uses V_draw-beta*D_draw, with beta=(D dot V)/(D dot D) held at its central value. It distinguishes common distance variation from changes that can separate a distance coefficient from an environment coefficient. This supplementary diagnostic was added after inspecting the original geometry results; it uses no redshift outcomes and is labeled as such in the protocol. The 64 draws provide finite sensitivity estimates, not precision claims about percentile convergence.

This confirms that the catalog proxy has some geometric variation beyond distance, but that alone does not prove the variation traces physical voids. Its distance trend, sampling bias and measurement error can mimic the proposed effect. No a or b was estimated from redshift in this audit.

## Verification and next action

An independent 20,000-point midpoint integration on nine directions checks all three radii; maximum discrepancy is 0.00763 Mpc, passing the predeclared 0.02-Mpc numerical tolerance. Bounds 0<=V<=D, target-ID exclusion, modulus reconstruction and the rank-deficient V=D control pass. Source and protocol hashes, all 774 selected tracers, and all 492 path/radius rows are retained. These controls verify calculations, not physical void accuracy.

The optimistic 11-Mpc radial boundary is now superseded by measured tracer positions extending across the target range; the coverage problem has become quantified sparsity and selection, not solved completeness. Next audit denser TF/FP inputs and their observer-level corrections rather than treating missing galaxies as empty space. A physically defined W, completeness and correlated position uncertainty are still required before a valid distance-versus-void likelihood. Meanwhile the surviving shared-messenger candidate still needs coupled clock, spectral-shape and brightness predictions; map development must not displace those goals.

No new target redshift was queried, no conversion or gravity parameter changed, and no fresh holdout was evaluated. The larger raw catalog was historically exposed; these extractions do not make it fresh. All four goal stages remain incomplete.

Artifacts: [protocol](protocol.md), [calculation](run.py), [results](results.json), [selected distance indicators](selected-tracers.csv), [path sensitivities](path-sensitivity.csv).
