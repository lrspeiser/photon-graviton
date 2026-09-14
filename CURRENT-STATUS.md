# Current status: photon–companion research

**Authoritative summary.** Updated 13 September 2026. Reviewed baseline: `main` at `3884b4f` (Codex handover, 13 September 17:23 PDT), plus the capture-to-orbit feasibility test and repository consolidation recorded here. Dated history is in [CHANGELOG.md](CHANGELOG.md); standing rules and the foundational index are in [research_plan/START-HERE.md](research_plan/START-HERE.md). All samples below were exposed during development. No genuinely untouched holdout has been opened.

## The programme in one paragraph

We study a hypothetical nonexpanding universe with our observed data. Photons transfer energy to traveling "companions". Captured companions form reservoirs whose stored energy gravitates. Exact one-third retention, η(X)=X^(1/3)/(1+X^(1/3)), is the empirical reference; its origin is under investigation. Dark matter, expansion and a Big Bang are not premises. Published distances are adopted as scenario inputs, and energy accounting is mandatory. The goal is a common mechanism predicting redshift, timing, brightness, rotation and lensing with shared parameters. Nothing in the repository achieves that yet.

## Active candidates

| Candidate | Definition | Standing |
|---|---|---|
| **Original exact-third reference** | Paper equations 5–8: attenuated isotropic companion bath, capture opacity k0/[1+(r/a)^2]^2, a=2.77 R_d, fitted C0, k0, qa | Reference. SPARC test RMSE 23.59 km/s; Milky Way fiducials 6.76 / 10.41 km/s |
| **MOND-guided fixed-inventory mixture** | Reference inventory redistributed toward the simple-MOND acceleration excess; f=0.923557 and a0=8.563e-11 m/s^2 frozen | SPARC test RMSE 16.24 (simple MOND 16.40; MOND better in log error and train/validation). Both Milky Way fiducials worsen (9.26 / 11.87). 21 galaxies need inventory caps. The shape is borrowed from MOND, not derived |
| **Coma Plummer support candidate** | n=5 polytrope, P=Kρ^(6/5), fitted to six weak-shear bins | χ² 3.727 versus NFW 3.855 (total gravity only). New analytic check: a universal K implies a∝M^2 for the isolated profile. A baryon-coupled solve is still needed |
| **Capture-to-orbit RB-1** (new, [report](research_work/results/capture-to-orbit/report.md)) | Receiver-assisted s-wave threshold production of one fixed-mass species on ordinary baryons, orbit-averaged in a fixed potential; nothing fitted | **Not promoted.** SPARC test RMSE 45.40, worse than the reference in every split. Milky Way fiducials 24.98 / 7.76. Supply multiplier about 4×10^9. Even at ideal efficiency, receivers would lose a median 25× their angular momentum |

## Frozen benchmarks

- **SPARC rotation.** 149 galaxies, 3,150 radii, original 89/29/31 train/validation/test split. Archived baryon model uses Υ_disk=0.5, Υ_bulge=0.7. Metrics are equal-galaxy RMSE (km/s) and log RMS.
- **Milky Way.** 38 Eilers circular-speed bins (5–25 kpc) and two ordinary-matter baselines. Fiducial R_d=2.6 kpc, luminosity factor 1, plus 18 declared sensitivities. These are Jeans-inferred speeds, not individual orbits.
- **Coma.** Six reconstructed Kubo weak-shear bins in dimensionless units: shape only, no matched gas/star model.
- **Lenses.** Six SLACS systems with stellar motions. Freely fitted companion profiles give total χ² 26.49 versus free NFW 49.71. These are fits, not a test of the one-third law.
- **Redshift.** 164 groups: cross-validated cz-residual RMS 448.38 km/s for constant conversion versus 448.09 for the linear control. Nearby fit c·α = 74.62 km/s/Mpc versus H0 = 74.92 for the expansion comparator on the same inputs. The supernova-brightness optimum is 70.48, with redshift-dependent residuals.

Current rotation scores on those benchmarks (RMSE km/s):

| Model | SPARC train / validation / test | Milky Way baseline I / II |
|---|---:|---:|
| Ordinary matter | 52.57 / 58.22 / 47.77 | 52.57 / 62.34 |
| Simple MOND, fitted a0 | 19.89 / 26.88 / 16.40 | 9.53 / 12.04 |
| Original exact-third | 29.03 / 32.50 / 23.59 | 6.76 / 10.41 |
| MOND-guided mixture | 20.47 / 27.14 / 16.24 | 9.26 / 11.87 |
| Capture-to-orbit RB-1 | 49.88 / 55.56 / 45.40 | 24.98 / 7.76 |

## Unresolved failures and open requirements

1. **Formation.** No local interaction yet produces the reservoir that the gravity fits need. RB-1 shows the obstruction is structural. Ordinary-matter receivers must supply the created mass's momentum (κ E v/c^2, κ between 4/3 and 1.74), and median inventories are 12.5× the baryons. Smooth spectra retain only (v_esc/c)^3 of the captured energy, and sharp lines react only in a ~0.1 km/s receiver-speed window.
2. **Supply.** C0 is fitted storage, not a demonstrated radiation budget. The recovered 26-galaxy comparison falls short by a median factor of about 5,400 for 10 Gyr at present luminosity. Age is free, but the required histories must be stated.
3. **Distribution.** The best rotation shape (MOND-guided) is borrowed. It worsens both Milky Way fiducials, and inventory is short in 21 SPARC galaxies (NGC3741 needs 39× its model inventory).
4. **Outer profile.** The MOND-guided construction completes its exterior arbitrarily beyond the last measured radius. RB-1 shows a forward model can instead return a convergent (r^-4) profile on its own radial domain, independent of which radii are sampled.
5. **Joint lensing.** No single response yet predicts motions and lensing together. Only two of six lens systems improve under the transferred profile.
6. **Clusters.** No matched gas-plus-star and fixed-inventory cluster transfer exists. The Coma fits cover total gravity only.
7. **Redshift and timing.** Conversion microphysics is missing. Whole-photon mixing does not redshift survivors, supernova event stretching is unexplained, and the 74.62 versus 70.48 tension stands.
8. **Stability.** Collective stability of any reservoir is untested.

## Next experiment

The binding constraint is now momentum, not gravity or shape. Proposed order:

1. **Self-illumination momentum budget.** Companions made from a rotating disk's own starlight carry that disk's aberration. Compute their mean momentum per energy at each receiver. A co-rotating bath would reduce the drag that kills RB-1; an external isotropic bath cannot.
2. **Receiver back-reaction.** Evolve disk angular momentum under the drag. This bounds the reservoir mass any ordinary-matter channel can create before the disk contracts measurably.
3. **Two-step absorb-and-emit** with a fixed internal Q-value. It removes the (v/c)^3 efficiency penalty but must carry item 1's momentum accounting.

The baryon-coupled Coma polytrope with one universal K stays queued as the support test. Its isolated limit, a∝M^2, is now a regression check. Keep the original and MOND-guided branches side by side, and freeze any new rule before touching an unexposed sample.

## Reproduce

```sh
python research_work/run_checks.py                          # 49 jobs, fresh output directory
python research_work/results/capture-to-orbit/runner.py     # full frozen comparison, about 7 min on 8 workers
```

The three latest canonical diagnostics (`mond-inventory.py`, `mond-cross-scale.py` and `coma-inverse.py` in `research_work/results/companion-extensions/`) now regenerate into a fresh directory. Each compares its numbers with the archived result and overwrites the archive only with `--canonical`. The comparison is exact, with one exception: the Coma flexible-mixture fits have non-unique weights and stop at slightly different points under different BLAS threading. Their χ² agrees to about 10^-5 and their derived masses to about 0.1%, so they use documented tolerances. The Coma LP extrema and single-Plummer fit are compared exactly. The working paper is v1.5 (`output/pdf/theory-basis.pdf`). Later results are in the supplement `papers/cumulative-time-companions/cross-scale-performance.md`.
