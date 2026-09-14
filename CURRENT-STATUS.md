# Current status: photon–companion research

**Authoritative summary.** Updated 13 September 2026. Reviewed baseline: `main` at `3884b4f` (Codex handover, 13 September 17:23 PDT). Added since then:
- the capture-to-orbit test and its consistency revision after review
- the self-illumination pilot
- repository consolidation

Dated history is in [CHANGELOG.md](CHANGELOG.md); standing rules and the foundational index are in [research_plan/START-HERE.md](research_plan/START-HERE.md). All samples below were exposed during development. No genuinely untouched holdout has been opened.

## The programme in one paragraph

We study a hypothetical nonexpanding universe with our observed data. Photons transfer energy to traveling "companions". Captured companions form reservoirs whose stored energy gravitates. Exact one-third retention, η(X)=X^(1/3)/(1+X^(1/3)), is the empirical reference; its origin is under investigation. Dark matter, expansion and a Big Bang are not premises. Published distances are adopted as scenario inputs, and energy accounting is mandatory. The goal is a common mechanism predicting redshift, timing, brightness, rotation and lensing with shared parameters. Nothing in the repository achieves that yet.

## Active candidates

| Candidate | Definition | Standing |
|---|---|---|
| **Original exact-third reference** | Paper equations 5–8: attenuated isotropic companion bath, capture opacity k0/[1+(r/a)^2]^2, a=2.77 R_d, fitted C0, k0, qa | Reference. SPARC test RMSE 23.59 km/s; Milky Way fiducials 6.76 / 10.41 km/s |
| **MOND-guided fixed-inventory mixture** | Reference inventory redistributed toward the simple-MOND acceleration excess; f=0.923557 and a0=8.563e-11 m/s^2 frozen | SPARC test RMSE 16.24 (simple MOND 16.40; MOND better in log error and train/validation). Both Milky Way fiducials worsen (9.26 / 11.87). 21 galaxies need inventory caps. The shape is borrowed from MOND, not derived |
| **Coma Plummer support candidate** | n=5 polytrope, P=Kρ^(6/5), fitted to six weak-shear bins | χ² 3.727 versus NFW 3.855 (total gravity only). A universal K implies a∝M^2 for the isolated profile (regression-checked). A baryon-coupled solve is still needed |
| **Capture-to-orbit RB-1** ([report](research_work/results/capture-to-orbit/report.md)), archived failed candidate | Receiver-assisted s-wave threshold production of one fixed-mass species on ordinary baryons, orbit-averaged in a fixed potential; nothing fitted | **Not promoted, revised after review.** SPARC test RMSE 45.40; Milky Way fiducials 24.98 / 7.76. Supply multiplier 4.1×10^9 (spectrum extending below threshold) or 8.2×10^9 (threshold-cut). Receivers would lose a median 34× their angular momentum even at ideal efficiency. The tested prescription fails at the reference inventory; other capture interactions are not excluded |
| **Self-illumination pilot** ([report](research_work/results/self-illumination/report.md)), diagnostic | RB-1 reaction under companion fields derived from declared rotating emitters; no co-rotation factor | **Does not relieve the RB-1 debt.** Receiver drag falls to κ=0.19–0.51 in a flat-rotation disk, and to zero inside a rigid ring. Emitters pay instead. With every companion absorbed internally, the combined baryonic loss is a median 19× the baryons' angular momentum at ideal efficiency (above 1 in 146 of 149). Threshold-cut efficiency is still halved |

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
| Capture-to-orbit RB-1 (archived) | 49.88 / 55.56 / 45.40 | 24.98 / 7.76 |

## Unresolved failures and open requirements

1. **Formation.** No local interaction yet produces the reservoir that the gravity fits need. The tested RB-1 prescription fails on supply and on receiver back-reaction. It was evaluated with one incident field per spectrum control, per-site anisotropic drag and a closed energy ledger. Self-illumination by rotating emitters does not change that: momentum conservation leaves the combined baryonic loss at a median 19× at ideal efficiency. Other capture interactions are not excluded. RB-1's 0.5 km/s quadrature gate still fails in 18 of 167 systems; passing verification jobs does not mean its full numerical model has converged there.
2. **Supply.** C0 is fitted storage, not a demonstrated radiation budget. The recovered 26-galaxy comparison falls short by a median factor of about 5,400 for 10 Gyr at present luminosity. Age is free, but the required histories must be stated.
3. **Distribution.** The best rotation shape (MOND-guided) is borrowed. It worsens both Milky Way fiducials, and inventory is short in 21 SPARC galaxies (NGC3741 needs 39× its model inventory).
4. **Outer profile.** The MOND-guided construction completes its exterior arbitrarily beyond the last measured radius. RB-1 shows a forward model can instead return a convergent (r^-4) profile on its own radial domain, independent of which radii are sampled.
5. **Joint lensing.** No single response yet predicts motions and lensing together. Only two of six lens systems improve under the transferred profile.
6. **Clusters.** No matched gas-plus-star and fixed-inventory cluster transfer exists. The Coma fits cover total gravity only.
7. **Redshift and timing.** Conversion microphysics is missing. Whole-photon mixing does not redshift survivors, supernova event stretching is unexplained, and the 74.62 versus 70.48 tension stands.
8. **Stability.** Collective stability of any reservoir is untested.

## Next experiment: a decision for the project owner

RB-1 and the self-illumination pilot pin down what a viable capture channel would have to do:

1. **Momentum.** Retained products must be born without drawing their motion from ordinary matter. Products born co-moving with rotating receivers or emitters make the baryons fund the reservoir's angular momentum, or the escaping flux's. The median reference inventory is 12.5× the baryonic mass.
2. **Energy.** A non-rotating receiver population avoids a net angular-momentum debt but not the drag: at ideal efficiency it drains about 2κ·M_ret/M_rec ≈ 40× the receivers' own kinetic energy (a first-order estimate, not yet computed per galaxy).
3. **Efficiency.** The retained fraction must be far above the (v_esc/c)^3 ≈ 10^-9–10^-10 of a threshold reaction with a smooth spectrum.

Candidates that could meet these, each needing its own declared protocol:
- a collective or field-mediated capture in which the reservoir, not ordinary matter, is the receiver, with its own drag and contraction accounting;
- a two-step absorb-and-emit process with an internal Q-value, which addresses item 3 but must still meet items 1 and 2.

The cheap next check is item 2's per-galaxy kinetic-energy budget. The baryon-coupled Coma polytrope with one universal K stays queued. Keep the original and MOND-guided branches side by side, and freeze any new rule before touching an unexposed sample.

## Reproduce

```sh
python research_work/run_checks.py                          # 49 jobs, fresh output directory
python research_work/results/capture-to-orbit/runner.py     # RB-1 populations and predictions, about 7 min on 8 workers
python research_work/results/capture-to-orbit/revision.py   # RB-1 consistency revision, about 20 s
python research_work/results/self-illumination/pilot.py     # self-illumination pilot, about 4 min
```

The three latest canonical diagnostics (`mond-inventory.py`, `mond-cross-scale.py` and `coma-inverse.py` in `research_work/results/companion-extensions/`) regenerate into a fresh directory. Each compares its numbers with the archived result and overwrites the archive only with `--canonical`. The comparison is exact, with one exception: the Coma flexible-mixture fits have non-unique weights and stop at slightly different points under different BLAS threading. Their χ² agrees to about 10^-5 and their derived masses to about 0.1%, so they use documented tolerances; the Coma LP extrema and single-Plummer fit are compared exactly. The working paper is v1.5 (`output/pdf/theory-basis.pdf`). Later results are in the supplement `papers/cumulative-time-companions/cross-scale-performance.md`.
