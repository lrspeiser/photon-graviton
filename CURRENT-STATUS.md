# Current status: photon–companion research

**Authoritative summary.** Updated 14 September 2026. Reviewed baseline: `main` at `3884b4f` (Codex handover, 13 September 17:23 PDT). Added since then:
- the capture-to-orbit test and its consistency revision after review
- the self-illumination pilot
- repository consolidation
- the three branches the project owner opened: the evolving propagation field PF-1, radiation-polarized gravity RPG-1 and the driven collective reservoir CR-1
- a [model contract](research_plan/model-contract.md) proposing one field that connects them, awaiting the owner's decision D1

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
| **Evolving propagation field PF-1** ([report](research_work/results/propagation-field/report.md)), diagnostic | One wave law ∂t(n∂tA) − (c²/n)∇²A = 0 for electromagnetic and gravitational waves, with a dynamical homogeneous index n(t); matter standards uncoupled by postulate | **Internally consistent, not promoted.** One driven simulation gives 1+z in carrier, envelope width and pulse spacing to 1.1×10^-4, with photon number conserved, and the radiation–index exchange closes to 10^-12. Its light sector is exactly coasting FLRW with a = n. Pantheon+ Δχ² = +35 against flat FLRW (Ω_m = 0.3) at equal freedom. Nonexpansion rests on a clock/ruler completion that no tested coupling supplies |
| **Radiation-polarized gravity RPG-1** ([report](research_work/results/radiation-polarized-gravity/report.md)), diagnostic | AQUAL field equation ∇·[μ(\|∇Φ\|/a*)∇Φ] = 4πGρ_b, μ = x/(1+x), archived a* frozen, solved on the frozen baryons; declared lensing Φ = Ψ in PF-1 geometry; link a* = ξc\|d ln n/dt\| | **Not promoted.** Rotation matches algebraic simple MOND on the same baryons, slightly worse in RMSE (SPARC 22.16 / 29.42 / 16.90; Milky Way 9.89 / 13.70). The Milky Way vertical force at 1.1 kpc worsens (RMS 42.0 against 15.6 Msun/pc² for baryons). SLACS Einstein radii are 0.35–0.48 of observed with Chabrier masses and 0.56–0.74 with Salpeter, failing the declared rule. With FLRW distances (post hoc) they are 0.48–0.89, and the remaining need matches the lenses' inner-dynamics masses. Field energy diverges logarithmically; ξ = 0.118. The declared convergence gate failed (bulge galaxies, up to 1.6 km/s) |
| **Driven collective reservoir CR-1** ([report](research_work/results/collective-reservoir/report.md)), diagnostic | One coherent spherical Gross–Pitaevskii–Poisson condensate with no self-interaction, for two constituent masses; a 10⁹ Msun counted seed in the Milky Way, fed by the conversion law's photon source collected inside r₉₉ | **Does not support a reservoir in this implementation.** In 10 Gyr the derived source adds 1.6×10³ Msun (light constituent) or 0.5 Msun (heavy), 10⁸ to 5×10¹¹ times short of the 2.1×10¹¹ Msun reference. At the required rate, the heavy constituent contracts below resolution within 0.4–1 Gyr (analytic endpoint: a 1.6 pc soliton; rotation RMSE 110 km/s). The light one grows 123× without overheating but reaches 59% of the inventory, and photon-shaped growth disperses. Zero self-interaction, the collection region and the 10 Gyr span are labeled assumptions |

## Frozen benchmarks

- **SPARC rotation.** 149 galaxies, 3,150 radii, original 89/29/31 train/validation/test split. Archived baryon model uses Υ_disk=0.5, Υ_bulge=0.7. Metrics are equal-galaxy RMSE (km/s) and log RMS.
- **Milky Way.** 38 Eilers circular-speed bins (5–25 kpc) and two ordinary-matter baselines. Fiducial R_d=2.6 kpc, luminosity factor 1, plus 18 declared sensitivities. These are Jeans-inferred speeds, not individual orbits. RPG-1 adds the 43 Bovy–Rix K_z/2πG values at |z| = 1.1 kpc (model-inferred).
- **Coma.** Six reconstructed Kubo weak-shear bins in dimensionless units: shape only, no matched gas/star model.
- **Lenses.** Six SLACS systems with stellar motions. Freely fitted companion profiles give total χ² 26.49 versus free NFW 49.71. These are fits, not a test of the one-third law.
- **Redshift.** 164 groups: cross-validated cz-residual RMS 448.38 km/s for constant conversion versus 448.09 for the linear control. Nearby fit c·α = 74.62 km/s/Mpc versus H0 = 74.92 for the expansion comparator on the same inputs. Supernova brightness prefers c·α = 70.48 (z < 0.3) or 69.76 (all z ≥ 0.1), with redshift-dependent residuals.

Current rotation scores on those benchmarks (RMSE km/s):

| Model | SPARC train / validation / test | Milky Way baseline I / II |
|---|---:|---:|
| Ordinary matter | 52.57 / 58.22 / 47.77 | 52.57 / 62.34 |
| Simple MOND, fitted a0 | 19.89 / 26.88 / 16.40 | 9.53 / 12.04 |
| Original exact-third | 29.03 / 32.50 / 23.59 | 6.76 / 10.41 |
| MOND-guided mixture | 20.47 / 27.14 / 16.24 | 9.26 / 11.87 |
| Capture-to-orbit RB-1 (archived) | 49.88 / 55.56 / 45.40 | 24.98 / 7.76 |
| RPG-1 field equation (surface-density baryons) | 22.16 / 29.42 / 16.90 | 9.89 / 13.70 |

## Unresolved failures and open requirements

1. **Formation.** No local interaction yet produces the reservoir that the gravity fits need.
   - **RB-1.** The tested prescription fails on supply and on receiver back-reaction. It was evaluated with one incident field per spectrum control, per-site anisotropic drag and a closed energy ledger. RB-1's 0.5 km/s quadrature gate still fails in 18 of 167 systems; passing verification jobs does not mean its full numerical model has converged there.
   - **Self-illumination.** Rotating emitters do not change the RB-1 verdict: momentum conservation leaves the combined baryonic loss at a median 19× at ideal efficiency.
   - **CR-1.** A single coherent condensate fed by the conversion law cannot grow.

   Other capture interactions are not excluded.
2. **Supply.** C0 is fitted storage, not a demonstrated radiation budget. The recovered 26-galaxy comparison falls short by a median factor of about 5,400 for 10 Gyr at present luminosity. Inside a reservoir's own volume the shortfall is far larger: CR-1's seed collects 0.05–160 Msun per Gyr. Age is free, but the required histories must be stated.
3. **Distribution.** The best rotation shape (MOND-guided) is borrowed. It worsens both Milky Way fiducials, and inventory is short in 21 SPARC galaxies (NGC3741 needs 39× its model inventory).
4. **Outer profile.** The MOND-guided construction completes its exterior arbitrarily beyond the last measured radius. RB-1 shows a forward model can instead return a convergent (r^-4) profile on its own radial domain, independent of which radii are sampled. RPG-1's field energy grows by M v_f²/3 per e-fold of radius, so a field response needs an outer boundary from outside its equation.
5. **Joint lensing.** No single response yet predicts motions and lensing together.
   - Only two of six lens systems improve under the transferred profile.
   - RPG-1's declared response in PF-1 geometry gives 0.35–0.74 of the SLACS Einstein radii from population stellar masses; at those radii g_N ≈ 10 a*.
   - With FLRW distances (post hoc) the lenses still need 2.1–2.7× Chabrier. That matches their inner-dynamics masses, so it is one inner-mass shortfall.
6. **Clusters.** No matched gas-plus-star and fixed-inventory cluster transfer exists. The Coma fits cover total gravity only.
7. **Redshift and timing.** Conversion microphysics is missing, and whole-photon mixing does not redshift survivors. PF-1 derives redshift and event stretching together from one wave law. It does so only with clocks and rulers that no tested matter completion provides, and its light sector is that of a coasting expanding universe. The brightness tension stands: 74.62 from galaxy groups against 69.8–70.5 from supernovae, with high-redshift residuals.
8. **Stability.** CR-1 evolved a coherent seed for 10 Gyr. Without a source it stays stationary. At the required rate the heavy constituent contracts below resolution; its analytic endpoint is a stable 1.6 pc soliton at 25% of the Kaup limit, not a black hole. Photon-shaped growth disperses. A supported reservoir, with self-interaction or a mode population, is untested.
9. **Vertical force.** Responses that repair Milky Way rotation over-predict K_z at 1.1 kpc: RPG-1 gives RMS 42.0 against 15.6 Msun/pc² for baryons alone, and the archived conservative completion went from 27.99 to 33.31.

## Next experiments

All three branches the owner opened on 13 September have run, and none is promoted. Each first implementation was narrower than its general hypothesis; section 1 of the [model contract](research_plan/model-contract.md) separates the two.

The model contract proposes connecting the branches through one clock field:
- PF-1's propagation field, completed by matter scales that co-evolve with it;
- its static response around galaxies, which is RPG-1's equation with a* tied to the redshift rate;
- its condensate, which supplies only the inner-core and cluster mass that an acceleration-keyed response cannot.

**Decision for the owner (D1).** Does "nonexpanding" allow matter's rulers to co-scale with the field? Galaxies stay fixed in the field's coordinates, but separations counted in rulers grow. Observationally that is identical to expansion with scale factor n(t).
- **If yes,** the queue is:
  1. CC-1, the clock completion with measured ratios;
  2. CG-1, clock-field gravity with derived lensing;
  3. CR-2, a supported condensate core.
- **If ruler counts must stay fixed,** the first test is whether a screened field bends light toward galaxies at all.

Standing constraints remain: retained products must not draw their motion from ordinary matter, and the retained fraction must far exceed (v_esc/c)^3. Keep the original and MOND-guided branches side by side, and freeze any new rule before touching an unexposed sample.

## Reproduce

```sh
python research_work/run_checks.py                                                # 53 jobs, fresh output directory
python research_work/results/capture-to-orbit/runner.py                           # RB-1 populations and predictions, about 7 min on 8 workers
python research_work/results/capture-to-orbit/revision.py                         # RB-1 consistency revision, about 20 s
python research_work/results/self-illumination/pilot.py                           # self-illumination pilot, about 4 min
python research_work/results/propagation-field/pf1.py                             # PF-1 T1-T3
python research_work/results/propagation-field/brightness.py                      # PF-1 T4 on exposed Pantheon+
python research_work/results/radiation-polarized-gravity/rpg1.py                  # RPG-1, about 90 s on 8 workers; exits 1 because its declared gate V1 failed
python research_work/results/radiation-polarized-gravity/convergence.py           # RPG-1 post-hoc convergence diagnostic, about 5 min
python research_work/results/radiation-polarized-gravity/geometry-diagnostic.py   # RPG-1 post-hoc lens-geometry diagnostic
python research_work/results/collective-reservoir/cr1.py                          # CR-1, fails fast and logs progress
```

The latest canonical diagnostics regenerate into a fresh directory:
- `mond-inventory.py`, `mond-cross-scale.py` and `coma-inverse.py` in `research_work/results/companion-extensions/`;
- the two PF-1 scripts;
- RPG-1's `rpg1.py`;
- CR-1's `cr1.py`.

Each compares its numbers with the archived result and overwrites the archive only with `--canonical`. The comparison is exact, with one exception: the Coma flexible-mixture fits have non-unique weights and stop at slightly different points under different BLAS threading. Their χ² agrees to about 10^-5 and their derived masses to about 0.1%, so they use documented tolerances; the Coma LP extrema and single-Plummer fit are compared exactly. The working paper is v1.5 (`output/pdf/theory-basis.pdf`). Later results are in the supplement `papers/cumulative-time-companions/cross-scale-performance.md`.
