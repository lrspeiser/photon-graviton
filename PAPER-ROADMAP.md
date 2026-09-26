# What it would take to publish this

**22 September 2026.** An honest assessment of the gap between where the work is and
a paper that survives peer review, with the tasks ordered by what blocks what.

---

## Status — 26 September 2026 (hot-companion rows updated during round 25)

| Task | Status | Result |
|---|---|---|
| **RULE** No MOND, no Newton, no dark matter; check every formula | **In force** | Recorded in [RULES.md](RULES.md). The automatic check lives in `research_work/tools/formula_guard.py`. It correctly flags MOND itself. The MOND law used until now, and everything built on it (T0.2's dynamics equation, the AeST home, the η values), is retired as our answer. |
| **NEW** Hot-companion gravity | **Works on galaxies, clusters and the Solar System** | Heat feeds the companion without cancelling, and strong fields hold it back. **Galaxies:** 15.93 km/s, against MOND's 16.13, better on train, validation and test. **Clusters:** error 0.223, against MOND's 1.062; held-out cross-check 0.234. MOND refit on clusters only reaches 0.338, and needs a constant 9.7 times its galaxy value. **Solar System:** zero anomalous pull. Three universal constants. **Round 2 (rows below):** MOND is derived as its cold limit, it has a field equation and a momentum-conserving action, and the KiDS, lens and slip problems are resolved or made testable. **Open:** merging clusters. Details in [research_work/results/hot-companion/](research_work/results/hot-companion/README.md). |
| **NEW-2** MOND from our law | **Done** | Cold matter emits in step, so companion flows add like Newton's field (identity to 1.8 × 10⁻¹⁵). Heat Doppler-scrambles them into a plain sum, shown numerically. Energy balance gives the deep law with **a₀ = 2ℓ/u** (ℓ = 2.9 × 10⁻⁵ W/kg). The resulting interpolating function, `1 + e^(−y/λ)/√y`, is new: 0.031 dex from "simple". MOND is an output, not an input. Test: the Sun's extra mass loss, 1.4 × 10⁻¹⁴ per year. `code/derive_mond.py`. |
| **T0.2′** Field equation and action for the new law | **Done (non-relativistic)** | Conservative potential form; the cold limit is QUMOND with our ν. The heat term as first written broke momentum (25%); an action adds a reaction force on hot matter and restores it (1.3% at the largest box, shrinking with box size). `code/field_equation.py`. Relativistic form open. |
| **KiDS** Early/late lensing offset | **Reproduced** | Hot gas haloes of about one stellar mass (0.6–1.0 keV) give +0.20 to +0.22 dex beyond 300 kpc, against the observed ≥ 0.2. MOND with the same haloes gives 0.14. Prediction: the offset grows with radius. `code/kids_haloes.py`. |
| **T3.5** Stellar-mass convention | **Done** | Standard distances plus published Chabrier masses remove 0.13–0.16 dex. Our law then needs stars 1.2–1.6× Salpeter in the six SLACS lenses; dark-matter models need about 1.0. Spectroscopic IMFs are the test. `code/lenses_t35.py`. |
| **Slip** Re-derived against the new law | **Not needed** | Lensing and resolved kinematics agree to +0.040 ± 0.026 dex. The earlier η ≠ 1 was an artefact of measuring against the MOND law. Cluster X-ray masses sit 22% below our prediction (b = 0.22), in the direction of the known hydrostatic bias. |
| **ROUND 3** Collision rule + flow direction | **Done; law locked in** | Gas collides, so its companion stays in step (Dicke narrowing): only stars and galaxies carry heat. The companion pulls along its net energy flow. Refit: u = 197 km/s; galaxies 15.85 km/s (MOND 16.13); clusters 0.227 (held-out 0.244; law-only stellar speeds 0.329). KiDS early/late from early-type stars: 0.17–0.27 dex (obs 0.17/0.27). SLACS: stars 1.05–1.35 × Salpeter; no slip. Momentum is carried by the companion; a matter-only action would push cluster stars outward by 30–70% of gravity, so it is ruled out. Scripts: `code/dicke_toy.py`, `code/run_v3.py`, `code/kids_v3.py`. |
| **T2.1** Bullet Cluster (round 6) | **Main cluster settled; subcluster lensing mass open** | Legacy Survey DR10 star count (1,652 cluster galaxies by photometric redshift): light rises to 2.5–3 Mpc with 23–28% beyond 1.5 Mpc; total 85–103% of the round-5 model inferred from the galaxy speed. Subcluster: the 1:3 solution's lost galaxies are not in Barrena et al.'s 71 non-core velocities (about 3σ for a compact population, 1.4σ for a wide fast one) or in the starlight within 500 kpc (0.46 of the main's central light; core alone 0.32; 1:3 would give 1.55). In our law its lensing needs its hot stars, so the 2× lensing gap is open. Next: the heat gained in the crossing, carried by the memory. Scripts: `code/bullet_light_v6.py`, `code/bullet_members_v6.py`. |
| **R8** Regression suite; three more colliding clusters (round 8) | **Done: suite in place; MACS J0025 passes; Abell 520's dark core comes from gas and heat (5 of 6 clumps pass); El Gordo points to its star masses** | **Suite** (`research_work/results/hot-companion/regression/`): every test the law has faced, graded the same way (pass within 2σ, close within 3σ), compared with a saved baseline. It reproduces every published number. Amendments tested: gradual release (L = 30,000 AU) fixes Cassini and costs nothing; a single external hold of 0.1 fixes Cassini and lifts every dwarf, but makes wide binaries 2.5× Newtonian, so the hold must depend on relative speed; g_d × 1.25–1.5 lifts the Sun to 215–218 km/s but pushes the vertical pull and the SLACS gap to 'close', so the disk's shape is the better lever. **MACS J0025.4−1222:** lensing inside 300 kpc 2.0 / 1.6 × 10¹⁴ against 2.5 (+1.0/−1.7) / 2.6 (+0.5/−1.4); peaks on the galaxies; galaxy speeds 770 against 835 ± 59 km/s. **Abell 520:** five of six clumps within the errors; the galaxy-poor 'dark core' P3 gets 3.06 × 10¹³ inside 150 kpc (measured 2.84–3.35) from gas and the surrounding galaxies' heat. **El Gordo:** with the published (Chabrier, SED) stars, lensing inside 1 Mpc 11.3 against 15.8 × 10¹⁴ and galaxy speeds 944 / 839 against 1,290 / 1,089 km/s; twice the stars recovers both. Scripts: `code/collisions_v8.py`, `regression/run_suite.py`. |
| **R9** Adopt the Cassini fix; what it opens; a plan for every open item (round 9) | **Done: release length adopted (62 pass, 8 close, 6 fail); the dwarfs' way out found, not yet derived** | **Adopted:** the companion is released over L = 0.15 pc (R = 1 − e^(−r/L)). Cassini's Q2 3.1 × 10⁻²⁶ → 4.6 × 10⁻²⁷ s⁻²; nothing else changes (full tier; new baseline `round9`); wide binaries 9% at 20,000 AU. **Opened:** the dwarfs are slowed by the Galaxy's pull inside their square root (external field effect). Without it, 6 of 10 agree (χ² 135 → 60; MOND 119). With L ≥ 0.5 pc the binaries stay in range (1.19–1.36). Not adopted until derived: a separate-roots rule would contradict the heat rule. Script: `code/release_hold_v9.py`. |
| **R10** First principles; our own distances; a locked forecast (round 10) | **Done: postulate 4 has a mechanism (locked-emitter recoil, simulated); of three combination rules only round 3's fits both SPARC and X-COP; far collisions and SLACS in static distances, now graded that way; wide-binary forecast locked** | **Derived:** a quarter-cycle-leading emitter is pulled by (q₀k/2) × the local amplitude at every distance (3D FDTD within 4%); locally ⟨F⟩ = (q₀/2)\|Ψ\|∇θ. Conditions: active medium, spectral locking, v_phase ≲ 0.05 u. **Tested:** independent waves (SPARC 19.57 km/s, X-COP 0.382) and one stream with nothing cancelling (18.97, 0.222) against round 3 (15.85, 0.227): ordered opposing flows cancel, scrambled companion adds in full, one stream at the full phase rate. **Static distances:** suite 58 / 11 / 7 (was 62 / 8 / 6); every moved grade is a far cluster with 20–29% lighter stars; 1.4× the stars gives 63 / 7 / 6. SLACS light = matter −0.012 dex, stars 1.44–1.95 × Salpeter. Forecast: γ = 1.039 (7,000 AU), 1.093 (20,000 AU). Scripts: `code/first_principles_v10.py`, `code/combination_rules_v10.py`, `code/collisions_v10.py`, `code/collisions_star_sweep_v10.py`, `code/forecast_wide_binaries_v10.py`. |
| **R11** Energy bookkeeping; the stars audited; everything in our own distances (round 11) | **Done: the ordered rule derived from energy conservation plus one stream; X-COP's stars corrected (projected → spherical) and every distance-dependent test in the static law; constants refitted and adopted (u = 162.6 km/s); suite 57 / 9 / 10** | **Derived:** with div J = ℓρ, curl J = 0 and one stream at the full phase rate, J = (ℓ/4πG)(−g_N), so the ordered intensity is \|g_N\|; "one stream of everything" carries 1.06–11.8× the power emitted inside spheres (Milky Way model 1.06–1.70×). Dials (400 bootstraps): w = 1 (SPARC), γ ≤ 0.1 (both), h = 1 (X-COP). **Audit:** X-COP's stellar profiles are projected (Ghizzardi et al. 2021, Sect. 4.1; deprojected fit reproduces Table 3 within 0–13%); MACS J0025's stars are Salpeter-based (Drory et al. 2004); distant cluster stars are younger, so no age-cap correction. **Recalibration:** a 6.547 × 10⁻¹¹, g_d 2.107 × 10⁻¹⁰ m/s², u 162.6 km/s (148–178 over 924 half-splits); SPARC 15.94 km/s, X-COP 0.221. **Far clusters** 12 / 4 / 1 of 17 checks (were 9 / 7 / 1); MACS NW peak slides to the gas at 0.5 Gyr (returns at ≤ 0.35 Gyr or 1.33× stars). **Bullet:** κ 0.715 / 0.259; sub mass 1.29 against 2.35–2.70 × 10¹⁴ (still half). **KiDS in our distances:** all / blue / disks +0.063 / +0.077 / +0.107 dex (level ∝ α; centred at 0.87 α), gap 0.234 against 0.153 (set by u). Two static geometries compared: clusters prefer D_A = D. Scripts: `code/companion_flow_v11.py`, `code/xcop_static_v11.py`, `code/kids_static_v11.py`, `code/bullet_static_v11.py`, `code/distance_variants_v11.py`, `code/macs_peak_scan_v11.py`; audit `literature/star_mass_audit_v11.md`. |
| **R12** Round 11's next steps, one at a time (round 12) | **Done: MACS J0025 dated by its shock fronts (both peaks on their galaxies); the lenses' heat measured (the early/late gap matches); the distance scale fitted jointly and adopted (α ×0.95; suite 59 / 11 / 7); a dynamical toy in which only guided streaming meets the three requirements** | **MACS J0025:** relics at the galaxies' distances; shock and separation clocks 0.08–0.36 Gyr in our distances; at 0.3 Gyr the NW peak is 72–79 kpc from its galaxies (pass). **Heat (SDSS DR17, 119,000 galaxies):** k = 2.14 / 0.55 / 2.15 / 0.30 (red / blue / bulges / disks) against 2.9 / 0.1 assumed; KiDS gaps 0.132 / 0.163 against 0.153 / 0.154; IMF ≤ 0.09 dex above 2 × 10¹¹ only. **Distance scale:** Pantheon+ ×0.955 ± 0.013 (drifting 0.967 → 0.944 with depth: the law's shape), SPARC's Hubble-flow galaxies ×0.87 (0.82–0.92); X-COP and KiDS flat; adopted ×0.95 with SPARC in the static law: a 6.298 × 10⁻¹¹, g_d 2.027 × 10⁻¹⁰ m/s², u 169.4 km/s. KiDS's level (+0.065 dex) follows the lenses' mass at fixed light; circumgalactic gas of 0.5–1 M* within 100 kpc closes it. **Toy:** free streaming f 0.83–0.87; scattering f 0.12; annihilation loses 19–24%; alignment f 0.98 but whirlpools (0.43–0.58) and 0.35–0.75 dex off J_N; guided f = 1, lossless, J_N to the grid's precision; a flux-tube argument makes it exact. Wide-binary forecast at the new constants: 1.032 / 1.076 at 7,000 / 20,000 AU. Scripts: `code/macs_timing_v12.py`, `code/lens_heat_sdss_v12.py`, `code/kids_heat_v12.py`, `code/sn_scale_v12.py`, `code/distance_scale_v12.py`, `code/kids_level_v12.py`, `code/companion_toy_v12.py`. |
| **R13** The pull and the scrambling in one experiment (round 13) | **Done: the joint experiment built and validated; with one fixed offset free motion halves the pull, for an exact reason; an offset that depends on the local wave's strength gives dark cold clouds that pull, with free motion strengthening the pull (×2), collisions suppressing the gain and rotation giving none; the data pin the heat exponent at p = 1.75–2** | **Experiment** (N = 100 emitters, one locking rule, 32 test bodies at 6λ feeling only the local wave; validation (k0/2)|E| to five digits): round 10's offset gives a chorus (4.39 × independent) whose pull free motion halves (0.36–0.49 of rest at q = 1–10) and collisions protect (0.64–0.96). **Exact reason:** F = −(P/c) k̂, and locked identical emitters share one phase to their local wave, so test bodies are pulled iff the cloud outshines independent emitters (16 offsets, no exception). The absorbing offset gives the heat pattern in the output (0.197 → 0.79 free, 0.40 colliding; 1/35 darker cloud 0.028 → 0.45) but pushes; free amplitudes are brighter (5–15×). **Strength-dependent offset** (−π/2 in weak waves, +π/2 in strong): dark and pulled for E_s = 0.1–0.3; densest cloud (1/87): free ×1.32 → 1.99 for q = 1/32 → 1/4, collisions ×1.07 → 1.39, rotation ×1.01 → 0.90; test-body wave ∝ √(radiated) (the law's √(|g_N| + S)); first doubling ×1.9 in the extra pull; beyond q ≈ 1/4 the test bodies lose step. **Data:** k ∝ σ^p, X-COP flat in p, KiDS gap p = 1.75–2 (doubling σ: ×1.8–2.0 extra pull); p = 1.75 removes KiDS's level (+0.011 dex). **Energy:** the heat's extra power drawn from motion would last u/a = 85 Myr, so it comes from matter's supply. The law and suite unchanged (59 / 11 / 7). Scripts: `code/coherent_force_v13.py`, `code/joint_checks_v13.py`, `code/strength_offset_v13.py`, `code/heat_exponent_v13.py`. |
| **R14** Round 13's next steps, one at a time (round 14) | **Done: the heat exponent as a switch, p = 1.75 tested on the full suite (60 / 10 / 7; not adopted); the suite's u refit fixed; round 13's timing rule derived from a power balance, with no memory; the singers' hold fits SPARC like the release factor; darker cold states widen the heat range bodies can follow** | **Exponent:** p = 1.75 (a 6.181e-11, u 132.2): KiDS all 0.065 → 0.012, red 0.078 → 0.018 (pass); Mistele ETG 2.6 → 5.2 (fail), SLACS −0.028 → −0.054, SPARC bulges 29.4 → 30.5. **Suite fix:** u refit on the graded static, deprojected X-COP (was the round-3 sample: 192 km/s for the adopted law); the control reproduces the baseline exactly. **Singer that absorbs** (s = a + i c0 E): dark and pulled (1/47, +0.075), but motion's gain is un-jamming (×3.6–4.0, memory: 2 of 4 stay bright); resonant, no release. **Power balance** (|E| sin(lead) = P_s − c_L|E|²): all twelve clouds dark and pulled; densest, free ×1.30–1.81 (q ≤ 1/8), colliding ×0.98–1.08, rotating none; no memory. **Hold on SPARC:** exp 15.87, absorber 15.86, balance 16.08 km/s; pushing or returning holds 16.19–16.78. **Release depth:** q* ≈ 4 × darkness; bodies keep step to q ≈ 1; after onset the release grows σ^0.6–0.8 (the law: σ²). Scripts: `code/singer_v14.py`, `code/release_shape_v14.py`, `code/release_depth_v14.py`, `regression/candidates/heat_p175.json`. |
| **R25** Flow frustration tested as registered (round 25) | **Not supported: the misses do not grow with the share of the companion's energy that is not flowing (pooled p = 0.44; passes median chi 0.092, misses 0.006; all three registered conditions fail), because the law already spends the hot glow's non-flowing energy in the pull's size; a companion that keeps turning (τ = u/a, 14.8 kpc) is excluded by lensing and flat rotation curves. Leads found after looking: the misses are one-sided (36 of 43 low), 8–16% even with nothing frustrated, and the dwarfs are ordered by the Galaxy's pull over their own (6 of 6 with η ≥ 0.78 miss, 4 of 4 with η ≤ 0.31 pass), the case `no_hold_r12` addresses** | **Registration:** `round25-flow-frustration.md` (commit e74f3f8), before any value. **Test** (`code/frustration_v25.py`, 3,435 points): SPARC ρ −0.03; X-COP +0.28, −0.56 at fixed radius; Milky Way +0.09; dwarfs −0.37; collisions +0.10 (all p > 0.01); KiDS and Mistele chi = 0 exactly (point lenses) with the law 16% and 13% low. **Exact:** Q = Q_hot + Q_dir; flicker m² + f² = 1 for two opposed coherent streams. **Turning** (`code/frustration_transport_v25.py`): reach 211 kpc in 13 Gyr (straight: 2,253), energy ×9 at 60–150 kpc (random walk agrees). **Exploratory** (`code/frustration_explore_v25.py`): sign test p = 9 × 10⁻⁶; dwarfs ρ(η) = +0.87. Next: a registered no-hold test on new dwarfs (M31's); weigh the lenses' circumgalactic gas. |
| **R24** The Casimir-EFT candidate reviewed and run on the data (round 24) | **The note's numbers all reproduce (bag model: 95% of a, strictly of a/u; its strong-field lead is the same scale again); its local field equation fails test A (47 / 13 / 17: the clusters and collisions), with no rescue from the law's memory or a refit (which also drops the bag match to 58%); one added term, the hot glow's non-flowing brightness S_ex = S − \|g_hot\|, gives 59 / 11 / 7, the law's tally, and 58 / 11 / 8 with no memory (registered, not adopted)** | **Suite (77 graded):** law 59/11/7; note's equation 47/13/17 (X-COP 0.486, trend 0.715; collisions 12/6/9; Bullet smaller half's peak 845 kpc; stack β 0.180); refit a 5.98 × 10⁻¹¹, u 98.6: 50/12/15 (ellipticals' lensing z 11.3); no hold 41/12/24 (planets 3 × 10⁻³); with memory 50/12/15; **repaired 59/11/7** (with memory), **58/11/8** (none; MACS J0025 NW peak 259 kpc); brightness only 59/10/8, 51/11/15. **Exact (curl) solution:** note's equation < 2% change (Bullet); repair, damped, settles: 58/11/8 on the full suite (Bullet peaks 27/50 kpc, κ 0.695/0.222, stack β 0.021; only MACS J0025 NW peak lost, 213 kpc). **Test B bound:** saturating heat factor safe by 5.6 × 10¹⁰; γv²/u² marginal (52 Mpc); (γv)²/u² fails by 2 × 10¹¹. Scripts: `code/eft_field_v24.py`, `code/eft_checks_v24.py`, `code/eft_suite_table_v24.py`. |
| **R23** The two clusters' companions interfere (round 23) | **The owner's question: interference adds no energy, but makes the combined companion flicker at each crossing star, with depth m = 2√(I_other I_own)/(I_other + I_own), far larger than the share p when the other companion is faint; driving the heat of crossing by m supplies about a third of the Bullet's missing lensing, keeps the lensing on the galaxies, leaves MACS J0025 and El Gordo within their errors, and cuts the fast glow's extra power from 3.3× to about 2.2× (registered, not adopted)** | **Bullet (22.5 kpc), smaller / main:** share at u 1.586 / 3.101; flicker at u 1.695 / 3.194; flicker 300 km/s 1.874 / 3.203; flicker 600 km/s 1.995 / 3.190; ×1.5 power 2.216 / 3.300; ×2 2.401 / 3.404; peaks 5–11 kpc. **MACS J0025 SE/NW:** 2.24/2.03 → 2.23/2.00 (u), 2.23/2.02 (600), peaks 13, 18 → 4, 9 kpc. **El Gordo 500/1,000 kpc:** 8.86/22.1 → 8.68/21.6, 8.61/21.7. **Averaging:** two equal interfering waves average 1.27 of one wave's height, against 1.41 without interference. Scripts: `code/beat_heat_v23.py`, `code/beat_heat_far_v23.py`. |
| **R22** The fast glow from collisions, taken as far as the model allows (round 22) | **Not a matter of frame: a glow left partly behind in the other system's frame adds no mass to the Bullet's smaller half and drags its lensing onto the gas, so the fast glow must spread evenly around its stars; its power (3.3× the heat rule, 5.8× the collision's kinetic energy) must come from matter's internal store, a quarter of what a cluster galaxy's settled heat already costs over 10 Gyr; in the medium it needs a two-regime dispersion (slow crests at the companion's frequency, group velocity rising 3.5-fold up to the crest-crossing frequency, 350–2,300 ω₀); the even glow passes the stack's β, the glow left behind fails it** | **Frame test (22.5 kpc; μ = share of the star's motion kept):** smaller half 1.49–1.68 for μ = 0–0.9 (round 16 1.60, even glow 2.54); its peak 93–189 kpc toward its gas and gas lensing 0.15–0.17 for μ = 0.8–0.9. **Energy (smaller half, inner stars):** 1.3 × 10¹³ J/kg (1.8× kinetic at 3,900 km/s); ×3.3 = 4.3 × 10¹³ (4.7 × 10⁻⁴ of rest energy; 0.24 of 10 Gyr of settled cluster heat); peak 5.9× (19×) a Sun-like star's light per kg. **Medium:** needed slope d ln v_g/d ln ω 0.16–0.22; a power law with v_g/v_phase = 20 gives ×263–1,563; two-regime curvature 0.12–0.16 u/k₀. **Stack (Bullet estimate):** Δβ −0.04 (u), −0.03 (even glow as proposed), −0.02 (energy booked), +0.17/+0.22/+0.36 (μ = 0.9/0.85/0.8). Scripts: `code/crossing_frame_v22.py`, `code/crossing_energy_v22.py`, `code/fast_glow_medium_v22.py`, `code/crossing_offsets_v22.py`. |
| **R21** Three proposals from the owner's feedback checked (round 21) | **A faster glow from collisions (v_h) reproduced: both Bullet halves in range for v_h ≈ 570–820 km/s, but only if the collision puts about 3.3× the power into it (with the energy booked the smaller half reaches 1.76); the heat of crossing added to MACS J0025 and El Gordo; the distance law's path factor √(1 + η z/(1 + z)) with η = 1/2 fits the supernovae as well as the standard model, needs photon conservation and an intergalactic property that changes in time, and brings KiDS into tolerance in the suite (62/8/7, registered, not adopted); the proposed relativistic action checked: sign, instantaneous constraint, reaction force, and χ in the metric gravitational waves feel** | **Bullet (22.5 kpc; 15 kpc in brackets), inside 250 kpc, 10¹⁴ suns:** no crossing heat 1.37/2.95; u 1.59/3.10; v_h 600 as proposed 2.54/3.37 (2.56/3.39), peaks 19–20/8–9 kpc, gas residuals 0.05/0.10; energy booked 1.75/3.10 (1.76/3.11); power into the fast glow ×2.5 / ×3: 2.26 / 2.40 (600 km/s), 2.27 / 2.41 (800 km/s), so about ×3.3 is needed. **Far collisions:** MACS J0025 SE/NW inside 300 kpc 2.05/1.85 (suite) → 2.24/2.03 (u), 2.58/2.50 (600 as proposed), 2.23/2.06 (energy booked), against 3.64/3.79 (large errors), NW peak 79 → 18–44 kpc; El Gordo inside 500/1,000 kpc 8.42/21.4 → 8.86/22.1 (u), 9.75/23.7 (as proposed), 8.84/22.1 (energy booked), against 9.45/24.3. **Supernovae (1,365):** η 0.449 (0.41–0.49), η = 1/2 Δχ² 1.09, χ² 1210.4 (flat ΛCDM 1208.3), η = 0 +100.9; energy-conserving stretch Δχ² +59.9; power of (1 + z) free 0.89. **Suite:** dist_fixed_eta 62/8/7, dist_metric_eta 62/8/7, dist_metric 57/13/7 (round 12: 59/11/7). **Action:** 105 vs 161 km/s at 10 kpc as written; λ mass 1.2 × 10⁻⁶; reaction 0.40–0.77 of round 2's; GW170817 delay from the Milky Way's χ 3.6 yr if waves ignored it. Scripts: `code/hot_mode_speed_v21.py`, `code/hot_mode_collisions_v21.py`, `code/distance_eta_v21.py`, `code/action_checks_v21.py`; `regression/candidates/dist_*.json`. |
| **R20** A second mechanism track: the energy and stress of a shared state (round 20) | **The owner's review of the step-back audit adopted as the next track (force = −⟨∂H/∂X⟩ in a shared matter–medium state; the rhythm models kept as the active-recoil control); the audit narrowed in four places; equilibrium media, stiffening or at their gap's closing, only change a length or give Newton's 1/r² (the review's stop rule: the long reach and √M need a non-equilibrium medium); in the linear medium a probe's pull is at best Newton-like in mass with a Yukawa range; the cold law's macroscopic target is AQUAL's deep equation; the cluster heat without X-ray input: X-COP 25% → 40%; no_hold on today's law 61/12/4 (registered, not adopted)** | **Supplied check reproduced:** rms force 1, 1.72, 2.15, 2.29 for 1–4 of 8 excited, 0 when all are; zero one-point amplitude. **Stiffening medium (32³, exact statistics):** range 3.18 → 0.93 as T 0.03 → 3 (β = 1), the harmonic estimate to 7%; force exponents ≥ 2.8 at r = 3; energized regions only steeper (minimum 2.9–3.2 against 2.2 cold); gap closed: source exponent 1 → 0, distance 2.1–2.8 → 2. **Probe (Ns 3–13, 1,260 cases, all attractive):** half filling: the kernel's Yukawa slope 4.9 → 8.5, mass exponent 1.1; one or two excitations: exponent 9.7 → 17, mass 1.0 / 0.35. **Clusters, stars' speeds from the law's own gravity:** rms 0.335 (0.306 with u 186); law ÷ hydrostatic gravity 0.83 → 1.51 at 0.1 → 1 R500 (X-COP's non-thermal pressure ≈ 6% at R500). **no_hold_r12:** dwarfs χ² 137.7 → 60.4, Draco/Ursa Minor unchanged, binaries +18% at 20,000 AU, Q₂ = 0. Scripts: `code/anharmonic_medium_v20.py`, `code/finite_population_probe_v20.py`, `code/xcop_selfconsistent_v20.py`; `energy-shift-v20/` (supplied); `regression/candidates/no_hold_r12.json`. |
| **R19** The medium tested against the data; the Milky Way refitted (round 19) | **The absorbing stream solved exactly (passive at every strength with wavelength-sized absorbers moving at the wave's speed; point absorbers drag; the ray form right at weak absorption); the data allow only weak absorption (the hot glow must reach inward: absorption length ≳ 300 kpc); there the models give the cold law (p = 0.58 ± 0.13, q = 1.11 ± 0.22) and the collision rule but not the heat gain; the Milky Way with independent matter follows the curves' shape, 3–6% slow (the law's SPARC offset); Cassini 0.2–0.4σ; a supplied blocker mechanism for screening and release, carried to the Solar System; a frozen lensing–dispersion prediction; NOVELTY.md** | **Exact medium (uniform stream):** backward transmission at 1 / 2 λ: 0.62 / 0.37 at κ = 0.5 (ray 0.61 / 0.37), 0.11 / 0.009 at κ = 2.3 for the cardioid (ray 0.10 / 0.010); floor a few % (cardioid) to 10–28% (inward rule) at strong κ; passivity min eigenvalue = the free medium's at every κ; self-force −0.015 to −0.13 P/c (SPARC allows ≈ 0.05 a); point absorbers +0.06 → +0.54 P/c with the cut-off (κ = 0.5). **Hot shell:** 69% of X-COP's S at 0.1 R500 from outside the receiver; strict one-way 0.41 (0.28 refitted; suite 32/6/11 refitted); a round source's inner shells heard exactly half by a fully absorbing stream; stream 3/Mpc refitted 35/8/6 (u 140), 10/Mpc 32/10/7, 30/Mpc 30/9/10. **Models (84 runs):** heat gain 0.90–0.97 at k = 2, push at k = 8; colliding 0.78–0.90; compact sources 2–15× weaker. **Milky Way (Feng 2026, 903 Cepheids; Eilers, Zhou, Ou):** χ² 35 / 34 / 35 / 62 for 12 / 38 / 34 / 35 points with K_z; scale 0.935–0.966; Newton 393; SPARC median +0.026 dex at g_N ≈ 10⁻¹⁰. **Cassini:** Q₂ 1.96–2.37 × 10⁻²⁷ at L = 0.15 pc; L ≥ 0.09–0.10 pc (1σ); binaries amended 1.024–1.028 / 1.057–1.068 (round-10 forecast kept). **Blockers:** release switch at L ln(load): 0.4–4 pc at 868 yr; same Q₂ at 6–60 yr, step to +12% at 5,800–7,400 AU. **Frozen:** −0.183 / −0.087 / 0 / +0.077 / +0.144 dex at σ_e = 100–300 km/s. Scripts: `code/full_wave_v19.py`, `code/self_force_sparc_v19.py`, `code/hot_shell_v19.py`, `code/emergent_law_v19.py`, `code/mw_joint_v19.py`, `code/blocker_release_v19.py`, `code/frozen_prediction_v19.py`; `screening-blockers-v19/` (supplied). |
| **R18** The review's checklist; the companion as a medium (round 18) | **Write-up matched to the evidence (postulates table, conditional collision rule, suite and scoreboard caveats); Cassini updated to the 2026 value (1.54σ). Step 1 started: waves carried by the companion's stream are one-way (derived) but push every emitter downstream (ruled out as is); waves the stream absorbs when they move against it are one-way, push-free and passive, and in the full model make the pull grow with heat at 73–84% of the square root of the model's own glow (about half the law's heat gain; corrected 25 September, STEP-BACK-AUDIT.md §5), without a single beat for warm sources** | **Carried waves (1D exact, checked by simulation to 1–4%):** push \|q\|²/(4(U² − c²)) = stream-frame power/U; P/F = (U² − c²)/U, galaxies need c ≥ U/√2. **Absorbing stream (full model, κ = 2.3 / 5):** R_net 0.62 / 0.73 (k = 8), 0.81 / 0.84 (k = 2); cold keeping step 0.64 / 0.72; absorbed 33–48%; passivity −8 × 10⁻⁴ / −8 × 10⁻³ of 4.8 (imposed mask −0.87 of 4.9); warm rhythm spread 3 × 10⁻³ (mask 4.6 × 10⁻⁵). **Heat factor:** 76–87% of 1 + k (dilute), 45–61% (dense). **Cassini 2026:** Q₂ 4.38 against (1.6 ± 1.8) × 10⁻²⁷; L ≥ 0.19 pc for 1σ. Scripts: `code/cassini_2026_v18.py`, `code/heat_factor_v18.py`, `code/flowing_medium_v18.py`, `code/absorbing_stream_v18.py`, `code/absorb_reduced_v18.py`. |
| **R17** Keeping warm matter in tune (round 17) | **Step A done: the damaging term found (the non-gradient half of the two-way tugs between warm pieces; removing it gives keeping step 0.88 against 0.63 cold). Step B: the review's pair reproduced, but its protection is lost in a shared wave; at equal glow reciprocity fixes the tugs in every internal structure; the velocity's time-reversal parity (a first-principles correction) and a widely detuned pair help (0.18 → 0.41 at k = 8). A one-way (outward-only) wave removes the damage completely in the reduced model (R = 1.0–1.5 at four distances, power ∝ mass); in the full, energy-booked model too: warm matter pulls distant matter harder than cold at four distances, R_net ≈ 1 on average, and colliding and stopped sources pull like cold ones; the glow is in proportion to mass and the pull grows as its square root. The one-way coupling is imposed; next, the companion as a flowing medium** | **Budget:** the full model's rhythms = the three budget terms to 1.5 × 10⁻⁷; C = C₀ + C₁ + C₂: own offsets 2.0 × 10⁻⁴ ⊕ 3.3 × 10⁻⁴ √k ⊕ 1.9 × 10⁻⁴ k, tugs 1.9 × 10⁻³ ⊕ 2.7 × 10⁻³ √k ⊕ 1.9 × 10⁻³ k. Diagnostics (k = 2 / 8): gradient part only 0.75 / 0.88, non-gradient only 0.15 / 0.00, no heat tugs 0.66 / 0.63. **Pair:** ×4.00000 per doubling at Δ = γ, shift 4.3 × 10⁻¹⁹, 1% → 0.0200, collisions 0.800 / 0.500; sensitivity at Δ = γ: separate channels 0.00, one shared channel 0.82 (1.67 / 1.32 / 0.32 at Δ = 0 / 0.5 / 2). **Structures** (12, equal glow): heat tugs 1.64–1.76 × 10⁻² at k = 8 in all; keeping step at k = 8: single 0.18, odd 0.33, pair Δ = 4(γ + γᵢ) odd 0.41, valve odd 0.21 (0.67 at k = 2; at k = 8, in one arrangement, it settles on a beat shifted by 3.1 × 10⁻³); sub-wavelength source ≤ 0.16. **One-way** (reduced): spread 6 × 10⁻⁵ at k = 8; keeping step 0.98 / 0.96 / 0.85 / 0.81 at r = 6 / 9 / 13.5 / 20; needs an inward leak ≲ 10%; intensity per piece constant for N = 24 / 48 / 96. **One-way, full model** (energy booked to 2 × 10⁻¹²): net pull at r = 6 / 9 / 13.5 / 20, k = 8: +2.36 / +2.12 / +1.75 / +1.07 × 10⁻⁴ against +1.61 / +1.14 / +0.41 / +0.53 cold; R_net 1.06 on average (single), 0.94 (odd); colliding (ν = 50) +1.53 / +0.99 / +0.38 / +0.50; stopped: back to cold; the sources' rhythm spread 4.6 × 10⁻⁵ at k = 8 (cold 1.7 × 10⁻⁴). With f = 0.05: the receivers' radiators take 8–11% of the pull; R_net 1.00 at k = 8, 1.12 at k = 2. **Two-way baseline** (full model): k = 8 pushes, −0.89 / −0.55 × 10⁻⁴ at r = 6 / 9. **Mass** (full model, 24 / 48 / 96 pieces): glow per piece 2.74 / 2.78 / 2.88 × 10⁻³ (two-way), 2.19 / 2.30 / 2.37 × 10⁻³ (one-way, f = 0.05); one-way pull ×1.71 at rest and ×2.11 at k = 8 for four times the mass (the law: ×2). Scripts: `code/rhythm_budget_v17.py`, `code/rhythm_protect_v17.py`, `code/paired_channel_v17.py`, `code/one_way_v17.py`, `code/one_matter_v17.py`, `code/one_matter_summary_v17.py`. |
| **R16** Round 15's next steps, one at a time (round 16) | **In progress. Step 1 done: the pull's energy bill uses the phase speed; fed back into the stream, the galaxies allow v_phase ≲ u/2. Step 2 done (conditional): the store opened by the velocity relative to the flow gives the cold leak and 1 + 3σ²/u² with one u; the crossing heat lifts the Bullet's smaller half 12–18%. Step 3a done: with the complete, energy-exact coupling the σ² and the collision rule survive (3.3–4.0 per doubling). Step 3b done: one kind of matter balances its energy, a cold source pulls, the recoil rule (1 − 2f) is exact, and a warm source's wave carries √ of its glow; but distant matter falls out of step with warm matter (the decisive test not yet passed)** | **Feeding feedback** (I = |g_N|(1 + 2βW/(aM)), β = v_phase/u; a, g_d refitted): SPARC 15.87 / 15.76 / 15.79 / 16.35 km/s at β = 0 / 0.2 / 0.5 / 1; outer residual slope on W/(aM) +0.074 ± 0.062 / −0.023 ± 0.034 / −0.054 ± 0.019 (β = 0 / 0.5 / 1); β fitted on 89 training galaxies = 0.95, worse out of sample (validation 19.88 against 19.10, test 13.80 against 12.40 km/s). X-COP: rms 0.222 / 0.223 / 0.228 with u = 170 / 142 / 113 km/s. So v_phase ≲ 85 km/s and the bill ≤ ℓ per kilogram; round 10's 0.05u was too strict. **Flow-opened store** (δ = χ(v − V_flow), γ₀ = 0): cold leak 2.0002 × 10⁻⁴ (expected 2 × 10⁻⁴); released/cold 1.16 / 1.70 / 3.88 / 12.7 / 48.6 / 198 at σ/u = 0.25 … 8 (law 1.19 … 193); collisions 5.90 / 2.96 / 1.07 / 0.383 (expected 6 / 3 / 1.09 / 0.387); drift (1 − v_r/u)² exact; rotation 1.000 with memory. **Crossing heat, Bullet:** M(<250 kpc) smaller half 1.385 → 1.556–1.641 × 10¹⁴ M☉ (target 2.47–2.85), main 2.969 → 3.073–3.160 (target 3.09–3.46). **Shared wave:** power from M = far-field flux (7.01166 / 7.01162), two monopoles λ/10: 1.9355 / 0.0645; dissipative-only 3.88 / 3.53 / 2.45 / 1.42 (t = 1000), complete 3.50 / 3.71 / 3.84 / 3.91 (t = 50, 20 pieces in 1λ), dilute 4.08 / 4.04 / 4.02 / 4.00; collisions 0.55 / 0.091 of free (γ/(γ+ν) = 0.5 / 0.091). **One kind of matter** (48 source pieces, 8 receivers of the same matter, f = 0.2, 3 arrangements): energy closes to 10⁻¹² (free) / 4 × 10⁻⁷ (collisions); recoil (1 − 2f) exact (one body: +1.050 / −1.050 / −1.126 × 10⁻⁴ at f = 1; net +0.833 × 10⁻⁴ at f = 0.1); cold: lead 0.55, net pull +1.21 × 10⁻⁴; warm k = 0.5 / 2 / 8 / 16: glow ×1.21 / 2.11 / 5.00 / 8.07, wave ×1.01 / 1.31 / 2.24 / 2.81, lead 0.26 / 0.27 / 0.27 / 0.16, net +0.43 / +0.25 / −1.04 / −2.92 × 10⁻⁴, sources' rhythm spread 7.7 × 10⁻⁴ … 3.9 × 10⁻³ (4.7 × 10⁻⁴ cold); collisions ν = 50: lead 0.44, net +1.07 × 10⁻⁴. f = 0.05: net +1.64 / +1.00 / +1.25 × 10⁻⁴ (k = 0 / 2 / 8), lead 0.40 / 0.19 / 0.21. Cold source detuned by hand (Δ₀ = 0.001 / 0.0025 / 0.004): lead 0.28 / 0.12 / 0.19 — the cause confirmed. Over 24,000 time units a cold source settles and its receivers lock at 0.91–0.92; a warm one (k = 2) never settles. Law and suite unchanged. Scripts: `code/feeding_feedback_v16.py`, `code/feeding_checks_v16.py`, `code/stream_store_v16.py`, `code/crossing_heat_v16.py`, `code/shared_wave_v16.py`, `code/one_matter_v16.py`. |
| **R15** An independent calculation joined to the local force (round 15) | **Done, narrowed in round 16 after an independent audit: the finite internal store (quiet mode opened by motion) reproduced byte for byte; the tested cloud gives its symmetries and Dicke narrowing but not its σ² (gapless collective darkness; an internal store is the candidate); round 10's quarter-cycle lead comes out of inverted self-sustained receivers' own equations; a first version of the chain, with δ = χw as input: pull ≈ √(cold + released), ×2 per doubling of σ where heat dominates. Open: the shared wave's energy budget (flux 0.6–1.4 × modal loss), the σ² with the feedback in, the full spatial law** | **Reproduction:** `independent-r15/` (4 of 5 files byte-identical, the fifth number for number): extra radiation ×4.00 … 3.01 per doubling, collisions 50.7 / 9.2 / 4.8% against γ/(γ + ν). **Wave only** (40 emitters, λ/4 ball, exact e^{ikR}/(kR) coupling): rotation and boosts release < 10⁻¹²; collisions suppress as ~1/(1 + ν/γ_eff); steady leak ×1.3–1.8 per doubling (σ^0.4–0.9) because the collective rates fill every decade 10⁻¹⁰–10. **Receivers** (one body at 6λ): passive pushed ∝ |E|²; amplifier pulled ∝ |E|²; the calculation's gain oscillator locks a quarter cycle behind, pushed ∝ |E| (−0.071 / −0.76 / −13.3); an inverted self-sustained (Bloch) emitter locks a quarter cycle ahead, pulled ∝ |E| (+0.0138 / +0.129 / +0.537), saturating at its pump rate (0.613); pull = P/v exactly. **Integration** (100 pieces, D monopole + B dipoles; 64 bodies at 6λ; 90 runs): released/cold = 1 + k to 3% for k ≤ 50; pull ∝ √(intensity) to 1–3% for locked and inverted bodies; extra-pull ratios per doubling 3.55, 2.95, 2.48, 2.20, 1.97 (law 3.60 … 2.15); collisions γ/(γ + ν) at the source to 1%, slow bodies lose step; no memory; pull/(P/v) 0.989–1.002; flux vs modal loss 0.6–1.4 per arrangement (mutual interference not fed back). **Cost:** F = P/v_phase (corrected in round 16; first written with the travel speed u): a·v_phase = 2βℓ per kilogram, ≤ ℓ with the galaxies' β ≲ 0.5; light-speed crests would cost 3,500ℓ (the Sun: 5 × 10⁻¹²/yr, against ephemerides' ~10⁻¹³/yr; as a simple budget v_phase ≲ 5,700 km/s; round 15's "2ℓ at v = u" was the identity a = 2ℓ/u). **Audit (round 16):** inverted bodies ×2.41 / 4.38 / 8.24 against √(1 + k) = 2.00 / 3.61 / 7.00; with energy-exact damping ż = −iHz − ½W†Wz a reduced 20-source test's doubling ratios fall from 4.13 … 3.91 to 3.80 / 3.28 / 2.02 / 1.64. **Constants:** γ₀ = ℓ/(2c²) = 3.0 × 10⁻²³ s⁻¹, γ/χ² = 9.7 × 10³² m²/s; the collision rule puts 1/γ between ~100 s and ~200 Myr. The law and suite unchanged (59 / 11 / 7). Scripts: `code/wave_dark_v15.py`, `code/receivers_v15.py`, `code/reservoir_force_v15.py`. |
| **R7** Full check against the Milky Way and the standard tests (round 7) | **Done: most pass; three shortfalls with candidate fixes** | **T1.1 Milky Way:** with McMillan (2017) matter our law gives 211 km/s at the Sun (Gaia: 229–234; MOND 223), and runs 3–13 km/s below the four Gaia curves from 15 to 27 kpc, where McMillan's dark-matter model is 12–23 km/s too high. The shortfall at the Sun matches our law's 0.03 dex offset in SPARC at the same pull; g_d × 1.5–2 gives 218–221 km/s, and SPARC's typical miss is then 15.80 km/s (now 15.85). Vertical pull at 1.1 kpc: 74 with McKee's counted local matter (68 ± 4, 74 ± 6; MOND 84). Mass inside 100 and 200 kpc: 6.6 and 12.5 × 10¹¹ (measured 6.1–7.3 and 11.0; MOND 8.0 and 15.6). Escape speed 507–522 km/s (measured 445–580). **T1.4 Dwarfs:** Fornax, Leo I, Leo II and Sculptor agree within 1–2σ; Draco, UMi, Sextans, Carina, Crater II and Antlia 2 are 1.5–5× too slow with the Galaxy's field (MOND, computed identically, is 10–20% higher). Without the external hold, Crater II gives 3.4 km/s (2.7 ± 0.3). **T2.3 Solar System:** planets, S2, pulsars and light bending are exactly GR. Cassini's Q2 = 2.4–3.1 × 10⁻²⁶ s⁻² (measured (3 ± 3) × 10⁻²⁷; MOND 2.8–3.1 × 10⁻²⁶); gradual release over L ≥ 0.15 pc brings it to ≤ 4.6 × 10⁻²⁷ and makes wide binaries 1–5%. **Lensing:** Einstein Cross and microlensing are stars only. KiDS-1000 absolute lensing RAR: spirals −0.005 dex, ellipticals +0.021 (MOND −0.073, +0.113). Scripts: `code/milky_way_v7.py`, `code/mw_dwarfs_v7.py`, `code/strong_field_v7.py`, `code/lensing_census_v7.py`, `code/transition_check_v7.py`. |
| **WB** Wide binaries (round 6) | **Predicted** | Exact orientation-averaged solution of our field equation in the Galaxy's field at the Sun (1.58 × 10⁻¹⁰ m/s², half the companion released): 19% more pull than Newton beyond 7,000 AU; 1.14–1.26 for ±20% in the Galaxy's field. MOND (simple) 43%. Data disputed (Chae ~1.4; Banik et al. Newton). `code/wide_binaries_v6.py`. |
| **T2.1** Bullet Cluster (round 5, subcluster size superseded by round 6) | **Pattern, speeds and lensing masses reproduced; subcluster size to test** | Main cluster: half of round 4's speed gap was the model cutting its galaxies at 1.5 Mpc; outer stars at star/gas 0.048 at R500 (X-COP 0.035–0.073) give 1,249 km/s and a lensing mass of 2.42–2.54 × 10¹⁴ inside 250 kpc (obs 2.5 ± 0.1, 2.8 ± 0.2). Clowe et al.'s κ are lower bounds, so round 4's "within 1σ" and its lighter stars are withdrawn. Subcluster: its strong-lensing mass (2.0–2.3) needs a pre-collision size of about 1:3 in visible matter (1.93 at 1:3; 1.04 at 1:8); its original galaxies then moved at about 800 km/s (Barrena et al.'s X-ray estimate: about 700). Predictions: about 7 × 10¹² M☉ of its stars now around it; star/gas about 0.05 in the main cluster's outskirts. `code/bullet_main_v5.py`. |
| **T2.1** Bullet Cluster (round 4, strength readings superseded by round 5) | **Pattern and strengths reproduced; main-cluster galaxy speeds open** | The companion is slow (u = 197 km/s) and keeps its emitter's velocity, so each cluster is still wrapped in the companion of the settled cluster it was before the collision, riding with the galaxies. Only a sphere of radius u·t (30 kpc) has been rebuilt around the stopped gas. No new rule; nothing fitted to the lensing. Subcluster: 0.175 (0.13–0.25 for pre-collision ratios 1:6–1:10), against 0.20 ± 0.05. At star M/L 1–1.5 all four strengths agree within about 1σ. Barrena et al. 2002 independently argued the subcluster is the stripped core of a cluster with a pre-merger ratio of about 1:6. Collision stack with memory: β = 0.03 (−0.01 to 0.12), against −0.04 ± 0.07. Prediction: the subcluster's stars move at 460–610 km/s (7 galaxies now give 212 ± 60). **Open:** the main cluster's galaxies come out about 20% slower than the 1,249 ± 100 km/s measured (2.5σ). Scripts: `code/bullet_v4.py`, `code/collisions_v4.py`, `code/bullet_speeds_v4.py`, `code/stream_tidal_v4.py`. |
| **T2.1** Bullet Cluster (round 3, strengths superseded by round 4) | **Pattern reproduced; strength open** | On Clowe et al. 2006's published masses, both lensing peaks sit on the galaxies (8 and 34 kpc from the BCGs) and the gas residuals match (0.04 vs 0.05 ± 0.06 and 0.02 ± 0.06). The old direction rule, hot gas, or no heat each put the peaks on the gas. Strengths: main 0.51 (obs 0.36 ± 0.06; 0.31–0.41 at M/L 1–1.5), sub 0.07 (obs 0.20 ± 0.05). **The subcluster strength is the open item.** 20 modelled collision pieces all keep lensing with the galaxies (Harvey et al. 2015: 5.8 ± 8.2 kpc). Scripts: `code/bullet_v3.py`, `code/collisions_v3.py`. |
| **T2.1** Bullet Cluster (round 2, superseded) | **Partly; top open problem** | Companions that keep their emitter's velocity move the lensing peak 42–75% of the way from gas toward galaxies in a toy merger. But Harvey et al. 2015 (72 collisions) keep lensing on the stars, which needs a much slower companion rebuild on shocked gas than u = 874 km/s gives. `code/bullet_toy.py`. |
| **T0.1** Resolve the contradiction | **Done** | Branch (a) fails on energy by 1.5 × 10⁴; the field's own energy is 1.7 × 10⁻⁴ of what it would need. Branch (b), a modified propagator, survives. τ splits into two objects: in lensing it is the **gravitational slip η** (legitimate, no energy cost); in cluster dynamics it has **no mechanism**. The cluster "closed within scatter" claim is **retracted**. |
| **T2.2** GW170817 | **Done** | Survives, with a construction constraint: the modification and the slip must live in the scalar sector. TeVeS-type completions are excluded; a surviving class exists. |
| **T0.2** Field equations | **Done**; relativistic home identified | Three equations. The dynamics equation reproduces the fitted law to 4.4 × 10⁻¹⁶ and is Milgrom's QUMOND, attributed. The slip equation is ours; JR-10's γ_χ is identically η. **Relativistic:** slip is exactly a radial stress of the extra gravity (the radial Einstein equation). Of four places the extra gravity could live, one passes GW170817 and our data: the aether-scalar-tensor theory (which gives η = 1, verified from the paper) plus a conformal coupling κ of matter to its scalar, giving η = (1−κ)/(1+κ). Clusters need κ ≈ −0.12 to −0.17. What sets κ is open. |
| **T1.2** Cluster lensing masses | **Done, via the published bias** | The X-COP hydrostatic bias implies η = 1.27–1.42, overlapping the under-bent SLACS lenses (1.36–1.54). Degenerate with non-thermal pressure in one comparison; a relaxed-versus-disturbed split at matched f separates them. A published slip measurement from galaxy orbits, which gas pressure does not affect (Pizzuti et al., MACS J1206), gives 1.01 +0.31/−0.28 against our predicted 1.21–1.37. That is consistent at about 1σ. |
| **T1.3** Groups | **Done, inconclusive** | Three SL2S groups with real weak-lensing masses point to η > 1 at face value, but a 30% dispersion bias the source paper itself reports erases it. **Correction:** an earlier "lensing masses ~50% above dynamical" attribution was not in the paper. |
| **T2.3** Solar system (slip part) | **Done** | Slip moves the Cassini γ by at most 3.7 × 10⁻¹², against a bound of 2.3 × 10⁻⁵. Still to check, separately: the simple interpolation function's high-acceleration tail against Cassini (Hees et al. 2016 rule out several popular choices). |

Full working in [THEORY.md](THEORY.md); relativistic slots in
[research_work/results/relativistic-slip/](research_work/results/relativistic-slip/README.md).

**Next, after round 24 (hot-companion gravity):** results README §35.4. For the field-theory route: (1) done in round 24: the repaired equation solved exactly (damped iteration) keeps 58 / 11 / 8; a faster, convergent solver (Newton-type) would let the exact form join the suite's standard tier; (2) derive the non-flowing brightness term from the companion's own dynamics (u U_X − \|J_X\| for the hot glow, and why the cold glow counts as pure flow); (3) put the strong-field hold into the field theory; (4) test B: the canonical pair-emission rate, which hinges on whether the heat factor saturates near light speed; (5) the owner's decision on RULES.md for an AQUAL-form action, or its derivation from the mechanism. Round 23's list (below) continues.

**Next, after round 23 (hot-companion gravity):** results README §34.4. The list after round 22 (below) continues, with the power the fast glow needs now about 2.2× (flicker-driven heat of crossing) and one new item: derive how the quiet store responds to a flickering companion (heat ∝ the flicker's depth m), and test the pattern it predicts: in unequal collisions the smaller cluster carries relatively more extra lensing.

**Next, after round 22 (hot-companion gravity):** results README §33.5; round 21's list (below) continues.
1. **The fast glow in a microscopic medium:** a medium with the two-regime dispersion (slow crests at the companion's
   frequency, group velocity rising about 3.5-fold up to the crest-crossing frequency), driven at that frequency by a
   crossing star; compute the hot glow's speed and power from first principles with every watt booked.
2. **The quiet store's capacity:** at least 4 × 10¹³ J/kg for the Bullet's stars (§33.2), in a store model.
3. **The rest of round 21's list:** the distance law's path factor derived; the relativistic action revised; the
   energy-shift track out of equilibrium; the cluster heat without X-ray input; the frozen lensing test.

**Next, after round 21 (hot-companion gravity):** results README §32.5; round 20's list (below) continues.
1. **The fast glow from collisions, derived:** v_h/u and the power the collision puts into that glow from the
   medium's dispersion, with every watt booked (§32.1 needs about 3.3× the rate the heat rule gives); then the
   72-collision stack with the heat of crossing.
2. **The distance law's path factor, derived:** η = 1/2 and a photon-conserving stretch from one intergalactic
   process, and what in intergalactic space changes with time (7.5 × 10⁻¹¹ per year, not in laboratories); the
   geometry (fixed or metric) tested with standard rulers and surface brightness.
3. **The relativistic action, revised:** the corrected sign, the constraint along the companion's world lines, I and
   s as the companion's own fields, and χ in the metric that gravitational waves feel.
4. **The energy-shift track, out of equilibrium** (round 20): a medium carrying the companion's steady outward
   energy flux, with every watt booked; exponents in source mass and distance against §31.5's target.
5. **The cluster heat without X-ray input:** independent velocity-dispersion profiles of cluster galaxies and
   weak-lensing masses (§31.6), with orbital anisotropy.
6. **The frozen lensing test** (§29.7).

**Next, after round 20 (hot-companion gravity):** results README §31.8.
1. **The energy-shift track, out of equilibrium:** a medium carrying the companion's steady outward energy flux; its
   energy density, flux, stress and correlations around a source; a probe's energy shift and force with every watt
   booked; exponents in source mass and distance against §31.5's target.
2. **Population dynamics:** how the shared state's occupation is established, maintained and depleted, with an energy
   account (not identified with heat by instruction).
3. **Motion with its recoil:** moving sources and receivers with forces, propagation and backreaction evolved
   together, cold, free, colliding and rigidly moving from the same internal resources.
4. **The cluster heat without X-ray input:** independent velocity-dispersion profiles of cluster galaxies and
   weak-lensing masses (§31.6), with orbital anisotropy.
5. **The frozen lensing test** (§29.7).
6. **Registered comparisons:** `no_hold_r12`; alternative static distance laws (brightness, angular size and event
   duration together).
7. **The owner's decisions** (audit §4.3–4.5): the field equation (either form derived from the energy balance is in
   the MOND family), the standard for no-hold, the distance formula.

**Next, after round 19 (hot-companion gravity):** results README §29.8 and BLOG §9.

1. **Keep a warm source in step without strong one-way absorption** (the data allow only weak absorption, and there the
   models give no heat gain): the reservoir-engineered one-way coupling inside a source with a two-way far field;
   receivers whose locking survives a detuned source; the models' collective states in compact sources.
2. **The radial medium exactly** (partial waves), κ from the stream's density; the self-force's scaling.
3. **The blocker candidate in the matter model:** a metastable blocking excitation, its formation rate and lifetime
   measured separately, the force rule; then wide binaries (a step against a gradual rise).
4. **The law's 3% at the Milky Way's pull** (SPARC's median +0.026 dex) with a robust fit statistic.
5. **The frozen prediction** tested on σ-split lensing; inertia; likelihoods; a reproducible release.

**Next, after round 18 (hot-companion gravity):** results README §28.6 and BLOG §9, following the review of
25 September 2026 (its items in dependency order).

1. **Strong absorption as a wave problem:** the absorbing stream's full wave equation between pieces closer than a
   wavelength, and whether a passive medium with strong absorption reaches the law's square root and a single beat.
2. **What sets the absorption:** a kinetic model of wave quanta scattering off the streaming companion, energy and
   momentum followed in both.
3. **The one-way law it implies,** tested before any fit is changed: S and the pull's direction counted over what the
   medium hears (the hot-shell benchmark), against the X-COP profiles and the collision maps, constants fixed.
4. **Mass scaling and distance with the medium:** the exponent with its uncertainty, physical mass and resolution
   varied separately.
5. **Inertia and universal free fall;** then the frozen law, likelihoods and fair baselines; a reproducible release.

**Next, after round 17 (hot-companion gravity):** results README §27.5 and BLOG §9.

1. **The companion as a flowing medium:** derive its wave in its own outflow (the crests' speed against the stream;
   whether a source's companion streams as one medium; the energy and momentum the stream exchanges with the pieces).
   This decides whether the wave is one-way, the property that kept warm matter in tune in round 17's tests, and gives
   the coupling to use in place of the imposed one.
2. **Step C with that wave, in the full model:** R(σ, r) at more distances; cold, free, colliding and stopped sources;
   mass scaling (physical mass and refinement separately); receivers of different kinds; the stream's share of energy
   and momentum booked explicitly.
3. **What inertia is** for a piece of this matter, so that "every kind of matter falls alike" can be tested; the
   radiators' share f common to all matter (or f ≪ 1/2).
4. A physical quiet store with a gap (χ, γ, and u from them); the release factor from drained receivers; the Doppler
   shifts of real motion.
5. Only then the astrophysical fits.

**Next, after round 15 (hot-companion gravity):** results README §25.5 and BLOG §9; reordered in round 16 after an
independent audit.

1. **One kind of matter, energy-exact:** every piece both source (a quiet store opened by motion) and receiver
   (inverted, self-sustained, powered by the same finite store), one local coupling for emission, loss, response and
   force, the mutual waves fed back from the start (ż = −iHz − ½W†Wz); cold, free and colliding motion; no
   strong-field switch.
2. **A physical quiet store with a gap:** an internal pair of oscillations that relative motion detunes (the motional
   mixing of a metastable atomic state is the laboratory analogue); derive χ and γ, and so u = √(γ₀γ)/χ.
3. **The release factor from the receivers:** whether a receiver drained by its neighbours' loud waves switches its
   pull off, rather than capping it.
4. Round 14's list: the two KiDS analyses; round 12's list.

**Next, after round 14 (hot-companion gravity):** results README §24.4 and BLOG §9.

1. **The release's steepness:** larger clouds, a loss rising faster than |E|², emitters of several strengths; the
   target is released/held ∝ σ² over a factor of ten in σ.
2. **The release factor from the loss:** the singers' hold in place of exp(−|g_N|/g_d) across the whole suite; why it
   stops at zero.
3. **The two KiDS analyses** (Brouwer et al. against Mistele et al.) and the heat they want for ellipticals.
4. Round 12's list: the guided companion with moving sources, the lenses' gas, the distance law's shape, MACS J0025's
   star masses, the faint dwarfs, the Sun's speed, Abell 1689.

**Next, after round 13 (hot-companion gravity):** results README §23.10 and BLOG §9.

1. **Derive the strength-dependent offset** from a saturable emitter (gain that saturates, loss that does not), and
   test whether the same saturation gives the law's release factor exp(−|g_N|/g_d).
2. **Bodies that keep in step** with a fast-changing wave: a pull drawn from the companion's energy flow (the guided
   stream), or a cold state dark enough that the whole heat-dominated regime lies at slow speeds; map u onto the
   locking rate and the darkness together.
3. **The release's steepness:** the toy's release flattens (σ^1.7 → σ^0.8); the data want about σ² over a factor of
   ten in σ.
4. **The heat exponent p = 1.75 on the full suite:** it removes KiDS's common level with the clusters unchanged.
5. Round 12's list, still open: the guided companion with moving sources, the lenses' gas, the distance law's shape,
   MACS J0025's star masses, the faint dwarfs, the Sun's speed, Abell 1689.

**Next, after round 12 (hot-companion gravity):** results README §22.5 and BLOG §9.

1. **The guided companion with moving sources:** does a companion carried with its matter keep the
   memory the collisions use, and does the Bullet's smaller half get its lensing back?
2. **Weigh the lenses' circumgalactic gas** (X-ray stacks for early types, UV absorption for late
   types) and predict KiDS's level from it; check the gas's effect on SPARC's outermost points.
3. **The distance law's shape:** the supernovae's best α drifts with depth (0.967 → 0.944); the
   bounded beam-area term removes it (η ≈ 0.45). Derive it, or find what else does.
4. **MACS J0025's star masses** from infrared light (its galaxy speeds prefer the published ones).
5. Items 5–8 of the round-11 list below, still open.

**Next, after round 11 (hot-companion gravity):** results README §21.7 and BLOG §9.

1. **The dynamical toy** of full-speed waves that cancel head-on while weakening away from the
   source: does a local rule give the curl-free flow that the energy argument needs? A
   self-scattering (diffusing) companion is one candidate.
2. **The distance law's scale from galaxy lensing:** refit α jointly on KiDS, SPARC's Hubble-flow
   distances and the light-transport fit that set it.
3. **The ellipticals' heat:** clusters want u = 163 km/s, the KiDS early/late gap about 200. Test
   a bottom-heavy IMF for massive ellipticals (spectra: about 1.85 × Chabrier) in X-COP and KiDS
   together, and the dispersions assumed for KiDS's red lenses.
4. **MACS J0025:** its collision age and star masses (the NW peak is a knife edge).
5. Items 1–8 of the round-9 list below, still open.

**Next, after round 9 (hot-companion gravity):** the full plan is in the results README §19.3
and BLOG §9. Every item ends in a regression-suite run.

1. **The dwarfs' hold.** Derive why the Galaxy's pull weighs less inside a separate system,
   consistent with the heat rule, and test it on the ten dwarfs. Check tides for Draco and
   Ursa Minor.
2. **The Bullet subcluster:** crossing heat as a collisions candidate. Also check item 1's rule
   on it.
3. **The Milky Way's matter:** the disk scale length and the bar mass as suite options.
4. **Abell 1689:** model it and add it to the suite.
5. **Cluster stars:** El Gordo's infrared stellar masses; Abell 520 with two main subclusters.
6. **Mistele's spirals:** grade against the lensing profiles directly.
7. **Wide binaries:** Gaia's next release measures the release length.
8. **Theory:** the field equation with the release length; what L is; u and g_d.

**Earlier list, after round 8 (hot-companion gravity):**

1. **A speed-dependent external hold** (loss of step between companions moving past each other
   at ≳ u), which the suite shows is needed to help the dwarfs without over-boosting wide binaries.
2. **Adopt gradual release** (costs nothing in the suite) after a full-tier run and refit. *Done in
   round 9; no refit needed.*
3. **The Milky Way disk's scale length** as the lever for the Sun's speed.
4. **Independent stellar masses** for El Gordo and Abell 520 (near-infrared, spectroscopy).
5. **Abell 1689**, and MACS J0025's elongated gas.

**Earlier list, after round 7 (hot-companion gravity):**

1. **Refit with the three candidate amendments** (gradual release L ≥ 0.15 pc; g_d × 1.5–2; a weaker external hold for
   bodies moving past each other at ≳ u), then re-run the coverage table (BLOG §7).
2. **Model the recorded collisions:** Abell 520 (the disputed dark core vs our fresh-sphere memory rule), MACS J0025,
   El Gordo, and the dense cluster Abell 1689.
3. **The Galaxy's hot flow near the Sun**, which sets the direction of the extra pull there (radial vs vertical).
4. Everything in the round-6 list below that is still open.

**Next, after round 6 (hot-companion gravity):**

1. **The Bullet subcluster's lensing mass** (twice our law's): model the heat its galaxies gained
   while crossing the main cluster, carried with them by the memory; look for diffuse
   starlight moving with it; a few hundred more velocities around it.
2. **A field theory for the companion,** with a momentum flux, its travel time (a retarded
   source) and a relativistic form, so lensing is derived.
3. **Microphysics:** derive `u` and `g_d`; pin the companion wavelength and locking time;
   check that cluster gas scatters faster than about 0.1–1 Myr.
4. **Wide binaries:** a proper external-field prediction.
5. **Tests others can run now:**
   * spectroscopic IMFs of the six SLACS lenses (1.05–1.35 × Salpeter);
   * the KiDS offset flat beyond 100 kpc;
   * galaxy and gas dynamics equal in relaxed clusters;
   * a larger spectroscopic sample of the Bullet's subcluster (460–610 km/s predicted);
   * extra lensing around the smaller clump in older collisions.

**Superseded list (22 September, built on the retired MOND law and slip), kept for the record.
Next, in order of value per effort:**

1. **Stack galaxy-orbit against lensing slip in CLASH-VLT clusters.** This is pressure-free,
   builds on a published single-cluster result, and tests η directly.
2. **Relaxed-versus-disturbed hydrostatic bias at matched f.** This separates slip from
   non-thermal pressure using public catalogues.
3. **The weak-lensing acceleration relation by galaxy type (KiDS-1000).** Brouwer et al. report a
   ≥6σ early/late difference that a universal law cannot make. A type-dependent κ can,
   and this is the most direct handle on what sets κ.
4. **T3.5, the stellar-mass convention**, before re-reading the three over-bent lenses.

**What changed about the paper.** The headline is no longer "τ solves the cluster
problem." It is: **the gravitational slip of the extra gravity is not 1**, measured
object by object with every constant frozen. That is smaller, contradiction-free, not a
rediscovery, and falsifiable. The cluster dynamical deficit stays open and is reported
as open.

---

## The strategic read, before the list

**The galaxy-scale result is not publishable as new.** Our `a₀ = 1.171 × 10⁻¹⁰ m/s²`
and our interpolation function are Milgrom's, and the acceleration relation we
reproduce was published by McGaugh, Lelli and Schombert in 2016 from the same SPARC
data. A paper whose headline is "we fit rotation curves with one constant" will be
desk-rejected as a rediscovery, correctly.

**The novelty is the slip η, and only η.** ~~The claim worth a paper is that the
residual factor of two blocking MOND in clusters is the Tolman active-mass factor of
the gravitational field.~~ *Superseded by T0.1 — see the status table above.* The
claim worth a paper is now: *the extra gravity bends light and moves matter through
potentials that differ, by a slip η ≠ 1 measurable object by object, where every
relativistic MOND theory and every dark-matter halo predicts η = 1.*

**Nature is the wrong first target.** Nature wants a result that is both novel and
settled. Ours is novel and *suggestive*: eleven clusters, a factor that needs a
hydrostatic-bias correction sitting at the edge of the measured range, and a theory
that is not yet written down as field equations. The right sequence is a strong
MNRAS or Physical Review D paper first, then Nature only if someone independent
confirms it. Aiming at Nature now costs a year and produces a rejection.

**One strategic amputation.** The "no expanding universe" position is not required by
any result we have — every test here is local, at z < 0.3, and none of them touches
cosmology. Carrying it into the paper means owing the referee a replacement account
of the CMB, big-bang nucleosynthesis and structure formation. That is a decade of
work attached to a result that does not need it. **Drop it from the paper.** Keep it
in the notebook if you like, but a paper that claims it will not be read past the
abstract.

---

## Tier 0 — Blockers. There is no paper until these are done.

### T0.1 Resolve the contradiction at the centre of the argument

This is the single most important task on the list, and a referee will find it in
five minutes.

Section 3.2 of the notebook argues the extra gravity is **not matter** — that there
is no reservoir, that the force is the source's own field with a longer reach, and
that this is why the energy catastrophe dissolves. Section 2 then uses **τ**, which
is the Tolman active mass of a **stress-energy tensor**.

Those cannot both be true as written. Either:

* **(a) the extra gravity carries a real T_μν** — in which case τ applies cleanly,
  but the 50,000× energy problem comes straight back and must be answered; or
* **(b) it is a modified propagator** — in which case the energy problem stays
  dissolved, but τ has to be *re-derived* as a property of the modified field
  equations rather than borrowed from general relativity's Tolman formula.

(b) is the more promising branch and the one consistent with the rest of the work.
It is real theory, not bookkeeping: what plays the role of "streaming versus
standing" when there is no substance to stream?

**Done when:** the paper can state in one paragraph what the extra gravity *is*,
and τ follows from that statement rather than from analogy.

### T0.2 Write down field equations

A covariant action, or at minimum a set of field equations, that reduces to:

* the interpolation function in the static weak-field limit;
* τ = 1 and τ = 2 as its two limiting stress states;
* a transition between them.

Without this the work is a fitting formula with a physical story attached. PRD will
not take it; MNRAS might, but the referee will ask anyway, and the answer "we don't
have one" caps the paper's reach permanently.

**Note:** T0.2 done properly may hand us T0.3 for free — if the field equations
contain the switch, we stop hunting for it empirically.

### T0.3 Derive what sets τ, or demote it honestly

Four candidates tested and ruled out: local starlight flux, galaxy mass/size/
luminosity/gas fraction, surface density, and support mode (blocked by a stellar-mass
convention). Only *scale* survives as a discriminator, and we cannot yet say why.

Two acceptable outcomes. Derive the switch from T0.2. Or present τ honestly as a
**measured two-valued parameter** — weaker, but still publishable if the measurement
itself is clean, which is what Tier 1 is for.

Not acceptable: fitting a crossover function to the gap between two populations and
presenting it as a prediction. We did that once already (notebook §3.5) and had to
retract it when raw data replaced a smoothed model.

---

## Tier 1 — Observational work the claim rests on

### T1.1 The Milky Way, properly — *done for hot-companion gravity in round 7 (see the status table); the text below is the earlier MOND-law era*

**Already run today, with nothing refitted:**

| | RMS vs 542-Cepheid rotation proxy |
|---|---:|
| Visible matter alone | 75.24 km/s |
| **This law, `a₀` from SPARC only** | **15.71 km/s** |
| Repo's earlier model, *fitted to this data* | 8.76 km/s |

Mean fractional error 6.4%. But **every one of the twelve bins is under-predicted**,
by 7 to 22 km/s. That is a systematic, not scatter, and it needs explaining before a
referee does it for us — most likely the baryonic mass model or the Jeans proxy,
possibly real.

**Still to do:** a modern Gaia DR3 rotation curve rather than a Cepheid proxy; the
**vertical force K_z(z)** at the solar radius, which is a genuinely different test of
the same law and one the literature takes seriously; and a baryonic model with
propagated uncertainties.

### T1.2 Cluster masses from weak lensing, not X-rays — *not optional*

The headline cluster result currently needs a hydrostatic bias of b ≈ 0.17–0.195
against a measured 0.10–0.15. That is the same size, but it sits at the edge, and the
whole claim leans on it. Weak-lensing masses for the same eleven clusters remove the
correction instead of applying it. **This is the task that most improves the paper per
unit effort.**

### T1.3 Groups — the decisive measurement

The only systems between the two regimes. With four candidate variables eliminated we
now know what to measure against.

### T1.4 Dwarf spheroidals — the test we might fail (*run in round 7: 4 of 10 agree; the Galaxy's hold is the issue*)

Deepest into the low-acceleration regime, and where MOND-like laws have known
tensions. If the law breaks here we need to know before the referee does. Run it even
though — especially though — it may hurt.

### T1.5 The reserved SLACS lenses

The KCWI release carries 14; only six are usable here. The rest are a genuine holdout
for the redshift trend, and the way to keep them one is to declare the prediction in
writing before computing a single angle.

### T1.6 SPARC with standard quality cuts

Our 16.25 km/s uses no inclination or quality flags, deliberately, so the comparison
against our own old model stayed like for like. For publication it has to be run the
way the field runs it, so the numbers are comparable to what is already in print.

---

## Tier 2 — The four things that kill no-dark-matter papers

Referees ask all of these. Silence on any one is fatal.

### T2.1 The Bullet Cluster
Mandatory. The canonical argument against modified gravity. The paper needs a
position, even if the position is "this is outside our scope and here is why."

### T2.2 Gravitational wave speed — GW170817
The simultaneous arrival of gravitational waves and light killed a large class of
modified-gravity theories overnight. Whatever comes out of T0.2 must be checked
against it. Do this *early* — it can invalidate a theory before it is written.

### T2.3 Solar system and binary pulsars (*run in round 7: GR exactly, except Cassini's Q2, which needs gradual release*)
The law must return to Newton and to standard general relativity where those are
tested to high precision. The interpolation function does this by construction; the
τ mechanism needs checking.

### T2.4 Cosmology — scope it out
See the strategic amputation above. The honest move is to state plainly that the
paper addresses galaxy and cluster dynamics and takes no position on cosmology.

---

## Tier 3 — Rigor the reviewers will check

* **T3.1** Full posteriors with propagated uncertainties, not point fits. Every
  headline number needs a credible interval.
* **T3.2** A blind analysis on reserved data, declared before unblinding.
* **T3.3** A formal prior-art section comparing against MOND, TeVeS, EMOND,
  emergent gravity and superfluid dark matter — stating exactly what is ours. The
  register at `research_plan/prior-art/` is the right home and needs updating.
* **T3.4** Independent reproduction: a clean public repository someone else can run.
* **T3.5** Fix the stellar-mass convention (notebook §3.9). It currently moves
  results by a factor of 2.4 and is the largest known systematic in the programme.

---

## Answers to the three questions asked

**Do we need first principles?** Yes, and it is the top of the list. T0.1 and T0.2
are the difference between "a formula that works" and "a theory." Without them the
ceiling is a solid MNRAS paper; with them, PRD and a real claim.

**Do we need more data?** Yes, but less than you might think, and specifically:
weak-lensing cluster masses (T1.2) and groups (T1.3). Everything else on the list is
confirmation rather than foundation.

**Do we need the Milky Way?** Yes — and it partly works already, with zero parameters
refitted. Finishing it properly with Gaia and the vertical force is one of the
cheaper high-value tasks here.

---

## Suggested sequencing

1. **T0.1** — resolve the contradiction. Everything downstream depends on which
   branch is taken, and it costs thinking rather than compute.
2. **T2.2** — check the gravitational wave constraint early, before investing in a
   theory it might rule out.
3. **T1.2 and T1.3** in parallel — the two measurements that most strengthen the
   cluster claim.
4. **T0.2** — field equations, informed by 1 and 2.
5. **Paper 1** to MNRAS or ApJ: the cluster τ result, scoped to galaxies and
   clusters, cosmology explicitly out of scope.
6. **Paper 2** to PRD: the theory, if T0.2 succeeds.
7. **Nature** only after someone independent reproduces the cluster result.
