# CL-1: the written-track response at cluster scale and for light — the first domain test

Declared before execution, 19 September 2026. Baseline: `main` at 3a80fec (RUT-1 stage 8). It runs beside the RUT-1 line and edits nothing in it; stages 0–8 stay frozen.

**The question.** Goal 8 of the [collective-memory goals](../../../research_plan/collective-memory-goals.md) asks whether the same physics behaves sensibly for other motions and for light, and requires an actual propagation rule for light: "an extra acceleration of massive bodies is not a lensing prediction without a derivation." Goal 5 asks for one source and one field to predict forces at many positions. CL-1 takes the written-track response exactly as RUT-1 left it — the two-stage equation, writing proportional to mass, the Gaussian footprint — and asks what its **steady state** does to the gravity of galaxy clusters and to light, with the data already in this repository. It is the first test of the model outside the planar orbital problem, and it is built so that it can fail.

**Why the steady state is the right object.** Every equilibrium system probes the mature field. For a source that does not change in time, both response equations give E = S and C = τ_keep·S whatever τ_form is (RUT-1 stage 3, declared and verified), and a memory of density is blind to the motion inside a steady distribution (PM-1 candidate A, PM-3). So a cluster in equilibrium, or a lens galaxy, tests the response law's steady state alone — the writing rule, the footprint and the retention — with no formation history entering. Formation and the response times are the RUT-1 line's business and are not tested here.

**What this is not.** No spectrum, decay law, energy budget or causal propagation is derived. The light rules below are *declared*, not derived, exactly as RPG-1's R4 was, and CL-1 tests what each declaration implies. The data are exposed (X-COP, SLACS and Coma have all been used before in this repository), so every observational reading is exploratory and is labelled so; nothing is frozen against an unexposed sample and no gate is an observational gate. "Exploration first, then frozen verification" is the rule this stage follows.

## 1. The physics, fixed

**The footprint** is the three-dimensional isotropic Gaussian K_w(x − x′) = exp[−|x − x′|²/(2w²)], whose restriction to a plane is RUT-1's footprint. **Writing is proportional to mass**: a continuum of density ρ writes S(x) = α·(ρ ∗ K_w)(x). **The steady state** of the two-stage response is C = τ_keep·S, with Φ_mem = −C and a_mem = ∇C. Two constants therefore enter, and only two:

    Λ ≡ α·τ_keep     the amplitude, in (km/s)² per M☉,
    w                 the footprint width, in kpc,

and it is convenient to write ℓ ≡ G/Λ, a length: a point mass M writes a well of depth ΛM = GM/ℓ, so **w/ℓ = Λw/G is the dimensionless writing strength**. In RUT-1's units (G = M = 1) the mature-ring control of stage 1 — 10.03% extra support at w = 0.1 — has w/ℓ = 0.501 in the narrow-track limit, and that number is recomputed here from `rut1.py`'s own functions so the two stages speak one language.

**A universal length.** w and Λ are universal constants in this stage. A width that scales with the system is a different law (RUT-1 stage 0: "a different width scaling … would be a declared mechanism to test"), and goal 8 forbids a label such as "galaxy" or "cluster" from acting as a switch. What each system *would* prefer is reported as a diagnostic reading, never adopted.

**The spherical source.** For a spherical density with cumulative mass M(<r), averaging the footprint over a shell gives the closed form

    C(r) = Λ ∫ k_w(r, r′) dM(r′),     k_w(r, r′) = exp[−(r − r′)²/(2w²)] · (1 − e^{−2u})/(2u),   u = r r′/w²,

which is exp[−(r² + r′²)/2w²]·sinh(u)/u written without overflow. Its radial derivative and Laplacian are taken analytically (the derivative of sinh u/u is (u cosh u − sinh u)/u², in the same scaled form, with series below u = 10⁻³). The extra inward acceleration is g_mem = −dC/dr, the **equivalent Newtonian mass** is M_eff(r) = r²g_mem/G, and the **equivalent density** is ρ_eff = −∇²C/(4πG).

**A theorem, to be verified numerically rather than assumed.** ∫∇²K_w d³x = 0 for the Gaussian, so ∫ρ_eff d³x = 0 and M_eff(r) → 0 outside the source: **the written field of a bounded source carries no net equivalent mass.** It redistributes force; it has no far field. Every cluster reading below is read in the light of this.

## 2. Light, declared

Two rules, both declared and neither derived, bracketing what a relativistic completion could say:

- **L0 — light sees ordinary matter only.** Deflection from Φ_N of the baryons alone: α̂ = (4/c²)∫∇⊥Φ_N dl. This is the case in which the written field acts on massive bodies and not on the light cone, as the conformal scalar of [gravity-response](../gravity-response/motion-and-lensing.md) does.
- **L1 — light sees the total potential with the factor two**, Φ = Ψ = Φ_N − C, exactly RPG-1's R4. Under L1 the memory field lenses as its equivalent density ρ_eff, including the negative shell.

Under L1 the memory part of the deflection has a closed two-dimensional form that needs no deprojection. Because ∫∂²_z C dz = 0, the projected equivalent density is Σ_eff = −(√(2π)·w·Λ/4πG)·∇²_2D F with F = Σ_b ∗ K_2D, the projected baryons convolved with the two-dimensional Gaussian — whose azimuthal average is **RUT-1's ring kernel** exp[−(R − R′)²/2w²]·I₀e(RR′/w²). Then

    M_proj,eff(<R) = −√(2π)·w·R·F′(R)/(2G),     α̂_mem(b) = −(2√(2π)·w·Λ/c²)·F′(b),

with F′ and F″ from the ring kernel's analytic derivatives (`equilibrium.py` already differentiates I₀e this way). Tangential shear is γ_t = [M_proj(<R)/πR² − Σ(R)]/Σ_crit with M_proj and Σ the sums of the baryonic and equivalent parts. This route shares no kernel with the spherical one, and their agreement is a gate.

**Geometry.** The six SLACS lenses use RPG-1's inputs unchanged: PF-1's static Euclidean distances from `lensing-data-readiness/conditional-geometry.json`, the archived Sérsic light profiles deprojected as spheres for the three-dimensional route and their analytic projected light for the two-dimensional one, and the population stellar masses (Chabrier and Salpeter) of the flux branch RPG-1 used. Coma's six Kubo bins are used as radius ratios with one adopted conversion, h = 0.7 — the same H₀ = 70 the X-COP masses carry — because a physical baryon model needs physical radii; the amplitude Σ_crit is a nuisance, so the Coma test is a **shape** test exactly as the [Coma inverse](../companion-extensions/coma-inverse-protocol.md) was. All of this is conditional model geometry, not measured distance, and the report says so.

## 3. Inputs, fixed

- **X-COP, 13 clusters** — [Ettori et al. 2019](https://arxiv.org/abs/1805.00035), Table 1, transcribed into `cl1-inputs-xcop.json` from the text already in `companion_wave_test/data/xcop.txt`: z, c₂₀₀, the backward-method NFW masses at 0.5, 1 and 1.5 Mpc, R₅₀₀, M₅₀₀, R₂₀₀, M₂₀₀ with their published statistical errors. [Eckert et al. 2019](https://arxiv.org/abs/1805.00034), Table 2 (`xcop_gas.txt`): M_HSE and f_gas at R₅₀₀ and R₂₀₀. **Status:** the hydrostatic masses are model-inferred (hydrostatic equilibrium, an NFW form, H₀ = 70, Ω_m = 0.3 distances), adopted as the *required* total acceleration g_req = GM(<r)/r² at the five published radii with their published errors; the five values of one cluster come from one two-parameter fit and are not independent, and the report says so. HydraA has no gas-fraction row and is excluded from every gas-dependent reading (12 clusters remain). A2029's M₅₀₀ differs between the two tables (8.82 against 8.65 × 10¹⁴ M☉); the required profile takes Ettori's, the gas mass takes Eckert's own f_gas·M_HSE, and the 2% difference is recorded. The non-thermal-pressure correction (Eckert's M_tot) is not applied; it raises the required mass, median 3% at R₅₀₀ and 8% at R₂₀₀, never lowers it.
- **The cluster baryons.** Gas: a β-model ρ_g ∝ [1 + (r/r_c)²]^{−3β/2} with r_c = 0.15·R₅₀₀ fixed and β solved so that the model's M_gas(<R₂₀₀)/M_gas(<R₅₀₀) equals the measured f_gas,200·M_HSE,200/(f_gas,500·M_HSE,500), then normalised to M_gas(<R₅₀₀). Stars: 0.09 of the gas, distributed like it (Ettori's median M_star/M_gas at R₅₀₀; the range 0.07–0.12 is a declared sensitivity). Declared sensitivities, each run and reported, none selected by fit: r_c/R₅₀₀ ∈ {0.10, 0.25}, M_star/M_gas ∈ {0.07, 0.12}. A stellar component distributed like the gas has no central cusp; a BCG would write a stronger local field than this model can show, and that is a stated limitation, not a knob.
- **Coma** — the two ends of CF-1's declared bracket, copied rather than imported so that CF-1's heavy modules are not loaded: β = 0.75, r_c = 296 kpc, μ_e = 1.17, gas truncated at 3 Mpc, n_e0 = 2.5×10⁻³ and 4.5×10⁻³ cm⁻³, stars 0.5×10¹³ and 2×10¹³ M☉ distributed like the gas. The six Kubo tangential-shear bins with plotted errors from `cluster-observation-readiness/kubo-figure-data.json`.
- **The six SLACS lenses** — J0037-0942, J1112+0826, J1204+0358, J1402+6321, J1621+3931, J1630+4520, loaded through RPG-1's `lensing.py` unchanged.
- **Reference laws, for comparison only:** Newtonian baryons; PM-1's local law g = g_N + √(a\*g_N) at the repository's archived a\* = 8.563×10⁻¹¹ m/s² and at PM-1's own fitted 6.54×10⁻¹¹ m/s². In spherical symmetry this is exactly what PM-2A's Completions I and II give, so it is also their cluster prediction. Not refitted.

## 4. Numerical verification: the gates, each with the control it must reject

Every tolerance below was measured on this configuration in a prototype before declaration; the measured value is quoted so a later reader can see the margin.

| gate | requirement | measured before declaration | the control it must reject |
|---|---|---|---|
| **G1 the shell kernel** | closed form against direct angular quadrature of the footprint at seven declared (r, r′, w) triples, including r ≪ w and r ≫ r′: relative difference < 10⁻¹² | 1.5×10⁻¹⁴ | RUT-1's *planar* ring kernel substituted for the shell kernel — the dimensionality mistake — must fail |
| **G2 derivatives** | analytic dC/dr and ∇²C against central differences on a Plummer sphere (M = 10¹² M☉, a = 5 kpc) at five radii and w ∈ {0.5, 2, 8} kpc: relative < 10⁻⁸ | 1.2×10⁻⁹ | the local-limit force −(2π)^{3/2}w³Λρ′ in place of the full gradient must fail at every width |
| **G3 the link to RUT-1** | (a) the ring kernel used here equals `rut1.phi_ring` at `rut1.gate_kernel`'s five cases to 10⁻¹⁴; (b) `rut1.gate_kernel()` and `rut1.gate_mature_ring()` rerun through import and pass; (c) the mature-ring control's w/ℓ, computed from `rut1.ring_depth`, equals 0.501 to three figures; (d) the ring kernel's analytic second derivative against a central difference of its analytic first derivative: relative < 10⁻⁸ | (a) identical formula; (c) 0.5013 | (d) with the cross term 2E′(I₀e)′ dropped must fail |
| **G4 no net equivalent mass** | for the Plummer sphere truncated at 60 kpc, M_eff(r_t + 8w)/(ΛM/G) < 10⁻⁸ at the three widths, and ∫₀^{r_t+8w} ρ_eff·4πr² dr < 10⁻⁶·ΛM/G | to be measured in the run; the estimate is 10⁻¹¹ | the Newtonian equivalent mass of the same source at the same radius, which is M, must fail |
| **G5 the two lensing routes** | the three-dimensional deflection integral of g_mem (shell kernel, adaptive quadrature) against the two-dimensional ring-kernel projected mass, on the Plummer sphere at b ∈ {1, 4, 10, 25} kpc and w ∈ {1, 3, 10} kpc: relative < 10⁻⁶, and the difference must fall when both source grids are doubled | 3.6×10⁻⁷ worst | the two-dimensional route with the shell kernel in place of the ring kernel must fail |
| **G6 the lens anchor** | with Λ = 0, the six Einstein radii for both IMFs through the three-dimensional route reproduce RPG-1's archived `newtonian_baryons` θ_E to 10⁻¹⁰; the two-dimensional route from the analytic projected Sérsic light reproduces them to 10⁻⁵ (it is limited by the archived deprojection's own quadrature) | 2×10⁻¹⁵; 8×10⁻⁷ | Λ = 10⁻⁶ (km/s)²/M☉ at w = 3 kpc must move every θ_E by more than 10⁻³ |
| **G7 X-COP ingestion** | NFW(c₂₀₀, R₂₀₀, M₂₀₀) reproduces the tabulated M(<0.5), M(<1), M(<1.5 Mpc) and M₅₀₀ to 0.5% for all 13, and R₅₀₀, R₂₀₀ from M₅₀₀, M₂₀₀ with the paper's H₀ = 70, Ω_m = 0.3 to 0.1%; the transcribed f_gas,500 agree exactly with the earlier hand transcription in `companion_wave_test/cluster_test.py` for the 11 clusters it holds | 0.40%; 0.04% | the same check with A85's and A2255's concentrations swapped must fail |
| **G8 Coma ingestion** | Σ(γ_t/σ)² over the six bins equals the paper's 23.33 to 0.1% | 23.336 | — (an ingestion check, as recorded in `kubo-figure-report.md`) |
| **G9 convergence of the exploratory quantities** | the cluster χ²(w) curve and Λ\*(w) unchanged to 10⁻⁵ relative when the radial grid is doubled (3,000 → 6,000 points); each lens's Λ(w) unchanged to 10⁻⁴ when its projected grid is doubled | 1.5×10⁻¹⁰ (clusters); lenses to be measured | — |

**Status rule.** `cl1_checks.py` is the suite job; it reruns G1–G9 and a regression anchor against `cl1-results.json`, and exits non-zero if any gate fails or any anchored number moves. Reproduction, numerical verification and the scientific reading are three statuses and are reported separately. A failed gate blocks every reading that depends on it.

## 5. The exploratory scan, declared

Nothing here is a gate. The parameter grid and every reading are fixed now; what the numbers turn out to be is the result.

**E1 — clusters.** Widths w on 41 logarithmic points from 10 kpc to 10 Mpc. At each w the amplitude Λ\*(w) is the weighted linear least-squares solution of g_N + Λ·s_w(r) = g_req over all 60 (cluster, radius) points, with s_w = −d(ρ ∗ k_w)/dr per unit Λ and weights from the published mass errors — linear because C is linear in Λ, so there is no optimiser to converge. Readings: χ²(w) and the log-rms in dex; the best w, its Λ, ℓ and w/ℓ; the same for each sensitivity; the per-cluster, per-radius ratio (model/required); the required/baryonic mass factor; and the three reference laws' χ². The prototype's values, which the run is expected to reproduce and which are quoted so no one mistakes them for a discovery of the run: best w ≈ 480 kpc, Λ ≈ 8.6×10⁻⁸, w/ℓ ≈ 9.6, χ² ≈ 4,030 for 60 points against 83,250 (baryons) and 24,450 (PM-1's law), with the model over the requirement at 1 Mpc and at half of it at R₂₀₀.

**E2 — lenses, under L1.** Widths w on 31 logarithmic points from 0.1 kpc to 10 Mpc. For each lens, IMF and w, the Λ_i(w) that gives the observed θ_E, by bracketing and Brent's method on the two-dimensional route (θ_E is monotone in Λ, checked at ½Λ and 2Λ). Readings: Λ_i(w), w/ℓ and w³/ℓ; the spread across the six lenses at each w (max/min); the width at which each lens's w/ℓ is smallest, in kpc and in units of its effective radius; and, under L0, the archived Newtonian ratios restated. The prototype: at w = 10 kpc the six lenses need w/ℓ between 5.1 and 9.0 (Chabrier).

**E3 — Coma, both rules.** For each bracket end: under L0 the shape χ² of ΔΣ_b against the six bins with Σ_crit fitted as one nuisance; under L1 at the cluster-best (Λ, w) of E1, the ratio ΔΣ_L1/ΔΣ_L0 at each bin and the shape χ² with the same single nuisance; and a free two-nuisance fit (amplitude and Λ) reported with the sign of its Λ. The reduced-shear-versus-shear difference in units of the plotted error is reported as the Coma inverse did. Reference: the inverse fits' 3.855 (NFW) and 3.727 (Plummer), which have more freedom and are not a benchmark to beat.

**E4 — the compatibility map.** A table over w of Λ_cluster(w) and the six Λ_i(w), their ratios, and the cluster-best pair applied to the lenses (θ_E ratios under L1) and the lens-preferred pairs applied to the clusters (χ²). The declared reading is a *statement of the region*: for which w, if any, do the clusters and all six lenses lie within a factor of two of one Λ. A factor of two is the exploratory criterion for "the same law could serve both", chosen now. RUT-1's supported populations (w/ℓ ≈ 0.5 at w/R = 0.1–0.2) are placed on the same map.

## 6. What would make this informative either way

If a common (Λ, w) served clusters and lenses within the criterion, the written-track steady state would have a domain worth taking to the SPARC and Milky Way comparison of PM-2A stages C and D with its own predicted fields. If the regions do not meet, the map says by how much and in which direction, and whether the obstruction is the amplitude, the width or the no-far-field theorem — which is the information goal 3 needs about what sets the width, and what goal 8 asked for. Either way the cluster residual pattern is a shape statement about the steady state that no formation run could change.

## 7. What CL-1 does not do

It does not fit galaxies or the Milky Way; it does not use the lenses' stellar kinematics (the joint motion-and-lensing test of CR-2's machinery is stage 2 if this stage earns it); it does not derive either light rule; it does not model a central galaxy in the clusters, non-thermal pressure, or the source geometry of the Coma shear; it does not evolve anything in time; and it opens no unexposed sample. Its exploratory readings are exploratory.

## Files

`steady_field.py` is the library (kernels, sources, deflection, shear); `cl1-inputs-xcop.json` the transcribed tables with provenance; `cl1.py` the driver, which runs the gates, then the scan, and writes `cl1-results.json` and compares it with the archive; `cl1_checks.py` the suite job, registered before the stage 8 job so the suite still ends on the job the owner is to rule on. Nothing here writes under another experiment's directory or edits any file of RUT-1's.
