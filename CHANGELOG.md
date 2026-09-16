# Research changelog

Dated notices that accumulated at the top of `README.md` and `research_plan/START-HERE.md` until 13 September 2026, moved here verbatim (newest first) so that [CURRENT-STATUS.md](CURRENT-STATUS.md) can remain the single authoritative summary. The date and commit are those that last touched each notice (`git blame`). Relative links were rebased to the repository root, and mis-encoded dashes in the START-HERE notices were repaired. Statements such as "latest" or "current" refer to the date of each entry, not to the present state.

## 2026-09-16

New entries are added here, newest first; [CURRENT-STATUS.md](CURRENT-STATUS.md) remains the authoritative summary.

- **PM-2A stage B: four field equations verified on a sphere, and the readout found to be the accuracy bottleneck** (protocol fa8b407). Plus the owner's two source-geometry corrections, neither of which changes a fitted parameter or the outer continuation.
  - **Correction 8.** `m_cyl` added a spherical bulge mass to a cylindrical integral. The shell projection M_cyl(R) = M_sph(R) + ∫M′(r)[1−√(1−R²/r²)]dr reproduces the Plummer closed form to 1.1×10⁻¹⁶; substituting M_sph understates the cylinder by 29.2893% at R = b and 55.3% at R = b/2. Impact on the published numbers: none, because every bulge's tabulated support ends at the outermost rotmod radius; the maximum change at any inner radius is 0.120% (UGC11914). It also replaces a false sentence: in spherical symmetry M_force equals M_sph, not M_cyl.
  - **Correction 9.** The thick-disk builder returns Σ(0)πr_in² — a cylinder's mass — for the solver's inner *sphere*, 15.57–15.60× the exact integral 2π∫R·Σ(R)[1−e^(−√(r_in²−R²)/h)]dR against a leading 3h/(2r_in) = 15. That inner region is 10⁻⁵ of the total and yet accounts for most of the published cell-total discrepancy, which falls from 3.9×10⁻⁵ to 5.3×10⁻⁶ once corrected — so that discrepancy was never mostly quadrature and truncation.
  - **Stage B's four equations** each carry their own μ, spherical inverse, relaxation weight and field functional, subclassing RPG-1's solver rather than editing it. μ(x)·x = y to 8.9×10⁻¹⁶ over 10⁻¹⁰ ≤ y ≤ 10¹⁰; every weight matches a numerical dlnμ/dlnx to 3.7×10⁻¹⁰; every F′ matches 2xμ(x) to 8.5×10⁻⁹; and the completions' functional densities agree identically on the matched spherical solution, F_I(y+√y) = y² + (2/3)y^{3/2}, to 1.3×10⁻¹⁵.
  - **Cancellation is real and now avoided.** The inherited simple-μ functional returns exactly 0.0 at x = 10⁻⁸ against a true 6.667×10⁻²⁵, and is wrong by 7.7×10⁻⁵ relative at x = 10⁻⁶; the series branch returns it correctly. Completion I's t³(t+2/3) form subtracts no nearly equal numbers at all.
  - **Solved on a Plummer sphere** at compactness 10⁻², 1 and 10²: the interior flux identity r²μ(g/a\*)g = GM(<r) holds to 4.4×10⁻⁷ with imposed faces excluded, the analytic law to 2.3×10⁻⁷, and the two completions agree with each other to 2.3×10⁻⁷ — as they must on a sphere. At g_N = a\* they give 2.000000 against simple-μ's 1.618034. Convergence is gated on the equation residual (1.8×10⁻¹²) as well as the iteration change (3.4×10⁻¹²); `aqual.py` computes a residual but sets its flag from the change alone.
  - **The readout is the bottleneck, by three to four orders.** `gradient_at` (cell-centre averages) and `midplane_speed` (linear in ln r) both carry 1.6×10⁻⁴ at nr = 300 and converge at second order, on a solution whose face gradients are accurate to 10⁻⁷. A cubic spline of the face gradients gives 1.3×10⁻⁸ and converges at about fourth order. Only the spline is gated; the inherited numbers are published as a measurement. RPG-1's curves are unaffected in any observationally meaningful sense, but differencing two completions is a different demand.
  - **Stage B is an equation implementation verified on a sphere.** It is not path memory, and the solver's relaxation time is never a physical memory time.

- **PM-2A stage A: the source passes its audit, and PM-1's mass is not the source's mass** (protocol d28b67a). One frozen ordinary-matter source per SPARC galaxy, built from the archived construction and integrated in closed form, gated before any field is solved.
  - **The outer source follows the owner's ruling.** The archived endpoint-anchored exponential continuation is primary, removal of the extrapolated stellar tail is the declared sensitivity, the gas and bulge prescriptions are unchanged, and neither prescription is selected by which fits better. The extrapolated tail is a median 1.06% of the total source mass, 8.0% at the 90th percentile and at most 28.1% (UGC06667); of the stellar mass alone, 2.77% / 15.9% / 57.1%.
  - **All four gates pass.** Refinement conserves the source mass to 2.2×10⁻¹⁶ with no clipping; the predicted force is invariant to velocity-row removal at 0.0 across the 131 galaxies with at least eight rows, every source hash unchanged; the closed forms agree with an independently converged Gauss–Legendre integration to 6.6×10⁻¹³, that quadrature's own 8-to-16-node refinement moving by 1.8×10⁻¹⁵; and the collective rule's 1/√(1−f) tail algebra holds to 2.2×10⁻¹⁶.
  - **The invariance gate carries a positive control,** because it is passed trivially by any pipeline that never consults the mask. Rebuilding the source from a shortened rotmod array — the mistake the gate exists to catch, which moves the stellar profile, the gas scale and the bulge at once, one galaxy's total mass by 8.6% — changes the hash every time and moves predictions by 1.4×10⁻⁴ to 2.3×10⁻².
  - **M_force exceeds the source's cylindrical mass at every galaxy:** median 1.381 at the last sampled radius, range 1.040–3.174, and 1.177 against the total. That is the expected sign for a flattened source, and the measured size of the gap stages B–D have to account for by solving a field rather than substituting a mass. The ratio mixes disk flattening with a second cause, in unknown proportions: `v_bar` is SPARC's tabulated contribution while M_cyl is the source reconstructed here, and RPG-1 documents that the two differ materially at some radii. Separating them needs the source's own Newtonian field on the same grid, which is stage C's recomputed baseline. M_sph is deliberately not quoted: it is not defined by a surface density alone.
  - **The archived quadrature is measured, not trusted.** RPG-1's `_cylinder_mass` differs from the exact totals by a median 2.4×10⁻⁷ (max 4.5×10⁻⁷) across all 149 galaxies, so it is kept as a regression diagnostic; two quantities produced by that route agreeing cannot establish 10⁻⁸ accuracy. The AQUAL cell construction's total differs by 0.9–3.9×10⁻⁵ and is reported rather than gated, because it decomposes a thickened disk over spherical shells. Gate 3's declared wording, which applied 10⁻⁸ to that comparison, was corrected in the protocol before the run rather than quietly relaxed. RPG-1's published results are untouched.
  - **Stage A solves nothing and scores nothing.** Passing says the source is well defined, not that it is right.

- **PM-1: a broadband path-memory force reaches MOND's fit and MOND's acceleration scale, and its feedback version is unstable** (protocol dcddb70). Six candidate families for "gravity remembers where matter has been", scored on SPARC's 149 galaxies and 3,152 radii and the Milky Way's 38 Eilers bins, one universal constant each, fitted on the 89 training galaxies only.
  - **Delay alone explains nothing.** An exponential memory of density returns αΦ_N for a stationary source, and a rigidly rotating axisymmetric density differs from a static one by 3.7×10⁻¹¹.
  - **One wavelength cannot work.** k J₀(kR)J₁(kR) changes sign 35 times across UGC02953's 115 radii.
  - **A broad spectrum can.** Equal potential weight per logarithmic wavenumber gives (C/2R)[J₀(k_min R)² − J₀(k_max R)²], so rings sum and the force at r comes from the matter inside r.
  - **The sourcing rule is what the data decide.** Amplitude ∝ M: 30.53/34.03/29.46 and a mass–speed slope of 0.412 against the observed 0.289, which fails. Amplitude ∝ √M: 20.21/27.74/18.07 and a slope of 0.280. One collective saturated mode, weight ∝ mass and no extra parameter: 21.66/29.67/18.71 and 0.2886, the closest slope of the three.
  - **That is MOND's fit, not better than it** (simple MOND 19.89/26.88/16.40; Milky Way 9.53/12.04 against 18.25/9.51), and 1.6%, 3.2% and 10.2% higher on the three splits with no paired uncertainty computed.
  - **The acceleration scale is the fitted amplitude re-expressed,** not a second prediction: g_mem = √(a\*·g_mono) identically with a\* = β²/G = 6.54×10⁻¹¹ m/s², and v⁴ ∝ M was built into the √M choice.
  - **The mass is force-equivalent, so the rule is a local acceleration law.** M = R·v_bar²/G = R²g_N/G, and substituting it gives identically g = g_N + √(a\*·g_N) for the disk inputs too (4×10⁻¹⁶ across 149 galaxies). The scores stand as an empirical acceleration law; the mass–speed regression is against force-equivalent mass, so the baryonic mass–speed test is still owed; the collective control's total depends on where each curve stops; and the √8 result is a mathematical warning rather than a physical annulus test. The archive is preserved under that label.
  - **The mechanism does not produce the fitted law.** Independent per-ring saturation gives Σ√mᵢ, not √(Σmᵢ): refining one galaxy 8× multiplies that field by 2.8284, exactly √8. The square root is derivable from a collective nonlinear field instead — the spherical limit of ∇·[(|∇ψ|/a\*)∇ψ] = 4πGρ_b is √(Ga\*M(<r))/r, matching the fitted amplitude identically.
  - **An audit withdrew one row.** The finite-band force had been computed with the on-ring expression at the observation radius, which moves the source: for a ring at R_s = 1 observed at r = 2 the source-integrated value is 0.49877 against the on-ring 0.24989. The wide-band enclosed-mass rule is confirmed as the correct limit (the ring integral is 1/r outside and ~0 inside).
  - **The feedback version is unstable, analytically.** Routh–Hurwitz on the owner's cubic requires q₀ < 0, so every q₀ > 0 at finite T > 0 has a growing mode; the root reproduces to five figures (0.15334 ± 1.28578i at q₀ = 0.9, T = 1). At real SPARC outer radii the median q₀ is 0.80, an e-folding of about one orbit, though the growth rate vanishes as T → 0 and T → ∞. It does not transfer to the force law: in the frozen potential −GM/r + K ln r, κ² = GM/r³ + 2K/r² > 0.
  - **Energy deposition is 1.7×10⁹ short** in the median galaxy, over a range of 1.4×10⁸ to 1.2×10¹¹.
  - Candidates D and F were declared and not run: D needs a thickness of tens of kiloparsecs that SPARC does not measure, and F needs a trajectory-history integrator. `path-memory/checks.py` becomes suite job 64.
- **BRIDGE-1B stage 2: a saturating mass removes the standing drain, and the brightness threshold with it** (protocol and corrections in 8119bcc). The companions' mass law becomes m_C² = m0² + M²tanh²[(n − n\*)/Δn] with g = M/Δn, so the crossing is stage 1's to leading order and the mass stops growing afterwards.
  - **The production event survives only if the plateau is far from the crossing.** The occupations match the crossing formula to 0.22% at M/k\* = 20 and 0.16% at 100, but to 3.07% at M/k\* = 5, which fails the declared 1%. The owner's independent four-mode test found 0.20% and 0.18% at M = 20k\*.
  - **Saturation ends the back-reaction.** With the coupling that traps stage 1's field at Δn = 124, the saturating field rolls to 568 and is still at 0.96 of its speed; beyond 5Δn the force is 4.6×10⁻⁴ of its peak, occupations hold to 3.6×10⁻⁴, and the energy delivered equals N(m_final − m0) to 0.24%.
  - **The brightness threshold becomes a condition on the epoch.** The drain is a pulse, so the history coasts on both sides of a step in ṅ of √(1 + 1/R_E) − 1. The Pantheon+ χ² is 871.6, indistinguishable from PF-1's coasting fit, for every crossing at z\* ≥ 3.2 at every stored energy tested including R_E = 1; z\* = 1 needs R_E = 30, and later crossings never qualify. Stage 1's "145 times the companions' rest energy" was a property of indefinite growth.
  - **The energy is still the field's, now as a pulse:** the companions take their rest energy in 1.78 hours, 1.4×10⁻¹¹ W/m³ against the photons' 1.0×10⁻³¹. Coldness is unaffected (2.1×10⁻⁹ km/s, 2×10⁻⁵ kpc travelled in 9.8 Gyr, now at constant speed). The abundance is scanned rather than fixed to 2B-F1's rate.
  - **An independent cross-check.** The owner's illustrative 1 GeV mapping, recomputed here: k\* = 1.0133×10⁻⁵ eV against their 1.01×10⁻⁵, 2.099×10⁻⁹ km/s against 2.10×10⁻⁹, 1.781 hours against 1.8.
  - **Two corrections.** The driven variant's definition had inverted a sign (â is the drain minus the drive, as coded). P2 as declared asked for energy conservation in every trial, but P1's trials prescribe the field, which is not a closed system; it is now P2a (closed runs, 3.2×10⁻¹¹) and P2b (the prescribed runs' gain against an independent quadrature of the driver's work, 3.7×10⁻⁹). Stage 1's V1 has the same non-closure, 0.36, which stage 1's V2 never covered.
  - The atomic-reference bound p ≤ 6.6×10⁻⁷ is unchanged: it depends on the coupling family, not the mass law.
- **BRIDGE-1B: the field that reddens light can make the companions, but it must already carry their energy** (protocol 31d143e). One interaction, calculated in a homogeneous closed volume: radiation gives energy to the matched-wave propagation field, the field produces companion excitations through m_C²(n) = m0² + g²(n − n\*)², and those react back on it. One accounting of the photon's loss.
  - **The source is no longer a set of knobs.** Production is one non-adiabatic crossing of the mass minimum, with the spectrum exp[−π(k² + m0²)/(g|ṅ\*|)] (verified against the exact mode integration to 3×10⁻⁴). Afterwards momenta are frozen and the mass grows as g(n − n\*), so the mass production rate is q = N g ṅ, tied to the redshift rate, and the companion density follows the redshift history, ρ_C ∝ n − n\*.
  - **Coldness is free.** The products' r.m.s. speed is 2.3×10⁻¹⁹ km/s a Gyr after the crossing; a companion travels 89 m before it is slower than 3 km/s, and 427 m in all. Mass growth with conserved momentum cools without expansion; the prototype's 0.054c is a code-unit artefact.
  - **The energy is the field's, not the light's.** Radiation delivers 1.5×10⁻⁶ of the companions' rest energy between the crossing and today (the owner's power ratio of 2.75 million, reproduced). It delivers everything only if the history begins at rest at a turnaround where the radiation holds it all: 6.3×10⁷ K at R_E = 10, 8.4×10⁸ K at the brightness threshold.
  - **Draining the field brightens distant sources.** The Pantheon+ χ² stays within 1 of PF-1's coasting fit only if the field carries at least 145 times the companions' rest energy: 3.0×10⁸ times the microwave background's energy, 1.4×10⁴ times the critical density. The owner's damping approximation, scored on the same 960 rows, gives 896.7, 1045.1 and 1298.3 against PF-1's 871.6.
  - **A driving potential moves it the right way.** A linear potential supplying 4–6 times the companions' drain improves the fit to 840.4 (flat FLRW 836.5 at one parameter). That is a two-parameter fit to exposed data, and the potential's energy is stored, not light.
  - **The atomic-reference factor, as a bound.** In the minimal coupling family, fixed rulers and a fixed fine-structure constant force the hyperfine-to-optical ratio to drift at exactly 2p·ṅ/n, where p is the surviving fraction of the coordinate shift. Yb⁺/Cs then gives p ≤ 6.6×10⁻⁷, so a homogeneous propagation field cannot supply the measured redshift in that family. None of the five tested completions passes the four optical gates.
  - **Handed to RC-2a:** one production epoch, the momentum spectrum, a mass that keeps growing (a drag of 0.1 per Gyr today) and the donor depletion that is the field's own energy loss. `shared-field-bridge/checks.py` becomes suite job 63.
- **Stage 2A's archive regenerated after the seed-order fix** (a53f45c). Stage 2A's verification rounds 2–4 had drawn their seeds in the order runs finished, so they depended on timing. The canonical run was repeated from scratch with the fix (8,462 s on 7 workers). Rounds 0 and 1 reproduced the archived run line for line, every validation still passes, and no conclusion changes. What moved:
  - **Shape.** The Milky Way's companion mass across 5–25 kpc grows as r^1.7–2.7 (was r^1.9–2.7), and across every system and tested radius as r^1.5–2.7 (was r^1.7–2.7).
  - **The bath's gravity.** With it on, B1 needs 0.16–0.67 of the gravity-omitted density (was 0.16–0.71), and the companions inside r_half are 0.02–0.9% of the baryons there (was 0.03–0.9%).
  - **B2.** The total confined mass is 490–670 times the baryons in the Milky Way (was 410–670). J1630's lower end, 87, had been printed as 88.
  - **Controls.** The counted seed changes M(T) by −1% to +5% (was +8%). A frozen potential forms 0.59–0.98 of the coupled mass (was 0.58–0.97). Ejections remove 0.01–0.75 of the confined mass in the galaxies (was 0.72).
  - **Labels.** One more combination is transparent (26 of 47). Ten runs rather than eight reached the collision cap, all at 100–1,000 cm²/g, outside the regime.
  - The current documents and 2B-F1's report now cite the new ranges; declared protocols keep the numbers they were written with.
- **RC-2a stage 1: the 300 kpc zone had chosen 2B-F1's fit** (protocol 1585c9f, Amendment 1 in 2a25cb3). Cold collisionless decay at 3 km/s, at 2B-F1's rate, with the donor counted and every daughter followed in an open region.
  - **Only at 300 kpc.** The fit exists in the 300 kpc region (RMSE 7.9 km/s), and there it matches 2B-F1's static-bath result. From 600 kpc outward the inner rotation overshoots by about 170 km/s.
  - **Convergence to the wrong galaxy.** The apertures converge from 600 kpc (statistically from 1,200) to 8×10¹² M☉ inside 200 kpc, with circular speeds of 410–446 km/s.
  - **Where the mass came from.** The companions inside 25 kpc were born at 140–570 kpc, so the measured infall horizon is 570 kpc. The 300 kpc runs' supply stopped at their edge.
  - **A transient.** Continuous, stopped and matched-budget production each cross near the observed profile, then overshoot.
  - **Controls.** The quiet null passes. Random sampling without baryons seeds central concentrations of up to 7.6×10¹⁰ M☉ at 16,000 births; they vanish at 64,000. The two resolutions agree without being resolved, so the declared gate closed and the search did not run.
  - **Validation passes.** `donor-companion/checks.py` becomes suite job 62.

  See the [report](research_work/results/donor-companion/report.md).
- **The 2B-F1 revision's results** (protocol 9232e07, with Amendment 1).
  - **The search.** In 22 of 26 combinations, the direct search (three seeds per rate) selects the rate 2B-F1 sampled best. The minimum is sharp: one step of 3^(1/8) either side raises the RMSE by about 4 km/s.
  - **Four robust passes.** 3 km/s collisionless, and 10 km/s at 0, 0.1 and 1 cm²/g, all with the bath's gravity on. Each passes on its seed means and again at four times the tracers, with RMSE 7.0–9.2 km/s, slope 1.07–1.26 and net cost 1.0–2.7.
  - **Two marginal passes fail.** At 3 km/s with 0.1 and 1 cm²/g, the seed-mean slopes are 0.99 and 0.77. The first was 2B-F1's best source.
  - **Numerical controls.** Halving the timestep raises the reference case's slope from 1.20 to 1.32, while the grid and the pools are converged. So the slope carries a numerical uncertainty of about 0.1.
  - **F3,** uncapped (1,024,000 births, no thinning), agrees with the quadrature to within −3.4% to +4.7%, at standard errors of 1.8–2.6%.
  - **The declared best source** is now 3 km/s, collisionless, bath gravity on. J1630 completes and Coma keeps no static bath. As the owner's review asks, this rule does not choose the next mechanism.
  - **RC-2a's amended code is in place.** Stage 2A's engine gains optional per-tracer tags, bit-identical when unused, and RC-2a's canonical run is under way.

  See the [report](research_work/results/companion-source/report.md).
- **Amendment 1 to RC-2a,** after the owner's review of 1585c9f. The first canonical attempt was stopped during its first round, before any result was read.
  - **Scope.** One decay speed, 3 km/s, run in three stages:
    1. controls and the reference at the fixed rate;
    2. the direct search, only if the controls are interpretable;
    3. the region-size comparison repeated at the selected rate, frozen.
  - **Sampling noise.**
    - A quiet representation for the reference: each shell and step gets its exact converted mass, with stratified radii and paired radial velocities.
    - C1a, a numerical null without baryons.
    - C1b, the same with random sampling at two resolutions.
    - C1c, the reference at two resolutions.
    - Noise excursions are measured as root-mean-squares over seeds.
  - **Convergence of the gravitational signal.**
    - Tested at every aperture (circular speed and net contrast) and at all 38 radii (rotation speed).
    - Each quantity must both agree and be resolved: its standard error must be small enough to detect a tolerance-sized difference.
    - The positive inventory is reported but no longer decides.
    - R_conv needs a comparison with a larger completed region.
  - **Where the mass came from.** Every tracer carries its birth radius, birth time and angular momentum. The report gives birth distributions for the companions inside 10–200 kpc, and measures the infall horizon instead of estimating it.
  - **A matched-budget control, C3b:** twice the rate for 5 Gyr, then nothing.
  - **The engine.** Stage 2A's engine gains optional per-tracer tags. They are absent in every earlier run, which stays bit-identical.

  See the [protocol](research_work/results/donor-companion/protocol.md).
- **Amendment 1 to the 2B-F1 revision,** after the owner's review of 9232e07. The first canonical attempt was stopped after about 18 minutes, before any results.
  - **F3's tolerance.** Three standard errors of 4% is 12%, so the 5% bound does not bind there. It controls only below 1.67%. F3 stays a statistical check and is no longer called a 5% convergence test.
  - **F3's cap.** The helper passed no tracer cap, so the engine's default of 40,000 applied. 2B-F1's F3 thinned twice and ended with 21,169 tracers, with 32–66 effective samples per bin. The revision caps F3 at four times its births and reports thinnings and each bin's effective sample size.
  - **Numerical controls (R7).** Halved timestep, doubled grid, and doubled pools for the reference case. Round D's extra tracers test tracer convergence only.
  - **The best-source rule** rewards cancellation within a passing tier, so it no longer chooses the mechanism investigated next.

  See the [protocol](research_work/results/companion-source/revision-protocol.md).
- **RC-2a declared: a donor-and-companion calculation in an open region,** the owner's highest-priority physical experiment after the reviews of ed96b00 and 9232e07.
  - **The question.** With the source law and its parameters held fixed, do the masses and velocities at fixed observational radii converge as the computational region grows, or does the boundary choose them?
  - **The model.** Cold collisionless decay at 3 km/s is the reference case, with 10 km/s as the comparison. The static bath and its subtraction of the incident density are replaced by daughters followed in time, bound and unbound, minus the donor's depletion. So the donor is counted once, and no second −qt is added.
    - Regions of 150–2,400 kpc run at 2B-F1's best sampled rates for the collisionless, bath-gravity-on combinations, never refitted per region.
  - **Scoring.** Seed means must agree from one doubling to the next: the 38-bin RMSE within 1 km/s, the 8–20 kpc slope within 0.05, and the companion mass inside each outer aperture within 5% (or two standard errors, whichever is larger).
    - A direct search at 1,200 kpc gives the normalization.
    - J1630 and Coma follow as predictions.
    - Inventory, mass contrast and source energy are reported separately.
  - **Controls.**
    - C1: the same source without baryons.
    - C2: births against depletion at every aperture.
    - C3: the source switched off after 5 Gyr.

  See the [protocol](research_work/results/donor-companion/protocol.md).
- **The 2B-F1 revision's code, stage 2A's seed order, and the F1 archive's profile radii** (a53f45c, and the archive refresh after it).
  - **The revision.** `companion-source/revision.py` implements the revision declared in 9232e07. Its first canonical run was stopped after the owner's review of 9232e07 found that F3's helper passes no tracer cap, so the engine's default of 40,000 thinned the 1,024,000 tracers requested. It restarts after an amendment.
  - **Stage 2A's seed order.** Stage 2A's `run_tasks` now returns results in the order the tasks were issued, so its rounds 2–4 no longer draw seeds in completion order. Three smoke runs on 4, 4 and 3 workers agree bit for bit. Its canonical rerun is under way.
  - **Profile radii.** Each profile row now stores its shell's lower edge, outer edge and centre, with the mean radial velocity and the radial dispersion about it. The regenerated F1 archive (37.7 min on 16 workers) matches ed96b00 value for value apart from the profile keys, and every row's centre equals its old radius exactly.
- **Corrections after the owner's review of ed96b00.**
  - **The outer profile.** Each profile row stores its shell's centre as the radius but its enclosed masses at the outer edge, 1.1926 times farther out. The first report paired them, so its outer circular speeds were 9.2% too high. The table and prose now use the outer edge: 166–195 km/s out to about 110 kpc, and 0.8–1.0×10¹² M☉ inside 148 kpc. The fit scores are unaffected.
  - **The cost gate.** With the bath's gravity on, the subtraction of the uniform density is, for a homogeneous donor, the donor's own depletion, so the G3 net is a mass contrast. Positive inventory, mass contrast and source energy are kept separate; the verdict does not make the source budget small.
  - **What is robust.** The inner shape holds without the bath's gravity and without scattering: 8.0 km/s and slope 1.32 at 3 km/s, collisionless.
  - **Narrowed claims.** The normalization search is not exhaustive, since the rescaling root is not in general the RMSE minimum; at 150 kpc no sampled rate passed. The Jeans number is a growth-time diagnostic. F3's bins differ from its quadrature by up to 21%, inside their 12–20% standard errors: a statistical check, not a 5% convergence test.
  - **Stage 2A's driver** draws the seeds of its rounds 2–4 in task completion order, as F1's did, so its archived canonical is not bit-reproducible. The fix and a rerun are queued with F1's revision.
  - **The queue** follows the review: F1's revision (a direct search around the best actual RMSE, seed repeats at higher resolution, F3 with more tracers); RC-2 as a time-dependent donor-and-companion calculation whose apertures converge as the region grows; an extended, baryon-triggered source that counts its donor; three-dimensional stability and a common field frame; finite source histories before Bose enhancement.

  See the [report](research_work/results/companion-source/report.md).
- **CC-2 stage 2B-F1, a decaying-field source** (protocol 2cf7ea1), requested in the owner's review of 68eb17c.
  - **The shape.** Cold decays (3–10 km/s) with the bath's gravity on pass the Milky Way's profile gate in five of six combinations. The rotation RMSE is 6.8–8.7 km/s over the 38 Eilers bins, and the enclosed-mass slope across 8–20 kpc is 1.03–1.31 against the required 1.31, at q\* = 1,160–1,440 M☉ kpc⁻³ Gyr⁻¹. Companions born throughout the zone fall in on radial orbits, so M(<r) ∝ r/v_esc.
  - **The cost passes on the net.** The confined companions are 12.6–16.2 baryon masses. The cold bath barely enters the zone, and its shortfall against the uniform density, which is assigned to the background, brings the total to 0.75–2.6. With the bath's gravity omitted, no cold combination passes.
  - **Faster decays fail.** At 30 km/s the slope is 2.3–2.6, and at 100 km/s the cost is 23–56. At 300 km/s only an omitted-gravity run just below the runaway passes, at a million times the cosmic mean.
  - **Supply and universality.** The field must make 280–360 times the cosmic mean density in companions (1.7–2.1 million times the microwave background's energy), and its uniform medium is Jeans-unstable within about 1.2 Gyr. With the same source Coma keeps no static bath, and a 150 kpc zone fails the gate.
  - **Numerics.** The first canonical run failed validations F2 and F4 because of the checks, not the engine: a coarse quadrature grid, and a prediction drawn from one fresh set of pools. Both checks were corrected. The rerun then showed that the refinement rounds took their seeds in the order runs finished. The driver now fixes that order, and the canonical run was repeated.

  See the [report](research_work/results/companion-source/report.md).
- **RC-1 narrowed; stage 2B redirected** (owner's review of 68eb17c).
  - **RC-1 is a local stock bound.** Its energy comparison shows that today's local radiation cannot become today's local bath. The global source requirement needs a power and flux ledger, now queued as RC-2 for open, three-sphere and three-torus backgrounds.
  - **Stage 2B tests phase space too.** The elastic mechanism builds an outer envelope (M ∝ r^1.7–2.7, against the Milky Way's r^1.3), so more supply alone does not help. Stage 2B becomes 2B-F1: a coherent field decaying into companion pairs, propagated through stage 2A's machinery and scored on energy supply and on a hard profile gate.
  - **The recording-transition item** now asks first whether one transition can leave both a blackbody and persistent correlated modes. The unification problem stays explicit in the queue.
- **RC-1, the source-side radiation check** (protocol cc02334), before stage 2B selects an interaction.
  - **Energy.** At stage 2A's required densities, the companions carry 1.1×10⁵–2.5×10⁸ times the microwave background's energy per unit volume. Present radiation holds the requirement in none of the 33 in-regime combinations.
  - **FIRAS.** Within the diagnostic Δχ² = 4, a frequency shift at fixed photon number can take at most 3.3×10⁻⁵ of the background's energy, and whole-photon removal 8.4×10⁻⁵.
  - **The Planck-preserving law** is the only radiative route left. It needs a past background at 50–340 K, 38–63 Gyr of transfer at α₀, and the removal of all but 10⁻⁶–10⁻⁴ of its photons. Applied to every beam, it dims sources by 0.81 magnitudes at 1 Gpc, so it needs a coupling about 74 times weaker on starlight.

  See the [report](research_work/results/radiation-budget/report.md).
- **CC-2 stage 2A, complete-channel formation** (protocol 1ccd762), requested in the owner's review of 299ba1f.
  - **Supply.** All three collision classes act together, with mass, energy, orbits, depletion and gravity evolving. Within the model's regime, the incident density that forms B1 is 29–580 times below CF-1's seeded-only requirement in the Milky Way, 2.8–44 times in J1630 and 77–1,600 times in Coma (bath gravity omitted). It is still 19–42,000 times the cosmic mean.
  - **An outer envelope.** At B1 the companions inside r_half are 0.03–3% of the baryons there, and the enclosed mass grows as r^1.7–2.7 across the tested radii. The inner benchmark needs 2.2–3.9 times the density and a total confined mass 16–670 times the baryons.
  - **The bath's own gravity.** A stable static bath exists only below a limit density. With its gravity on, B1 lies below the limit in nine galaxy combinations, where it needs 0.16–0.71 of the gravity-omitted density. Coma keeps no static bath with its gravity on: the boundary model's own limit.
  - **Controls.** The result is seed-independent in all 33 in-regime combinations, and self-gravity matters in 29.
  - **Numerics.** A first complete canonical run was discarded because its seed controls were wrong. Two defects caused it: tracer splitting with random-removal thinning in dense cores, and production shells switched off for the whole run when their first pool held no event. Both were fixed and the canonical run was repeated. The report lists every deviation from the protocol.

  See the [report](research_work/results/companion-formation/report.md).
- **Shared constraints with the microwave background** (owner's direction after 2c9a110). The [persistent-background plan](research_plan/persistent-microwave-background.md) now states what the companion mechanism shares with it: the radiation budget and spectrum on the source side, a separately specified global background with an open-boundary, three-sphere and flat three-torus controls, the recording transition first, and one gravity for motions, galaxy lensing and CMB lensing. The queue gains a source-side radiation check before stage 2B, and those three later items.

## 2026-09-14

- **Stage-1 corrections after the owner's review of 299ba1f, and the stage-2A protocol.**
  - **No premise decision is needed.** The stage-1 report had asked whether a dominant, self-interacting companion medium fits the premises. The universe contract already prohibits an inserted halo or bath and permits a gravitating population derived from the permitted supply and interactions. CC-1's clock field stays a separately labeled comparison; the active objective is a source compatible with operational nonexpansion.
  - **The old density is not a lower bound.** Stage 1 found that seedless formation dominates at CF-1's seeded-only density but did not solve again for the density the complete mechanism needs. Seedless formation grows as the density squared, so the supply is undetermined until it is recomputed.
  - **Narrowed labels.** Jeans results are instability diagnostics, boundary-volume masses are environmental estimates, and uT is a transport benchmark. The velocity thresholds and the 0.34% fast fraction belong to CF-1's seeded coefficient and the tested distributions.
  - **Stage 2A** ([protocol](research_work/results/companion-formation/protocol.md)) recomputes the supply with all three collision classes together, evolving mass, energy, orbits, bath depletion and gravity on the Milky Way, J1630 and both Coma brackets, with zero-seed and counted-seed controls and no prescribed profile.

- **CC-2 stage 1, a regime and supply audit** (protocol 5381c68), requested in the owner's review of f13c09f. Is CF-1's incoming bath consistent at the supply CF-1 requires?
  - **No.** Across the Milky Way's and Coma's boundary regions its optical depth is 4–8. For every trial σ/m from 0.1 to 1,000 cm²/g, its own mass there exceeds the baryons by 10¹–10⁷, and in Coma it is Jeans-unstable at every trial scattering strength. It would be the dominant self-gravitating, self-interacting medium, not a trace bath.
  - **The seed stops mattering.** Under the same law, capture from collisions between two incoming companions beats seeded capture by 10²–10⁷, and it is depth-keyed.
  - **Retention windows.** Retention needs a Maxwellian bath dispersion below about 210 km/s (Milky Way) or 510–700 km/s (Coma). A 0.34% tail at 3,000 km/s reverses the Milky Way's.
  - **Coma's supply region.** A 300 km/s population reaches Coma from only about 3 Mpc in 10 Gyr.
  - The CF-1 report is qualified to match: its survival statements are initial statistics, and its runaway belongs to the fixed-shape rule.

  See the [report](research_work/results/companion-supply/report.md).
- **Focusing and capture CF-1** (protocol 012f6ca). A conditional transport-and-capture diagnostic, as the owner directed. One slow companion bath is placed in potentials built from ordinary matter plus a counted seed: the Milky Way, six SLACS lens hosts, and Coma over bracketed gas-and-star inputs, with a full NFW well only as a labeled inverse diagnostic.
  - **Transport verified.** Focusing raises the entry rate and the unbound density exactly as predicted, to 10⁻⁵ and 10⁻⁸. With capture off, nothing is retained.
  - **One interaction: elastic companion scattering.** Retention needs a bath slower than the host.
    - At 300 km/s, Coma retains 140–290 times faster than the Milky Way with favorable initial energy statistics.
    - The Milky Way faces a large heating diagnostic, and the lens hosts are marginal.
    - (Qualified after the owner's review of f13c09f. These are initial rate coefficients. The fixed-shape growth rule does not evolve the heating, so neither survival nor overheating is an evolved result.)
    - At the trial speed of 3,000 km/s every system erodes.
  - **Supply.** It needs 10⁴–10⁶ times the cosmic mean density at 1 cm²/g.
  - **Correction.** From baryons alone, Coma's escape speed is 1,100–1,500 km/s, not the 4,500 km/s of its inferred well that the earlier focusing estimates used.

  See the [report](research_work/results/gravitational-focusing/report.md).
- **CR-2 corrections after the owner's review of e5f6fc3.** The report had said the condensate's edge near 75 kpc costs it the Milky Way. Material beyond the measured 5–25 kpc cannot change those speeds, so that explanation was wrong.
  - A post-hoc interior diagnostic (`milky-way-interior.py`) finds the actual cause. The condensate is nearly uniform across the disk, so its enclosed mass grows as r^2.7 where the data need r^1.3. It is too thin inside about 18 kpc and too heavy beyond.
  - The one-scale NFW comparison is relabeled a restricted benchmark. NFW's own universal law lets the scale follow mass and history.
  - z ≈ 5,200 is relabeled a required-budget estimate.
  - The queue now starts with CF-1, a conditional transport-and-capture diagnostic, as the owner directed.
- **Supported reservoir CR-2** (protocol a5b2c41). A Thomas–Fermi condensate with one shared constant, in equilibrium with each lens's stars, is fitted jointly to the KCWI stellar motions and exact lensing of six SLACS lenses, with identical geometry and mass conventions.
  - **Result: it fails the declared rule in both geometries.**
    - The best shared core size is about 80 kpc in both.
    - χ² is 113.6 against free NFW's 85.3 (FLRW) and 89.3 against 59.1 (co-scaling), beyond the declared allowance of 10.
    - The same core size cuts the Milky Way rotation RMSE from 52.6 to 22.5 km/s, short of 20.
  - **Separate results.** Supply falls short by 10⁶–10⁹. Chabrier population masses cannot fit; the lenses need 1.5–3 times Chabrier.
  - **Validation.** The archive's benchmark fits are reproduced to 10⁻¹³. The co-scaling benchmarks are new.
  - **Post hoc, universal against universal** (the owner prefers universal settings to per-system tuning). At equal freedom, the condensate beats an NFW halo with one shared scale (a restricted benchmark) on the lenses (113.6 against 128.6; 89.3 against 105.7). It loses in the Milky Way (22.5 against 14.1 km/s).

  See the [report](research_work/results/supported-reservoir/report.md). Its leads section proposes a two-phase support law and a cosmological supply test.
- **Model contract revised after the owner's review.** Co-scaling is now a separately labeled branch, not a declaration that the nonexpanding requirement is met; the fixed-ruler branch stays open. The contract changes these claims:
  - "the only remaining alternative" becomes "one completion compatible with the kinematics of the tested family; other completions have not been exhaustively classified";
  - "one field does all three jobs" is marked as a goal, with no complete action yet;
  - the energy-density argument for baryon dominance is replaced by a requirement for coupling equations;
  - the claim that lensing and inner dynamics agree becomes a reason for a joint fit with identical conventions;
  - the claim of no Big Bang singularity is withdrawn.

  See [research_plan/model-contract.md](research_plan/model-contract.md).
- **Co-scaling completion CC-1** (separately labeled branch; protocol 6ce4c7a). A consistency calculation with matter and light coupled universally to g_m = −c²dt² + n(t)²dx².
  - **The six declared checks pass.** The coupling determines the clock factor. The measured redshift is n_o/n_e in both frames to 7×10⁻⁹. Atomic clocks and their ratios stay constant to 2×10⁻⁸. Cavities keep 2ℓ/c. Light and gravitational waves arrive together. D_L = (1+z)²D_A holds.
  - **The energy ledger closes to 9×10⁻¹⁵.** Radiation and free motion pay the field; bound systems and rest mass do not.
  - **History.** Ruler counts grow. With V = 0 the history is coasting with a turnaround. If the field's energy gravitates, a recombination-hot turnaround needs field energy of at least 5.5% of 3H²c²/8πG.

  See the [report](research_work/results/clock-completion/report.md).
- **Clock-gradient check CG-0** (protocol bc4a64b). One scalar's gradient cannot be both the rolling frame and the galaxy force:
  - carrying the total acceleration, X turns negative at every SLACS Einstein radius, inside 2.0 kpc in the Milky Way, in 23 SPARC galaxies and in the Solar System;
  - a limiting gradient leaves a constant 7×10⁻¹⁰ m/s² force;
  - only an excess-only, two-field form keeps the frame defined.

  So clock-field gravity needs a vector or a second field. See the [report](research_work/results/clock-gradient/report.md).
- **CR-1 ground-reference check (post hoc, requested in review).** The excess energies are recomputed against ground states relaxed afresh on each run's grid and on refined grids. See the addendum to the [CR-1 report](research_work/results/collective-reservoir/report.md).
- **Driven collective reservoir CR-1.** This is the third of the owner's branches: a coherent Gross–Pitaevskii–Poisson condensate, seeded with 10⁹ Msun in the Milky Way potential and fed by the conversion law's photon source.
  - **Supply.** The derived source falls short of the reference inventory in 10 Gyr by 1.3×10⁸ for the light constituent and 4.5×10¹¹ for the heavy one.
  - **Dynamics at the required rate.** Growth into the occupied mode is adiabatic, but the condensate contracts: it reaches 59% of the inventory at 7.4 kpc for the light constituent, and drops below the grid within 1 Gyr for the heavy one. The heavy constituent's analytic endpoint is a stable 1.6 pc soliton. Photon-shaped growth overheats or runs away.
  - **Numerics.** The first execution took 9.8 hours and hid its NaNs. The driver now fails fast, logs its progress and calibrates its time step on energy conservation. That calibration removed a 1.4% energy artifact.

  No channel supports a reservoir in this implementation. Zero self-interaction, the collection region and the 10 Gyr span are labeled assumptions. See the [report](research_work/results/collective-reservoir/report.md).
- **RPG-1 lens-geometry diagnostic (post hoc).** With flat FLRW angular distances and the published population masses, RPG-1's predicted SLACS Einstein radii rise to 0.48–0.58 of the observed values with Chabrier masses and 0.75–0.89 with Salpeter. The declared rule still fails. The remaining need, 2.1–2.7 times the Chabrier mass, is comparable to the archive's masses fitted to the same lenses' inner stellar motions. Those masses were fitted in a different geometry, so the similarity calls for a joint fit with identical conventions. (Corrected after review; this entry first said lensing and inner dynamics share one inner-mass shortfall.) See the addendum to the [RPG-1 report](research_work/results/radiation-polarized-gravity/report.md).
- `15bd1f6` **Model contract.** A proposal that one clock field does all three jobs:
  - PF-1's propagation field, completed by matter scales that co-evolve with it, so that the measured redshift survives with every clock ratio fixed;
  - its static galaxy response, which gives RPG-1's equation with a* tied to the redshift rate;
  - the companion condensate, reassigned to inner-core and cluster mass.

  The contract separates each branch's general hypothesis from its first implementation, and asks the owner to decide what "nonexpanding" means (decision D1). See [research_plan/model-contract.md](research_plan/model-contract.md). (Superseded the same day: the owner's review led to the revision entered above.)

## 2026-09-13, after the Codex handover

New entries are added here, newest first; [CURRENT-STATUS.md](CURRENT-STATUS.md) remains the authoritative summary.

- **Radiation-polarized gravity RPG-1.** This is the second of the owner's three branches. It solves the AQUAL equation with μ = x/(1+x) and the archived a* on the frozen SPARC and Milky Way baryons, using an axisymmetric solver validated against exact solutions.
  - Rotation matches algebraic simple MOND, with slightly worse RMSE.
  - The Milky Way vertical force gets worse, with RMS rising from 15.6 to 42.0 Msun/pc².
  - The declared lensing response predicts only 0.35–0.74 of the SLACS Einstein radii.
  - The field energy diverges logarithmically, and ξ = 0.118.
  - The declared convergence gate failed in one of its three galaxies.

  Not promoted. See the [report](research_work/results/radiation-polarized-gravity/report.md).
- **Evolving propagation field PF-1.** This is the first of three branches opened by the project owner.
  - One driven simulation reproduces the redshift, envelope stretch and pulse spacing, all equal to 1+z within 1.1×10⁻⁴, with photon number conserved.
  - The exchange between radiation and the index closes to 10⁻¹².
  - The light sector is exactly a coasting expanding universe with a = n. The claim of nonexpansion therefore rests on clock and ruler couplings that no tested completion supplies.
  - On Pantheon+, Δχ² = +35 against flat FLRW (Ω_m = 0.3) at equal freedom.

  Not promoted. See the [report](research_work/results/propagation-field/report.md).
- **Self-illumination pilot.** Companion fields derived from declared rotating emitters, with no co-rotation factor. Receiver drag falls to κ=0.19–0.51 in a flat-rotation disk, but the emitters pay instead. With every companion absorbed internally, the combined baryonic loss at ideal efficiency is a median 19× the baryons' angular momentum (above 1 in 146 of 149). This does not relieve the RB-1 debt. See the [report](research_work/results/self-illumination/report.md).
- `70cebef` **RB-1 consistency revision.** Review found three problems:
  - mixed incident spectra
  - an isotropic drag assumption
  - an energy ledger that omitted the receivers' energy change

  One incident field per labeled spectrum control, per-site pressure-tensor drag and a closed ledger now replace the first evaluation. The threshold-cut spectrum's supply multiplier doubles to 8.2×10^9. The ideal-case angular-momentum debt is a median 33.7. Population shapes are unchanged. The conclusion is narrowed to the tested prescription. See the [report](research_work/results/capture-to-orbit/report.md).
- `21147a1` **Repository consolidation.** Added CURRENT-STATUS.md and moved the dated notices here. The latest diagnostics gained entry points and archive regression checks.
- `e04d013` **Capture-to-orbit RB-1.** Receiver-assisted threshold production tested with orbit-averaged populations. Not promoted.

## 2026-09-13

- 17:23 · `3884b4f4` · from `research_plan/START-HERE.md`

Coma inverse diagnostic: positive finite Plummer profiles fit six saved shear points (single profile chi2 3.727 versus earlier NFW 3.855; flexible mixtures about 3.56). Outer mass is highly nonunique. These are total-gravity shapes, not baryon-subtracted companion predictions. Known n=5 polytropic support supplies a stationary candidate, not a derived capture law, stable state or shared K. See coma-inverse-report.md; PDF v1.5 predates this result.

- 16:51 · `1967ff8a` · from `research_plan/START-HERE.md`

Cross-scale transfer: frozen MOND-guided redistribution worsens both fiducial Milky Way baselines (6.761 to 9.259; 10.410 to 11.870 km/s), despite improving 13/18 sensitivity cases. Keep candidate and original side by side. Cluster transfer requires matched extended baryons, fixed inventory and an exterior lensing profile; no new cluster score. Conversion-only redshift unchanged at fixed alpha/path; no new timing mechanism. See mond-cross-scale-report.md. PDF v1.5 predates this update.

- 16:44 · `16d6100f` · from `research_plan/START-HERE.md`

MOND-guided fixed-inventory redistribution improves SPARC rotation predictions: trained shared mixture f=0.923557 gives train/validation/test RMSE 20.468/27.136/16.241 km/s versus original 29.025/32.495/23.591. No negative sampled target shells; 21 galaxies require inventory caps. This borrows MOND radial guidance, is not derived capture physics, and does not establish superiority to MOND or lensing success. See mond-inventory-report.md. PDF v1.5 predates this result.

- 16:10 · `47b4d4f8` · from `research_plan/START-HERE.md`

Shared deposition scaling tested: common NFW-shaped amplitude/size laws from population mass, half-light size and exact-third eta fail frozen six-system transfer. Omitted-system motion chi-squared is 12667/8821/25236 for 2/3/4 shared parameters, with lens RMS 41/41/66 percent. Keep original reference; these are inverse-profile failures, not gravity-solver failures. See shared-deposition-report.md. PDF v1.5 predates this test.

- 15:58 · `342af88c` · from `research_plan/START-HERE.md`

Halo equivalence passes: frozen NFW density interpreted as companions reproduces 40 stellar-motion bins across six systems within 0.000201 km/s and nonzero halo lensing within 5.71e-14 relative. Production axisymmetric solver spherical-NFW control also passes. No refit or capture-origin claim; see halo-equivalence-report.md. PDF v1.5 predates this diagnostic.

- 14:19 · `ca3cd25c` · from `research_plan/START-HERE.md`

Well-depth capture tested: three bounded laws across 149 galaxies all select zero modification under the declared training proportional-error criterion. Weak shallow suppression improves refined test RMSE from 23.587 to 23.307 km/s but worsens training/validation log error and middle-radius underprediction. Boosting deep-well capture worsens the speed scores. Keep exact-third reference. See research_work/results/companion-extensions/depth-capture-report.md. Sphericalized depth and stationary feedback only; PDF v1.5 predates this result.

- 13:43 · `5c75ebc5` · from `research_plan/START-HERE.md`

Directional capture tested in four oblate geometries, each at raw or fixed-mass normalization. Refined thin/raw mixture gives test RMSE 23.334 vs reference 23.591 km/s, but same-inventory spherical control gives 23.243. Fixed-mass flattening raises inner and outer speeds; no variant fixes both radial signs. Keep reference. User well-depth feedback suggestion is distinct from this prescribed-opacity test; review earlier feedback branches before extending it directionally. See directional-capture-report.md. PDF v1.5 predates this result.

- 13:27 · `daec3d88` · from `research_plan/START-HERE.md`

Homology change tested: fitting one shared size/retention parameter gives s=0.999614, effectively the unchanged reference. Larger reductions worsen rotation errors; at s=0.90 outer mean error grows from +10.41 to +12.75 km/s. Keep exact-third reference; no lensing or formation success claimed. See homology-retention-report.md. PDF v1.5 predates this result.

- 13:02 · `b1595ef5` · from `research_plan/START-HERE.md`

Quick residual diagnostic completed: unchanged companion predictions underpredict inner/middle speeds by about 7 km/s and overpredict outer speeds by about 10 km/s, with the same mean signs in each old split. Only two of six transferred lens fits improve the matched stellar-only score. Radial shape deserves a bounded follow-up; a uniform amplitude increase does not fix both signs. See shared-residual-map-report.md. No new fit or blind validation; PDF v1.5 predates this diagnostic.

- 12:36 · `5a4109d6` · from `research_plan/START-HERE.md`

Current paper v1.5 incorporates the capture-recoil result through 2929cc8. Section 8 now prioritizes three bounded tasks: consistent existing-data propagation comparison; capture-to-orbit viability; and a shared-law residual map before adding parameters. These tasks are not yet executed. Earlier PDF notices below are historical.

- 12:13 · `2929cc89` · from `research_plan/START-HERE.md`

New capture diagnostic: receiver-assisted absorption or separate-particle production can conserve momentum with a small recoil budget; two-packet production needs favorable energies and directions. These are conditional kinematic branches, not capture rates or orbit formation. See research_work/results/companion-extensions/capture-recoil-threshold-report.md. The v1.4 PDF predates this calculation.

- 12:04 · `4f2a9c4b` · from `research_plan/START-HERE.md`

Current paper: v1.4 (13 September 2026), evidence through 0756dec. The review now includes matched rotation comparisons with ordinary matter, MOND and specified NFW models; qualified Coma and lens comparisons; conversion-rate versus FLRW and published H0 values; timing/brightness limitations; and the latest reservoir scale and orbital-support results. Older notices below describe historical paper versions. No new fitting was performed for this update.

- 11:53 · `0756decf` · from `research_plan/START-HERE.md`

Circular support result: randomly oriented circular particle orbits can maintain the tested profiles in the frozen Newtonian potential, with individual radial stability. This is not collective stability or a capture mechanism. The extreme reservoir needs about 1.05e60 J in orbital kinetic energy, small relative to its mass-energy but explicitly accounted. See circular-reservoir-support-report.md; PDF v1.3 predates this result.

- 11:48 · `d19c4cd9` · from `research_plan/START-HERE.md`

New support constraint: outward-rising deposited density excludes a spherical bound collisionless isotropic particle completion for 4/6 free and 5/6 transferred Chabrier profiles. It does not exclude all companion interpretations. Tangential support, wave stresses or evolution need explicit calculations. See reservoir-isotropic-support-report.md; PDF v1.3 predates this result.

- 11:45 · `45ab1dea` · from `research_plan/START-HERE.md`

Frozen branch comparison: the lower-mass J1621 ac/Re=10 branch predicts stronger gravity nearer the galaxy (843 km/s circular diagnostic and 0.03823 reduced shear at 300 kpc), unlike the almost-equivalent inner fit at ac/Re=100 (116 km/s and 0.000984). No outer observations were fitted. See reservoir-branch-outer-report.md; PDF v1.3 predates this calculation.

- 11:42 · `2b8d8bb3` · from `research_plan/START-HERE.md`

New fixed-scale result: J1621+3931 can store about 486 times less mass than the extreme fit with only 0.000850 higher chi-squared. Stellar-only differs by just 0.332. The huge reservoir is not uniquely required by these inner data. See reservoir-scale-scan-report.md; v1.3 PDF predates this scan.

- 11:38 · `e42733fa` · from `research_plan/START-HERE.md`

Current paper: v1.3 (13 September 2026) incorporates evidence through 2054afb, including matched free-companion/NFW fits, the extreme-reservoir inventory and frozen outer predictions. No new fitting was done for this update. Notices below referring to v1.2 are historical; the current manuscript and PDF include those results.

- 11:30 · `2054afb7` · from `research_plan/START-HERE.md`

Frozen outer predictions: the extended J1621+3931 free-fit reservoir implies a circular-speed peak about 4246 km/s near 9.17 Mpc and radial shear at sampled 1 and 3 Mpc, returning tangential by 10 Mpc. This is an extrapolated test target, not new observed data. No additional critical ring is found in the finite scan. Numerical projection checks pass after resolving the extended envelope. The v1.2 PDF predates this experiment.

- 11:19 · `eb8a6d5a` · from `research_plan/START-HERE.md`

Matched six-parameter companion fit: total chi-squared 26.49 versus free NFW 49.71. This is an inverse fit, not validation of one-third retention: free density absorbs that factor. Some density normalizations grow by hundreds/thousands; J1621+3931 reaches the largest capture scale and implies about 8.42e16 solar masses of stored effective mass. Supply and large-radius gravity remain unproved. The v1.2 PDF predates this experiment.

- 11:11 · `fe68306e` · from `research_plan/START-HERE.md`

Broader NFW comparison changes the assessment: fitting halo size and strength gives total chi-squared 49.71 versus companion 111.27-113.56. It adds two fitted parameters per galaxy and retains many boundary solutions. Earlier lens-matched NFW scores were restricted shape tests, not a general advantage over dark matter. See free-nfw-report.md. The v1.2 PDF predates this experiment.

- 11:04 · `a8649dd6` · from `research_plan/START-HERE.md`

NFW shape control: at the same added bending at the lens radius, three prescribed NFW scales give total chi-squared 162.62-194.57 versus companion 111.27-113.56. This restricted inherited-normalization comparison favors the current profile, not a general superiority claim over dark matter. Best NFW score lies at the largest tested scale. The v1.2 PDF predates this experiment.

- 10:59 · `025b3453` · from `research_plan/START-HERE.md`

Matched stellar-only control: total motion chi-squared is 175.09 versus 113.56/111.27 with companions, a 35%-36% aggregate improvement. Only J0037-0942 and J1402+6321 improve with the added profile; the other four slightly prefer stellar-only. This supports conditional profile usefulness, not photon origin or a complete solution. The v1.2 PDF predates this experiment.

- 10:54 · `7622bdf3` · from `research_plan/START-HERE.md`

Free orbital transition radius: total all-motion chi-squared improves about 30%, to 113.56/111.27, but J1402+6321 still contributes roughly 57% of the residual and multiple stellar parameter boundaries remain. All bins are fitted; no independent prediction or physical solution is established. See orbit-transition-report.md. The v1.2 PDF predates this experiment.

- 10:50 · `f50d3d3a` · from `research_plan/START-HERE.md`

All-motion diagnostic: fitting all 40 bins reduces total chi-squared from 193.62/190.91 to 161.81/159.83, but J1402+6321 remains the dominant discrepancy. The outer bins are now fitted, not predicted; orbital boundary dependence remains. See all-motion-orbits-report.md. The v1.2 PDF predates this experiment.

- 10:47 · `ecdc078f` · from `research_plan/START-HERE.md`

New constrained orbital fit: inner chi-squared 87.58/86.79 and outer residual-square sums 106.03/104.12. All twelve fits satisfy the imposed central and sampled slope/anisotropy conditions, with six cases on the central bound. Positivity, stability and joint agreement remain unresolved. See constrained-gradient-orbits-report.md. The v1.2 PDF predates this experiment.

- 10:42 · `2a0a8b46` · from `research_plan/START-HERE.md`

New result after paper v1.2: combining stellar gradients and radial orbits lowers inner chi-squared to 82.40/81.67, but outer residual-square sums remain 103.20/101.32 and six of twelve fits fail the necessary orbital condition near extrapolated centers. See the combined-fit report; no physical solution is adopted. The v1.2 PDF predates this experiment.

- 10:38 · `20fcf40d` · from `research_plan/START-HERE.md`

Paper updated to v1.2 (13 September 2026): theory-basis.md and output/pdf/theory-basis.pdf now incorporate capacity/recycling, joint lensing, stellar-gradient and radial-orbit results, and the necessary orbit-slope check. No joint physical solution is claimed. Earlier snapshot notices below are historical.

- 10:38 · `20fcf40d` · from `research_plan/START-HERE.md`

**Necessary orbital-consistency check.** All 18 tested population/system/orbit cases satisfy sampled gamma>=2 beta over 1e-6..100 Re, with minimum margins 0.324-0.325. This necessary condition applies to a separable augmented-density completion with beta0<=1/2; it neither constructs a positive distribution function nor establishes stability. Outer-motion discrepancies remain. [Report](research_work/results/companion-extensions/orbit-slope-check-report.md).

- 10:28 · `3f57bcc3` · from `research_plan/START-HERE.md`

**Radial-orbit diagnostic.** A known varying-beta profile reduces exact-lens inner chi-squared to 448-454 under previous bounds. Expanding the outer bound for the three boundary cases lowers combined totals to 106-107, but all hit beta_infinity=0.95 and outer residual sums remain 114-116 versus free fits at 53-55. Orbital restrictions materially affect the tension; physical distribution functions and a joint fit remain unestablished. [Report](research_work/results/companion-extensions/radial-orbits-report.md).

- 10:22 · `c0a5dff4` · from `research_plan/START-HERE.md`

**Stellar M/L gradient diagnostic.** A bounded fixed-Re gradient reduces exact-lens inner chi-squared from 611-614 to 472-477, still far above the free inner-motion fits near 45. Two gradients and two anisotropies hit bounds per population. Large-gradient endpoint checks give modest further gains, not reconciliation. These are fitted nuisance profiles, not observed population gradients or an adopted gravity change. [Report](research_work/results/companion-extensions/stellar-gradient-report.md).

- 10:16 · `ae6aa6af` · from `research_plan/START-HERE.md`

**Exact-angle orbital compatibility.** Fixing stellar mass to the catalogue lens angle and refitting constant anisotropy raises inner-motion chi-squared totals from 41-46 to 585-625; outer residuals worsen too. All fitted anisotropies remain interior. This conditional incompatibility persists across the four profiles and both population proxies; the lens angle is consumed as a constraint, not predicted. [Report](research_work/results/companion-extensions/exact-lens-orbits-report.md).

- 10:12 · `0fea706f` · from `research_plan/START-HERE.md`

**Lens bending budget.** At the current inner-star fits, stars alone exceed the catalogue-angle bending in J0037-0942, J1204+0358 and J1402+6321, under all four profiles and both population assumptions. Positive companion changes alone cannot resolve these three while stars and geometry stay fixed. This is a conditional stellar-model/geometry incompatibility, not a measured excess stellar mass. [Report](research_work/results/companion-extensions/lens-bending-budget-report.md).

- 10:07 · `976e5b1f` · from `research_plan/START-HERE.md`

**Capacity profiles transferred to six lenses.** Galaxy-trained amplitudes give lens-angle RMS 13.07%-13.19% for the matched reference, 13.29%-13.46% for local filling and 13.05%-13.16% for refined recycling. Recycling helps only 3/6 systems per population and performs worse on galaxy test motions; outer stellar discrepancies remain. No joint solution is adopted. Conditional geometry and luminosity proxies are retained. [Report](research_work/results/companion-extensions/capacity-lensing-report.md).

- 10:01 · `24534de4` · from `research_plan/START-HERE.md`

**Normalization correction.** Recent absolute capacity and conservative-recycling mass outputs used C0 instead of the source amplitude A=2 C0. Those masses are corrected upward by two; area/capacity halves. Rotation predictions and scores are verified unchanged for all 149 galaxies because mass ratios cancel the factor. Main-paper equation (7) and the original lens code were already correct. New-profile lens transfer remains pending. [Audit](research_work/results/companion-extensions/capacity-normalization-report.md).

- 09:57 · `9abdd909` · from `research_plan/START-HERE.md`

**Finite local formation.** Empty local stores under the full-opacity field were tested at six common durations plus equilibrium. Training selection chooses u=10, giving validation/test RMSE 29.41/25.93 km/s versus matched reference 30.82/21.41. The common-history grid is not an improvement across partitions; this is a fixed-field diagnostic, not coupled time-dependent transport or a physical age estimate. [Report](research_work/results/companion-extensions/local-formation-report.md).

- 09:53 · `eabf07e9` · from `research_plan/START-HERE.md`

**Conservative-recycling endpoint.** Isotropic same-channel re-emission admits stationary J=1 under an isotropic boundary bath, giving density Cg eta(X). Training-adjusted validation/test RMSE is 28.82/30.88 km/s versus matched reference 30.82/21.41. Stationary recycling avoids a separate escape-channel drain but does not establish formation energy, post-source retention or a physical return mechanism. [Report](research_work/results/companion-extensions/conservative-recycling-report.md).

- 09:50 · `5e768634` · from `research_plan/START-HERE.md`

**Occupancy-dependent transport.** Solving empty-site capture and local filling together across 149 inputs increases deposits further. Training-adjusted validation/test RMSE is 28.64/26.70 km/s versus matched reference 30.82/21.41. Stationary absorbed power is balanced by an explicitly assumed non-recaptured release channel; its microscopic cause and local support remain open. [Report](research_work/results/companion-extensions/occupancy-transport-report.md).

- 09:44 · `6d4744e0` · from `research_plan/START-HERE.md`

**Local-capacity comparison.** Replacing CgJ eta(X) by Cg eta(XJ) with fixed transport greatly increases deposits in weakly illuminated regions. Training-only amplitude adjustment yields mixed validation/test RMSE 28.63/25.19 km/s versus matched control 30.82/21.41. No replacement is adopted; occupancy-dependent transport and the physical capacity remain unresolved. [Report](research_work/results/companion-extensions/local-capacity-report.md).

- 09:40 · `2ea6d50d` · from `research_plan/START-HERE.md`

**Capture/capacity audit.** Interpreting the reference eta=1 profile as capacity gives M_cap=(C/k0) sigma. Incoming power and capacity therefore scale with the same absorption area: size alone cannot change the normalized background rate. An alternative unattenuated volume capacity yields a factor 0.495-0.981 across existing scales, but changes the storage/profile mapping and needs its own local solution. [Report](research_work/results/companion-extensions/capture-capacity-report.md).

- 09:38 · `49c90996` · from `research_plan/START-HERE.md`

**Uniform-background comparison.** A common additive input strong enough to preserve 90% of every equilibrium deposit worsens frozen galaxy rotation scores. A training-only amplitude adjustment partly compensates, but validation/test RMSE remains 31.54/22.52 km/s versus 30.82/21.41 for the identically adjusted no-background control. This tests a uniform rate proxy with fixed geometry, not actual distant-source transport. [Report](research_work/results/companion-extensions/background-retention-report.md).

- 09:34 · `214d6af7` · from `research_plan/START-HERE.md`

**Continuing-source checkpoint.** With initially equilibrated stores, maintaining 90% of the deposit indefinitely after fading requires a constant input floor of 34.09%-68.62% across the existing galaxy proxies. This is a requirement on external supply, not a measurement of it. Finite-threshold kinetics requires continuing energy turnover; actual sources and rates remain unspecified. [Report](research_work/results/companion-extensions/threshold-floor-report.md).

- 09:31 · `65c8d1c0` · from `research_plan/START-HERE.md`

**Main paper updated to v1.1.** The manuscript and verified PDF now incorporate findings through `6e6b022`: coupled Milky Way endpoints, transfer and scattering constraints, reciprocal storage, and exact-third threshold formation and retention. These remain partial, conditional results. [Manuscript](papers/cumulative-time-companions/theory-basis.md) | [PDF](output/pdf/theory-basis.pdf).

- 09:26 · `6e6b0224` · from `research_plan/START-HERE.md`

**Latest checkpoint: formation versus retention.** In the stipulated threshold kinetics, individually filling to 90% equilibrium and then removing the source loses the first 10% in 4.92%-41.79% of the formation time. The ratio is independent of the unknown rate. A long release tail is not permanent retention; actual histories and other protection laws remain open. See [report](research_work/results/companion-extensions/threshold-memory-report.md).

- 09:25 · `b76b60fe` · from `research_plan/START-HERE.md`

**Latest checkpoint: finite threshold formation histories.** Starting empty, 90%-equilibrium times vary by about 1209 across the existing galaxy sample under one rate scale. Common-duration scans give mixed score changes; none is selected. Equilibrium requires continuing energy turnover and is not automatically permanent storage. See [report](research_work/results/companion-extensions/threshold-history-report.md).

- 09:23 · `9d94da33` · from `research_plan/START-HERE.md`

**Latest checkpoint: finite thresholds on 149 galaxies.** With all original parameters fixed, a 12-decade threshold range changes speeds by at most 0.640 km/s; a six-decade range reaches 7.03 km/s and worsens aggregate validation/test velocity RMSE. This is an exposed-sample sensitivity comparison, not a new blind test or a physical derivation. See [report](research_work/results/companion-extensions/threshold-galaxies-report.md).

- 09:21 · `92c7f35a` · from `research_plan/START-HERE.md`

**Latest checkpoint: exact-third threshold mixture.** A known positive mixture of ordinary saturation responses represents the reference exactly. The distribution is inverse-designed, not derived. Finite cutoffs alter the curve; stipulated release kinetics produce slow depletion rather than permanent storage. See [report](research_work/results/companion-extensions/threshold-mixture-report.md).

- 09:18 · `d2fcbc94` · from `research_plan/START-HERE.md`

**Latest checkpoint: finite reciprocal protected reservoir.** Collective transfer can replenish leakage faster, but no-decay equilibrium protected occupancy is (n+1)/(n+1+r n), independent of domain size. Returning waves and finite capacity remain essential. Fifty-four stationary cases include the reverse export and close the energy ledger. See [report](research_work/results/companion-extensions/collective-reservoir-report.md).

- 09:17 · `9c31e959` · from `research_plan/START-HERE.md`

**Latest checkpoint: collective-state normalization.** An empty domain has first-excitation rate factor N, not N^2. Partly excited symmetric states can enhance rates more strongly, but require preparation energy and enhance reverse transitions too. The tested closed ladders reach zero net current. A physically specified protected export remains necessary. See [report](research_work/results/companion-extensions/collective-normalization-report.md).

- 09:15 · `3088e3a8` · from `research_plan/START-HERE.md`

**Latest checkpoint: angular-rate consistency.** Matching the source-size-tracking angle with the existing Gaussian overlap suppresses the transfer rate as s^-2 exp[-(s/s_c)^2]. At fixed coupling and target density, redshift accumulation saturates near the source. The favorable observer moments cannot be combined with an unchanged long-path rate in this branch. See [report](research_work/results/companion-extensions/kernel-rate-closure-report.md).

- 09:13 · `3ac30254` · from `research_plan/START-HERE.md`

**Latest checkpoint: observer geometry and delay.** Accounting for kick location gives about 0.520 arcseconds and 5025 seconds mean delay for the fixed-angle optical example over 30.66 Mpc. A source-size-tracking kernel yields much smaller moments but lacks a derived local cause. Stationary scattering still does not explain event time dilation. See [report](research_work/results/companion-extensions/observer-transport-report.md).

- 09:11 · `fb101ddf` · from `research_plan/START-HERE.md`

**Latest checkpoint: angular escape geometry.** In an isolated-source toy geometry, small direction changes can reach less populated modes and favor loading. The illustrative cumulative direction spread can be small for distant sources, so this is not a universal image-blurring exclusion. Real source coverage, outgoing-mode filling, kernel strength and observer geometry remain unmodeled. See [report](research_work/results/companion-extensions/angular-outlet-report.md).

- 09:09 · `604b7016` · from `research_plan/START-HERE.md`

**Latest checkpoint: self-consistent source and outlet.** A replenished high-energy field and escaping low-energy field can sustain loading with occupations solved together. Slow escape throttles the rate. The terminal protected sink remains imposed; a one-step 2 eV to 1e-8 eV loading route requires processed light power 2e8 times net assembly power. See [report](research_work/results/companion-extensions/source-outlet-report.md).

- 09:07 · `df8adc0c` · from `research_plan/START-HERE.md`

**Latest checkpoint: occupation-driven assembly bias.** A reciprocal two-mode interaction favors loading only when the input mode is more occupied than its paired output mode. Smooth thermal examples do not provide that preference; a closed reservoir consumes an imposed contrast. Sustained loading needs specified source and escape channels, with momentum and field depletion included. See [report](research_work/results/companion-extensions/mode-driven-bias-report.md).

- 09:06 · `c6aca10c` · from `research_plan/START-HERE.md`

**Latest checkpoint: restarted assembly and recycling.** Repeated attempts give finite mean completion times for the constant-rate ladder, but unbiased assembly can circulate about 200,000 times its net energy in the representative packet. Keeping replacement input below net assembly then requires returned-energy losses below about five parts per million. A forward bias helps, but its physical source remains unspecified. See [report](research_work/results/companion-extensions/assembly-recycling-report.md).

- 09:04 · `2c246711` · from `research_plan/START-HERE.md`

**Latest checkpoint: collective release and assembly.** Larger release quanta conditionally reduce the finite-site lifetime bound as K^-4; mere grouping with unchanged radiation does not. Twenty-four designs quantify packet sizes, while a reversible assembly ladder shows strong sensitivity to forward versus reverse rates. A microscopic accumulation mechanism and release spectrum remain required. See [report](research_work/results/companion-extensions/collective-assembly-report.md).

- 09:02 · `dd7bae74` · from `research_plan/START-HERE.md`

**Latest checkpoint: physical site/mode budget.** Connecting finite protected excitations to the fitted Milky Way energy inventory yields a necessary stationary lifetime exceeding 10^32 to 10^41 years in the 120 kpc, 90%-protected examples. This assumes all extra mass-energy is in tiny excitations and a homogeneous coupled escape channel; it is not a cosmic-age exclusion. Collective storage, stable states or a separately justified release channel remain alternatives. See [report](research_work/results/companion-extensions/site-mode-budget-report.md).

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Latest checkpoint — escaping-wave feedback:** An energy-conserving finite-site model now includes stimulated return while radiation escapes. Eighteen illustrative histories show delayed protection and extra loss through the bright state. Full stationary protection requires a stipulated non-decaying state and continued illumination; storage saturates at one excitation per site. Site/mode counts and physical lifetimes remain un-derived. See [report](research_work/results/companion-extensions/escape-feedback-report.md).

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Resonance ensemble candidate:** An inverse-designed ensemble of independent damped resonances with weight proportional to resonance energy squared gives response approximately proportional to 1/E, canceling the Gaussian spatial model's residual E scaling. Across 0.01–100 eV, broad-cutoff p=2 cases vary by 0.0063–0.1891% after normalization. This is a conditional rate-shape success, not an absolute redshift fit: mode origin, damping-energy channels, extinction, image preservation and event timing remain unresolved. [Calculation](research_work/results/companion-extensions/resonant-spectrum-report.md).

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Saturation profile tested:** A constant donor-density saturation closure was calibrated to the old mean contraction in baseline I and frozen for II. All six coupled endpoints converge with positive release, but all worsen their uniform-contraction counterparts: all-bin RMS 7.21–7.23/8.19–9.47 km/s versus 6.27–6.28/8.10–9.26. Mean settling alone does not preserve the radial gravity profile. This candidate is recorded, not adopted as the new reference. [Calculation](research_work/results/companion-extensions/saturation-report.md).

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Coupled redistribution:** Coupled circular endpoints now converge for both matter baselines and all three receiver bands. Receivers expand 5.4–8.2%, leaving 3.65–4.80e50 J to release. All-bin speed RMS is 6.27–6.28 km/s for I and 8.10–9.26 for II, versus no-settling 6.34/10.91. This preserves a conditional internal-exchange route; it does not derive the torque, contraction, formation rate or collective stability. [Calculation and limitations](research_work/results/companion-extensions/coupled-torque-report.md).

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Shared conservative settling tested:** [Capture and cooling ledger](research_work/results/companion-extensions/settling-report.md) uses one phase mass/contraction pair for both Milky Way matter baselines. Inner RMS improves 3.68/11.83 → 2.74/6.70 km/s; outer RMS remains 8.36/9.78. Conditional virial accounting requires about 4.5e50 J of released energy to leave the bound subsystem. All grid cases permit nonnegative release, so conservation does not derive the fitted contraction. Capture recoil can be balanced using known relativistic kinematics; rates, local support, feedback and phase formation remain unresolved. All goals remain open.

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Retention lens transfer and one-third test:** [Executed results](research_work/results/isotropic-galaxy-transfer/retention-lensing-report.md) improve outer stellar motion but worsen lens RMS to 13.94–14.11% under two declared photometric-to-3.6-micron proxy mappings. Exact q=1/3 produces nearly unchanged galaxy performance after training-only refits; simplicity is supported, a fundamental exponent is not. The lens run retains fitted q and is not a separate one-third test.

- 09:00 · `f1e4b8f3` · from `research_plan/START-HERE.md`

**Three reservoir options now formulated and compared:** [Sun–Earth–Moon specification](research_plan/sun-earth-moon-three-branches.md) and [executed transfer results](research_work/results/isotropic-galaxy-transfer/three-reservoirs-report.md). A common normalization calibrated on the Sun fails to predict terrestrial/lunar mass under the declared source/path proxies. Static stored and steady-renewal versions are degenerate; a freely escaping cloud lacks inventory, while independently self-bound states still require a support law.

- 08:56 · `eb04a066` · from `research_plan/START-HERE.md`

**Nonthermal escape occupancy:** Free escape of nonthermal low-energy relaxation waves does not justify empty coupled modes in the tested uniform galaxy-scale source. The 120-kpc orbital-heat example over 1 Gyr gives mean occupation 4.67e16 for a narrow band, so negligible thermal return cannot be assumed solely from escape. This is an optically thin occupation diagnostic, not a computed reabsorption rate; actual coupling, geometry and population feedback remain open. [Analysis](research_work/results/companion-extensions/escape-occupation-report.md).

- 08:54 · `71fee377` · from `research_plan/START-HERE.md`

**Physical bath budget:** A physical two-polarization, c-speed thermal bath has far too little local cold capacity for the tested heat budgets. In the 120-kpc larger-gap benchmark, prior orbital cooling exceeds available capacity by 5.58e20; a cold thermal outlet would need 2.91e26 years. These are conditional local-volume and throughput requirements, with no imposed universe age. Nonthermal escape and distinct slow/massive baths remain alternatives requiring derivation. [Calculation](research_work/results/companion-extensions/bath-budget-report.md).

- 08:52 · `3fdf06d4` · from `research_plan/START-HERE.md`

**Finite-bath feedback:** A finite-temperature bath now warms under protected-state relaxation, with heat capacity and cooling-outlet energy explicitly included. Small capacity substantially weakens protection; larger capacity or exported heat can preserve it in the conditional model. Thirty-six cases and six refinements close the energy ledger. Physical bath degrees of freedom, volume, heat capacity and cooling channels remain to be supplied. [Analysis](research_work/results/companion-extensions/finite-bath-report.md).

- 08:50 · `4b5b9660` · from `research_plan/START-HERE.md`

**Protected-state candidate:** A finite-site illuminated-to-protected-state model now includes thermal backflow, bright-state radiation and a signed bath-energy ledger. Twenty-seven cases show stronger thermal isolation for deeper relaxation, at the cost of less retained energy per site. Direct protected-state decay is omitted only as an optimistic assumption. No absolute lifetime, site capacity or microscopic selection rule is derived. [Calculation and energy costs](research_work/results/companion-extensions/protected-state-report.md).

- 08:47 · `43877762` · from `research_plan/START-HERE.md`

**Reciprocal storage check:** An ideal flat forward fractional-loss law, combined with reciprocal heavy-store scattering, fixes inverse/forward ratio to (1+gap/E)^3. In the inherited bosonic pair-mode model this yields finite stationary storage and no net continued photon-energy deposition when spontaneous loss is absent. The small-gap E/3 per-mode limit is unrelated to the astrophysical one-third retention exponent. Protected-state transport, mode count and actual decay remain open. [Derivation](research_work/results/companion-extensions/reciprocal-storage-report.md).

- 08:45 · `d5844942` · from `research_plan/START-HERE.md`

**Joint width and survival:** An exact discrete energy-transfer/removal process now couples spectral width and photon survival. For illustrative width 1e-5 and 90% survival at z=1, 2 eV light requires gap <=2e-10 eV and removal/useful ratio <=2.11e-11. The previous small-gap drift survival estimates remain accurate, but their brightness allowance alone did not guarantee this linewidth. These are model design conditions, not measured exclusions or derived branching ratios. [Exact calculation](research_work/results/companion-extensions/discrete-transfer-report.md).

- 08:43 · `d3a46ea7` · from `research_plan/START-HERE.md`

**Resonant loss channels:** The resonance candidate now has an explicit photon-survival requirement and drift energy ledger. At 2 eV, gap 1e-8 eV and z=0.00766048, an illustrative 90% survival requires removal/useful rate ratio below 6.93e-8; at z=1 it is 1.05e-9. Removed radiation energy is separately counted. Damping is not automatically absorption, and the actual common-coupling branching ratios remain unknown. [Analysis](research_work/results/companion-extensions/resonance-loss-report.md).

- 08:39 · `93104861` · from `research_plan/START-HERE.md`

**Finite spatial response:** A stipulated Gaussian receiving overlap makes scattering forward-directed and changes the high-energy fractional-loss scaling from cubic to approximately linear. It still fails common fractional redshift. Increasing size also suppresses the integrated rate; illustrative small deflections require compensating coupling or mode abundance that has not been derived. The angular diagnostic is not an observed image-width prediction. [Calculation](research_work/results/companion-extensions/spatial-response-report.md).

- 08:37 · `9ec9d403` · from `research_plan/START-HERE.md`

**Continuum color bound:** For the previously derived forward point interaction, any positive incident-energy-independent gap spectrum obeys d ln(alpha)/d ln(E)>=3. A continuum of soft states alone therefore cannot produce common fractional redshift. Eighty-four numerical examples confirm the analytic bound; spatial response, inverse populations or a coherent evolving field must change the assumptions. This is a conditional interaction result, not an exclusion of the overall hypothesis. [Derivation and consequences](research_work/results/companion-extensions/continuum-color-report.md).

- 08:34 · `8a3bd6ac` · from `research_plan/START-HERE.md`

**Microscopic scale bridge:** The redshift-to-bound-state bridge now separates traveling quantum energy, creation gap and fitted phase mass. Under the illustrative Poisson linewidth case, an optical transfer supplies at most 2.62e-8 eV versus a hypothetical 17.8 eV creation gap. Accumulated storage or soft collective modes remain possible; one redshift step cannot be equated to creation of one fitted-mass particle. The width allowance is a design case, not a measured exclusion. [Calculation and alternatives](research_work/results/companion-extensions/quantum-bridge-report.md).

- 08:30 · `bdac724f` · from `research_plan/START-HERE.md`

**Stopping mechanism audit:** The coupled endpoint family has no sampled energy stop between contraction s=1 and 0.4; the prior fitted s=0.803 remains downhill toward further settling. Inverse restoring terms can create a local minimum but consume most of the release budget and need a new force/support solution. Density-dependent transport shutdown and finite receiver capacity remain separate candidates, with no independently derived threshold yet. [Results and alternatives](research_work/results/companion-extensions/stopping-report.md).

- 08:19 · `427effaa` · from `research_plan/START-HERE.md`

**Torque and cooling separated:** [Budget calculation](research_work/results/companion-extensions/torque-report.md) finds direct orbital ray-emission braking needs 1,394/1,504 times the settling energy. Internal receivers at 15-120 kpc need roughly 3-4% specific angular-momentum increases and 4-8% outward motion in frozen-potential scenarios. This motivates internal torque plus separate cooling; coupled self-gravity, orientation matching and rates remain unresolved. The ray bound does not exclude coherent wave torque or extra reservoir depletion.

- 08:17 · `54c95031` · from `research_plan/START-HERE.md`

**Orbital support checked:** [Existence and transport test](research_work/results/companion-extensions/orbital-support-report.md) rejects isotropic energy-only support for the settled profiles through a necessary positivity condition. Nonnegative circular-orbit shells preserve the profiles in the spherical approximation and pass individual radial stability checks. Moving the selected component between circular states requires about 17-18% specific angular-momentum reduction. This is not formation, collective stability or a new speed prediction.

- 08:14 · `c00ad96a` · from `research_plan/START-HERE.md`

**Shared pressure-law equilibria tested:** [Results](research_work/results/companion-extensions/interaction-support-report.md) preserve exact-third mass but do not preserve the speed-fit improvement. Gamma=5/3 gives all-bin RMS 26.93/19.47 km/s; gamma=2 gives 43.73/35.98 versus original 6.34/10.91. Gamma=4/3 boundary-limited candidates are rejected. Outward-rising density intervals also rule against a single positive-compressibility barotropic fluid reproducing the exact fitted shape. Next: orbital-distribution positivity and separately supported components.

- 08:11 · `15389646` · from `research_plan/START-HERE.md`

**Local support audit:** [Frozen phase-profile check](research_work/results/companion-extensions/local-support-report.md) finds the assigned compact fraction incompatible with ideal thermal condensate support. At 5 kpc, assigned 82-85% contrasts with pressure-implied 2.5-8.3%; retaining those fractions needs about 94-95% additional pressure. The empirical rotation fit remains, but interactions, orbital/wave support or revised phase physics are needed. No refit or new observational success.

- 08:02 · `f3a17892` · from `research_plan/START-HERE.md`

**Theory-basis paper v1.0, 13 September 2026:** [Editable paper](papers/cumulative-time-companions/theory-basis.md) and [PDF](output/pdf/theory-basis.pdf) consolidate the core postulates, twelve provenance-labeled equations, current successes/failures and eight requirements for a competitive theory. Nuclear emission is retained as an optional track alongside settling, conversion/time, lensing/field response and cross-system transfer. No superiority claim or newly completed physical goal is made.

- 07:56 · `35a7ea79` · from `research_plan/START-HERE.md`

**Fictional nuclear-release branch:** [Energy requirements](research_work/results/companion-extensions/nuclear-release.md) reinterpret Sgr A*/M87 core light as companion-powered. Luminosity constrains throughput, not a unique capture law. The Milky Way settling budget divided by the faint Sgr A* reference luminosity gives 1.43e14 years, exposing the need for channel/location/history assumptions. No Hawking detection or brightness prediction is claimed.

- 07:46 · `3846e2f5` · from `research_plan/START-HERE.md`

**Companion extensions analyzed:** [Refraction, phase selection and conversion](research_work/results/companion-extensions/report.md) retain exact-third capture. No selected refraction variant improves outer Milky Way predictions. An ideal-boson phase selector improves inner RMS from 3.68/11.83 to 1.81/2.85 km/s; outer RMS is 8.36/9.29 versus matched controls 8.36/9.78. These use a gas-monopole pilot, inner-bin fitting and previously seen data. Cross-baseline transfer exposes sensitivity; massive bound-state formation and self-consistent support remain missing. Whole-photon mixing is not survivor redshift; a Poisson inelastic-loss calculation quantifies illustrative spectral broadening. No lensing or energy-supply solution, new-galaxy validation, or closed goal is claimed.

- 07:32 · `a2b9519e` · from `research_plan/START-HERE.md`

**Milky Way circular-speed transfer completed:** [Current frozen predictions](research_work/results/isotropic-galaxy-transfer/milky-way-current-report.md) compare 38 Eilers bins with two archived stars+gas baselines. Exact-third RMS is 6.76/10.41 km/s; released-binding feedback gives 12.03/5.75, versus ordinary-only 52.57/62.34. No Milky Way speed fit was performed. Disk scale and infrared luminosity remain conditional proxies; all sensitivity cases are retained. Individual star orbits, bar/vertical dynamics and photon supply remain unresolved.

- 07:15 · `b17ac705` · from `research_plan/START-HERE.md`

**Companion-field feedback tested:** [Coupled exchange closure](research_work/results/isotropic-galaxy-transfer/feedback-report.md) includes the current deposits in the binding potential and density penalty. Frozen parameters worsen; shared refits recover about 20.04-point retained-geometry omitted error, essentially the prior 20.01 benchmark. Tested fixed points agree from opposite initial seeds. Traveling-energy scenarios are explicit Newtonian sensitivities, not a measured radiation field or complete relativistic evolution. All goals remain open.

- 06:55 · `a3f6e39f` · from `research_plan/START-HERE.md`

**Compact/extended exchange tested:** [Three local-rate rules](research_work/results/isotropic-galaxy-transfer/reservoir-exchange-report.md) preserve one-third capture. The released-binding candidate predicts shared-fit compact fractions of 30.9-37.6% from initial stellar conditions and gives 20.15-point omitted-target RMS versus 20.01 for prior partial migration. It is a possible local rationale, not an accuracy gain or self-consistent equilibrium. Extended-target deficits remain; no new lens/motion success is claimed.

- 06:47 · `a606f52c` · from `research_plan/START-HERE.md`

**Wave/particle support compared:** [Equilibrium candidates](research_work/results/isotropic-galaxy-transfer/wave-particle-report.md) retain the original inventory and capture law. Gaussian wave, repulsive wave and circular-particle omitted-target RMS are 39.08, 33.25 and 25.78 points under retained geometry, versus 20.01 for earlier partial migration. Gaussian radial energy minima and individual circular-orbit stability pass, but full stability and formation are not established. No candidate improves the existing profile benchmark or earns a new lens/motion success claim.

- 06:38 · `787b905d` · from `research_plan/START-HERE.md`

**Saturation tested:** [Three bound-capacity rules](research_work/results/isotropic-galaxy-transfer/saturation-report.md) preserve the one-third capture law and count unbound overflow. Constant, stellar-density and baryonic-gravity capacities do not beat prior partial migration in retained-geometry transfer. Best distant shared fits have inactive saturation. A nearby-source gravity-linked case helps its own baseline but remains worse than partial migration. Some allocations depend strongly on packet order; support, source history and the traveling component gravity remain unsolved.

- 06:31 · `07203be8` · from `research_plan/START-HERE.md`

**Four inward-migration variants tested:** [Results](research_work/results/isotropic-galaxy-transfer/migration-variants-report.md) keep capture and one-third retention fixed. Partial mobility performs best: retained-geometry distant-input cumulative RMS changes from 26.74 to 20.01 percentage points in omitted-galaxy transfer. It improves concentrated targets but worsens three others; standard-geometry distant-input transfer does not improve. The fitted 35.8% mobile fraction is a halo-shape diagnostic, not the earlier rotation-trained fraction. Support, gravitational work and absolute supply remain unresolved.

- 06:26 · `451ad23d` · from `research_plan/START-HERE.md`

**Straight-stream capture executed:** [Fixed-law radial test](research_work/results/isotropic-galaxy-transfer/stream-capture-report.md) traces nearby and distant external inputs through the original capture profile within diagnostic 5 Re. Target-specific distance mixtures improve three extended standard-geometry halo shapes but leave central concentration deficits in three others; none improves the retained-geometry targets in this basis. Improved cases favor the nearest tested source shell. This is normalized placement, not verified source populations, absolute supply or new lensing success.

- 06:22 · `21abe0be` · from `research_plan/START-HERE.md`

**Near/far companion streams:** [Source-flow calculations](research_work/results/isotropic-galaxy-transfer/companion-streams-report.md) retain the exact one-third law. Source distance relative to halo radius controls illumination gradients; a uniform distant population can collectively outweigh nearby sources. The cached M87 catalog supplies directional and distance-bin diagnostics only, extending to about 72 Mpc. No billion-year source history, steering solution, or new halo fit is established. All six goals remain open.

- 06:11 · `4b5564e4` · from `research_plan/START-HERE.md`

**Halo configuration and deposition map:** [Exact fitted parameters and finite-region targets](research_work/results/isotropic-galaxy-transfer/halo-deposition-map-report.md) give r_s, rho_s, stellar masses, boundary flags and radial comparisons with the fixed one-third deposit inventory. Four standard-geometry targets inside 5 Re have enough fitted inventory in principle; two do not. This is conditional mass bookkeeping, not demonstrated capture, realistic halo inference or a verified photon-energy budget.

- 06:06 · `80db170a` · from `research_plan/START-HERE.md`

**Same-data NFW/geometry comparison complete:** [Results](research_work/results/isotropic-galaxy-transfer/nfw-geometry-report.md) fit stars plus NFW to all six lens galaxies under our optical distances and a separately labeled standard FLRW benchmark. Standard geometry improves every NFW stellar score, but J1402 remains mismatched and parameter boundaries persist. This identifies an optical-geometry contribution to the current tension; no universal light-bending multiplier or expansion premise was adopted for our model.

- 06:00 · `d01b2c3c` · from `research_plan/START-HERE.md`

**Lens/profile compatibility analyzed:** [Inverse diagnostic](research_work/results/isotropic-galaxy-transfer/lens-profile-compatibility-report.md) allows per-target mixtures of seven mass-conserving dilations while keeping exact-one-third inventory fixed. Enforcing catalog lens angles leaves severe motion mismatches in J0037, J1204 and J1402; J1621 admits a low-residual fit at the orbital bound. Several profile scales reach the basis edge. This constrains the tested fixed-input family, not all possible redistribution. Geometry/stellar-orbit mapping remains a concrete next diagnostic; no bending multiplier was adopted.

- 00:30 · `b26caeb6` · from `research_plan/START-HERE.md`

**One-third redistribution iteration complete:** [Four-variant results](research_work/results/isotropic-galaxy-transfer/redistribution-report.md). Keeping capture constants and eta fixed, a compact 0.37% component reduces validation/test speed RMS to 31.67/22.59 km/s; lens errors remain near 14%. It is a two-parameter rotation improvement candidate, not a joint solution. Shared spreading and retention-conditioned dilation do not help materially. The compact follow-up has no bound hit, and all lens transfers used training-frozen parameters.

- 00:22 · `54b30063` · from `research_plan/START-HERE.md`

**Running paper assessment and new test:** [Cross-scale performance](papers/cumulative-time-companions/cross-scale-performance.md) records successes and failures with formula provenance. [Exact-one-third lens transfer](research_work/results/isotropic-galaxy-transfer/cross-test-audit-report.md) retains about 14% lens errors. Paired galaxy speed errors improve over baryons for 56/60 comparison galaxies but beat tested MOND for only 20/60. No cluster normalization-transfer success is claimed.

- 00:12 · `b7516552` · from `research_plan/START-HERE.md`

**NFW cluster follow-up:** [Full-profile and omitted-bin Coma checks](research_work/results/isotropic-galaxy-transfer/cluster-comparison-detail-report.md) reproduce the published NFW fit approximately (3.855 versus 3.87). NFW and transparent companions are essentially tied in omitted-bin prediction (6.82/6.86); point mass deteriorates to 32.17 and is only a control. The initial outer-only ranking is not robust evidence for companions. [Next sample specification](research_work/results/isotropic-galaxy-transfer/cluster-comparison-next-sample.md).

- 00:08 · `687f2001` · from `research_plan/START-HERE.md`

**Cluster comparison:** [Coma results](research_work/results/isotropic-galaxy-transfer/cluster-model-comparison-report.md) compare companion, NFW and compact-baryon MOND shear shapes using the same inner/outer bins. Several shapes remain viable; this does not yet test the frozen one-third normalization or multiple clusters. See the linked protocol for the physical inputs needed to do that.

- 00:04 · `61b21b9a` · from `research_plan/START-HERE.md`

**Matched MOND/dark-halo comparison:** [Rotation benchmarks](research_work/results/isotropic-galaxy-transfer/model-comparison-report.md) show simple MOND outperforming our exact-one-third candidate on frozen galaxies with fewer shared parameters. A restricted shared NFW mapping performs worse and reaches a scale bound; target-inner-fitted NFW gives much better outer predictions with extra target information. No full-theory ranking or dark-matter exclusion is claimed.

## 2026-09-12

- 23:39 · `82f24792` · from `research_plan/START-HERE.md`

**Improved galaxy candidate:** [Bounded radiation retention](research_work/results/isotropic-galaxy-transfer/bounded-retention-report.md) changes radiation's role from interception to storage probability eta=X^q/(1+X^q). Training-only q=0.33999 improves both frozen log and km/s errors; depth-dependent retention stays near constant. This candidate still needs motion/lensing transfer and a Solar System source bridge; it is not a unified solution.

- 23:36 · `2c05c483` · from `research_plan/START-HERE.md`

**Capture-law variety evaluated:** [Solar System family comparison](research_work/results/isotropic-galaxy-transfer/capture-law-family-report.md) finds no unified full-mass scaling among six formulas (two equivalent). [Radiation-conditioned galaxy capture](research_work/results/isotropic-galaxy-transfer/radiative-flux-report.md) improves frozen fractional errors but worsens km/s errors; it remains a mixed candidate rather than replacing the reference. Full-mass and extra-gravity targets remain explicitly distinct.

- 23:24 · `aca766d7` · from `research_plan/START-HERE.md`

**Earth conversion source decision:** [Outgoing EM budget](research_plan/earth-em-conversion-budget.md) uses already emitted terrestrial IR with the retained alpha, not a new ground-state electron drain. It predicts about 0.373 W converted across all outgoing directions before lunar-orbit distance; local retention is separate. Static source energy is distinguished from a hypothetical decaying well requiring replenishment.

- 23:01 · `cbb62efb` · from `research_plan/START-HERE.md`

**Common planetary well law:** [Energy to potential across Earth, Moon and planets](research_plan/planetary-well-energy-law.md) uses W=G(M+E_retained/c^2)/R and derives solar-fed growth. The same solar history/common geometric capture predicts Earth/Moon added-depth ratio 3.67, versus existing-depth ratio 22.28; it cannot explain their entire gravity on those assumptions. Ordinary-gravity plus a small retained-energy increment remains the stated conditional interpretation.

- 22:58 · `2c85c43f` · from `research_plan/START-HERE.md`

**Solar-system specialization:** [Photon-fed well formulas](research_plan/solar-system-photon-feeding.md) calculates solar photon output, conversion at the retained alpha, capture-dependent stored energy and its gravitational effect. Perfect capture within 100 AU bounds solar-generated deposits at about 16,219 kg/year; the solar source loses the emitted energy. No graviton count or capture efficiency is inferred from force strength.

- 22:55 · `c89de20e` · from `README.md`

**Theory extension:** [Persistent microwave background and collective modes](research_plan/persistent-microwave-background.md) and [academic addendum](papers/cumulative-time-companions/persistent-background-addendum.md). Three alternative mechanisms are specified for testing, with optional spherical geometry and no fixed cosmic age or radius. No CMB match or energy-supply solution is claimed. See the [current checkpoint](research_work/results/RESEARCH-CHECKPOINT.md) for executed results.

- 22:55 · `c89de20e` · from `research_plan/START-HERE.md`

**New theory extension:** [Persistent microwave background and collective spatial modes](research_plan/persistent-microwave-background.md) adds recording-transition, driven-resonance and spatial-pattern alternatives, with optional boundary-free spherical geometry. Known mode mathematics is distinguished from proposed physics. No radius/age is fixed, no CMB fit is claimed, and energy supply/retention constraints remain. Use the [current checkpoint](research_work/results/RESEARCH-CHECKPOINT.md) for the latest executed results.

- 17:23 · `eb202a82` · from `README.md`

Latest predictive diagnostic: [outer stellar motions with mass and orbit uncertainty](research_work/results/slacs-outer-predictive/report.md). Three of six exposed outer measurements remain outside the extra-force model's conditional 95% intervals (four for the baryonic baseline), under both mass priors. All 24 integration checks pass. Fixed structural assumptions remain; this is not the reserved independent prediction test.

- 17:23 · `eb202a82` · from `README.md`

Previous predictive diagnostic: [outer radial bins omitted from fits](research_work/results/slacs-outer-bin-check/report.md). Best-fit-only predictions motivated the parameter-uncertainty calculation above.

- 17:16 · `172f4771` · from `README.md`

Latest orbit-boundary check: [six resolved training lenses and comparison chart](research_work/results/slacs-orbit-boundary-extension/report.md). Wider tangential range removes baseline boundary fits; radial chi-square is 152.04 baseline versus 85.87 extra force, with extra force better in five of six. Significant discrepancies remain; no causal validation.

- 17:13 · `5589f379` · from `README.md`

Latest profile refit: [six training lenses with fixed published stellar components](research_work/results/slacs-component-refit/report.md). Extra-force radial chi-square 95.86 to 85.87; lensing RMS 10.60% to 11.60%. Four baseline orbit-boundary fits retained. No unified solution or reserved score claimed.

- 17:10 · `1b017cb3` · from `README.md`

Latest stellar-profile audit: [six published component models reconstructed](research_work/results/slacs-light-profile-audit/report.md). Half-light radii match release within 0.029%; one is 52% larger than our older pilot input. Seventh profile missing. No gravity refit yet; preserve this uncertainty in interpreting earlier results.

- 17:08 · `d1152e10` · from `README.md`

Latest resolved fit: [seven training lenses, 47 radial bins](research_work/results/slacs-resolved-fit/report.md). Fixed extra force lowers summed radial chi-square 205.74 to 111.22, but residuals remain and lensing RMS does not improve. Training-only conditional result; all six objectives open.

- 17:05 · `877721d9` · from `README.md`

Latest resolved-input audit: [54 training radial measurements](research_work/results/slacs-resolved-input-audit/report.md). All covariance matrices numerically valid; correlations reach 0.843. Release PSF FWHM is 0.8 arcsec; J0330-0020 use flag is zero and remains unresolved before fitting. No reserved score opened.

- 17:04 · `40b27e55` · from `README.md`

Latest data acquisition: [KCWI radial profiles and covariances for 14 SLACS lenses](research_work/results/slacs-resolved-data/report.md). Eight overlap training; three validation, two test and one unassigned retained. Pinned source hashes verified; no resolved fit or reserved evaluation performed.

- 17:02 · `58da91aa` · from `README.md`

Latest lensing sensitivity: [common orbital anisotropy on 33 training lenses](research_work/results/slacs-orbit-sensitivity/report.md). Median extra-force angle excess ranges 3.7-15.5%; orbital constraints are needed before a gravity-law conclusion. No beta selected from these scores, reserved roles untouched.

- 17:00 · `5c70481e` · from `README.md`

Latest observational transfer: [33 SLACS training lenses](research_work/results/slacs-motion-lensing-pilot/report.md). Mass inferred from stellar dispersion alone; fixed extra-force template predicts median lens angles 9-10% high and does not improve descriptive RMS. Conditional spherical empirical pilot, not companion-source validation; reserved scores untouched.

- 16:57 · `b06c6476` · from `README.md`

Latest numerical evidence: [641-node timing refinement](research_work/results/timing-scatter-grid-refinement/report.md). Zero-scatter b changes by 0.0019 and 0.0023 on two exposed artificial samples; revised uncertainty calibration remains necessary. Original 160-case experiment unchanged.

- 16:54 · `02f5ad5a` · from `README.md`

Latest timing development: [continuous source-scatter integration](research_work/results/timing-continuous-scatter/report.md). Exact integration includes zero scatter; exposed-case grid and optimizer limitations remain. Separate revision only: frozen 160-case calibration unchanged, all six demonstrations open.

- 16:52 · `7970edf8` · from `README.md`

Latest source-clock evidence: [SN 2011fe observed tail audit](research_work/results/sn2011fe-clock-audit/report.md). All 14 published LBT R observations retained; a single exponential has diagonal chi-square 1142/12. Its naive apparent stretch is not a propagation measurement. Multiband source physics remains necessary; all six demonstrations stay open.

- 16:48 · `377eb8a3` · from `README.md`

**Physical source-clock option assessed:** [Radioactive-tail calculation](research_work/results/radioactive-clock-feasibility/report.md) identifies a possible nuclear-decay anchor but shows that ignoring deposition evolution biases the inferred stretch. A conditional slope bound is derived and checked; real bolometric coverage and source-physics constraints remain necessary. No calibrated independent clock or new astronomical success is claimed.

- 16:46 · `c5fc6f40` · from `README.md`

**Joint pipeline recovery checked:** [Two artificial end-to-end controls](research_work/results/joint-light-recovery/report.md) recover b=-0.021/0.982 for injected 0/1, together with a shared rate, luminosity and duration. Independent analytic injection agrees with numerical real-filter predictions. Identical sources and exact distances are restrictive control assumptions; this does not validate the theory or override the broader timing-calibration failures.

- 16:43 · `2e72c01c` · from `README.md`

**Joint flux likelihood checked:** [Shared-calibration likelihood](research_work/results/joint-light-likelihood/report.md) retains model-dependent covariance and normalization, with independent dense and scalar-integration checks. It must be evaluated jointly across events sharing zero points. Source, selection, nonlinear-systematic and redshift-likelihood components remain required; no observed fit or theory validation is claimed.

- 16:42 · `9f6729a6` · from `README.md`

**Joint light forward calculator assembled:** [Verified distance/SED-to-detector calculation](research_work/results/joint-light-forward/report.md) predicts redshift, time stretch and calibrated DES band flux from one conditional transport rule. Absolute-reference, inverse-square and 24 pulse checks pass. The source model, calibration likelihood and observational validation remain required; no causal or astronomical success is claimed.

- 16:40 · `e7aadc50` · from `README.md`

**DES light-curve versions match:** [Full paired release audit](research_work/results/des-release-compatibility/report.md) finds exact agreement in all 1,779,030 aligned measurement rows for 19,706 events, including flux and errors. Later host-galaxy metadata differ and remain pinned to the corrected release. This closes measurement-version ambiguity for the historical calibration workflow; source and remaining systematic requirements stay open.

- 16:35 · `006d91c0` · from `README.md`

**Measured zero-point covariance acquired:** [Verified Fragilistic DES5YR block](research_work/results/des-calibration-covariance/report.md) provides correlated g/r/i/z uncertainties of about 0.0057-0.0061 mag. Full cross-survey covariance is retained; this is only the zero-point component, not the entire brightness error budget. Compatibility, passband/source-systematic treatment and physical source constraints remain open.

- 16:34 · `be584249` · from `README.md`

**Historical DES offset convention checked:** [Source/FITS audit](research_work/results/des-calibration-convention/report.md) establishes the synthetic-model offset sign and implements native/AB-equivalent flux transforms with shared-band covariance propagation. Four reference checks and a finite-difference covariance check pass. Actual observations remain untouched; release compatibility and measured covariance are still required.

- 16:31 · `f0d8ffad` · from `README.md`

**Real-filter timing bridge checked:** [Controlled spectral-duration experiment](research_work/results/bandpass-timing-coupling/report.md) shows the fixed-filter width exponent is b-q for the declared self-similar source; matching emitted response recovers b. All 320 numerical cases match their analytic scalings. This is not a real-supernova result; it defines a necessary spectral/passband component of the joint-light inference.

- 16:29 · `8d98788b` · from `README.md`

**Faint timing boundary diagnosed:** [First two flagged controls](research_work/results/timing-boundary-diagnosis/report.md) converged but reached the imposed intrinsic-scatter floor 0.03. Timing estimates remain near injected values; nominal truth inclusion does not override the frozen validity screen. The zero-invalid-fits criterion has already failed in at least one cell. The full batch continues unchanged, and uncertainty treatment needs further work before observational inference.

- 16:28 · `39d86fbf` · from `README.md`

**DES calibration acquired:** [Checksum-verified original resources](research_work/results/des-photometric-calibration/report.md) supply the four real passbands and calibration inputs. Independent AB-reference photon integrals agree with the archived program log within 0.0114%. Version compatibility, offset application, covariance and source-evolution constraints remain unresolved; no observed brightness fit is claimed.

- 16:23 · `1d0aa6af` · from `README.md`

**Joint-light input requirement derived:** [Executed identifiability audit](research_work/results/joint-light-identifiability/report.md) shows that free duration and peak-luminosity evolution can exactly offset a changed timing exponent, even with fixed distances; fluence does not add an independent constraint in the fixed-shape bolometric model. Source-evolution constraints and the DES external photometric calibration are needed for a predictive joint test. This is a conditional deduction, not an observational failure or completed mechanism.

- 16:20 · `0b753e3d` · from `README.md`

**Expanded timing calibration launched:** [Frozen 160-case batch](research_work/results/timing-coverage-calibration/README.md) spans two source shapes, two signal levels and both timing laws, with 20 seeds per cell and eight later numerical refinements. Completed cases are hash-checked and checkpointed; the execution session, not this historical note, determines whether the job is live. Scientific interpretation requires the completed report and full per-cell failures/uncertainties. All six demonstrations remain open.

- 16:18 · `a6e90725` · from `README.md`

**Timing uncertainty method:** [Nominal profile intervals](research_work/results/timing-profile-uncertainty/report.md) contain the injected truth in 6/6 exposed base controls at nominal 95%, and 3/6 at nominal 68.27%. The controls have only three paired seeds, so these counts do not establish coverage. Fixed-b optimization and interval-root checks pass; broader calibration and joint physical brightness remain outstanding. All six demonstrations remain open.

- 16:15 · `77652e76` · from `README.md`

**Timing null/stretch controls:** [Six artificial-data controls](research_work/results/timing-null-stretch-controls/report.md) pass the frozen small-control test: median recovered exponent -0.006 for injected 0 and 1.022 for injected 1. Refining the worst-error case changes the exponent by 0.00134. Every base estimate is slightly above its true sample slope, so residual bias remains to investigate. This is not observed time dilation, calibrated interval coverage or a joint brightness result. All six scientific demonstrations remain open.

- 16:08 · `1b419d5e` · from `README.md`

**Joint-light timing prerequisite:** [Refined synthetic-cadence calculation](research_work/results/timing-quadrature-refinement/report.md) passes the unchanged mean-based numerical gate: shape integration changes 0.0531/0.0356 (limit 0.1), and the separate duration-grid comparison changes 0.0346. Some individual curves exceed 0.1 and remain reported. This is not a real-supernova result or a brightness test; repeated-trial bias/coverage, source/filter/selection calibration and joint physical flux predictions remain required. The six-demonstration evidence program remains active.

- 14:11 · `4b0497b6` · from `README.md`

**Age and size constraints released:** [Executed assumption audit](research_work/results/assumption-release/report.md) replaces fixed cosmic age/size with free histories and analytic infinite limits. The old ideal full-capture energy benchmark requires a median 54.3 trillion years at constant present luminosity; this is not a measured age or a funded history. Infinite lossless traveling supply and permanent capture have separate convergence/growth problems. Conditional checks pass; no complete mechanism or improved observational fit is claimed.

- 13:59 · `01f44d16` · from `README.md`

**Adaptive first force comparison complete:** [Fixed-parent-mass audit](research_work/results/bar-cell-geometry/adaptive-report.md) completes all 218 new trajectories and passes all five age-layer checks. Four of five spatial force gates pass, but the central target still changes by 32.2%, exceeding 5%. Source mass and sampled roundoff checks pass; geometry, capture-quadrature sensitivity and other failures remain unresolved. No physical/observational closure follows; all nine goals remain open.

- 13:43 · `ea9e1ebc` · from `README.md`

**Angular refinement targets identified:** [Complete contribution audit](research_work/results/bar-cell-geometry/angular-attribution-report.md) reproduces the problematic age-bin force and finds that 66/560 parent regions supply 90% of summed difference magnitudes. All directions and their masses remain included. The ranking guides actual refinement; it does not resolve the 14 force failures or establish physical/observational success. All nine goals remain open.

- 13:42 · `4d27e008` · from `README.md`

**Largest force discrepancy localized in age:** [Mass-preserving attribution](research_work/results/bar-cell-geometry/age-attribution-report.md) reproduces both archived forces and traces 92.3% of their signed difference to the 61–73 Myr source-age bin at the selected failing target. No ages are removed or renormalized. Next resolve angular structure during that passage; other failures and all nine physical/observational goals remain open.

- 13:40 · `685949eb` · from `README.md`

**Curved interpolation does not resolve the late map:** [Quadratic audit](research_work/results/bar-cell-geometry/quadratic-report.md) uses the same actual vertex/midpoint paths with known six-node shape functions. Early RMS improves, but final errors worsen to 25.63% and 40.52% of launch radius. Retain this failed remedy; no new force result or mass reweighting is adopted. The model still needs resolved late-time source geometry and all physical/observational closures. All nine goals remain open.

- 13:39 · `1ede0057` · from `README.md`

**Second full-sphere gravity comparison complete:** [Completed audit](research_work/results/bar-cell-geometry/nested2-report.md) passes all 30 age-layer potential and force gates but fails 6/30 spatial potential and 14/30 spatial force gates; the largest force change is 135.2%. Source mass and sampled high-precision checks pass, while the source map remains unconverged. Next address curved/folded cells without physical retuning. No stellar/lensing fit, energy closure or full-model validation is established; all nine goals remain open.

- 13:05 · `4c0b1f41` · from `README.md`

**Spatial coverage gap confirmed:** [New interior-trajectory audit](research_work/results/bar-cell-geometry/local-refinement-report.md) finds that prior probes occupy cells containing only about 24% of source mass. Eight actual trajectories in previously unchecked cells have final interpolation errors of 23–78% of launch radius despite passing orbit checks. Errors also span many already sampled cells. Next refine across the sphere with actual trajectories, preserving mass and force-convergence requirements; do not infer that unchecked cells pass. All nine goals remain incomplete.

- 13:03 · `e8543bca` · from `README.md`

**Nested mesh improves geometry but remains unconverged:** [Completed nested comparison](research_work/results/bar-cell-geometry/nested-report.md) adds 160 actual trajectories and reduces final interpolation RMS from 66.7/70.2 to 34.6/49.7 percent of launch radius. Mass is preserved; a 1024-layer supplement resolves the remaining age gate, but matched 512 angular comparisons fail 15/30 potential and 14/30 force gates. Nine of 10 stricter orbit checks pass in the original domain; the polar cache-boundary failure is retained, with separate shrinking-core sensitivity tests passing at orbit tolerances. Next refine the late source map before stellar predictions. Physical supply, self-gravity, redshift/lensing and all nine goals remain incomplete.

- 13:03 · `e8543bca` · from `README.md`

**Source-cell geometry needs refinement:** [Direct trajectory interpolation audit](research_work/results/bar-cell-geometry/report.md) finds final RMS errors of66.7/70.2 percent of launch radius when source6 cells predict actual source8 trajectories. This prevents treating the coarse affine density as resolved, even with exact cell gravity. The nested mesh has now been evaluated; no refined-force success is claimed. All nine goals remain incomplete.

- 12:44 · `a83f829c` · from `README.md`

**Mass-preserving volume gravity implemented:** [Trajectory-volume audit](research_work/results/bar-volume-gravity/report.md) replaces point-stream source quadrature with exact unsoftened tetrahedral gravity and an affine source-age mesh. Cube, independent exterior, thin-cell, cache and sampled60-digit checks pass; normalized mass error is below1.7e-15. At256/512 age layers,59/60 force gates pass; source6/source8 still fails15/30 force and8/30 potential gates. Angular interpolation/folds and one late central age check remain unresolved. The method changes numerical representation, not the force law, and is not yet a reliable stellar prediction. Next refine the curved source geometry while preserving mass, then close physical supply/self-gravity and the remaining redshift/lensing requirements. All nine goals remain incomplete.

- 12:32 · `4ab7c9bd` · from `README.md`

**Finer force grid still fails; integration remedy benchmarked:** [Source8/source12 comparison](research_work/results/bar-population-gravity/source-refinement12.md) completes both partitions and passes all208 combined orbit/age status checks, but2/30 potential and19/30 force gates fail. [Unsoftened straight-stream benchmark](research_work/results/stream-force-quadrature/report.md) shows a direct96x96 grid with0.00343 percent potential error still has8.69 percent force error; exact radial integration plus angular quadrature passes reference/gradient checks. Next adapt close-passage integration to curved trajectories instead of relying on larger uniform source grids. No modified force law, fit, physical supply closure or collective model is established; all nine goals remain incomplete.

- 12:32 · `4ab7c9bd` · from `README.md`

**Moving population now produces three-dimensional gravity:** [Linear-response coefficients](research_work/results/bar-population-gravity/report.md) integrates the source-age distribution into potential and acceleration at five positions in and above the plane. Analytic ring/Kepler checks, all36 orbit checks and108 age-quadrature gates pass. These are coefficients per unit injected rest mass in a fixed ordinary bar; source-direction convergence, physical normalization, self-gravity, ordinary-matter response and lensing remain unverified. Both angular runs are now complete; their later comparisons are recorded above. All nine goals remain incomplete.

- 12:24 · `d3a28fec` · from `README.md`

**Source refinement exposes force uncertainty:** [Angular population refinement](research_work/results/bar-angular-refinement/report.md) completes64 and144 new trajectories, with all individual numerical checks passing but10/36 encounter/residence comparison gates failing; final residence changes are19.6 and33.2 percent relative to the finer values. [Gravity source refinement](research_work/results/bar-population-gravity/source-refinement.md) passes all100 combined orbit/age status checks and all30 potential comparisons within2 percent, but21/30 vector-force comparisons fail5 percent. Well depth is better resolved than local acceleration. Next refine source integration for forces before a stellar/bulge fit or self-gravity evolution. Physical supply, joint redshift, lensing and all nine goals remain incomplete.

- 12:10 · `35c6615a` · from `README.md`

**Bar source and repeated passages:** [Central residence audit](research_work/results/bar-source-persistence/report.md) follows source-weighted angular grids at 1 and 3 kpc through about 244 Myr. Later central entries occur, while actual central residence is much smaller than the fraction ever visiting. All 88 orbit checks and field-order comparisons pass after repairing missed short crossings; 10 angular encounter gates fail. At 1 kpc the residence estimates differ by a factor of 2.5 across angular grids, despite passing the loose absolute gate. Source precision, self-consistent formation and all nine goals remain incomplete.

- 11:54 · `14afc973` · from `README.md`

**Nonspherical-torque checkpoint:** [Ordinary bar first-pass test](research_work/results/bar-capture-torques/report.md) generates angular momentum from purely radial probes without an extra force law. Eight selected probes enter0.1kpc in the spherical control; all turn outside it on the first approach in the bar cases. Energy/Jacobi and field-order checks pass. This conditional bar-component experiment is not a formed population or a measured capture fraction.

- 11:48 · `9c485fde` · from `README.md`

**Coupled formation checkpoint:** [Moving deposits with capture feedback](research_work/results/coupled-radial-formation/report.md) starts from zero deposits and evolves their gravity and depleted supply together. The lower rate passes the tested slow-motion numerical checks; the higher rate leaves the approximation. Central potential remains sensitive to numerical core size even where outer mass fractions are stable.

- 11:41 · `ecabf28a` · from `README.md`

**Moving-deposit checkpoint:** [Radial capture response](research_work/results/radial-capture-response/report.md) follows a negligible bound source in the frozen well. The innermost added-mass fraction grows about96-fold relative to in-place capture, while the projected fraction grows10.34-fold. Orbital energy is conserved, but the source does not build the assumed stationary profile. Changing self-gravity and source feedback remain necessary.

- 11:36 · `b18d5204` · from `README.md`

**Capture-momentum checkpoint:** [Newly captured cohort audit](research_work/results/capture-injection/report.md) derives rest mass and inward motion for one local merger rule. Weak absorption can inject bound radial material; stronger absorption mostly injects unbound material at the tested scales. The rule creates no circular angular momentum, so snapshot support is not a formation solution.

- 11:32 · `80aa37ee` · from `README.md`

**Retention checkpoint:** [Bound-particle support audit](research_work/results/deposit-retention/report.md) rules out isotropic bound-particle support for a strongly shielded deposit profile. Randomly oriented circular orbits provide a conditional stationary alternative in the frozen spherical potential, with kinetic and potential energies calculated. Formation, collective stability and continuing capture remain unresolved.

- 11:27 · `a5904e66` · from `README.md`

**Depleted-supply checkpoint:** [Attenuated companions with gravity feedback](research_work/results/depleted-capture-feedback/report.md) separates incident fluence from opacity and removes captured energy from onward rays. The absorption ledger balances, but several models remain boundary-dependent; unresolved numerical gates are retained explicitly. Mechanical retention and a photon-funded source still require closure.

- 11:13 · `50ed8802` · from `README.md`

**Integrable-tail checkpoint:** [Steeper capture with feedback](research_work/results/integrable-capture-feedback/report.md) finds finite-range mass convergence in weak-source p4/p6 examples, but domain-filling capture at higher supply. Inner force convergence alone misses the divergent mass and depth. Supply depletion and retention dynamics remain necessary before another galaxy calibration.

- 11:10 · `a23f39a3` · from `README.md`

**Capture-boundary checkpoint:** [Boundary and finite-mass audit](research_work/results/coupled-capture-boundary/report.md) finds that refitting the amplitude masks a more than 400-fold change in deposited mass across assumed boundaries. The best tested boundary still misses the radial shape. A physical outer capture/retention law is required before further tuning.

- 11:05 · `e7dd2ed5` · from `README.md`

**Coupled-calibration checkpoint:** [Feedback-aware amplitude fit](research_work/results/coupled-capture-calibration/report.md) reduces the large overprediction, but both capture-law variants remain too weak inside and too strong outside. The full curves and required masses show that spatial/retention assumptions need testing beyond an overall amplitude adjustment.

- 11:02 · `d8bac559` · from `README.md`

**Capture-feedback checkpoint:** [Chronological deposited-gravity test](research_work/results/milky-way-capture-feedback/report.md) shows that the previously fitted capture amount overproduces stellar speeds once deposits deepen the capture well. The coupled growth equations must be calibrated directly; the earlier no-feedback agreement is not a validated complete model.

- 10:58 · `9bf04f46` · from `README.md`

**Milky Way calibration checkpoint:** [Depth-dependent capture against training stars](research_work/results/milky-way-depth-capture/report.md) substantially improves radial agreement with one fitted amplitude, but requires a deposit potential large enough to change the capture rule itself. This is a conditional calibration, not a new-data validation or a complete photon-to-gravity mechanism.

- 10:53 · `8ea467dc` · from `README.md`

**Capture-law checkpoint:** [Depth-law comparison](research_work/results/capture-depth-law/report.md) finds a nearly flat finite-range force contribution for one thin-capture case, but the same profile declines farther out. Full curves, interception efficiency and lensing remain to be checked against actual galaxy geometry and data.

- 10:50 · `44182bac` · from `README.md`

**All-direction capture checkpoint:** [Isotropic deposition profile](research_work/results/isotropic-capture-profile/report.md) computes where companions are intercepted, then derives gravity and lensing. Strong capture can load outer layers while starving the center; that pattern does not automatically fit stellar motion. Supply and nonspherical tests remain open.

- 10:47 · `a3f617d2` · from `README.md`

**Deposit/gravity checkpoint:** [Growing-well clock and lensing test](research_work/results/deposit-clock-gravity/report.md) links a prescribed deposit profile to stellar forces, clocks and light bending. Growing wells can redshift crossing signals, through a known evolving-potential effect, but the strength is tied to the gravitational potential and its growth. Capture-funded profiles and full energy accounting remain unresolved.

- 10:44 · `74923192` · from `README.md`

**Response-history checkpoint:** [Slow-response and steady-flow comparison](research_work/results/slow-response-history/report.md) shows that a longer response time does not sustain stretching in the tested open channel. Independent stationary equations reproduce the late results. Deposits grow, but their gravity/clock feedback remains to be added and tested.

- 10:41 · `0fff7d1b` · from `README.md`

**Initiation checkpoint:** [Radiation-powered startup](research_work/results/radiation-powered-initiation/report.md) removes the external seed and funds deposits, but mainly by removing whole photons. Early redshift does not persist after the field settles. Surviving-light frequency, timing and attenuation are reported separately in the [working paper addendum](papers/cumulative-time-companions/temporal-field-addendum.md).

- 10:38 · `f1725ab8` · from `README.md`

**Spatial energy-assignment checkpoint:** [Field-funded transport](research_work/results/field-funded-transport/report.md) keeps reservoirs nonnegative and energy moving in a tested finite domain, with consistent measured clock/frequency changes. It still gives path-dependent redshift/blueshift and net energy gain by radiation; its deposits are funded by the supplied seed overall. A common radiation-powered history remains to be derived.

- 10:34 · `cfb1c6b6` · from `README.md`

**Finite-response checkpoint:** [Moving response-energy audit](research_work/results/finite-response-transport/report.md) improves local propagation in an explicit slower-field variant, but a spatial gradient can demand energy from an empty companion reservoir. Uniform-cell conservation passes do not resolve that failure. The next constitutive rule must specify spatial energy exchange as well as temporal response.

- 10:31 · `4deeac47` · from `README.md`

**Traveling field-energy checkpoint:** [Propagation and capture audit](research_work/results/advected-temporal-state/report.md) tests tying the temporal state directly to traveling companion energy. Its reference-energy budget closes, but capture reaches an unstable equal-speed boundary from initially regular examples. The tested instantaneous feedback requires revision before a joint redshift/gravity claim; energy conservation alone does not establish sensible signal propagation.

- 10:15 · `f170bbfb` · from `README.md`

**Operational clock result:** [Evolving nonuniform temporal field](research_work/results/inhomogeneous-clock-completion/report.md) retains measurable frequency/event stretching with local c and fixed spatial rods in a prescribed lapse geometry. Uniform changes cancel; generic endpoints and time-field forces matter. The radiation-fed origin of that spatial field and its gravity/energy completion remain unproved.

- 10:12 · `cbb9c35a` · from `README.md`

**Execution roadmap:** [Postulate-based analysis goals](research_plan/postulate-based-analysis-goals.md) define the sequence, purpose and completion evidence for developing and testing the two selected branches. We can postulate operational rules without explaining the ultimate origin of time; their observable consequences and internal consistency still need testing.

- 09:32 · `da10abe8` · from `README.md`

**Selected direction: shared temporal transport, plus energy-fed temporal feedback.** [Current branch specification](research_plan/shared-time-and-energy-feedback.md) makes light and gravitational waves share an evolving propagation field, and adds a receiving-energy sector that can travel onward or be captured in wells. [Initial feedback calculation](research_work/results/temporal-energy-feedback/report.md) checks energy accounting and signal stretching in four nondimensional cells. Seed dependence, strong radiation feedback, clock/local-speed consistency, tensor dynamics and the field's full energy/momentum action remain open. These are explicit hypotheses, not a completed nonexpanding theory or an observational validation.

## 2026-09-10

- 08:46 · `07f43e72` · from `README.md`

**Small transfers can reproduce the mean redshift:** [Fixed-energy transfer audit](research_work/results/fixed-energy-transfer/report.md) gives an explicit statistical rule with the existing exponential mean and calculable spectral broadening. The reused radio-line width limits its step size conditionally. This is a postulated energy process, not a derived graviton interaction; event stretching, coherence, clocks and momentum remain unresolved.

- 08:43 · `4045c03e` · from `README.md`

**Conversion is not automatically redshift:** [Direct spectral discriminator](research_work/results/conversion-versus-redshift/report.md) gives stationary photon/graviton mixing exactly the energy-transfer fraction required by the existing calculator. Surviving frequencies and event-copy separation remain unchanged. This specialization can be a source channel but cannot alone explain the desired redshift; inelastic or time-dependent extensions require separate derivation. The calculator itself remains unchanged.

- 08:38 · `4fb544dd` · from `README.md`

**Brightness-distance consistency:** [Fixed-rate exposed-data sensitivity](research_work/results/brightness-distance-consistency/report.md) finds that ideal shared-stretching dimming changes median inferred path distance by -0.714% and overall redshift RMS from 447.05 to 444.45 km/s. The large scatter remains, and training RMS worsens. Published distances are not overwritten; actual passband/calibration inference and a genuinely withheld test remain required.

- 08:33 · `baed4053` · from `README.md`

**Shared-transport photometry:** [One-rate brightness and spectrum predictions](research_work/results/shared-transport-photometry/report.md) link spectral/event stretching to bolometric dimming, photon arrival rate and calibration bias. Numerical identities pass; ordinary endpoint clocks remain a stated assumption requiring physical completion. The gamma/GW arrival lag includes source timing and does not prove unequal propagation speeds; the earlier conditional photon-only exclusion did not assume simultaneous emission.

- 08:30 · `7bff3c8f` · from `README.md`

**Extended environment geometry:** [Method-specific CF4 coordinates](research_work/results/extended-environment-geometry/report.md) provide 774 selected distance-indicator tracers out to 131.8 Mpc, but only 16 beyond 50 Mpc. Distance-error propagation and three proximity definitions quantify sparsity and distance/environment degeneracy. The numerical geometry passes its independent check; the result is not a physical void map or a redshift fit. Selection and denser-tracer calibration remain required.

- 08:27 · `02eefa65` · from `README.md`

**Independent environment coverage:** [869-galaxy catalog audit](research_work/results/independent-environment-readiness/report.md) retains 335 individual distance-indicator tracer candidates inside 11 Mpc while excluding redshift-inferred coordinates. That nominal region covers only 37.93% of the median path in the 164-group redshift sample. Missing coverage is retained as unknown; no void coefficient is fitted to a falsely empty map. Numerical proximity geometry passes its independent check, but physical density, selection and full-path coverage remain unresolved.

- 08:23 · `5efadf85` · from `README.md`

**Joint propagation test:** [Timing, spectra and messenger arrivals](research_work/results/joint-propagation-audit/report.md) distinguish candidates that share a redshift formula. With the prior conversion rate, the specified uniform photon-only evolving-time law predicts about 675,000 years of extra delay at GW170817's independently estimated host distance. Its associated gamma/GW lag is about 1.7 seconds including source timing. Shared messenger transport removes that particular relative-delay conflict but remains physically incomplete. This is a conditional exclusion of one specialization, not a general nonexpansion exclusion. [Forward test contract](research_plan/joint-observation-test-contract.md) tracks independent void paths, shared spectrum/clock/brightness predictions and genuinely withheld evaluation.

- 04:45 · `8510c09f` · from `README.md`

**Nearby test discrimination:** [Post-evaluation sensitivity](research_work/results/nearby-redshift-discrimination/report.md) finds that the eight-galaxy exponential/linear prediction differences are small compared with quoted distance-error projections. Even the independent-error, zero-motion scenario has expected log likelihood separation only 0.0513 natural-log units. These are assumed-covariance test-design calculations, not an observed significance or inferred motion model. No new holdouts opened or prior scores changed.

- 04:44 · `bbab05c9` · from `README.md`

**Eight additional redshift predictions:** [Frozen ELVES-distance/ALFALFA-velocity crosscheck](research_work/results/nearby-redshift-crosscheck/report.md), protocol commit 427fa9f, yields RMS 271.88 km/s for the existing exponential and 270.11 km/s for the fixed linear control. Seven objects lie below the original distance calibration range. No rates or individual motions were fitted; this is a conditional cross-catalog result, not certified independent validation. The eight outcomes are now exposed, with 16 previously pending targets unscored by this test; see its exposure update.

- 04:38 · `a3f783e8` · from `README.md`

**Shared-speed structural restriction:** [Conditional Hamiltonian derivation](research_work/results/shared-speed-obstruction/report.md) shows that a single canonical companion sector cannot preserve its reference energy at every momentum while sharing an evolving nondispersive speed c0/n. Offset, momentum-dependent compensation and common-clock controls illustrate the limitation. This rules out parameter-only repair under those assumptions, not all companion theories. Distinct conversion/internal-state dynamics must be explicit before adoption; no holdouts opened.

- 04:37 · `6f8629ed` · from `README.md`

**Equal-speed retention failure:** [Autonomous wave test](research_work/results/companion-retention-test/report.md) evolves the proposed logarithmic-gradient field without photons or external driving. In rolling backgrounds the perturbation retains about 76.9% or 62.5% of its energy, with the loss accounted for by background kinetic-energy gain. Grid and doubled-domain checks pass. This candidate therefore does not satisfy the intended lossless transport condition in the tested reference-energy convention; total conservation alone is insufficient. No holdouts opened.

- 04:33 · `9bf898cc` · from `README.md`

**Companion speed and energy retention:** [Speed-consistency test](research_work/results/companion-speed-consistency/report.md) shows that fixed companion speed 1 does not equal photon speed 1/n in the evolving optical field. A proposed positive logarithmic-gradient Hamiltonian passes limited small-wave equal-speed checks, but its full photon interaction is untested and equal-speed dispersion does not establish lossless travel. These remain separate requirements; no new branch is adopted as a completed theory.

- 04:30 · `0512dd5f` · from `README.md`

**Smooth-pulse feedback:** [Fixed-profile refinement](research_work/results/continuous-signal-feedback/report.md) resolves the same finite pulse with 8/16/32 integration samples, checks coincident packet splitting and doubles the field-mode resolution. The brightness-dependent redshift/duration effect persists, with conserved energy and momentum. The resolved result is not an arbitrary packet-count effect, but still lacks physical normalization and a derived detector-clock law. No astronomical values or holdouts were fitted.

- 04:26 · `0af2f262` · from `README.md`

**Finite signal feedback:** [Closed packet/timing calculation](research_work/results/finite-signal-feedback/report.md) includes the measured signal packets in the shared field's energy and momentum exchange. Faint signals approach the weak-background result, while finite loading changes carrier redshift and successive event-duration stretches. Both mode resolutions, matched source-time controls, fixed-energy color control and the analytic homogeneous clock check pass. This is a dimensionless consistency experiment; physical normalization, matter clocks, astrophysical source independence, capture and gravity remain unresolved. No holdouts opened.

- 04:22 · `e72eff31` · from `README.md`

**Bulge velocity-component control:** [Training decomposition](research_work/results/bulge-velocity-components/report.md) separates spectroscopic and proper-motion/distance contributions using one exact coordinate transform and the full variance/covariance identity. For the selected bulge-plane sample, the flag-removal sensitivity changes the proper-motion vertical spread from 64.50 to 58.49 km/s, while raw APOGEE line-of-sight spread changes from 112.27 to 111.22 km/s. Different sight-line projections prevent treating these as interchangeable gravity measurements. No associations, gravity parameters or holdout outcomes were changed.

- 04:19 · `13cfed71` · from `README.md`

**Bulge association sensitivity:** [Training position/epoch audit](research_work/results/stellar-association-window/report.md) checks 77,819 applicable 2MASS identifiers among 77,927 training stars. All 140 matches flagged above 0.5 arcsec remain flagged using downloaded actual observing dates. In the selected bulge-plane sample, 54/299 are flagged, versus 0/577 off-plane; temporarily omitting flags reverses the small descriptive vertical-dispersion contrast. This is a selection/association warning, not an adopted deletion or gravity result. Original catalogs and holdouts remain unchanged.

- 04:10 · `4f83fe53` · from `README.md`

**Reserved Cepheid validation:** [Frozen-pipeline evaluation](research_work/results/cepheid-validation/report.md) was committed at e5cf484 before inferred validation outcomes. Fixed cuts retain 167 stars; nine one-kpc bins are scored and three sparse outer bins are listed unscored. Validation RMS is 67.74 km/s for original ordinary matter, 15.36 for the primary original completion, and 12.36 for the frozen balanced-component completion. The primary is low in all nine bins; the secondary in eight. Input hashes, role separation and numerical refinement pass. This is conditional held-out-star evidence with shared Gaia/calibration assumptions, not full theory validation. These validation outcomes are now exposed; the final test role remains untouched.

- 04:05 · `f7249507` · from `README.md`

**Ordinary-component response:** [Separate stellar, gas and central masses](research_work/results/baryon-component-response/report.md) re-solves the unchanged nonlinear field for each perturbation. A balanced local proposal with stellar disks x0.962, gas x0.7 and central stars x1.3 improves training rotation RMS from 18.79 to 16.36 km/s and exposed vertical RMS from 33.31 to 29.11 units. An alternate combined-error proposal improves vertical but worsens rotation, preserving the tradeoff. The +/-30 percent box is illustrative, not an allowed mass posterior; both proposals use full refined field checks. No reserved outcomes opened, field coefficients changed, or photon mechanism validated.

- 03:59 · `485f0e93` · from `README.md`

**Shared potential-shape response:** [Exploratory two-parameter audit](research_work/results/potential-shape-response/report.md) finds that vertical stretching q=1.8465 reduces vertical RMS from 33.31 to 15.18 while leaving Cepheid rotation unchanged. Rotation fitting prefers 22.60 percent more ordinary matter, but ordinary matter alone then exceeds 27/43 vertical estimates. Extending the joint search beyond q=4 drives the chosen loss toward the infinite-height cylindrical limit (rotation RMS 7.80 km/s, vertical RMS 14.89), not a finite companion source. Finite tested shapes have no negative effective-density grid values; global source admissibility, capture and lensing remain unproved. No reserved outcomes opened; this fitted geometry is not adopted as the main theory.

- 03:53 · `44461faf` · from `README.md`

**Cepheid orbital requirements:** [Jeans-condition audit](research_work/results/cepheid-jeans-requirements/report.md) shows the frozen additional field below the nonincreasing-radial-pressure, zero-tilt speed floor in all twelve training bins, by 5.83-25.87 km/s. Thus changing a declining radial profile alone cannot repair the field under steady axisymmetry. An exploratory radial/vertical cross-motion fit in four regions does not supply the required correction, but is not a selection-corrected midplane measurement or a formal exclusion. Algebra, source integrity, and training-only checks pass; no field parameter or reserved outcome changed.

- 03:50 · `9ad59fab` · from `README.md`

**Individual Cepheid training comparison:** [Common-frame reduction](research_work/results/cepheid-common-frame/report.md) uses a published Gaia-calibrated period-Wesenheit relation and 542 training stars. With fixed approximate Jeans assumptions, twelve-bin RMS is 75.24 km/s for ordinary matter and 18.79 for the frozen additional field, which remains below every bin. A fixed-membership +/-7 percent distance-scale check shifts the inferred proxy by at most 2.13 km/s. Numerical refinement and Monte Carlo propagation checks pass. There are 433 validation and 443 test measurement/mode candidates whose inferred motions remain unopened. This is conditional training evidence, not a fresh validation or a photon-derived mechanism.

- 03:43 · `219a4632` · from `README.md`

**Cepheid individual measurements:** [Parent-catalog readiness](research_work/results/cepheid-star-readiness/report.md) archives 13,077 Gaia DR3 classical Cepheids with 57 fields and an independent archive-count check. Initial cumulative measurement flags leave 2,304 candidates; none of the parent IDs overlap the 140,407 prepared red giants. All rows remain available. Distances, spatial selection and orbital inference are pending; this is not the published 903-star reconstruction or a new gravity score.

- 03:36 · `ac194f05` · from `README.md`

**Rotation-frame diagnostic:** [Published-frame sensitivity](research_work/results/rotation-frame-sensitivity/report.md) finds a 0.153 kpc solar-radius and 4.4 km/s solar-rotation difference between the two reductions. A central-ray, unchanged-drift approximation reduces curve-to-curve RMS disagreement from 12.95 to 8.68 km/s. This is not an exact catalog correction or revised model score; all twelve differences remain positive. Original frozen Cepheid scores are unchanged. A common-frame stellar-level inference is still required.

- 03:34 · `fc38572f` · from `README.md`

**New Cepheid catalog check:** [Predeclared twelve-bin evaluation](research_work/results/cepheid-rotation-check/report.md) freezes the field and mass variants at commit f7b2980 before acquiring Feng et al. 2026 Table 1. RMS errors are 76.84 km/s for ordinary matter, 20.46 for the primary extra-field model, and 14.82 for the previously Eilers-calibrated mass. Both extra-field variants underpredict every bin; no Cepheid fitting or residual exclusions occurred. Shared Gaia and inference assumptions limit independence. This table is now exposed and cannot be reused as fresh validation after tuning.

- 03:29 · `405a2da9` · from `README.md`

**Ordinary-mass sensitivity:** [Shared normalization test](research_work/results/milky-way-mass-response/report.md) finds that the conservative completion prefers 7.78% more ordinary mass for rotation but 21.62% less for vertical pull. Rotation RMS at the vertical-selected setting worsens to 24.59 km/s; vertical RMS at the rotation-selected setting worsens to 42.33 surface-equivalent units. An independent rescaled field solve verifies the exact homogeneity shortcut. This is an exposed-data global normalization diagnostic, not full component-mass uncertainty or a fresh holdout; a single mass correction does not repair the geometric mismatch.

- 03:25 · `2ad387d9` · from `README.md`

**Conservative field-equation comparison:** [Frozen Milky Way evaluation](research_work/results/conservative-field-completion/report.md) implements a known QUMOND-style potential equation with the archived empirical coefficients. Rotation RMS improves from 67.44 to 8.61 km/s, but vertical RMS worsens from 27.99 to 33.31 surface-equivalent units; median vertical overprediction is 26.9%. Global disk, analytic-limit, refinement and exterior-domain checks are archived. These are exposed model-dependent force summaries, not raw-star holdouts or a photon-derived response. The mathematical completion alone does not supply the required field geometry.

- 03:16 · `f581da36` · from `README.md`

**Production-to-storage matching:** [Same-field kinematics](research_work/results/companion-source-matching/report.md) applies the frozen 1e-24 eV/c^2 lens mass to photon conversion. Optical two-photon production of one such scalar requires an opening angle around 1e-24 rad and yields an ultra-relativistic state, not the fitted slow source. Retaining 1 eV in slow excitations of that same mass requires occupation of order 1e24, or a different receiver. This quantifies a missing state-conversion mechanism; it is not a photon-supply budget or an exclusion of stored geometric deformation. No holdout scores changed.

- 03:13 · `390ee6b3` · from `README.md`

**Carrier-resolved conversion test:** [Outgoing-spectrum calculation](research_work/results/gauge-kinetic-spectrum/report.md) finds that a resonant scalar-photon collision transfers 3.611% of its energy into the scalar while the surviving light has a 0.2208% higher energy-weighted mean frequency. Refinement and conservation checks pass. This dimensionless interaction reshapes/depletes light rather than deriving the required uniform redshift; no observational holdouts were rescored. Spectrum, event timing and matter-clock response must follow from a common completed mechanism before its gravity fits count as photon-origin evidence.

- 03:02 · `16c472f3` · from `README.md`

**Photon-to-scalar mechanism check:** [Closed gauge-kinetic pulse calculation](research_work/results/gauge-kinetic-transfer/report.md) confirms that the existing phi F²-type coupling does not source a scalar from an isolated ideal null light pulse in its vacuum. Colliding pulses do generate scalar energy with conserved total energy/momentum, but part returns to light; the dimensionless example is not a redshift, capture or astrophysical-rate calculation. This identifies the need for a specified non-null environment or other interaction and a derived spectrum/clock response before the galaxy source can be connected to photon loss.

- 02:56 · `122d1a74` · from `README.md`

**Disk-transfer bound:** [149-galaxy necessary-mass check](research_work/results/wave-disk-transfer-bound/report.md) finds the lens fraction f=0.6 insufficient at the outermost measured radius in 147/149 SPARC galaxies, even assigning all stellar light the larger adopted mass-to-light ratio and putting all source mass inside that radius. Increasing stellar masses by 1.645 still leaves 134 failures. Median necessary source fractions are at least 5.58 and 3.04 respectively. This conditional spherical-Newtonian bound shows the lens normalization cannot be universal; it is not a photon-supply calculation or a fresh SPARC holdout. A derived environment/history-dependent accumulation law or a consistently changed field response is now a central missing requirement.

- 02:51 · `750c2749` · from `README.md`

**First final-role test completed:** [Frozen 12-system evaluation](research_work/results/wave-lens-test/report.md) gives wave dispersion RMS 36.14 km/s versus 50.27 for the shared stellar-mass benchmark, but worse lensing RMS, 0.2839 versus 0.2381 arcsec. The wave joint score is only slightly lower (0.03828 versus 0.04226), and a post-test influence diagnostic reverses that ranking when J0935-0003 is omitted; all 12 remain in the primary score. The protocol was committed before predictions. This is not robust predictive superiority or a complete photon theory. Test residuals are now exposed; future model changes require a fresh independent evaluation.

- 02:47 · `d0c7d4d2` · from `README.md`

**Shared orbital sensitivity:** [Tangential/isotropic/radial comparisons](research_work/results/stellar-orbit-degeneracy/report.md) retain all mass calibrations and apply beta=-0.3, 0 and +0.3 to both hypotheses. Wave joint errors remain lower at every tested beta; reused-validation dispersion RMS spans 34.8–42.0 km/s for waves versus 40.0–50.5 for the stellar-mass benchmark. Lensing is unchanged. Independent Jeans integrations agree within 7.4e-6 in variance. This is a limited sensitivity result, not a measured anisotropy, a full mass/orbit posterior or fresh validation; the original calibration and unopened final test remain intact.

- 02:41 · `5698148e` · from `README.md`

**Stellar-mass degeneracy:** [One shared mass-correction benchmark](research_work/results/stellar-mass-degeneracy/report.md) fits lambda=1.645 on training and freezes it for reused validation. Its validation errors are 45.85 km/s and 0.3541 arcsec, versus 38.22 and 0.2699 for frozen waves. Much of the gain over unadjusted ordinary matter can therefore arise from a coherent stellar-mass shift; the wave advantage is smaller but remains in both validation RMS measures. The shift is not independently established as population-allowed, and wave/grid versus mass-benchmark flexibility differs. No final-test scores were opened.

- 02:38 · `4d4b128d` · from `README.md`

**Frozen-wave validation:** [Seven-system reused validation comparison](research_work/results/wave-lens-validation/report.md) keeps m=1e-24 eV/c^2, source fraction 0.6 and photometric masses fixed. Dispersion RMS improves from 92.12 to 38.22 km/s and lens-angle RMS from 0.6437 to 0.2699 arcsec versus the paired ordinary-matter benchmark; all seven absolute errors improve for both observables. The wave dispersion median is still 11% low and J0946+1006 remains a large mismatch. This is conditional reused validation, not pristine blinding or a complete photon theory. Calibration hashes are unchanged; final test predictions remain unopened.

- 02:34 · `8b5f6d20` · from `README.md`

**Shared wave calibration:** [Complete nine-setting training grid](research_work/results/wave-lens-training/grid-report.md) selects m=1e-24 eV/c^2 and source/stellar mass ratio 0.6 using a declared joint log-error score. The same setting minimizes both RMS errors: 31.99 km/s for stellar dispersion and 0.1440 arcsec for lens angle, compared with 62.77 and 0.4628 for the fixed-photometric-mass ordinary benchmark. Median ratios are 0.960 and 1.044. This is training improvement, not holdout validation or photon origin; all nine settings and residuals are retained. The shared calibration is frozen for a subsequent conditional validation comparison, with test predictions still unopened.

- 02:27 · `463030f1` · from `README.md`

**Supported-wave observational pilot:** [Fixed-photometric-mass comparison](research_work/results/wave-lens-training/report.md) tests 32 training lenses with shared m=1e-24 eV/c^2 and source/stellar mass ratio 1. Compared with the same Salpeter-based ordinary-matter masses, dispersion RMS changes from 62.77 to 52.50 km/s and Einstein-angle RMS from 0.4628 to 0.4400 arcsec, but wave predictions are median 11.1% and 28.7% high respectively. This is not a satisfactory joint match or holdout validation. Independent aperture, mass-scaling and projected-lensing checks show numerical errors are negligible relative to the discrepancies. Source normalization, population uncertainty and capture physics remain unresolved.

- 02:22 · `8880115b` · from `README.md`

**Supported wave profiles:** [Self-consistent stationary solver](research_work/results/self-consistent-wave/report.md) constructs nine nodeless Schrodinger-Poisson configurations in a fixed Hernquist ordinary galaxy, including source self-gravity. Enlarged-domain/inner-radius/tolerance checks give enclosed-mass agreement within 1.45e-7 of source mass; virial residuals are below 4e-9. Motion and conditional weak lensing follow from the same density. These are numerical equilibria, not galaxy fits, capture formation or collective-stability results. Shared physical field parameters and new observational likelihoods remain required.

- 02:18 · `02d6eff9` · from `README.md`

**Oscillating-field profile test:** [Single-state wave equilibrium](research_work/results/oscillating-field-profile/report.md) requires a particle mass varying by factors 75–90 between 0.1 and 4 effective radii to support the existing empirical source profiles. This rules out their exact completion as a free, stationary, no-flow single scalar state, even away from the cutoff. It does not exclude self-consistent wave profiles, self-interactions or mixed states; each needs a specified action and new predictions. No observational or holdout scores were changed.

- 02:16 · `ec2e8796` · from `README.md`

**Stored-field stress:** [Canonical-field comparison](research_work/results/stored-field-stress/report.md) shows that a literal frozen, minimally coupled canonical scalar cannot supply the positive, nearly pressureless extra source assumed by the lens pilot. Its tangential tension equals minus its energy density; for nonnegative potential energy its local temporal-potential source has the wrong sign. This is a restricted candidate failure, not a rejection of modified gravity or all stored deformation. A local oscillating-field example has pressureless average stress, but no supported galaxy or capture mechanism is established. All previous observational scores remain unchanged.

- 02:13 · `1c4a192d` · from `README.md`

**Companion-source support:** [Boundary and orbit-support audit](research_work/results/deposit-boundary-admissibility/report.md) finds that the sharp source boundary in all 99 training lens configurations cannot describe an isotropic collisionless particle reservoir. A necessary positivity condition fails in an outer region containing median 24–28% of the equivalent source mass, depending on cutoff. A nonnegative circular-shell construction preserves the profile but requires an unexplained angular-momentum supply; it is not a capture or collective-stability proof. This restriction does not apply automatically to stored spacetime deformation, which needs its own stress equation. Stellar-tracer admissibility and all previous lens residuals remain unchanged; no held-out scores were opened.

- 02:06 · `47ea8c02` · from `README.md`

**Stellar orbit admissibility:** [Known DF inversion applied to the pilot](research_work/results/lens-orbit-admissibility/report.md) establishes nonnegative stellar tracer distribution functions for beta=-0.3, 0 and +0.3 in the stipulated Hernquist-plus-truncated-extra-force family. Analytic curvature bounds cover all radii and positive amplitudes; numerical transforms and 54 field configurations support the implementation. This removes a tracer-consistency gap, not stability, observed anisotropy, deposit support or the missing companion mechanism. No new held-out scores were opened.

- 02:02 · `052c28fe` · from `README.md`

**Orbital sensitivity after validation:** [Training-only constant-beta test](research_work/results/lens-training-pilot/anisotropy-report.md) shows radial preference beta=0.3 reduces companion median angle overprediction from 9.7% to 4.2%, with RMS 0.1942 arcsec versus 0.1904 for the same-beta ordinary-matter benchmark. Tangential preference worsens the companion result. No beta is fitted per lens or promoted as measured; the old negative isotropic validation remains. Test predictions stay unopened. Independent orbital/population constraints and physical distribution-function checks are needed before a revised predictive claim.

- 01:58 · `d763b69e` · from `README.md`

**Reserved lens validation evaluated:** [Frozen seven-system comparison](research_work/results/lens-training-pilot/validation-report.md) gives Einstein-angle RMS 0.2505 arcsec for the empirical companion pilot versus 0.1712 for the ordinary-matter benchmark. All seven companion primary residuals are positive and larger in absolute value than the benchmark. Shared parameters, primary case and eligibility rules were fixed before execution; all sensitivity cases are retained. Validation is now exposed, test predictions remain unopened. This is a negative result for the restricted pilot, not a completed physical-theory test or proof of photometrically allowed ordinary-matter masses.

- 01:55 · `63a86c8d` · from `README.md`

**Stellar-population cross-check:** [Multiband photometry and updated-size rerun](research_work/results/lens-photometric-audit/report.md) cover 33 training lenses, 32 with published population masses. Under fixed-population normalization and event stretching, the updated companion pilot still requires median ordinary-matter masses 2.37 times Chabrier or 1.35 times Salpeter estimates. Updated lens-angle RMS is 0.1840 arcsec for the ordinary-matter benchmark versus 0.2083 for companions. No full population refit or holdout score is claimed. Propagation-dependent luminosity, population priors, realistic profiles and photometric errors remain necessary before model selection.

- 01:50 · `943b1eec` · from `README.md`

**First observed lens training comparison:** [33-system pilot](research_work/results/lens-training-pilot/report.md) infers mass from aperture dispersion and then predicts Einstein angle. The spherical isotropic Hernquist ordinary-matter benchmark has RMS 0.206 arcsec versus 0.243 with the frozen empirical extra-gravity relation for the declared representative geometry/seeing. Companion predictions are median 9.7% high. These are conditional training diagnostics with masses not independently constrained by starlight; no held-out score or theory confirmation is claimed. Next require photometric mass constraints, orbital/seeing treatment and a shared companion metric/source before selecting a final model.

- 01:45 · `20ed4bdf` · from `README.md`

**Measured lensing inputs:** [SLACS acquisition](research_work/results/lensing-data-readiness/report.md) archives 131 candidates and 63 modeled lenses, with 58 measured stellar dispersions. Deterministic system roles are 41 training, 9 validation and 13 test before lensing scoring. Observed redshifts and spectroscopy are separated from image-model summaries and cosmology-corrected luminosities. Conditional static geometry is calculated and labeled; no observational lensing fit or score is claimed. Next specify the common metric/source, ordinary-matter photometry and aperture-dispersion/imaging likelihood.

- 01:42 · `c073479d` · from `README.md`

**Lensing boundary audit:** [1,788 conditional predictions](research_work/results/lensing-boundary-audit/report.md) across 149 archived galaxy templates show that changing the outer deposited-source boundary changes extra light bending while preserving inner rotation. Median change is about 23% at half the last measured radius for the declared boundary sweep. Equal metric potentials are an additional assumption, not a companion derivation; no observed lensing data were scored. Three spherical interpolants fail positive-source eligibility and remain flagged. The capture boundary and metric response must be specified before a joint lensing claim.

- 01:39 · `33e818d8` · from `README.md`

**Deposit-source consistency:** [Poisson reconstruction](research_work/results/gravity-geometry-audit/source-density-report.md) finds positive equivalent companion densities at 63 field probes for both conservative geometries, with finite-step, bar-order and declared-source cross-checks. A small negative reconstructed ordinary-matter density is recorded as numerical error, not new physics. Indefinite power-law extrapolation requires unbounded equivalent mass; a finite spherical source continuation exists, but its transition radius and capture cause remain unpredicted. No energy normalization or held-out observational test is claimed.

- 01:33 · `3c3a7dd5` · from `README.md`

**Three-dimensional gravity ambiguity:** [Executed comparison](research_work/results/gravity-geometry-audit/report.md) freezes the existing empirical coefficients and compares two conservative well shapes. Both give 8.20 km/s RMS error on 38 exposed Milky Way rotation bins; their vertical-force RMS errors differ, 15.76 versus 36.13 surface-equivalent units on 43 model-dependent estimates. Directly scaling the local force fails closed-path work consistency. Numerical refinement does not erase this distinction. No new coefficients, raw-star repairs or held-out likelihood scores were introduced. A shared capture/response law must determine geometry; rotation agreement alone cannot do so.

- 01:29 · `e8b3f2f4` · from `README.md`

**Association audit:** [Epoch-aware positional check](research_work/results/stellar-orbit-support/association-report.md) finds 2.2–2.9 arcsec residuals for four extreme-speed training examples. One alternate Gaia source agrees within 0.045 arcsec; no replacement is adopted. Resolve source associations before combining parallax likelihoods or interpreting provisional high-speed launches as gravity evidence. Original parent and holdout roles remain unchanged.

- 01:22 · `fa2dde95` · from `README.md`

**Latest real-star test and correction:** [Training orbit support and input provenance](research_work/results/stellar-orbit-support/report.md) analyzes 27,884 training stars and identifies an omitted StarHorse input flag. Some distances did not use Gaia parallax; their treatment cannot follow the earlier blanket reuse assumption. The recovered metadata and newly downloaded fidelity scores preserve all original stars and holdout roles while flagging disagreements and different processing needs. This must be resolved before interpreting stellar fits as evidence for extra gravity.

- 00:53 · `9cbfb03c` · from `README.md`

**Latest Milky Way forward calculation:** [Rotating-bar orbits](research_work/results/rotating-bar-orbits/report.md) compares 30 synthetic trajectories under ordinary matter, a separate halo, and three deposit geometries. Integration, coordinate and sampled-field checks pass; higher-order bar/disk calculations quantify remaining numerical sensitivity. Added well depth changes complete orbital paths rather than imposing one speedup. This is not a stellar population fit. The [likelihood contract](research_work/results/rotating-bar-orbits/likelihood-contract.md) specifies the next training-data comparison and preserves stellar holdouts.

- 00:22 · `ab1b8097` · from `README.md`

**Latest joint clock/energy test:** [Matter and the autonomous receiver](research_work/results/matter-clock-closure/report.md). Including universal matter and clock response in the homogeneous receiver conserves energy but cancels the locally measured redshift and event stretch. Weakening the clock response changes local light speed under the tested fixed-ruler prescription. This is a specific candidate failure, not an exclusion of all localized time fields or photon conversion; it makes explicit what the next source–propagation–detector model must repair.

- 00:11 · `e1a25cb8` · from `README.md`

**Latest core-mechanism test:** [Energy conservation and the receiving state](research_work/results/coherent-receiver-audit/report.md). A finite energy-preserving quantum interaction can lower photon energies and retain high two-frequency beat visibility when its receiver is coherent and nonstationary. A stationary or dephased receiver does not retain that pattern in the tested construction. This identifies a physical state/resource that the propagation model must specify; it is not yet a supernova timing solution, a derived graviton interaction or an astronomical fit.

- 00:01 · `e6f44629` · from `README.md`

**Current theory-level priority:** [Broad review and executed coupling test](research_work/results/theory-priority-review/report.md). The gravity amplitude can exactly compensate changes in the photon-loss scale across all 3,150 retained rotation rows. Thus the useful rotation fit does not independently establish photon conversion as its cause. A shared timing/energy/receiving-sector law and a constrained gravitational response now take priority over further numerical tuning. A [provisional 3D field foundation](research_work/results/bar-field-foundation/report.md) is saved for the later orbital test; no new stellar holdout scores have been opened.

## 2026-09-09

- 23:49 · `e3009415` · from `README.md`

**Spatial holdouts and measurement uncertainty:** [Stellar holdout audit](research_work/results/stellar-holdout-audit/report.md) allocates 77,927 training / 24,755 validation / 26,090 test stars, propagates Gaia/StarHorse uncertainties in four explicit sensitivity modes, and evaluates the frozen photon-loss formula. About 5% of target assignments change under distance draws. The apparent bulge height trend changes sign after coarse population matching; no gravity inference follows from raw averages. Shared-field and prior-summary exposure limits are recorded. The bar-aware orbital likelihood remains unfinished.

- 23:40 · `e817974e` · from `README.md`

**Stellar data acquired:** [Bulge/disk data-readiness report](research_work/results/stellar-data-readiness/report.md). Recovered local BRAVA/Gaia observations; downloaded official APOGEE DR17 and StarHorse v2, matched 554,724 distinct sources, prepared 140,407 giants and obtained Gaia DR3 covariance/quality fields for all of them. Additional Gaia cuts retain 128,772 candidates. Includes reproducible preparation, hashes, coverage and an audit of AGAMA's bar/disk/halo components. This prepares the bulge test; it does not yet fit a dynamical model.

- 23:25 · `471d4245` · from `README.md`

**Latest joint analysis and user clarification:** [Galaxy/redshift report](research_work/results/joint-galaxy-audit/report.md), [observed/predicted explorer](research_work/results/joint-galaxy-audit/comparison.html), and [local well deepening beneath/above the bulge](research_work/results/joint-galaxy-audit/bulge-local-well.md). The suite recalculates 164 redshift groups and 149 rotation-curve galaxies, compares gas-inclusive Milky Way models, and retains both successful fits and remaining conflicts. A deposit may reshape a scalar well without an arrival-directed push; spherical cold-source fits are comparison models, not the selected mechanism. Fifty-seven bulge field summaries are imported but not yet dynamically predicted. No complete all-observation theory is claimed.

- 23:25 · `471d4245` · from `research_plan/START-HERE.md`

**Latest:** [Joint galaxy/redshift analysis](research_work/results/joint-galaxy-audit/report.md) and [local well response/bulge comparison](research_work/results/joint-galaxy-audit/bulge-local-well.md). A companion deposit may deepen or reshape a local well; incoming direction is not an imposed force direction. Keep ordinary cold-source fits as conditional benchmarks. New bulge field moments have no derived dynamical prediction yet. The actual/predicted explorer marks these missing predictions explicitly.

- 22:59 · `1223295d` · from `README.md`

**Milky Way capture update:** [Radial and vertical diagnostic](research_work/results/milky-way-capture/report.md), [data and calculations](research_work/results/milky-way-capture/results.json), and [figure](research_work/results/milky-way-capture/summary.png). Imported 38 rotation bins and 43 provisional vertical-force estimates; recalculated the user's six speed-excess rows and tested 48 all-direction capture geometries. Rising percentage speed excess does not imply increasing extra acceleration. Outer capture can also load the disk faces. A common capture/storage/gravity model and an assumption-audited vertical-tracer fit remain required; no joint fit or energy-supply success is claimed.

- 22:59 · `1223295d` · from `research_plan/START-HERE.md`

**Milky Way priority:** [Executed radial/vertical capture diagnostic and next goals](research_work/results/milky-way-capture/report.md). We now have 38 rotation bins, 43 provisional vertical-force estimates, and an all-direction ellipsoid capture control. Before a joint fit, recover the stellar/gas baseline and connect incident intensity, one capture rule, supported storage, and both force components. Preserve the published vertical inference's halo-model dependence pending raw-tracer reanalysis. No cosmological expansion or dark matter source is adopted by this test.

- 22:50 · `d02ceef6` · from `README.md`

**Current formulas and audit:** [Electromagnetic transfer specification](research_plan/electromagnetic-transfer-specification.md), [executed results](research_work/results/electromagnetic-audit/report.md), [comparison figure](research_work/results/electromagnetic-audit/summary.png), and [coverage/missing inputs](research_work/results/electromagnetic-audit/coverage.md). The suite recalculates 164 group predictions, six maser galaxies, 35 spectral-aging summaries, two independent radio centroid groups, and 43 FIRAS channels. Seventy synthetic band/path predictions and 60 prescribed-field cases are separately labeled. A new inverse-affine time candidate links wavelength and event stretch exactly, but source, companion loss and background consistency remain unresolved. Comparable performance to expansion is sufficient; no global parity is claimed.

- 22:50 · `d02ceef6` · from `README.md`

**Baseline comparison:** [Direct Photon Conversion and Bound Companion Fields](research_plan/direct-conversion-bound-companions.md) and its [earlier benchmark](research_work/results/direct-conversion/report.md) remain preserved. The user has reopened altered time as an alternative cause; neither a static energy-loss rule nor the prescribed time field has passed the full observational program.

- 22:50 · `d02ceef6` · from `research_plan/START-HERE.md`

**Latest scope:** evaluate both energy conversion and altered time against the electromagnetic observations. Read the [current formulas](research_plan/electromagnetic-transfer-specification.md), [executed electromagnetic audit](research_work/results/electromagnetic-audit/report.md), and [coverage register](research_work/results/electromagnetic-audit/coverage.md). Direct conversion is the baseline, not an exclusive branch choice. The new inverse-affine time profile repairs wavelength/event kinematics but is not sourced or compatible with no-loss same-metric companions. No complete all-data solution or fresh validation is claimed. The historical direct-conversion decision below is superseded on the question of which causes may be investigated.

- 18:55 · `3b831632` · from `README.md`

**Earlier conceptual specification:** [Environmental Time and Bound Companion Fields](research_plan/environmental-time-bound-companions.md) preserves the environmental-time branch and the development of the bound-level picture. Its time and screening prescriptions are not silently included in the active direct-conversion branch. The existing paper PDF remains a separately versioned artifact.

- 18:55 · `3b831632` · from `research_plan/START-HERE.md`

The latest integrated step is the [direct-conversion redshift benchmark](research_work/results/direct-conversion/report.md): z=0.0076605 at 100 million light-years with the reused rate, and very similar recovered-data redshift errors to explicit expansion comparisons. This does not yet establish complete parity, a microscopic interaction, transient time stretching or a formed halo. Earlier [pair production and loss balance](research_work/results/pair-production-balance/report.md) calculations remain relevant conditional constraints, not demonstrations of unlimited stable growth. Retain all 20 tasks/32 requirements.

- 13:59 · `455a6190` · from `README.md`

**New working paper:** [Environmental Time Stretching and Companion-Energy Deposition](papers/cumulative-time-companions/manuscript.md), with a [PDF](papers/cumulative-time-companions/manuscript.pdf) and [evidence map](papers/cumulative-time-companions/evidence-map.md). This current synthesis develops the cumulative-time premise, conditional derivations and remaining research objectives; the historical version 9 below is preserved separately.
