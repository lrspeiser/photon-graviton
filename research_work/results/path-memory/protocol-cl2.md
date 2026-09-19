# CL-2, stage 1: one written-field response tested jointly on galaxies, galaxy lenses and clusters

Declared before execution, 19 September 2026. Baseline: this branch at 40144ef, which is CL-1 (ece306a) merged with `main` at 271def5, so the owner's stage-8 corrections and the prior-art audit are both in the tree. CL-1's equations, archive and readings stay frozen; this stage adds files and edits none of CL-1's.

**The objective, stated by the owner:** one explanation of gravitational behaviour that works for galaxies and clusters without inserting an unexplained dark-matter distribution. Reproducing a galaxy law while leaving clusters unresolved does not meet it. **The deliverable of this stage** is the reviewer's first milestone: given independently specified ordinary matter, what field does the written-track mechanism predict for a galaxy, a lens galaxy and a cluster, using one law and common parameters — and does any member of the family serve all three?

## 0. The reviewer's plan: what is adopted, what is corrected, what is deferred

The plan (five work packages, an implementation order, and the rule that a static cross-scale screening precedes any formation campaign) is adopted as the structure of CL-2. Six points are corrected or qualified before anything runs, each with its reason.

1. **The point-source spectrum formula is reproduced**, and one more analytic member is recorded. For amplitude dΛ = A·d ln w the Gaussian forces sum to g = (AM/r)[e^{−r²/2w_max²} − e^{−r²/2w_min²}], numerically 1 ± 10⁻⁷. And with weight dΛ = A·w⁻¹·d ln w the potential sums to A·M·√(π/2)/d exactly: **the Newtonian Green's function is a member of the family.** So a rescaled G is inside the space of universal spectra, which is what a cluster's flat boost wants, and a linear law with v² ∝ M is all the far field can give a galaxy, which is what PM-1's B-linear rule already showed fails the mass–speed relation. Both are known before the run; they are what make the diagnostic decisive rather than open-ended.
2. **Widths narrower than the measured density shells are excluded from the cluster block.** The pressure operators for w < 10 kpc are not grid-converged (they change by 11–17% when the radial grid is doubled, because the electron density is measured in shells of that size and a 0.1-kpc footprint reads the interpolation's kinks); above 10 kpc they change by 3.6×10⁻⁵. Those widths stay in the family and are constrained by the galaxies and lenses only. This is a measured numerical limit, not a physical switch.
3. **"No far field" is measured, not asserted.** The reviewer is right that CL-1's gas extended to 20 Mpc and the deficit at R₂₀₀ lies inside it; the exact statement is that a local response to a declining density falls faster than the enclosed-mass force. The shell decomposition below says which shells the force at each radius comes from, with the gas continued or truncated.
4. **The cluster comparison is the forward pressure**, as the plan asks: with the measured electron density, a model gravity predicts the pressure profile up to one boundary value per cluster, and that is compared with the measured X-ray and SZ pressure points. The X-COP release makes this possible; CL-1's five NFW numbers per cluster stay archived as the regression case only.
5. **Coma stays shape-only and outside the joint fit.** Calibrating its lensing amplitude needs the source population under the adopted geometry, which this repository does not hold; CL-1's Coma record stands.
6. **No validation or test SLACS lens has prepared inputs**, so the reserved-lens forecast the earlier summary suggested cannot be run in this stage. The reserved information here is: the SPARC validation and test splits (the spectrum is fitted on the training split only), lensing-only fits predicting stellar motions and motion-only fits predicting Einstein radii, and the seventh KCWI lens, J1538+5817, which CL-1 never used and which is scored, never fitted.

CL-2E (formation, energy and propagation at the inferred strengths) is deferred until a static response survives this screening, as the plan says.

## 1. Inputs and sources (CL-2A)

**Clusters: the X-COP data release**, acquired by `cl2_xcop_acquire.py` (archive SHA-256 pinned, 315 MB cached and ignored, only the profile tables read) into `cl2-inputs-xcop-profiles.json`, committed. Per cluster: the electron-density shells with the release's bounds (57–69 shells to 1.6–2.7 Mpc), the X-ray pressure points (n_e kT, 8–19 per cluster, to ≈0.9 R₅₀₀) and the SZ pressure points (Planck, 9–11 per cluster, to 1.8–2.6 R₅₀₀), the release's gas-mass profile, its hydrostatic mass reconstructions (forward, NFW, Einasto) and, for seven clusters, the cumulative stellar-mass profiles of Ghizzardi et al. Every quantity is a release product (deprojection, spectral fits, the release's H₀ = 70 distances adopted as scenario inputs); none is a raw observation, and the hydrostatic masses are used only as a consistency floor.

- **The forward pressure.** For electron pressure P_e and gas density ρ_gas = μ_e m_p n_e, hydrostatic balance reads dP_e/dr = −μ m_p n_e g with the release's μ = 0.6 and μ_e = 1.14, so P_e(r) = P_out + ∫_r^{r_out} μ m_p n_e(s) g(s) ds. r_out is the smaller of the last density shell and the last SZ point; P_out is one nonnegative nuisance per cluster; the points used are all X-ray and SZ pressure points inside r_out (246 in all), with the release's symmetric errors, diagonal. The predicted pressure is **linear in the force**, so it is linear in the amplitudes of any spectrum and the boundary values, and the joint solve stays a nonnegative least-squares problem. Thermal pressure is taken equal to total pressure in the primary run; the declared sensitivity applies Eckert et al.'s non-thermal fractions (Table 2, α at R₅₀₀ and R₂₀₀, interpolated linearly in r/R₅₀₀ from zero at the centre, capped at R₂₀₀'s value) to the model's total pressure. μ ∈ {0.59, 0.61} is the second sensitivity.
- **The baryons.** Gas mass from the release's own profile; stars from the release's cumulative profiles where they exist (7 clusters) and, for the other five, the mean stellar-to-gas ratio profile of those seven as a function of r/R₅₀₀ (declared, reported as a limitation; a central galaxy's cusp is present only in the seven).
- **The consistency floor.** The release's own NFW mass profile run through the same pressure operator: it is the best any model can be asked to do on these points with this pipeline (prototype: χ²/N = 5.82 over 246 points, per cluster 2.4 to 12.4, the merging A2319 worst). A model is "described" in the cluster block if its χ²/N is within a factor 2 of that floor.
- **Diagnostics (E4):** the contribution of each source shell to the model force at each pressure radius for the best spectrum; the same with the gas truncated at 1, 1.5, 2 and 3 R₅₀₀ instead of continued; per-cluster χ² and residual patterns; the CL-1 five-point comparison recomputed with the new baryons as the regression case.

**Galaxies: SPARC and the Milky Way.** The 149 galaxies with the frozen 89/29/31 split; sources from RPG-1's `baryons.sparc_components` (stellar and gas surface densities, spherical bulge), the written field of a razor-thin disk in its plane through RUT-1's ring kernel and of the bulge through the shell kernel; the Newtonian baseline is SPARC's tabulated ordinary-matter contribution at Υ = 0.5/0.7 as PM-1 used it, with PM-2A's caveat recorded (the reconstructed source and the tabulated contributions differ at some radii; recomputing the disk's Newtonian field from the source is PM-2A stage C and is not done here). Errors: σ_g = 2 v e_v/R, floored at 0.1% of g. Scored by the repository's equal-galaxy RMSE and by χ²; the spectrum is fitted on the training split only; validation and test are reported unchanged. The Milky Way's 38 Eilers bins (baselines I and II) are transfer targets, never fitted.

**Lenses: six SLACS systems plus one held out.** J0037-0942, J1112+0826, J1204+0358, J1402+6321, J1621+3931, J1630+4520 through CR-2's measurement interface — the deprojected Sérsic light, the KCWI apertures, PSF and covariance, `AnnularModel`'s second moments — anchored to CR-2's archive; J1538+5817 (light profile, masses and KCWI profile all prepared, never used in CL-1) is scored under every fitted law and enters no fit. **Geometry registry** `cl2-geometry.json`: one record per lens with three labelled scenarios — G1 flat FLRW (CR-2's, the primary here because the kinematics anchors live in it), PF-1 static Euclidean (CL-1's) and G2 co-scaling coasting — each fixing D_l, D_ls/D_s, the aperture and PSF scales and the stellar-mass rescaling (D_L/D_L,FLRW)². The Einstein constraint enters as one row per lens with a declared 3% uncertainty on the required bend; the second moments enter with the release covariance, linearised in V²_rms about the measured values. Stellar masses: population values (Chabrier primary, Salpeter reported); the lens-determined mass M1 is CR-2's benchmark and is reported beside. Anisotropy β constant per lens within CR-2's bounds [−2, 0.45], refitted for every trial law. Light rule L1 throughout; L0 is the archived Newtonian bracket.

## 2. The response family (CL-2B)

    C(x) = Σ_j Λ_j ∫ ρ_b(y) exp[−|x−y|²/(2w_j²)] d³y,    Λ_j ≥ 0,   w_j on 26 logarithmic points from 0.1 kpc to 10 Mpc (5 per decade),

the same widths and amplitudes for every galaxy, lens and cluster, with the unnormalised Gaussian convention of RUT-1 and CL-1 kept so the amplitudes keep their meaning (w/ℓ = Λw/G). `cl2_response.py` sums CL-1's verified per-width operators; nothing about a single width is recomputed. The objective is the sum over blocks of χ²/N_block, so no block's size buys it weight; each block is also solved alone. There is no smoothing penalty in the primary solve — the amplitudes span ten decades and a linear penalty has no natural scale — and the sensitivity is the width grid itself: 4 and 6 points per decade, and endpoints 0.03 kpc–30 Mpc.

**The decision rule, fixed now.** A block is *described* by a solve if: clusters, χ²/N ≤ 2 × the release-NFW floor; galaxies, training RMSE ≤ 1.1 × simple MOND's 19.89 km/s (that is ≤ 21.9) and the mass–speed slope within 0.05 of the observed 0.289; lenses, every Einstein bend within 3% and the kinematics χ² ≤ 1.5 × CR-2's free-NFW benchmark summed over the six. Then: **(a)** all three described by the joint solve — a shared spectrum works provisionally, reduce it and predict; **(b)** each described alone but not jointly — a universal linear convolution is missing source or environmental physics; **(c)** some block not described even alone — within the tested sources, bandwidth and light rule, changing widths is not enough, and the next work is nonlinear response, transport or the light coupling, not more widths. The prototype (quoted so it is not mistaken for a discovery of the run): clusters alone χ²/N ≈ 9.1 against the floor 5.8 — described; galaxies alone 30.95 km/s with slope 0.388 — not described; the two blocks use disjoint widths and the joint solve barely moves either. Outcome (c), driven by the galaxies, is the expectation; the run will say.

## 3. The lens tests (CL-2C)

With the six lenses' Einstein rows and second moments, four comparisons, all under L1 with population stellar masses:

| fit | what is predicted, not used in the fit |
|---|---|
| lensing only (Einstein rows) | the stellar velocity profiles |
| stellar motions only | the Einstein radii |
| joint, common amplitudes | J1538+5817's bend and motions; the Einstein residuals and motion χ² per lens |
| CR-2's benchmarks, reproduced | stars-only (lens-determined mass) and free NFW, for scale |

The regularity CL-1 found — six lenses agreeing on one amplitude at w ≈ 4.6 kpc to a factor 1.39 — is tested here by the stellar motions: if the amplitude the Einstein radii want at that width is not the amplitude the second moments want, the regularity is a coincidence of one observable.

## 4. Transfer (CL-2D)

Three calibrations, each carried without retuning: on the training galaxies, predicting lenses and clusters; on the clusters, predicting galaxies and lenses; on galaxies and clusters together, predicting the lenses, J1538+5817, the SPARC validation and test splits and the Milky Way. Residuals are reported per observable and per system, never as one aggregate. These are exploratory transfer tests on exposed data; the SPARC splits and J1538+5817 are the only role-respecting reservations.

## 5. Numerical verification: the gates, each with the control it must reject

| gate | requirement | measured before declaration | control |
|---|---|---|---|
| G1 analytic members | log-flat spectrum against (AM/r)[e^{−r²/2w_max²} − e^{−r²/2w_min²}] and 1/w spectrum against √(π/2)·M/r², point source, r ≤ 3 kpc inside the band: < 10⁻⁶ | 1.0×10⁻⁷, 1×10⁻⁸ | the 1/w member with √π in place of √(π/2) must fail |
| G2 single-Gaussian limit | the spectrum machinery with one width at CL-1's archived best reproduces CL-1's E1 χ² and Λ to 10⁻⁹ | same operators | a second width with nonzero amplitude must change it |
| G3 mixture lensing routes | three-dimensional deflection of a three-width mixture on the Plummer sphere against the two-dimensional ring-kernel route: < 10⁻⁶ | per width 3.6×10⁻⁷ (CL-1) | the mixture with the shell kernel in the 2D route must fail |
| G4 grid convergence | cluster pressure operators (w ≥ 10 kpc) under a doubled radial grid: < 10⁻⁴; galaxy operators under a doubled source grid: < 10⁻³ (all 89 training galaxies); the excluded widths' non-convergence reported | 3.6×10⁻⁵; 1.1×10⁻⁴ on three galaxies | — |
| G5 galaxy anchor | baryons-only, PM-1's law at its fitted a\* and simple MOND on the training split reproduce the archived 52.57, 20.21 and 19.89 km/s to 0.02 | 52.56, 20.21, 19.89 | Υ_disk = 1 must move the baryon score by more than 5 km/s |
| G6 lens anchor | CR-2's archived stars-only χ² and β for the six lenses in G1 reproduced to 10⁻⁶ | 10⁻¹² on two | a PSF doubled must move χ² by more than 1% |
| G7 X-COP ingestion | the extract's hash pinned; the gas mass integrated from the density shells against the release's own gas mass at R₅₀₀ within 5% for all 12; the release's NFW through the pressure operator χ²/N < 10 | 0.959–0.990; 5.82 | μ = 1.2 must raise the consistency χ²/N above 10 × its value at 0.6 (measured 595 against 7.9 on three clusters) |
| G8 the solve | Karush–Kuhn–Tucker: gradient on active amplitudes below 10⁻¹⁰ of the objective's scale, nonnegative on inactive ones | 10⁻¹⁶ | — |

`cl2_checks.py` is the suite job; it reruns G1–G8 (G4 at reduced size) and anchors the archived readings. A failed gate blocks every reading that depends on it; the three statuses are reported separately.

## 6. What stage 1 does not do

No formation or energy accounting at the inferred strengths (CL-2E); no Coma amplitude; no reserved-lens forecast beyond J1538+5817; no Milky Way vertical force; no recomputation of the disks' Newtonian field from the reconstructed source; no derivation of the light rule. Every observational reading is exploratory on exposed data, and "described" is a declared threshold, not a claim of correctness.

## Files

`cl2_xcop_acquire.py` and `cl2-inputs-xcop-profiles.json` (the acquisition and the extract, committed with this declaration as inputs); `cl2_response.py` (the family's operators and the nonnegative solve); `cl2_sources.py` (cluster, galaxy and lens blocks with their anchors); `cl2.py` (gates, then E1–E4, writing `cl2-results.json`); `cl2_checks.py` (the suite job, registered before stage 8's); `cl2-geometry.json`; `report-cl2.md`.
