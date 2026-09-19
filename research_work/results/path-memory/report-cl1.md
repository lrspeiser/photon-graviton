# CL-1: the written field has no far field — it can carry a cluster's interior but not its outskirts, and no one footprint width serves both clusters and lens galaxies

Protocol [protocol-cl1.md](protocol-cl1.md), declared in ccde494 before any run. Library `steady_field.py`, inputs `cl1-inputs-xcop.json`, driver `cl1.py`, archive `cl1-results.json`, suite job `cl1_checks.py`. RUT-1's files are imported and unchanged; nothing of stages 0–8 is edited.

This is the first test of the RUT-1 response outside the planar orbital problem, on goal 8's terms: the steady state of the two-stage equation — writing proportional to mass, the Gaussian footprint, retention — applied to spherical sources, with two declared light rules. Two constants enter and only two: the amplitude Λ = α·τ_keep, in (km/s)² per M☉, and the width w. With ℓ ≡ G/Λ, the dimensionless writing strength w/ℓ puts RUT-1's supported orbits (w/ℓ = 0.50 at 10% support) on the same map as the clusters and the lenses.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | **passes.** `cl1_checks.py` reruns the nine gates with their controls (G9 at reduced size) and recomputes fourteen anchored numbers from scratch; all match the archive to 10⁻⁹. 41 s. |
| **numerical verification** | **passes: 9 of 9 gates, every negative control rejected.** One thing to record: the archive is from the *second* run. On the first run both of G4's declared criteria passed with the numbers below, but a third criterion I had coded into that gate without declaring it — Gauss's law between the gradient kernel and the Laplacian kernel at interior radii — read 1.8×10⁻⁴ to 1.2×10⁻³, because its quadrature ended at the grid point below the target radius instead of at the radius. The quadrature was corrected to end at the radius (it now reads 6.5×10⁻⁹ to 7.1×10⁻⁷), the first archive was discarded before being committed, and no declared tolerance changed. A gate failed on a first run, whatever the cause, and that is written down. |
| **scientific outcome** | **exploratory, on exposed data, as declared.** A domain map, not a fit: (i) one universal (Λ, w) cuts the X-COP cluster misfit to χ² 4,383 against 83,267 for the baryons and 24,469 for PM-1's law, but with a fixed pattern — right at 0.5–1 Mpc, half the required force at R₂₀₀ — that follows from a theorem verified here: the written field of a bounded source carries **no net equivalent mass**; (ii) under the rule that light sees the total potential, each of the six SLACS lenses can be brought to its observed Einstein radius, and at w ≈ 5 kpc the six need the same amplitude to within 39% (Chabrier) — but the cluster-fitted pair leaves every lens exactly where the Newtonian baryons leave it (0.30–0.65 of observed); (iii) the same pair steepens Coma's shear profile and worsens its shape (χ² 4.28 → 10.2); (iv) **no width brings the clusters and the six lenses within a factor of two of one amplitude**. The closest approach is a factor 2.6–2.8 at w ≈ 70 kpc, where neither system is described. |

In plain words: the written field is a smoothed copy of where the matter is. It pulls toward where matter is densest and fades a few footprint widths outside it. So it can add pull *inside* a cluster but cannot supply the pull the outskirts need, and a footprint the size of a cluster does nothing to a galaxy while one the size of a galaxy does nothing to a cluster.

## What was computed

For a spherical density the footprint's shell average has a closed form, k_w(r, r′) = exp[−(r − r′)²/2w²]·(1 − e^{−2u})/(2u) with u = rr′/w², so the steady field C = Λ∫k_w dM, its gradient and its Laplacian are one-dimensional integrals with analytic kernels. The extra inward acceleration is −ΛdC/dr; the equivalent Newtonian mass M_eff = r²g_mem/G; the equivalent density ρ_eff = −∇²C/4πG. Under the light rule **L1** (light sees Φ_N − C with the factor two, RPG-1's R4) the memory field lenses as ρ_eff, and because ∫∂²_zC dz = 0 its projected equivalent density is a two-dimensional Laplacian of the projected baryons convolved with the two-dimensional Gaussian — whose azimuthal average is exactly RUT-1's ring kernel. So the lensing of the written field is computed from the observed light profile with `rut1.py`'s own kernel and its analytic derivatives, and needs no deprojection. Under **L0** light sees the ordinary matter only, which is the archived Newtonian bracket.

## Numerical verification, as it fell

| gate | requirement | measured | the control it had to reject |
|---|---|---|---|
| G1 shell kernel | closed form vs direct angular quadrature, seven triples, < 10⁻¹² | **1.5×10⁻¹⁴** | RUT-1's planar ring kernel in its place: off by a factor 8 — rejected |
| G2 derivatives | analytic dC/dr and ∇²C vs fourth-order differences, Plummer sphere, < 10⁻⁸ (scale-relative, since ∇²C crosses zero) | **2.7×10⁻¹²**, **2.7×10⁻¹³** | the local-limit force −(2π)^{3/2}w³Λρ′: off by 12%, 253% and 8,990% at w = 0.5, 2, 8 kpc — rejected |
| G3 link to RUT-1 | ring kernel = `rut1.phi_ring`; `rut1.gate_kernel` and `gate_mature_ring` pass through import; the mature-ring control's w/ℓ = 0.501; second derivative < 10⁻⁸ | **identical (0.0)**; both pass; **w/ℓ = 0.5007**; **3.0×10⁻¹²** | the second derivative without its cross term: 5.0×10⁻⁴ — rejected |
| G4 no net equivalent mass | M_eff(r_t + 8w)/(ΛM/G) < 10⁻⁸; ∫ρ_eff 4πr²dr < 10⁻⁶ ΛM/G, at w = 0.5, 2, 8 kpc | **4.5×10⁻¹⁹, 1.1×10⁻¹⁷, 3.1×10⁻¹⁶**; **−6.5×10⁻¹³, −1.0×10⁻¹¹, −7.8×10⁻¹²**; (the undeclared interior Gauss check, above: 6.5×10⁻⁹ to 7.1×10⁻⁷ against 10⁻⁶, the thinnest margin in the stage) | the Newtonian equivalent mass of the same source at the same radius, 1.00 — rejected |
| G5 two lensing routes | three-dimensional deflection integral (shell kernel) vs two-dimensional ring-kernel projection, 12 (w, b) cases, < 10⁻⁶ and falling under grid doubling | **3.6×10⁻⁷**, **9.1×10⁻⁸** doubled, falls in all 12 | the two-dimensional route with the shell kernel: 4.4% — rejected |
| G6 lens anchor | Λ = 0 reproduces RPG-1's archived Newtonian θ_E, six lenses, two IMFs: < 10⁻¹⁰ (3D route), < 10⁻⁵ (2D route) | **5.1×10⁻¹⁴**; **1.8×10⁻⁶** | Λ = 10⁻⁶ at w = 3 kpc moves every θ_E by at least 25% — rejected |
| G7 X-COP ingestion | NFW(c₂₀₀, R₂₀₀, M₂₀₀) reproduces the four tabulated masses < 0.5% and both radii from the paper's cosmology < 0.1%; f_gas,500 agrees with the earlier transcription | **0.40%**; **0.038%**; exact for all 11 | A85 and A2255 concentrations swapped: 57% — rejected |
| G8 Coma ingestion | Σ(γ_t/σ)² = 23.33 to 0.1% | **23.336** | — |
| G9 convergence | cluster χ²(w) and Λ\*(w) under a doubled radial grid < 10⁻⁵; every lens's Λ(w) under a doubled projected grid < 10⁻⁴ | **5.9×10⁻⁷**; **1.0×10⁻⁷** | — |

The two lensing routes share no kernel — one is the shell average of the footprint, the other RUT-1's ring average — and their agreement to 3.6×10⁻⁷, improving fourfold under refinement, is the substantive check that the written field's lensing has been computed rather than assumed.

## E1 — twelve X-COP clusters, sixty published radii

Baryons: a β-model gas matched to the measured M_gas(<R₅₀₀) and M_gas(<R₂₀₀) with r_c = 0.15 R₅₀₀ (fitted β from 0.61 to 0.90), stars at 0.09 of the gas and distributed like it (baryon fractions at R₅₀₀ from 0.116 to 0.206). Required: the published backward-method NFW masses at 0.5, 1, 1.5 Mpc, R₅₀₀ and R₂₀₀ with their published errors — five values per cluster from one two-parameter fit, so **not independent**, and hydrostatic, so model-inferred. The clusters need 3.9 to 10.0 times their baryonic mass (median 6.5).

| law | free constants | χ² (60 points) | log-rms |
|---|---|---|---|
| Newtonian baryons | 0 | 83,267 | 0.817 dex |
| PM-1's law g_N + √(a\*g_N), a\* = 8.56×10⁻¹¹ m/s² frozen (= Completions I and II in spherical symmetry) | 0 | 24,469 | 0.212 dex |
| the same at PM-1's own fitted 6.54×10⁻¹¹ m/s² | 0 | 29,663 | 0.252 dex |
| **the written field, one (Λ, w) for all twelve** | 2 | **4,383** | **0.159 dex** |

The best pair is **w = 531 kpc, Λ = 7.77×10⁻⁸ (km/s)²/M☉, ℓ = 55 kpc, w/ℓ = 9.6**. The sensitivities move it little and are not selected by fit: r_c = 0.10 R₅₀₀ gives w = 531, w/ℓ = 9.2, χ² 4,228; r_c = 0.25 gives w = 447, w/ℓ = 10.8, χ² 3,882; M_star/M_gas = 0.07 and 0.12 give χ² 4,402 and 4,355 at the same w. χ²/point is 73 against published errors of 1–10%, so this is not a fit that describes the data; it is a factor 5.6 below the PM-1 law and 19 below the baryons, with a **fixed pattern**. The median of model/required across the twelve clusters is

| 0.5 Mpc | 1 Mpc | 1.5 Mpc | R₅₀₀ | R₂₀₀ |
|---|---|---|---|---|
| 0.90 | 1.24 | 0.84 | 1.05 | **0.56** (0.27 to 0.82) |

too much force at 1 Mpc and half the required force at R₂₀₀ in every cluster. That is the theorem of G4 at work: the written field's equivalent mass integrates to zero, so the force it adds must fall away outside the matter, and a footprint of half a megaparsec has run out by two. The χ²(w) curve says the same from the other side: for w ≤ 100 kpc the field is local, Λ\*(w) scales as w⁻³ (w/ℓ = 4,400 at 10 kpc) and χ² sits on a plateau at 24,300, fitting the tightly measured inner points and nothing else; the minimum at 450–530 kpc is sharp; beyond 1 Mpc the field is a blob wider than the cluster and χ² rises to 51,800 at 10 Mpc.

## E2 — the six SLACS lenses, under L1

θ_E rises monotonically with Λ at every width (checked at ½Λ and 2Λ), so each lens has one amplitude Λ_i(w) that gives its observed Einstein radius. Under L0 the archived Newtonian ratios stand: 0.30–0.42 of observed with Chabrier masses, 0.49–0.65 with Salpeter.

| lens | L0 ratio (Chab / Salp) | smallest w/ℓ (Chab / Salp) | at w | w/R_e | Λ at w = 4.64 kpc (Chab / Salp) |
|---|---|---|---|---|---|
| J0037-0942 | 0.328 / 0.534 | 4.49 / 1.88 | 4.64 kpc | 0.75 | 4.16 / 1.74 ×10⁻⁶ |
| J1112+0826 | 0.299 / 0.487 | 5.32 / 2.35 | 4.64 | 0.61 | 4.93 / 2.18 |
| J1204+0358 | 0.415 / 0.649 | 3.18 / 1.15 | 3.16 | 1.18 | 3.54 / 1.29 |
| J1402+6321 | 0.329 / 0.554 | 3.93 / 1.60 | 4.64 | 0.47 | 3.64 / 1.48 |
| J1621+3931 | 0.310 / 0.503 | 4.99 / 2.16 | 4.64 | 0.92 | 4.62 / 2.00 |
| J1630+4520 | 0.360 / 0.578 | 3.97 / 1.59 | 6.81 | 0.78 | 3.89 / 1.56 |

Two things are worth stating plainly. The amplitude each lens needs is smallest at a width comparable to its effective radius (0.5–1.2 R_e), and **at w = 4.64 kpc the six lenses need the same Λ to within a factor 1.39 (Chabrier) or 1.69 (Salpeter)** — the spread is 5 to 7 at w < 2 kpc and about 3 above 100 kpc. Six early-type lenses with Einstein radii from 1.29″ to 1.78″ agreeing on one amplitude at one width is a regularity, not a prediction: the data are exposed, the number of lenses is six, and the width was found by looking. It is recorded so that a later stage can test it on lenses this repository has not touched. And the cluster-fitted pair applied to the lenses gives θ_E ratios of 0.328–0.415 (Chabrier) — **identical to L0 to three decimals**: a 531-kpc footprint smooths a galaxy into a uniform potential offset and bends no light.

## E3 — Coma, both rules

Six figure-reconstructed Kubo bins at 1.7 to 13.8 Mpc (h = 0.7 adopted), CF-1's two baryon brackets, Σ_crit a nuisance. Under L0 the baryons' ΔΣ shape gives χ² 4.28, the same for both brackets (the stars follow the gas, so the brackets differ only in amplitude; the implied Σ_crit is 8.0×10⁸ and 1.5×10⁹ M☉/kpc²). Under L1 at the cluster-best pair the written field multiplies ΔΣ by **6.1 at the innermost bin, 1.7 at the second and 1.00 beyond**, and the shape χ² rises to **10.2**: the profile steepens where the data are flat. A free two-nuisance fit prefers a *negative* Λ (−2.8×10⁻⁹) for a χ² of 4.12 — the data would rather the written field removed shear inside 2 Mpc than added it. The reduced-shear correction is at most 0.018 plotted standard errors. The inverse fits' 3.855 (NFW) and 3.727 (Plummer) have more freedom and are not a benchmark.

## E4 — the compatibility map

No width, for either IMF, brings the clusters and all six lenses within a factor of two of one amplitude — the criterion declared before the run. The closest approach is at w = 68 kpc, a factor 2.8 (Chabrier) or 2.6 (Salpeter); but there the clusters' own best χ² is 24,356, the local-regime plateau, so the crossing describes nothing. At the cluster-best width the lenses need 385 to 1,180 times the cluster amplitude (Chabrier); at the lens-preferred widths of 3–7 kpc the clusters would need 2,400 to 13,300 times the lens amplitude, and the lens-preferred pairs applied to the clusters give χ² 83,218–83,258 against the Newtonian 83,267 — no effect at all. On the same map, RUT-1's supported orbits sit at w/ℓ = 0.50 at w/R = 0.1: the clusters want w/ℓ ≈ 10 at w ≈ 0.4 R₅₀₀, the lenses w/ℓ ≈ 3–5 (Chabrier) at w ≈ R_e. Each system prefers a width near its own size and a strength ten to twenty times RUT-1's; **a universal length does not serve both, and the obstruction is the width, not the amplitude alone.**

## What this does and does not establish

**Shown.** The steady written field of a spherical source, its equivalent mass and density, and its lensing under a declared light rule, computed two independent ways that agree to 3.6×10⁻⁷, anchored to RPG-1's archived lens numbers to 5×10⁻¹⁴ and to RUT-1's kernel exactly. A theorem verified numerically to 10⁻¹⁶: the written field of a bounded source has no far field. And a map, on exposed data, of what a universal (Λ, w) can and cannot do: it can carry a cluster's interior and any single lens, and it cannot carry a cluster's outskirts, Coma's outer shear, or a cluster and a lens together.

**Not shown.** No galaxy or Milky Way comparison (PM-2A C/D); no use of the lenses' stellar kinematics, so no joint motion-and-lensing test; neither light rule derived; no central galaxy in the clusters, no non-thermal pressure, no source geometry for Coma; nothing evolved in time — the response times and formation are RUT-1's business and do not enter a steady state; no unexposed sample opened. The five cluster masses per system are one NFW fit; the Coma bins are a figure reconstruction; the lens masses are population-synthesis values in PF-1's conditional geometry. The 39% agreement among the lenses at 4.64 kpc was found after looking and is not a prediction.

## Next

- **Whether a width that scales with the system is worth declaring** is now a concrete question for goal 3 — what sets the width — and it is the owner's to pose, because RUT-1 stage 0 named it as a different law and goal 8 forbids a label as a switch. The map says what such a law would have to do: w ≈ R_e in a lens, w ≈ 0.4 R₅₀₀ in a cluster, and even then a cluster's outskirts stay unserved unless the field acquires a far field.
- **The joint motion-and-lensing test on the six lenses**, with the KCWI kinematics through CR-2's machinery under L0 and L1, would say whether the amplitude the Einstein radii want at w ≈ R_e is the one the stellar motions want — the test that decides whether the lens regularity means anything.
- **The galaxy-scale comparison** of the written field's own predicted fields with PM-2A stages C and D, where w must be one number for 149 disks and the Milky Way's vertical force.
- **A stellar component with a cusp** for the clusters, since a BCG would write a stronger local field than a gas-shaped stellar model can show.

## Reproduce

```sh
python research_work/results/path-memory/cl1.py           # the nine gates, then E1-E4; about four minutes on one core
python research_work/results/path-memory/cl1_checks.py    # the suite job, about 40 s
```
