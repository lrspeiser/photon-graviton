# PM-1: a √M, 1/r force law is a competitive target; one explicit delayed-response realization fails; and the law is not derived by the mechanism that motivated it

Six candidate families for "gravity remembers where matter has been", scored against the data already in this repository: SPARC's 149 galaxies on the archived baryon model and the frozen 89/29/31 split, and the Milky Way's 38 Eilers bins on both fiducials. One universal constant per rule, fitted on the training galaxies only. The owner's review of the first run supplied six corrections, all carried in [protocol.md](protocol.md) and applied here.

> **Correction 7, after the owner's review of the landed run.** The mass this experiment uses is `M_force = R·v_bar²/G = R²g_N/G`, the mass a spherical Newtonian source would need to produce the baryon model's radial force — an enclosed mass only in spherical symmetry, and these are disks. Substituting it into the fitted rule gives, identically and for the disk inputs too, **g = g_N + √(a\*·g_N)** with a\* = β²/G, verified to 4×10⁻¹⁶ across all 149 galaxies. **PM-1 therefore fitted a local acceleration law, a pointwise function of the local Newtonian field, not a force sourced by the matter inside r.** No observed speed entered the predictor — `v_bar` is the ordinary-matter model's own prediction — so the scores stand, as scores of a one-parameter empirical acceleration law. Three labels change with it: the mass–speed regression is against force-equivalent mass rather than an independently integrated stellar-plus-gas mass, so the advertised baryonic mass–speed test has not yet been done; the collective control's "total" is `M_force` at the last sampled radius, which depends on where the curve stops; and the √8 demonstration is a mathematical warning about per-ring saturation, not a physical annulus-convergence test. The source-mass version is PM-2A stage A; this archive is preserved as the force-proxy version.

- **Delay alone explains nothing (A).** An exponential memory of density returns αΦ_N for a stationary source, and a rigidly rotating axisymmetric density is indistinguishable from a static one to 4×10⁻¹¹. A memory of density cannot see an orbit.
- **A single remembered wavelength cannot work (B, control).** On UGC02953's 115 radii, k J₀(kR)J₁(kR) changes sign 35 times.
- **A broad spectrum gives a competitive force law.** With the amplitude ∝ √M_b(<r): **20.21 / 27.74 / 18.07** on SPARC against simple MOND's 19.89 / 26.88 / 16.40 — that is 1.6%, 3.2% and 10.2% higher, and whether those gaps are meaningful needs paired resampling that this run does not do. On the Milky Way it is 18.25 / 9.51 against 9.53 / 12.04: better on fiducial II, much worse on I, so the law is sensitive to the adopted ordinary-matter model.
- **The mass–speed relation separates the sourcing rules.** Amplitude ∝ M gives a slope of 0.412 against the observed 0.289 and fails; ∝ √M gives 0.280; the collective rule gives **0.2886**, the closest of the three.
- **The acceleration scale is a re-expression of the fitted amplitude, not a second prediction.** g_mem = √(a\*·g_mono) identically with a\* = β²/G, so fitting β *is* fitting a\*, and v_f⁴ ∝ M was built into the √M choice. The fitted value, 6.54×10⁻¹¹ m/s² (0.54 of MOND's canonical scale, 0.76 of this repository's own fit), is a consistency, not a discovery.
- **The mechanism does not produce the law that was fitted.** Independent per-ring saturation gives Σ√mᵢ, not √(Σmᵢ). Refining the same galaxy 8× multiplies that field by **2.8284 = √8 exactly**: a representation artefact, not physics. The cumulative rule is a fitting law.
- **A collective field does produce it.** The spherical limit of ∇·[(|∇ψ|/a\*)∇ψ] = 4πGρ_b gives g_ψ = √(G a\* M_b(<r))/r, exactly the fitted form, with β = √(G a\*) matching the fitted β to zero relative difference. The square root comes from one shared nonlinear field, not from a square-root charge per ring.
- **The delayed-response realization fails analytically (C).** Routh–Hurwitz on Tλ³ + (1+q₀)λ² + Tλ + (1+3q₀) requires q₀ < 0, so every q₀ > 0 at finite T > 0 has a growing mode. The owner's root reproduces to five figures (0.15334 ± 1.28578i). At real SPARC outer radii the median q₀ is 0.80, an e-folding of about one orbit — but the growth rate vanishes as T → 0 and T → ∞, so the survival question needs the actual memory time.
- **C's failure does not transfer to B-root**, whose frozen potential −GM/r + K ln r has κ² = GM/r³ + 2K/r² > 0: its circular orbits are radially stable. Whether they survive when matter and field evolve together is the next experiment's question.
- **Energy deposition is far short (E),** by a median 1.7×10⁹ (range 1.4×10⁸–1.2×10¹¹) — conditional on converting the catalogue's 3.6 μm band luminosity with the bolometric solar constant, which is an assumed spectral conversion and is labelled as one.

Protocol: [protocol.md](protocol.md), declared in dcddb70 before execution, with the owner's corrections declared before the corrected run. Code: [pm1.py](pm1.py). Results: [pm1-results.json](pm1-results.json).

## The three sourcing rules, scored the same way

| Rule | one constant | SPARC train / validation / test | Milky Way I / II | mass–speed slope (obs 0.2894) |
|---|---|---|---|---|
| Baryons alone | — | 52.57 / 58.22 / 47.77 | 52.57 / 62.34 | — |
| Simple MOND, fitted a₀ | a₀ | 19.89 / 26.88 / 16.40 | 9.53 / 12.04 | — |
| **B-root**, amplitude ∝ √M(<r) | β = 9.3147×10⁻² | **20.21 / 27.74 / 18.07** | 18.25 / **9.51** | 0.2804 |
| **Collective**, one shared mode ∝ M(<r)/√M_tot | β = 1.0399×10⁻¹ | 21.66 / 29.67 / 18.71 | 20.90 / 12.85 | **0.2886** |
| B-linear, amplitude ∝ M(<r) | κ = 3.2746×10⁻⁷ | 30.53 / 34.03 / 29.46 | 20.45 / 10.30 | 0.4120 |

The collective rule is the one with a derivation behind it — a single saturated mode whose weight follows the mass — and it costs no extra parameter. It fits the curves slightly worse than the cumulative rule and the mass–speed relation slightly better. Neither is preferred by the data at the precision of this test.

## The two audits the owner's review required

**The ring force.** The on-ring expression (C/2R)[J₀(q_min R)² − J₀(q_max R)²] is exact only at r = R_s; the field of a ring at R_s observed at r is C∫J₁(qr)J₀(qR_s)dq. At r = 2, R_s = 1 over [0.01, 100] those are 0.49877 and 0.24989 — a factor of two, reproducing the owner's example to five figures. The first run used the on-ring form at the observation radius, which silently moves the source; that row is withdrawn and replaced by this audit. The wide-band limit is verified directly rather than assumed: the ring integral is 1/r outside (0.4988 at r = 2, 0.2001 at r = 5) and near zero inside (0.011 at r = 0.5, −0.005 at r = 0.2). So the enclosed-mass rule *is* the correct wide-band limit; it is the finite-band claim that was wrong.

**The representation test.** Splitting each annulus of UGC02953 into 1, 2, 4 and 8 equal sub-rings and summing √mᵢ gives ratios to √(Σmᵢ) of 8.83, 12.49, 17.66 and 24.98 — a growth of 2.8284 over the 8× refinement, exactly √8. Independent ring saturation therefore cannot be the mechanism behind the fitted law, whatever its merits as a fit.

## What this does and does not establish

**Shown.**
- A one-amplitude, √M, 1/r force law is a competitive description of 149 rotation curves under the frozen baryon model, with no per-galaxy freedom.
- The baryonic mass–speed relation, not the curves, discriminates the sourcing rules.
- The square root is derivable — from a collective nonlinear field equation, not from per-ring saturation, which is representation dependent by √N.
- The specific delayed-receiver realization has a growing mode for every positive memory time, analytically.

**Not shown.**
- **That the force has path memory.** The fitted equation contains enclosed mass, radius and a constant. No memory time, history or phase appears in it. A successful static fit is not evidence of memory.
- **An independent prediction of the acceleration scale.** It is the fitted amplitude, re-expressed.
- **That B-root beats MOND.** It is 1.6–10.2% worse on the three SPARC splits and much worse on Milky Way fiducial I, and no paired uncertainty was computed.
- **Anything about a disk.** B-root's enclosed-mass prescription is not the Newtonian radial field of a disk, and the nonlinear field's disk solution will not equal it. The repository's AQUAL solver (RPG-1) is the tool for that comparison, and it was not run here.
- **Candidates D and F**, declared and not run: D needs a thickness of tens of kiloparsecs SPARC does not measure; F needs a trajectory-history integrator.
- **No blindness.** These data are exposed and the split is reused; per [[holdout-provenance]] this is validation, not confirmation.

## Reproduce

```sh
python research_work/results/path-memory/pm1.py       # about a minute on one core
python research_work/results/path-memory/checks.py    # suite job: the control, the three rules, both audits and the anchor
```
