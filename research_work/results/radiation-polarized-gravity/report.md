# Radiation-polarized gravity (RPG-1)

**Result.** RPG-1 solved the AQUAL field equation, with μ = x/(1+x) and the archived a*, on the frozen baryons. It is not promoted.

- **Rotation.** It reproduces algebraic simple MOND, slightly worse in RMSE.
- **Lensing.** It fails the declared test: stars would need 2.7–3.9 times their Chabrier mass.
- **Milky Way.** Its vertical force gets worse.
- **Field energy.** It diverges logarithmically.
- **The link to PF-1.** It is a coincidence-level statement that this repository's data cannot test.
- **Validation.** The declared convergence gate failed in one of its three galaxies (0.56 km/s against a 0.5 km/s tolerance). By the protocol's own rule, V1 did not pass.

**Protocol.** [protocol.md](protocol.md), declared at `main` 1e9b9d6 and committed as 693f7c5 before any science output existed.

**Code.**
- [aqual.py](aqual.py): the solver.
- [baryons.py](baryons.py) and [lensing.py](lensing.py): baryon and lens models.
- [rpg1.py](rpg1.py): the driver.
- [checks.py](checks.py): the suite job.
- [convergence.py](convergence.py): the post-hoc diagnostic.

**Results.** [rpg1-results.json](rpg1-results.json) and [convergence-results.json](convergence-results.json).

## Provenance

- **Established mathematics and existing proposals:**
  - the AQUAL equation (Bekenstein and Milgrom 1984) and the simple μ (Famaey and Binney 2005);
  - the exact AQUAL solution for a razor-thin Kuzmin disk (Brada and Milgrom 1995);
  - the exactness of the algebraic relation g = ν(g_N/a*)g_N in spherical symmetry.
- **Project hypotheses:** the field response replaces the reservoir (R1); a* frozen at the archived value (R2); the link to PF-1 (R3); the declared lensing response Φ = Ψ (R4); PF-1's Euclidean geometry (R5).
- **Measured:**
  - SPARC rotation speeds;
  - SLACS Einstein radii (bSIE) and redshifts.
- **Model-inferred:**
  - SPARC baryon surface densities (Υ_disk = 0.5, Υ_bulge = 0.7, exponential gas);
  - the declared disk thickness;
  - Milky Way Eilers speeds (Jeans-inferred);
  - Bovy–Rix K_z (from stellar kinematics);
  - SLACS population stellar masses (Auger et al. 2009, rescaled to PF-1 luminosity distances);
  - all lens distances (PF-1).
- **New here:**
  - the field-equation solution on the frozen baryons;
  - its lensing and vertical-force consequences;
  - the growth of the field energy;
  - the radiation-radius diagnostic.

## V1. Solver validation (declared)

| Check | Tolerance | Result | |
|---|---|---:|---|
| Plummer sphere, exact AQUAL | 1e-6 relative | 1.5×10⁻⁸ | pass |
| Razor-thin Kuzmin disk, exact AQUAL: speed / K_z | 0.5% | 0.14% / 0.05% | pass |
| Miyamoto–Nagai, Newtonian: speed / K_z | 0.5% | 0.16% / 0.06% | pass |
| NGC3741 (gas-richest): doubled resolution / r_out × 10 | 0.5 / 0.1 km/s | 0.002 / 0.0004 km/s | pass |
| UGC05005 (median mass) | 0.5 / 0.1 km/s | 0.007 / 0.002 km/s | pass |
| NGC4217 (most bulge-dominated) | 0.5 / 0.1 km/s | **0.557 / 0.265 km/s** | **fail** |

- Nonlinear residuals are below 9×10⁻¹⁰ everywhere, and grid masses agree with the direct integrals to 2.1×10⁻⁵.
- A rerun reproduced every number exactly.
- By the declared rule, V1 did not pass.

**Post-hoc diagnostic (not declared).** [convergence.py](convergence.py) repeated the gate for all 149 galaxies.

- **Doubled resolution.** 12 galaxies change by more than 0.5 km/s, and all 12 have bulges (12 of 31). The 118 bulgeless galaxies change by at most 0.11 km/s (median 0.011). The largest change is 1.59 km/s (0.55%), in NGC7814 at 0.63 kpc.
- **Ten times r_out.** 15 galaxies change by more than 0.1 km/s (largest 0.99 km/s).
- **Smoothing the bulge helps only partly.** The diagnostic replaced the piecewise-linear bulge profile with a monotone cubic through the same SPARC points. Doubled-resolution failures fall only to 7 of 31, so kinks between tabulated radii are not the main cause. The likely remaining sources, not tested here, are the tabulated bulge's constant-density core inside its first radius and its abrupt edge beyond its last. The interpolation choice alone moves speeds by up to 1.55 km/s.
- **The conclusion survives on the converged galaxies.** In bulge galaxies, this numerical and representation uncertainty (up to about 1.6 km/s) is comparable to the field-equation effect. The 118 bulgeless galaxies converge to 0.11 km/s and give the same conclusion:
  - AQUAL RMSE is 19.12 / 16.73 / 12.39 km/s, against 17.99 / 15.66 / 11.69 for the algebraic relation on the same baryons;
  - the median difference is −2.5 km/s;
  - AQUAL is better in 44 of 118.

## T1. SPARC rotation at frozen a* (all 149 galaxies)

| Model and baryons | RMSE, km/s (train / validation / test) | log RMS |
|---|---:|---:|
| Newtonian, this surface-density model | 54.94 / 60.45 / 49.61 | 0.292 / 0.269 / 0.258 |
| Algebraic simple MOND, same baryons | 21.28 / 28.32 / 16.30 | 0.112 / 0.099 / 0.078 |
| **AQUAL (RPG-1), same baryons** | **22.16 / 29.42 / 16.90** | **0.109 / 0.099 / 0.072** |
| Archived algebraic, SPARC rotation contributions | 19.89 / 26.88 / 16.40 | 0.109 / 0.095 / 0.078 |

- **The field equation lowers the midplane speed** relative to the algebraic relation on the same baryons.
  - Median −2.0 km/s (−1.4%); 5–95% range −6.1 to +0.4 km/s; largest −11.7 km/s.
  - By acceleration: −2.9% median below g_N/a* = 0.1, −1.4% for 0.1–3, and −0.7% above 3.
  - AQUAL fits better than the algebraic relation in 57 of 149 galaxies (30/89, 12/29 and 15/31).
  - The Kuzmin test reproduces the exact thin-disk solution to 0.14%, so the shift is the field equation, not numerics.
  - a* was fitted in the algebraic limit and is deliberately not refitted here, so the field solution under-rotates slightly at the same a*.
- **The baryon model matters as much.** The surface-density baryons and the archived rotation contributions differ by a median of only −0.25 km/s in the algebraic prediction. But the 5–95% range is −9.8 to +3.9 km/s, individual radii differ by up to 80 km/s, and the training RMSE rises from 19.89 to 21.28. That is why the field equation is compared on the same baryons.
- **Thickness.** h_z = 0.2 R_d lowers the AQUAL speeds by a median 0.86 km/s, giving RMSE 22.99 / 30.05 / 17.44.

## T2. Milky Way

The table gives rotation RMSE over the 38 Eilers bins, and K_z/2πG over the 43 Bovy–Rix values at |z| = 1.1 kpc (RMS, bias and χ² with the quoted errors).

| Variant | Model | Rotation RMSE (km/s) | K_z RMS (Msun/pc²) | K_z bias | K_z χ² (43 values) |
|---|---|---:|---:|---:|---:|
| I | Newtonian baryons | 52.56 | 15.6 | +1.2 | 70.6 |
| I | Algebraic | 9.54 | 42.1 | +40.1 | 671.6 |
| I | **AQUAL** | **9.89** | **42.0** | **+40.0** | **667.2** |
| II | Newtonian baryons | 62.32 | 19.1 | −8.9 | 60.6 |
| II | Algebraic | 12.04 | 33.4 | +30.4 | 446.2 |
| II | **AQUAL** | **13.70** | **33.3** | **+30.4** | **443.8** |

- **Consistency.** The Newtonian and algebraic rotation scores reproduce the archived 52.57 / 62.34 and 9.53 / 12.04.
- **The vertical force is over-predicted.** The response that repairs rotation over-predicts the vertical force at 1.1 kpc by 30–40 Msun/pc², and χ² rises from 61–71 to 444–667. The archived conservative-field completion, on a different baryon baseline, moved the same way (27.99 → 33.31).
- **Caveats.** The K_z values are model-inferred in an R₀ = 8 kpc frame. The algebraic K_z column applies ν(|g_N|) point by point and serves only as a reference.

## T3. Lensing: declared response in PF-1 geometry

Each cell gives the predicted Einstein radius divided by the observed one, for Newtonian baryons / RPG-1 / RPG-1 with linked a*. The last column is the stellar-mass factor RPG-1 needs to reproduce the observed radius.

| Lens | z_l | Chabrier | Salpeter | Mass factor, Chabrier / Salpeter |
|---|---:|---|---|---:|
| J0037-0942 | 0.196 | 0.328 / 0.384 / 0.393 | 0.534 / 0.614 / 0.628 | 3.41 / 1.92 |
| J1112+0826 | 0.273 | 0.299 / 0.352 / 0.365 | 0.487 / 0.562 / 0.580 | 3.89 / 2.19 |
| J1204+0358 | 0.164 | 0.415 / 0.473 / 0.481 | 0.649 / 0.731 / 0.743 | 2.74 / 1.54 |
| J1402+6321 | 0.205 | 0.329 / 0.386 / 0.397 | 0.554 / 0.636 / 0.651 | 3.02 / 1.74 |
| J1621+3931 | 0.245 | 0.310 / 0.365 / 0.376 | 0.503 / 0.581 / 0.597 | 3.66 / 2.06 |
| J1630+4520 | 0.248 | 0.360 / 0.428 / 0.443 | 0.578 / 0.674 / 0.694 | 3.06 / 1.72 |

- **The declared rule fails for both a* choices.** It required all six ratios within 0.9–1.1 for at least one IMF.
- **The response is weak where the lenses are.** At the Einstein radii g_N/a* = 9–14, so the response adds only 15–18% to θ_E. Without it, stars would need 3.2–4.6 times their Chabrier mass (1.8–2.6 times Salpeter).
- **PF-1 geometry raises the requirement.** It is larger than in the archived bending budget because PF-1's Euclidean distances make D_l larger and D_ls/D_s smaller than the archive's regular-optics convention. The archive's statement that "stars alone supply 114–118%" used masses fitted to stellar motions in that convention, 2.3–2.7 times the Chabrier values.

## T4. Field energy and the outer boundary

- **The field energy grows without limit.** Between 10 and 100 MOND radii (r_M = √(GM/a*)), dE_F/d ln r is a median 1.031 × M v_f²/3 (range 0.979–1.031; Milky Way 1.016), where v_f⁴ = GMa*.
  - The field energy therefore grows by M v_f²/3 per e-fold of radius.
  - The potential rises by a median 1.034 v_f² ln 10 per decade.
  - Inside 10 r_M the field energy is a median 0.67 M v_f² (range 0.20–2.95). The Newtonian gradient energy of the same baryons converges, at a median 0.27 M v_f².
- **Consequence.** An isolated RPG-1 galaxy has no finite total energy and no escape speed. Its outer boundary must be supplied from outside the equation:
  - an external field;
  - a cosmological cutoff;
  - or the owner's radiation-state argument s.

  The numerical solutions do not depend on where the grid ends, because Gauss's law fixes the flux.
- **The radiation radius falls inside the data.**
  - R_u is the radius at which a galaxy's own starlight energy density falls to the microwave background's 4.17×10⁻¹⁴ J/m³, using a point source with L_bol = L[3.6] (an upper bound; R_u scales as √(L_bol/L[3.6])).
  - R_u lies inside the last measured radius in 140 of 149 galaxies, at a median 0.30 R_last. A response that switched off where starlight falls below the background is therefore excluded by the measured flat outer curves.
  - R_u sits near the MOND radius, at a median 1.31 r_M (5–95%: 0.68–1.74). But R_u/r_M ∝ √(L/M_b): it is 1.63 for the 55 galaxies with gas fraction below 0.3 and 0.76 for the 34 above 0.7 (Spearman −0.98).
  - A transition keyed to starlight would therefore fall at g ≈ 0.38 a* in gas-poor galaxies and ≈ 1.7 a* in gas-rich ones, about 4.6 times apart. One a* fits the SPARC rotation curves.
  - This is a design constraint on any dependence on s, not support for one.

## T5. The link to PF-1

- **The value of ξ.** ξ = a*/(c²α) = 0.1181 with the galaxy-group α, or 0.1263 with the supernova α. This restates the long-known coincidence a0 ~ cH0/2π (1/2π = 0.159); it is not a derivation.
- **The predicted evolution.** A linear index predicts a*(z) = (1+z) a*(0). SPARC and the Milky Way lie at z < 0.03 and cannot test it. For the SLACS lenses (z_l = 0.16–0.27) it raises θ_E by only 1.7–3.6%, far below the lensing deficit.
- **A test needs high-redshift kinematics.** Rotation curves or the radial-acceleration relation at z ≳ 1 (measured kinematics, model-inferred baryons) would test it. None are in the repository.

## Decision

RPG-1 is not promoted, for five reasons:
1. The declared lensing response fails. It predicts 0.35–0.48 of the SLACS Einstein radii with Chabrier masses and 0.56–0.74 with Salpeter.
2. The Milky Way vertical force worsens, from K_z RMS 15.6 to 42.0 Msun/pc² in variant I.
3. Rotation matches algebraic simple MOND and is slightly worse in RMSE. RPG-1 is the field form of that model.
4. V1 failed its declared convergence gate in one galaxy.
5. The field energy diverges, so the outer boundary must come from outside the equation.

**What survives.**
- A field solver validated against exact solutions.
- The link's ξ ≈ 0.12, with an evolution prediction that needs high-redshift data.
- Two constraints on the owner's s-dependence:
  - it cannot switch the response off inside the measured curves;
  - a starlight-keyed transition would differ between gas-rich and gas-poor galaxies by about 4.6 in acceleration.

**What the lensing deficit leaves open.** SLACS Einstein radii lie at g ≫ a*, where any acceleration-keyed response is weakest. The factor of about 3 in Chabrier mass therefore needs one of:
- more stellar mass (a bottom-heavy IMF);
- a lensing response stronger than the dynamical one;
- a different geometry.

## Addendum: how much of the lens deficit is geometry (post hoc, not declared)

[geometry-diagnostic.py](geometry-diagnostic.py) recomputes the Einstein radii with flat FLRW angular distances (H0 = 70, Ω_m = 0.3) and the published population masses. Those are the conventions the masses were derived in. Results are in [geometry-diagnostic.json](geometry-diagnostic.json).

The FLRW history here is an adopted comparison. Co-scaling rulers ([CC-1](../clock-completion/report.md)) give FLRW-form distances for whatever history n(t) the field's dynamics supply; they do not select this one.

- **The ratios rise, but the rule still fails.** RPG-1's ratios go from 0.35–0.47 to 0.48–0.58 with Chabrier masses, and from 0.56–0.73 to 0.75–0.89 with Salpeter.
- **A reason for a joint fit, not agreement.** In that geometry RPG-1 needs 2.1–2.7 times the Chabrier mass (1.2–1.5 times Salpeter). The archive's masses fitted to the same lenses' inner stellar motions are 2.3–2.7 times Chabrier, but they were derived in a regular-optics geometry with its own mass conventions. The similarity is therefore a reason to fit stellar motions and lensing jointly, with identical geometry and mass conventions; it is not evidence that they agree. (Revised 14 September 2026 after review; the first version of this addendum said they agree.)
- **In these conventions the deficit is an inner mass budget.** It needs a heavier IMF or an extra compact component. PF-1's static Euclidean distances make it about 30% worse in Einstein radius.

## Reproduce

```sh
python research_work/results/radiation-polarized-gravity/rpg1.py                  # about 90 s on 8 workers; exits 1 because V1 failed
python research_work/results/radiation-polarized-gravity/checks.py                # suite job: analytic V1 and three galaxies against the archive
python research_work/results/radiation-polarized-gravity/convergence.py           # post-hoc convergence diagnostic
python research_work/results/radiation-polarized-gravity/geometry-diagnostic.py   # post-hoc lens-geometry diagnostic
```
