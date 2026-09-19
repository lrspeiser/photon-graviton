# CL-2 stage 1: a universal spectrum of footprints describes the clusters' pressure profiles nearly as well as their own NFW fits, cannot describe the galaxies or the lenses, and no member serves all three

Protocol [protocol-cl2.md](protocol-cl2.md) (14e4b34, declared before any run); [amendment 1](protocol-cl2-amendment-1.md) (80c4bef, the held-out lens withdrawn before any run); [amendment 2](protocol-cl2-amendment-2.md) (b932d78, after gate G4 failed as declared and before its replacement was run). Inputs `cl2-inputs-xcop-profiles.json` from `cl2_xcop_acquire.py`; operators `cl2_response.py`; blocks `cl2_sources.py`; driver `cl2.py`; archives `cl2-results-as-declared.json` (the run as declared, numerical verification FAILED) and `cl2-results.json` (the run as amended); geometry registry `cl2-geometry.json`; suite job `cl2_checks.py`. CL-1's files are imported and unchanged.

This is the reviewer's first milestone run: given independently specified ordinary matter, what field does the written-track mechanism predict for a galaxy, a lens galaxy and a cluster with one law and common parameters, and does any member of a universal family of footprint widths serve all three. The family is the nonnegative spectrum C = Σ_j Λ_j (ρ_b ∗ K_{w_j}) over 26 widths from 0.1 kpc to 10 Mpc, the same for every system, built from CL-1's verified per-width operators. The clusters are compared through the **measured pressure**: with the X-COP electron density, a model gravity predicts the pressure profile up to one boundary value per cluster, against 246 X-ray and SZ pressure points in 12 clusters, with the release's own NFW mass as the floor no model can be asked to beat.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | **passes.** `cl2_checks.py` reruns the eight gates with their controls (G4 and G6 at reduced size) and recomputes ten anchored numbers from scratch; all match the archive to 10⁻⁹. 83 s. |
| **numerical verification** | **FAILED as declared; passes as amended, and both are quoted wherever it is.** G4 required every galaxy operator to change by less than 10⁻³ when the source grid is doubled; the narrowest width, 0.1 kpc, changes by **3.4×10⁻²** on UGC12506 and exceeds 10⁻³ on eight galaxies, because a 0.1-kpc footprint reads the slope discontinuities of the interpolated photometry rather than a resolved field, the same limit that excluded widths below 10 kpc from the clusters. The declared run is archived as it fell. Amendment 2, declared before the rerun, excludes widths below 0.15 kpc from the galaxy block (they stay in the family, constrained by the lenses, whose analytic light profiles are converged to 3×10⁻⁶); G4b then passes with the retained galaxy columns at **5.9×10⁻⁴** and the clusters at 5.2×10⁻⁵. Every other gate passes with its control rejected: the two analytic members to 1.0×10⁻⁷ and 6.9×10⁻⁹ (the √π control off by 29%); CL-1's archived best reproduced exactly by the one-width spectrum (a second width moves χ² from 4,383 to 3,213); the two lensing routes on a three-width mixture to 1.4×10⁻⁷ (the shell-kernel control off by 6%); the archived galaxy scores 52.57, 20.21 and 19.89 km/s reproduced to 0.01 (Υ = 1 moves the baryon score to 38.4); CR-2's stars-only χ² and β reproduced to 1.7×10⁻⁷ on all six lenses (the PSF doubled moves χ² by 93%); the release's gas mass reproduced to 1.0–4.1% and its NFW floor at 5.75 per point (μ doubled sends it to 349); the solve's optimality residual 2×10⁻¹⁶. |
| **scientific outcome** | **(c) by the declared rule: even flexible positive spectra fail — within the tested sources, bandwidth and light rule, changing widths is not enough.** The clusters are described (8.55 per point against the floor's 5.75, a factor 1.49 where 2 was allowed); the galaxies are not (31.0 km/s where 21.9 was required, and a mass–speed slope of 0.388 against the observed 0.289); the lenses are not (Einstein radii off by up to 13% where 3% was allowed, and stellar-motion χ² 408 where 128 was allowed). Jointly the same three verdicts hold. The prototype's expectation, quoted in the protocol, is what the run found. |

In plain words: give the mechanism every footprint size at once and let the data choose. The clusters choose five sizes between 60 and 630 kpc and their pressure profiles come out nearly as well as with a dark halo fitted to each. The galaxies choose sizes below 100 kpc and still cannot get the rotation curves right, because a field that is proportional to the matter that writes it can never make speed grow as the fourth root of mass, which is what galaxies do. The lens galaxies cannot satisfy their light bending and their stellar motions with one field. The next step is therefore not another width; it is a response that is not proportional to its source.

## E1 — the universal spectrum, block by block and jointly

| block | points | reference laws | the spectrum alone | described? |
|---|---|---|---|---|
| clusters, forward pressure | 246 | release NFW floor **5.75** per point; baryons 211; PM-1's law 90.7 | **8.55** per point, five widths 63–631 kpc (w/ℓ 2.3–9.3) | **yes** (≤ 2 × floor) |
| SPARC training galaxies | 1,989 | baryons 52.56 km/s; PM-1's law 20.21; simple MOND 19.89 (slope 0.272) | **30.99 km/s**, slope **0.388**; validation 33.6, test 28.0 | no (needed ≤ 21.9 and a slope within 0.05 of 0.289) |
| six SLACS lenses, Einstein + motions | 6 + 40 | CR-2 free NFW motions χ² 85.3; stars-only 883 | Einstein residuals −9.3% to +13.1%; motions χ² 407.7 | no (needed 3% and ≤ 128) |
| **joint** | 2,281 | — | clusters 8.83; galaxies 31.8 km/s (val. 35.0, test 30.7); Einstein to −22%, motions 1,306 | clusters yes; galaxies no; lenses no |

**The clusters.** One spectrum shared by twelve clusters, five active amplitudes plus one boundary pressure each, sits a factor 1.49 above twelve individually fitted NFW halos. The declared sensitivities move it little: μ = 0.59 or 0.61, 8.55; the non-thermal pressure fractions applied, 8.59; the four- and six-per-decade and the wider width grids, 8.89–8.93; no stars at all, 9.15. Per cluster it is within 1.3 × the floor for nine of the twelve, and worst for RXC1825 (23.8 against 4.1), one of the five without a measured stellar profile. **Where the force comes from** (E4): at every pressure radius, 72–77% of the written force comes from source shells between r/2 and r, 29–43% from inside r/2, −5% to −15% from r to 2r, and nothing from beyond 2r; truncating the gas at R₅₀₀ moves the score from 8.55 to 8.84 at fixed amplitudes and 8.50 refitted. So the reviewer's question is answered by measurement: the response is local-to-interior, and a cluster's outer pressure is carried here because the pressure is an integral of the force from the outside in, not because the field has a far field. On CL-1's own five-point NFW comparison, with the release's baryons, the cluster spectrum scores 2,135 against CL-1's single width 4,383, PM-1's law 28,309 and the baryons' 86,192.

**The galaxies.** A nonnegative combination of widths from 0.16 to 100 kpc halves the baryons' misfit (52.6 → 31.0 km/s) and stops far short of the MOND-like law (19.9), and its mass–speed slope, 0.388, is PM-1's B-linear failure again (0.412 there; the observed 0.289; the √M law 0.280). This is not a numerical accident and the width grid does not change it (31.0 on every grid): a field linear in its source gives v² ∝ M, and the data want v⁴ ∝ M. The declared run, which still held the 0.1-kpc column, scored 30.95: the unresolved column had taken the role the 0.16-kpc column now takes, and the reading is the same.

**The lenses.** No positive spectrum, with population stellar masses, satisfies both observables. The solver's best compromise leaves Einstein radii off by 13% and the motions at 408 (CR-2's free NFW: 85). The amplitude CL-1 found the six lenses agreeing on at w = 4.64 kpc, to a factor 1.39, receives **no weight at all** in any lens fit; the fits put their weight at 0.4–0.6 kpc and at 25–40 kpc instead.

## E2 — the lens tests: light bending against stellar motions

| fit (Chabrier masses) | Einstein radii | stellar motions, χ² over six lenses | what it predicts |
|---|---|---|---|
| lensing only | within 9.6% (widths 4 and 100 kpc) | 1,923 (per lens 627, 109, 488, 369, 116, 214; CR-2's stars-only 25, 388, 6, 69, 50, 344) | a field that bends the light spoils the motions for four of six |
| stellar motions only | **6.6 to 21 times too large** (the fit puts w/ℓ of 10⁵–10⁶ at 1.6–2.5 Mpc: a harmonic force across the galaxy) | 323 | a field that fits the motions bends light absurdly |
| joint | within 13% | 408 | the compromise above |

With Salpeter masses the same three rows read 10%/767, 25–81×/339, and 14%/424. The held-out lens is not available (amendment 1). **The verdict on CL-1's regularity:** it was a regularity of one observable at one width; once the stellar motions are included the width carries nothing, and the two observables do not want the same field. That is what the reviewer said the joint test would decide, and it decided against.

## E3 — transfer without retuning

| calibrated on | clusters (per point; floor 5.75) | galaxies, train / validation / test (km/s) | lenses, Einstein residuals | Milky Way I / II (km/s; baryons 52.6 / 62.3) |
|---|---|---|---|---|
| galaxies | 157 (baryons 211) | — | −16% to −39% | 25.9 / 12.5 |
| clusters | — | 42.5 / 45.9 / 38.2 (baryons 52.6 / 58.2 / 47.8) | −50% to −60% (the Newtonian bracket) | 41.2 / 51.9 |
| galaxies and clusters | 8.96 | 31.0 / 33.6 / 28.0 | −16% to −38% | 25.8 / 12.4 |

The two scales do not interfere and do not help each other: the widths the clusters use (60–630 kpc) are invisible in a galaxy, and the widths the galaxies use (≤ 100 kpc) do nothing for a cluster's pressure. The one transfer worth a second look is the Milky Way on baseline II, 12.4 km/s from a spectrum that never saw it (simple MOND: 12.0; PM-1's √M law: 9.5), against 25.8 on baseline I (MOND 9.5).

## What was corrected on the way, and what it cost

- The held-out lens J1538+5817 has no published light-profile components; the audit says so and forbids a fallback. I declared it without checking. Withdrawn before any run (amendment 1); the stage has no held-out system.
- Gate G4 failed as declared on the 0.1-kpc galaxy width. The declared run is archived unchanged; the replacement was declared before it was run; the cause is a resolution limit of the source, identical in kind to the one already declared for the clusters, and the galaxy reading it blocked (30.95 km/s) is reproduced by the amended run (30.99) with the resolved widths.
- The pressure comparison uses the release's electron density and its symmetric errors, diagonal; the SZ covariance in the release is not used yet. Non-thermal pressure and μ enter as sensitivities, not as parameters.

## What this does and does not establish

**Shown.** The reviewer's structure works and runs end to end: one family, one solve, three classes of system, each block anchored to something already archived. The clusters' pressure profiles are described by five universal footprints nearly as well as by their own halos, with the force measured to come from inside each radius. The galaxies cannot be described by any universal linear convolution, for the structural reason PM-1 already exposed. The lenses' two observables want different fields. Therefore, by the rule fixed before the run, changing widths is not enough: the next ingredient must be a response that is not proportional to its source — nonlinear, or transported, or a different light coupling.

**Not shown.** No formation, energy or propagation at the inferred strengths (CL-2E); no Coma amplitude; no held-out system; no Milky Way vertical force; no recomputation of the disks' Newtonian field from the reconstructed source; neither light rule derived. The cluster success is exposed-data description with 17 parameters against 246 points, on hydrostatic pressures with the release's modelling inside them, and "described" is a declared threshold, not a claim that the mechanism is right. The clusters' five widths and the galaxies' eight are disjoint; a universal law that used both would be two laws with one name.

## Next

The decision rule points at three ingredients and the map says what each must do: a nonlinear response, which can make v⁴ ∝ M in galaxies (PM-2A's completions do, as the reproduced Bekenstein form) while leaving the clusters' flat boost to be earned by something clusters have and galaxies lack; transport, which would give a far field the local footprint cannot; and the light coupling, which the lens test now constrains directly, since the motions and the bending disagree under L1 with population masses. Each belongs in its own protocol with universal parameters fitted on a declared subset and carried, and the machinery here runs any of them without change. Formation and the response times stay deferred, as the plan says, until a static response survives.

## Reproduce

```sh
python research_work/results/path-memory/cl2_xcop_acquire.py    # rebuilds the extract from the cached or downloaded release (needs astropy)
python research_work/results/path-memory/cl2.py                 # the gates, then E1-E4; about twelve minutes
python research_work/results/path-memory/cl2_checks.py          # the suite job, about 80 s
```

## Corrections after review (19 September 2026; text only, no archived number changed)

The owner's review of bddb6df was checked line by line against the code and the archive; every point held, and two are more consequential than the review stated. The archive `cl2-results.json` is unchanged; the sentences below correct this report's reading of it.

1. **The "4.64-kpc amplitude" is the 10-kpc column.** `Lambda_at_4p64_kpc` reads `amps[10]`, and the 26-width grid has no 4.64-kpc column (its neighbours are 3.98 and 6.31 kpc). The fits in fact put weight at 3.98 kpc: 2.5×10⁻⁶ in the lensing-only Chabrier fit, against CL-1's 3.5–4.9×10⁻⁶ at 4.64 kpc, and 1.1×10⁻⁶ once the motions are included. The sentence "CL-1's regularity at 4.64 kpc receives no weight" is withdrawn: the neighbouring width carried weight in every fit that includes the bending, halved when the motions are added.
2. **The joint fit did not test one law.** Cluster columns below 10 kpc and galaxy columns below 0.15 kpc are zero while the lenses see every width. The joint fit's largest amplitude, 5.0×10⁻⁴ at 0.1 kpc, sat in a column two of the three blocks could not see; had the galaxies seen it, the added force would have been 18% of Newtonian at the median point and a 10 km/s speed change per galaxy at the median. The clusters' omitted narrow terms are at most 0.09 of one error bar, so the cluster verdict stands. The sentence "a universal law that used both would be two laws with one name" is withdrawn: one spectrum can hold several scales; the defect was the system-dependent zeroing.
3. **Outcome (c) is the output of the declared procedure, not a certified exclusion.** The solver minimised error-weighted acceleration residuals; the rule judges the equal-galaxy speed error, a different objective whose minimum was not computed. The mass–speed slope used the force-equivalent mass R²g_N/G at the last point, not the integrated baryonic mass, and the Newtonian force came from the SPARC tables while the written force came from the reconstructed components. The standalone galaxy and lens failures are therefore failures of the returned solutions under these choices; a certified minimum of the acceptance metric is CL-2 stage 2's task.
4. **The non-thermal sensitivity's reference is wrong.** `fit_boundary` multiplies its model by the thermal fraction that `pressure_of` has already applied, so the reference models in that sensitivity carry the fraction twice; its "floor 10.16" is not a floor. The primary comparison (fraction one) and the spectrum's 8.59 are unaffected. The corrected reference is recomputed in stage 2.
5. **Wording.** The lens "Einstein residuals" are deflection mismatches at the observed Einstein radius relative to the needed deflection, not ratios of predicted to observed Einstein radii; and 5.75 per point is the release's NFW reference under an approximate pressure comparison, not a floor no model can beat: "within a factor 1.49 of the reference" is the supported statement.
