# CR-2 report: one supported condensate does not fit the lenses and the Milky Way together

**Result.** Tested: a Thomas–Fermi condensate with pressure P = Kρ²/2 and one shared constant K, in equilibrium with each lens's stars. It is fitted jointly to KCWI stellar motions and exact lensing, with identical geometry and mass conventions. It fails the declared rule in both geometries:
- **Shared core.** The best shared size is R_TF ≈ 80 kpc in both geometries: 79.6 kpc in flat FLRW, 81.0 kpc in the co-scaling history.
- **Lenses.** The six-lens χ² is 113.6 (FLRW) and 89.3 (co-scaling). Free NFW halos score 85.3 and 59.1 on the same data. The condensate is 28 and 30 higher, with five fewer parameters; the declared allowance was 10.
- **Milky Way.** The same R_TF cuts the rotation RMSE from 52.6 km/s (baryons alone) to 22.5 km/s, and the inner 20 bins from 45.7 to 25.4. That misses the declared limit of 20 km/s.
- **Supply fails separately.** Each lens's starlight could supply at most 10⁻⁶–10⁻⁹ of the condensate mass it needs.
- **Population stellar masses do not fit.** With Chabrier masses the best χ² is 1,608 (FLRW) and 2,070 (co-scaling). The stellar masses the lens fits require are 1.5–3.1 times Chabrier (0.8–1.7 times Salpeter).

Protocol: [protocol.md](protocol.md), declared in a5b2c41 before execution. Code: [tf.py](tf.py) (equilibrium) and [cr2.py](cr2.py) (fits). Results: [cr2-results.json](cr2-results.json). The suite job is [checks.py](checks.py).

## Validation (declared tolerances)

| Check | Result | Tolerance |
|---|---|---|
| V1: without a host the profile is sin(x)/x | edge π to 3×10⁻¹¹, mass π to 6×10⁻¹¹ | 10⁻⁸ |
| V2: the archive's FLRW stars-only and NFW fits, all six lenses | reproduced to 1.9×10⁻¹³ | 10⁻⁴ |
| V3: lens constraint at every reported fit, independent quadrature | closes to 2.1×10⁻⁹ | 10⁻⁸ |

V2 shows the data, geometry and likelihood are identical to the archive's. The co-scaling benchmarks are new and come from the same code.

## Joint fits (M1: stellar mass set by the lens)

| Geometry | Stars only | Free NFW (18 parameters) | Condensate (13 parameters) | Best R_TF | Rule |
|---|---|---|---|---|---|
| G1, flat FLRW (adopted comparison) | 882.9 | 85.3 | 113.6 | 79.6 kpc | fails: 113.6 > 95.3 |
| G2, co-scaling V = 0 history | 1,090.4 | 59.1 | 89.3 | 81.0 kpc | fails: 89.3 > 69.1 |

How the joint χ² varies with the shared R_TF:

| R_TF (kpc) | 0.3–3 | 5.3 | 9.5 | 17 | 30 | 53 | 95 | 169 | 300 |
|---|---|---|---|---|---|---|---|---|---|
| G1 | 847–868 | 876 | 371 | 231 | 153 | 116 | 124 | 135 | 155 |
| G2 | 1,066–1,085 | 840 | 343 | 241 | 169 | 105 | 91 | 114 | 135 |

- **Compact cores barely help.** With R_TF below a few kpc the six lenses stay near their stars-only totals.
- **An extended condensate helps.** Its edge falls at 38–78 kpc, compressed from R_TF by the stars, and it carries up to 43% of the lens deflection.

Per lens at the best shared R_TF, flat FLRW. The co-scaling history is similar; see the results file.

| Lens | Condensate χ² | Stars only | NFW | Condensate's share of the deflection | M*/M_Chabrier | Condensate mass (Msun) | Edge (kpc) |
|---|---|---|---|---|---|---|---|
| J0037-0942 | 18.3 | 25.0 | 20.2 | 0.07 | 2.8 | 2.4×10¹¹ | 38 |
| J1112+0826 | 13.2 | 388.4 | 11.6 | 0.40 | 1.9 | 8.1×10¹² | 75 |
| J1204+0358 | 5.8 | 5.8 | 4.1 | 0.001 | 2.4 | 1.0×10⁸ | 2.3 |
| J1402+6321 | 55.7 | 69.4 | 43.1 | 0.08 | 2.3 | 5.2×10¹¹ | 46 |
| J1621+3931 | 13.7 | 49.8 | 4.8 | 0.20 | 2.4 | 3.0×10¹² | 68 |
| J1630+4520 | 6.9 | 344.4 | 1.6 | 0.42 | 1.5 | 7.5×10¹² | 75 |

**What drives the failure.**
- **Three lenses carry almost all the excess over NFW:** J1402 (+12.7), J1621 (+8.9) and J1630 (+5.3). Their NFW fits want scale radii from 0.01 to 100 times the effective radius. One shared polytrope size cannot give all of those profiles at once.
- **J1204 wants no condensate at all.** J0037 does slightly better with the condensate than with NFW.

## Population masses (M2, reported, not scored)

- **With Chabrier masses** the stars supply only 28–41% of each lens's deflection. The condensate that supplies the rest cannot fit the motions: χ² 1,608 (G1) and 2,070 (G2).
- **With Salpeter masses** the stars supply 49–74%, and the χ² is 396 (G1) and 463 (G2). That is still far above the lens-determined fits.
- **In both geometries** the motions and the lens together want stellar masses of 1.5–3.1 times Chabrier. The archive found the same in its own geometry. With identical conventions, then, this is a stellar-mass (IMF) requirement, whatever the extra component is.

## Milky Way with the same R_TF

- **Setup.** The condensate mass is fitted to the 38 Eilers speeds, with spherical baryon model I for the equilibrium and the archive's baseline-I baryon speeds.
- **Result.** At R_TF = 79.6 (81.0) kpc it needs 1.6 (1.7)×10¹² Msun, with its edge at 75 (77) kpc. RMSE falls from 52.6 to 22.5 (22.6) km/s, and the inner 20 bins from 45.7 to 25.4 (25.5).
- **Comparison.** That is better than baryons alone, but worse than the archived references: 6.8 exact-third, 9.3 MOND-guided and 9.5 simple MOND. It misses the declared 20 km/s.

## Supply, formation and validity (separate results)

- **Supply.** The upper bound on what each lens's starlight could deliver through the conversion law in 10 Gyr is M_supply ≤ αLR_edgeT/c², with L ≤ M*_Chabrier L☉/M☉. It falls short of the condensate mass needed by 1.6×10⁶ to 2.1×10⁹ in G1 and 6×10⁶ to 3.4×10⁹ in G2. If a supported condensate is the reservoir, it cannot come from the lens's own starlight.
- **Formation history.** Not modeled.
- **Thomas–Fermi validity.** The healing length over the condensate edge is at most 0.0065 for a 10⁻²² eV constituent, so the limit holds. For 1.34×10⁻²⁴ eV it reaches 0.48 in J1204's small condensate and 0.01–0.03 elsewhere, so the limit fails there.

## What this means

- **One supported law with one shared size does better than stars alone.** It needs the lens-determined, roughly Salpeter-like stellar masses.
- **It does worse than free halos.** Free NFW halos with a scale per lens fit these lenses better, even after the Akaike penalty for their five extra parameters.
- **The same size improves the Milky Way, but not enough.** The declared rule fails in both geometries.
- **What stands out is that the best size is about 80 kpc in both geometries.** The lenses do constrain an extended scale, consistent across two distance conventions, but this law's profile shape is not flexible enough.

Not tested:
- other repulsive laws, such as a polytropic index other than 1 or a finite-temperature component;
- anisotropy beyond constant β;
- clusters.

## Post hoc: universal against universal

Declared after the results, at the project owner's direction. The same laws and constants should hold everywhere, with differences coming from each system's specified inputs and environment, not from per-system corrections.

The declared rule measured the condensate against NFW halos tuned lens by lens; its allowance of 10 was an Akaike penalty for their extra parameters. The comparison below matches freedoms instead. Both comparisons stand, because they answer different questions.

An NFW halo with one scale radius for every lens is a restricted benchmark, not NFW's own universal law, in which the scale follows each object's mass and assembly history. That version is not tested here.

[universal-halo.py](universal-halo.py) makes the like-for-like comparison ([universal-halo.json](universal-halo.json)):
- an NFW halo with one scale radius shared by all six lenses;
- the same freedom as the condensate: a per-lens share and constant β;
- exact lensing, and the same data, geometry and likelihood;
- the lens-fixed scale then carried to the Milky Way with only the mass free.

| | Condensate, one shared R_TF | NFW, one shared r_s (restricted benchmark) | Free NFW, tuned per lens |
|---|---|---|---|
| Parameters | 13 | 13 | 18 |
| FLRW: six-lens χ² | **113.6** (R_TF = 79.6 kpc) | 128.6 (r_s = 251 kpc) | 85.3 |
| Co-scaling: six-lens χ² | **89.3** (R_TF = 81.0 kpc) | 105.7 (r_s = 305 kpc) | 59.1 |
| Milky Way RMSE, 38 bins (inner 20), km/s | 22.5 (25.4) | **14.1 (13.4)** | not a universal law |

- **Lenses.** At equal freedom the condensate fits better, by 15.0 (FLRW) and 16.4 (co-scaling). It wins in J0037, J1112 and J1402, by 3–14. The shared-scale NFW wins in J1621 and J1630, by 2–5. J1204 is a tie.
- **Milky Way.** The shared-scale NFW does better. Its advantage comes from the interior profile inside the measured 5–25 kpc, analyzed in the next section. The first version of this report blamed the condensate's edge near 75 kpc. Material beyond 25 kpc cannot change these speeds, because the model adds GM(<r)/r to v², so that explanation was wrong.
- **Neither is as good as the archive's frozen references** in the Milky Way (6.8, 9.3 and 9.5 km/s). Those references use their own universal settings, not lens-fixed ones.

The declared verdict is unchanged. At equal freedom, the supported law describes the six lenses better than the restricted NFW benchmark. The χ² values measure the stellar-motion discrepancy with the lens angles imposed exactly, so this is a better joint fit, not an independent prediction of lensing.

## Post hoc: what the Milky Way comparison measures

Requested in the owner's review. [milky-way-interior.py](milky-way-interior.py) compares, inside the measured range, the extra enclosed mass each fitted model supplies with the extra mass the observed speeds require, M_req(<r) = r(v_obs² − v_b²)/G. Results: [milky-way-interior.json](milky-way-interior.json). The FLRW-fixed scales are shown; the co-scaling ones agree to 1%.

| Across 5–25 kpc | Extra mass supplied/required at 6, 10, 15, 20, 24 kpc | Enclosed-mass slope, 8–20 kpc | Mean signed residual (km/s) at 5–10, 10–15, 15–20, 20–25 kpc |
|---|---|---|---|
| Required by the data | 1 | 1.31 | — |
| Condensate (R_TF 79.6 kpc) | 0.19, 0.38, 0.65, 1.16, 1.94 | 2.71 | −27, −24, −9, +21 |
| NFW with one shared r_s (251 kpc) | 0.51, 0.67, 0.83, 1.19, 1.73 | 1.93 | −15, −12, −3, +18 |
| Baryons only | 0 | — | −38, −52, −61, −56 |

- **The condensate is too uniform across the disk.** Its core (~80 kpc) is much larger than 25 kpc, so its density is nearly constant there, and its enclosed mass grows as r^2.7. The data need about r^1.3. It supplies too little inside about 18 kpc and too much beyond.
- **Its edge plays no role.** Only 11–12% of its mass lies inside 25 kpc, and nothing outside 25 kpc affects these speeds.
- **The shared NFW has the same problem, more mildly.** Its density falls as about r⁻¹ inside its scale, so its enclosed mass grows as r^1.9.
- **What the Milky Way asks for** is extra mass concentrated like the requirement, with density falling roughly as r^−1.7 across the disk. A Thomas–Fermi condensate with the lenses' constant cannot be both extended enough for the lenses and concentrated enough for the Milky Way. That tension belongs to this law; it is not an edge effect.

## Leads (hypotheses for declared tests, not results)

- **A universal scale.**
  - One constant takes the six lenses from χ² 883 (stars only) to 114 in FLRW and from 1,090 to 89 in co-scaling. That is about 96–97% of the χ² improvement over stars alone that free per-lens NFW halos give, with five fewer parameters.
  - Its size, about 80 kpc, is the same in both geometries, and it also improves the Milky Way.
  - A law with one shared scale is worth pursuing further.
- **The Milky Way failure is the interior slope, not the edge.**
  - A revised law must concentrate mass inside the Milky Way's disk while staying extended in the lenses.
  - A two-phase law, with a Thomas–Fermi core inside a thermal component, as in superfluid dark-matter models (Berezhiani & Khoury, arXiv:1507.01019), could do this only if its thermal part is concentrated across the measured region. That must be calculated, not assumed.
  - A universal interaction can produce different equilibrium dispersions in different hosts. One identical dispersion everywhere is not what universality requires.
- **Supply must be cosmological if it exists.**
  - Local starlight falls short by 10⁶–10⁹ here, as in CR-1.
  - [CC-1](../clock-completion/report.md)'s ledger stores all the energy light has lost since the turnaround in the clock field. With V = 0, that store is z_max times today's radiation energy. Matching about 26% of 3H²c²/8πG would require z_max ≈ 5,200. That is a required-budget estimate for the toy, not evidence that such a turnaround occurred or that the energy can form halos.
  - Whether that energy can clump into companions is untested. No mechanism for it exists yet, and z_max is free.
- **The stellar mass is part of the story.** With identical conventions, the lenses need 1.5–3 times Chabrier stellar mass in both geometries. That part of the inner shortfall is an IMF question, not a reservoir question.

## Reproduce

```sh
python research_work/results/supported-reservoir/cr2.py      # about 5 minutes on 6 workers; fails fast and logs progress
python research_work/results/supported-reservoir/checks.py   # suite job: V1, V2 for all six lenses, one archived fit per geometry
python research_work/results/supported-reservoir/universal-halo.py   # post-hoc universal comparison, about 1 minute on 6 workers
python research_work/results/supported-reservoir/milky-way-interior.py   # post-hoc Milky Way interior diagnostic, seconds
```
