# CC-2 stage 2B-F1 report: cold decays pass the Milky Way's profile gate, as a compensated collapse from a field of about 300 times the cosmic mean density

- **The profile gate passes for cold decays.** At 3 and 10 km/s with the bath's gravity on, five of the six combinations pass G1, G2 and G3.
  - Rotation RMSE: 6.8–8.7 km/s over the 38 Eilers bins, and 2.3–6.8 over the inner 20.
  - Enclosed-mass slope across 8–20 kpc: 1.03–1.31, against the required 1.31. Two of the five lie within the run-to-run scatter (about ±0.1) of G2's lower edge.
  - Best rate: q\* = 1,160–1,440 M☉ kpc⁻³ Gyr⁻¹, nearly the same in every combination.
  - Stage 2A's elastic envelope grew as r^1.7–2.7, and CR-2's condensate as r^2.7. The same machinery, fed by a cold source, gives the Milky Way's shape.
- **The mechanism is cold radial infall, as the protocol expected.** At these speeds the field binds 98–100% of what it makes inside R_b. The companions fall in on nearly radial orbits (anisotropy β = 0.97–0.997 at 250 kpc), so their density goes as 1/(r²·v_esc).
- **G3 passes on the net, not on the confined mass.**
  - The confined companions are 12.6–16.2 baryon masses, close to the protocol's estimate of 16. That is above G3's 10 and inside the lenient 20.
  - A cold incident bath barely enters the zone. Inside R_b it falls 10.5–15.0 baryon masses short of the uniform density, which is assigned to the background. That shortfall brings the net cost to 0.75–2.6.
  - In effect the zone has collapsed inward and left its outer part depleted: a compensated structure. With the bath's gravity omitted, no cold combination passes.
- **Beyond the fitted range** the net extra mass reaches 7.0–9.4 baryon masses near 124 kpc, a total of 0.8–1.0×10¹² M☉, and the circular speed stays at 180–212 km/s out to about 90 kpc. Toward the zone's edge the depletion takes over.
- **Faster decays fail.**
  - At 30 km/s the slope is 2.3–2.6; at 100 km/s the cost is 23–56.
  - At 300 km/s one combination passes, with the bath's gravity omitted, just below the runaway. It needs q\* = 4.6×10⁶ (a million times the cosmic mean), while the omitted bath excess inside r_half equals the baryons there. Every 300 km/s run with the bath's gravity on ran away or lost its static bath.
- **The supply is enormous, and it is not radiation.** The passing combinations need the field to make 280–360 times the cosmic mean density in companions by T. That takes a field energy at least 1.7–2.1 million times the microwave background's, turned into companion rest mass: about two million times today's radiation, in a local stock comparison like RC-1's (the flow budget is RC-2).
- **The medium it implies is Jeans-unstable** above 12–46 kpc, with growth times of 1.1–1.2 Gyr, about a tenth of the span. The smooth bath the model assumes is therefore not self-consistent for these decays.
- **It is not universal as modeled.**
  - With the Milky Way's best source unchanged, Coma keeps no static bath in either bracket.
  - J1630 completes, with companions only 0.15–0.25 of the baryons at 3–21 kpc.
  - With a 150 kpc zone, no rate passes the Milky Way's gate.

Protocol: [protocol.md](protocol.md), declared in 2cf7ea1 before execution (baseline 6a3b14b). Code:
- [field.py](field.py): the decay kinematics, the field as a source channel of stage 2A's engine, and the validation quadratures;
- [f1.py](f1.py): the driver;
- [checks.py](checks.py): the suite job;
- stage 2A's [mc.py](../companion-formation/mc.py), with two added hooks that do nothing by default.

Results: [f1-results.json](f1-results.json).

**Why it was run.** It follows the owner's review of 68eb17c. Supply stopped being the principal obstacle; where the mass goes became it.
- **Stage 2A's result.** Its elastic capture mechanism builds an outer envelope: M(<r) ∝ r^1.7–2.7, while the Milky Way needs about r^1.3.
- **Why more supply won't do.** A larger source with the same phase space would only build a heavier envelope.
- **So stage 2B tests a phase-space distribution as well as an energy supply.** The owner's first toy was a coherent field decaying into companion pairs, which also meets RC-1's momentum objection.

## The source

- **The decay.** A homogeneous coherent field decays into companion pairs, χ → C + C, at a constant rate q per unit volume. Its lifetime is taken to be much longer than the 10 Gyr span.
  - Each companion leaves with the decay speed v_d = c√(2ε) in an isotropic direction. The pair carries no momentum.
  - Per unit of field energy, companion rest mass takes 1/(1+ε) and kinetic energy ε/(1+ε). The field loses exactly what it gives.
- **Trial parameters.** Decay speeds of 3, 10, 30, 100 and 300 km/s (ε = 5×10⁻¹¹ to 5×10⁻⁷), each the same in every system, and σ/m = 0, 0.1 and 1 cm²/g under CF-1's elastic law.

## How the population is built

- **Stage 2A's machinery:**
  - exact orbits in the evolving spherical potential;
  - all three collision classes;
  - self-gravity;
  - the bath's focused excess, gravitating or omitted.
- **One added source channel.** Pairs are born throughout each zone at rate q, and a born companion that is confined becomes a tracer.
- **A growing incident bath.** Companions produced outside arrive at R_b with speed v_d, at a density growing as qt.
- **The zone of influence.** R_b is 300 kpc for the Milky Way and J1630 and 30 Mpc for Coma. The Milky Way is repeated with a 150 kpc zone.

## Validation (declared tolerances)

| Check | Declared tolerance | Result | Passed |
|---|---|---|---|
| F1, decay kinematics | rounding; isotropy within 3 standard errors | energy 0, momentum 0, mass shell 2.2×10⁻¹⁶; largest isotropy deviation 1.9 standard errors | yes |
| F2, born-bound mass (edge-refined quadrature) | 3 standard errors | 100 km/s: −0.6; 300 km/s: +0.2 | yes |
| F3, cold-birth density at 5–25 kpc | 3 standard errors or 5% | worst bin 0.39 of its tolerance | yes |
| F4, seedless production with a growing bath (prediction from the pools the run used) | 3 standard errors | +1.1 standard errors | yes |
| V5, collisions | rounding | kinetic 1.2×10⁻¹⁵, momentum 3.1×10⁻¹⁶ | yes |
| V6, ledgers | energy 10⁻⁹; mass exact | energy 3.5×10⁻¹⁴, mass 8.4×10⁻¹⁵ | yes |

- **F2 and F4 are the corrected checks.** The first canonical run failed both, for reasons in the checks and not in the engine (deviation 9). The corrections draw no random numbers, and every ladder run and first verification came out identical in all three canonical runs (deviation 10 covers the refinements).
- **Stage 2A's engine checks carry over.** Its suite job reruns V2, V3, V5, V6 and V7 with the hooks in place. V5 and V6 are recomputed here across every run.
- **The passing collisional combinations are transparent:** τ(r_half) ≤ 0.011, against the regime's limit of 0.3.

## The Milky Way's profile gate

The table gives, for each combination, the gates at the best-fitting rate q\*. It lists:
- the rotation RMSE over the 38 Eilers bins, and over the inner 20;
- the companions' enclosed-mass slope across 8–20 kpc (the required slope is 1.31);
- the cost: all gravitating non-baryonic mass inside R_b, over the baryons.

G1 is CR-2's fit rule, RMSE ≤ 20 km/s without worsening the inner 20. G2 is the slope within 0.3 of 1.31. G3 is a cost of at most 10, with 20 as a lenient variant.

| Decay speed (km/s) | σ/m (cm²/g) | Bath gravity | q\* (M☉ kpc⁻³ Gyr⁻¹) | RMSE, 38 bins (inner 20), km/s | slope 8–20 kpc | cost: non-baryonic mass in R_b over baryons | G1 | G2 | G3 | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 0 | on | 1,438 (bracketed) | 6.8 (2.4) | 1.3 | 1.1 (16 confined, −15 bath) | yes | yes | yes | **yes** |
| 3 | 0 | omitted | 1,361 (bracketed) | 8.0 (2.2) | 1.3 | 15 (15 confined) | yes | yes | no | **no** |
| 3 | 0.1 | on | 1,200 (bracketed) | 7.5 (5.5) | 1.1 | 0.75 (13 confined, −12 bath) | yes | yes | yes | **yes** |
| 3 | 0.1 | omitted | 1,102 (bracketed) | 8.2 (7.7) | 0.95 | 12 (12 confined) | yes | no | no | **no** |
| 3 | 1 | on | 685 (bracketed) | 10.9 (8.1) | 0.75 | 0.17 (7.2 confined, −7 bath) | yes | no | yes | **no** |
| 3 | 1 | omitted | 696 (bracketed) | 10.3 (10.8) | 0.74 | 7.3 (7.3 confined) | yes | no | yes | **no** |
| 10 | 0 | on | 1,437 (bracketed) | 8.7 (6.8) | 1.3 | 2.6 (16 confined, −13 bath) | yes | yes | yes | **yes** |
| 10 | 0 | omitted | 1,373 (bracketed) | 7.7 (4.1) | 1.3 | 15 (15 confined) | yes | yes | no | **no** |
| 10 | 0.1 | on | 1,412 (bracketed) | 6.9 (2.3) | 1.2 | 2.5 (16 confined, −13 bath) | yes | yes | yes | **yes** |
| 10 | 0.1 | omitted | 1,367 (bracketed) | 7.5 (4.3) | 1.2 | 15 (15 confined) | yes | yes | no | **no** |
| 10 | 1 | on | 1,158 (bracketed) | 6.8 (4.0) | 1 | 2 (13 confined, −11 bath) | yes | yes | yes | **yes** |
| 10 | 1 | omitted | 1,128 (bracketed) | 7.0 (3.5) | 1.1 | 12 (12 confined) | yes | yes | no | **no** |
| 30 | 0 | on | 1,885 (bracketed) | 18.4 (20.0) | 2.4 | 6.4 (20 confined, −13 bath) | yes | no | yes | **no** |
| 30 | 0 | omitted | 2,009 (bracketed) | 17.4 (19.9) | 2.3 | 21 (21 confined) | yes | no | no | **no** |
| 30 | 0.1 | on | 1,925 (bracketed) | 19.7 (20.3) | 2.6 | 6.2 (20 confined, −14 bath) | yes | no | yes | **no** |
| 30 | 0.1 | omitted | 1,854 (bracketed) | 19.6 (23.0) | 2.5 | 20 (20 confined) | yes | no | no | **no** |
| 30 | 1 | on | 1,811 (bracketed) | 17.6 (19.0) | 2.5 | 6 (19 confined, −13 bath) | yes | no | yes | **no** |
| 30 | 1 | omitted | 1,809 (bracketed) | 17.7 (20.5) | 2.4 | 19 (19 confined) | yes | no | no | **no** |
| 100 | 0 | on | 7,435 (bracketed) | 22.8 (24.7) | 2.4 | 25 (56 confined, −30 bath) | no | no | no | **no** |
| 100 | 0 | omitted | 1.1×10⁴ (bracketed) | 26.2 (31.7) | 3.4 | 56 (56 confined) | no | no | no | **no** |
| 100 | 0.1 | on | 6,685 (bracketed) | 27.7 (30.7) | 3.1 | 23 (48 confined, −25 bath) | no | no | no | **no** |
| 100 | 0.1 | omitted | 9,888 (bracketed) | 28.0 (31.9) | 2.7 | 44 (44 confined) | no | no | no | **no** |
| 100 | 1 | on | 6,562 (bracketed) | 24.4 (12.8) | 2.2 | 23 (48 confined, −25 bath) | no | no | no | **no** |
| 100 | 1 | omitted | 10⁴ (bracketed) | 19.1 (20.2) | 1.9 | 50 (50 confined) | yes | no | no | **no** |
| 300 | 0 | on | – (none completed) | no verification | – | – | – | – | – | – |
| 300 | 0 | omitted | 4.6×10⁶ (high edge) | 8.2 (3.7) | 1.4 | 5.4 (5.4 confined) | yes | yes | yes | **yes** |
| 300 | 0.1 | on | – (none completed) | no verification | – | – | – | – | – | – |
| 300 | 0.1 | omitted | – (none completed) | no verification | – | – | – | – | – | – |
| 300 | 1 | on | – (none completed) | no verification | – | – | – | – | – | – |
| 300 | 1 | omitted | – (none completed) | no verification | – | – | – | – | – | – |

The extra mass the model supplies at five radii, over what the observed speeds require:

| Decay speed | σ/m | Bath gravity | 6 kpc | 10 kpc | 15 kpc | 20 kpc | 24 kpc |
|---|---|---|---|---|---|---|---|
| 3 | 0 | on | 1.1 | 1 | 0.95 | 1.1 | 1.4 |
| 3 | 0 | omitted | 1.1 | 1 | 1 | 1.2 | 1.5 |
| 3 | 0.1 | on | 1.2 | 0.99 | 0.87 | 0.99 | 1.3 |
| 3 | 0.1 | omitted | 1.4 | 1.1 | 0.95 | 1 | 1.3 |
| 3 | 1 | on | 1.4 | 1 | 0.79 | 0.81 | 0.96 |
| 3 | 1 | omitted | 1.6 | 1.2 | 0.91 | 0.91 | 1.1 |
| 10 | 0 | on | 0.87 | 0.84 | 0.8 | 0.98 | 1.3 |
| 10 | 0 | omitted | 0.9 | 0.9 | 0.88 | 1.1 | 1.4 |
| 10 | 0.1 | on | 1 | 0.99 | 0.91 | 1.1 | 1.4 |
| 10 | 0.1 | omitted | 0.88 | 0.94 | 0.85 | 1 | 1.3 |
| 10 | 1 | on | 1.1 | 1 | 0.89 | 0.99 | 1.2 |
| 10 | 1 | omitted | 1.1 | 1.1 | 0.95 | 1.1 | 1.4 |
| 30 | 0 | on | 0.28 | 0.51 | 0.77 | 1.3 | 1.8 |
| 30 | 0 | omitted | 0.23 | 0.53 | 0.78 | 1.2 | 1.6 |
| 30 | 0.1 | on | 0.25 | 0.49 | 0.86 | 1.4 | 1.9 |
| 30 | 0.1 | omitted | 0.22 | 0.43 | 0.69 | 1.1 | 1.7 |
| 30 | 1 | on | 0.24 | 0.53 | 0.88 | 1.3 | 1.8 |
| 30 | 1 | omitted | 0.29 | 0.49 | 0.75 | 1.1 | 1.6 |
| 100 | 0 | on | 0.1 | 0.53 | 0.56 | 1.2 | 2.1 |
| 100 | 0 | omitted | 0.22 | 0.26 | 0.57 | 1.3 | 1.9 |
| 100 | 0.1 | on | 0.062 | 0.31 | 0.47 | 1.3 | 2.2 |
| 100 | 0.1 | omitted | 0.22 | 0.21 | 0.39 | 0.77 | 1.6 |
| 100 | 1 | on | 0.47 | 0.69 | 0.91 | 1.7 | 2.6 |
| 100 | 1 | omitted | 0.094 | 0.63 | 0.76 | 1.2 | 1.9 |
| 300 | 0 | omitted | 0.82 | 0.93 | 0.95 | 1.1 | 1.5 |

Beyond the fitted range, the net extra mass over the baryons (the confined companions plus, where it gravitates, the bath's excess), with the circular speed from all the mass in brackets (km/s). The profile-gate passes are listed, with the cold collisionless runs without bath gravity for contrast:

| Decay speed | σ/m | Bath gravity | 43 kpc | 87 kpc | 124 kpc | 177 kpc | 252 kpc | profile gate |
|---|---|---|---|---|---|---|---|---|
| 3 | 0 | on | 3.5 (210) | 7.3 (202) | 9 (186) | 8.4 (152) | 1.1 (61) | yes |
| 3 | 0 | omitted | 3.6 (213) | 7.8 (208) | 11 (202) | 13 (188) | 15 (167) | no |
| 3 | 0.1 | on | 2.8 (194) | 5.8 (183) | 7 (167) | 6.7 (137) | 0.75 (55) | yes |
| 10 | 0 | on | 3.1 (201) | 7.5 (204) | 9.4 (190) | 8.9 (156) | 2.6 (79) | yes |
| 10 | 0 | omitted | 3.4 (208) | 7.5 (204) | 11 (201) | 14 (191) | 15 (167) | no |
| 10 | 0.1 | on | 3.2 (205) | 7.1 (200) | 8.8 (185) | 8.4 (152) | 2.5 (78) | yes |
| 10 | 1 | on | 2.5 (187) | 5.7 (181) | 7.1 (167) | 7 (140) | 2 (72) | yes |
| 300 | 0 | omitted | 3.5 (211) | 5 (172) | 5.3 (148) | 5.4 (125) | 5.4 (105) | yes |

**Why cold decays work.**
- **Nearly every birth is bound.** A companion born at 3–10 km/s is bound nearly wherever it is born inside 300 kpc. The field binds 98–100% of what it makes in the zone, and 96–100% of the mass made is still confined at T.
- **The radial-infall profile.** Each companion falls from almost rest on a nearly radial orbit and spends time at radius r in proportion to 1/v_esc(r). The time-averaged density is ρ ∝ 1/(r²·v_esc(r)), so M(<r) ∝ r/v_esc(r). The protocol expected this before running: about r^1.4 in the baryonic potential, and a confined mass of about 16 baryon masses at R_b = 300 kpc.
- **Against the requirement.** In all five cold passes the extra mass is 0.80–1.24 of what is required from 6 to 20 kpc, and 1.24–1.40 at 24 kpc.
- **Why the rate barely moves between combinations.** The mass inside 25 kpc is set by how much of the zone has fallen in by T, not by the interaction. The confined mass is still growing in proportion to time (d ln M/d ln t = 0.86–1.01 at T).
- **Scatter.** Across the three canonical runs, seven cold refinement runs were repeated at the same rate with a different seed (deviation 10). Their slopes moved by 0.01–0.19 and their RMSE by up to 2.5 km/s, so a slope carries about ±0.1.
  - The collisionless and 0.1 cm²/g slopes differ by 0.12–0.21, comparable to that scatter.
  - At 1 cm²/g, collisions flatten the 3 km/s profile to r^0.72–0.75 in every run, which fails G2. At 10 km/s the slope is 1.03.
  - Two of the passes lie within the scatter of G2's lower edge, 1.01: 3 km/s at 0.1 cm²/g (1.05) and 10 km/s at 1 cm²/g (1.03).

**What G3 counts.** The protocol's cost is all gravitating non-baryonic mass inside R_b: the confined companions plus, where it gravitates, the bath's focused excess. The uniform part is assigned to the background, a labeled assumption.
- **For a cold bath that excess is negative.** Only orbits that reach R_b are populated, and at 3–10 km/s they are nearly radial, so the incident bath fills little of the zone. Inside R_b it falls short of the uniform density by 10.5–15.0 baryon masses, nearly the mass the field made there and bound.
- **So G3 passes on a compensated structure.** The companions born in the zone have in effect collapsed toward the galaxy and left the zone's outer part depleted. The net gravitating mass inside R_b is therefore small.
- **Without that accounting nothing passes.** Counted alone, the confined companions (12.6–16.2 baryon masses) fail G3 and pass the lenient 20. With the bath's gravity omitted, no cold combination passes.
- **The verdict rests on a labeled assumption:** that the uniform density belongs to the background. RC-2 and the global background must settle that.

**The outer profile.**
- **Inside about 125 kpc.** Out to about 90 kpc the circular speed stays at 180–212 km/s. The total mass reaches 0.8–1.0×10¹² M☉ near 124 kpc. That is the order of the Milky Way's mass that its halo stars and satellites indicate, from which G3's benchmark was declared.
- **Beyond about 150 kpc** the depleted zone dominates; by 252 kpc the circular speed has fallen to 55–79 km/s. That fall comes from the zone's boundary model and from the unmodeled infall beyond R_b, so it is not a prediction.

**Faster decays.**
- **30 km/s.** Pericenters keep companions out of the inner kiloparsecs. The extra mass is 0.22–0.29 of the requirement at 6 kpc and 1.6–1.9 times it at 24 kpc (slope 2.3–2.6). The fits still meet G1 (17.4–19.7 km/s), which is why G2 was declared.
- **100 km/s.** The slope is 1.9–3.4 and the cost 23–56.
- **300 km/s.** Almost nothing is born bound in the baryons' potential. The population that does form deepens the potential and binds more births.
  - **Bath gravity on.** Every ladder run ran away (1,000 baryon masses within 0.03–0.59 Gyr) or lost its static bath.
  - **Bath gravity omitted, no collisions.** Rates up to 4.3×10⁶ complete, and 4.9×10⁶ runs away at 9.5 Gyr. The verification at 4.6×10⁶ passes the gate just below that runaway, but only because the bath's gravity is omitted; its excess inside r_half equals the baryons there.

## Energy supply

| Decay speed | σ/m | Bath gravity | companion mass made by T (M☉/kpc³) | over the cosmic mean | field energy over the microwave background | field energy loss (W/m³) | Jeans length (kpc) | Jeans growth time (Gyr) |
|---|---|---|---|---|---|---|---|---|
| 3 | 0 | on | 1.4×10⁴ | 352 | 2.1×10⁶ | 2.8×10⁻²⁵ | 12 | 1.1 |
| 3 | 0 | omitted | 1.4×10⁴ | 334 | 2×10⁶ | 2.6×10⁻²⁵ | 13 | 1.1 |
| 3 | 0.1 | on | 1.2×10⁴ | 294 | 1.7×10⁶ | 2.3×10⁻²⁵ | 14 | 1.2 |
| 3 | 0.1 | omitted | 1.1×10⁴ | 270 | 1.6×10⁶ | 2.1×10⁻²⁵ | 14 | 1.3 |
| 3 | 1 | on | 6,847 | 168 | 10⁶ | 1.3×10⁻²⁵ | 18 | 1.6 |
| 3 | 1 | omitted | 6,963 | 171 | 10⁶ | 1.3×10⁻²⁵ | 18 | 1.6 |
| 10 | 0 | on | 1.4×10⁴ | 352 | 2.1×10⁶ | 2.8×10⁻²⁵ | 41 | 1.1 |
| 10 | 0 | omitted | 1.4×10⁴ | 337 | 2×10⁶ | 2.6×10⁻²⁵ | 42 | 1.1 |
| 10 | 0.1 | on | 1.4×10⁴ | 346 | 2.1×10⁶ | 2.7×10⁻²⁵ | 42 | 1.1 |
| 10 | 0.1 | omitted | 1.4×10⁴ | 335 | 2×10⁶ | 2.6×10⁻²⁵ | 42 | 1.1 |
| 10 | 1 | on | 1.2×10⁴ | 284 | 1.7×10⁶ | 2.2×10⁻²⁵ | 46 | 1.2 |
| 10 | 1 | omitted | 1.1×10⁴ | 276 | 1.6×10⁶ | 2.2×10⁻²⁵ | 46 | 1.3 |
| 30 | 0 | on | 1.9×10⁴ | 462 | 2.7×10⁶ | 3.6×10⁻²⁵ | 108 | 0.97 |
| 30 | 0 | omitted | 2×10⁴ | 493 | 2.9×10⁶ | 3.9×10⁻²⁵ | 104 | 0.94 |
| 30 | 0.1 | on | 1.9×10⁴ | 472 | 2.8×10⁶ | 3.7×10⁻²⁵ | 107 | 0.96 |
| 30 | 0.1 | omitted | 1.9×10⁴ | 454 | 2.7×10⁶ | 3.6×10⁻²⁵ | 109 | 0.98 |
| 30 | 1 | on | 1.8×10⁴ | 444 | 2.6×10⁶ | 3.5×10⁻²⁵ | 110 | 0.99 |
| 30 | 1 | omitted | 1.8×10⁴ | 443 | 2.6×10⁶ | 3.5×10⁻²⁵ | 110 | 0.99 |
| 100 | 0 | on | 7.4×10⁴ | 1,823 | 1.1×10⁷ | 1.4×10⁻²⁴ | 181 | 0.49 |
| 100 | 0 | omitted | 1.1×10⁵ | 2,668 | 1.6×10⁷ | 2.1×10⁻²⁴ | 150 | 0.4 |
| 100 | 0.1 | on | 6.7×10⁴ | 1,639 | 9.7×10⁶ | 1.3×10⁻²⁴ | 191 | 0.51 |
| 100 | 0.1 | omitted | 9.9×10⁴ | 2,424 | 1.4×10⁷ | 1.9×10⁻²⁴ | 157 | 0.42 |
| 100 | 1 | on | 6.6×10⁴ | 1,608 | 9.6×10⁶ | 1.3×10⁻²⁴ | 193 | 0.52 |
| 100 | 1 | omitted | 10⁵ | 2,513 | 1.5×10⁷ | 2×10⁻²⁴ | 154 | 0.42 |
| 300 | 0 | omitted | 4.6×10⁷ | 1.1×10⁶ | 6.7×10⁹ | 8.9×10⁻²² | 22 | 0.02 |

**What the field must be.**
- **Nearly all of its energy becomes companion rest mass.** Per unit of field energy the companions' kinetic share is ε = 5×10⁻¹¹ at 3 km/s and 5.6×10⁻¹⁰ at 10 km/s.
- **How much the gate needs.** To pass it the field must deliver 1.2–1.4×10⁴ M☉ kpc⁻³ of companions by T, 280–360 times the cosmic mean density (a comparison unit only). Its own energy density is at least 1.7–2.1 million times the microwave background's, larger still since ΓT ≪ 1 is assumed, and it loses 2.2–2.8×10⁻²⁵ W/m³.
- **It is not photon energy.** Today's radiation holds about two million times less energy than the field must, the same kind of local stock shortfall RC-1 found for stage 2A's bath (10⁵–10⁸). Converting today's light cannot supply it; a flow budget is RC-2's.
  - In this toy the companions come from a separate dark sector, and the programme's photon–companion link plays no part in them.
  - The unification problem is sharper, not smaller.
- **The medium is unstable.** A uniform medium of density qT at 3–10 km/s is Jeans-unstable above 12–46 kpc and grows in 1.1–1.2 Gyr.
  - So the model's smooth incident bath and uniform background are not self-consistent for the decays that pass: the field's own companions would clump on galaxy scales long before T.
  - What that collapse does, and whether a static background can hold such a field, are questions for RC-2 and the global background.

## Universality and the zone

The best-scoring source, run unchanged, in J1630 and Coma (the same field, so the same q per unit volume):

| System | status | total companion mass over baryons | companions over baryons at the test radii (enclosed; projected) |
|---|---|---|---|
| J1630 | completed | 3.3 | 3 kpc 0.15 (0.17); 10 kpc 0.18 (0.25); 21 kpc 0.25 (0.39) |
| Coma low | no static bath | – | – |
| Coma high | no static bath | – | – |

With a 150 kpc zone the best rate gives RMSE 30.0 km/s (inner 20: 35.9), slope 0.75 and cost −0.2: G1 no, G2 no, G3 yes, gate no (best point: high edge).

- **Coma.** As in stage 2A, no static bath exists in either bracket: with the bath's gravity on, the net enclosed mass is not positive. A uniform density qT over the 30 Mpc zone is about 1.4×10¹⁸ M☉, and a cold bath that barely enters leaves most of it as a shortfall.
- **J1630.** Its companions total 3.3 times its baryons, but inside 21 kpc they are 0.15–0.25 of the baryons (0.17–0.39 in projection): most of them lie farther out. These are predictions, not gates.
- **The zone.** It sets the cost, since the field binds nearly everything it makes there, and here it also sets the profile.
  - **With R_b = 150 kpc, no rate passes.** The best-fitting completed rung (3.2×10³) fits to 17.8 km/s, with slope 0.98. The verification at 4.8×10³ gives 30.0 km/s and slope 0.75, with a net shortfall inside the zone (cost −0.2). Rates of 9.6×10³ and above lose the static bath.
  - **The protocol had expected** a 150 kpc zone to halve the cost.
  - **The zone is not derived.** It comes from CF-1's boundary model, and infall from beyond it is not modeled. A result this sensitive to it needs the global background to set it.

## Deviations from the protocol, and implementation choices

1. **Engine hooks.** Stage 2A's engine (`mc.py`) gained two hooks that do nothing by default:
   - births from an added source, inside each step's energy-closure window;
   - an update of the incident density before each step's bath tables.

   With them in place, stage 2A's suite job reproduces its archived regression anchor exactly (5,303,548.224133587 M☉). `formation.summarize` counts the field's births in its mass ledger (zero in stage 2A), and leaves the bath's self-collision time undefined at σ/m = 0.
2. **The field channel.**
   - **Pools.** Births are drawn from pools, one per logarithmic shell (40 shells): positions uniform in volume, speed v_d in isotropic directions, kept when confined. Pool sizes adapt as in stage 2A, aiming at 100 confined samples in 1,000–20,000 draws. The pools are redrawn with stage 2A's every ten steps.
   - **Tracer masses.** Masses per shell follow stage 2A's square-root rule within a factor of 10, from the confined production rate at the start of the run. A shell with no confined birth yet takes its nearest producing neighbour's mass.
   - **Pairs.** The two members of a pair are drawn as independent isotropic companions. The pair correlation does not affect a collisionless population, and a partner born at the same point with the opposite velocity has the same energy and |J|.
3. **The growing bath.** The incident density is q times the middle of each step (qΔ/2 in the first). The seedless channel's tracer masses come from the end-of-span density qT, so its births total about a third of the target over the span.
4. **Regime diagnostics.** Stage 2A's summary assumes a fixed 300 km/s bath of constant density. Here the Jeans length, Jeans growth time and capture efficiency are recomputed for a bath at speed v_d growing as qt, with an inflow over the span of πR_b²v_d·qT²/2.
5. **Starting rates.** Each ladder starts at the rate whose born-bound mass alone, in the frozen baryons-only potential, reaches ten baryon masses by T (F2's quadrature on its plain grid). The ladder spans 1/27 to 3 times that rate.
6. **Early stops, added after the first smoke run.**
   - **Why.** Hot decays (300 km/s) at the higher starting rates ran away: the population born bound deepened the potential, which bound ever more births. They reached 3,000–100,000 baryon masses and took up to 15 minutes per smoke run. Collisional runs with hot decays were far into the opaque regime.
   - **The stops.** A run stops as "not completed: runaway" once its confined mass exceeds 1,000 baryon masses, 50 times the lenient cost gate, so no gate outcome can change. A collisional run stops as "not completed: opaque" once the bath's optical depth to r_half exceeds 3, ten times the regime's limit of 0.3.
   - Stopped runs are excluded from the normalization.
7. **The smoke span** was 2 Gyr, too short for companions born in the outer zone to fall in, so smoke-size fits were not representative. The canonical span is 10 Gyr.
8. **How the normalization is found.** The protocol said only that a ladder brackets the RMSE-minimizing rate and a verification run sits at the interpolated q.
   - **The fitted quantity.** Each completed run gives s(q): the factor by which its extra mass, shape held fixed, would have to be scaled to minimize the RMSE, the way CR-2 fits a profile's mass.
   - **The rate** is the root of ln s = 0, interpolated in log q between the runs that bracket it, as stage 2A interpolates its benchmarks. Ladders extend until the root is bracketed: up to three times, by at most a factor of 3, and never past a stopped run.
   - **Refinements.** After the verification run, up to three refinements repeat the interpolation with every run so far. They keep the verification whose own s is closest to 1 and stop once it is within 10%.
   - **What failed first.** A parabola in RMSE against log q was tried in the smoke runs and failed. The RMSE is V-shaped near its minimum, and self-gravity makes the extra mass grow faster than q; one refinement jumped from 29 to 460 km/s.
9. **Validation corrections after the first canonical run.** The first canonical run (33 minutes) failed F2 and F4. Both failures were in the checks, not in the engine. The corrected checks draw no random numbers, and in the rerun every ladder run and first verification came out identical.
   - **F2** gave z = −1.01 at 100 km/s and −3.80 at 300 km/s.
     - **The cause** was the quadrature's fixed grid of 600 log-radius points, which misplaced the confinement edge. At 300 km/s the confined fraction falls from 1 to 0 within a sliver of radius (angular momentum barely moves the escape energy). There the plain trapezoid was 2.6% high, and 0.25% at 100 km/s.
     - **The fix** is an edge-refined quadrature, which splits every interval whose ends differ 64 ways. It agrees with uniform grids of 2,400 and 9,600 points and gives z = −0.61 and +0.18. The engine's own expectation from its pools agrees with it to 0.03% and 0.02%.
     - The ladders' starting rates keep the plain grid, since they need only a rough value, so the ladders are unchanged.
   - **F4** gave z = −6.56.
     - **The cause** was a prediction from one set of pools drawn after the run, while the engine redraws its pools every ten steps (201 sets over the span). One set's sampling error, 3.9% (the scatter of 40 fresh sets), was not in the standard error. The check's set lay 4.7% above the run; with that error included, the first check's z would have been −1.2.
     - **The standard error** born/√births also ignored unequal tracer masses and events confining both partners: 5.2×10³ against 5.7×10³.
     - **The corrected check** integrates the production of the pools the run used, each weighted by the exact integral of q²t² over the steps it served, and takes the births' variance from the same pools: z = +1.13. The engine's midpoint density matches the exact t² weighting to 6×10⁻⁸. The one-set value is kept in the results for the record.
10. **Seeds in the refinement rounds.** The rerun still differed from the first canonical run in two refinement runs and in one run that followed from them, all at 10 km/s without collisions.
    - **The cause.** The driver drew each refinement round's seeds in the order the previous round's runs had finished, so two runs could exchange seeds between executions.
    - **The fix.** The driver now restores the order in which the tasks were issued before anything reads the results. Every seed is drawn while iterating over issued tasks or a fixed list of combinations, so every seed is now drawn in a fixed order.
    - **The check.** Two concurrent smoke runs, whose runs finished in different orders in rounds 1 and 4, gave bit-identical results. The canonical run was then repeated a third time.
    - **The third run.** Its ladders, first verifications, universality runs and 150 kpc runs are identical to both earlier runs. Its refinements drew their seeds in the fixed order, so they differ from both. The results reported here are from this run.
    - **A direct measure of the scatter.** Across the three runs, 14 refinement runs were repeated at the same rate with a different seed.
      - **Cold decays.** The seven cold-decay runs moved by up to 0.19 in slope and 2.5 km/s in RMSE. For example, 10 km/s, collisionless, bath gravity on, at q = 1.39×10³ gave slope 1.10 and RMSE 9.6 km/s with one seed, and 1.29 and 8.6 km/s with another.
      - **The largest change** was at 100 km/s, 1 cm²/g, bath gravity on: RMSE 17.4 km/s and slope 1.56 with one seed, and 24.4 and 2.16 with another.

## What this does and does not show

**Shown.**
- **The phase space was the problem.** Stage 2A's machinery and CF-1's interaction, fed by a cold source, give the Milky Way's required enclosed-mass shape across 5–25 kpc. The elastic envelope came from the source's phase space, as the owner's review argued.
- **Cold radial infall is the mechanism.** The normalization depends on the zone and the span, not on σ/m below 1 cm²/g.
- **The cost gate passes only on the net,** through the cold bath's shortfall against the background. The confined companions alone are 12.6–16.2 baryon masses.
- **Faster decays fail the gate,** apart from one 300 km/s pass that depends on ignoring the bath's gravity at a million times the cosmic mean.
- **The field needed is far beyond today's radiation:** 280–360 times the cosmic mean density, and 1.7–2.1 million times the microwave background's energy (a local stock comparison; the flow budget is RC-2). Its uniform medium is Jeans-unstable within about a gigayear.
- **The passing configuration is not yet universal.** Coma keeps no static bath, and a 150 kpc zone fails.

**Not shown.**
- **Any fit to lensing or clusters.**
- **A global energy or flux budget.** That is RC-2.
- **The collapse of the uniform companion medium.** A cold medium at these densities is Jeans-unstable within about a gigayear.
- **Infall from beyond R_b, and so the outer halo.** Beyond about 150 kpc the profile reflects the zone's boundary model.
- **More than one field toy.**
- **A microscopic model of χ.** The small mass defect is stated, not explained.

## Next

- **RC-2 comes next, in the owner's order.** For this source it must say three things:
  - whether a static background can hold a field of at least 280–360 times the cosmic mean density;
  - what the companions do across the whole volume, including the uniform medium's collapse;
  - what sets the zone of influence that the profile depends on.
- **The pre-declared follow-ups** were declared for a failed gate. The gate passed, but two of them answer the weaknesses it exposed, each with its own protocol:
  - **a field sourced by the baryons.** Production follows the baryons, so there is no uniform medium to deplete or collapse, and the rate differs between systems, which bears on Coma and J1630;
  - **Bose-enhanced decay.** Production concentrates where companions already are, which also avoids a uniform medium.
  - A bulk flow through the field frame addresses neither and can wait.
- **The optimistic reading.** Cold phase space solves the shape problem that defeated stage 2A, with the interaction and machinery unchanged. What remains is a cold source that does not fill the universe uniformly. The next toy should keep the cold infall and localize the production.
- **Unification stays explicit.** The field is a separate dark sector. Nothing here connects its decay to redshift, and its energy is not photon energy.
- **Then the recording-transition toy,** as queued.

## Reproduce

```sh
python research_work/results/companion-source/f1.py       # about 34 minutes on 16 workers
python research_work/results/companion-source/checks.py   # suite job
```

The driver compares its output with the archived results and overwrites them only with `--canonical`. Wall-clock timings are excluded from the comparison.
