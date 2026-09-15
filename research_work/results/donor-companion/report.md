# RC-2a report: with the donor counted and the boundary removed, cold decay at 2B-F1's rate overshoots the Milky Way; the 300 kpc zone had chosen the fit

- **The question.** With the source law and its parameters held fixed, do the masses and velocities at fixed observational radii converge as the computational region grows, or does the boundary choose them? The reference case is cold collisionless decay at 3 km/s, at 2B-F1's rate of 1,438 M☉ kpc⁻³ Gyr⁻¹. The donor is counted once, and every daughter is followed, bound or unbound.
- **The fit exists only at the old zone's size.**
  - In a 300 kpc region the Milky Way fits to 7.9 ± 0.6 km/s, with slope 1.25. That reproduces 2B-F1's static-bath result for the same case.
  - At 150 kpc too little falls in: RMSE 43 km/s.
  - From 600 kpc outward the inner rotation overshoots by about 170 km/s: RMSE 176–189 km/s over 600–2,400 kpc. The slope, 1.15–1.17, stays near the required 1.31.
- **The apertures converge, to the wrong galaxy.**
  - From 600 to 2,400 kpc the net contrast inside each aperture out to 200 kpc agrees to within about 10%: 19–21 baryon masses inside 50 kpc, and 78–79 inside 200 kpc. That is 8×10¹² M☉, against the roughly 10¹² M☉ the Milky Way's halo tracers indicate.
  - The circular speed is 410–446 km/s at every aperture out to 200 kpc.
  - Under the declared statistical rule, agreement holds from 1,200 kpc. Convergence is demonstrated for no doubling, because the seeds of these runs scatter by up to 35 km/s in RMSE.
- **Where the mass came from.**
  - In the open region, the companions inside 25 kpc at T were born between 140 and 570 kpc (10th and 99th percentiles), mostly in the first 7 Gyr.
  - So the infall horizon, measured, is 570 kpc. It is the same in the 1,200 and 2,400 kpc regions, against the 332 kpc the protocol estimated for the baryons alone.
  - In the 300 kpc region, the companions inside 25 kpc came from all the way out to its edge (99th percentile 296 kpc): the zone cut the supply off.
- **The inner profile is a transient, not a reservoir.**
  - Continuous production crosses near the Milky Way's profile at 8 Gyr (RMSE 27 km/s) and overshoots to 188 km/s by 10 Gyr.
  - Stopping production at 5 Gyr gives 20.3 km/s at 8 Gyr, then 96 km/s at 10 Gyr.
  - Converting the same total mass in the first 5 Gyr gives 21 km/s at 5 Gyr and 981 km/s at 10 Gyr.
  - No history passes G1 and G2 at any snapshot, so the declared readings are empty.
- **The controls.**
  - **The quiet representation's null (C1a) passes.** Without baryons, it puts 0.02% of the reference's mass inside 25 kpc.
  - **Random sampling (C1b) seeds structure at 2B-F1's resolution.** Without baryons, at 16,000 births, noise built up to 7.6×10¹⁰ M☉ inside 25 kpc. At 64,000 births that fell to under 6×10⁸.
  - **The reference agrees across the two resolutions (C1c) without being resolved.** The declared gate to the search is therefore closed, and stages 2 and 3 did not run.
- **Validation passes:** the ledgers to 10⁻¹³, the depletion's potential to 4×10⁻⁵, bit-identity with 2B-F1's field model, and births against depletion at every aperture.
- **The supply is unchanged:** 352 times the cosmic mean density, and 2.1 million times the microwave background's energy as a minimum stock.

Protocol: [protocol.md](protocol.md), declared in 1585c9f and amended in 2a25cb3 before execution. Code:
- [open_region.py](open_region.py): the model, 2B-F1's field channel in stage 2A's engine with the donor's depletion, an open boundary, the quiet representation and birth tags;
- [rc2a.py](rc2a.py): the driver;
- [checks.py](checks.py): the suite job;
- stage 2A's [mc.py](../companion-formation/mc.py), with optional per-tracer tags that are absent in every earlier run.

Results: [rc2a-results.json](rc2a-results.json).

## What was run

Stage 1 of Amendment 1, at the fixed rate, in 78 minutes on 4 workers:
- the reference, in regions of 150, 300, 600, 1,200 and 2,400 kpc, three seeds each, using the quiet representation;
- C1a, the quiet null without baryons, at 2,400 kpc;
- C1b, random sampling without baryons at 16,000 and 64,000 births, at 2,400 kpc;
- C1c, the reference at 64,000 births, at 1,200 kpc;
- C3 and C3b at 1,200 kpc;
- D4.

C2 is evaluated in every run. Stages 2 and 3, the search and the check at the selected rate, did not run: the declared gate needs C1a to pass and C1c to demonstrate convergence.

## The region-size comparison at the fixed rate

Seed means; masses are in baryon masses (1.004×10¹¹ M☉), with the circular speed of baryons, companions and depletion in brackets.

| R_comp (kpc) | RMSE, km/s (inner 20) | slope 8–20 kpc | net inside 50 kpc, M_b (v_circ, km/s) | net inside 100 kpc, M_b (v_circ, km/s) | net inside 150 kpc, M_b (v_circ, km/s) | net inside 200 kpc, M_b (v_circ, km/s) | net inside 300 kpc, M_b (v_circ, km/s) | escaped, M_b | tracers at T |
|---|---|---|---|---|---|---|---|---|---|
| 150 | 43.3 ± 0.8 (38.4) | 1.41 ± 0.01 | 0.6 (115) | 1.0 (94) | – | – | – | 0.0 | 4,973 |
| 300 | 7.9 ± 0.6 (4.2) | 1.25 ± 0.07 | 3.4 (193) | 7.2 (187) | 8.5 (166) | 8.2 (141) | – | 0.1 | 1.5×10⁴ |
| 600 | 189.4 ± 11.2 (176.5) | 1.15 ± 0.05 | 20.8 (433) | 45.1 (446) | 62.9 (429) | 78.2 (414) | 80.5 (343) | 0.4 | 4.4×10⁴ |
| 1200 | 187.5 ± 35.0 (174.3) | 1.17 ± 0.02 | 19.5 (419) | 40.8 (424) | 61.4 (424) | 79.4 (417) | 90.3 (362) | 9.5 | 1.3×10⁵ |
| 2400 | 175.5 ± 29.5 (162.0) | 1.17 ± 0.03 | 18.6 (411) | 40.4 (423) | 61.7 (425) | 78.3 (414) | 93.5 (369) | 112.0 | 3.6×10⁵ |

| Doubling | agrees (declared statistical rule) | demonstrated (agrees and resolved) | signal quantities that disagree | agree but unresolved |
|---|---|---|---|---|
| 150 → 300 kpc | no | no | 44: net 100 kpc, net 50 kpc, rmse_38, slope_8_20, … | 0 |
| 300 → 600 kpc | no | no | 48: net 100 kpc, net 150 kpc, net 200 kpc, net 50 kpc, … | 0 |
| 600 → 1200 kpc | no | no | 2: net 300 kpc, v_circ 300 kpc | 48 |
| 1200 → 2400 kpc | yes | no | 0 | 49 |

- **The net contrast is the inflow.** Inside every aperture the births match the donor's depletion (control C2, below). So the net contrast equals the companion mass carried in across the aperture, the identity in the owner's review.
- **Escapes.** Companions near the region's edge leave, which is why the whole region's net is negative. That does not act on the interior.

## Where the inner mass came from

Mass-weighted quantiles over the companions inside each radius at T; seed means.

| R_comp (kpc) | inside (kpc) | companion mass | birth radius, kpc: 10% / 50% / 90% / 99% | birth time, Gyr: 10% / 50% / 90% | median birth angular momentum, kpc km/s | reached 25 kpc |
|---|---|---|---|---|---|---|
| 300 | 10 | 4×10¹⁰ M☉ | 83 / 196 / 276 / 294 | 0.6 / 3.0 / 7.5 | 453 | 1.00 |
| 300 | 25 | 1.3×10¹¹ M☉ | 82 / 191 / 274 / 296 | 0.5 / 2.9 / 7.5 | 430 | 1.00 |
| 300 | 100 | 7.8×10¹¹ M☉ | 111 / 207 / 277 / 297 | 0.5 / 3.1 / 7.7 | 469 | 0.93 |
| 300 | 200 | 1.3×10¹² M☉ | 130 / 224 / 284 / 298 | 0.8 / 4.4 / 8.2 | 515 | 0.86 |
| 1200 | 10 | 3×10¹¹ M☉ | 142 / 352 / 535 / 570 | 0.6 / 3.2 / 6.9 | 796 | 1.00 |
| 1200 | 25 | 8.8×10¹¹ M☉ | 144 / 337 / 534 / 569 | 0.6 / 3.0 / 6.9 | 752 | 1.00 |
| 1200 | 100 | 4.2×10¹² M☉ | 196 / 371 / 537 / 572 | 0.6 / 3.3 / 7.1 | 851 | 0.84 |
| 1200 | 200 | 8.5×10¹² M☉ | 242 / 431 / 551 / 593 | 0.7 / 3.5 / 7.3 | 984 | 0.77 |
| 2400 | 10 | 2.8×10¹¹ M☉ | 143 / 352 / 540 / 570 | 0.5 / 3.1 / 7.1 | 774 | 1.00 |
| 2400 | 25 | 8.2×10¹¹ M☉ | 141 / 325 / 537 / 568 | 0.5 / 3.0 / 7.0 | 730 | 1.00 |
| 2400 | 100 | 4.1×10¹² M☉ | 196 / 366 / 535 / 571 | 0.6 / 3.3 / 7.1 | 843 | 0.86 |
| 2400 | 200 | 8.3×10¹² M☉ | 241 / 428 / 549 / 590 | 0.7 / 3.6 / 7.3 | 980 | 0.77 |

## The controls

**C1a, the quiet null** (no baryons, 2,400 kpc):

| Aperture (kpc) | net contrast without baryons, rms over seeds | companions with baryons | ratio |
|---|---|---|---|
| 50 | 6×10⁹ M☉ | 1.9×10¹² M☉ | 0.0032 |
| 100 | 2.8×10¹⁰ M☉ | 4.1×10¹² M☉ | 0.0068 |
| 150 | 3.6×10¹⁰ M☉ | 6.4×10¹² M☉ | 0.0057 |
| 200 | 8.7×10¹⁰ M☉ | 8.3×10¹² M☉ | 0.0104 |
| 300 | 4.5×10¹¹ M☉ | 1.1×10¹³ M☉ | 0.0412 |

Inside 25 kpc: 1.8×10⁸ M☉ without baryons, against 8.2×10¹¹ with them (ratio 0.0002). The limit is 0.1; C1a passes: yes.

**C1b, random sampling without baryons** (2,400 kpc):

| Aperture (kpc) | net contrast, rms over seeds: 16,000 births | 64,000 births | ratio |
|---|---|---|---|
| 50 | 5.1×10¹⁰ M☉ | 4.4×10⁹ M☉ | 0.09 |
| 100 | 2×10¹¹ M☉ | 5.4×10¹⁰ M☉ | 0.26 |
| 150 | 1.9×10¹¹ M☉ | 6.5×10¹⁰ M☉ | 0.34 |
| 200 | 5.5×10¹¹ M☉ | 5.6×10¹⁰ M☉ | 0.10 |
| 300 | 2.2×10¹¹ M☉ | 4.2×10¹¹ M☉ | 1.88 |

The excess inside 25 kpc over the uniform expectation, rms over seeds: 4.4×10¹⁰ M☉ at 16,000 births and 6.4×10⁸ M☉ at 64,000 (ratio 0.015). Per seed, the companions inside 25 kpc: 2.7×10⁸, 1.3×10¹⁰, 7.6×10¹⁰ M☉ at 16,000 births; 5×10⁸, 2.9×10⁸, 1.6×10⁸ M☉ at 64,000.

Random draws alone formed central concentrations comparable to the baryons at 16,000 births. 2B-F1 and stage 2A used random draws, but with the baryons as the dominant perturbation. The reference at 300 kpc, in the quiet representation, reproduces 2B-F1's result for the same case, so that result is not a sampling artifact.

**C1c, the reference at two resolutions** (1,200 kpc):

| Quantity | 16,000 births | 64,000 births | difference ± standard error | tolerance | agrees | resolved |
|---|---|---|---|---|---|---|
| RMSE, km/s | 187.5 | 173.6 | −14.0 ± 21.5 | 1 | yes | no |
| slope 8–20 kpc | 1.17 | 1.20 | +0.03 ± 0.03 | 0.05 | yes | no |

Across all 50 signal quantities, 0 disagree and 50 agree without being resolved. The largest rotation-speed difference is −18.0 ± 19.6 km/s at 5.7 kpc, against a tolerance of 2 km/s. Demonstrated: no; statistical: yes.

**C3 and C3b, the source switched off** (1,200 kpc; seed means; "pass" marks G1 and G2 both passing):

| History | 2 Gyr: RMSE (slope) | 5 Gyr: RMSE (slope) | 8 Gyr: RMSE (slope) | 10 Gyr: RMSE (slope) |
|---|---|---|---|---|
| continuous (the reference) | 51.7 (1.50) | 42.8 (1.23) | 26.7 (1.26) | 187.5 (1.17) |
| stopped at 5 Gyr (C3) | 51.8 (1.48) | 41.7 (1.29) | 20.3 (1.25) | 96.4 (1.16) |
| twice the rate until 5 Gyr (C3b) | 51.0 (1.56) | 21.1 (1.28) | 183.1 (1.18) | 980.6 (0.53) |

Declared readings: none (no history passes G1 and G2 at any snapshot).

## Validation

| Check | result | declared limit | passed |
|---|---|---|---|
| D1: companions plus escaped against births | 1.3e-14 | 10⁻⁹ | yes |
| D1: depletion against the mass converted | 8.2e-14 | 10⁻⁹ | yes |
| D1: births drawn against expectation | 1.86 standard errors | 4 | yes |
| D1: energy ledger closure | 2.6e-15 | 10⁻⁹ | yes |
| D1: tags follow the tracers | yes | every run | yes |
| D2: depletion potential against the uniform sphere | 3.7e-05 | 10⁻⁴ | yes |
| D4: bit-identical to 2B-F1's field model | state yes, ledger yes, tags yes | identical | yes |
| C2: births inside each aperture against its depletion | 2.54 standard errors | 4 | yes |

## Energy

At the fixed rate (1,438 M☉ kpc⁻³ Gyr⁻¹), the field makes 1.4×10⁴ M☉ kpc⁻³ of companions by T: 352 times the cosmic mean density, and 2.1×10⁶ times the microwave background's energy as a minimum stock (ΓT = 1). At ΓT = 10% and 1% the parent stock would be 2.1×10⁷ and 2.1×10⁸ times it. The field loses 2.8×10⁻²⁵ W m⁻³.

## Deviations from the protocol, and implementation choices

1. **Three starts.**
   - **The first canonical attempt** began at 12:48 under the protocol as declared in 1585c9f. It was stopped during its first round, before any result was read, to apply Amendment 1.
   - **The second,** under Amendment 1, was stopped a minute in. With 8 workers beside the other canonical runs, the machine's commit memory fell to 7.9 GB, and stage 2A's rerun had already lost runs to memory errors once.
   - **The run reported here** used 4 workers. Seeds are drawn in a fixed order, so the worker count does not change the results.
2. **The engine's tags.** Stage 2A's engine gained an optional per-tracer array: absent by default, and refused with bound-bound collisions, whose resampling replaces tracers. With it absent, F1's and stage 2A's short runs, collisions included, are bit-identical to the unmodified engine.
3. **Pool regeneration.** With every draw kept, the field's pools do not depend on the potential, so the engine's regeneration on a 1% change of the potential at r_half is turned off. That trigger divides by the potential, which is zero without baryons. The pools are still regenerated every 10 steps.
4. **The circular speed's tolerance** sets the net contrast's absolute tolerance: the mass that moves the aperture's circular speed by 2 km/s, using the mean speed of the two runs compared.

## What this does and does not show

**Shown.**
- **2B-F1's Milky Way fit depended on its 300 kpc zone.** With the donor counted and the boundary removed, the same cold source at the same rate delivers infall from out to about 570 kpc. The inner mass then overshoots by about a factor of seven inside 25 kpc.
- **The profile does not settle.** In every production history it passes near the observed profile on the way up and keeps growing.
- **The explicit donor-and-companion calculation agrees with 2B-F1's bath accounting** when both stop at 300 kpc (7.9 against 7.5 km/s), so the difference is the boundary, not the bookkeeping.
- **Random sampling can seed central structure** in the cold medium at 16,000 births. The quiet representation does not.

**Not shown.**
- **The rate that would fit at 1,200 kpc.** The gate closed, so the search did not run. A fitted rate would be tuned to a transient in any case.
- **Demonstrated convergence at 2 km/s.** The collapsing state's seeds scatter by up to 35 km/s. Demonstrating convergence at 2 km/s would take hundreds of seeds per region.
- **Three-dimensional dynamics.** Spherical symmetry keeps only the radial mode of the daughters' instability, which C1b shows is active.
- **A localized or finite source, the donor's own dynamics, and compact topologies.**

## What the protocol says next

- **The declared rule.** RC-2a stops after stage 1, and the representation is fixed before any search.
- **What the controls show.** They point elsewhere. The representation passes its null, and the two resolutions agree. The gate closed on precision, because the physical state's outcomes scatter between seeds.
- **The owner's table for this outcome.** "The solution converges but the inner/outer profile is wrong" leads to the birth-radius diagnostics above and a source-response feasibility calculation, to design an extended or finite source law.
- **The target that diagnostic sets.** A production history that puts about 10¹¹ M☉ inside 25 kpc by the present, and keeps it there, without the infall from 300–570 kpc that overshoots.
- **No search yet.** Whether to run stage 2 under a further amendment is the owner's decision. It is not taken here.

## Reproduce

```sh
python research_work/results/donor-companion/rc2a.py      # about 80 minutes on 4 workers (stage 1)
python research_work/results/donor-companion/checks.py    # suite job
```

The driver compares its output with the archived results and overwrites them only with `--canonical`. Wall-clock timings are excluded from the comparison.
