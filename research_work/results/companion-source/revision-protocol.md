# CC-2 stage 2B-F1 revision: a direct normalization search, seed repeats at higher resolution, and reproducible drivers

Declared before execution, 15 September 2026. Baseline: `main` at e3eb14a.

The source, interaction, systems, zones of influence, span, gates and validation tolerances are those of the 2B-F1 protocol (2cf7ea1). Only what is stated here changes.

**Why it is run.** The owner's review of ed96b00 asks for four things before any wider campaign.
- **The search.** The rescaling root is not in general the rate that minimizes the actual RMSE, because q changes the profile's shape as well as its normalization. Each combination's best actual sampled RMSE is to be kept, and its neighbourhood bracketed and evaluated directly.
- **Marginal passes.** Selected cases are to be repeated at fixed parameters, with independent seeds and more tracers, before a marginal pass is read.
- **F3.** Its 2 kpc bins hold 50–106 tracers, so it cannot show 5% convergence. It is to be repeated with more tracers.
- **The profile rows.** They are to store each shell's edges and centre separately, and the anisotropy is to separate coherent infall from random radial motion.

The seed-order defect found in F1's driver affects stage 2A's too: its rounds 2–4 draw seeds in the order runs finish, so its archived canonical run (7620668) is not bit-reproducible.

## R1. The direct normalization search

It covers every combination with at least one completed run in the canonical 2B-F1 results, and the 150 kpc sensitivity.
1. **The best sampled rate, q_s.** This is the completed run with the lowest actual RMSE over the 38 Eilers bins, among the ladder, extension, verification and refinement runs. Its scores and gates are reported as they are.
2. **The bracket.** Five rates, q_s·3^(k/4) for k = −2, …, 2, each with three independent seeds at the canonical tracer numbers.
   - q_s itself is rerun with new seeds, so the selected value does not inherit the favourable noise that made it the best sample.
   - Rates at or beyond one where a run of the combination stopped early (runaway, opacity or a lost static bath), on the far side of q_s, are not run.
3. **Extension.** If the lowest seed-mean RMSE lies at an end of the completed rates, the bracket extends one step of 3^(1/4) outward, with three seeds. It does so at most twice per side, and never past a rate where a run stopped early.
4. **Refinement.** Two more rates, a step of 3^(1/8) either side of the rate with the lowest seed-mean RMSE, three seeds each.
5. **The selected rate, q_sel,** has the lowest seed-mean RMSE among the evaluated rates at which every seed completed.
6. **Scores at q_sel.** The 38-bin RMSE, the inner-20 RMSE, the 8–20 kpc slope and the cost are each reported as a seed mean with its standard deviation. So is each seed's own gate outcome. The gates G1–G3 are those of 2B-F1, applied to the seed means.
7. **The rescaling-root verification** of 2B-F1 stays in the results beside it, as a diagnostic.

**The revision's profile-gate verdict is the one at q_sel.**

## R2. Repeats at higher resolution

At q_sel, every combination whose seed means pass G1 and G2, whatever G3 says, is repeated with four more independent seeds at four times the tracer numbers.
- **Reported:** each high-resolution seed's scores and gates; their means and standard deviations; and how far the means move from the canonical-size seeds.
- **Declared label "robust pass":** the high-resolution seed means pass G1, G2 and G3, and at least three of the four seeds pass all three individually.

## R3. F3 at higher resolution

The cold-birth density check is repeated with 16 times F3's tracers (1,024,000), under the same frozen potential and at v_d = 1 km/s.
- **The tolerance** stays three standard errors or 5% per bin, whichever is larger. With per-bin standard errors near 4%, the 5% bound now binds.
- **It passes** when every 2 kpc bin across 5–25 kpc meets it.

## R4. Profile rows and the anisotropy

- **Radii.** `formation.summarize` stores each shell's `r_lo_kpc`, `r_hi_kpc` and `r_center_kpc`. `r_center_kpc` is the geometric centre, previously `r_kpc`.
  - The enclosed-mass fields stay at `r_hi_kpc`.
  - Shell density and velocity statistics stay at the centre.
- **Anisotropy.** Each shell adds its mass-weighted mean radial velocity, the radial dispersion about that mean, and the anisotropy computed from that dispersion, alongside the raw second moments.
- **F1's outer table** is regenerated from `r_hi_kpc`.

## R5. Reproducible drivers

Stage 2A's `run_tasks` returns results in the order the tasks were issued, as F1's now does.
- **Two concurrent stage 2A smoke runs** must agree exactly.
- **Stage 2A's canonical run is repeated.** The report compares the new archive with 7620668's and lists every declared benchmark, label and headline range that changes. The stage 2A report is updated from the new archive.
- **F1's canonical run is repeated,** for the new profile fields. Every score must reproduce ed96b00's exactly; that is the fix's end-to-end determinism check.

## R6. Universality

If the best-scoring source changes under the revised verdicts, J1630 and Coma are rerun, unchanged, with the new best source at its q_sel. Otherwise the 2B-F1 universality runs stand.

## Declared labels and reporting rules

- **"Profile gate passed (revision)"** at q_sel, on the seed means, and **"robust pass"** as in R2.
- **Three quantities stay separate:** positive inventory, gravitational mass contrast and the source's energy expenditure.
- **Negative search results** read "no sampled rate passed".
- **The Jeans growth time** is reported as a growth-timescale diagnostic.
- **Outer radii** come from `r_hi_kpc`.

## Not changed

- the source and its kinematics;
- CF-1's elastic law;
- the three systems and their zones of influence;
- the 10 Gyr span;
- the gates and their thresholds;
- the validation tolerances F1–F4 and V5–V6;
- the early stops.

## Not claimed

The time-dependent donor-and-companion calculation, boundary convergence, three-dimensional stability and finite source histories. Those are the queue's next items, each with its own protocol.

## Files

- **In `companion-source/`:** `revision.py` (R1–R3 and R6), `revision-results.json`, a revision section in `report.md`, and a refreshed `f1-results.json`.
- **In `companion-formation/`:** `formation.py` (R4, R5), with a refreshed `formation-results.json` and report.

## Amendment 1 (2026-09-15), after the owner's review of 9232e07

Declared before any revision run completed. The first canonical attempt started at 11:57 and was stopped at about 12:15 on 2026-09-15, after the review found the F3 cap below. It wrote nothing. Its seeds come from the revision's own counter in a fixed order, so the restarted run draws the same seeds for rounds A–D.

### A1. F3's tolerance, stated correctly (replaces R3's last sentence)

The tolerance stays three standard errors or 5% per bin, whichever is larger. That is a statistical check, not a 5% convergence test.
- **Where the 5% term would control.** Only where the relative standard error falls below 5%/3 = 1.67%.
- **At the expected resolution.** 2B-F1's bins had standard errors of 12–20% of the quadrature. An uncapped pilot before this amendment gave the following:
  - It ran 128,000 births, with no thinning.
  - It gave 178–366 effective samples per bin, which is relative standard errors of 5.2–7.5%.
  - Scaled to the revision's 1,024,000 births, those errors become 1.9–2.7%. The allowance is then about 6–8%, so the three-standard-error term still controls.
- **A genuine 5% test** needs trajectories sampled preferentially where they reach the measured bins at T, with correct weights. It is a separate item, not part of this revision.

### A2. F3's tracer cap

2B-F1's F3 helper passed no cap, so the engine's default of 40,000 applied. The engine thins by merging pairs whenever the population exceeds its cap. So 2B-F1's F3 requested 64,000 tracers and ran at no more than 40,000, and the revision's 1,024,000 would have ended the same way.
- **The fix.** The revision passes a cap of four times the requested births, so nothing thins.
- **What it reports:** the thinnings (required: none), the birth-mass doublings, the cap used, and each bin's effective sample size, N_eff = (Σw)²/Σw².
- **F3's budget** is 12 hours, in its own process.
- **The disclosure.** A diagnostic rerun of 2B-F1's configuration, made before this amendment, measured how much it thinned.
  - The population passed 40,000 twice and was halved by pair-merging each time. It ended with 21,169 tracers.
  - Its 2 kpc bins held 49–94 tracers, with effective sample sizes of 32–66 and relative standard errors of 12–18%. That agrees with the archived 12–20%.
  - The run took 382 s. An uncapped run at 1,024,000 births is estimated at two to three hours.

### A3. Numerical controls beyond tracer count (new R7)

Four times the tracers tests one source of numerical error. For the reference case (3 km/s, collisionless, bath gravity omitted), three seeds are run at its selected rate under each of three changes, one at a time:
- the timestep halved;
- the potential grid doubled (2,048 nodes);
- the source pools regenerated every 5 steps instead of 10, with twice as many draws per pool.

Each control's seed means are compared with the canonical-size seed means at the same rate. They agree when they differ by less than two standard errors of the difference, or by less than 1 km/s in RMSE and 0.05 in slope, whichever is larger. Round D's result is reported as tracer convergence only.

### A4. The best-source rule

Within a passing tier, round E's rule prefers the lowest net cost. It is kept for this declared experiment. But net cost includes the negative background contrast, so the rule can reward stronger cancellation rather than a smaller source.
- **So it does not choose the mechanism investigated next.**
- **The report compares the credible candidates separately:** numerical stability, positive inventory, source energy and transferred predictions.
- **The reference case.** The report treats cold collisionless decay (3 km/s) as the main physical reference, with the other combinations as comparisons.
