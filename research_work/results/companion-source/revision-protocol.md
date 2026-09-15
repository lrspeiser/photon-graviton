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
