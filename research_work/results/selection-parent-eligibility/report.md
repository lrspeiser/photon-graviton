# Pilot parent eligibility: reconstructed colors, unresolved infrared cut

The previous turn recovered drilled plate membership. This analysis applies published targeting criteria to the 189,661-row photometric parent and audits them against the 1,271 observed main targets. No stellar velocities or gravity parameters enter this calculation.

## Established agreement

The conventional extinction correction is `(J-K)_0 = J-K-1.5 A_K`, obtained from `A_J=2.5 A_K` in the [APOGEE catalog reader](https://github.com/jobovy/apogee/blob/main/apogee/tools/read.py). This is a known survey convention, not a new formula in our hypothesis. All 1,271 observed stars have RJCE_IRAC extinction labels, and every reconstructed color lies within its recorded targeting bounds (tolerance 1e-6 mag).

All observed targets pass the 0.8-degree design radius, finite photometry/extinction, JHK errors, quality/read/confusion flags, galaxy-contamination flag, negative extkey sentinel, neighbor-separation and overall H/color cuts used here. The generic instrument footprint is not substituted for the actual design radius.

Three observed stars have stored H exactly equal to the float32 value of 12.8 and belong to the medium cohort. Using lower-exclusive, upper-inclusive H bins reproduces the six recorded observed cohort/color counts when the problematic infrared cut is omitted. This is an explicit endpoint convention, supported by these stored cases; original unrounded selection values are unavailable.

## A failed eligibility reconstruction

[Zasowski et al. Table 2](https://arxiv.org/html/1708.00155#S4.T2) gives a mid-infrared uncertainty limit of 0.1 mag. The [final southern targeting paper](https://arxiv.org/html/2108.11908) refers back to those quality requirements. Applying this threshold to TARG_4_5_ERR rejects 42 actual observed main targets. Their IRAC_4_5_ERR values also exceed the threshold: simply switching these columns does not repair the discrepancy. The remaining cuts reject none of the observed main targets individually.

| Candidate rule | Parent survivors | Observed main survivors |
|---|---:|---:|
| Including the 0.1 mag infrared-error cut | 12,150 | 1,229 |
| Omitting only that cut, for sensitivity | 13,193 | 1,271 |

The second row is not silently adopted as the true selection. It changes the candidate denominator by 1,043 stars (about 8.6% of the stricter count), with unequal effects across cohorts. In particular, the strict cut rejects one of the two actual long-cohort stars. A universal completeness ratio from either row would hide these differences.

## Consequence and remaining work

The dust-color conversion now reproduces the observed targeting intervals. The parent pool is far smaller than the raw catalog, but its exact eligibility is not settled. Trace the historical science input lists and targeting implementation to determine how infrared errors were used. Central/off-axis camera exclusions, fiber competition, repeated cohort versions and project-specific selection also remain to be handled. Camera geometry cannot be replaced by a circular radius alone.

We have produced no selection weights, intrinsic density estimate or gravity likelihood. Matching the observed targets is a necessary check, not proof that all and only eligible unobserved targets are included. Full three-dimensional tracer modeling and all six scientific goals remain open.

Run `python research_work/results/selection-parent-eligibility/run.py`. Results record input hashes, every cut's individual and cumulative effects, all 42 rejected identifiers, endpoint cases, and six bin counts for both sensitivity variants. This run changes no training/holdout assignments and reads no held-out kinematics.
