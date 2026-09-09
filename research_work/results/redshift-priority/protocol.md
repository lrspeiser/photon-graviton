# Frozen exploratory comparison, 9 September 2026

Written before running this comparison. Source: unchanged all_164_groups.csv, SHA-256 8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0. All 164 groups and all historical partitions have been exposed. No blind-validation claim.

Compare exactly three candidates, without choosing a winner inside this run:

1. Linear mathematical control: z=k D/c, k in [0,150].
2. Constant conversion: z=expm1(k D/c), k in [0,150].
3. Smooth rate: alpha(s)=[k0+(k100-k0)s/100]/c over 0-100 Mpc, both endpoints in [0,150] km/s/Mpc. Thus z=expm1([k0 D+(k100-k0)D^2/200]/c). This is a two-parameter empirical path-age/distance diagnostic, not a local environmental law. No extrapolation beyond 100 Mpc or claim of derived time physics.

Fit squared c*z residuals. Leave each entire existing sky tile out once, with all parameter fitting restricted to other tiles. Keep all rows; no outlier clipping, per-object rate, intercept, distance adjustment, environmental proxy from residuals, or hyperparameter selection. Compare every candidate's pooled out-of-fold RMS, MAE, median absolute error, bias, and per-tile metrics. Use 2000 paired sky-tile bootstrap samples, seed 2026090903, for descriptive RMS-difference uncertainty. It does not erase prior sample exposure.

Audit frozen baseline residuals by tile and quartiles of adopted distance, and their Pearson correlations with distance and the three Cartesian sky-direction components. These correlations are descriptive, not independent significance tests or evidence of void physics. No independently measured path-environment columns exist in the recovered sample; environmental tuning is deferred.

Record source selection and frame. A later fresh-data protocol must independently constrain nuisance scatter and set interval/precision criteria before opening outcomes. This run does not certify uncertainty coverage or select a final physical model.
