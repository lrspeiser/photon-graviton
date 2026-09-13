# Finite formation history on the frozen galaxy sample

13 September 2026. Exposed sample sensitivity analysis, no optimization.

## Result

The threshold-mixture interpretation requires time to approach the equilibrium used in the reference comparison. With one common microscopic time scale, the duration to reach 90% equilibrium spans 0.08478 to 102.462 dimensionless units across the 149 galaxies: a factor of about 1209. Low source intensity can delay accumulation substantially. These units are not years; the microscopic rate has not been derived.

## Hypothesis and calculation

Retain the prior threshold distribution truncated to [10^-6,10^6], its normalized capacity weighting, and the ordinary local rate equation. Starting empty under constant X gives

\[
f_t(u)=\frac{X}{X+t}[1-e^{-(X+t)u}],\qquad u=\lambda\,\Delta time.
\]

This is the standard solution of the stipulated linear capture/release equation. Averaging over thresholds gives the time-dependent occupied fraction. The equilibrium reference is approached as u grows. Since every local rate is at least X in these units, u_90 is no larger than ln(10)/X; the executable solves the actual mixture crossing.

The same durations 0.1, 1, 10, 100 and 1000 are applied to every galaxy. They are exploratory common histories, not fitted galaxy ages. The equilibrium case is also evaluated. All original ordinary-matter inputs, companion shape and normalization remain fixed; only the occupied-amplitude multiplier changes using the earlier squared-speed rescaling. Interpreting occupancy as that empirical multiplier remains a hypothesis.

## Rotation comparison

Scores retain the original equal-galaxy weighting and already-examined 89/29/31 partitions.

| Common duration u | Train RMSE (km/s) | Validation RMSE (km/s) | Test RMSE (km/s) |
|---:|---:|---:|---:|
| 0.1 | 33.8215 | 38.3154 | 30.9692 |
| 1 | 28.0350 | 34.1063 | 25.2714 |
| 10 | 28.9719 | 32.6056 | 23.5027 |
| 100 | 29.0345 | 32.5404 | 23.6186 |
| 1000 | 29.0325 | 32.5406 | 23.6186 |
| Finite-mixture equilibrium | 29.0325 | 32.5406 | 23.6186 |

For comparison, the original exact-third validation/test RMSE values were 32.4949 and 23.5908 km/s. The u=10 test RMSE is slightly lower, but its validation RMSE and logarithmic errors do not establish an improvement. No duration is selected from this scan. These are not new blind tests, and no parameter has been optimized on any split.

## Energy bookkeeping

For each threshold, k=X+t and f_eq=X/k. Cumulative capture and release, in units of energy capacity, are

\[
C_t=X[(1-f_{eq})u+f_{eq}(1-e^{-ku})/k],
\]
\[
R_t=t f_{eq}[u-(1-e^{-ku})/k],\qquad C_t-R_t=f_t(u).
\]

The same identity holds after averaging over thresholds. Maintaining equilibrium can therefore consume and release energy even while the stored amount remains constant. A good equilibrium curve is not evidence for a sufficient radiation supply or permanent storage. These formulas are derived from the stipulated local dynamics, not independently measured energy flows.

## Implications and verification

A future physical version needs a common rate scale and independently constrained source histories; otherwise age can become another fitting freedom. The actual sample spans enough intensity that assuming every galaxy reached equilibrium needs justification. The exact-third empirical reference remains unchanged. No age of the universe is fixed or excluded.

All original input hashes are verified. The equilibrium limit reproduces the prior finite-threshold split scores within 1e-8 km/s. The executable records 894 galaxy/history cases, 149 numerical 90%-formation crossings and six independent integrated capture-minus-release checks. Both velocity and logarithmic residual metrics are retained. These checks validate this conditional calculation, not its microscopic interpretation.

Files: `threshold-history.py`, `threshold-history-results.json`. Related: [finite-threshold galaxy comparison](threshold-galaxies-report.md).
