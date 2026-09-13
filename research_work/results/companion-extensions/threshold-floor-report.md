# Continuing illumination after source fading

13 September 2026. Conditional extension of the threshold-retention branch.

## Finding

Continuing external supply can support an enduring deposit in the stipulated capture/release kinetics. Maintaining 90% of an initially equilibrated deposit requires the input to remain at 34.09%-68.62% of its original value across the existing 149 galaxy intensity proxies. Maintaining 99% requires 89.09%-96.43%. Thus an external background is a possible maintenance mechanism, but its required strength is calculable and is not automatically negligible.

These are source requirements, not measured external radiation levels, new galaxy fits, or evidence that distant stars satisfy them. X is the existing dimensionless source proxy. Identifying its change with a physical incident companion flux requires an explicit source-to-rate mapping.

## Protocol and derivation

Unlike the previous individually stopped 90%-formation experiments, begin here at full equilibrium under X. At dimensionless time u=0, reduce the source to a constant Y=beta X, with 0<=beta<=1. Preserve the normalized threshold distribution on 1e-6<=t<=1e6, the shared rate scale lambda, and the equation

\[
\frac{df_t}{du}=Y(1-f_t)-t f_t,\qquad u=\lambda\Delta t.
\]

The known linear rate-equation solution is

\[
f_t(u)=\frac{Y}{Y+t}+\left[\frac{X}{X+t}-\frac{Y}{Y+t}\right]e^{-(Y+t)u}.
\]

Every threshold population decreases monotonically to the new equilibrium. Therefore retaining at least a fraction p of the initial total occupancy at every subsequent time is equivalent to eta(Y)>=p eta(X). No cosmic age or fitted absolute rate is needed for this perpetual-floor criterion. A finite maintenance interval can require less supply and needs a specified duration.

For the infinite-support exact-third reference eta(X)=X^(1/3)/(1+X^(1/3)), algebra gives

\[
\boxed{\beta_{\min}=\frac{p^3}{[1+(1-p)X^{1/3}]^3}.}
\]

Provenance: this is an algebraic consequence of the known Hill response, not a new fundamental interaction. Applying it to companion storage is our hypothesis. The executable separately inverts the finite-cutoff integral used in the recent galaxy calculations rather than silently replacing that response by the exact infinite-support expression.

## Results

| Retained deposit | Required input fraction, finite thresholds | Infinite-support comparison | Continuing turnover power / original turnover power |
|---|---:|---:|---:|
| 50% | 0.00879-0.09889 | 0.00834-0.08672 | 0.02195-0.11133 |
| 90% | 0.34089-0.68624 | 0.33710-0.67512 | 0.44303-0.70349 |
| 99% | 0.89087-0.96426 | 0.88973-0.96279 | 0.91756-0.96668 |

The ranges span existing inputs X=0.017419 to 25.1995; they are not uncertainty intervals. A highly saturated system can lose more input while retaining a given fraction of its deposit. There is no parameter optimization and no new observational validation.

## Energy requirement

If E_cap is the energy represented by full capacity, the late-time capture and release powers are equal:

\[
P_{\rm cap}=P_{\rm release}=E_{\rm cap}\lambda Y[1-\eta(Y)].
\]

The equality follows by integrating t Y/(Y+t)=Y[1-Y/(Y+t)] over normalized capacity weights. The table divides this by the original equilibrium turnover E_cap lambda X[1-eta(X)]. The stored energy remains steady while energy continues flowing through it. These expressions retain the branch's energy-capacity convention; without E_cap and lambda they do not give watts. Release need not be electromagnetic heat, but any proposed outgoing channel must carry this energy. A static gravitational field itself is not assumed to require continuous power: the turnover follows from this particular leaking-reservoir model.

If initial input is X_local+X_external and only the local component vanishes, beta is the external fraction only when those input rates add linearly and the external source remains fixed. This provides a target for a transport calculation of distant stars. It cannot justify importing an arbitrary background or ignoring depletion, geometry, source spectra, or backreaction.

## What changes next

The earlier source-removal result remains valid. This calculation resolves the partial-fading loophole for a constant floor and initially equilibrated stores: supply can prevent further depletion, at the stated cost. Alternatives include physically protected states or altered kinetics, which must retain their own reverse transitions and energy accounting. None is established here. Actual source histories, a physical rate scale, external transport, and joint motions/lensing predictions remain outstanding.

## Reproduction and checks

Run `python research_work/results/companion-extensions/threshold-floor.py`. It records the history-input SHA-256 and evaluates 447 cases (149 inputs, three retained fractions). Checks recover the existing finite equilibrium values, solve each target occupancy, verify equality of capture and release, and evaluate monotonic trajectories at six durations. The full trajectories follow analytically from the expression above. These checks validate the calculation under its stated assumptions; they do not validate those assumptions.

Files: `threshold-floor.py`, `threshold-floor-results.json`. Related: [complete source removal](threshold-memory-report.md), [formation histories](threshold-history-report.md), [threshold representation](threshold-mixture-report.md).
