# Formation versus retention after source removal

13 September 2026. Controlled-history consequence of the threshold interpretation.

## Result

For all 149 source-intensity values in the existing sample, the first 10% of stored energy is lost after source removal in only 4.92% to 41.79% of the time used to form it. A much slower remainder persists: a long release tail does not imply near-permanent retention of the full deposit. This conclusion follows from the stipulated threshold kinetics, not from observed fading histories or a constraint on the universe's age.

## Controlled protocol and formulas

Use the same 12-decade capacity-weighted threshold distribution as the previous galaxy comparison. Start empty, maintain constant X, and stop the source individually when the occupied fraction reaches 90% of its equilibrium value. If that duration is u_90=lambda T_build, the subsequent occupancy at v=lambda T_dark is

\[
F(v)=\int \rho(t)\frac{X}{X+t}[1-e^{-(X+t)u_{90}}]e^{-tv}\,dt.
\]

Solve F(v_90)=0.9 F(0) and F(v_half)=0.5 F(0). These exponentials are known solutions of the local capture/release equation; their application and the threshold distribution are hypotheses.

Because both times use the same lambda, v_90/u_90 is independent of the unknown absolute rate. Changing lambda cannot change this ratio. It can change the elapsed years, which remain unspecified.

## Sample examples

| Source proxy from galaxy | X | Formation u_90 | First 10% loss v_90 | Half-loss duration | v_90/u_90 |
|---|---:|---:|---:|---:|---:|
| UGCA444 | 0.017419 | 102.462 | 5.03822 | 168.021 | 0.04917 |
| UGC11914 | 25.1995 | 0.084782 | 0.035429 | 2.17917 | 0.41789 |

The names identify the catalogue inputs used, not measurements of these galaxies' actual storage or source shutoffs. No rotation curve is refitted here. Across all sample inputs the individual first-10%-loss/formation ratio lies between those two values. Half-loss times can be substantially longer, consistent with the distribution's long tail.

## Rate-window interpretation

For a chosen maximum formation duration T_b and a desired post-source interval T_d retaining at least 90% of the formed deposit, this controlled protocol requires

\[
\lambda\ge u_{90}/T_b,\qquad \lambda\le v_{90}/T_d,
\]

hence T_d/T_b<=v_90/u_90. If the same lambda and deadline bounds must work for all these separate, individually stopped experiments, the intersection requires T_d/T_b<=0.000345779. The tightest formation and early-loss bounds come from different inputs.

This last statement does not describe galaxies illuminated for one common duration and switched off simultaneously. Their final threshold populations would differ and need a new calculation. It also does not apply unchanged to partial fading, ongoing external companion supply, initially loaded storage, or a different release law.

## Consequences

The equilibrium one-third representation alone cannot justify permanent storage. This particular local kinetics needs continued supply to maintain occupancy, or an additional physical protection mechanism. Changing release kinetics while preserving the equilibrium curve remains possible, but is a new assumption whose preparation and energy accounting must be tested.

No fixed cosmic age is used. Actual source histories remain unknown here. The reference formula, existing rotation predictions and all other theory branches are unchanged.

## Verification

The executable evaluates all 149 existing intensity inputs and uses the prior numerical formation crossings. It checks each initial occupancy against 90% equilibrium, solves both post-source thresholds, and verifies the resulting fractions within 1e-8. All times are explicitly dimensionless. The results establish a conditional history tradeoff, not observational validation.

Files: `threshold-memory.py`, `threshold-memory-results.json`. Related: [formation-history calculation](threshold-history-report.md).
