# Exact-third response as a threshold mixture

13 September 2026. Inverse representation, not a derived microscopic law.

The reference eta(X)=X^(1/3)/(1+X^(1/3)) can be represented exactly by a positive mixture of ordinary first-power saturation responses. This gives an interpretation involving varied capture thresholds, but the distribution is chosen to reproduce the exponent. It is not independent evidence for our theory.

## Local rule and known mixture identity

Assume a site with threshold t>0 obeys df_t/dtime=lambda[X(1-f_t)-t f_t]. This familiar capture/release equation gives equilibrium f_t=X/(X+t). X remains the empirical dimensionless intensity proxy; its connection to a local loading rate is still hypothetical.

For 0<q<1, choose the capacity-weighted density

\[
\rho_q(t)=\frac{\sin(\pi q)}\pi\frac{t^{q-1}}{1+2t^q\cos(\pi q)+t^{2q}}.
\]

The known positive-mixture identity is

\[
\int_0^\infty\rho_q(t)dt=1,\qquad
\int_0^\infty\frac{X}{X+t}\rho_q(t)dt=\frac{X^q}{1+X^q}.
\]

Setting q=1/3 reproduces our reference. In y=ln(t), its density is

\[
w_q(y)=\frac{\sin(\pi q)}{2\pi[\cosh(qy)+\cos(\pi q)]}.
\]

These formulas belong to known fractional-relaxation spectral representations; see [the distribution discussion in this research paper](https://www.mdpi.com/1996-1944/3/1/585). Their use for companion thresholds is our hypothesis. The mathematics is not unique to this project. One derivation takes the negative imaginary part of the Stieltjes function X^(q-1)/(1+X^q) at X=-t+i0, divided by pi.

The mixture describes equilibrium occupied capacity, not automatically the fraction of incident energy permanently captured. At equilibrium, capture and release balance. Mapping this occupancy to the empirical halo multiplier requires capacity, history and energy-flow accounting. The weighting need not equal a number distribution if individual capacities differ.

## Finite thresholds

Truncating the distribution to [10^-W,10^W] and renormalizing changes the curve. Over the tested X range 10^-6 to 10^6:

| W | Retained distribution mass | Maximum absolute response error |
|---:|---:|---:|
| 3 | 0.842833 | 0.063282 |
| 6 | 0.983543 | 0.007816 |
| 9 | 0.998347 | 0.000811 |

These are numerical cutoffs, not measured threshold ranges. No galaxy fit has been rerun or improved. A narrower application range may require less threshold breadth.

## Conditional release prediction

Starting from X=1 equilibrium and removing the source, the stipulated local dynamics gives f_t(u)=exp(-tu)/(1+t), where u=lambda times elapsed time. The ideal unbounded mixture has

\[
\bar f(u)\sim\frac{\Gamma(q)\sin(\pi q)}\pi u^{-q}.
\]

Thus q=1/3 produces slow power-law depletion, not permanent retention. Initial mean occupancy is 0.5. At u=10^3, 10^6 and 10^9, about 14.02%, 1.4695% and 0.1476% of that initial occupancy remains. Lambda is undetermined; these are not years.

This prediction assumes release rate lambda t and the same lambda across thresholds. The equilibrium curve alone does not determine kinetics. A nonzero lower threshold eventually cuts off the power-law tail. No cosmic age is imposed.

## Interpretation and verification

A physical explanation must predict the threshold distribution from environments or an interaction before applying it to a new galaxy. Choosing it from the desired curve is inverse design. Source-removal behavior and finite cutoffs offer ways to distinguish this interpretation from protected-state alternatives. The exact-third reference remains unchanged.

Twenty-five direct log-threshold integrals match the reference within absolute tolerance 1e-10. Three finite-cutoff mixtures are evaluated on the same grid. Four release calculations check the long-time coefficient, with the final case within 0.1% of its asymptotic value. These checks validate the representation and stipulated dynamics only.

Files: `threshold-mixture.py`, `threshold-mixture-results.json`.
