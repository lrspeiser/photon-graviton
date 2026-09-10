# The pilot's stellar orbit populations are mathematically admissible

The three constant-anisotropy stellar populations used in the lens pilot admit nonnegative distribution functions in the stipulated potentials. This strengthens the previous Jeans-moment calculation: it is possible to populate the model with stars having the assumed density and anisotropy, rather than merely obtaining positive-looking velocity averages.

This is a **conditional mathematical consistency result**, not new observational agreement. It does not establish stability, a distribution of companion deposits, observed stellar anisotropy, or the photon-production/capture/metric equations. The earlier unfavorable isotropic validation remains unchanged, and no additional held-out scores were opened.

## Exact scope

The tracer has the spherical Hernquist profile j(r) proportional to 1/[r(r+a)^3]. The gravitational field has positive ordinary-matter amplitude with g_b proportional to (r+a)^(-2), plus either no companion term or a positive extra acceleration g_c proportional to g_b^p inside a finite radius. Outside that radius, the extra force keeps its enclosed equivalent mass fixed and falls as r^(-2).

The empirical p=0.4624587420 is unchanged. The tested beta values are −0.3, 0 and +0.3, and the existing source radii are 5, 20 and 100 effective radii, with R_e=1.8153a. The analytic sign bounds below hold for every positive mass/coupling normalization in this family, not just the 33 fitted training masses. This does not apply automatically to a bar, arbitrary tracer profile, evolving field or a different truncation prescription.

## Known inversion criterion

Let Ψ(r)=Φ(infinity)−Φ(r) be relative potential and E=Ψ−v²/2 the binding energy per unit mass. For constant beta, use the **known spherical distribution-function form**

\[
f(E,L)=L^{-2\beta}f_E(E),\qquad h(\Psi)=r^{2\beta}j(r).
\]

L is specific angular momentum. For −1/2<beta<1/2 and the outer conditions h(0)=h'(0)=0, a known inversion is

\[
f_E(E)=\frac{1}{2^{3/2-\beta}\pi^{3/2}\Gamma(1-\beta)\Gamma(1/2+\beta)}
\int_0^E h''(\Psi)(E-\Psi)^{\beta-1/2}\,d\Psi.
\]

The prefactor and kernel are positive. Therefore h''(Ψ)>=0 is sufficient for a nonnegative stellar distribution function in this range. This is established distribution-function mathematics, associated with generalized Eddington/Cuddeford inversion; it is not a new photon-companion law. See [Cuddeford (1991)](https://adsabs.harvard.edu/pdf/1991MNRAS.253..414C) and [Ciotti & Morganti's consistency criteria](https://academic.oup.com/mnras/article/401/2/1091/1152199).

The present potentials have finite depth relative to infinity because their exterior accelerations fall as r^(-2). At large radius h is proportional to r^(-(4−2beta)), hence to Ψ^(4−2beta), which supplies the required outer conditions. The central tracer cusp can make f diverge near maximum binding energy; nonnegativity does not require a finite phase-space density at that endpoint. The spatial stellar mass remains finite.

## Applying the criterion to the actual pilot

Set x=r/a and k=1−2beta. Apart from positive constants, h=x^(-k)(1+x)^(-3). Write

\[
\ell=\frac{d\ln h}{dx}=-\frac{k}{x}-\frac3{1+x}.
\]

Since dΨ/dx is a negative positive-factor multiple of g, the sign of h''(Ψ) is the sign of

\[
Q=\ell^2+\ell'-\ell\frac{g'}g.
\]

Primes on g here denote x derivatives. The following bounds are direct algebraic deductions for this stated model, not a novelty claim or an inferred physical mechanism.

Inside the companion boundary, g'/g=−s/(1+x), where 2p<=s<=2 because g is a positive sum of powers with those slopes. Consequently

\[
Q\ge\frac{k(k+1)}{x^2}+\frac{4k}{x(1+x)}+\frac6{(1+x)^2}>0.
\]

Outside, the baryonic slope is −2/(1+x) and the companion slope is −2/x. A sufficient lower bound is

\[
Q\ge\frac{P_k(x)}{x^2(1+x)^2},
\quad P_k(x)=(k+2)(k+3)x^2+2(k+3)(k-1)x+k(k-1).
\]

For beta=−0.3 and 0, this polynomial is positive for x>0. For beta=+0.3 it becomes 8.16x²−4.08x−0.24 and is positive for x>0.55317. The smallest pilot boundary is x_t=9.0765, well beyond that sufficient threshold. Thus the inner and exterior bounds cover all radii for every tested case.

The extra force is continuous at the source boundary. Its derivative may jump, but h and its first potential derivative remain continuous, so there is no unaccounted negative delta contribution to the inversion. The nonnegative curvature holds on both sides. A smooth source boundary would need its own check; it is not silently assumed to inherit this proof.

## Numerical verification

The code checks 54 force configurations: three beta values, three boundaries, and six relative companion strengths from zero through 10^6, including 10^-6. Each has 4,003 radial probes, with explicit samples immediately to both sides of the boundary. The analytic lower bounds and direct augmented-density curvature agree in sign and ordering.

The inverse-transform normalization is checked independently using h(Ψ)=Ψ^n with a known power-law inverse. Nine energy/beta cases are inverted numerically and re-integrated to recover h. Maximum inverse relative error is below 2.5×10^-12, and maximum density-recovery relative error below 7.1×10^-11. These power-law checks validate the transform normalization; they are not additional galaxy fits. The global admissibility conclusion rests on the analytic bounds, rather than extrapolating positivity from a finite radial grid.

## What is still missing

The result establishes existence of the **stellar tracer** distribution function in the specified potential. It does not establish that this distribution is stable, that observations favor beta=+0.3, or that the extra source has a physically supported distribution. In particular, the companion source remains prescribed; its capture, transport, pressure/support and energy accounting are not supplied by the stellar inversion.

We should therefore keep radial and tangential orbital alternatives available for a properly constrained population analysis, without declaring either measured or using an independent beta for each galaxy to force agreement. The next decisive evidence is independent stellar-population/orbital information together with a shared companion field equation. More successful moment fits alone would not complete the theory.

Reproduce with run.py. It reads the archived p, records its input hash, verifies the polynomial bounds and transform identities, and saves all curvature checks. No training coefficients, source catalog entries, original validation outcomes or test-role assignments are changed.
