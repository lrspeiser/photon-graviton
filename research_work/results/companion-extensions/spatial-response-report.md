# Finite receiving size: forward scattering and residual color dependence

13 September 2026. Conditional spatial extension of the earlier scalar interaction.

## Result

A Gaussian receiving overlap can soften the point interaction's color dependence and make individual deflections small. It does not make the fractional loss achromatic: at high photon energy relative to the gap and inverse response size, the loss rate grows approximately linearly with photon energy. Large forward-concentrating regions also suppress the integrated event rate. The unknown coupling cannot be increased silently to compensate.

This is a tested spatial hypothesis, not a derived bound-companion wavefunction or a successful redshift fit. It does not establish that real companion stores have the specified sizes.

## Proposed spatial overlap and known scattering mathematics

Use natural units hbar=c=1 for energies and inverse lengths. **Hypothesis:** the transition density has a normalized isotropic Gaussian profile with spatial standard deviation a. **Known Fourier-transform consequence:** its amplitude is F(q)=exp(-a^2 q^2/2), and the squared amplitude is exp(-a^2 q^2).

For an incident photon E, outgoing photon E'=E-Delta, and cosine of scattering angle mu, **known momentum geometry** gives

\[
q^2=E^2+E'^2-2EE'\mu=\Delta^2+2EE'(1-\mu).
\]

Keep the prior unpolarized point-response angular factor 1+mu^2. Its [operator and rate derivation](../pair-production-balance/derivation.md) remain the underlying conditional model. The normalized angular suppression is

\[
S(E,a)=\frac38 e^{-a^2\Delta^2}
\int_0^2(2-2u+u^2)e^{-xu}\,du,
\qquad x=2a^2EE',\quad u=1-\mu.
\]

The forward fractional-loss coefficient becomes alpha(E,a) proportional to Delta*(E-Delta)^3*S(E,a). This is a consequence of the stipulated Gaussian extension, not an empirical law. The heavy receiving store absorbs recoil momentum; recoil energy is neglected as in the parent leading-order rate.

## Why the size does not fully remove color dependence

For x much greater than one, the integral is approximately 2/x. Therefore

\[
S\simeq\frac{3e^{-a^2\Delta^2}}{8a^2E(E-\Delta)},\qquad
\alpha\propto\frac{3\Delta e^{-a^2\Delta^2}}{8a^2}
\frac{(E-\Delta)^2}{E}.
\]

At E much greater than Delta, alpha is proportional to E. These are asymptotic derivations for this Gaussian model, not universal bounds on all finite targets. The interaction preferentially allows a narrowing forward cone; integrating over that cone removes two of the point model's three energy powers, leaving one. Additional frequency dependence or different dynamics is still required.

We evaluate five response sizes and six photon energies (0.01–100 eV) at an illustrative gap 1e-8 eV. Rate curves are normalized at 1 eV only to compare color dependence; no absolute coupling is fitted. The point limit is recovered. Analytic incomplete-gamma angular moments agree with an independently rescaled quadrature to better than 1e-10. The rescaling resolves extremely narrow forward cones without mistaking an unresolved numerical integral for zero scattering.

## Conditional accumulated deflection

For independent axisymmetric scattering events, the mean direction projected along the original direction after Poisson mean N is exp[-N*<1-mu>]. This is a known angular-random-walk identity. When deflections remain small, the accumulated direction RMS is approximately sqrt(2N*<1-mu>). It is not automatically the angular size of an observed source: a source-observer calculation needs scatterer locations, ray selection, path lengths and telescope geometry.

As a design example, use 2 eV photons, Delta=1e-8 eV, and the previous alpha*D=0.0076313 over approximately 100 million light-years. At leading small-loss order N=(alpha*D)*E/Delta is 1.53 million. The same fixed-gap Poisson approximation gives relative energy width about 6.18e-6. This is an illustrative calculation, not a measured linewidth constraint. Energy and rates evolving along the path are neglected at this order.

| Gaussian spatial standard deviation | Integrated rate / point rate at 2 eV | Conditional accumulated direction RMS |
|---|---:|---:|
| 0, point limit | 1 | Small-angle approximation invalid |
| 0.197 micrometers | 0.0835 | Small-angle approximation invalid |
| 0.197 millimeters | 9.37e-8 | Small-angle approximation invalid |
| 0.197 meters | 9.37e-14 | 127 arcseconds |
| 19.7 meters | 3.45e-18 | 1.27 arcseconds |

The event count is held at the redshift target for this comparison. With the same original coupling and target abundance, suppressed cross sections would instead produce fewer events and less redshift. Recovering the count requires a changed coupling, overlap normalization or number of receiving modes; none has been independently supplied. Large parameter changes could also invalidate weak-coupling or heavy-store approximations.

The e^[-a^2 Delta^2] factor expresses another tradeoff: even perfectly forward transfer changes photon momentum by Delta, so a very extended smooth overlap can suppress the energy transfer itself. Larger size is not a cost-free way to obtain both sharp images and sufficient conversion.

## Decision and remaining possibilities

Do not adopt this Gaussian spatial extension as a complete redshift solution. It improves directional concentration and softens the color problem, but leaves a substantial color dependence and an absolute-rate problem. It has not calculated event-duration stretching or a common clock law.

A useful next comparison is an explicitly energy-dependent polarizability or a coherent evolving background derived from the same source interaction. Such a replacement must predict both its dispersive and dissipative behavior; dividing out the remaining color dependence by hand would only impose the desired answer. The existing reverse-transfer and storage constraints must also remain attached to the interaction. No change is made to the exact-third retention reference or the gravity-fit records.

Reproduce with `python research_work/results/companion-extensions/spatial-response.py`. [Source](spatial-response.py), [all calculated values](spatial-response-results.json).
