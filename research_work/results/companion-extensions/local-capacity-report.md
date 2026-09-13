# Local threshold filling with illumination-independent capacity

13 September 2026. Exposed-sample diagnostic; the reference is not replaced.

**Normalization correction (13 September 2026).** In the density, capacity and rate equations in this report, C denotes the pre-retention amplitude A=2 C0=9.457178483e7 Msun/kpc^3, where C0=4.728589242e7 is the stored fit parameter. C/original C multipliers are unchanged because the factor of two cancels. See [normalization audit](capacity-normalization-report.md).

## Outcome

Moving attenuation inside the one-third response produces much more extra gravity in weakly illuminated regions. Keeping all reference constants gives validation/test velocity RMSE 67.14/66.01 km/s, versus 32.49/23.59 for the reference. Reducing the shared amplitude using training galaxies gives a mixed result: validation improves to 28.63 from the matched control's 30.82, but test worsens to 25.19 from 21.41 km/s. This is not a general improvement or fresh observational evidence.

## What is being changed

The preceding capture/capacity audit identified a missing distinction: is available capacity itself attenuated, or is only its filling rate attenuated? Define the existing assumed profile g(r)=[1+(r/a)^2]^-2 and incident-field attenuation J(r). The reference density is

\[
\rho_{\rm ref}(r)=C g(r) J(r)\eta(X),\qquad \eta(X)=\frac{X^{1/3}}{1+X^{1/3}}.
\]

Here X is the existing luminosity-density proxy, not an independently measured companion bath. For this diagnostic assume instead that the energy capacity density is Cg, independent of illumination, while the local dimensionless loading input is XJ. Equilibrium filling then gives

\[
\rho_{\rm local}(r)=C g(r)\eta[XJ(r)].
\]

The profile g is the hypothetical capture profile, not an observed ordinary-matter density. The local loading rule and capacity interpretation are hypotheses. The one-third Hill response is known mathematics, previously represented by an inverse-chosen mixture of threshold responses; this substitution does not derive the exponent.

For each threshold t, the local rate equation is df_t/du=XJ(1-f_t)-t f_t. Its known stationary solution f_t=XJ/(XJ+t), averaged over the chosen threshold distribution, gives eta(XJ). No universe age is assumed; stationarity is a stipulated endpoint, not established formation history.

## Why the profile changes

Ordinary algebra gives

\[
\frac{\rho_{\rm local}}{\rho_{\rm ref}}=
J^{-2/3}\frac{1+X^{1/3}}{1+(XJ)^{1/3}}\geq1\quad(0<J\leq1).
\]

Thus a region receiving much less incoming intensity loses density more slowly under the local occupancy rule than under the original factor J. Both versions agree when J=1, so the transparent limit is preserved. The inequality concerns deposited density and extra gravity, not the entire ordinary-matter-plus-companion gravitational field.

The calculation retains the existing spherical transport field J computed from kappa=k0 g. It integrates both densities inside every accepted measured radius and uses

\[
v_{\rm local}^2(R)=v_b^2(R)+
\frac{M_{\rm local}(<R)}{M_{\rm ref}(<R)}[v_0^2(R)-v_b^2(R)].
\]

This uses the known spherical mass/acceleration relation and the reference's geometry. Ratios use the same numerical quadrature for both mass integrals. Measured baryonic components, distances, reference k0 and scale-to-disk remain unchanged.

## Numerical comparison

All scores below are equal-galaxy velocity RMSE in km/s across the existing 89/29/31 training/validation/test partitions. These partitions have been inspected repeatedly and are not blind evidence. No observational covariance or baryonic uncertainty is propagated in this diagnostic.

| Branch | C/original C | Train | Validation | Test |
|---|---:|---:|---:|---:|
| Frozen reference | 1 | 29.03 | 32.49 | 23.59 |
| Frozen local occupancy | 1 | 57.71 | 67.14 | 66.01 |
| Training-adjusted reference | 0.747109 | 26.87 | 30.82 | 21.41 |
| Training-adjusted local occupancy | 0.352688 | 28.75 | 28.63 | 25.19 |

The adjustments each fit one shared amplitude C on training velocity MSE only, with identical bounds 0<C/original C<4; both solutions are interior. The other two reference constants are frozen. Using this velocity objective changes the calibration criterion relative to the original fit, which is why an identically adjusted reference is included. There are no per-galaxy fitted parameters and neither amplitude change is adopted into the reference.

Before amplitude adjustment, the largest predicted speed change is 339.94 km/s for NGC6195. Enclosed extra-mass ratios over all sampled radii range from 1.020 to 476.129; large ratios can arise where the old attenuated component is small and must not be described as equivalent factors in total gravity. The pronounced effect is a formula consequence, not a small numerical perturbation.

## Transport and energy limitations

This is a fixed-incident-field calculation, not a self-consistent absorptive equilibrium. If occupied sites stop capturing, the effective capture opacity contains an empty-site factor, and the originally computed J generally changes. Alternatively, an independent attenuating channel could maintain the specified J, but its outgoing energy and source budget would have to be included explicitly. Neither closure is established here.

For a prescribed local input XJ, threshold equilibrium has equal capture and release: per capacity and per dimensionless time, both are XJ[1-eta(XJ)]. The overall calculation still needs a source-to-rate mapping, a physical capacity and lambda, and transport with returning or released companions. Local rate balance alone does not close that global ledger. Likewise, the modified stationary density does not demonstrate orbital support, stability, or lensing consistency.

## Interpretation and next discriminator

The unattenuated capacity alternative is mathematically usable, and its transparent limit agrees with the reference, but direct insertion changes rotation curves markedly. A training amplitude adjustment yields conflicting validation and test changes. Therefore retain this as a diagnostic branch and keep the original reference intact.

The next physical discriminator is occupancy-dependent transport: does reduced capture in filled regions change illumination enough to improve or further worsen the stationary mass profile? That requires a coupled transfer/occupancy solution, with conservation and appropriate outgoing channels, rather than further free adjustment of the present fixed field. Deriving actual capacity from a specified microscopic state model remains a separate requirement.

## Reproduction

Run `python research_work/results/companion-extensions/local-capacity.py`. It verifies catalogue/archive hashes and every accepted radius and observed velocity (149 galaxies, 3150 points), reproduces the baseline scores, and compares radial/angular Gauss-Legendre quadrature orders 64 and 128, with a programmed 256 fallback. All cases converge at 128; the largest 64-to-128 predicted-speed difference is 1.28e-9 km/s. This establishes integration convergence for the stated field, not physical correctness. Model-input hashes, per-radius predictions, extra-mass ratios and fitted scores are recorded in `local-capacity-results.json`.

Related: [capture/capacity audit](capture-capacity-report.md), [uniform-background comparison](background-retention-report.md).
