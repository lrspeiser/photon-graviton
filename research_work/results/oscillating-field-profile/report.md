# A single stationary wave cannot support the current fitted profile

The simplest oscillating-field completion fails a necessary equilibrium test in the interior of the archived lens models. Between 0.1 and 4 effective radii, the field mass required to support the prescribed density changes by a factor of **75.5 to 90.2**, with median **82.0**, across 99 training configurations. A single species must have one constant mass.

This is not an observational exclusion of all wave reservoirs. It shows that our empirical gravity profile is not already a solution of the simplest free stationary wave equations. An oscillating field's approximately pressureless average stress was a useful possibility, but it did not establish spatial support.

## Candidate and provenance

Assume a weak, nonrelativistic, minimally coupled massive scalar envelope described by the **known Schrodinger-Poisson system**:

\[
i\hbar\partial_t\psi=-\frac{\hbar^2}{2m}\nabla^2\psi+m\Phi\psi,
\qquad \nabla^2\Phi=4\pi G(\rho_b+m|\psi|^2).
\]

The envelope can describe a slowly varying oscillating real scalar to leading order or a complex scalar in the nonrelativistic limit. Here there are no self-interactions, capture sources during the equilibrium test, external confinement or modified gravitational coupling. We assume one stationary spherical state with positive amplitude, no nodes and no bulk flow. This is considerably narrower than arbitrary evolving or multistate fields.

The wave equations and their fluid representation are established mathematics; see [Construction and Evolution of Equilibrium Configurations of the Schrodinger-Poisson System](https://doi.org/10.3390/universe8080432). We borrow those equations without adopting expansion, a cosmological formation history or a dark-matter interpretation. Identifying their field with photon-fed deposits remains a hypothesis without a derived capture interaction.

## Invert the equilibrium requirement

Write psi=sqrt(rho/m) exp(-i mu t/hbar). The stationary equation becomes

\[
\frac\mu m=\Phi-\frac{\hbar^2}{2m^2}\frac{\nabla^2\sqrt\rho}{\sqrt\rho}.
\]

The left side is spatially constant. Differentiating eliminates the unknown eigenvalue and additive potential constant:

\[
g(r)=\frac{\hbar^2}{2m^2}\frac{d}{dr}
\left(\frac{\nabla^2\sqrt\rho}{\sqrt\rho}\right),
\qquad g=\Phi'>0.
\]

Consequently the following **conditional inverse diagnostic**, derived from the known equation, must be one positive constant:

\[
m_{\rm required}^2(r)=\frac{\hbar^2}{2g(r)}\frac{d}{dr}
\left(\frac{\nabla^2\sqrt\rho}{\sqrt\rho}\right).
\]

A negative right side would fail immediately. A positive but strongly radius-dependent value also fails. This does not measure a particle mass; it tests whether the stipulated density and total gravity could belong to this candidate. A constant value would still require boundary conditions, normalization, self-consistency and stability checks before acceptance.

## Apply it to the existing source, without refitting data

Inside the cutoff, the lens pilot specifies

\[
g=g_b+g_c,\quad g_b=GM/(r+a)^2,\quad
g_c=Aa_* (g_b/a_*)^p,
\]

\[
\rho_c=\frac{g_c}{4\pi G}\left[\frac2r-\frac{2p}{r+a}\right],
\qquad p=0.4624587420.
\]

The extra-force relation is an empirical fit. The density follows conditionally from ordinary Poisson gravity; neither is a microscopic wave law. The calculation uses rho=rho_c and the total baryon-plus-companion g. Using only the extra force in the equilibrium denominator would incorrectly omit the ordinary galaxy's pull on the wave.

With x=r/a, the density shape is proportional to

\[
f(x)=(1+x)^{-2p}\left[\frac2x-\frac{2p}{1+x}\right].
\]

Set ell=(ln f)'/2. Then the dimensionless Laplacian ratio is ell'+ell^2+2ell/x. Its derivative is evaluated symbolically and independently checked by high-precision differentiation of sqrt(f). The factor converting that derivative into physical units is a^(-3).

The inputs are the 33 training lens geometries and updated I-band sizes, with all three archived cutoffs for the 1.5-arcsecond seeing case. The ordinary-matter acceleration normalization is reconstructed from each previously saved dimensionless companion strength. These masses were inferred from aperture dispersions and are not independently measured stellar masses.

Eighty radii per configuration span 0.1 to 4 R_e, below the smallest cutoff of 5 R_e. This deliberately avoids the previously identified sharp-edge obstruction. The inverse mass is real and positive at every probe, but its largest/smallest ratio within each configuration is:

| Statistic across 99 configurations | Required mass variation |
|---|---:|
| Minimum | 75.49-fold |
| Median | 82.00-fold |
| Maximum | 90.19-fold |

These are 99 model configurations of 33 systems, not 99 independent observations. The result is an exact-profile compatibility test; observational and ordinary-matter uncertainties have not been propagated into a model-selection likelihood.

For an additional scale diagnostic, the code chooses one constant mass per configuration minimizing unweighted squared log acceleration residuals over these same radii. Even granting this otherwise unjustified per-galaxy freedom, predicted wave support divided by required acceleration ranges from about 0.0104 to 84.5 over all probes. The sampling and logarithmic objective are computational choices, not measurement weights. No such masses are adopted or used to update galaxy predictions.

## Why the mismatch is structural

For a power-law density rho proportional to r^(-2k), known differentiation gives

\[
\frac{\nabla^2\sqrt\rho}{\sqrt\rho}=\frac{k(k-1)}{r^2},
\qquad g_{\rm wave}=\frac{\hbar^2 k(1-k)}{m^2 r^3}.
\]

For 0<k<1 this support has the required sign, but falls as r^(-3). In the outer power-law approximation to the empirical interior profile, k=p+1/2=0.96246, while companion-dominated gravity falls approximately as r^(-2p)=r^(-0.92492). Thus m_required is proportional to r^(p-3/2), rather than constant. This asymptotic explanation is conditional on that radial regime; the quoted numerical results use the full profiles including baryons.

The prescribed central cusp and abrupt cutoff raise additional wave-regularity problems, but this interior result does not depend on reaching either endpoint.

## Consequence for the next calculation

We should not assign a radius-dependent particle mass or an independently adjustable support term to force the empirical profile into the wave equation. The scientifically useful alternatives are to solve for the density from a chosen action, or specify and justify an additional term before fitting:

1. A self-consistent stationary wave solution produces its own core and outer shape; compare its predicted motions and lensing with the data instead of insisting on the empirical profile exactly.
2. A self-interaction adds a density-dependent contribution to the equilibrium equation. Its coupling must be common and its stability checked, rather than chosen separately at each radius.
3. A mixture of occupied states or time-dependent wave motions can support a different average density, but rho is then a sum or time average; inserting sqrt(total rho) into a single-state equation is generally invalid. State occupations and formation need an explicit model.

All three routes still need a photon-to-companion capture mechanism, conservation and the same metric response for stars and light. No new observational agreement or holdout result is claimed. Existing negative lens validation remains part of the evidence.

Reproduce with `python research_work/results/oscillating-field-profile/run.py`. Dependencies are NumPy, SymPy and mpmath. The symbolic expression, input hashes, per-configuration diagnostics and derivative checks are archived alongside this report. Independent 45-digit amplitude differentiation agrees with the symbolic derivative to better than 1e-15 relative error at six separated radii.
