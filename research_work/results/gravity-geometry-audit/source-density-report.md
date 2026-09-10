# Can accumulated deposits produce these well shapes?

Both conservative companion-well constructions from the previous comparison require a positive equivalent source density at all 63 sampled locations. This removes one possible objection on that grid: neither needs local subtraction of a Newtonian source there. It does not derive the deposition process or prove positivity everywhere.

The main additional requirement is an outer boundary or transition if the model is to describe a finite isolated deposited reservoir. Extending the fitted power law indefinitely would require unbounded equivalent mass and would not admit the usual finite potential referenced to infinity. This is an extrapolation diagnosis, not a failure measured at a galaxy's edge.

## What was tested and why

We retained the same empirical coefficients, ordinary-matter model, and two conservative geometries. The grid covers R=0.5, 1, 1.5, 2, 3, 5, 8, 12 and 20 kpc, each at z=0, 0.1, 0.3, 0.7, 1.1, 2 and 4 kpc. Reflection symmetry supplies the equivalent negative-z result. Those are numerical field probes, not 63 observed stars.

**Known Newtonian Poisson relation, used conditionally:**

\[
\rho_{c,\mathrm{eff}}=\frac{\nabla^2\Phi_c}{4\pi G}
=-\frac{\nabla\cdot\mathbf a_c}{4\pi G}.
\]

This density tells us what an ordinary Newtonian source would need to look like to produce the proposed extra well. It is not automatically the deposited photon energy divided by c². The user's more general possibility—deposits altering a gravity response—would require its own field equation, and this equivalent density need not be a material density. We therefore do not use this test to impose a particle-mass interpretation on all versions of the hypothesis.

**Known radial differentiation applied to the spherical geometry postulate:**

\[
\rho_c^{\rm sph}(r)=\frac{1}{4\pi G}
\left[\frac{2g_c(r)}r+\frac{dg_c}{dr}\right].
\]

**Known chain rule applied to the proposed potential composition:**

\[
\Phi_c=F(\Phi_b),\qquad
\rho_c^{\rm shape}=F'(\Phi_b)\rho_b+
\frac{F''(\Phi_b)|\nabla\Phi_b|^2}{4\pi G}.
\]

The second term matters: multiplying ordinary-matter density by the force multiplier would not reconstruct the source of this well. These equations are established mathematics applied to the project's effective models, not new fundamental physics.

At the finest tested finite-difference step and L64 bar representation:

| Equivalent source | Smallest sampled density (solar masses/pc³) | Largest sampled density | Negative sampled locations |
|---|---:|---:|---:|
| Spherical extra well | 0.0011184 | 0.44950 | 0 |
| Potential-composition extra well | 0.0010676 | 0.63770 | 0 |

The underlying ordinary-matter source is explicitly nonnegative, but differentiating its finite-resolution cached potential gives a small negative result at R=0.5, z=4 kpc: approximately −0.0000273 solar masses/pc³. This is a numerical reconstruction artifact, not negative matter in the input model. It is retained in the result files, not clamped away. We separately evaluate the analytic bar, nuclear, disk and softened-central-source densities and substitute that density into the chain rule as a sensitivity check. That substitution is not silently presented as the exact Laplacian of the cached potential.

We compare steps of 0.004, 0.002 and 0.001 kpc and bar orders L40/L64, using 64 azimuth samples and 960 equatorial interpolation nodes. Independent acceleration-divergence and chain-rule calculations are compared in source-density-verification.json. These are finite numerical checks. Sparse-grid positivity does not exclude a negative pocket elsewhere; it does not establish a globally physical source, equilibrium or stability. The previous force-accuracy checks alone were insufficient to certify second spatial derivatives, which is why the explicit source comparison was necessary.

## What indefinite extrapolation would require

**Conditional deduction from the frozen empirical power law, using known Newtonian exterior gravity.** Suppose the ordinary matter has finite total mass M_b, so sufficiently far away g_b≈GM_b/r². With the fitted p=0.4624587420:

\[
g_c\propto r^{-2p}=r^{-0.92491748},
\]
\[
M_{c,\mathrm{eff}}(<r)=\frac{r^2g_c}{G}
\propto r^{2-2p}=r^{1.07508252},
\]
\[
\rho_{c,\mathrm{eff}}\propto r^{-(1+2p)}=r^{-1.92491748}.
\]

The potential difference grows as r^(1−2p)=r^0.07508252. Thus its integral to infinity diverges, although slowly. The two geometries approach the same radial limit around a finite isolated ordinary-matter source.

In plain language, extending the formula farther keeps adding required deposits without ever finishing the reservoir. A finite isolated reservoir needs a turnover, truncation or a different exterior response. A model embedded in a cosmic companion environment could instead specify a non-isolated boundary condition, but then that environment and its field must be calculated. Neither solution follows from the existing rotation fit. This does not adopt an expanding universe, assign a cosmic age, or reopen the deferred total photon-source energy calculation.

## Consequences for the next physical model

A mathematically consistent finite-source continuation is available for the spherical benchmark. At an as-yet-unspecified radius r_t, stop adding source density while retaining the mass already enclosed. Outside it, use g_c(r)=g_c(r_t)(r_t/r)². The force is continuous, the exterior density is zero, and the exterior potential integral is finite. A smooth taper of source density can replace the abrupt edge. **This is known spherical Newtonian construction, not a prediction of r_t or of capture.** No cutoff radius is selected or fitted here.

Simply multiplying the acceleration by exp(−r/r_t) would not accomplish the same thing. In the asymptotic power-law example its inferred density changes sign when r/r_t exceeds 2−2p≈1.0751. That would remove enclosed equivalent mass as radius increases. The dimensionless calculation in exterior_check.py preserves this failed shortcut and verifies the positive-source continuation. This distinction explains why an outer transition must be specified at the source/field-equation level rather than attached arbitrarily to the acceleration.

The positive sampled densities are compatible with exploring accumulation; they do not make either geometry correct. The earlier vertical-force comparison still disfavors the particular shape-preserving completion under its fixed baseline. Capture, subsequent redistribution/support and the gravity response must together predict both the interior shape and the exterior transition. Picking an arbitrary radius after looking at the data would add a fitted parameter, not explain the mechanism.

No new astronomical data were fitted, no raw-star identities or distances changed, and no reserved stellar likelihood scores opened. Source positivity and a sensible exterior are necessary checks for a simple accumulated-source interpretation; successful redshift, event timing, lensing, transport, conservation and holdout predictions remain separate requirements.

Reproduce with source_density.py followed by verify_source_density.py, using the existing field caches and archived empirical fit. The full 378 numerical rows, summaries and resolution checks are stored beside this report. The [previous geometry comparison](report.md) contains observational provenance and references; the formulas here follow directly from the explicitly stated Poisson and power-law assumptions.
