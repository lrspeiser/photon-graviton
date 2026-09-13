# Does capture geometry change the background loading rate?

13 September 2026. Audit of the existing exact-third density profile and its possible storage interpretation.

**Normalization correction (13 September 2026).** In the density, capacity and rate equations in this report, C denotes the pre-retention amplitude A=2 C0=9.457178483e7 Msun/kpc^3, where C0=4.728589242e7 is the stored fit parameter. C/original C multipliers are unchanged because the factor of two cancels. See [normalization audit](capacity-normalization-report.md).

## Finding

If full storage capacity is identified with the existing fitted deposition profile at retention eta=1, its integrated capacity is proportional to its absorption cross section. Consequently, the incoming energy per capacity in a common isotropic bath is independent of galaxy size. Enlarging the capture area alone cannot repair the uniform-background problem under this capacity interpretation.

This is a conditional cancellation, not a proof that physical capacities must follow fitted profiles. The original profile is a deposition calculation; calling its eta=1 limit a capacity is an added storage interpretation. A genuinely independent capacity law could behave differently and must be derived and tested.

## Derivation and provenance

Retain the project's postulated opacity

\[
\kappa(r)=k_0 g(r),\qquad g(r)=[1+(r/a)^2]^{-2}.
\]

Let J(r) be the angular average of the attenuated incident intensity, normalized to the external isotropic bath. The reference deposited density has the form rho_d=C eta J(r)g(r). If eta is interpreted as the occupied fraction of an energy capacity, the corresponding full capacity density is C Jg (in mass-equivalent units).

Known straight-ray absorption gives the chord optical depth and effective cross section

\[
\tau(b)=T[1+(b/a)^2]^{-3/2},\qquad T=\pi k_0a/2,
\]

\[
\sigma=2\pi\int_0^\infty b(1-e^{-\tau(b)})\,db
=\pi a^2\{T^{2/3}\gamma(1/3,T)-1+e^{-T}\}.
\]

Here gamma is the lower incomplete gamma function. These are known integration and radiative-transfer mathematics applied to our assumed opacity, not newly discovered gravity laws. This profile has an infinite declining tail, so sigma is not the area of a fixed opaque sphere.

The transfer identity integral kappa J dV=sigma implies

\[
M_{\rm cap}=\int CJg\,dV=\frac{C}{k_0}\sigma.
\]

For a uniform isotropic energy density u_c and light-speed companions, the absorbed incident power before additional storage-availability factors is P_in=c u_c sigma. With E_cap=M_cap c^2,

\[
\boxed{\frac{P_{\rm in}}{E_{\rm cap}}=\frac{k_0 u_c}{C c}.}
\]

Use consistent SI units in this last expression: C is a mass density and k_0 is inverse length. The ratio has units of inverse time. No surface-area factor remains. The full absorption cross section already integrates directions; multiplying again by four times the projected area would double count. This is the gross loading normalization, not a claim that all arriving energy remains stored. Occupancy factors, outgoing energy and reverse processes still belong in the ledger.

Identifying this rate divided by a shared lambda with a dimensionless B is a further kinetic mapping. Under that mapping, geometry produces the same B across galaxies when C, k_0, u_c and lambda are common. Changing their fitted common normalization does not create environmental variation. Actual source fields, occupancy-dependent opacity and redistribution could invalidate the fixed-J assumptions and require solving transfer again.

## Numerical audit using current reference constants

Use C=2 C0=9.457178483e7 Msun/kpc^3, k_0=0.2039029004 kpc^-1 and a=2.770766589 R_d from the exact-third reference. This deliberately updates the constants relative to the older pre-third `energy.py` audit without altering that historical calculation.

| Quantity across 149 inputs | Minimum | Maximum |
|---|---:|---:|
| Capture scale a (kpc) | 0.49874 | 51.97958 |
| Central chord optical depth | 0.15974 | 16.64854 |
| Absorption cross section (kpc^2) | 0.24482 | 139773.47 |
| Inferred full capacity (Msun equivalent) | 1.1355e8 | 6.4828e13 |
| Area/capacity (kpc^2/Msun) | 2.15606e-9 | 2.15606e-9 |

These capacities are the eta=1 model integrals over the entire infinite profile, not measured galaxy masses or actual occupied deposits. The common area/capacity is exact algebraically; numerical values agree within floating-point precision. Both capture area and capacity can vary enormously while their ratio stays fixed.

## Alternative volume capacity and its cost

One alternative is rho_cap=Cg, independent of attenuation. Its total capacity is M_cap,vol=C pi^2 a^3. The loading rate now acquires a geometry factor

\[
\frac{P_{\rm in}}{E_{\rm cap,vol}}=\frac{k_0u_c}{Cc}\,S(a),\qquad
S(a)=\frac{\sigma}{k_0\pi^2a^3}.
\]

For the current scales, S ranges from 0.49454 to 0.98062. In the optically thin limit S tends to one, not an unrestricted inverse-radius enhancement. This alternative is a capacity hypothesis; it changes the relation between occupancy and the fitted spatial mass profile. It cannot be inserted as a background modifier while claiming that the original capacity and deposition equations are unchanged. A local occupancy solution would need to establish whether the previous J factor is recovered. No new rotation fit is performed or claimed for this alternative.

## Consequence for the research direction

The earlier common-background comparison was not resolved by invoking larger capture area. Under the capacity interpretation that preserves the reference profile, capture and capacity share that same area. Geometry can help only through independently specified capacities, nonuniform illumination, changed opacity or a spatially resolved storage solution. This audit selects the missing physical relationship to derive, rather than introducing a free per-galaxy correction.

The absolute external supply, microscopic rate, storage degrees of freedom, redshift/timing mechanism, and motions/lensing consistency remain unresolved. No universe age or size is assumed. The reference is not replaced.

## Reproducibility

Run `python research_work/results/companion-extensions/capture-capacity.py`. The calculation verifies original input hashes, records reference-model hashes, evaluates all 149 catalogue scales, and compares the incomplete-gamma cross section against independent impact-parameter quadrature to relative error below 1e-7. The capacity identity is an analytic transfer consequence, not an independent observation or a numerical simulation of storage. Files: `capture-capacity.py`, `capture-capacity-results.json`.

Related: [uniform-background comparison](background-retention-report.md), [historical capture-energy calculation](../isotropic-galaxy-transfer/energy-report.md).
