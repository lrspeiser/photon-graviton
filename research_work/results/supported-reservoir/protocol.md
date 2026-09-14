# Supported reservoir CR-2: one repulsive support law against joint lensing and stellar motions

Declared before execution, 14 September 2026. Baseline: `main` at 5dc018c.

**Why it is run.** The project owner's review requested it. CR-1's failures stand with their stated scope; CR-2 is a different model.

## Question

Can one declared repulsive support law, with one shared constant, supply the extra mass six SLACS lenses need? The test fits stellar motions and lensing jointly, in identical geometry and with identical mass conventions.

Three things are reported as separate results:
- support, meaning whether the law's equilibrium profiles fit;
- supply, meaning whether anything could deliver the mass;
- formation history.

## Declared law (support)

- **Matter.** A condensate of bosons with a repulsive contact interaction, in the Thomas–Fermi limit: P = (K/2)ρ², an n = 1 polytrope. K = g/m² is one constant shared by every system.
- **Equilibrium.** Hydrostatic equilibrium in the combined Newtonian potential of the host and the condensate: K∇ρ = −∇Φ, so ρ = (μ − Φ)/K where positive. In spherical symmetry,
  - dψ/dr = −G[M_h(r) + M_c(r)]/r²,
  - dM_c/dr = 4πr²ψ/K,
  - with ψ = μ − Φ, integrated from the centre to the edge where ψ = 0 ([tf.py](tf.py)).
- **The shared constant.** A condensate alone has ρ ∝ sin(kr)/(kr), with k² = 4πG/K and edge R_TF = π/k, whatever its mass. R_TF is quoted in place of K.
- **Per system,** the condensate mass is free. It is the reservoir inventory, not a prediction of supply.
- **Validity.** Thomas–Fermi needs the healing length ħ/(m√(2Kρ)) to be small compared with R_TF. The constituent mass is not fitted; for CR-1's two masses, the report states whether the condition holds.
- **Gravity** is Newtonian for both probes, and light is bent by twice the Newtonian deflection. No modified response is included.

## Data and conventions (identical for both probes)

- **Lenses.** The six SLACS systems of the archive: J0037-0942, J1112+0826, J1204+0358, J1402+6321, J1621+3931 and J1630+4520. For each:
  - the catalogue Einstein radius bSIE;
  - the published Sérsic light profile, spherically Abel-deprojected;
  - its KCWI Vrms bins with their full released covariance (40 bins in all).
- **Geometry.** One set of angular-diameter distances per run, used both for lensing and to convert every angle to kpc.
  - **G1, an adopted comparison.** Flat FLRW with H0 = 70 and Ω_m = 0.3, the convention in which the published stellar masses were derived. It is adopted, not derived.
  - **G2, the co-scaling branch's own history.** CC-1's V = 0 coasting history, with the repository's α:
    - D_M = ln(1+z)/α;
    - D_A = D_M/(1+z);
    - D_ls/D_s = 1 − D_M(z_l)/D_M(z_s).
- **Stellar mass.** Two conventions, each applied identically to lensing and dynamics:
  - **M1, lens-determined.** Stellar mass follows from the exact lens constraint, as in the archive's exact-lens fits. The share f of the deflection supplied by the condensate is fitted per lens.
  - **M2, population masses.** Auger et al. (2009) Chabrier and Salpeter masses, rescaled to each geometry's luminosity distance (M ∝ D_L²). The condensate supplies exactly the remaining deflection. A lens whose stars alone over-deflect has no admissible condensate.
- **Dynamics.** Spherical Jeans modeling with the archive's engine (aperture and PSF kernels, stellar mass following light) and a constant anisotropy β per lens in [−2, 0.45].
- **Likelihood.** χ² with the full covariance over all 40 bins. Lensing is imposed exactly.

## Freedoms and benchmarks

- **Shared:** R_TF, a single value for all six lenses, per geometry and mass convention. It is searched from 0.3 to 300 kpc and then refined.
- **Per lens:** β, and under M1 also the condensate share.
- **Benchmarks,** using the same data, likelihood and code in each geometry:
  - stars only, with β per lens;
  - free NFW, with share, scale and β per lens.
- **Provenance of the benchmarks.** G1's benchmarks are in the archive ([nfw-geometry-results.json](../isotropic-galaxy-transfer/nfw-geometry-results.json)); G2's are computed here with the same code.

## Validation (tolerances declared)

- **V1.** With no host, the profile is sin(x)/x: the edge is at π and the dimensionless mass is π, both to 10⁻⁸.
- **V2.** This code reproduces the archive's G1 stars-only and free-NFW exact-lens χ² for every lens to 10⁻⁴ relative. This establishes that the data, geometry and likelihood are identical.
- **V3.** At every reported fit, the lens constraint closes to 10⁻⁸.

## Declared outputs

- R_TF with its profile likelihood. Per lens: the share, the stellar mass, the condensate mass, β, the condensate edge, and M*/M_pop.
- Total χ² per geometry and convention, against the benchmarks.
- **Rotation curves.** The same R_TF in the Milky Way, using spherically averaged baryon model I for the equilibrium and the archive's baseline speeds.
  - The condensate mass is fitted to the 38 Eilers speeds.
  - Its RMSE is compared with baryons alone, over all 38 bins and over the inner 20.
- **Supply (separate).** For each lens, the condensate mass it needs against an upper bound on what its starlight could supply through the conversion law in 10 Gyr: M_supply ≤ α L R_edge T/c². L is bounded above by M*_Chabrier × (L☉/M☉).
- **Formation history.** Not modeled; labeled.

## Pass rule (declared)

The support law passes in a geometry, under M1, if all three hold:
1. **Admissible.** Every lens has an admissible solution: a share between 0 and 1, a stellar mass within the archive's bounds (10⁷–10¹⁴ Msun), and the condensate's edge inside the grid.
2. **As good as NFW at equal Akaike score.** Its total χ² is at most the free-NFW total in the same geometry plus 10. CR-2 has 13 parameters against NFW's 18.
3. **Milky Way.** With that R_TF, the fitted condensate reaches an RMSE of at most 20 km/s on baseline I, without worsening the inner 20 bins relative to baryons alone.

M2 is reported, not scored. Supply and formation history are reported separately and do not enter the rule.

## Not claimed

- That a supported condensate forms, or is supplied, by any mechanism.
- Any modified gravity response.
- Any cluster result.
- Calibrated significance. The released covariance is used as given, and the benchmarks share the same data.
