# What source would produce our full-bar extra gravity?

**All 240 sampled locations have positive equivalent source density**, both in the fine/finer reconstructed potentials and in a direct evaluation of the requested field-equation source. This is a limited compatibility result for a positive Newtonian source. It does not show that photons produce the source, establish positivity everywhere, or validate its detailed density distribution.

This audit applies to the current conservative three-dimensional bar field used in the orbit integrations. The earlier source audit concerned different potential constructions; its positivity results could not automatically be transferred to this field.

## Equations and their status

The known Newtonian Poisson relation gives an equivalent source:

`rho_extra = Laplacian(Phi_extra)/(4*pi*G) = -div(a_extra)/(4*pi*G)`.

This is the density that would reproduce the extra potential if ordinary Newtonian sourcing applied. It is not automatically deposited energy divided by c squared. That identification would be an additional assumption about how companion energy gravitates; a spacetime-response interpretation needs its own field equation and stress/energy accounting.

The implemented field uses the known [QUMOND construction (Milgrom 2010)](https://arxiv.org/abs/0911.5464) with the project's empirical power-law coefficients. Here Q denotes only the extra response; ordinary gravity is added separately:

`Q = A (|grad(Phi_b)|/a_star)^(p-1) grad(Phi_b)`
`Laplacian(Phi_extra) = div(Q)`.

The direct div(Q) calculation tests the intended source, while the Laplacian of the cached potential tests the source actually represented by our numerical solver. These need not agree perfectly at finite resolution. Neither expression derives a photon-conversion rate, capture cross section or occupied storage capacity.

For each real spherical harmonic of degree l, the known Laplacian identity used here is:

`Laplacian[f_lm(log r) Y_lm] = [f_lm'' + f_lm' - l(l+1) f_lm] Y_lm / r^2`.

The primes refer to log-radius derivatives. This evaluates the second derivatives of the same spline potential used for the force. An independent finite difference of its acceleration checks the sign and normalization. These are established mathematical identities, not claimed novel formulas.

## Grid and results

The grid uses R = 0.5, 1, 2, 3, 5, 8, 12 and 20 kpc; z = 0, 0.1, 0.5, 1, 2 and 4 kpc; and five bar angles from 0 to pi/2. These are numerical probes, not observed stars. The retained even harmonics impose reflection symmetries; this grid cannot test real-galaxy asymmetries. No dark halo, new parameter adjustment or holdout observation is introduced.

The reconstructed equivalent density ranges from 0.001162 to 0.792340 solar masses per cubic parsec. The directly requested source ranges from 0.001159 to 0.792336. These are conditional model outputs, not measurements of companions.

| Numerical comparison | Largest fractional difference | Median fractional difference |
|---|---:|---:|
| fine to finer source | 5.41425% | 0.05625% |
| analytic laplacian vs acceleration divergence | 0.00136% | 0.00011% |
| requested divQ vs reconstructed laplacian | 1.09069% | 0.10885% |
| requested divQ step change | 0.00414% | 0.00033% |

The largest fine/finer density change is 5.414% at R=1 kpc, z=4 kpc, on the bar axis. The largest intended-versus-reconstructed source discrepancy is 1.091% at R=12 kpc, z=0.5 kpc, on the bar axis. The latter comparison uses the same cached ordinary-matter field and is not an independent assessment of its mass accuracy. The axisymmetric reference, ordinary-matter model and outer boundary were not further refined in this audit.

Finite differences use steps 0.004, 0.002 and 0.001 kpc. The final requested-source step change is at most 0.00414%; the cached-potential Laplacian and acceleration-divergence calculations agree to 0.00136%. Thus the larger source discrepancies are not removed merely by reducing this finite-difference step. More angular/radial field resolution and ordinary-matter uncertainty would be needed before claiming percent-level source accuracy throughout the grid.

## What this changes—and what remains to be derived

The sign check does not reveal a negative-source obstruction on this grid. It therefore leaves open the conditional possibility that a positive source produces the proposed extra field. However, a positive source is only a target: we still need a common transport/capture/support law that predicts this spatial distribution, rather than assigning the required density after fitting the gravity.

The previously demonstrated amplitude degeneracy also remains: changing the photon-loss scale and compensating with a free gravity amplitude preserves the rotation predictions. Positive reconstructed density does not remove that freedom. A photon-origin claim requires an independently constrained response or capture law.

Practical next steps are to keep the density-resolution limitation visible, specify the physical response that maps deposited energy to the scalar potential, and test whether that one law predicts radial, off-plane and lensing behavior. The total cosmic photon-supply budget remains deferred, not passed. No new observational agreement or global source admissibility is claimed here.

## Reproduction

Run `run.py`, then `report.py`. Results retain every sampled location, all derivative steps, both potential resolutions and the direct source. Input hashes are checked before and after calculation and again during reporting. A quadratic-potential control verifies divergence sign and normalization. No catalog, empirical coefficient, orbit cache or sample role is changed.
