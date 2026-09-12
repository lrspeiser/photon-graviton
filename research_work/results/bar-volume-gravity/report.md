# Gravity from mass-preserving trajectory volumes

The new solver integrates the gravity of finite volume cells instead of treating a few sampled streams as point contributions. It retains all source mass and passes independent field checks. It does not yet produce converged Galactic forces: one final time-resolution comparison and fifteen direction-resolution comparisons still fail. No new observation fit or successful redshift mechanism is claimed.

## Physical assumptions and formula provenance

The continuous source law, initial50km/s radial motion, two separate launch spheres and fixed ordinary bar are unchanged. Each source is normalized per unit total injected rest mass. Physical injection rates, their energy funding, other Galactic components and deposited self-gravity remain absent. This is still the project cold-particle interpretation, not a derived graviton interaction. The retained bar normalization comes from a model-dependent reconstruction; omitting its upstream halo component does not make it an independent ordinary-mass measurement.

The gravitational volume-to-surface identities are established mathematics, described for example by [Pearl and Hitt, MNRAS](https://academic.oup.com/mnras/article/492/1/420/5673493). We use analytic planar-face integration rather than their approximate face quadrature. The source transport and cell discretization here are project numerical choices; the gravity equations are not claimed novel.

For a cell of mass m, volume V and constant density rho=m/V, with outward face normal n_f, face distance d_f=n_f dot(x_f-p), and I_f=integral_face dS/|x-p|,

\[\Phi(p)=-\frac{G\rho}{2}\sum_f d_f I_f,\qquad \boldsymbol a(p)=-G\rho\sum_f\boldsymbol n_f I_f.\]

Output coefficients use G=1 and total source mass1: potential has units1/kpc and acceleration1/kpc^2. Multiply by G times the physical total injected mass when that mass is established.

The face integral is evaluated using edge logarithms and a signed solid angle. Code and sign conventions are in tetra.py. No softened force or excluded central volume is introduced. A degenerate cell causes an error rather than loss of its mass.

## How the moving source is represented

The source6/source8 angular nodes give72/128 full-sphere directions per launch sphere after symmetry reconstruction. Their convex-hull triangles have140/252 faces. Each triangle receives its solid-angle fraction times the mean vertex capture factor; this is a different numerical quadrature of the same continuous law from the previous tensor Gauss weights. Its total source normalization is recorded, not silently equated with the earlier approximation.

Adjacent trajectory ages form triangular prisms. Each prism is divided into three tetrahedra, each receiving one third of the prism source mass. Density is mass divided by the mapped cell volume. Folded or overlapping cells contribute positive mass independently. Shared-face diagonals are consistent, and an explicit symmetry average prevents arbitrary hull diagonals from breaking plane reflection or half-turn symmetry.

Straight-sided cells approximate curved trajectories and the angular source map. This can introduce substantial interpolation bias on a coarse grid. Mass preservation alone does not establish a correct density or force. It is a numerical approximation to refine, not a new physical smoothing length.

## Verification

The cube interior potential and force agree with the independent stream benchmark to8.14e-15 and1.05e-14 relative error. A regular tetrahedron has zero center force; exterior fields agree with independent positive volume quadrature. Flattened-tetrahedron checks through a height factor1e-6 remain below4.8e-9 relative error.

All100 cached trajectories pass Jacobi checks. Full-versus-thinned Hermite interpolation changes position by at most1.99e-08kpc and rotating velocity by0.000976km/s. Cache files and field inputs have recorded hashes. These are numerical checks of this fixed model, not energy-supply closure.

For180 sampled actual-cell/point combinations,60-digit arithmetic changes potential by at most3.15e-12 and force by1.06e-11 per G total source mass. This samples the thinnest retained cells; it is not a global error bound on every cell.

Across48 source/age configurations, maximum normalized mass discrepancy is1.67e-15. Initial64/128 and128/256 failures are preserved. At256/512, all60 potential checks pass and59/60 force checks pass. The remaining force failure is source8, R3kpc, T0.25, point(0.1,0,0), with17.8 percent change. At matched512 layers, source6/source8 gives8/30 failed potential gates and15/30 failed force gates. Thresholds remain2 percent potential and5 percent vector force, with the documented0.01 floors; these are not observational uncertainties.

## Direction refinement at the final epoch

| Launch R (kpc) | Position (kpc) | Potential change (%) | Vector-force change (%) |
|---:|---|---:|---:|
| 1 | (0.1, 0, 0) | 6.271 | 43.481 |
| 1 | (1, 0, 0) | 1.357 | 3.888 |
| 1 | (1, 0, 1) | 0.096 | 1.068 |
| 1 | (0, 1, 1) | 0.416 | 1.401 |
| 1 | (3, 0, 0) | 0.165 | 0.488 |
| 3 | (0.1, 0, 0) | 5.305 | 40.555 |
| 3 | (1, 0, 0) | 2.052 | 10.214 |
| 3 | (1, 0, 1) | 1.848 | 9.268 |
| 3 | (0, 1, 1) | 2.922 | 13.473 |
| 3 | (3, 0, 0) | 0.592 | 5.351 |

## What still needs work

The source map and cell geometry need angular refinement or local refinement where they bend and fold. The late central time-resolution failure also remains open. These errors cannot be removed by changing the physical source amplitude. Do not equate the new coefficients with a measured stellar force or claim that the volume method has already solved source convergence.

After numerical convergence, the model must still fund injection and evolve deposited gravity and the ordinary-matter response together. Lensing requires a specified metric/stress response, and the photon/time redshift branch remains unresolved. All nine project goals remain active.

Reproduce: verify.py; prepare.py6 and prepare.py8 (numbers are separate arguments); cache_check.py; volumes.py with each source number and64,128,256,512 layers; roundoff.py; export.py. Generated trajectory caches live under research_work/data-cache/bar-volumes and are hash-listed in prepared6.json/prepared8.json. All mass, field and interpolation assumptions are stated in protocol.md.
