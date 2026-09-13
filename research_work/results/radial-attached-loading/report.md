# Can arbitrary radial loading reproduce the three-dimensional target?

**Even an independently adjustable loading at every sampled radius cannot reproduce the existing target source everywhere.** Within R<=8 kpc and |z|<=1 kpc, some source-density mismatch must be at least 43.13% under the specified relative-error measure. This is a comparison of model shapes, not a measured discrepancy in stellar motion or a rejection of the broader companion hypothesis.

## Optional closure tested

Project postulate, not claimed novel: rho_extra(R,z,phi) = eta(R) rho_b(R,z,phi), with nonnegative eta allowed to vary freely at every sampled radius. This allows any radial loading pattern, including enhanced outskirts, but assumes identical loading per unit ordinary matter across height and angle at a given radius. It also assumes ordinary Poisson sourcing, no redistribution away from receiving matter, and no separate gravitational response kernel.

The target is the previously reconstructed empirical extra-field source. It is not proven to fit all stellar and lensing observations. Ordinary-matter uncertainty remains unpropagated. No new stars or observational holdouts were opened.

## Known minimax algebra

At each radius let q=rho_target/rho_b, with extrema q_min and q_max over sampled heights and angles. Balancing opposite relative errors gives eta_opt=2 q_min q_max/(q_min+q_max) and minimum worst relative error (q_max-q_min)/(q_max+q_min). The global bound is the largest radial-group bound because eta is free for every radius. For multiplicative errors, the optimal worst factor is sqrt(q_max/q_min); this uses a different loss and generally a different optimal eta. These are known algebraic results applied conditionally, not new laws. Independent scalar minimization verifies every fractional-error bound.

| Subset | Points | Unavoidable worst relative error | Unavoidable worst multiplicative factor |
|---|---:|---:|---:|
| full | 240 | 99.85% | 36.5356 |
| inner_R_le_8_z_le_1 | 120 | 43.13% | 1.5865 |
| midplane | 40 | 10.47% | 1.1108 |

The inner worst group occurs at R=0.5 kpc. Required loading spans 0.198146 to 0.498703 at that radius. Exact locations and every group are retained in results.json.

## Consequence for capture modeling

Compared with the previous constant-loading inner bound of 71.54%, arbitrary radial freedom improves the fit but leaves at least 43.13% mismatch in some inner sampled source density. Even the midplane has nonzero residual variation with bar angle. These are not percentage errors in acceleration or velocity: gravity integrates source density nonlocally.

To retain this target under attached storage, capture/exposure must vary with height or angle as well as radius, or deposits must redistribute, or the sourcing equation must change. We have not selected which explanation works. A local binding-depth rule can depend on height and angle and is not ruled out by this test. Setting eta separately to the target ratio at every point would be reconstruction, not prediction.

Next physical work must calculate that spatial dependence from shared capture and transport rules and then forward-predict motion and lensing. This bound does not determine photon supply, universe age or deposit stability. All six research objectives remain open.
