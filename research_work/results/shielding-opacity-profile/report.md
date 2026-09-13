# Shared-opacity calibration of the shielding model

We vary one common mass absorption coefficient over 1e-13 to 1e-7 kpc^2/Msun and analytically optimize one exposure normalization. All original 240 target probes remain exposed development inputs. The full and inner subsets are calibrated separately for diagnosis, not combined into a single galaxy solution.

## Assumptions and formula provenance

The optional postulate remains rho_extra=C rho_b T(x;kappa), with fixed ordinary-matter geometry, isotropic boundary illumination at 30 kpc, straight-ray absorption beta=kappa rho_b, and local attached storage. No new capture dependence, transport history, focusing or self-gravity feedback is introduced.

Known absorption integration gives log T=logsumexp(log angular_weight-kappa*column). Known minimax algebra optimizes C using q/T, where q is the target source per ordinary mass. We minimize the range of log(q/T); this avoids numerical saturation of the equivalent relative-error bound tanh(range/2). These are standard mathematical constructions, not newly derived physical laws.

At each resolution a 121-point scan locates candidate minima; bounded scalar optimization refines every sampled local minimum. Both search endpoints are retained. This does not prove global optimality between all scan points or outside the declared range. Original pair-calibrated transmissions reproduce to relative 1e-12. Ray columns are cached with hashes.

| Subset | Angular mu nodes | Fitted kappa (kpc^2/Msun) | Lowest found worst source-density mismatch |
|---|---:|---:|---:|
| full | 16 | 3.6104478e-09 | 99.518% |
| inner | 16 | 7.6220066e-10 | 46.543% |
| full | 32 | 3.6202802e-09 | 99.524% |
| inner | 32 | 7.5871279e-10 | 46.541% |

These are source-density shape errors for a conditional target, not observed velocity errors, confidence levels or evidence against every companion model. The normalization absorbs unspecified incoming intensity and exposure duration; there is no supply prediction or fixed universe age. A remaining mismatch indicates this particular geometry/capture/storage combination does not reproduce the target within the explored parameter range. Target and ordinary-matter uncertainties also remain open.

All scans, optimizer outcomes and resolution results are retained. A parameter fit to this exposed target cannot count as an independent prediction. Motion, lensing, energy supply and receiver dynamics still require a common physical model. All six scientific goals remain open.
