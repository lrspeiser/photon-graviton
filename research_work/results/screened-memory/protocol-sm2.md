# SM-2: finish the source-based benchmark; do not relabel AQUAL as a new theory

Declared 2026-09-21 against main `abd026f324c6bf90f71a68f0b0d6a4f8b889c6c3`, before the SM-2 galaxy runs. This is exposed-data validation of an effective candidate, not a blind discovery experiment.

## Novelty boundary

The static SM-1 PDE is exactly in the Bekenstein-Milgrom AQUAL class, with a particular interpolation function. The unscreened PM law and rationalized mu are exactly the published Bekenstein toy function printed by Famaey and Binney (2005), equation 5. These antecedents were already documented in `research_plan/prior-art/audit.md`; the SM-1 summary should have made this equivalence more prominent. A chosen screening denominator, a variational energy, and a dynamical scalar extension do not establish a new fundamental theory. Static agreement with a known theory does not by itself preclude a new microscopic derivation, but none has been supplied for SM-1. This experiment retains SM-1 as a benchmark rather than promoting it as unique photon-graviton physics.

## Fixed source and equations

Use all 149 names and the frozen 89/29/31 split loaded by capture-to-orbit `inputs.sparc_galaxies`. Use all positive-radius raw rotmod rows, exactly as SM-1. Build each source once from the FULL source rows with PM-2A `build_source` and `sigma_of`; retain the archived endpoint-anchored stellar continuation, gas prescription, bulge, disk mass-to-light 0.5 and bulge 0.7. No velocity residual enters the source. Exponential vertical scale is 0.1*rd for both stellar and gas disk, fixed before runs. These are declared inferred density models, not assumption-free observations. Retain hashes and mass totals.

Solve four total-potential equations on the SAME cell masses: Newtonian; unscreened Bekenstein-root AQUAL at a0=6.54e-11 m/s^2; frozen SM-1 screened AQUAL at that same a0; simple-mu AQUAL comparison at the previously selected 8.563e-11 m/s^2. No new fit and no coefficient/crossover/exponent scan. The known AQUAL and Bekenstein forms are controls, not claimed discoveries.

Inner radius = .002*min(rd, smallest positive observed radius). Outer radius = max(80*rd, 20*h_gas, 10*largest observed radius, 80*sqrt(G*M_total/a0)). The grid is geometric radial with the existing angular clustering. Primary grids are (192,32) and (384,64); a third (768,128) is NOT automatically substituted after an unfavorable score. For EVERY source, the screened boundary control doubles outer radius, uses enough radial cells to preserve or improve the fine grid log spacing, and retains 64 angular cells. Source quadrature uses 8 radial and 12 vertical nodes; record cell-total discrepancy against independent PM-2A analytic total. Use the spherical unresolved-core integral rather than the inherited thin-cylinder approximation. Radial forces are read with PM-2A's face-gradient cubic spline, without extrapolation.

## Accuracy gates and reporting

Require every solve to meet the inherited residual <1e-8 and iteration-change <1e-11 gates. Require mass discrepancy <1e-3. For EACH equation and each galaxy record max pointwise force refinement change, dividing by max(abs(fine force at the point), 0.01*max(abs(fine force))). The common matched numerical certification requires every equation below 2% refinement and screened boundary change below 0.5%. Failure is an unresolved numerical prediction, not a rejected force law. Preserve every failure and all scores; report all-source provisional and common-certified-subset scores separately. Do not silently replace the full sample with passing galaxies.

Positive inward radial acceleration predicts v=sqrt(R*g). Nonpositive force gives no real circular orbit at that row: record it explicitly, do not hide it by clipping to a plausible speed, and do not assign a complete per-galaxy RMSE when a model lacks a prediction. The observational error bars are retained, but the primary metric remains the unweighted mean of per-galaxy velocity RMSE in km/s. Also report median and pooled velocity RMSE. Report the reconstructed Newtonian baseline versus the original SPARC baryonic-force baseline, so source reconstruction differences cannot masquerade as a theory improvement.

Keep per-row model predictions, source fingerprints, split counts, all refinement/boundary results, source and input hashes, executed git SHA, and software versions in the result archive. Paired bootstrap differences use 10000 galaxy-level resamples and a fixed seed (210926); they quantify variation over these exposed galaxies under the fixed source model, not astrophysical systematic errors or new blind confirmation.

## Independent dynamics/novelty calculation

The SM-1 state (U,Ut) is a local Markovian scalar system. A transient after changing the source is not automatically a new history-dependent law or an irreversible deposit. Derive its linear response around a constant nonzero field without changing its kinetic coefficient, and distinguish static susceptibility, propagation speed and relaxation time. A physical memory kernel must eventually be calculated from an explicit microscopic coupling and initial state. Do not graft a fitted tau onto the galaxy law and call it derived.

The source-based galaxy benchmark does not certify Solar-system external-field quadrupoles, lensing, tensor waves, scalar radiation, microscopic derivation, long-lived deposits or cosmology. These remain separately open. In particular do not use a QUMOND quadrupole as if it were the prediction of this AQUAL PDE.

## Primary references to verify directly

Bekenstein & Milgrom (1984), ApJ 286, 7-14, equations 2b-3.
Famaey & Binney (2005), MNRAS 363, 603-608, arXiv:astro-ph/0506723, equation 5.
The original PM-2A source and field provenance in this repository remains unchanged.
