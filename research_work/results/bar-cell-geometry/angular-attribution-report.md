# Angular attribution within the problematic passage

The previous age decomposition localized most of the largest spatial force change to the 61–73 Myr source-age interval. This calculation compares common angular regions within that interval, at the same R=3 kpc source, T=0.1 and target (0.1,0,0).

## Verified decomposition

The exact tetrahedral kernel now optionally returns individual cell contributions before summation. Its default output and gravity law are unchanged; existing cube, exterior-volume, thin-cell and central-symmetry regression checks pass. Three tetrahedra per source-age prism are grouped, the same fourfold target symmetry is applied, and each four-child group of the second mesh is aggregated into its first-mesh parent region.

The decompositions reproduce the archived age-bin forces within 3.74e-17 absolute units per G times total source mass; potentials agree within 2.78e-17. All 560 common parent regions remain represented, and their masses sum to 1/8 of the original source mass. Differences include the previously used capture-weight quadrature at each resolution; this does not isolate trajectory interpolation from angular-weight error.

## Where the difference is concentrated

Rank regions by the magnitude of their force-difference vector. Of 560 parent regions, 18 account for 50%, 66 for 90%, and 179 for 99% of the sum of those magnitudes. The top 66 contain 0.02243 of the original total source mass, or 17.94% of the selected age bin. These are numerical error-attribution statistics, not observed matter or gravity fractions.

The largest four region IDs are 367, 58, 362 and 63. No single region dominates: the largest contributes 5.99% of the summed difference magnitudes. Signed projections along the net difference are retained separately, including negative contributions, because opposing directions can cancel. The summary JSON retains every region, both mass weights, force differences and rankings.

## Next decision

The result supplies a finite set of angular targets for actual trajectory refinement during this passage. It does not validate the other cells, justify dropping them, or solve the other failing target positions. A refinement must retain the full source, preserve symmetry and mass, and compare the total force again. The existing 14 spatial force failures remain unresolved.

The formulas used are known linear superposition, vector norms and projection identities; no capture law, coupling, force law or physical source amplitude changed. No stellar/lensing fit, energy closure or full-model validation is established. All nine goals remain open.

Run `angular_force_attribution.py 1`, `angular_force_attribution.py 2`, and `summarize_angular_attribution.py` from this directory. Both decomposition runs finished successfully; no process from this audit remains active.
