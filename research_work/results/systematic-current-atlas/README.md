# JR-7: systematic current atlas and observational characteristics

**21 September 2026 (America/Los_Angeles). Executed exploratory research; not a validated gravity theory.**

## Theory first

Ordinary matter and radiation may generate a companion state whose spatial organization, currents and history supply a connected gravitational response. Gas can promote conversion or return, and its overlap with illumination can change the resulting state. The present positive two-channel transport model is one implementation, not a restriction of the wider branch portfolio to diffusion, local Markovian laws or one metric. No observation, R10 coefficient, or earlier result was altered.

## Completed scope

- 576 primary 3D finite-volume stationary cases in a full factorial design.
- Eight representative cases at three resolutions (24 solves), plus seven tilt, source-offset, two-source and mirror controls.
- Seven fixed statistical predictors of the existing R10 galaxy-level residual across 149 SPARC objects; no refit of the gravitational formula.
- Three real MaNGA DR17 spatial-map products downloaded, with projected stellar/gas velocity and light-overlap features extracted.
- A descriptive census of the existing 12 X-COP cluster profile sets.
- 2,903 local numerical/readback checks and 72 cluster-value checks, zero failures. These verify arithmetic and implementation, not astronomical validity.

No new joint real-galaxy/cluster gravitational fit was performed. The map prototype is not matched to the SPARC sample, and cluster thermodynamic profiles are not full 3D current maps.

## Spatial construction

In the frame of fixed gas/source patterns:

    P_t + Omega_P P_phi = D_P lap(P) - k_plus P + k_minus C - lambda_P P + J
    C_t + Omega_C C_phi = div(D_C grad C) + k_plus P - k_minus C - lambda_C C.

P and C are coarse-grained energy densities. Positive rates are

    h=rho_g/(rho_g+0.15)
    k_plus=0.5 exp(a h), k_minus=0.5 exp(b h^2).

The four (a,b) pairs are (0,0), (2,0), (0,4), (2,4): constant, enhanced forward, enhanced return, and mixed. The factorial additionally varies source/gas axis angle 0/45/90 degrees, three companion angular drifts, smooth/structured gas, two gas thicknesses, two companion escape rates, and two vertical/horizontal transport ratios. Gas inventory and luminosity are equal across geometries on each grid. The discretized reference gas profile also changes slightly on refinement; absolute grid comparisons include that quadrature effect.

Exchange cancels locally. At stationarity, source power equals outgoing radiation plus outgoing companion power. Prescribed rotation and diffusion are not a microscopic entrainment derivation, real photon propagation at c, gas dynamics, recoil closure, or a complete gravitational binding-energy ledger. Source/gas phase locking is a declared preparation, not a demonstrated permanent physical alignment.

## Direction now changes the mean

JR-6's axisymmetric rates left the angular mean independent of angular currents. Structured gas in JR-7 couples the angular pattern to mean exchange. The exact identity is

    mean(Q)=mean(k_plus)mean(P)-mean(k_minus)mean(C)
            +Cov(k_plus,P)-Cov(k_minus,C).

The covariances are computed from fields, not fitted independently for each galaxy.

Across the primary scan, aligning gas/source axes rather than making them perpendicular changes integrated companion energy by:

| Rate family | Range across 24 structured matched pairs | Signs |
|---|---:|---|
| Constant | 0 | No phase effect |
| Forward-enhanced | +2.04% to +2.71% | All positive |
| Return-enhanced | -2.12% to -0.46% | All negative |
| Mixed | -1.35% to +1.08% | 12 positive, 12 negative |

Selected matched pairs were refined before interpreting their magnitudes. The following is the finest-grid result, aligned relative to perpendicular:

| Selected case | Companion energy | Mean inward force, companion only | Mean inward force, both channels | Mean face-on deflection, both channels |
|---|---:|---:|---:|---:|
| Forward-enhanced | +2.535% | +2.676% | +1.684% | +1.408% |
| Return-enhanced | -2.119% | -3.808% | -1.058% | -0.811% |
| Mixed: thick gas, strong vertical transport | +1.138% | +0.662% | +0.276% | +0.358% |
| Mixed: thin gas, weak vertical transport | -1.042% | -1.855% | -0.763% | -0.510% |

The two mixed configurations also differ in companion drift; this is not one-factor attribution to thickness. Every non-phase parameter is fixed within each aligned/perpendicular pair. Source energy is not itself gravitational force, and a percentage force change is not the same percentage circular-speed change.

Force and lensing use one fixed softened 3D potential (G=c=1, softening 0.35), evaluated at radius/impact 2.5. They are conditional finite-domain toy readouts, not measured galaxy effects. Counting both energy channels retains the sign changes, so this is not only counting a relabeled channel while ignoring the rest of the energy.

### Numerical boundaries

All primary states are nonnegative. Maximum stationary energy error is 1.78e-15. Axisymmetric and constant-rate controls have zero phase effect; reversing currents and reflecting the pattern mirrors the solution. A separate source-removal evolution conserves active plus escaped energy while active energy decreases.

From the medium to fine selected grids, integrated companion energy changes by at most 0.262%, mean radial force by 0.718%, and mean face-on lensing by 0.518%. Matched contrast changes are smaller and preserve both mixed signs. This is not an infinite-domain or full-continuum certificate: all grids retain the same reflecting domain. Seven geometry controls include a two-source configuration with displaced broad gas components; it is not an observed cluster, a shock simulation, or a dynamically derived merger offset.

## Static galaxy characteristics: modest, uncertain signal

The target is one signed residual per galaxy, mean ln(v_R10/v_observed), not a new force or observation. Predictors are surface brightness, inferred gas fraction, stellar size, gas-to-light extent, luminosity, morphology, light concentration and mean gas surface density. One candidate additionally includes observational nuisance metadata. Observed rotation amplitude is not a feature.

Seven predictors were fixed before training. Training uses 89 objects, selection 29, and the final comparison 31; all objects were historically exposed. Training-only imputation and five-fold predictions are retained. The selected regularized pair-interaction predictor gives comparison residual RMSE 0.11148 versus the constant-control 0.12693, or 22.86% less squared error. But the paired-bootstrap 95% interval for selected-minus-control MSE is [-0.011415,+0.002882], including no improvement. Residual-sign accuracy is 58.1% for both selected and constant predictors.

The largest single training association is the quality flag (Spearman 0.267); gas fraction is only 0.096. Correlated descriptors and prior source-model choices complicate attribution. These associations do not establish a new physical law, a 22.86% reduction in actual rotation-curve errors, or the values of unmeasured currents.

## Real projected maps, not invented 3D orientations

The three MaNGA IDs below are exact examples in the official DR17 access documentation, not a representative or counterrotation-selected sample. Weighted planar stellar and Halpha velocity gradients use a common valid mask and Halpha signal-to-noise at least five.

| ID | Common valid pixels | Projected gradient mismatch | Raw continuum/Halpha fractional covariance |
|---|---:|---:|---:|
| 7443-12703 | 2533 | 14.17 degrees | 4.129 |
| 8082-12704 | 618 | 12.79 degrees | 0.215 |
| 8081-1901 | 404 | 11.25 degrees | 1.810 |

These are not full kinematic-position-angle analyses or 3D spin measurements. Velocity-plane residuals, beam effects, inclination and bin covariance remain relevant. Halpha brightness is not total gas density; continuum is not bolometric source power. No companion rate is calibrated from these three examples.

Removing a coarse circular-annulus brightness trend reduces the three fractional covariances to 1.040, 0.077, and 0.181. This demonstrates why raw correlation cannot be equated directly with conversion overlap; it is not itself a deprojection. Original FITS, extracted arrays, masks and hashes are in the delivered archive. An independent normal-equation fit reproduces the original weighted least-squares angles.

## Cluster characteristic census

The pinned X-COP extract provides 12 clusters, 778 electron-density shells, 150 X-ray temperature points, 150 X-ray pressure points, 122 SZ pressure points, and seven cumulative stellar profiles. Missing stellar profiles remain missing. Median reported cluster temperatures span 4.06-8.76 keV; the density ratio in the nearest shells to 100 and 500 kpc spans 2.75-10.73. Exact shell midpoints are saved, not replaced with unperformed interpolation.

These are deprojected/spectrally modeled release products retaining the release's distances and normalizations, not assumption-free raw data. Hydrostatic and NFW/Einasto mass estimates are excluded from companion-source inference. The extract has no full gas velocity direction, stellar velocity field, angular light/gas overlap, or matched lens-image likelihood. No cluster gravity fit was performed.

The initial helper's substring search for 'star' missed the actual key 'stellar_mass'. A separate explicit-key readback corrected that metadata to seven stellar profiles. Original inputs and the original summary are retained, with all 72 numerical value checks passing. It is not a changed stellar-mass measurement.

## How the full exploration should proceed

Use the observed source geometry to constrain a common family of positive exchange laws, anisotropic transport, relative-current response, nonlinear capacity, and memory. Intrinsic variables should include gas phase/temperature, clumping, source spectrum, relative velocities and velocity-dispersion tensors, overlap, thickness, and physical crossing/exchange/escape ratios. Viewing inclination belongs in the observation model, not in the intrinsic law.

For example, a positive rate can be an exponential of a small basis of intrinsic scalar descriptors; an anisotropic diffusion tensor can be LL^T. These are controlled parameterizations, not derivations. Keep zero-effect limits, reciprocal energy accounting, whole-object comparison roles, uncertainty and model complexity visible. Do not create a free invisible source for each object or tune missing orientations from the residuals being explained.

Compare complete stellar/gas velocity maps and lens-image patterns. Spectral cubes have two sky dimensions plus wavelength, not measured depth. Stars, gas and the companion need not share velocities. Individual stellar or planetary spin enters only through a specified weighted interaction; those spins were not measured here. Magnetic effects on plasma do not automatically imply a direct neutral-companion coupling.

Cluster comparisons should combine thermodynamics, stellar distributions and actual image constraints for the same targets. Hubble Frontier Fields provides multiple model-dependent reconstructions for six clusters, not an automatic matched extension of X-COP. Prefer the underlying image constraints; do not input reconstructed dark-halo mass as the proposed companion source. A static two-Gaussian analogy does not reproduce a merger's history.

This atlas does not exhaust all formulas, thermal phases, dust, magnetic couplings, spin, nonlinear self-consistency or memory histories. The new concrete result is a sign-changing connection from overlap and transport to mean response. The new observational capability is measured-map ingestion. The source-to-gravity coupling, energy supply, nested planetary behavior and joint astrophysical prediction remain to be solved.

## Provenance and reproduction

Baseline: f55c07df0f2da42b6a2d4d49c0b73975ad17d4ec. Protocol: e094f685f3be8badccf425afbb23e4db95b52d04. Local scientific environment: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, scikit-learn 1.8.0.

MaNGA acquisition run 35686945515 completed; artifact SHA-256 8dd8f8e35c1859d34fb8c102311b7b19d7225083d4988b7dda6df54351c16043. Cluster census run 35687241997 completed; artifact SHA-256 47a1a89b8448551a08cd38b3bdfd1e78d16dbc62f614190817c204ea61bbc214. These workflows did not run the local spatial atlas or demographic regression.

The delivered complete archive contains all inputs, code, primary/fine arrays, feature outputs, audits and detailed reports. The core atlas.py and refine.py are committed. The attempted characteristics.py write was blocked, so that exact executed code and local audit utilities are provided in the archive rather than claimed present here.

For the primary atlas use `python atlas.py --output fresh --skip-refinement`. Write the preselected representative IDs [168,264,316,412,445,456,541,552] to fresh/extended-selection-before-refinement.json as an ids array, then use `python refine.py --results fresh`. The delivered select_refinements.py reconstructs that deterministic selection. Use fresh paths; earlier evidence is not overwritten.

## Primary data sources

- SDSS DR17: https://www.sdss4.org/dr17/manga/manga-data/data-access/
- Co/counterrotation comparison: Beom et al. 2024, https://arxiv.org/html/2410.06256v1 . Their matched observational groups are not evidence for this companion theory and were not the three demonstration objects.
- THINGS: https://arxiv.org/abs/0810.2125 ; PHANGS: https://sites.google.com/view/phangs/home . Complementary maps, not newly ingested here.
- X-COP release: https://dominiqueeckert.wixsite.com/xcop/data ; original provenance is retained in cl2-inputs-xcop-profiles.json.
- Hubble Frontier Fields: https://archive.stsci.edu/prepds/frontier/lensmodels/ . No new HFF image/model array was ingested.
