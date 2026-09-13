# Different companion loading for different matter components

**Component-dependent loading does not reproduce the fixed three-dimensional source target.** With seven freely fitted nonnegative component coefficients, the inner sampled region still has a minimum worst source-density mismatch of 42.47%. This is conditional model geometry, not an observed velocity error or rejection of all companion storage.

## Candidate and provenance

Optional postulate, not claimed novel: rho_extra(x)=sum_j eta_j rho_b,j(x). Each eta_j is constant across its matter component and nonnegative. Components are the existing bar, nuclear component, two stellar disks, two gas disks and softened central source. Thin, uniform illumination with different receiver efficiencies could motivate this shape, but opacity, exposure and storage are not derived. Shielding from the preceding study is not combined with these coefficients: component-dependent absorption would require recomputing the attenuation consistently.

We permit unbounded positive loadings to give the shape hypothesis generous freedom. The full-grid optimizer uses a softened-center coefficient around 3.25e7; this is a diagnostic mathematical optimum, not a plausible energy supply or fitted central mass measurement. Zero coefficients are optimizer choices, not physical evidence that components cannot capture companions.

Known linear programming minimizes the maximum relative source-density error over the selected subset. The primal residual and dual objective agree within 1e-7, verifying the finite-grid optimum for this specified linear model. Coefficients and geometry predictions are retained. No new gravity formula or microscopic cross section follows.

| Target resolution | Fitted subset | Best worst source-density mismatch | Worst mismatch on excluded locations |
|---|---|---:|---:|
| fine | full | 99.352% | — |
| fine | inner | 42.474% | 99.666% |
| fine | midplane | 11.576% | 99.656% |
| finer | full | 99.350% | — |
| finer | inner | 42.473% | 99.655% |
| finer | midplane | 11.707% | 99.659% |

The inner subset is R<=8 kpc and |z|<=1 kpc (120 points), the midplane has 40 points, and the full grid has 240. Each subset is fitted separately to diagnose geometry; coefficients from different rows do not constitute one model. Excluded locations have already been examined in earlier work and are not a blind validation sample. Two target reconstruction resolutions give closely similar inner and full-grid mismatch bounds; this does not propagate ordinary-matter or observational uncertainty.

## What the result changes

Uniform loading previously failed; allowing separate receivers also leaves substantial mismatch. A midplane-only fit can conceal severe off-plane errors, strengthening the need for three-dimensional tests. These are source-density discrepancies; accelerations and lensing require spatial integration, so the percentages must not be reported as stellar-speed errors.

The test does not include spatially varying histories within a component, combined receiver-dependent shielding, movement after capture, or independently supported field energy. Those remain explicit physical options. It also does not establish that the empirical target itself is the right one. Progress toward a predictive mechanism now requires a consistent transport/storage equation or forward observable calculation, not assigning more arbitrary coefficients to individual target locations. Photon supply, timing and joint motion/lensing remain open. All six research objectives remain open.
