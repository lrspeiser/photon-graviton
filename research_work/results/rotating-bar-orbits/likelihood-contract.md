# From a gravity field to the stellar comparison

This document specifies the next observational step. It is not an executed fit. The rotating-bar tests establish a numerical orbit operator; six illustrative trajectories cannot specify the distribution of real stars.

## The quantity to predict

Use one nonnegative orbital population in each declared chemical/tracer class, in a shared candidate gravitational potential. Predict the measured three-component velocity distribution conditional on position and tracer properties. Do not assign every bulge star the circular speed or derive an independent correction in each observed cell.

**Known statistical mixture construction, not a new gravity formula:** for potential parameters theta and nonnegative orbit weights w_k,

\[
f_\theta(\mathbf x,\mathbf v\mid C)=\sum_k w_{k,C} K_{k,\theta}(\mathbf x,\mathbf v),
\qquad w_{k,C}\geq0.
\]

K_k describes time occupancy of an orbit, with a declared finite-resolution representation. Weights and orbit sampling are model choices, not predictions of the companion theory. Orbit mixing duration, coverage, regularization and the resulting effective complexity must be assessed. Equal weights for the six demonstration probes are not a defensible Galactic distribution function.

The observable probability includes the survey's selection S and the measurement response p(d_i|x,v,C):

\[
p(d_i\mid\theta,w,\mathrm{selected})=
\frac{\int p(d_i\mid\mathbf x,\mathbf v,C)\,S(\mathbf x,\mathbf v,C)\,f_\theta(\mathbf x,\mathbf v\mid C)\,d\mathbf x\,d\mathbf v\,dC}
{\int S(\mathbf x,\mathbf v,C)\,f_\theta(\mathbf x,\mathbf v\mid C)\,d\mathbf x\,d\mathbf v\,dC}.
\]

Survey-field dependence and the observed tracer class must be included in S and the conditioning. Conditioning on exact position and properties can cancel a velocity-independent selection factor, but uncertain distances, population mixing and velocity-sensitive quality cuts prevent assuming this cancellation automatically. Document and test any conditional-likelihood approximation.

The existing StarHorse files contain posterior distance summaries, not a standalone measurement likelihood. A hierarchical fit must use appropriate posterior-recycling prior corrections where known, or label its uncertainty integration as approximate. **Route by the actual StarHorse input flag:** the [training input audit](../stellar-orbit-support/report.md) found many distances without recorded PARALLAX use. For parallax-used distances, never count the same Gaia parallax likelihood twice. For parallax-omitted distances, a reliable independent parallax likelihood may instead be added once, after checking association and astrometric quality; disagreement is not repaired by blindly imposing Gaussian correlations. The existing four uncertainty reconstructions are sensitivity inputs, not an exact survey likelihood. Use the recovered provenance sidecar and retain the original sample/holdout history when revising the model.

## Frozen comparison rules

1. Keep the existing 77,927/24,755/26,090 training/validation/test assignments and their documented nonblind history. Shared sky or APOGEE-field systematics remain a limitation.
2. Use the same orbit-population flexibility and selection/error treatment for ordinary matter, the separate standard-halo comparison and companion candidates. Share ordinary-matter nuisance assumptions across the comparisons.
3. Vary or marginalize ordinary-matter normalization, bar orientation/pattern speed and solar-frame parameters within explicitly sourced alternatives. The published baseline parameters are not independent proofs of the ordinary-matter distribution.
4. Restrict companion parameters to a common spatial response/capture rule. The current ring/cap/shell placements and amplitude are synthetic probes and must not become fitted observational claims without a physical specification and frozen parameter count.
5. Fit weights and gravitational parameters on training data. Use validation for model/regularization choices. Freeze the complete procedure before evaluating the test likelihood and calibration by location and tracer class.
6. Report all original and revised variants, interval calibration, region coverage, and failures. Compare predictions of radial, rotational and vertical velocity distributions jointly. A difference in average rotation alone is insufficient to identify an extra gravitational field.

The likelihood implementation, a sufficiently complete orbit library, survey selection, and held-out scores are still outstanding. Even successful stellar prediction would not establish photon origin while the gravity/transfer amplitude degeneracy, propagation/clock consistency, capture, lensing and deferred total source-energy budget remain unresolved.

## Full-bar empirical-response field available for orbit checks

The [full-bar completion](../full-bar-completion/report.md) supplies an additional conservative potential using the unchanged empirical response, with separate three-dimensional-grid, axisymmetric-reference and outer-boundary checks. It differs from the synthetic ring/cap/shell potentials used in the original demonstration trajectories. Use its single potential and derivative force together when testing the empirical-response candidate; do not substitute an algebraic local acceleration multiplier or splice separate radial and vertical corrections. First verify the actual orbit domain, interpolation and rotating-frame Jacobi invariant. The potential has not yet supplied the shared orbital population or selected stellar likelihood specified above, and the photon/capture interpretation remains underived.
