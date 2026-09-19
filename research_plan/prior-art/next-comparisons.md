# Prior-art comparison programme for the local research agent

This is a proposed work queue, not a report of completed experiments. Source IDs refer to [sources.json](sources.json). Preserve all current numerical failures and historical claims in their existing archives. Each new physical benchmark needs its own declaration before its verification runs; an exploratory scan must remain labeled exploratory.

## 1. Apply the already-settled attribution in new writing

Use “the Bekenstein toy interpolation reproduced in Famaey and Binney (2005), equation 5” for the PM function. Explain the source's spherical/regime qualifications. Label its fitted coefficient as calibrated here. Identify the AQUAL architecture and standard functional construction. Do not change frozen `fields.py` or old result hashes merely to edit its prose; use the new provenance document and apply source labels in the next manuscript revision.

**Done when:** every new claim about PM's force law distinguishes the known function, our implementation, our calibration and the unproven connection to RUT. No implication of full TeVeS equivalence is introduced.

## 2. Reproduce the nearest active-lattice control [TDR2020]

Implement the published constrained lattice with its stated kernel, drag and memory normalization. Reproduce a reported linear eigenvalue/instability boundary before changing those ingredients. Then change one assumption at a time: attractive rather than oscillatory kernel; absence of explicit particle drag; release of radial confinement; central gravitational potential; one versus two response states.

Keep units, base equilibrium and source gain matched where the comparison permits it. Compute the full coupled spectrum, not just one root guessed from a simulation. Specify whether a difference arises from geometry, dissipation, kernel sign or additional memory dynamics.

**Done when:** a mapping table and benchmark distinguish an exact limit, a parameter adaptation and a genuinely additional prediction. “Memory changes stability” and “geometric instability is distinct from memory instability” cannot be submitted as new broad results.

## 3. Establish the inertial-chemotaxis limit [C2010, CS2007, CS2008]

Write an explicit map between particle positions, the chemical field and the RUT depth field. State which terms are set to zero and how point production becomes a Gaussian source. Check the cold/memoryless limit of the response matrix against the known attractive-interaction problem.

Separate a pressure-supported hydrodynamic closure from a collisionless rotating distribution function; their thresholds need not agree. Reproduce one simple pressure/damping threshold before testing the actual warm RUT annulus. This provides an independent physics benchmark beyond agreement between two implementations of our own equations.

**Done when:** the familiar stabilizing role of velocity spread is separated from a verified, model-specific quantitative effect of the two-stage field. Do not transplant a chemical threshold as a galaxy prediction.

## 4. Test whether two-stage maturation supplies a new result, not only a known filter [SGH2011, Z1973]

Hold the stationary response fixed and compare the full coupled system under instantaneous, one-stage and two-stage responses. Match initial field histories, not just initial particle positions. State the exact class of allowed kernels and source laws for any claimed bound or scaling.

A useful target is a derived relation between support, torque and mode growth that predicts a new configuration without retuning. Explore whether another kernel in the existing literature reproduces the same effect after a parameter transformation. If it does, label the outcome an equivalence or replication rather than a new mechanism.

**Done when:** the claimed result survives the closest established response family and its mathematical scope is stated. Two cascaded relaxations, a damped second-order equation and the corresponding transfer function are not three independent discoveries.

## 5. Use phase feedback and collective experiments constructively [O2017, LV2026]

Compare the impact-phase adaptation in orbiting droplets with our strictly linear maturation rule. It suggests a physically different nonlinear response worth investigating, but an external bath's available energy must not be silently imported.

For the 2026 cluster preprint, retain its short-range repulsion and drive in the reference reproduction. Determine which parts of collective survival require those ingredients. The isolated droplet loses bounce synchronization and coalesces; that is not our ring's m = 2 mode. Only compare matched outcomes after distinguishing the two instabilities.

**Done when:** the project can state a specific difference from both precedents. A stable common-field cluster or a restored orbit, without more, is not a defensible first-ever claim.

## 6. Finish the exact-identity and geometry ledger [FB2005, BM1984, M2010]

For each field completion, list its action or equation, boundary conditions, source and observable acceleration. Verify spherical equivalence, then test a nonspherical source to expose differences. Compare a known algebraic acceleration law with a true field solve without defining source mass from the field being fitted.

The scope must distinguish equation identity, spherical-limit equivalence and observational degeneracy. The existing PM score does not become a RUT score until the dynamical model derives and predicts the same field.

## 7. Keep established numerical methods distinct from the physics [P2024, Beyn2012, Liu2026]

Use published stellar-response benchmarks to validate orbital quadrature and response matrices. Use known nonlinear spectra for contour solvers, including root-free domains and difficult near-boundary cases. Consider region partitioning only after understanding its declared assumptions.

For the current stage-8 discrepancy, first complete population-specific numerical verification and full-state seeding. Prior-art attribution cannot repair an under-resolved basis or account for a mismatched growth rate. The direct methods comparison and the physical novelty comparison need separate outcomes.

## 8. Define a paper-sized contribution that remains

For every candidate contribution, complete this sentence with citations:

> Previous work established X. Under the explicitly different assumption Y, we establish Z, which the cited work did not establish in that form or regime.

Possible targets are a general support-torque bound, a verified change in the kinetic mode spectrum at fixed static attraction, a reproducible supported-population regime, or a distinctive finite-history prediction. Do not preselect any as original before the comparison is finished. A new, correct result can be built from known ingredients; a new name cannot supply missing novelty.

The novel-physics question and the observed-gravity question remain separate. The latter additionally requires a universal dimensional source law, energy/momentum and causal consistency in the claimed domain, and predictions against measurements. Neither requires adopting an expanding cosmology merely because some antecedent models did.

## Work products and preservation

Create new comparison files outside frozen experiment directories until their interfaces and protocols are declared. Keep one table per comparator with equations, parameters, state/history initialization, drive/reservoir, boundary conditions, matched quantities and outcomes. Save negative outcomes too.

Do not rerun or re-pin historical science archives as part of a citation edit. The current audit's offline identity checks can run independently:

```sh
python -B research_plan/prior-art/verify.py
```

They are not replacements for the research suite or stage 8's unresolved numerical verification. Publication-ready originality requires a review of the resulting comparisons, not a successful exit code from this audit.
