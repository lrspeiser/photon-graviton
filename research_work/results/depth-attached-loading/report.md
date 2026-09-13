# Monotone depth-only loading of ordinary matter

**Increasing attached loading with well depth alone does not reproduce the existing target source under common exposure.** Both ordinary-field and full-field potential orderings give a minimum worst relative source-density error of 71.54% in the inner sampled region. This is a conditional comparison of modeled source shapes, not a measured velocity error or an exclusion of the broader companion mechanism.

## Optional closure and provenance

Project postulate, not claimed unique: rho_extra(x)=eta(W(x)) rho_b(x), where eta is any nonnegative nondecreasing function and W=-Phi. Common illumination/exposure, common receiver properties, no redistribution and monotone capture can motivate this optional closure. They are assumptions, not consequences of energy conservation. We do not impose a power law or fit a chosen exponent; the test permits every monotone shape on the sampled points.

We use both the ordinary-matter potential and the existing full empirical potential. Only their ordering is used. The full stored potential can have a different additive reference and negative W values; W here is a signed potential-ordering coordinate, not a measured escape energy. Adding a constant changes no ordering or bound. No clipping at zero or inferred escape speed is used. A nondecreasing threshold law would still belong to this ordering class.

## Conditional algebraic result

Known minimax reasoning: for W_i <= W_j but required loading q_i > q_j, any monotone eta must incur worst relative error at least (q_i-q_j)/(q_i+q_j). Taking the largest such inversion gives the optimum on this ordered finite grid. An independent linear program minimizes the common error subject to all per-point relative-error intervals and monotonicity constraints, confirming each bound to 1e-7. This is mathematical auditing, not a new gravity formula.

| Potential ordering | Sample | Points | Unavoidable worst source-density error |
|---|---|---:|---:|
| ordinary | full | 240 | 99.86% |
| ordinary | inner | 120 | 71.54% |
| full | full | 240 | 99.86% |
| full | inner | 120 | 71.54% |

Inner critical shallower: R=8 kpc, z=1 kpc, phi=0; required loading 1.19439, signed -Phi=36642.1 (km/s)^2.

Inner critical deeper: R=0.5 kpc, z=0.1 kpc, phi=0.785398; required loading 0.198146, signed -Phi=116164 (km/s)^2.

## Meaning for the research

This is stronger than rejecting one fitted depth exponent, but narrower than rejecting depth-dependent capture. Unequal incoming intensity, shielding, receiver composition, different exposure histories, motion and redistribution can all break the assumed final-loading ordering. A capture rate increasing with depth does not guarantee accumulated loading increases with present depth.

The earlier radial-only audit allowed eta(R), whereas this test allows eta(-Phi); neither fits the fixed target exactly. A next transport calculation must supply exposure and capture jointly, rather than prescribe each location's final loading from the desired source. Alternatively, a separately supported reservoir or a changed gravity response requires its own dynamics.

The empirical target itself remains conditional and has not passed all motion/lensing tests; ordinary-matter uncertainties and alternative source reconstructions are not propagated here. No radiation budget or new observation is supplied. Source files and field arrays are hashed, critical pairs and all fitted monotone values are retained. All six scientific goals remain open.
