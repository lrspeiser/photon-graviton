# Distance adoption and selection: implications for a fresh test

## Source findings

The [ELVES-Field paper](https://arxiv.org/abs/2602.16766), pages 12 and 15-16, describes querying SIMBAD but does not specify the per-object query result, match radius, retrieval epoch or bibliographic chain in those method passages. Those passages do not resolve the two LSBC associations.

For confirmed dwarfs, the distance priority is TRGB, then SBF, then redshift-derived distance. The paper does not describe blending the redshift distance into an adopted TRGB/SBF value. It assigns TRGB a typical 5% error. It averages upper/lower SBF errors for a symmetric reported error bar. Confirmation uses a two-sigma distance lower bound within 10 Mpc; SBF confirmation additionally requires trustworthiness, signal-to-noise above five and visual review.

The publisher header defines its three SBF flags as default mask/annulus use, ambiguous mask choice and failed visual check. **The first exported flag is not the separate trustworthiness flag.** Thus the initial parser's decision not to require a True first flag was appropriate; manual mask adjustment alone is not a failure.

Local paper SHA256: a418469836986225b2d2a21e91dbd96b13dad9c0ff3376584cb25267235a7c3c. Publisher table SHA256: 890b1d178feb505cb2d9105926d13b77bc639d882e18776a3b3561b0650d8903.

## What this establishes, and what it does not

The current 24 candidates consist of nine TRGB and 15 SBF entries. The distance-priority procedure supports retaining them as potential non-redshift-distance tests; it does not certify freshness, source identity or fully independent selection. A typical error assigned by the source is not a newly measured uncertainty for each target. A symmetric published error bar does not prove the underlying measurement likelihood is symmetric or Gaussian.

We preserve the user-adopted distances. The diagnostic below does not revise any galaxy's distance or infer a per-object correction. It demonstrates why a selection-aware uncertainty analysis may matter if interpreting the catalog's quoted errors as a measurement model.

## Synthetic calculation: a lower-bound cut changes the accepted errors

**Provenance: established Gaussian truncation mathematics applied to a deliberately simplified example, not a new physical law or a fitted model.** Let the measured distance be D_hat = D + epsilon, with epsilon normally distributed with mean zero and a fixed sigma. Accept the synthetic object when D_hat - 2 sigma <= 10 Mpc. This is a simplified representation of one source selection condition, not the entire survey pipeline.

Define a = (10 Mpc + 2 sigma - D)/sigma. With phi and Phi the standard-normal density and cumulative distribution:

P(accepted | D) = Phi(a),

E[D_hat - D | accepted, D] = -sigma phi(a)/Phi(a).

These are forward sampling statements at a fixed true distance, not a posterior distance correction for an observed candidate. No population distribution or galaxy-motion model is specified.

For a declared toy error sigma = 0.5 Mpc, the measured-distance cutoff is 11 Mpc:

| Synthetic true distance (Mpc) | Accepted fraction | Mean distance error among accepted objects (Mpc) | Conditional coverage of nominal 95% distance intervals |
|---:|---:|---:|---:|
| 8 | approximately 100% | approximately 0 | 95.00% |
| 10 | 97.72% | -0.0276 | 97.21% |
| 11 | 50.00% | -0.3989 | 95.00% |
| 12 | 2.28% | -1.1866 | 0.00% |

The final row is a rare-selection example: only measurements fluctuating down to 11 Mpc or below are admitted, and an interval extending 0.98 Mpc above such a measurement cannot reach 12 Mpc. This does not predict zero coverage in the real sample. The toy omits variable errors, population density, photometric selection, SBF quality and spectroscopic availability.

In plain language, a farther galaxy can enter a nearby catalog when its measured distance happens to be too small. The galaxies that enter this way therefore do not have a balanced mixture of overestimates and underestimates. Ignoring that can create a distance-dependent prediction residual even if the preselection measurement errors were unbiased. Its size in ELVES is not established here.

## Verification and consequence

Run `python research_work/results/redshift-priority/elves_selection_diagnostic.py`. The four distances and fixed error are synthetic constants declared in code, not candidate values. Independent numerical quadrature verifies the analytic accepted probabilities and conditional mean errors. Output retains the assumptions, full toy results and hashes of the feature-only inputs used for method counts. No target distance, redshift or residual is accessed by this calculation.

A final predictive interval model must account for the actual source likelihoods and relevant selection, or state a defensible limited conditional claim. It must also handle independent motion/endpoint uncertainty, which this example does not supply. Do not substitute the source's flow-model residual scatter as our nuisance distribution. Continue original spectral-source association checks and investigate independent inputs rather than adding flexible distance terms to absorb these uncertainties. Five candidates remain excluded and 24 pending; no fresh predictive improvement is claimed.
