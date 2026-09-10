# Can deposits simply add the same extra source to all ordinary matter?

**A position-independent loading per unit ordinary matter cannot reproduce our fixed extra-field source.** The required loading varies by about sixfold even on the inner sampled grid, and roughly 1,460-fold over the full source grid. This is a conditional failure of a simple attachment rule, not a rejection of all companion deposits or an observed exclusion of the full hypothesis.

The comparison uses the extra source reconstructed in the preceding full-bar audit and the declared analytic ordinary-matter model that feeds the field calculation. It includes the bar, nuclear components, two stellar disks, two gas disks and the softened central source. Disk parameters are read directly from the field-builder code. There is no dark-halo component and no new measured stellar dataset.

## The branch being tested

Hypothetical postulate: deposits remain attached to ordinary matter, and each unit of ordinary mass carries the same extra effective source. In the Newtonian source interpretation:

`rho_extra(x) = eta * rho_b(x)`, with one nonnegative constant `eta`.

This is a familiar proportional-density model applied as a candidate companion-storage rule, not a novel fundamental formula. The existing hypothesis did not require uniform loading; this audit tests that specific optional closure rather than attributing it to every version of the theory.

For example, standard optically thin transport with capture coefficient per length `beta = kappa * rho_b` gives deposited power density `Q_cap = c * kappa * rho_b * u_c`. If every location experiences the same integrated companion energy density, the deposited energy remains attached without redistribution, and the deposited source is its energy divided by c squared, then `eta = (kappa/c) * integral(u_c dt)` is spatially constant. Equal fully occupied storage capacity per unit ordinary mass would also lead to the same proportional shape. These are conditional deductions using known transport and mass-energy bookkeeping. No capture constant, exposure history or energy supply is measured here.

Uniform illumination alone does not guarantee these conditions: shielding, differing ages, varying capture coefficients, receiver capacity and matter motion can all change loading. Those effects require a quantitative law before they can explain the needed variation.

## Measured mismatch between the two model shapes

For the fixed target define `q(x) = rho_extra,target(x)/rho_b(x)`. Exact proportional loading requires q to be constant. The ratios below are model outputs, not measured extra mass attached to individual stars.

| Grid subset | Points | Smallest q | Largest q | Range factor |
|---|---:|---:|---:|---:|
| full grid | 240 | 0.1981 | 290.0226 | 1463.68 |
| inner sampled region R le 8 z le 1 | 120 | 0.1981 | 1.1944 | 6.03 |
| midplane only | 40 | 0.2002 | 1.0303 | 5.15 |

The full-grid minimum is at R=0.5 kpc, z=0.1 kpc, bar angle pi/4. The maximum is at R=20 kpc, z=4 kpc on the bar axis. The latter region lies outside the current red-giant orbit selection, so the table separately reports R<=8 kpc and z<=1 kpc. The effect is not limited to that distant off-plane extreme: midplane ratios alone span a factor of 5.15.

Changing to the lower-resolution reconstructed source gives a full-grid range factor of 1,464.09; the direct requested div(Q) source gives 1,459.36; the finer reconstruction gives 1,463.68. Thus the order-of-magnitude full-grid mismatch is not removed by these existing source-resolution alternatives. This does not incorporate uncertainty in the ordinary-matter mass model or prove the target extra field fits all observations.

## A normalization-independent bound

Known minimax algebra gives the smallest possible worst fractional error of any constant eta on this grid:

`min_eta max_x |eta/q(x)-1| = (q_max-q_min)/(q_max+q_min)`
`eta_opt = 2*q_min*q_max/(q_min+q_max)`.

This follows by balancing the opposite-sign relative errors at the smallest and largest q. An independent numerical scalar minimization checks the result. It is a deterministic mismatch bound for this target and chosen fractional loss, not a chi-squared statistic or a confidence level.

The bound is 99.86% over the full grid and 71.54% over the inner subset. Under a symmetric multiplicative-error criterion, even the best constant must miss some full-grid target by a factor of at least 38.26, or 2.46 on the inner subset. Different loss functions change the preferred constant; none makes unequal q values equal.

## Consequences for the hypothesis

This nonproportionality is expected from the nonlinear empirical field template; it is not a new observational discovery. The calculation quantifies how much spatial variation a simple attachment interpretation would have to explain in the actual three-dimensional model.

If we retain this target gravity field, the next physical closure must explain at least one of: nonuniform accumulated loading on ordinary matter, a deposited reservoir with its own spatial support, or a gravitational response that is not simply the local deposited energy as an ordinary Poisson source. A specified combination is possible. Merely setting eta(x)=q(x) would reproduce the target by definition and would not be an independent prediction.

For uneven loading, the required ratio is a target for a shared transport/capture/history calculation. For a separate reservoir, the earlier mechanical-support and long-lived-capture constraints apply. For modified sourcing, an explicit common field equation must determine stellar acceleration and lensing while accounting for energy and momentum. Naming the carriers gravitons does not choose among these equations.

The results do not establish a new gravitational field fit, determine cosmic ages, calculate the deferred total photon budget, or open observational holdouts. They remove one overly simple source-shape closure and identify what a nontrivial physical connection has to explain.

## Reproduction

Run `run.py` followed by `report.py`. Every probe, analytic baryon component and required loading is retained in `results.json`; the restricted-grid checks are in `assessment.json`. Source hashes and the previous source-audit hashes are verified. Existing catalogs, force coefficients, source grids and sample roles remain unchanged.
