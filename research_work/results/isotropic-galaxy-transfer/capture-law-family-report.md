# No unified capture equation emerges from this small family comparison

Six declared effective-retention formulas were evaluated with common solar normalization, and five fitted one exponent on Earth before predicting Moon and the remaining planets/Pluto. The constant model uses the Sun alone. No candidate reproduces the nine prediction masses well; these formulas are not promoted to universal laws. See the [protocol and implicit prediction equation](../../../research_plan/capture-law-family-protocol.md).

| Predictor for effective retention | Fitted exponent | Prediction log10 RMS | Prediction bodies within factor 2 |
|---|---:|---:|---:|
| Constant | 0 | 6.8224 | 0/9 |
| Radius | -2.87464 | 2.7426 | 0/9 |
| Well-depth proxy GM/R | -1.68095 | 0.8856 | 2/9 |
| Surface acceleration GM/R^2 | -4.04806 | 0.4286 | 3/9 |
| Mean-density proxy M/R^3 | +9.91693 | 0.7007 | 4/9 |
| Radiative surface flux L/R^2 | -1.11215 | 0.7007 | 4/9 |

These are descriptive errors with equal body weight, not uncertainty-based likelihoods. Two-parameter candidates reproduce Sun/Earth by construction. Rounded dynamically inferred masses and stipulated output-power/path proxies limit their interpretation. Source and retention histories are not observed or fitted independently.

## What the predictions say

The least-bad aggregate candidate, inverse surface acceleration, predicts mass ratios relative to catalog values:

| Body | Predicted / catalog mass |
|---|---:|
| Mercury | 3.187 |
| Venus | 1.242 |
| Moon | 4.691 |
| Mars | 1.969 |
| Jupiter | 0.3318 |
| Saturn | 0.5984 |
| Uranus | 0.4513 |
| Neptune | 0.2993 |
| Pluto | 2.662 |

This is not close agreement across the objects. Its negative exponent also means weaker acceleration requires more effective retention, contrary to a simple positive preference for deep/strong wells. The density candidate does better for the Moon (1.653) but predicts giant-planet masses 6.4–14.7 times their targets. A selected successful subset would misrepresent the result.

## Important algebraic checks

Mass-dependent predictors were solved implicitly for the predicted mass. The calculation never inserts a prediction object's observed mass into a supposed prediction of that same mass. None of the fitted implicit exponents is singular; each computed prediction satisfies its implicit equation numerically.

Density and radiative-flux candidates yield the same predictions in this setup. In the thin-conversion regime with D=R, P_c is proportional to L*R. Solving M=P_c*K*(M/R^3)^p gives powers of L and R equivalent to an appropriate power law in L/R^2. Their matching errors are therefore a reparameterization, not two independent pieces of support. Six formula labels represent only five distinct predicted families here.

## Relation to the three storage scenarios

These formulas model K=eta*T for constant, lossless accumulation or K=eta*tau for steady regeneration. They do not distinguish those mechanisms, determine eta separately, establish an energy history or give a graviton frequency. An independently bound reservoir still needs a confinement/stability law; fitting K does not supply one. Free-streaming escape remains constrained by the prior residence calculation.

The full-Solar-System-mass target also differs from galaxy extra-gravity residuals. The separate radiation-conditioned galaxy trial tests a related environmental idea under the existing galaxy transport calculation; it cannot be called the same fitted universal equation without a source/target bridge.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/capture-law-family.py`. All observed/predicted masses are saved in `capture-law-family-results.json`. No previously unseen observations are claimed. No fixed cosmic age, arbitrary per-body normalization, or observational confirmation of gravitons is introduced.
