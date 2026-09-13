# Does the proposed opacity describe the specified receiver?

The preceding goal turn made progress by establishing a local capture energy/momentum ledger. This continuation connects receiver identity to absorption geometry. All six objectives remain open.

## Ordinary-matter receiver comparison

Suppose ordinary matter receives and retains companion energy with constant absorption cross section per unit original receiver mass K. Then kappa=K rho_b. K has units area/mass; kappa has units inverse length. This is an optional comparison law, not a selected microscopic mechanism. The calculation neglects changes in receiver state, motion and mass after absorption, and assumes straight rays.

For ordinary Plummer matter of total mass M and scale a, known projection gives

    Sigma_b(b)=M/(pi a^2) [1+(b/a)^2]^-2,
    tau(b)=T [1+(b/a)^2]^-2,  T=K M/(pi a^2).

Applying the known absorption cross-section integral yields

    sigma_eff=pi a^2 [sqrt(pi T) erf(sqrt(T))+exp(-T)-1].

This is a derived closed form for the chosen comparison profile, not a novelty claim. For any nonnegative finite ordinary-matter distribution obeying kappa=K rho_b under the same ray assumptions,

    sigma_eff <= integral tau(b) d^2b = K M,

because 1-exp(-tau)<=tau. In weak absorption, sigma_eff approaches K M: there is little overlap or shielding between receivers. Increasing the region's size at fixed mass and K cannot exceed the sum of its receivers' cross sections.

## What this means for a bigger collecting region

The earlier pi R^2 c u_c f_capture formula remains correct for a defined collecting sphere. Its factor-of-four increase when R doubles requires f_capture to remain constant. It does not guarantee that a fixed population of receivers maintains that fraction when diluted across a larger volume. Here a is a smooth density scale, not a hard outer boundary; the Plummer distribution extends to infinity.

| Central optical depth T | sigma_eff/(K M) | Area increase when a doubles at fixed M and K |
|---:|---:|---:|
| 0.001 | 0.999833 | 1.000125 times |
| 0.1 | 0.983661 | 1.012396 times |
| 1 | 0.861528 | 1.114679 times |
| 10 | 0.460499 | 1.575286 times |
| 1000 | 0.055050 | 1.963669 times |

For this smooth profile, the optically thick asymptote is sigma_eff approximately pi a sqrt(K M), while the thin asymptote is K M. A more massive cluster can still have more collecting area. This test does not estimate the actual external energy supply or rule out its dominance. It shows why both receiver abundance and capture probability must accompany a geometric-area argument.

## Receiver implication of our depth-to-fourth law

For the previous optional law kappa=chi*(-Phi)^4 in the same Plummer well, the required cross section per unit ordinary mass would be

    K_required(r)=kappa/rho_b
                =4 pi chi G^4 M^3/(3a) sqrt(1+(r/a)^2).

Therefore it cannot simultaneously represent identical ordinary receivers with a universal constant K. Relative to the same well's center, K_required is 1.414 times larger at r=a, 3.162 times at 3a and 10.050 times at 10a. Across homologous wells it also scales as M^3/a at fixed r/a. The total opacity still declines outward; it is absorption strength per unit ordinary matter that increases.

This is not a contradiction in a field-based hypothesis: a distinct field reservoir could have a different spatial distribution, or microscopic capture could depend on the environment. But either needs an explicit law. In a region with no ordinary receivers, a nonzero depth opacity cannot be attributed to those absent receivers. Conversely, a constant ordinary-matter opacity is a materially different spatial model and cannot reuse the previous depth-profile predictions unchanged.

## Verification and decision consequence

receiver-consistency.py and receiver-consistency-results.json preserve five optical-depth cases, four independent line-of-sight column checks, and twelve checks of the implied K_required formula across three synthetic wells. Independent numerical cross-section integration agrees with the closed form within 1e-8; column checks pass at that tolerance, and algebraic receiver ratios agree within 1e-13. Parameters are synthetic; no observations or final holdouts were fitted or opened.

The next receiver model must link its abundance, per-receiver interaction, energy retention and spatial support. The existing depth law can remain a phenomenological candidate, but its receiver cannot be silently identified with constant-cross-section ordinary matter. This constraint advances the shared rotation/lensing problem: a single physical receiver model must supply the opacity and the persistent gravitating distribution together.
