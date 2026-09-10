# Momentum implied by the fitted shielding model

**An isotropic external bath does not give balanced local capture after shielding.** The computed directional imbalance ranges from 0.0079 to 0.6408 on the 240-point grid. In all sampled locations, its projection along the reference gravitational acceleration is positive. This is a consequence of the specified transport model, not a measured force on Galactic stars.

The calculation fixes the finer shielding fit's capture coefficient and normalization. It uses its actual angular transmission arrays instead of assuming that opposite directions cancel. No new coefficient is fitted, no target density is substituted for the model's own loading, and no holdout is opened.

## Formula provenance and physical assumptions

Let n be the direction in which a companion travels toward the receiving point. The earlier rays were backtracked as x-s*n. Standard radiation-momentum bookkeeping for full absorption of massless companions with momentum E/c gives:

`xi(x) = integral I(x,n)*n dOmega / integral I(x,n) dOmega`
`F_capture = (P_abs/c) * xi`.

These are known energy-flux and momentum-flux relations applied to the hypothetical companions. A different dispersion relation or a channel that emits compensating momentum needs a different calculation. The result does not establish that companions are ordinary gravitons.

Assume the receiving ordinary matter initially has mass M_b, and retained deposit energy adds inertial mass q*M_b under the Newtonian source interpretation. For constant loading over an unknown duration T, the snapshot mass is M_b*(1+q), while `P_abs = q*M_b*c^2/T`. The leading stationary-receiver acceleration is therefore:

`a_capture = [q/(1+q)] * (c/T) * xi`.

This is a conditional derivation from the preceding assumptions, not a new fundamental force law. q is the shielding model's fitted loading, which was already shown not to match the target source well. T is not measured or fitted here. Nonuniform loading histories replace q/T with the current deposited-energy rate per original rest energy.

For comparison with the magnitude g_ref of our frozen full-bar target acceleration, define:

`T_critical = [q/(1+q)] * c*|xi|/g_ref`
`|a_capture|/g_ref = T_critical/T`.

T_critical is the loading duration for equal instantaneous force magnitudes under the constant-rate snapshot assumption. T>10*T_critical makes this leading capture term less than 10% of the reference gravity. The reference field is NOT recomputed from the mismatching shielding density, so this is a diagnostic force scale, not a self-consistent orbit solution.

## Results

T_critical ranges from 0.030 to 3.477 billion years. Keeping the leading capture term below one tenth of the reference gravity at every sampled point would require T greater than about 34.77 billion years under these assumptions. This does not establish a cosmic age or rule out an older fictional universe.

The following examples lie on the bar axis. z denotes height above the disk plane; the imposed reflection symmetry supplies the corresponding reversed vertical direction below it.

| R (kpc) | z (kpc) | Model loading q | Directional imbalance | T_critical (Gyr) |
|---:|---:|---:|---:|---:|
| 1 | 0 | 0.141 | 0.2120 | 0.269 |
| 1 | 1 | 1.207 | 0.4752 | 2.401 |
| 1 | 4 | 2.577 | 0.1657 | 3.432 |
| 8 | 0 | 2.496 | 0.0561 | 1.814 |
| 8 | 1 | 2.643 | 0.0880 | 2.853 |
| 8 | 4 | 2.886 | 0.0753 | 3.194 |
| 20 | 0 | 3.121 | 0.0082 | 0.922 |
| 20 | 1 | 3.109 | 0.0162 | 1.833 |
| 20 | 4 | 3.131 | 0.0147 | 1.750 |

The cosine of the angle between capture momentum and reference gravity ranges from 0.645 to 1.000. Capture generally pushes inward or toward the plane in this calculation, rather than acting as an isotropic pressure with zero mean force. This does not prove that the matter collapses; support, motion and the evolving field must be solved together.

The force is directly exerted on an absorbing receiver. It cannot simply be counted as extra spacetime curvature: lensing light would need a separately demonstrated response. Any use of capture forces to explain stellar motions must remain consistent with the gravitational field used to predict lensing.

## Checks and limits

The angular weights sum to one, their unattenuated vector moment vanishes, and the transmitted mean direction has magnitude no greater than one. An independent optically thin uniform-sphere calculation approaches the analytic limit `xi = -beta*x/3`, verifying the travel-direction sign and first angular moment.

Holding the fitted coefficient fixed, coarse/fine quadrature changes xi by at most 0.01283 in absolute vector norm. Moving the finer boundary from 30 to 60 kpc changes it by at most 0.000280. These are absolute changes: small anisotropies need not have correspondingly small relative errors. No precision bound at every ray or point is claimed.

Moving receivers experience Doppler-dependent sampling and mass-growth drag in addition to this leading stationary term. The earlier all-direction capture calculation addresses that distinction. Capture recoil, changing optical depths, receiver trajectories, possible outgoing radiation and self-consistent gravity remain absent from the static shielding model.

Thus this audit supplies a momentum requirement for the candidate, not a successful evolution model. The previous source-shape mismatch remains. Long loading times might reduce the present instantaneous force, but do not automatically eliminate cumulative changes in inertia and orbits. No total source-energy budget has been evaluated.

## Reproduction

Run `run.py`, then `report.py`. The result retains vectors, reference accelerations, predicted loadings and critical times at every location for all three quadrature/boundary cases. Cached transmission arrays and field inputs are checked by hashes. No original model, catalog or sample role is changed.
