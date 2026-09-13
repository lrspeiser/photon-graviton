# Observer-regular area response improves farther brightness predictions

The previous bounded-area fit improved brightness transfer, but left its optical interpretation unspecified. This run imposes one necessary local condition for a metric-geodesic interpretation and tests a replacement with the same number of free parameters. It improves the farther standardized-magnitude score while retaining the redshift rate and event exponent. It is not yet a global metric or a completed interaction theory.

## Conditional optical derivation

Assume metric null geodesics, photon conservation, and D as accumulated local path length. Normalize the past-directed affine parameter by observer photon energy. Then dD/dlambda=S=exp(alpha D), so

    lambda=(1-exp(-alpha D))/alpha=f/alpha.

Known reciprocity gives D_A=D_G/S. For a regular observer vertex and finite optical tidal field, the known Jacobi initial conditions require D_A=lambda+O(lambda^3); see [the Jacobi-map framework](https://arxiv.org/abs/1308.4935). This is not an expansion assumption, but it is also not automatically applicable to non-geodesic energy-conversion rays. In particular, the path length, affine normalization and local photon energy identification are additional assumptions of this branch.

With B(f)=D_G^2/D^2, the previous B=1+eta f yields

    D_A=lambda+(eta-1) alpha lambda^2/2+O(lambda^3).

Thus eta=0.61512 does not satisfy regularity under these assumptions. This does not reject its algebraic brightness fit or every nonmetric optical mechanism. It shows that it cannot simply be relabeled as this smooth metric-geodesic branch.

The new **postulated** bounded response is

    B(f)=1+f/(1+q f), f=1-1/(1+z), q>=0.

The slope B'(0)=1 is fixed by the local condition; only saturation parameter q is fitted. Expansion gives

    D_A=lambda-alpha^2(q/2+13/24)lambda^3+O(lambda^4).

For D_A''=-R_opt D_A this implies finite observer R_opt(0)=alpha^2(3q+13/4). Symbolic expansion verifies these coefficients. These are known optical identities and a project response ansatz, not a novel field equation, proof of a global metric, or a derivation of companion production. Inserting a desired R_opt along one ray does not solve the gravitational field equations or consistency across observers.

## Actual brightness fit and transfer

Fit q on the same 466 nearer rows and freeze it for the 494 farther rows, recalibrating the common source magnitude consistently using the 77 calibrator rows. Alpha remains 0.0002488993286/Mpc, b remains 1, and photon number is retained. The full inherited covariance and frame convention remain unchanged. The data and earlier failures are exposed, so this is exploratory model development, not blind discovery.

Fitted q=2.260139, away from its declared bounds. At z=1, B=1.234734; at arbitrarily large z, B approaches 1.306735. The conditional observer optical coefficient is 6.21393e-7 Mpc^-2. It is not an independently measured matter density or stress-energy source.

| Group | No area correction | Earlier bounded area | Observer-regular area |
|---|---:|---:|---:|
| Training, 466 rows | 468.030 | 445.421 | 444.544 |
| Farther prediction, 494 rows | 424.482 | 402.335 | 389.177 |

These are point-prediction GLS scores on the same covariance, not significance levels. Each area model has one fitted shape parameter. No farther rows are used to fit q, but repeated exploration of this split limits claims of generalization.

Correlated mean residuals for the new model are 0.00999 mag at z=0.1-0.3, 0.02559 at 0.3-0.6, 0.00182 at 0.6-1, and 0.11591 at 1-3. The highest-redshift discrepancy is not eliminated; a lower aggregate score does not mean every observation or bin is matched. Released standardized magnitudes also retain source/color/selection assumptions and have not been reprocessed with a new physical propagation model.

## Consequence

A locally admissible optical response can improve this conditional brightness prediction without increasing event duration beyond the adopted timing relation. That is a stronger candidate than an unconstrained area correction, but only one local consistency condition has been supplied. The missing work is a common spacetime or interaction model that generates the redshift, arrival mapping and beam response with energy accounting and consistent angular sizes/lensing. The current deposited-gravity fits cannot simply be combined with this new optical law without that derivation.

All six goals remain open. Reproduce with `python research_work/results/brightness-distance-consistency/bounded-area.py --regular`. Separate regular-area result/prediction files preserve the previous bounded-area artifacts.
