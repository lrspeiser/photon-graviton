# Stage 5C: Nonlinear Constraint Closure and the Finite-Regulator Obstruction

**Status:** **PASS for the nonlinear smooth/continuum constraint algebra; explicit failure of exact closure for the current fixed central-difference regulator.** The derived scalar constraint closes with the spatial-diffeomorphism generator, preserves nontrivial constrained evolution without projection, and retains two principal physical modes. The local lattice errors fall quadratically but are nonzero at every finite spacing.  
**Date:** 2026-09-21  
**Issue:** #5  
**Executable:** [`check_nonlinear_constraint_closure.py`](check_nonlinear_constraint_closure.py)  
**Frozen output:** [`results/nonlinear_constraint_closure_full.json`](results/nonlinear_constraint_closure_full.json)

## 1. Coefficients and nonlinear scalar constraint

Stage 5B supplies the curvature coefficient

\[
b=\frac{\kappa}{4}
=7.135031174451665\times10^{-7}.
\]

The inherited frame kinetic residue is

\[
U_{\rm frame}=0.006789913753821769,
\]

so the canonical coefficient is

\[
a=\frac{U_{\rm frame}}{2}
=0.0033949568769108844.
\]

The tested nonlinear scalar density is

\[
\boxed{
\mathcal H
=
\frac{a}{\sqrt g}
\left(\pi^{ij}\pi_{ij}-\frac12\pi^2\right)
-b\sqrt g\,R.
}
\]

The curvature term includes the torsion-free connection energy after elimination. No second connection-energy term is added.

The spatial generator is

\[
D[\xi]=\int d^3x\,\pi^{ij}\mathcal L_\xi g_{ij},
\]

with \(\pi^{ij}\) treated as a weight-one contravariant density.

## 2. Cubic kinetic term is fixed by covariance

Starting only from

\[
T_2=\operatorname{tr}(\pi^2)-\frac12(\operatorname{tr}\pi)^2
\]

and the declared metric/density transformation laws, the four parity-even ultralocal \(h\pi\pi\) coefficients are uniquely fixed to

\[
\boxed{
(-\tfrac12,\tfrac14,2,-1).
}
\]

Equivalently,

\[
\boxed{
T_3
=2\operatorname{tr}(h\pi^2)
-\operatorname{tr}(\pi)\operatorname{tr}(h\pi)
-\frac12\operatorname{tr}(h)\operatorname{tr}(\pi^2)
+\frac14\operatorname{tr}(h)(\operatorname{tr}\pi)^2.
}
\]

The coefficient system has rank four and the held-out covariance residual is

\[
1.10\times10^{-13}.
\]

This is a bootstrap from spatial covariance; microscopic derivation of the kinetic residue remains a separate open question.

## 3. Exact smooth constraint algebra

For smooth periodic fields evaluated with spectral derivatives, the executable verifies

\[
\boxed{
\{D[\xi],D[\eta]\}=D[[\xi,\eta]],
}
\]

\[
\boxed{
\{D[\xi],H[N]\}=H[\mathcal L_\xi N],
}
\]

and

\[
\boxed{
\{H[N],H[M]\}
=ab\,D\!\left[g^{-1}(N\,dM-M\,dN)\right].
}
\]

At \(L=16\), the relative residuals are:

| Bracket | Relative residual |
|---|---:|
| \(D-D\) | \(1.85\times10^{-16}\) |
| \(D-H\) | \(1.10\times10^{-15}\) |
| \(H-H\) | \(0\) at reported precision |

The absolute \(H-H\) bracket is small because the coherent curvature amplitude is small, but the independently calculated structure-function target agrees exactly.

## 4. Local finite regulator result

The same brackets were evaluated with a local centered derivative on successively finer cubic grids.

The relative errors scale as

\[
D-D:\ L^{-1.97480},
\]

\[
D-H:\ L^{-1.97214},
\]

\[
H-H:\ L^{-1.97098}.
\]

At \(L=48\), they are still nonzero:

| Bracket | Relative error |
|---|---:|
| \(D-D\) | \(3.97\times10^{-5}\) |
| \(D-H\) | \(1.33\times10^{-6}\) |
| \(H-H\) | \(4.11\times10^{-6}\) |

Thus the present local central-difference algebra is **not exactly first class at finite spacing**. It converges toward the continuum algebra with the expected \(O(a^2)\) product-rule error.

No repeated projection or tolerance adjustment is used to convert that failure into a pass.

## 5. Unprojected nonlinear evolution

Start from flat metric and zero momentum, then evolve with a spatially varying lapse. The lapse is rescaled by the common coefficient product only to expose nontrivial evolution in a short coordinate-time test.

For the smooth spectral regulator, after 60 RK4 steps:

- manual projections: `0`;
- maximum metric change: \(4.36\times10^{-4}\);
- maximum momentum: \(1.48\times10^{-4}\);
- maximum scalar-constraint violation: \(3.55\times10^{-20}\);
- maximum vector-constraint violation: \(1.43\times10^{-17}\).

The evolution is therefore nontrivial and remains on the nonlinear constraint surface without manual projection.

For the local central regulator, the vector-constraint drift decreases with refinement:

| \(L\) | Maximum vector drift |
|---:|---:|
| 8 | \(1.63\times10^{-10}\) |
| 12 | \(9.72\times10^{-11}\) |
| 16 | \(6.29\times10^{-11}\) |
| 24 | \(3.03\times10^{-11}\) |

This is convergence, not exact finite-lattice preservation.

## 6. Nonlinear principal mode count

For 200 random positive-definite determinant-one backgrounds and nonzero covectors, the principal scalar/vector constraint Jacobian has:

- rank four;
- vanishing first-class bracket matrix to a maximum residual \(5.87\times10^{-14}\);
- physical configuration count \(\boxed{2}\).

No scalar mode is regenerated at the nonlinear principal level.

### Local determinant/trace conditions

Adding the local conditions

\[
\log\det g=0,
\qquad
g_{ij}\pi^{ij}=0
\]

does not create another independent reduction on top of four first-class constraints. The classification becomes

\[
\boxed{
4\ \text{second-class}+2\ \text{first-class},
}
\]

and still leaves

\[
\boxed{2\ \text{physical configurations}.}
\]

The local volume pair therefore acts as a gauge fixing of two old directions. The global homogeneous volume pair remains a separate second-class construction.

## 7. What has been completed

The requested local nonlinear program has now been executed to the following extent:

1. **Local determinant/trace construction:** classified consistently as gauge fixing at nonzero momentum.
2. **Scalar/vector constraints:** retained and shown to close in the smooth nonlinear theory.
3. **Quadratic-order brackets:** surpassed by testing the full nonlinear metric expressions on finite-amplitude smooth fields.
4. **Frame and connection energy:** included through the nonlinear DeWitt kinetic term and the coherent Palatini curvature term, with the auxiliary connection eliminated once.
5. **Constraint preservation:** demonstrated during nontrivial unprojected nonlinear evolution.
6. **Mode count:** two configurations; no principal scalar mode.
7. **Nonlinear interaction:** the spatial curvature term comes from the coherent frame-area junction, while the kinetic cubic coefficients are fixed by covariance rather than inserted as a fitted Einstein–Hilbert vertex.

## 8. What has failed or remains open

The strict finite-local version of the requested program is **not complete**, because the present centered derivative lacks an exact Leibniz rule. Its constraint algebra is only emergent as the regulator is removed.

Also still open:

- a unique microscopic derivation of the temporal kinetic coefficient;
- an exact finite quantum nonlinear constraint algebra and anomaly analysis;
- strong-field global solutions and collapse;
- local vacuum-energy response or sequestering;
- regulator universality beyond spectral versus central controls;
- Newton’s constant and empirical predictions.

## 9. Decision

The project now has a closed **classical nonlinear continuum candidate** with coefficients inherited from the finite coherent junction and the frame kinetic sector. It does not yet have an exact finite local nonlinear gravity theory.

The next branch point is explicit:

1. construct a finite move/derivative algebra with exact discrete product rules; or
2. accept diffeomorphism symmetry as emergent and establish regulator-independent continuum observables.

In either case, the next physical calculations should test local vacuum sources and strong-field configurations using the closed continuum constraints, while preserving the finite-lattice closure failure as an unresolved gate.

## 10. Reproduction

```sh
python phase_junction_network/validation/check_nonlinear_constraint_closure.py \
  --output phase_junction_network/validation/results/nonlinear_constraint_closure_full.json
```

Use `--quick` for the smaller regression set.
