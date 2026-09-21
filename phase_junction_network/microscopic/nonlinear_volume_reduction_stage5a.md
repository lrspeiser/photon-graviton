# Stage 5A: Exact Nonlinear Homogeneous Volume Reduction

**Status:** **PASS at nonlinear homogeneous scope.** The fixed-volume condition and its Hamiltonian-preservation partner form an exact nonlinear second-class pair. They remove the one negative homogeneous conformal canonical direction and leave five positive tracefree shape directions.  
**Date:** 2026-09-20  
**Issue:** #5  
**Executable:** [`check_nonlinear_volume_reduction.py`](check_nonlinear_volume_reduction.py)  
**Frozen output:** [`nonlinear_volume_reduction_results.json`](nonlinear_volume_reduction_results.json)

## 1. Why this stage is necessary

The earlier periodic calculation controlled the homogeneous conformal instability only in the linearized frame theory. A proper gravity theory cannot rely on a condition that fails as soon as the metric is treated nonlinearly.

Stage 5A tests the exact finite-dimensional homogeneous geometry rather than another linear perturbation.

It does **not** assume that the local Phase Junction theory is already Einstein gravity. It uses the nonlinear DeWitt kinetic form as the target kinetic geometry that any successful massless spin-2 completion must reproduce in its homogeneous sector.

## 2. Nonlinear canonical variables

Let

\[
g_{ij}=g_{ji}>0
\]

be a finite positive-definite homogeneous three-metric and

\[
\pi^{ij}=\pi^{ji}
\]

its canonical momentum.

The two constraints are

\[
\boxed{
C_V=\log\det g=0
}
\]

and

\[
\boxed{
C_P=g_{ij}\pi^{ij}=0.
}
\]

The first fixes total homogeneous volume. The second removes its conjugate trace momentum.

Their exact Poisson bracket is

\[
\begin{aligned}
\{C_V,C_P\}
&=
\frac{\partial C_V}{\partial g_{ij}}
\frac{\partial C_P}{\partial \pi^{ij}}\\
&=(g^{-1})_{ij}g_{ij}\\
&=3.
\end{aligned}
\]

Therefore

\[
\boxed{
\operatorname{rank}\{C_A,C_B\}=2
}
\]

for every positive-definite metric. The pair is exactly second class, not merely at first order around the identity.

Across 500 random determinant-one metrics, the maximum numerical bracket error from three is

\[
1.33\times10^{-15}.
\]

## 3. Exact nonlinear conformal decomposition

The homogeneous kinetic term is

\[
T
=
\frac{1}{\sqrt{\det g}}
\left(
\pi^{ij}g_{ik}g_{jl}\pi^{kl}
-\frac12(g_{ij}\pi^{ij})^2
\right).
\]

Before reduction, its six-dimensional momentum-space inertia is

\[
\boxed{(1\text{ negative},5\text{ positive}).}
\]

The negative vector is the exact conformal momentum direction

\[
\pi^{ij}\propto g^{ij}.
\]

Indeed, for

\[
\pi^{ij}=\alpha g^{ij}
\]

at \(\det g=1\),

\[
T=-\frac32\alpha^2.
\]

Project an arbitrary momentum to

\[
\boxed{
\pi_T^{ij}
=
\pi^{ij}
-\frac{1}{3}(g_{kl}\pi^{kl})g^{ij}.
}
\]

Then

\[
g_{ij}\pi_T^{ij}=0.
\]

On that five-dimensional subspace,

\[
T[\pi_T]
=
\pi_T^{ij}g_{ik}g_{jl}\pi_T^{kl}>0
\]

for every nonzero \(\pi_T\).

The finite tests found:

- unreduced inertia `(1 negative, 5 positive)` in all 500 samples;
- reduced inertia `(0 negative, 5 positive)` in all samples;
- minimum reduced kinetic eigenvalue `0.0375141`;
- minimum sampled projected kinetic energy `0.261045`;
- maximum trace-projection residual `1.78e-15`.

Thus the nonlinear reduction removes exactly the conformal canonical pair rather than hiding its negative sign.

## 4. Hamiltonian preservation generates the second constraint

For the generalized DeWitt form

\[
T_\lambda
=
\frac{1}{\sqrt g}
\left(
\pi^{ij}\pi_{ij}
-\lambda\pi^2
\right),
\]

Hamiltonian evolution of the volume condition gives

\[
\dot C_V
=
\{C_V,H\}
=
\frac{2(1-3\lambda)}{\sqrt g}\,C_P.
\]

At the required spin-2 value

\[
\lambda=\frac12,
\]

this becomes

\[
\boxed{
\dot C_V=-\frac{C_P}{\sqrt g}.
}
\]

Therefore preserving fixed volume dynamically generates

\[
C_P=0.
\]

The direct derivative calculation and the analytic formula agree to

\[
7.11\times10^{-15}
\]

or better across the finite sample.

This extends the earlier linear observation to the exact nonlinear homogeneous metric.

## 5. Pure volume energy on the reduced shape manifold

Consider a vacuum or volume term

\[
U_\Lambda=\Lambda\sqrt{\det g}.
\]

A determinant-preserving shape path can be written as

\[
g(s)=g^{1/2}e^{sA}g^{1/2},
\qquad
\operatorname{tr}A=0.
\]

Because

\[
\det e^{sA}=e^{s\operatorname{tr}A}=1,
\]

we have

\[
\det g(s)=\det g.
\]

The finite checks give:

- maximum determinant error along the tracefree path:
  \[
  7.77\times10^{-15};
  \]
- maximum first derivative of \(U_\Lambda\) along shape directions:
  \[
  1.67\times10^{-12};
  \]
- maximum second derivative:
  \[
  3.89\times10^{-8},
  \]
  consistent with the finite-difference floor.

By contrast, along the conformal path

\[
g(s)=e^{2s}g,
\]

we have

\[
\sqrt{\det g(s)}=e^{3s}\sqrt{\det g},
\]

and the measured first derivative is at least

\[
3.00000018.
\]

Thus a pure volume term has a conformal force before reduction but no force along the exact reduced homogeneous shape manifold.

## 6. What this accomplishes

Stage 5A establishes that the fixed-volume proposal is not merely a linear trick:

1. the volume and trace-momentum conditions form an exact nonlinear second-class pair;
2. the pair removes the exact negative homogeneous conformal momentum;
3. the surviving five global shape directions have positive kinetic energy;
4. preservation of fixed volume generates the trace constraint dynamically;
5. a pure volume energy becomes a constant on the reduced homogeneous shape space.

This is meaningful progress on the nonlinear gravity pillar.

## 7. What it does not solve

Stage 5A is **not** a cosmological-constant solution.

It does not yet establish:

- a local fixed-volume rule at every lattice cell;
- closure of the inhomogeneous scalar and vector constraints;
- suppression of local vacuum-energy curvature;
- the microscopic reason that the junction network selects fixed total volume;
- gravitational self-energy as a universal source;
- nonlinear tensor interactions or strong-field solutions;
- absence of an extra scalar away from the homogeneous sector.

A global constraint can make the integrated volume term nondynamical while local stress differences and local curvature remain physical. Those local equations still have to be derived.

## 8. Relation to prior work and originality

Fixed-volume and unimodular gravity have extensive prior art. The use of

\[
\det g=1
\]

or a volume constraint is not claimed as a Phase Junction invention.

The potentially project-specific result would be a microscopic derivation in which the finite junction state space and its constraint algebra force the volume pair, rather than adding it to an otherwise standard gravitational action.

Stage 5A supplies the exact target that such a derivation must reproduce.

## 9. Next gate

Promote the determinant/trace pair to the local finite frame lattice and test the nonlinear constraint algebra through the first nontrivial order.

The next calculation must:

1. define the local nonlinear frame determinant and trace momentum;
2. retain the existing scalar and vector first-class constraints;
3. compute their brackets through quadratic order in frame amplitude;
4. include the energy carried by the frame and connection themselves;
5. verify time preservation without repeated manual projection;
6. count physical modes and reject any regenerated scalar;
7. derive the nonlinear terms from the shared junction algebra rather than inserting an Einstein–Hilbert vertex as an unexplained answer.

## 10. Reproduction

```sh
python phase_junction_network/microscopic/check_nonlinear_volume_reduction.py \
  --samples 500 \
  --seed 57 \
  --output phase_junction_network/microscopic/nonlinear_volume_reduction_results.json
```
