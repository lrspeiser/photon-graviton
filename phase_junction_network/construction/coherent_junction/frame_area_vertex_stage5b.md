# Stage 5B: Coherent Frame–Area Curvature Vertex

**Status:** **PASS at finite coherent-parent and classical geometry scope.** The current four-stage coherent junction produces a first-order frame–area/curvature interaction with a fixed absolute coefficient. The same finite Clifford area satisfies the three-dimensional Palatini identity.  
**Date:** 2026-09-21  
**Issue:** #5  
**Executable:** [`check_frame_area_vertex.py`](check_frame_area_vertex.py)  
**Frozen output:** [`frame_area_vertex_results.json`](frame_area_vertex_results.json)

## 1. Missing structure addressed

The previous local nonlinear audit established two related facts:

1. holonomy multiplication supplies the connection commutator at quadratic order;
2. a Hermitian completed **unitary** plaquette by itself begins at curvature squared.

It therefore did not produce the first-order frame/area insertion required by a local gravitational scalar constraint.

The current coherent parent already fixes the complete return resolvent of a four-stage junction and its absolute gain. Stage 5B couples its entrance port to a Hermitian frame-area matrix rather than completing the already-eliminated transfer a second time.

## 2. Exact finite parent

For holonomy \(W\), the heavy chiral ring has the exact return block

\[
G_{00}(W)=\alpha^{-1/2}(I-\rho W)^{-1},
\]

with

\[
x=0.13554178509861228,
\qquad
\rho=0.0003133951506985591.
\]

The physical and identity-reference rings are retained. The left low port carries a Hermitian area matrix \(B\); the right port carries the identity. Exact heavy-state elimination gives

\[
\boxed{
H_{LR}(0)
=-\frac{x^2}{2}
B\,[G_{00}(W)-G_{00}(I)]^\dagger.
}
\]

The diagonal low-port Schur blocks vanish by the chiral structure. The reference cancellation is inside the Hamiltonian and is not a subtraction of an observed vacuum energy.

Project the low ports onto the transport quadrature

\[
|+y\rangle=\frac{|L\rangle+i|R\rangle}{\sqrt2}.
\]

The resulting internal operator is

\[
H_{\rm area}=\frac{i}{2}(H_{LR}-H_{LR}^\dagger).
\]

For

\[
W=e^{i\epsilon F},
\]

its first derivative is

\[
\boxed{
\left.\frac{dH_{\rm area}}{d\epsilon}\right|_{0}
=-\frac{\kappa}{2}\{B,F\},
}
\]

where the absolute gain is

\[
\boxed{
\kappa
=\frac{x^2\rho}{2\sqrt\alpha(1-\rho)^2}
=2.854012469780666\times10^{-6}.
}
\]

This coefficient includes the solved return amplitude. It is not divided out to normalize the curvature response.

## 3. Finite checks

Across 20 random four-component holonomies and Hermitian areas:

- maximum exact Schur-block residual: \(1.11\times10^{-17}\);
- maximum diagonal low-block norm: \(0\);
- flat-holonomy residual: \(0\);
- local conjugation-covariance residual: \(3.81\times10^{-18}\).

A central derivative of the finite parent converges to the anticommutator vertex with error power

\[
2.000010.
\]

Thus the first-order vertex is an actual derivative of the finite coherent Hamiltonian, not a fitted continuum interaction.

## 4. Palatini identity from the same finite area

Let \(e_i{}^a\) be a positive coframe, \(A_i{}^a\) its torsion-free spin connection, and

\[
F_{ij}{}^a
=\partial_iA_j{}^a-\partial_jA_i{}^a
-\epsilon^{abc}A_i{}^bA_j{}^c.
\]

The finite Clifford area is

\[
B^{ij}
=\sqrt g\,\frac{[\gamma^i,\gamma^j]}{2i}.
\]

Using the committed spin generators, the executable verifies

\[
\boxed{
\sum_{i<j}\operatorname{tr}(B^{ij}F_{ij})
=\epsilon^{ijk}e_i{}^aF_{jk}{}^a
=\sqrt g\,R.
}
\]

The smooth spectral control closes with relative errors \(2.38\times10^{-14}\) for the Palatini/metric identity and \(2.27\times10^{-16}\) for the finite Clifford-area trace.

A local central-difference regulator approaches the identity as

\[
L^{-1.98339},
\]

consistent with the expected second-order lattice error.

## 5. Scalar-curvature coefficient

The area block is four dimensional and the declared local energy convention is \(\operatorname{tr}H/4\). Therefore the coherent parent supplies the candidate scalar-curvature density

\[
\boxed{
\mathcal H_R=-b\sqrt g\,R,
\qquad
b=\frac{\kappa}{4}
=7.135031174451665\times10^{-7}.
}
\]

No Einstein–Hilbert coefficient was inserted to repair a bracket or match a measured speed.

## 6. What Stage 5B establishes

Stage 5B establishes:

1. an explicit finite parent for the missing first-order frame/area vertex;
2. exact heavy-state elimination and flat-reference cancellation;
3. a fixed absolute curvature amplitude and its frequency-resolvent ancestry;
4. local frame conjugation covariance;
5. the Palatini scalar-curvature identity for the same finite Clifford area.

## 7. Claim boundary

Stage 5B does **not** establish:

- uniqueness of the area-port layout;
- a complete finite quantum nonlinear gravity Hamiltonian;
- exact finite-lattice diffeomorphism closure;
- the temporal kinetic residue from the same coherent parent;
- a quantum anomaly-free scalar constraint;
- strong-field solutions, Newton’s constant, or empirical validity.

The area insertion and equal physical/reference splitter remain explicit architecture choices. The calculation shows what follows from them; it does not prove nature selects them.

## 8. Next gate

Use \(a=U_{\rm frame}/2\) and \(b=\kappa/4\) in the nonlinear scalar constraint and test:

1. the full \(D-D\), \(D-H\), and \(H-H\) brackets;
2. unprojected nonlinear time evolution;
3. physical mode counting at a nonlinear background;
4. the classification of local determinant/trace conditions;
5. exact versus emergent closure of the finite spatial regulator.

That gate is implemented by `validation/check_nonlinear_constraint_closure.py`.
