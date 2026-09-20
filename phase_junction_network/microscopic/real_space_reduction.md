# Full Real-Space Constrained Reduction

**Status:** propagating real-space sector passed; homogeneous conformal mode remains open  
**Date:** 2026-09-20

## 1. Purpose

The previous reduction was performed one momentum block at a time. It did not use a transverse-traceless projector, but it still used translation symmetry to diagonalize the lattice.

This calculation performs the constraint reduction on the **full periodic real-space matrices** before any spectral comparison. It asks whether the local constraints themselves remove the unwanted scalar and vector directions and leave two positive propagating modes per nonzero lattice momentum.

## 2. Real-space regulator

Use an odd \(L^3\) periodic cubic lattice and the central-difference derivative

\[
(D_i f)(x)=\frac12\left[f(x+\hat i)-f(x-\hat i)\right].
\]

Its Fourier symbol is

\[
i\widetilde k_i=i\sin\left(\frac{2\pi n_i}{L}\right).
\]

Odd lattice sizes avoid a separate Nyquist zero of the central derivative.

At each site the frame has six symmetric components and six conjugate momenta. The independent torsion-free connection has eighteen components and eighteen primary momenta.

## 3. Local real-space maps

The spatial gauge map is constructed directly from

\[
\delta h_{ij}=D_i\xi_j+D_j\xi_i.
\]

The scalar constraint is

\[
C(x)=D_iD_jh_{ij}(x)-D^2h(x).
\]

The local connection map is

\[
\Gamma^a{}_{bc}(x)
=\frac12\left[
D_bh^a{}_c(x)+D_ch^a{}_b(x)-D^ah_{bc}(x)
\right].
\]

All matrices are assembled in site space. No momentum-space polarization basis is used during the reduction.

The discrete identity

\[
C\,R=0
\]

holds to roundoff, so the scalar constraint is invariant under the vector gauge map.

## 4. Constraint count

For \(N=L^3\) sites:

- frame phase-space dimension: \(12N\);
- connection second-class constraints: \(36N\);
- vector plus scalar first-class rank on a periodic lattice:
  \[
  3N-3+N-1=4N-4.
  \]

The four rank deficiencies are the summed periodic redundancies.

The full periodic frame system therefore has

\[
12N-2(4N-4)=4N+8
\]

physical phase-space dimensions after the connection has been eliminated. This includes the twelve-dimensional homogeneous frame sector at \(k=0\).

## 5. Isolating the propagating sector

To test the finite-momentum gravitons, remove the six constant frame coordinates and six constant frame momenta. This is a boundary/global-sector choice, not a transverse-traceless projection.

The mean-zero phase space has dimension

\[
12(N-1).
\]

The first-class rank is

\[
4(N-1).
\]

After imposing the constraints and quotienting their gauge directions, the propagating physical phase space has

\[
\boxed{4(N-1)}
\]

dimensions, or

\[
\boxed{2(N-1)}
\]

configuration modes. That is exactly two modes for every nonzero lattice momentum.

## 6. Real-space symplectic reduction

Let \(F\) be the full real-space constraint matrix and \(Q_0\) the mean-zero embedding. The reduction is:

1. form the mean-zero constraints
   \[
   F_0=FQ_0;
   \]
2. extract an independent row basis of \(F_0\);
3. find a basis \(N_0\) for the full real-space constraint surface;
4. pull back the symplectic form,
   \[
   \Omega_0=N_0^TQ_0^T\Omega Q_0N_0;
   \]
5. quotient the null directions of \(\Omega_0\);
6. pull back the real-space Hamiltonian and diagonalize the resulting physical dynamical matrix.

The Fourier symbols are used only after this reduction to compare the resulting unordered frequency spectrum with the analytic lattice expectation.

## 7. Results

### \(L=3\)

\[
N=27.
\]

- mean-zero phase dimension: 312;
- first-class rank: 104;
- constraint-surface dimension: 208;
- physical phase dimension: 104;
- positive frequencies: 52;
- minimum physical Hessian eigenvalue: 0.75;
- maximum relative frequency error: \(3.63\times10^{-15}\);
- maximum growth rate: \(1.67\times10^{-15}\).

### \(L=5\)

\[
N=125.
\]

- mean-zero phase dimension: 1488;
- first-class rank: 496;
- constraint-surface dimension: 992;
- physical phase dimension: 496;
- positive frequencies: 248;
- minimum physical Hessian eigenvalue: 0.34549;
- maximum relative frequency error: \(2.20\times10^{-14}\);
- maximum growth rate: \(7.61\times10^{-15}\).

Both lattices have zero count, sign, or spectrum failures.

## 8. Homogeneous global sector

The real-space calculation exposes a point hidden by the nonzero-momentum analysis.

At \(k=0\):

- all derivative constraints vanish;
- all six constant frame coordinates have zero potential;
- the frame kinetic matrix has eigenvalues
  \[
  \left\{-\frac12,1,1,1,1,1\right\}.
  \]

Therefore the twelve-dimensional periodic global sector has Hamiltonian inertia

\[
\boxed{(1\text{ negative},6\text{ zero},5\text{ positive}).}
\]

The negative direction is the homogeneous trace momentum, the finite-volume conformal mode.

This mode is not a propagating finite-momentum ghost—the nonzero-momentum constrained sector is positive—but a fundamental periodic theory cannot simply ignore it. It needs one of:

- a fixed-volume or unimodular constraint;
- an additional global gauge redundancy;
- boundary conditions that fix the homogeneous frame;
- a nonlinear potential or junction rule that lifts it;
- a demonstration that it is not part of the physical Hilbert space.

This is directly connected to the unresolved vacuum-volume problem.

## 9. Decision

The real-space calculation closes another part of the linear bridge:

\[
\boxed{
\text{local real-space constraints}
\longrightarrow
2\text{ positive modes per nonzero momentum}
}
\]

without Fourier reduction or a TT projector.

The next work should proceed on two tracks:

1. implement the second-class connection constraints in a finite local Hilbert space and derive the reduced transfer matrix from those local variables;
2. add and test a principled global-volume rule that removes or controls the homogeneous conformal direction.

The second item is no longer optional: it is an explicit instability of the unrestricted periodic global sector.