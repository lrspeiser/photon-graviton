# Reduced Physical TT Transfer-Matrix Benchmark

**Status:** finite physical-sector benchmark passed; local microscopic embedding remains open  
**Date:** 2026-09-20

## 1. Purpose

The auxiliary-stability result shows that the independent connection cannot be quantized as 18 ordinary positive-energy oscillators with a zero bare frame block. The next logically prior step is therefore to solve the linear constraints and ask whether the **reduced physical sector itself** admits a controlled finite Hilbert-space transfer matrix.

This benchmark addresses that narrower question.

It does not derive the physical tensor modes from the local microscopic frame and connection variables. Instead, it establishes the spectrum and correlation-function target that any full local construction must reproduce.

## 2. Reduced physical Hamiltonian

For every nonzero lattice momentum, the four constraints leave two transverse-traceless coordinates,

\[
q_+,\qquad q_\times,
\]

with conjugate momenta

\[
p_+,\qquad p_\times.
\]

In units where the tensor speed and kinetic normalization are one, the physical Hamiltonian is

\[
H_{\rm TT}
=
\frac12\sum_{\lambda=+,\times}
\left[p_\lambda^2+\omega_{\mathbf k}^2q_\lambda^2\right],
\]

where the prototype cubic-lattice frequency is

\[
\omega_{\mathbf k}
=
\sqrt{4\sum_i\sin^2\left(\frac{k_i}{2}\right)}.
\]

For the lowest axial mode on an \(L^3\) periodic lattice,

\[
\omega_L=2\sin\left(\frac{\pi}{L}\right).
\]

## 3. Finite oscillator representation

Each polarization is represented by a \(d\)-level truncated oscillator. The canonical operators obey

\[
[Q,P]=i\left(I-d|d-1\rangle\langle d-1|\right).
\]

A boundary penalty

\[
H_{\rm boundary}
=\lambda\omega|d-1\rangle\langle d-1|
\]

pushes the anomalous top state above the low physical spectrum.

The two-polarization Hilbert space has dimension

\[
d^2.
\]

The benchmark used

\[
d\in\{3,4,5,6,8\},\qquad\lambda=10.
\]

## 4. Spectrum and transfer matrix

The Euclidean transfer matrix is

\[
T=e^{-a_tH_{\rm TT}}.
\]

Its first physical gap can be obtained either from the Hamiltonian spectrum or from

\[
\Delta E=-\frac1{a_t}\log\left(\frac{t_1}{t_0}\right),
\]

where \(t_0,t_1\) are the two largest transfer eigenvalues.

Across all tested lattice sizes and oscillator dimensions:

- maximum relative Hamiltonian-gap error: \(4.29\times10^{-15}\);
- maximum relative transfer-gap error: \(6.93\times10^{-15}\);
- first excitation degeneracy: exactly two;
- maximum boundary occupation among the four lowest states: \(4.75\times10^{-31}\).

The lowest gap fits

\[
\omega_L\propto L^{-0.993857},
\]

consistent with a massless linear mode.

## 5. Euclidean tensor correlator

For either polarization, the ground-state correlator is

\[
C(\tau)=\langle0|q(\tau)q(0)|0\rangle
=\frac{e^{-\omega\tau}}{2\omega}.
\]

The finite calculation evaluates the correlator spectrally and compares it with the analytic result at four positive Euclidean times.

The reflection kernel

\[
K_{ij}=C(\tau_i+\tau_j)
\]

must be positive semidefinite. The minimum eigenvalue over all tested sizes and dimensions was

\[
-2.64\times10^{-15},
\]

consistent with zero at floating-point precision. The nonzero reflection eigenvalue is positive.

## 6. What this establishes

The calculation establishes that:

1. after the constraints are solved, the two physical tensor modes have a positive finite transfer matrix;
2. a finite oscillator truncation can preserve the low tensor spectrum with negligible boundary leakage;
3. the twofold polarization degeneracy survives;
4. the finite-size gap closes as \(1/L\);
5. the Euclidean two-point function has the correct pole and a positive reflection kernel in the reduced free theory.

## 7. What it does not establish

This benchmark does not show that:

- a local finite frame–connection Hamiltonian reduces to these oscillators;
- the reduction is local;
- nonlinear constraints close;
- the connection has no hidden finite-state excitation;
- the physical transfer matrix follows from a positive local microscopic transfer matrix;
- the electromagnetic and gravitational stiffnesses share one microscopic origin.

The use of a momentum-space TT reduction is intentionally nonlocal. It is a target test, not the final microscopic construction.

## 8. Next step

Construct the same reduced transfer matrix **from the local constrained variables**, rather than inserting the TT oscillators directly.

The next implementation should:

1. classify the linear frame–connection constraints canonically;
2. eliminate the auxiliary connection and gauge variables at finite \(L\);
3. compute the reduced symplectic form without a pre-imposed TT projector;
4. truncate only after reduction;
5. compare the resulting finite spectrum and correlator with this benchmark;
6. quantify any mismatch and boundary leakage.

A successful match would establish the linear finite quantum phase. The nonlinear and parameter-closure problems would remain.
