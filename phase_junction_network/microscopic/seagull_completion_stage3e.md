# Stage 3E: Symmetry-Complete Local Seagull Candidate

**Status:** **PASS at candidate scope.** A unit finite-link seagull cancels the dangerous second-order isolated-move self-energy while leaving the photon ring and frame superexchange intact. One relativistic shared-cone root survives finite tests. The microscopic symmetry enforcing this completion and its matter/companion extension remain open.  
**Date:** 2026-09-20  
**Issue:** #3  
**Executable:** [`check_seagull_completion.py`](check_seagull_completion.py)  
**Frozen output:** [`seagull_completion_results.json`](seagull_completion_results.json)

## 1. Motivation

Stage 3D rejected the pure-swap shared graph because the same elementary moves generated

\[
\delta U=2x^2\Delta+\cdots
\]

while the photon ring began only at

\[
K_A=20x^4\Delta+\cdots.
\]

Deleting that self-energy by hand produced a spurious common-cone root. The next acceptable candidate had to cancel it through a local algebraic partner shared by all sectors, not through a photon-only counterterm.

## 2. Candidate local completion

For every finite spin-1 link move, pair the swap with its finite-link anticommutator:

\[
\boxed{
H_A/\Delta
=
\frac12\sum_vG_v^2
-x\sum_\ell(U_\ell+U_\ell^\dagger)
+x^2\sum_\ell
(U_\ell^\dagger U_\ell+U_\ell U_\ell^\dagger),
}
\]

where

\[
x=t/\Delta.
\]

For one isolated finite link,

\[
U^\dagger U+UU^\dagger
=
\operatorname{diag}(1,2,1).
\]

This is exactly the number of allowed up/down elementary moves from each electric state. The coefficient one is the unique local value that cancels the second-order self-energy of each individual move.

The frame mediator receives the same completion:

\[
\boxed{
H_g/\Delta
=
\begin{pmatrix}
x^2&-x&0\\
-x&1&-x\\
0&-x&x^2
\end{pmatrix}.
}
\]

The endpoint seagulls cancel the two diagonal self-energies, while the cross-mediator process between the two frame states remains.

This structure is analogous to a local finite graph-Laplacian or seagull completion. That methodology is not claimed as original. The project-specific test is whether the same completion can normalize the finite photon and frame sectors together.

## 3. Why the coefficient is not being fitted

At

\[
x=0.005,
\]

the residual electric correction divided by \(x^2\) is:

| Seagull coefficient | \(\delta U/x^2\) |
|---:|---:|
| 0.8 | +0.399958 |
| 0.9 | +0.199973 |
| 1.0 | -0.0000125 |
| 1.1 | -0.199998 |
| 1.2 | -0.399983 |

Only the unit coefficient cancels the leading second-order term. The ring amplitude remains unchanged to the displayed precision.

Thus the candidate is not chosen by fitting the common-cone result. It is fixed by the local move algebra and then tested.

## 4. Low-energy coefficients

Exact diagonalization of the full \(81\)-state spin-1 plaquette gives

\[
K_A\propto x^{3.99087},
\qquad
\delta U\propto x^{3.99530},
\qquad
v\propto x^{3.99587}.
\]

The asymptotic forms are

\[
K_A=20x^4\Delta+\cdots,
\]

\[
\delta U=-\frac12x^4\Delta+\cdots,
\]

\[
v=-\frac83x^4\Delta+\cdots.
\]

The dangerous \(O(x^2)\) term is absent, but the linked photon ring and diagonal plaquette operators survive.

The exact frame coupling scales as

\[
K_g\propto x^{1.99906},
\]

so the desired second-order frame superexchange also survives.

## 5. Two speed crossings and the relativistic gate

Use the Stage-3C common collective coefficient

\[
U_{\rm bare}=\Delta/144
\]

plus the small residual \(\delta U(x)\). The exact Stage-6B photon blocks and the completed frame mediator produce two speed crossings.

### Nonrelativistic control

The first crossing occurs at

\[
\frac{t}{\Delta}=0.092027515.
\]

It has

\[
\text{photon gap power}=0.456107
\]

and fit residual

\[
2.92\times10^{-2}.
\]

It is rejected. Matching velocities alone is not sufficient.

### Physical candidate root

The second crossing is

\[
\boxed{
\frac{t}{\Delta}=0.1355417851.
}
\]

At this point:

\[
\frac{K_A}{\Delta}=0.00560619757,
\]

\[
\frac{K_g}{\Delta}=0.01805149687,
\]

\[
\frac{U}{\Delta}=0.00678991375,
\]

\[
rac{v}{K_A}=-0.129245131,
\]

and

\[
\frac{J_2}{K_A}=0.01385122.
\]

The finite photon fit gives

\[
\boxed{
\text{gap power}=0.977345,
}
\]

with maximum linear-plus-cubic residual

\[
2.91\times10^{-4}.
\]

The photon and harmonic frame speeds are

\[
c_\gamma=0.0110710481371\,\Delta,
\]

\[
c_g=0.0110710481372\,\Delta,
\]

with fractional mismatch below numerical root tolerance.

The microscopic gauge band is separated from the next band by

\[
126.68K_A.
\]

## 6. Exact periodic-cube anchor

The complete spin-1 zero-charge, zero-winding \(2^3\) gauge sector has

\[
146{,}327
\]

states. At the candidate root, the lowest gaps are

\[
0,
3.59702,
3.59702,
3.97913,
4.46820,
4.73153,
4.89235,
5.16249
\]

in plaquette units.

The first gap is twofold within the computed window. The longitudinal electric operator has exactly zero norm, while the scalar electric-energy operator has only

\[
2.05\times10^{-22}
\]

of its spectral weight in that low multiplet.

The ground-state overlap with the RK equal-amplitude state is

\[
0.295422.
\]

Thus the candidate is substantially detuned from the RK anchor while remaining in the same finite gauge component.

## 7. Finite gravity regulator at the same coefficients

The compact frame clock is rerun using the same physical

\[
U=0.00678991375\Delta
\]

and

\[
K_g=0.01805149687\Delta.
\]

Results include:

| Local prime \(p\) | gap power | fitted \(c_g\) | error from harmonic target | max fit residual |
|---:|---:|---:|---:|---:|
| 5 | 0.400004 | 0.0196888 | +77.84% | 0.02348 |
| 7 | 1.062506 | 0.0084716 | -23.48% | 0.00126 |
| 11 | 1.032903 | 0.0091169 | -17.65% | 0.000439 |
| 13 | 1.028623 | 0.0093414 | -15.62% | 0.000396 |
| 17 | 1.024818 | 0.0095606 | -13.64% | 0.000355 |
| 23 | 1.022659 | 0.0096935 | -12.44% | 0.000332 |

The coarse \(p=5\) regulator fails. For \(p\ge7\), the branch is linear and the speed moves monotonically toward the harmonic common-cone target as the local dimension increases.

Finite-\(p\) equality has not been demonstrated; the evidence is convergence toward the shared continuum normalization.

## 8. What Stage 3E establishes

The seagull-completed candidate passes the current finite architecture gates:

1. the leading photon self-energy is cancelled by a coefficient fixed before spectral testing;
2. photon ring exchange survives;
3. frame superexchange survives;
4. a finite relativistic photon branch and frame target share one speed;
5. the exact full cube rejects longitudinal and scalar contamination;
6. the compact frame regulator remains linear at the same physical coefficients for adequate local dimension.

This repairs the specific Stage-3D failure without adding separate photon and gravity clock factors.

## 9. What remains open

The central theoretical issue is now the origin of the completion.

Stage 3E does not yet prove:

- a fundamental symmetry that requires the seagull partner;
- that the same factorization works for the full dressed-frame, matter, and companion operators;
- uniqueness among all symmetry-complete local Hamiltonians;
- a stable infinite-volume interacting phase;
- finite Ward identities or vacuum polarization;
- nonlinear gravitational closure;
- radiative stability of the cancellation;
- observed \(\alpha\), \(G\), particle masses, or formal novelty.

A cancellation that is not symmetry-protected could be destroyed by interactions or renormalization.

## 10. Next gate

Construct the completed moves directly in the already committed finite Hilbert spaces:

1. replace each dressed-frame elementary move by its finite algebraic completion;
2. construct the analogous charged-endpoint and companion terms;
3. verify exact Gauss and frame constraints;
4. derive the complete low-band operator inventory rather than retaining only the desired terms;
5. test whether one factorization symmetry protects the cancellation in all sectors;
6. compute the first interacting Ward identity and one-loop/finite-volume stability test;
7. reject the candidate if sector-specific seagulls or counterterms are required.

The next claim threshold is not another matched free spectrum. It is one shared finite Hamiltonian whose algebra fixes every off-diagonal move and diagonal partner across photon, frame, matter, and companion sectors.

## 11. Reproduction

```sh
python phase_junction_network/microscopic/check_seagull_completion.py \
  --output phase_junction_network/microscopic/seagull_completion_results.json
```

A reduced regression run is available with `--quick`.
