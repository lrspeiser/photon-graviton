# Stage 3D: Exact Plaquette Diagonal Matching

**Status:** **PASS as a rejection gate.** The current pure-swap shared graph does not produce a counterterm-free relativistic common cone once all induced low-band operators are retained.  
**Date:** 2026-09-20  
**Issue:** #3  
**Executable:** [`check_plaquette_diagonal_matching.py`](check_plaquette_diagonal_matching.py)  
**Frozen output:** [`plaquette_diagonal_matching_results.json`](plaquette_diagonal_matching_results.json)

## 1. Why Stage 3D was required

Stage 3A derived only the off-diagonal electromagnetic ring process,

\[
K_A\sim20\frac{t^4}{\Delta^3}.
\]

Stage 3C then introduced a collective charging coefficient

\[
U_A=U_g=\Delta/144.
\]

However, the same microscopic swaps that generate the ring process also generate diagonal self-energies and additional operators when the Gauss-violating states are eliminated. A common-cone claim is not meaningful unless those terms are included.

Stage 3D performs that exact matching in the actual four-link spin-1 Hilbert space.

## 2. Exact microscopic plaquette

Use four spin-1 links with normalized raising operator \(U\) and electric field \(E=S_z\). The square Gauss generators are

\[
G_0=E_0+E_3,
\quad
G_1=E_1-E_0,
\quad
G_2=-E_1-E_2,
\quad
G_3=E_2-E_3.
\]

The microscopic Hamiltonian is

\[
H
=\frac{\Delta}{2}\sum_{v=0}^{3}G_v^2
-t\sum_{\ell=0}^{3}(U_\ell+U_\ell^\dagger).
\]

Set

\[
x=t/\Delta.
\]

The exact Gauss-law subspace contains three loop-flux states,

\[
|m\rangle
=|m,m,-m,-m\rangle,
\qquad
m=-1,0,+1.
\]

For every \(x\), the full \(81\times81\) Hamiltonian is diagonalized. The three lowest eigenstates are projected into the exact gauge basis, and the polar decomposition gives a canonical Hermitian three-state effective Hamiltonian.

No perturbative coefficient is inserted into the matching.

## 3. Operator basis and linked-cluster prescription

The intended Stage-6B plaquette operators act in this gauge band as

\[
H_{\rm eff}
=2U\,m^2
-K(W+W^\dagger)
+vD
+J_2(W^2+W^{\dagger2})+\cdots,
\]

where

\[
D=\operatorname{diag}(1,2,1).
\]

The exact low-band Hamiltonian determines:

- nearest-neighbor flux tunnelling \(K\);
- edge-center diagonal contrast;
- direct \(m=-1\leftrightarrow+1\) tunnelling \(J_2\);
- separation from the next microscopic band.

To distinguish the extensive electric renormalization from the plaquette-local diagonal term, the calculation independently diagonalizes one link with fixed endpoint Gauss charges. Four copies of that exact one-link self-energy define \(\delta U\). The remaining connected edge-center contrast is matched to \(vD\).

This is a linked-cluster matching convention: disconnected one-link diagrams renormalize the electric coefficient; the connected remainder belongs to the plaquette interaction.

## 4. Small-move result

For

\[
0.005\le x\le0.05,
\]

the exact finite spectra give

\[
K
\propto x^{3.99090},
\]

\[
\delta U
\propto x^{1.99695},
\]

and

\[
|v|
\propto x^{3.95847}.
\]

The asymptotic coefficients are

\[
\boxed{
K=20x^4\Delta+O(x^6\Delta),
}
\]

\[
\boxed{
\delta U=2x^2\Delta+O(x^4\Delta),
}
\]

and

\[
\boxed{
v=-\frac83x^4\Delta+O(x^6\Delta).
}
\]

Therefore

\[
\boxed{
\frac{v}{K}\longrightarrow-\frac{2}{15}.
}
\]

At \(x=0.005\), the exact ratios to these leading forms are:

| Quantity | Exact / leading |
|---|---:|
| \(K/(20x^4\Delta)\) | 0.999760 |
| \(\delta U/(2x^2\Delta)\) | 0.999919 |
| \(v/[ -(8/3)x^4\Delta]\) | 0.998979 |

The central structural result is the order mismatch:

\[
\delta U=O(x^2),
\qquad
K=O(x^4).
\]

The generated electric self-energy is parametrically larger than the photon loop stiffness.

## 5. Counterterm-free common-cone test

Retain the Stage-3C collective bare coefficient,

\[
U_{\rm bare}=\Delta/144,
\]

and include the exact swap-generated correction,

\[
U_A=\Delta/144+\delta U(x).
\]

Use the matched

\[
K(x),
\qquad
v(x),
\]

in the actual Stage-6B spin-1 transverse Hamiltonian. Compare its fitted physical speed with the current shared-frame target

\[
c_g=\sqrt{(\Delta/144)K_g(x)}.
\]

A 36-point scan over

\[
0.05\le x\le0.40
\]

finds no common-cone crossing.

The smallest speed difference is still positive:

\[
\boxed{
\min(c_\gamma-c_g)
=0.01776797\,\Delta
}
\]

at \(x=0.05\).

Within the 19 scan points that independently pass the declared approximately linear photon criteria, the minimum difference is

\[
\boxed{
0.04234674\,\Delta.
}
\]

Thus the current shared move graph does not produce a counterterm-free relativistic common cone over the controlled low-band range.

## 6. Negative control: delete the predicted self-energy

If the exact \(O(x^2)\) electric correction is simply omitted, a superficially successful root appears at

\[
\frac{t}{\Delta}=0.1395938898.
\]

At that point:

\[
\frac{v}{K}=-0.01820925,
\]

\[
\text{photon gap power}=0.989096,
\]

and the fitted photon and gravity speeds agree at

\[
0.0114203163\,\Delta.
\]

But the deleted correction is

\[
\delta U=0.0367029121\,\Delta,
\]

more than five times the Stage-3C bare collective coefficient \(\Delta/144\).

Therefore this root requires a large, unearned cancellation. It is retained as a negative control, not accepted as a solution.

## 7. Negative control: copy the self-energy into gravity

A second possible shortcut is to add the same \(\delta U\) to the gravity charging coefficient by declaration.

That produces a speed crossing near

\[
\frac{t}{\Delta}=0.199959568,
\]

but the photon branch has

\[
\text{gap power}=0.904084
\]

and linear-plus-cubic residual

\[
6.79\times10^{-3}.
\]

It fails the Stage-6B relativistic branch criteria. Matching one velocity is not enough.

Moreover, the present frame mediator has not independently derived this gravity self-energy, so copying it would already be an unjustified sector identification.

## 8. Additional generated operator

The exact low band also contains direct

\[
|{-1}\rangle\leftrightarrow|{+1}\rangle
\]

mixing, equivalent to

\[
J_2(W^2+W^{\dagger2}).
\]

It is eighth order at small \(x\), but grows to

\[
J_2/K=0.07353
\]

at \(x=0.22\), and

\[
J_2/K=0.25422
\]

at \(x=0.40\).

At the upper scan edge, the three-state gauge band remains separated from the next microscopic band by

\[
3.91K,
\]

so the extra operator is not merely an artifact of band collapse. A complete shared microscopic theory must retain or dynamically suppress it.

## 9. Decision

The present pure-swap microscopic graph is rejected as the complete electromagnetic–gravity unification.

It correctly generates ring exchange, but it simultaneously predicts:

1. an \(O(t^2/\Delta)\) electric self-energy;
2. a linked diagonal plaquette term with negative small-\(x\) coefficient;
3. a higher flux-jump operator;
4. no counterterm-free relativistic common cone in the controlled scan.

The failure is specific. It does not refute the broader Phase Junction idea, the finite photon phase, or the finite frame regulator. It refutes the claim that a Hamiltonian containing only the currently declared off-diagonal elementary swaps and Gauss penalties is already the shared microscopic theory.

## 10. Originality significance

This is progress toward an alternative theory rather than away from it.

Known ring-exchange and emergent-gravity methods can be assembled with independently selected coefficients. Stage 3D prevents us from mistaking that assembly for a new theory. The project-specific mechanism must now explain the cancellation or sharing of the second-order diagonal terms as well as the curvature terms.

Any acceptable replacement must derive its diagonal partners from symmetry. Adding a fitted electric counterterm only in the photon sector would destroy the shared-mechanism claim.

## 11. Next gate

Enumerate **symmetry-complete elementary junction Hamiltonians**. Each candidate must include all diagonal or seagull terms required by its local algebra, not only the off-diagonal swap.

For every candidate:

1. derive the photon and frame low-band Hamiltonians through at least fourth order;
2. include all generated electric, diagonal, double-flux, and frame operators;
3. identify a symmetry that cancels or relates the \(O(x^2)\) charging corrections;
4. apply the same symmetry to photon, frame, matter, and companion channels;
5. reject any candidate requiring a photon-only or gravity-only counterterm;
6. rerun the common-cone and finite-spectrum gates.

The preferred next candidate is a positive local junction generator whose diagonal term is algebraically tied to each allowed swap, analogous to a finite graph-Laplacian move. Whether that structure can retain the Stage-6B \(z\simeq1\) region must be calculated, not assumed.

## 12. Reproduction

```sh
python phase_junction_network/microscopic/check_plaquette_diagonal_matching.py \
  --output phase_junction_network/microscopic/plaquette_diagonal_matching_results.json
```

A reduced regression run is available with `--quick`.
