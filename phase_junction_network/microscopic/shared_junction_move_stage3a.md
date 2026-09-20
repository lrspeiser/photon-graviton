# Stage 3A: Shared Junction Move for Electromagnetism and Gravity

**Status:** finite microscopic prototype pass; full issue #3 remains open  
**Date:** 2026-09-20  
**Issue:** #3  
**Executable:** [`check_shared_junction_move.py`](check_shared_junction_move.py)  
**Frozen output:** [`shared_junction_move_results.json`](shared_junction_move_results.json)

## 1. Why this calculation matters

Quantum links, ring exchange, finite helicity-2 models, first-order frame/connection gravity, and enlarged gauge–gravity formulations all have substantial prior art. The Phase Junction program is not original merely because it uses those tools.

The project-specific question is narrower:

> Can the same elementary matter–geometry junction move and the same virtual-state penalty generate both the electromagnetic and gravitational stiffnesses, so their relative normalization is calculated rather than selected independently?

Stage 3A supplies the first executable model aimed directly at that question.

## 2. Shared microscopic inputs

The prototype has only two microscopic energy quantities:

\[
 t=\text{one elementary junction-swap amplitude},
\qquad
 \Delta=\text{one virtual defect penalty}.
\]

Define

\[
x=\frac{t}{\Delta}.
\]

The local imbalance penalty is assumed to be one common junction Casimir,

\[
H_{\rm charge}
=\frac{\Delta}{2}
\left(E_A^2+\Pi_g^2\right),
\]

with equal generator normalization in the electromagnetic and frame sectors. This equality is an explicit hypothesis of the prototype, not a derived fact.

## 3. Electromagnetic virtual process

Use one four-link plaquette. The two low-energy states are

\[
|0000\rangle,
\qquad
|1111\rangle,
\]

representing no loop and one completed oriented loop. A partial flipped string costs \(\Delta\) for every connected component.

One elementary move toggles one link with matrix element \(-t\). A completed plaquette requires four moves. The 24 possible orderings divide into:

- 16 orderings with denominators \((\Delta,\Delta,\Delta)\);
- 8 orderings with denominators \((\Delta,2\Delta,\Delta)\).

Therefore the leading off-diagonal plaquette amplitude is

\[
\boxed{
K_A
=20\frac{t^4}{\Delta^3}
+O\!\left(\frac{t^6}{\Delta^5}\right).
}
\]

This is standard degenerate perturbation methodology. The project-specific content is that the same \(t\) and \(\Delta\) are also used below for the frame sector.

## 4. Frame–connection virtual process

For one physical frame mode, use two low-energy frame configurations connected through one gapped connection mediator:

\[
|L\rangle
\xleftrightarrow{\ t\ }
|C\rangle
\xleftrightarrow{\ t\ }
|R\rangle,
\qquad
E_C=\Delta.
\]

The leading effective frame tunnelling is

\[
\boxed{
K_g
=\frac{t^2}{\Delta}
+O\!\left(\frac{t^4}{\Delta^3}\right).
}
\]

The exact three-state Hamiltonian is diagonalized rather than relying only on the second-order formula.

## 5. Exact scaling checks

For

\[
x=0.01,0.015,0.02,0.03,0.04,0.05,
\]

the exact finite spectra give

\[
K_A\propto x^{3.99197},
\qquad
K_g\propto x^{1.99719}.
\]

The fitted powers reproduce the expected fourth- and second-order processes.

At \(x=0.01\):

\[
\frac{K_A^{\rm exact}}{20x^4}=0.999440,
\qquad
\frac{K_g^{\rm exact}}{x^2}=0.999800.
\]

The executable also enumerates all 24 electromagnetic orderings and reproduces the coefficient 20 exactly.

## 6. Reduced parameter count

Before the shared-move hypothesis, the effective theory contains four independent stiffnesses:

\[
U_A,\ K_A,\ U_g,\ K_g.
\]

Under the explicit common-Casimir hypothesis,

\[
U_A=U_g=\Delta,
\]

while the exact finite Hamiltonians make both curvature stiffnesses functions of the same ratio \(x=t/\Delta\).

Thus the four effective quantities become functions of:

- one overall scale \(\Delta\);
- one dimensionless ratio \(x\).

At leading order,

\[
\frac{c_A}{c_g}
=\sqrt{\frac{K_A}{K_g}}
\simeq\sqrt{20}\,x,
\]

and, in this common normalization,

\[
\frac{Z_g}{Z_A}
=\sqrt{\frac{K_A}{K_g}}
\simeq\sqrt{20}\,x.
\]

This is a derived correlation within the prototype, not yet a physical prediction.

## 7. Exact common-cone solution

The exact finite spectra, rather than the leading perturbative formulas, are used to solve

\[
K_A(x)=K_g(x).
\]

The result is

\[
\boxed{
\frac{t}{\Delta}
=0.27909563384302527.
}
\]

For comparison, leading perturbation theory gives

\[
\frac{t}{\Delta}=\frac1{\sqrt{20}}
=0.22360679774997896.
\]

The exact root is 24.8% higher, showing why the finite calculation matters at the common-cone point.

At the exact root,

\[
\frac{K_A}{\Delta}
=
\frac{K_g}{\Delta}
=0.0685077493800112.
\]

With the shared charging normalization,

\[
U_A=U_g,
\qquad
K_A=K_g,
\]

so the prototype has equal bare cones and

\[
Z_g/Z_A=1
\]

in its own normalization.

The next unwanted band remains separated by

\[
9.93K
\]

in the electromagnetic graph and

\[
16.60K
\]

in the frame graph.

## 8. Negative controls

The claimed correlation disappears if any shared assumption is removed.

### Independent swap amplitudes

With \(t_A\) and \(t_g\),

\[
\frac{K_A}{K_g}
\simeq
20\frac{t_A^4}{\Delta^2t_g^2},
\]

which is arbitrary.

### Independent virtual penalties

With \(\Delta_A\) and \(\Delta_g\),

\[
\frac{K_A}{K_g}
\simeq
20\frac{t_A^4\Delta_g}{\Delta_A^3t_g^2},
\]

which is again arbitrary.

### Unequal charging Casimirs

Even if \(K_A=K_g\), allowing \(U_A/U_g\) to float gives

\[
\frac{Z_g}{Z_A}
=\sqrt{\frac{U_g}{U_A}},
\]

so the impedance relation is lost.

These controls are essential: the result comes from the shared microscopic assumptions, not algebraic relabeling.

## 9. What is original here—and what is not

### Borrowed methods

- finite quantum links and plaquette moves;
- degenerate perturbation theory and virtual ring exchange;
- an auxiliary connection mediator;
- low-energy elimination of gapped states;
- a common-cone matching condition.

Those methods are credited in [`../novelty_boundary.md`](../novelty_boundary.md).

### Project-specific result

The project-specific result is the explicit finite model in which:

1. one amplitude \(t\) and one penalty \(\Delta\) generate both sectors;
2. the sectors arise at different perturbative orders;
3. exact diagonalization fixes the common-cone ratio;
4. the parameter count falls from four effective stiffnesses to one scale plus one ratio;
5. negative controls show that the relation disappears when the shared junction assumption is removed.

A limited literature audit did not identify this exact combined finite graph and relation. That is not proof of formal novelty.

## 10. Claim boundary

Stage 3A establishes a useful **mechanism prototype**, not issue-#3 closure.

It does not yet show that:

- the actual three-dimensional spin-1 photon Hamiltonian and dressed-frame Hamiltonian descend from this exact move graph;
- the common Casimir is enforced by a fundamental symmetry;
- physical normalization gives \(Z_g/Z_A=1\);
- the observed \(\alpha\) or \(G\) follows;
- the matter-localization parameters, companion hopping, conversion, capture, or binding use the same \(t\) and \(\Delta\);
- the relation is radiatively stable;
- nonlinear gravity preserves it;
- no closer prior construction exists.

## 11. Next gate

Embed the shared mediator into the actual finite structures already committed:

1. the Stage-6B spin-1 detuned photon link;
2. the finite dressed-frame Weyl algebra;
3. the finite connection lock;
4. the minimal \(|q|=1\) matter endpoint;
5. the neutral companion qudit.

Then perform a full Schrieffer–Wolff or exact low-band calculation and ask whether the common-cone root and impedance correlation survive **without** sector-specific move amplitudes, penalties, or clock rescalings.

Failure of that embedding would reject this shared-move prototype while preserving the separate finite electromagnetic and gravity regulators.