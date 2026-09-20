# Stage 3B: Embedding the Shared Move in the Actual Photon Phase

**Status:** diagnostic pass; the naïve one-hinge shared-Casimir embedding is rejected  
**Date:** 2026-09-20  
**Issue:** #3  
**Executable:** [`check_shared_move_embedding.py`](check_shared_move_embedding.py)  
**Frozen output:** [`shared_move_embedding_results.json`](shared_move_embedding_results.json)

## 1. Why this stage is necessary

Stage 3A produced a promising finite toy relation from one swap amplitude \(t\) and one virtual penalty \(\Delta\):

\[
K_A\sim20\frac{t^4}{\Delta^3},
\qquad
K_g\sim\frac{t^2}{\Delta}.
\]

Under the extra hypothesis

\[
U_A=U_g=\Delta,
\]

the toy common-cone equation fixed \(t/\Delta\).

That result is not useful unless it embeds in the actual finite photon and gravity regulators already committed. Stage 3B performs the first such compatibility test. It deliberately treats incompatibility as a scientific result rather than adding a fitting coefficient after the fact.

## 2. Mapping Stage 3A into Stage 6B

At the exact Stage-3A common-cone root,

\[
\frac{t}{\Delta}=0.27909563384302527,
\]

and

\[
\frac{K_A}{\Delta}
=0.06850774938001117.
\]

The Stage-6B electromagnetic Hamiltonian uses

\[
H_A/t_{\rm plaq}
=\frac{u/t_{\rm plaq}}2\sum E^2
-\sum(W+W^\dagger)
+\frac{v}{t_{\rm plaq}}\sum D.
\]

The naïve Stage-3A identifications are

\[
t_{\rm plaq}=K_A,
\qquad
u=U_A=\Delta.
\]

They imply

\[
\boxed{
\frac{u}{t_{\rm plaq}}
=rac{\Delta}{K_A}
=14.5968888053.
}
\]

The Stage-6B linearly dispersing region was demonstrated only for

\[
0\le\frac{u}{t_{\rm plaq}}\le0.2.
\]

Thus the naïve common-Casimir point lies a factor

\[
72.9844
\]

above the demonstrated upper edge.

## 3. Exact finite photon test at the naïve ratio

The spin-1 axial finite blocks were diagonalized at

\[
\frac{u}{t_{\rm plaq}}=14.5968888053,
\qquad
\frac{v}{t_{\rm plaq}}=0.2,0.4,0.6,
\]

on \(L=6,8,10\).

| \(v/t_{\rm plaq}\) | fitted gap power | \(L=10\) first gap |
|---:|---:|---:|
| 0.2 | 0.138893 | 10.601739 |
| 0.4 | 0.154351 | 9.835648 |
| 0.6 | 0.171816 | 9.071740 |

The gaps remain order one in plaquette units and barely decrease with momentum. They do not show

\[
\Delta_\gamma\propto|\widehat k|.
\]

Therefore

\[
\boxed{
U_A=U_g=\Delta
\text{ is incompatible with the demonstrated Stage-6B photon phase.}
}
\]

This rejects the simplest Stage-3A normalization while preserving the more general shared-move idea.

## 4. What the failure teaches us

A viable shared microscopic model needs a **collective charging or capacitance factor** \(r\):

\[
U_A=U_g=r\Delta,
\]

while the virtual curvature processes retain the same elementary \(t\) and \(\Delta\).

To reach the upper edge of the demonstrated photon region,

\[
\boxed{r\le0.01370155.}
\]

To reach the representative Stage-6B point \(u/t_{\rm plaq}=0.1\),

\[
\boxed{r=0.006850775.}
\]

If—and only if—a derived collective model later gives \(r=1/N_{\rm eff}\), these numbers correspond to

\[
N_{\rm eff}\gtrsim72.98
\]

and

\[
N_{\rm eff}\simeq145.97,
\]

respectively. These are target scales, not a derivation. Merely declaring that many hinges does not produce the required capacitance matrix.

## 5. The existing gravity normalization also fails the shared test

The current finite gravity regulator uses a self-dual clock normalization with

\[
U_g/K_g=1.
\]

At the Stage-3A curvature value, that corresponds to

\[
r_g=K_g/\Delta=0.0685077494.
\]

The representative photon point requires

\[
r_A=0.00685077494.
\]

Their ratio is exactly ten in the current conventions.

Therefore the photon and gravity regulators cannot both be retained with their independently chosen clock normalizations and still be called a shared microscopic theory. The gravity regulator must be rerun at the same **derived** charging ratio as the photon sector.

Independent \(r_A\) and \(r_g\) would trivially repair both spectra, but would restore the arbitrary impedance ratio that issue #3 exists to eliminate.

## 6. Originality significance

This negative result is progress toward the alternative gravity theory because it prevents an apparent unification by notation.

Known methods already allow us to build:

- a finite quantum-link photon sector;
- a finite helicity-2 regulator;
- a first-order frame/connection system.

The theory becomes genuinely distinct only if the same microscopic matter–geometry junction predicts their relative charging and curvature scales. Stage 3B shows that the simplest one-hinge equality does not work.

The remaining candidate original mechanism is more specific:

> a finite collective junction bundle whose capacitance/constraint matrix suppresses the shared charging mode while the same elementary local moves generate both electromagnetic loops and frame curvature.

That mechanism is not yet derived and should not be claimed as a result.

## 7. Acceptance accounting

| Diagnostic | Result |
|---|---|
| Map Stage-3A coefficients into Stage 6B | Pass |
| Test the exact finite photon blocks at the implied ratio | Pass |
| Naïve ratio lies outside demonstrated photon region | Yes |
| Naïve ratio has \(z\approx1\) | No; powers 0.139–0.172 |
| Naïve ratio remains gapped at \(L=10\) | Yes; gaps 9.07–10.60 |
| Quantify required common suppression | Pass |
| Test current gravity clock compatibility | Fail; factor-ten mismatch at representative photon point |
| Preserve one shared ratio rather than introduce sector-specific fits | Required next gate |

## 8. Claim boundary

Stage 3B establishes:

1. the naïve common-Casimir Stage-3A embedding is incompatible with the demonstrated photon phase;
2. a common charging suppression of order \(10^{-2}\) or smaller is required;
3. the current photon and gravity regulator clocks were chosen independently and cannot both survive unchanged;
4. any successful shared theory must derive one collective charging matrix and apply it to both sectors.

It does not establish:

- the microscopic origin of \(r\);
- that the full three-dimensional photon phase ends at \(u/t=0.2\);
- that \(r=1/N\);
- a revised shared model that passes both sectors;
- physical \(\alpha\), \(G\), or a common interacting cone.

## 9. Next gate

Construct an explicit finite collective junction bundle with a capacitance/constraint matrix \(C\). The next calculation must:

1. derive the uniform and relative charging eigenvalues of \(C^{-1}\);
2. obtain \(r\) from those eigenvalues without fitting photon or gravity data;
3. retain the same elementary swap amplitude and virtual penalty in both sectors;
4. rerun the Stage-6B photon spectrum at the derived \(r\);
5. rerun the finite dressed-frame Hamiltonian at the same \(r\);
6. reject the branch if separate temporal or charging normalizations remain necessary;
7. extend the same matrix to the companion and matter endpoint rather than introducing new sector clocks.

That collective-capacitance derivation, not another standalone photon or graviton calculation, is now the critical path for the alternative theory.