# Stage 3C: Collective Junction Capacitance and a Shared Cone

**Status:** finite structural prototype pass; one common charging coefficient is derived, while the plaquette diagonal coefficient and full microscopic embedding remain open  
**Date:** 2026-09-20  
**Issue:** #3  
**Executable:** [`check_collective_capacitance.py`](check_collective_capacitance.py)  
**Frozen output:** [`collective_capacitance_results.json`](collective_capacitance_results.json)

## 1. Why this stage exists

Stage 3B rejected the isolated-hinge normalization

\[
U_A=U_g=\Delta.
\]

Mapped into the actual Stage-6B photon regulator, that choice gave

\[
u/t_{\rm plaq}=14.5969,
\]

and a strongly gapped rather than linearly dispersing photon branch.

The next allowed move was not to assign separate photon and gravity clock factors. It was to derive one **collective** charging coefficient from the finite matter–geometry architecture and apply it to both sectors.

## 2. Structural channel count

The elementary periodic \(2^3\) cube already used in the finite electromagnetic calculation has two independent counts of the same microscopic channel number:

\[
8\ \text{vertices}
\times
18\ \text{torsion-free connection components}
=144,
\]

and

\[
24\ \text{positive-oriented links}
\times
6\ \text{symmetric frame components}
=144.
\]

This equality is not obtained from the desired photon or gravity speed. It follows from the existing finite cube and frame/connection bookkeeping.

The Stage-3C hypothesis is that these 144 microscopic matter–geometry hinge channels share one locked charging phase while retaining their orientation and tensor labels in the curvature operators.

This is a candidate architecture, not a proof that 144 is unique or fundamental.

## 3. Explicit constraint matrix

Index the channels by a cube vertex \(v\in Q_3\) and an internal connection channel \(a=1,\ldots,18\).

The relative-phase lock is the Cartesian-product graph Laplacian

\[
L_{\rm bundle}
=
L_{Q_3}\otimes I_{18}
+
I_8\otimes L_{K_{18}}.
\]

The first term locks the same internal channel across neighboring cube vertices. The second locks the 18 channels inside one local junction.

Its spectrum is

\[
\begin{array}{c|rrrrrrrr}
\lambda&0&2&4&6&18&20&22&24\\
\text{multiplicity}&1&3&3&1&17&51&51&17.
\end{array}
\]

Therefore

\[
\operatorname{rank}L_{\rm bundle}=143,
\]

and there is exactly one uniform mode.

The relative-mode gap is 2 in the declared lock normalization.

### Topology control

Replacing the complete internal graph by an 18-cycle retains:

- rank 143;
- one uniform zero mode;
- the same common charging coefficient;

while changing the relative-mode gap to

\[
0.1206147584.
\]

Thus the \(1/N\) charging result requires a connected lock but not the all-to-all internal topology.

### Negative controls

- spatial locking alone leaves 18 independent clocks;
- internal locking alone leaves 8 independent clocks;
- both together leave one.

The shared normalization is therefore lost when either part of the finite constraint architecture is removed.

## 4. Derived charging coefficient

Give every microscopic hinge the same capacitance

\[
C_0=1/\Delta.
\]

On the locked subspace,

\[
\theta_i=\Theta
\qquad(i=1,\ldots,N),
\]

so the kinetic term is

\[
\frac12\sum_i C_0\dot\theta_i^2
=
\frac{NC_0}{2}\dot\Theta^2.
\]

The conjugate collective charge therefore has Hamiltonian

\[
H_{\rm charge}
=
\frac{\Delta}{2N}P_\Theta^2.
\]

For \(N=144\),

\[
\boxed{
U_A=U_g=r\Delta,
\qquad
r=\frac1{144}=0.00694444444444.
}
\]

This is the first common charging suppression in the project that is derived from a declared finite matrix rather than selected from either spectrum.

## 5. Common-cone calculation using the actual photon spectrum

Stage 3A provides exact finite functions

\[
K_A(x),\qquad K_g(x),
\qquad x=t/\Delta.
\]

Stage 3C no longer sets \(K_A=K_g\). Instead, it compares:

- the physical photon speed extracted from the actual finite Stage-6B spin-1 transverse spectrum;
- the harmonic frame speed
  \[
  c_g=\sqrt{rK_g}.
  \]

The diagonal plaquette ratio \(v/t_{\rm plaq}\) is left explicit because the shared virtual graph has not yet derived it.

For each tested \(v/t_{\rm plaq}\), the common-cone equation fixes one \(x=t/\Delta\):

| \(v/t_{\rm plaq}\) | \(x=t/\Delta\) | \(u/t_{\rm plaq}\) | photon gap power | common speed in \(\Delta\) units |
|---:|---:|---:|---:|---:|
| 0.2 | 0.143652630 | 0.928700804 | 1.009630 | 0.0117403 |
| 0.4 | 0.153016758 | 0.735698834 | 1.055388 | 0.0124749 |
| 0.6 | 0.167443752 | 0.530591435 | 1.135849 | 0.0135964 |

The maximum photon linear-plus-cubic fit residual is

\[
2.98\times10^{-4}.
\]

The photon and harmonic gravity speeds agree at the roots to better than

\[
4.2\times10^{-6}
\]

fractionally.

These roots use the actual finite photon eigenvalues, not the bare relation \(K_A=K_g\).

## 6. Finite compact gravity at the same charging coefficient

For a \(p\)-state clock variable, the physical angle is

\[
\theta=\frac{2\pi q}{p}.
\]

To represent the same physical charging coefficient \(U_g=r\Delta\), the finite shift coefficient must be

\[
U_{\rm clock}
=U_g\left(\frac{p}{2\pi}\right)^2.
\]

This factor is canonical normalization, not an independently adjusted gravity clock.

At all three common-cone roots, the exact compact gravity spectra for

\[
p=11,13,17
\]

have:

- gap powers between 1.0217 and 1.0333;
- maximum linear-plus-cubic fit residual below \(4.46\times10^{-4}\);
- exactly degenerate signed momenta;
- speeds that move monotonically toward the harmonic common-cone target as \(p\) increases.

For example, at \(v/t=0.4\):

| \(p\) | fitted power | finite \(c_g\) | error from harmonic target |
|---:|---:|---:|---:|
| 11 | 1.032659 | 0.0102807 | -17.59% |
| 13 | 1.027414 | 0.0105776 | -15.21% |
| 17 | 1.022926 | 0.0108649 | -12.91% |

The finite regulator is converging toward, but has not reached, the harmonic speed at these local dimensions.

### Coarse-clock negative control

At \(p=5\), the representative spectrum has fitted power

\[
0.250188,
\]

not a linear branch. Thus the common charging matrix does not make every finite clock representation acceptable. Sufficient local dimension remains a real regulator requirement.

## 7. What Stage 3C accomplishes

Stage 3C repairs the specific failure found in Stage 3B without introducing separate photon and gravity charging factors.

It establishes one candidate chain:

\[
\begin{aligned}
&8\times18=24\times6=144\ \text{hinge channels}\\
&\Downarrow\\
&\text{one connected collective charging mode}\\
&\Downarrow\\
&U_A=U_g=\Delta/144\\
&\Downarrow\\
&\text{common-cone roots in the actual photon spectrum}\\
&\Downarrow\\
&\text{linear compact gravity spectra at the same }r.
\end{aligned}
\]

This is progress toward a genuinely shared alternative-gravity mechanism rather than two independently normalized regulators.

## 8. What remains unfixed

The result is still a family, not a unique theory.

The roots depend on

\[
\frac{v}{t_{\rm plaq}},
\]

the diagonal plaquette/flippability coefficient. The current shared-move graph derived the off-diagonal ring exchange but did not derive this diagonal coefficient.

Consequently:

- \(r\) is now structurally fixed;
- the common cone fixes \(x=t/\Delta\) only after \(v/t\) is specified;
- \(v/t\) remains one independent dimensionless input;
- the full three-dimensional photon phase has not been checked at the new roots, whose \(u/t\) values exceed the original Stage-6B scan;
- finite-p gravity has not yet reached exact speed equality;
- matter and companion couplings have not been derived from the same matrix.

## 9. Originality and prior-art boundary

Parallel capacitances adding and graph-Laplacian locking are standard methods. They are not claimed as inventions.

The candidate project-specific contribution is the structural identification

\[
8\times18=24\times6=144
\]

inside the already committed matter–geometry cube, combined with the requirement that this one finite charging mode normalize both the photon and frame sectors generated by the same elementary swaps.

This identification may still have undiscovered prior art and is not asserted to be unique.

## 10. Next gate

Derive the diagonal plaquette coefficient from the same virtual move graph.

The next stage must compute both:

\[
-t_{\rm plaq}(W+W^\dagger)
\]

and

\[
+vD
\]

from one low-energy elimination, including all return paths and counterterms allowed by the microscopic junction symmetry.

No spectral observable may be used to choose \(v/t\). Once \(v/t\) is derived:

1. the common charging matrix fixes \(r=1/144\);
2. the common-cone equation fixes \(x=t/\Delta\);
3. the full three-dimensional photon model must be rerun at that point;
4. the compact gravity regulator must be rerun at the same point and increasing \(p\);
5. matter and companion sectors must use the same capacitance matrix;
6. the branch must be rejected if sector-specific clocks or penalties are again required.

## 11. Reproduction

```sh
python phase_junction_network/microscopic/check_collective_capacitance.py \
  --output phase_junction_network/microscopic/collective_capacitance_results.json
```

A reduced CI run is available with `--quick`.
