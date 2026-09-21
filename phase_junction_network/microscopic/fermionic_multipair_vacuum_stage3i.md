# Stage 3I: Antisymmetric Multipair Vacuum

**Status:** **PASS at finite antisymmetric multipair scope.** The Stage-3H zero/one-pair truncation has been replaced by the complete neutral spinless-fermion Fock space on the periodic `2 x 2` gauge patch. Exact fermionic signs, Pauli blocking, local continuity, nonzero-momentum Ward identities, longitudinal gauge invariance, and transverse virtual polarization survive through all four allowed pairs.  
**Date:** 2026-09-20  
**Issues:** #3, #4, #6, and #8  
**Executable:** [`check_fermionic_multipair_vacuum.py`](check_fermionic_multipair_vacuum.py)  
**Frozen output:** [`fermionic_multipair_vacuum_results.json`](fermionic_multipair_vacuum_results.json)

## 1. Purpose

Stage 3H allowed only the vacuum and one neutral endpoint pair. That was sufficient to demonstrate sector-changing Ward identities and a finite virtual-pair photon self-energy, but it did not establish fermionic many-pair consistency or convergence of the pair expansion.

Stage 3I adds the complete neutral Fock space of four positive and four negative spinless endpoint modes on the same four-site torus. It asks:

1. whether exact fermionic antisymmetry survives the gauge constraint;
2. whether pair creation remains Pauli blocked and nilpotent;
3. whether Ward identities survive across all pair-number sectors;
4. whether the one-pair photon shift converges when two-, three-, and four-pair states are admitted;
5. how the response scales with the pair gap and common pair-creation amplitude;
6. how large a finite photon-cone shift remains before a matching frame-loop calculation.

## 2. Finite Hilbert space

The spatial patch remains:

- four sites;
- eight positive-oriented spin-1 links;
- four plaquettes;
- four positive-fermion modes;
- four negative-fermion modes.

Only neutral sectors with equal positive and negative particle number are retained. Every basis state also satisfies the exact link Gauss law.

The sector dimensions are:

| Pair number | Physical states |
|---:|---:|
| 0 | 115 |
| 1 | 1,484 |
| 2 | 3,138 |
| 3 | 1,484 |
| 4 | 115 |
| **Total** | **6,336** |

The particle–hole-symmetric counting is an output of the exact finite construction.

## 3. Fermionic algebra

The matter modes are ordered as four positive modes followed by four negative modes. Every hopping and pair-creation matrix element includes the exact Jordan–Wigner parity of the occupied modes below it.

On the complete unconstrained eight-mode matter Fock space,

\[
\{c_i,c_j^\dagger\}=\delta_{ij},
\qquad
\{c_i,c_j\}=0
\]

with maximum matrix residual

\[
0.
\]

Exchanging the order of two creation operators gives phase

\[
\boxed{-1}.
\]

Every oriented pair-creation move is Pauli nilpotent,

\[
P_\ell^2=0,
\]

with residual

\[
0.
\]

Thus the multipair extension is an antisymmetric fermionic Fock construction rather than a distinguishable-pair enlargement.

## 4. Signed translation symmetry

Translating a fermionic configuration requires both moving its occupied sites and including the parity of the induced permutation of the ordered creation operators.

With that sign included, the full Hamiltonian commutes with both lattice translations. The maximum translation residual is

\[
0.
\]

The transverse electric operator has characters

\[
T_x=-1,
\qquad
T_y=+1,
\]

and therefore remains the finite \((\pi,0)\) transverse probe.

## 5. Multipair ground state

At the shared Stage-3E/3F coefficients,

\[
M_{\rm pair}=\Delta,
\qquad
K_{\rm pair}=0.0338029748033\,\Delta,
\]

the complete Fock-space ground-state probabilities are:

| Pair number | Probability |
|---:|---:|
| 0 | 0.9917488633 |
| 1 | 0.00823075324 |
| 2 | \(2.03673\times10^{-5}\) |
| 3 | \(1.61925\times10^{-8}\) |
| 4 | \(5.28991\times10^{-12}\) |

The mean pair number is

\[
\langle N_{\rm pair}\rangle
=0.00827153639.
\]

The ground state remains strongly vacuum dominated, and every higher pair sector is hierarchically suppressed.

## 6. Convergence of the pair truncation

The dominant transverse photon pole converges rapidly as the maximum allowed pair number is increased:

| Maximum pairs | Hilbert dimension | Photon gap | Gap shift from pure gauge | Photon residue |
|---:|---:|---:|---:|---:|
| 1 | 1,599 | 0.01788196595 | \(2.25098\times10^{-5}\) | 0.98765585 |
| 2 | 4,737 | 0.01788165199 | \(2.21958\times10^{-5}\) | 0.98764958 |
| 3 | 6,221 | 0.01788165160 | \(2.21954\times10^{-5}\) | 0.98764958 |
| 4 | 6,336 | 0.01788165160 | \(2.21954\times10^{-5}\) | 0.98764958 |

The complete Fock result differs from the one-pair photon shift by

\[
-3.14350\times10^{-7}\,\Delta,
\]

or

\[
\boxed{-1.3965\%}
\]

of the one-pair correction.

The one-pair transverse static response is

\[
1.41073842\times10^{-4}\,\Delta,
\]

while the complete-Fock result is

\[
1.39868471\times10^{-4}\,\Delta.
\]

That is a change of less than one percent. The one-pair truncation was therefore quantitatively close on this finite patch, but Stage 3I no longer relies on that assumption.

## 7. Exact multipair Ward identities

The local charge density is

\[
\rho_r=n_{+,r}-n_{-,r}.
\]

The current contains:

- positive-endpoint hopping;
- negative-endpoint hopping with the opposite charge sign;
- pair creation and annihilation.

The complete 6,336-state Hamiltonian satisfies

\[
\boxed{
i[H,\rho_r]+\nabla\cdot j_r=0
}
\]

with residual

\[
0.
\]

At

\[
k=(\pi,0),
\]

the operator Ward residual is

\[
8.28\times10^{-18},
\]

and the maximum computed low-energy spectral residual is

\[
2.31\times10^{-17}.
\]

The exact Ward structure therefore survives the entire finite antisymmetric Fock space.

## 8. Longitudinal and transverse response

A longitudinal external mode is exactly related to the zero-field Hamiltonian by the finite gauge transformation of every occupied fermion mode. The matrix-equivalence residual is

\[
4.97\times10^{-16}.
\]

The full-Fock static curvatures are:

\[
\Pi_L=-3.23\times10^{-10},
\]

\[
\boxed{
\Pi_T=1.39868471\times10^{-4}\,\Delta,
}
\]

and

\[
\Pi_{LT}=-5.07\times10^{-11}.
\]

The transversality ratio is

\[
2.31\times10^{-6}.
\]

When pair creation is disabled, the transverse response falls to

\[
-4.86\times10^{-11},
\]

which is the numerical floor.

## 9. Pair-gap dependence

The complete Fock calculation was repeated for pair gaps between

\[
0.75\Delta
\quad\text{and}\quad
1.5\Delta.
\]

The mean pair number scales approximately as

\[
\langle N_{\rm pair}\rangle
\propto M_{\rm pair}^{-2.124}.
\]

The transverse polarization scales approximately as

\[
\boxed{
\Pi_T\propto M_{\rm pair}^{-3.146}.
}
\]

Representative values are:

| Pair gap | Mean pairs | \(\Pi_T/\Delta\) |
|---:|---:|---:|
| 0.75 | 0.0154831 | \(3.54673\times10^{-4}\) |
| 1.00 | 0.00827154 | \(1.39868\times10^{-4}\) |
| 1.50 | 0.00354552 | \(3.99488\times10^{-5}\) |

The first physical doubler gap from the earlier finite \(q=-1\) domain-wall calculation,

\[
1.167449903\Delta,
\]

was also used as a mass control. It gives mean pair number

\[
0.00596453
\]

and photon residue

\[
0.98861417.
\]

This is only a gap control; the full fifth-dimensional domain-wall multiplet is not yet present in the interacting patch.

## 10. Pair-coupling dependence

Scaling the common pair-creation amplitude gives:

| Coupling scale | Mean pairs | Photon gap shift | Photon residue |
|---:|---:|---:|---:|
| 0.25 | 0.00051955 | \(1.47717\times10^{-6}\) | 0.98921320 |
| 0.50 | 0.00207614 | \(5.83712\times10^{-6}\) | 0.98889991 |
| 0.75 | 0.00466359 | \(1.28641\times10^{-5}\) | 0.98837835 |
| 1.00 | 0.00827154 | \(2.21954\times10^{-5}\) | 0.98764958 |

The virtual occupation and pole shift increase monotonically without introducing a new counterterm.

## 11. Finite common-cone diagnostic

Relative to the pure-gauge pole, the complete multipair photon gap changes by

\[
\boxed{
\frac{\Delta_\gamma^{(\rm Fock)}}{\Delta_\gamma^{(0)}}-1
=0.00124278.
}
\]

Thus the finite photon cone receives a small but nonzero matter-loop correction of approximately \(0.124\%\) at this momentum.

This does **not** establish radiative stability of the shared photon–frame cone. That claim requires the frame sector to receive the matter loop through the same microscopic stress-energy coupling and to produce the corresponding renormalization without an independently chosen coefficient.

## 12. What Stage 3I establishes

At finite `2 x 2` scope, Stage 3I establishes:

1. the complete neutral antisymmetric endpoint Fock space through four pairs;
2. exact canonical anticommutation, exchange phase, and Pauli blocking;
3. exact signed translation symmetry;
4. exact local and nonzero-momentum Ward identities;
5. longitudinal pure-gauge invariance;
6. nonzero transverse multipair polarization;
7. quantitative convergence of the one-pair truncation;
8. controlled mass- and coupling-dependence with the same unit completion coefficient;
9. no matter-loop-specific counterterm.

## 13. Claim boundary

Stage 3I does not establish:

- the complete fifth-dimensional domain-wall light and mirror multiplet in the interacting patch;
- mirror-wall removal or symmetric mass generation;
- volume scaling beyond the `2 x 2` torus;
- an infinite-volume fermion determinant or charge-renormalization flow;
- the matching matter-loop correction to the frame cone;
- a `3+1`-dimensional interacting-QED continuum;
- nonlinear gravity or observed matter.

## 14. Next gate

Couple the finite domain-wall light and mirror modes to this exact antisymmetric Fock construction. Then compute the photon and frame self-energies from the same matter loop and ask whether

\[
\delta c_\gamma=\delta c_g
\]

follows without a new counterterm or sector-specific completion coefficient.

The shared candidate must be rejected if matter dressing preserves the electromagnetic Ward identity but independently shifts the photon and frame cones.

## 15. Reproduction

```sh
python phase_junction_network/microscopic/check_fermionic_multipair_vacuum.py \
  --output phase_junction_network/microscopic/fermionic_multipair_vacuum_results.json
```

A reduced two-pair regression run is available with `--quick`.
