# Stage 3H: Virtual-Pair Photon Dressing in a Finite Vacuum

**Status:** **PASS at finite one-pair-truncation scope.** A no-matter vacuum block and a gauge-invariant neutral-pair block coexist with the Stage-3E photon Hamiltonian. The same unit completed-move rule generates pair creation and annihilation, exact Ward identities survive across the sector-changing moves, and virtual-pair admixture produces a small transverse photon self-energy without a new counterterm.  
**Date:** 2026-09-20  
**Issues:** #3 and #6  
**Executable:** [`check_virtual_pair_vacuum.py`](check_virtual_pair_vacuum.py)  
**Frozen output:** [`virtual_pair_vacuum_results.json`](virtual_pair_vacuum_results.json)

## 1. Why this is different from Stage 3G

Stage 3G measured polarization in a sector containing a permanently occupied positive/negative endpoint pair. Stage 3H adds a genuine vacuum block and allows the Hamiltonian to create and annihilate one neutral pair virtually.

The finite model still truncates the matter Fock space to zero or one pair. It is therefore not a fermion determinant or an infinite-volume vacuum-polarization calculation.

## 2. Hilbert space

The periodic `2 x 2` spin-1 patch contains:

- vacuum gauge sector: `115` states;
- one neutral-pair sector: `1484` states;
- combined physical dimension: `1599`.

Every basis state satisfies the exact lattice Gauss law.

## 3. Pair creation rule

For an oriented link from tail to head, pair creation produces:

- a `q=+1` endpoint at the head;
- a `q=-1` endpoint at the tail;
- one unit decrease of the oriented link flux.

That mapping preserves Gauss law exactly.

The sector-changing move uses the same completion rule as the photon, frame, endpoint hopping, and companion sectors:

\[
H_P=-g(P+P^\dagger)+g^2(P^\dagger P+PP^\dagger),
\]

with

\[
g=0.0338029748033\,\Delta
\]

and completion coefficient one.

The one-pair block is assigned the shared virtual defect gap

\[
M_{\rm pair}=\Delta.
\]

This gap is a declared finite-stage input. It is not yet the observed mass of a relativistic fermion.

## 4. Vacuum composition

The interacting ground state remains vacuum dominated:

\[
\boxed{
\langle N_{\rm pair}\rangle=0.00823338709.
}
\]

Thus approximately `0.823%` of the finite ground-state norm lies in the one-pair block.

## 5. Photon pole dressing

The pure-gauge transverse photon pole is

\[
\Delta_\gamma^{(0)}=0.0178594561978367\,\Delta,
\]

with residue

\[
Z_T^{(0)}=0.9893176943.
\]

After virtual-pair mixing:

\[
\Delta_\gamma^{(\rm vp)}=0.0178819659525319\,\Delta,
\]

and

\[
Z_T^{(\rm vp)}=0.9876558472.
\]

Therefore

\[
\boxed{
\delta\Delta_\gamma
=2.25097546952\times10^{-5}\,\Delta
}
\]

and

\[
\boxed{
\delta\omega_\gamma^2
=8.04530645062\times10^{-7}\,\Delta^2.
}
\]

The dressed-to-bare residue ratio is

\[
0.9983202088.
\]

The finite photon remains sharply identifiable.

## 6. Ward identities across the vacuum/pair blocks

The full Hamiltonian, including pair creation and annihilation, satisfies exact local continuity:

\[
i[H,\rho_r]+\nabla\cdot j_r=0
\]

with maximum matrix residual

\[
0.
\]

At

\[
k=(\pi,0),
\]

the nonzero-momentum operator Ward residual is

\[
8.28\times10^{-18},
\]

and the maximum low-energy spectral residual is

\[
1.59\times10^{-17}.
\]

These identities include matrix elements that connect the vacuum and one-pair sectors.

## 7. Static transversality

A longitudinal external mode is related to the zero-field Hamiltonian by an exact finite gauge transformation. The unitary-equivalence residual is

\[
2.29\times10^{-16}.
\]

The longitudinal curvature is

\[
-3.47\times10^{-12},
\]

while the transverse virtual-pair response is

\[
\boxed{
\Pi_T=1.41073867865\times10^{-4}\,\Delta.
}
\]

The transversality ratio is

\[
2.46\times10^{-8}.
\]

With pair creation disabled, the transverse response falls to the numerical floor,

\[
-2.78\times10^{-11}.
\]

## 8. Pair-coupling scan

| Pair-coupling scale | Pair occupation | Photon gap shift | Photon residue |
|---:|---:|---:|---:|
| 0.25 | 0.000519 | \(1.48\times10^{-6}\) | 0.989213 |
| 0.50 | 0.002074 | \(5.86\times10^{-6}\) | 0.988900 |
| 0.75 | 0.004651 | \(1.30\times10^{-5}\) | 0.988380 |
| 1.00 | 0.008233 | \(2.25\times10^{-5}\) | 0.987656 |

Both virtual-pair occupation and photon dressing grow monotonically with the common coupling.

## 9. What Stage 3H establishes

At finite one-pair scope:

1. one Hamiltonian contains a gauge vacuum and gauge-invariant pair sector;
2. pair creation and annihilation use the same unit completion coefficient;
3. local and nonzero-momentum Ward identities survive sector-changing moves;
4. longitudinal backgrounds remain pure gauge;
5. virtual-pair admixture produces a nonzero transverse photon self-energy;
6. no pair-specific or photon-specific counterterm is added.

## 10. Claim boundary

Stage 3H does not establish:

- relativistic fermionic statistics in a many-pair Fock space;
- a fermion determinant;
- infinite-volume vacuum polarization or charge renormalization;
- a 3+1-dimensional interacting-QED fixed point;
- radiative stability beyond the one-pair truncation;
- mirror-wall completion;
- nonlinear gravity.

The declared pair gap remains an input at this stage.

## 11. Next gate

Replace the one-pair truncation with the finite domain-wall endpoint spectrum and multiple pair sectors. Then test:

1. fermionic antisymmetry and pair statistics;
2. scaling of the transverse polarization with volume and matter gap;
3. charge renormalization from the finite photon pole and static response;
4. stability of the shared photon/frame cone;
5. the Ward–Takahashi identity in the enlarged Fock space;
6. whether the unit completion remains sufficient without a new matter-loop counterterm.

## 12. Reproduction

```sh
python phase_junction_network/microscopic/check_virtual_pair_vacuum.py \
  --output phase_junction_network/microscopic/virtual_pair_vacuum_results.json
```

A reduced CI run is available with `--quick`.
