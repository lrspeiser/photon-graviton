# Phase Junction Network

**Status:** exploratory candidate foundation. The repository now contains a finite pure-gauge photon phase, a finite linear dressed-frame regulator, a finite chiral charged-endpoint prototype, a reversible photon–companion architecture, one shared completed-move rule across four finite sectors, an interacting real-pair Ward patch, and a first vacuum-plus-virtual-pair photon-dressing calculation. Multi-pair fermionic loops, interacting QED, nonlinear gravity, continuum universality, observed matter, many-body deposits, and empirical predictions remain open.  
**Created:** 2026-09-20  
**Prior-art boundary:** [`novelty_boundary.md`](novelty_boundary.md)

## Working ontology

| Object | Role |
|---|---|
| Photon | Transverse quantum of the compact U(1) link connection |
| Companion `chi=varphi` | Neutral matter–geometry relative-phase excitation |
| Deposit | Bound state of the same `chi` sector |
| Graviton | Transverse dressed-frame excitation |
| Frame connection | Constrained auxiliary comparison variable with a gapped relative sector |
| Matter | Finite charged or neutral endpoint defect |

The photon, companion, graviton, and connection-relative mode are distinct. Energy transferred among them must be carried by an explicit reversible interaction and counted once.

## Originality boundary

Finite quantum links, RK Hamiltonians, graph-Laplacian/seagull completions, relative-phase modes, domain-wall fermions, first-order frame/connection gravity, emergent helicity-2 spin models, finite Ward identities, and broad gauge–gravity unification all have prior art. They are not claimed as Phase Junction inventions.

The potentially original contribution is narrower: one finite matter–geometry junction Hilbert space and one symmetry-complete move/capacitance algebra must generate electromagnetic, frame, matter, and companion sectors while deriving cross-sector coefficients that would otherwise be independent.

## Finite electromagnetism

The finite integer-spin link Hamiltonian is

\[
H_A=\frac{U}{2}\sum_\ell E_\ell^2
-K\sum_p(W_p+W_p^\dagger)
+v\sum_pD_p
+J_2\sum_p(W_p^2+W_p^{\dagger2})+\cdots.
\]

### Stage 6A — static Coulomb anchor

At the RK point, the finite spin-1 model has exact Gauss symmetry, Coulomb winding response, a rank-two equal-time tensor, perimeter-favored Wilson shifts, and a fixed-charge response whose effective string slope decreases with size. The exact periodic `2^3` gauge component has `146,327` states and `1,236,144` plaquette-transition edges.

### Stage 6B — detuned pure-gauge dynamics

A finite detuned region supports two transverse photon-like branches. Across spin-1 and spin-2 controls:

- gap powers: `0.978964`–`1.152567`;
- maximum linear-plus-cubic residual: `4.31e-4`;
- minimum electric residue: `0.992529`;
- scalar residue below `5.47e-31` in the photon branch;
- scalar gap at least `2.8496` times the photon gap;
- cubic tensor symbol `{0,k_hat^2,k_hat^2}`;
- anisotropy recovery approximately `L^-2.0111`.

The RK controls remain approximately `z=2`, so the static RK point is not relabeled as the relativistic phase.

See [`microscopic/finite_em_dynamics.md`](microscopic/finite_em_dynamics.md).

## Finite linear gravity regulator

The committed continuum, exact-arithmetic, finite-state, and full-real-space checks establish under linear assumptions:

- three vector and one scalar frame constraints;
- 36 second-class connection constraints and four first-class frame constraints;
- two positive linearly dispersing configuration modes per nonzero momentum;
- no inserted transverse-traceless projector;
- a finite Weyl connection lock and exact dressed-frame algebra;
- linear tensor gaps with cubic cutoff corrections and vanishing lattice anisotropy.

This remains a linear regulator, not nonlinear quantum gravity.

## Shared electromagnetic/frame mechanism

### Stages 3A–3D

One elementary swap and one virtual penalty generate

\[
K_A=20t^4/\Delta^3+\cdots,
\qquad
K_g=t^2/\Delta+\cdots.
\]

A `144 x 144` collective lock matrix gives one uniform charging mode and

\[
U_A=U_g=\Delta/144.
\]

Exact low-band matching showed that a pure off-diagonal swap also generates a dominant `O(t^2/Delta)` electric self-energy and a double-flux operator. The pure-swap graph was rejected rather than patched with a sector-specific counterterm.

### Stage 3E — local completed move

Every elementary move is paired with its algebraic diagonal partner:

\[
\boxed{
H_M=-x(M+M^\dagger)+x^2(M^\dagger M+MM^\dagger)
}
\]

with unit coefficient. The completed photon/frame candidate has

\[
\frac{t}{\Delta}=0.1355417851,
\]

\[
K_A/\Delta=0.00560619757,
\quad
K_g/\Delta=0.01805149687,
\quad
U/\Delta=0.00678991375,
\]

and photon gap power `0.977345`.

See [`microscopic/seagull_completion_stage3e.md`](microscopic/seagull_completion_stage3e.md).

### Stage 3F — one rule in four finite sectors

The same unit completion was applied to photon, dressed-frame, charged-endpoint, and neutral-companion representatives.

Key results:

- exact positive-square identity for every completed move;
- dressed-frame lock commutator below `2.3e-16` and leakage below `1.2e-15`;
- complete frame harmonic inventory at roundoff;
- exact charged-endpoint Gauss law and continuity equation;
- two-site spectral Ward residual `4.65e-16`;
- pure-gauge diamagnetic/paramagnetic cancellation `5.55e-17`;
- charged hopping power `1.99533` and residual self-energy power `3.99129`;
- companion conversion power `1.99719` and residual self-energy power `3.99531`;
- companion reciprocity residual `2.78e-16`.

See [`microscopic/all_sector_completed_move_stage3f.md`](microscopic/all_sector_completed_move_stage3f.md).

### Stage 3G — finite spatial real-pair Ward patch

A periodic `2 x 2` spin-1 patch contains the Stage-3E photon operators and one dynamical `q=+1`, `q=-1` endpoint pair in one exact Gauss-law Hilbert space.

At momentum `(pi,0)`:

- local continuity residual: `0`;
- operator Ward residual: `8.28e-18`;
- low spectral Ward residual: `2.02e-15`;
- exact longitudinal unitary-equivalence residual: `6.94e-18`;
- transverse polarization curvature: `0.040161632 Delta`;
- transversality ratio: `3.46e-9`.

The photon pole moves from `0.0178594562 Delta` to `0.0354774349 Delta` while retaining residue `0.965396`. No new counterterm is introduced.

See [`microscopic/spatial_gauge_matter_stage3g.md`](microscopic/spatial_gauge_matter_stage3g.md).

### Stage 3H — finite virtual-pair vacuum

The same `2 x 2` photon patch now contains a no-matter vacuum block and a gauge-invariant one-pair block, coupled by the same unit completed move.

The combined dimension is `1599`. The interacting ground state remains vacuum dominated, with pair occupation

\[
0.00823338709.
\]

The photon pole shifts from

\[
0.0178594561978367\Delta
\]

to

\[
0.0178819659525319\Delta,
\]

giving

\[
\delta\omega_\gamma^2
=8.04530645062\times10^{-7}\Delta^2.
\]

The photon residue ratio is `0.998320209`.

Across the vacuum/pair blocks:

- local continuity residual: `0`;
- nonzero-momentum operator Ward residual: `8.28e-18`;
- spectral Ward residual: `1.59e-17`;
- longitudinal unitary-equivalence residual: `2.29e-16`;
- transverse virtual-pair response: `1.41074e-4 Delta`;
- no new pair- or photon-specific counterterm.

This is a finite one-pair-truncation vacuum effect, not an infinite-volume fermion determinant or interacting-QED continuum.

See [`microscopic/virtual_pair_vacuum_stage3h.md`](microscopic/virtual_pair_vacuum_stage3h.md).

## Finite charged matter prototype

The wider issue-#4 endpoint construction supplies a finite anomaly-free same-chirality search result,

\[
q=(-11,-5,-1,-1,9,9),
\qquad
\sum q=\sum q^3=0,
\]

plus a finite domain-wall slab with a light Weyl cone, gapped doublers, a remote mirror wall, exponentially protected finite-width gaps, and one bare frame operator for every species.

This is not the Standard Model. Mirror completion, observed charges and masses, generations, bound-state clocks, and interacting radiative stability remain open.

## Finite photon–companion bridge

A reversible finite source-to-receiver architecture contains photons, neutral `chi`, recoil, a computed bound `chi`, a receiver, and one shared frame. Its ledgers close, but its current propagation benchmark has

```text
maximum chi speed / maximum photon speed = 0.3571428571
```

so common-cone recovery for the companion remains open.

## Current work program

| Issue | Pillar | Status |
|---:|---|---|
| #2 | Finite dressed frame Hamiltonian | Linear finite-regulator scope closed |
| #3 | Shared microscopic coefficients | Stages 3E–3H pass finite all-sector, real-pair, and virtual-pair gates; many-pair stability open |
| #4 | Chiral matter defects | Finite free-regulator prototype closed; interacting domain-wall completion open |
| #5 | Nonlinear gravity, self-coupling, many-body binding, volume/vacuum | Open |
| #6 | Finite electromagnetic/QED phase | Static, pure-gauge, real-pair, and one-virtual-pair tests pass; multi-pair QED open |
| #7 | Photon–companion integration | Finite architecture scope closed; common microscopic propagation open |
| #8 | Continuum consistency, mirror completion, Lorentz/common-cone recovery | Open |
| #9 | Parameter ledger, frozen predictions, empirical tests | Open |

## Immediate next gates

1. Replace the one-pair truncation with the finite domain-wall endpoint spectrum and antisymmetric multi-pair sectors.
2. Test scaling of transverse polarization and photon-pole renormalization with volume, matter gap, and pair number.
3. Verify Ward–Takahashi identities in the enlarged Fock space.
4. Test whether the shared photon/frame cone remains stable after matter dressing.
5. Reject the candidate if a matter-loop-specific completion coefficient or counterterm is required.

## Stop rules

Do not claim originality for established component techniques. Do not call Stage 3H an interacting-QED continuum or its one-pair result an infinite-volume fermion determinant. Do not introduce sector-specific completion coefficients, ignore generated higher operators, identify the finite matter spectrum with observed particles, interpret finite companion probabilities as astrophysical rates, call fixed volume a cosmological-constant solution, or call the linear frame regulator nonlinear quantum gravity. Derive one symmetry-complete interacting Hamiltonian across all sectors—or reject the combined branch.
