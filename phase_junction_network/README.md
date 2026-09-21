# Phase Junction Network

**Status:** exploratory candidate foundation. The repository now contains a finite pure-gauge photon phase, a finite linear dressed-frame regulator, a finite chiral charged-endpoint prototype, a reversible photon–companion architecture, one shared completed-move rule across photon/frame/matter/companion sectors, exact real- and virtual-pair Ward tests, and a complete antisymmetric neutral Fock-space calculation through four pairs on the finite spatial patch. The interacting domain-wall multiplet, matching frame-loop renormalization, infinite-volume QED, nonlinear gravity, continuum universality, observed matter, many-body deposits, and empirical predictions remain open.  
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

Finite quantum links, RK Hamiltonians, graph-Laplacian/seagull completions, relative-phase modes, domain-wall fermions, first-order frame/connection gravity, emergent helicity-2 spin models, finite Ward identities, finite vacuum-polarization calculations, and broad gauge–gravity unification all have prior art. They are not claimed as Phase Junction inventions.

The potentially original contribution is narrower:

> one finite matter–geometry junction algebra applies the same completed move and collective normalization to photon, frame, matter, companion, and pair-creation sectors, derives cross-sector coefficients, and survives interacting and nonlinear consistency tests with fewer inputs than outputs.

## Finite electromagnetic branch

The finite link Hamiltonian has the low-band form

\[
H_A=\frac{U}{2}\sum_\ell E_\ell^2
-K\sum_p(W_p+W_p^\dagger)
+v\sum_pD_p
+J_2\sum_p(W_p^2+W_p^{\dagger2})+\cdots.
\]

### Stage 6A — static Coulomb anchor

The exact spin-1 RK parent has exact Gauss symmetry, Coulomb winding response, a rank-two equal-time tensor, perimeter-favored Wilson shifts, and a fixed-charge response whose effective string slope decreases with size. The exact periodic `2^3` component has `146,327` states and `1,236,144` plaquette-transition edges.

### Stage 6B — detuned pure-gauge dynamics

A finite detuned region supports two transverse photon-like branches. Across spin-1 and spin-2 controls:

- gap powers: `0.978964`–`1.152567`;
- maximum linear-plus-cubic residual: `4.31e-4`;
- minimum electric residue: `0.992529`;
- scalar residue below `5.47e-31` in the photon branch;
- scalar gap at least `2.8496` times the photon gap;
- cubic tensor symbol `{0,k_hat^2,k_hat^2}`;
- anisotropy recovery approximately `L^-2.0111`.

The RK controls remain approximately `z=2`, so the static RK point is not relabeled as the relativistic phase. See [`microscopic/finite_em_dynamics.md`](microscopic/finite_em_dynamics.md).

## Finite linear gravity regulator

The committed continuum, exact-arithmetic, finite-state, and full-real-space checks establish under linear assumptions:

- three vector and one scalar frame constraints;
- 36 second-class connection constraints and four first-class frame constraints;
- two positive linearly dispersing configuration modes per nonzero momentum;
- no inserted transverse-traceless projector;
- a finite Weyl connection lock and exact dressed-frame algebra;
- linear tensor gaps with cubic cutoff corrections and vanishing lattice anisotropy.

This remains a linear regulator, not nonlinear quantum gravity.

## Shared microscopic move program

The candidate completed move is

\[
\boxed{
H_M=-x(M+M^\dagger)+x^2(M^\dagger M+MM^\dagger),
}
\]

with one unit completion coefficient in every sector.

| Stage | Main result | Boundary |
|---|---|---|
| 3A | One elementary swap gives `K_A ~ 20 t^4/Delta^3` and `K_g ~ t^2/Delta` | Toy mediator graph only |
| 3B | `U_A=U_g=Delta` maps to a strongly gapped photon point and is rejected | Requires collective charging |
| 3C | A `144 x 144` lock leaves one collective mode and gives `U_A=U_g=Delta/144` | Bundle not proven unique |
| 3D | Exact matching finds a dominant `O(t^2/Delta)` electric self-energy and rejects pure swap | No fitted deletion allowed |
| 3E | The unit completed move cancels the dangerous self-energy and yields a finite shared-cone root | Symmetry origin still open |
| 3F | The same rule passes finite photon, dressed-frame, charged-endpoint, and companion representatives | Separate finite representatives |
| 3G | One spatial Hamiltonian with a real charged pair passes nonzero-momentum Ward and photon-dressing gates | Real-pair response |
| 3H | Vacuum plus one virtual pair passes sector-changing Ward and transverse-response gates | One-pair truncation |
| 3I | Complete antisymmetric neutral Fock space through four pairs passes and converges | Finite `2 x 2` light-endpoint model |

### Stage 3E common finite normalization

The accepted finite candidate root is

\[
\frac{t}{\Delta}=0.1355417851,
\]

with

\[
K_A/\Delta=0.00560619757,
\quad
K_g/\Delta=0.01805149687,
\quad
U/\Delta=0.00678991375,
\]

and photon gap power `0.977345`.

See [`microscopic/seagull_completion_stage3e.md`](microscopic/seagull_completion_stage3e.md).

## Stage 3I — complete finite antisymmetric Fock space

The periodic `2 x 2` spin-1 patch now contains four positive and four negative fermion modes with exact Jordan–Wigner signs. All neutral sectors are included:

| Pair number | Physical states |
|---:|---:|
| 0 | 115 |
| 1 | 1,484 |
| 2 | 3,138 |
| 3 | 1,484 |
| 4 | 115 |
| **Total** | **6,336** |

The finite matter algebra has:

- canonical anticommutator residual `0`;
- two-creation exchange phase `-1`;
- pair-creation nilpotency residual `0`;
- signed translation commutator residual `0`.

The complete-Fock ground-state pair distribution is

```text
N=0: 0.9917488633
N=1: 0.00823075324
N=2: 2.03673e-5
N=3: 1.61925e-8
N=4: 5.28991e-12
```

The photon pole converges rapidly with pair-number cutoff:

| Maximum pairs | Photon gap | Gap shift from pure gauge | Residue |
|---:|---:|---:|---:|
| 1 | 0.01788196595 | `2.25098e-5` | 0.98765585 |
| 2 | 0.01788165199 | `2.21958e-5` | 0.98764958 |
| 3 | 0.01788165160 | `2.21954e-5` | 0.98764958 |
| 4 | 0.01788165160 | `2.21954e-5` | 0.98764958 |

The full result changes the one-pair photon correction by only `-1.3965%`.

Across the complete 6,336-state Hamiltonian:

- local continuity residual: `0`;
- nonzero-momentum operator Ward residual: `8.28e-18`;
- low spectral Ward residual: `2.31e-17`;
- longitudinal gauge-equivalence residual: `4.97e-16`;
- transverse polarization: `1.39868e-4 Delta`;
- transversality ratio: `2.31e-6`;
- no new matter-loop counterterm.

The finite photon pole is shifted by `0.124278%` relative to the pure-gauge pole. That is small but nonzero. A matching frame-loop calculation is required before claiming that the photon–gravity cone is radiatively stable.

The pair-gap scan gives approximately

\[
\langle N_{\rm pair}\rangle\propto M_{\rm pair}^{-2.124},
\qquad
\Pi_T\propto M_{\rm pair}^{-3.146}.
\]

See [`microscopic/fermionic_multipair_vacuum_stage3i.md`](microscopic/fermionic_multipair_vacuum_stage3i.md).

## Finite charged matter prototype

The wider issue-#4 construction supplies an anomaly-free same-chirality finite search result,

\[
q=(-11,-5,-1,-1,9,9),
\qquad
\sum q=\sum q^3=0,
\]

plus a finite domain-wall slab with a light Weyl cone, gapped doublers, a remote mirror wall, exponentially protected finite-width gaps, and one bare frame operator for every species.

The Stage-3I interacting patch currently contains the light `|q|=1` endpoint algebra, not the full fifth-dimensional wall and mirror multiplet.

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
| #3 | Shared microscopic coefficients | Stages 3E–3I pass finite all-sector, Ward, and complete finite-Fock gates; matching frame-loop renormalization open |
| #4 | Chiral matter defects | Free domain-wall prototype and interacting light endpoint pass; full wall/mirror Fock embedding open |
| #5 | Nonlinear gravity, self-coupling, many-body binding, volume/vacuum | Open |
| #6 | Finite electromagnetic/QED phase | Static, pure-gauge, real-pair, virtual-pair, and finite multipair tests pass; volume scaling and QED continuum open |
| #7 | Photon–companion integration | Finite architecture scope closed; common microscopic propagation open |
| #8 | Continuum consistency, mirror completion, Lorentz/common-cone recovery | Open |
| #9 | Parameter ledger, frozen predictions, empirical tests | Open |

## Immediate next gates

1. Couple the finite domain-wall light and mirror modes to the exact antisymmetric Fock construction.
2. Compute photon and frame self-energies from the same matter loop.
3. Test whether the `0.124%` finite photon-cone correction has the matching frame correction without a new coefficient.
4. Add spatial-volume scaling beyond the `2 x 2` torus.
5. Reject the candidate if matter loops require a photon-only or frame-only counterterm.

## Stop rules

Do not claim originality for established component techniques. Do not call Stage 3I an infinite-volume fermion determinant, charge-renormalization flow, or interacting-QED continuum. Do not introduce sector-specific completion coefficients, ignore generated higher operators, identify the finite matter spectrum with observed particles, interpret finite companion probabilities as astrophysical rates, call fixed volume a cosmological-constant solution, or call the linear frame regulator nonlinear quantum gravity. Derive one symmetry-complete interacting Hamiltonian across all sectors—or reject the combined branch.
