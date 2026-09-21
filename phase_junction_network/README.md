# Phase Junction Network

**Status:** exploratory candidate foundation. The repository contains a finite pure-gauge photon phase, a finite linear dressed-frame regulator, a finite chiral charged-endpoint prototype, a reversible photon–companion architecture, and a structurally derived collective charging mode. Stage 3D rejected the pure-swap shared graph; Stage 3E now identifies a symmetry-complete local seagull candidate that cancels the dangerous self-energy and restores one finite relativistic shared-cone root. Its fundamental symmetry, full finite embedding, interactions, nonlinear gravity, and empirical predictions remain open.  
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

Finite quantum links, RK Hamiltonians, graph-Laplacian/seagull completions, relative-phase modes, domain-wall fermions, first-order frame/connection gravity, emergent helicity-2 spin models, and broad gauge–gravity unification all have prior art. They are credited in [`novelty_boundary.md`](novelty_boundary.md) and are not claimed as Phase Junction inventions.

The potentially original contribution is narrower: one finite matter–geometry junction Hilbert space and one symmetry-complete local move/capacitance structure must generate the electromagnetic, frame, matter, and companion sectors while deriving cross-sector coefficients that would otherwise be independent.

## Finite electromagnetism

The finite integer-spin link Hamiltonian is

\[
H_A(u,t,v)=\frac{u}{2}\sum_\ell E_\ell^2
-t\sum_p(W_p+W_p^\dagger)+v\sum_pD_p.
\]

### Stage 6A — static Coulomb anchor

At `u/t=0, v/t=1`, the finite spin-1 model has an exact positive RK parent, exact Gauss symmetry, Coulomb winding response, a rank-two equal-time tensor, perimeter-favored Wilson shifts, and a fixed-charge response whose effective string slope decreases with size. The exact periodic `2^3` gauge component has `146,327` states and `1,236,144` plaquette-transition edges.

See [`microscopic/finite_em_coulomb_phase.md`](microscopic/finite_em_coulomb_phase.md).

### Stage 6B — detuned pure-gauge dynamics

A scan over

\[
u/t\in\{0,0.1,0.2\},
\qquad
v/t\in\{0.2,0.4,0.6\}
\]

finds a finite detuned region with two transverse photon-like branches:

- gap powers `0.978964`–`1.152567`;
- maximum linear-plus-cubic residual `4.31e-4`;
- minimum electric residue `0.992529`;
- scalar residue below `5.47e-31` in the photon branch;
- scalar gap at least `2.8496` times the photon gap;
- full cubic symbol `{0,k_hat^2,k_hat^2}`;
- anisotropy recovery approximately `L^-2.0111`.

The RK controls remain approximately `z=2`, so the static RK point is not being relabeled as the relativistic phase.

See [`microscopic/finite_em_dynamics.md`](microscopic/finite_em_dynamics.md).

Issue #6 remains open for dynamical matter, Ward identities, vacuum polarization, charge renormalization, a regulator-universal QED limit, and precision Lorentz tests.

## Finite linear gravity regulator

The committed continuum, exact-arithmetic, finite-state, and full-real-space checks establish under linear assumptions:

- three vector and one scalar frame constraints;
- 36 second-class connection constraints and four first-class frame constraints;
- two positive linearly dispersing configuration modes per nonzero momentum;
- no inserted transverse-traceless projector;
- a finite Weyl connection lock and exact dressed-frame algebra;
- two finite logical Weyl pairs carrying positive compact clock/Villain Hamiltonians;
- linear tensor gaps with cubic cutoff corrections and vanishing lattice anisotropy.

See [`microscopic/finite_dressed_frame_hamiltonian.md`](microscopic/finite_dressed_frame_hamiltonian.md).

This is a linear finite regulator, not nonlinear quantum gravity. Its old independently chosen self-dual clock normalization is not treated as fundamental.

## Shared electromagnetic/frame mechanism

### Stage 3A — one move, one penalty

A toy finite mediator uses one elementary swap amplitude `t` and one virtual penalty `Delta`. Exact path counting and diagonalization give

\[
K_A=20\frac{t^4}{\Delta^3}+O(t^6/\Delta^5),
\qquad
K_g=\frac{t^2}{\Delta}+O(t^4/\Delta^3).
\]

See [`microscopic/shared_junction_move_stage3a.md`](microscopic/shared_junction_move_stage3a.md).

### Stage 3B — naïve normalization rejected

The assumption

\[
U_A=U_g=\Delta
\]

maps into the real photon regulator as

\[
u/t_{\rm plaq}=14.5969,
\]

where the exact finite branch is strongly gapped rather than linearly dispersing.

See [`microscopic/shared_move_embedding_stage3b.md`](microscopic/shared_move_embedding_stage3b.md).

### Stage 3C — collective capacitance

The elementary cube has two independent counts of the same proposed microscopic bundle:

\[
8\times18=24\times6=144.
\]

An explicit `144 x 144` connected lock matrix has rank 143 and one uniform charging mode. With identical hinge capacitances,

\[
U_A=U_g=\Delta/144.
\]

See [`microscopic/collective_capacitance_stage3c.md`](microscopic/collective_capacitance_stage3c.md).

### Stage 3D — pure swap rejected

Exact three-state matching of the spin-1 plaquette showed that the same swaps generate

\[
K=20x^4\Delta+\cdots,
\qquad
\delta U=2x^2\Delta+\cdots,
\qquad
v=-\frac83x^4\Delta+\cdots.
\]

The second-order electric self-energy dominates the fourth-order loop term and eliminates the counterterm-free common cone. A direct `W^2+W^dagger^2` operator is also generated. The pure-swap graph is therefore rejected as the complete microscopic theory.

See [`microscopic/plaquette_diagonal_matching_stage3d.md`](microscopic/plaquette_diagonal_matching_stage3d.md).

### Stage 3E — local seagull completion

Pair every finite-link move with its algebraic diagonal partner:

\[
H_A/\Delta
=
\frac12\sum G^2
-x\sum(U+U^\dagger)
+x^2\sum(U^\dagger U+UU^\dagger).
\]

Use the same completion for the frame mediator. The unit coefficient is fixed by cancellation of each isolated second-order move self-energy; varying it above or below one leaves an uncancelled correction with the corresponding sign.

Exact matching gives

\[
K_A\propto x^{3.99087},
\qquad
\delta U\propto x^{3.99530},
\qquad
v\propto x^{3.99587},
\qquad
K_g\propto x^{1.99906}.
\]

The dangerous `O(x^2)` electric correction is gone while the photon ring and frame superexchange remain.

At the shared `r=1/144`, one relativistic common-cone root is

\[
\boxed{t/\Delta=0.1355417851.}
\]

There,

\[
K_A/\Delta=0.00560619757,
\qquad
K_g/\Delta=0.01805149687,
\]

\[
U/\Delta=0.00678991375,
\qquad
v/K_A=-0.129245131,
\]

and the finite photon gap power is

\[
0.977345
\]

with fit residual `2.91e-4`. The photon and harmonic frame speeds agree at `0.01107104814 Delta` within numerical root tolerance.

The complete `2^3` gauge-sector anchor has two lowest equal gaps, exactly zero longitudinal electric norm, and scalar contamination `2.05e-22` in the low multiplet. Compact frame regulators with `p>=7` remain linear at the same physical coefficients and move toward the harmonic speed as `p` grows; `p=5` is a failure control.

See [`microscopic/seagull_completion_stage3e.md`](microscopic/seagull_completion_stage3e.md).

**Claim boundary:** Stage 3E advances a viable candidate, not a completed fundamental theory. A symmetry that requires the completion, its matter/companion analogues, radiative stability, nonlinear gravity, and the full interacting continuum remain open.

## Finite charged matter prototype

Issue #4 constructs finite odd-strand endpoint defects with exact flux-bundle Gauss covariance. Within the declared search class, the first primitive non-vectorlike anomaly-free spectrum is

\[
q=(-11,-5,-1,-1,9,9),
\qquad
\sum q=\sum q^3=0.
\]

A finite domain-wall slab supplies a light Weyl cone, gapped doublers, a remote opposite-chirality wall, exponentially protected finite-width gaps, and one bare frame operator for every species.

This is not the Standard Model. Mirror completion, observed charges and masses, generations, bound-state clocks, and interacting radiative stability remain open.

See [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md).

## Finite photon–companion bridge

Issue #7 supplies a reversible finite source-to-receiver architecture containing photons, neutral `chi`, recoil, a computed bound `chi`, a receiver, and one shared frame. The finite ledgers close, and no independent motion or lensing multiplier is introduced.

The benchmark also finds

```text
maximum chi speed / maximum photon speed = 0.3571428571
```

so common-cone recovery is not yet established.

See [`companion_integration.md`](companion_integration.md).

## Current work program

| Issue | Pillar | Status |
|---:|---|---|
| #2 | Finite dressed frame Hamiltonian | Linear finite-regulator scope closed |
| #3 | Shared microscopic coefficients | Stage 3D rejects pure swap; Stage 3E provides a viable seagull-completed candidate; full cross-sector embedding and symmetry protection open |
| #4 | Chiral matter defects | Finite free-regulator prototype closed |
| #5 | Nonlinear gravity, self-coupling, many-body binding, volume/vacuum | Open |
| #6 | Finite electromagnetic/QED phase | Static and pure-gauge dynamics pass; interacting matter/QED open |
| #7 | Photon–companion integration | Finite architecture scope closed |
| #8 | Continuum consistency, mirror completion, Lorentz/common-cone recovery | Open |
| #9 | Parameter ledger, frozen predictions, empirical tests | Open |

## Immediate next gates

1. Construct the seagull-completed moves directly in the dressed-frame, charged-endpoint, and companion finite Hilbert spaces.
2. Verify exact Gauss and frame constraints and derive the full generated operator inventory in every sector.
3. Reject the candidate if photon, gravity, matter, or companion requires a different completion coefficient.
4. Compute the first finite interacting Ward identity and test whether the cancellation is radiatively stable.
5. Continue the nonlinear gravity and many-body companion program only after the shared finite Hamiltonian survives these gates.

## Stop rules

Do not claim originality for established component techniques. Do not call Stage 6B interacting QED, Stage 3A a derivation of `alpha` or `G`, Stage 3C a unique microscopic bundle, Stage 3D's counterterm root a solution, or Stage 3E a fundamental symmetry. Do not introduce sector-specific seagulls, ignore generated operators, identify the finite matter set with observed particles, interpret finite companion probabilities as astrophysical rates, call fixed volume a cosmological-constant solution, or call the linear frame regulator nonlinear quantum gravity. Derive one symmetry-complete finite Hamiltonian across all sectors—or reject the combined branch.
