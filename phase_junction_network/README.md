# Phase Junction Network

**Status:** exploratory candidate foundation. Issue #2 supplies a finite linear dressed-frame regulator, issue #4 supplies a finite chiral charged-endpoint prototype, issue #7 supplies a reversible photon–companion integration architecture, and issue #6 now has both a finite static Coulomb anchor (**Stage 6A**) and a detuned pure-gauge two-photon dynamical region (**Stage 6B**). Interacting QED, observed matter, shared microscopic coefficients, nonlinear gravity, common-cone recovery, many-body deposits, and empirical predictions remain open.  
**Created:** 2026-09-20

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

## Finite electromagnetism

The finite integer-spin link Hamiltonian is

\[
H_A(u,t,v)=\frac{u}{2}\sum_\ell E_\ell^2
-t\sum_p(W_p+W_p^\dagger)+v\sum_pD_p,
\qquad
D_p=\sqrt{W_p^\dagger W_p}+\sqrt{W_pW_p^\dagger}.
\]

### Stage 6A — static Coulomb anchor

At `u/t=0, v/t=1`, every local term is a positive weighted graph Laplacian and commutes exactly with Gauss law. The committed finite checks establish:

- an exactly enumerated periodic `2^3` spin-1 gauge component with `146,327` states and `1,236,144` undirected plaquette transitions;
- Coulomb winding scaling over `L=4,6,8`;
- a rank-two equal-time transverse tensor with no longitudinal component;
- perimeter-favored finite Wilson shifts;
- a fixed-charge response consistent with the periodic lattice Green function and a decreasing effective string slope;
- persistence of the static diagnostics in a spin-2 control.

See [`microscopic/finite_em_coulomb_phase.md`](microscopic/finite_em_coulomb_phase.md).

### Stage 6B — detuned pure-gauge dynamics

The same finite link operators were scanned over

\[
u/t\in\{0,0.1,0.2\},
\qquad
v/t\in\{0.2,0.4,0.6\},
\]

for spin 1 and spin 2. Exact principal-axis Gauss-reduced transverse blocks give two identical photon branches. All 18 detuned spin/coupling fits pass:

- finite-size gap powers: `0.978964`–`1.152567`;
- maximum linear-plus-cubic fit residual: `4.31e-4`;
- minimum electric-field residue in the first photon branch: `0.992529`;
- maximum scalar residue in that branch: `5.47e-31`;
- minimum scalar-gap/photon-gap ratio: `2.8496`;
- exact `+k/-k` degeneracy at reported precision;
- full cubic symbol `{0,k_hat^2,k_hat^2}` and anisotropy recovery `L^-2.0111`.

The RK controls remain nonrelativistic, with fitted powers `1.96037` and `1.82408`. Thus the static RK point is not being relabeled as the photon phase; the approximately linear branch appears in a finite detuned region.

An exact `2^3` full-cube anchor retains continuous overlap with the RK ground state throughout the tested detuning and has an exactly null longitudinal electric operator.

See [`microscopic/finite_em_dynamics.md`](microscopic/finite_em_dynamics.md), [`microscopic/check_finite_em_dynamics.py`](microscopic/check_finite_em_dynamics.py), and [`microscopic/finite_em_dynamics_results.json`](microscopic/finite_em_dynamics_results.json).

**Claim boundary:** Stage 6B passes the pure-gauge dynamical gate. Issue #6 remains open for dynamical charged matter, finite Ward identities, vacuum polarization, charge renormalization, regulator-universal QED scaling, and precision Lorentz tests.

## Finite linear gravity regulator

The continuum, exact-arithmetic, finite-state, full-real-space, and finite-size checks establish under their linear assumptions:

- three vector and one scalar frame constraints;
- 36 second-class connection constraints and four first-class frame constraints;
- four physical phase-space dimensions, or two positive linearly dispersing modes, per nonzero momentum;
- no inserted transverse-traceless projector;
- a finite odd-prime Weyl connection lock and exact dressed-frame algebra;
- two finite logical Weyl pairs at `p=5,7,11` carrying positive compact clock/Villain Hamiltonians;
- linearly closing tensor gaps with cubic cutoff corrections;
- connection-lock leakage at roundoff and rapidly decreasing compact-boundary occupation.

See [`microscopic/finite_dressed_frame_hamiltonian.md`](microscopic/finite_dressed_frame_hamiltonian.md).

This closes the **linear finite constrained-regulator** gate, not nonlinear quantum gravity. The gravity scale and clock ratio remain inputs until issue #3 derives shared coefficients. Nonlinear closure and the volume/vacuum sector remain issue #5; universality remains issue #8.

## Finite charged matter prototype

Issue #4 constructs finite odd-strand endpoint defects with exact flux-bundle Gauss covariance. Within the declared finite search class, the first primitive non-vectorlike anomaly-free spectrum is

\[
q=(-11,-5,-1,-1,9,9),
\qquad
\sum q=\sum q^3=0.
\]

A finite Wilson/domain-wall slab supplies one light Weyl corner on each wall, gaps the other seven physical corners, and produces exponentially protected finite-width gaps through a shared charge-to-localization rule. Every species uses the same bare frame derivative.

This is a regulator prototype, not the Standard Model. The remote mirror wall, observed charges and masses, generations, bound-state clocks, and interacting radiative stability remain open. See [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md).

## Finite photon–companion bridge

Issue #7 fixes the ontology and supplies a finite reversible source-to-receiver architecture containing a high-frequency photon, lower-frequency photon, neutral `chi`, converter recoil, dynamically computed bound `chi`, capture recoil, receiver, and one constrained frame. Every interaction has its Hermitian reverse, and the energy/component ledgers close to roundoff.

The finite benchmark uses the same bound source for synthetic clock response, massive-body acceleration, and light bending, with no independent motion or lensing multipliers. It also records an unresolved speed ratio:

```text
maximum chi speed / maximum photon speed = 0.3571428571
```

Thus common-cone recovery must be derived under issue #8. The large finite transition probability is an architecture stress test, not an astrophysical conversion rate. See [`companion_integration.md`](companion_integration.md).

## Current work program

| Issue | Pillar | Status |
|---:|---|---|
| #2 | Finite dressed frame Hamiltonian | Linear finite-regulator scope closed |
| #3 | Shared microscopic coefficients and `Z_g/Z_A` | Open |
| #4 | Chiral matter defects and protected gaps | Finite free-regulator prototype closed |
| #5 | Nonlinear gravity, self-coupling, many-body binding, and volume/vacuum | Open |
| #6 | Finite electromagnetic Coulomb/QED phase | Stage 6A static and Stage 6B pure-gauge dynamics pass; interacting matter/QED open |
| #7 | Photon–companion integration | Finite architecture scope closed |
| #8 | Continuum consistency, mirror completion, Lorentz/common-cone recovery | Open |
| #9 | Parameter ledger, frozen predictions, and empirical tests | Open |

The immediate next gate is to embed the minimal `|q|=1` issue-#4 endpoint in the Stage-6B detuned photon region and test exact current continuity, a finite Ward identity, photon dressing, vacuum polarization, and the issue-#7 conversion vertex using actual photon eigenmodes.

## Relationship to the wider repository

The root program remains governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, `research_plan/universe-contract.md`, and the solution-goal ledger. It requires fixed published distances, an operationally nonexpanding universe, complete energy accounting, no inserted dark halo, and joint predictions for redshift, timing, brightness, motion, lensing, gravitational waves, and background observations.

The finite constructions here supply architecture pieces. They do not yet derive a physical redshift coefficient, continuum conversion rate, many-body galactic deposit, universal weak-field correction, or observational fit.

## Stop rules

Do not call the RK equal-time tensor a photon spectrum or its winding free energy a photon gap. Do not call the Stage-6B principal-axis/full-symbol result a completed interacting QED continuum. Do not treat the finite matter charge set as observed particles, the companion transition probability as an astrophysical rate, one-particle binding as a halo, fixed volume as a cosmological-constant solution, or the finite gravity regulator as nonlinear quantum gravity. Derive shared coefficients, interacting Ward identities, nonlinear closure, common-cone recovery, regulator universality, and frozen predictions—or reject the combined branch.
