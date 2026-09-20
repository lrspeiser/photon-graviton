# Phase Junction / Photon–Companion Integration Contract

**Date:** 2026-09-20  
**Issue:** #7  
**Decision:** the finite integration-architecture gate is closed; continuum, many-body, nonlinear, and empirical gates remain open in their owning issues  
**Executable evidence:** [`microscopic/check_companion_bridge.py`](microscopic/check_companion_bridge.py)  
**Frozen output:** [`microscopic/companion_bridge_results.json`](microscopic/companion_bridge_results.json)

## Conclusion and claim boundary

The Phase Junction module and the active photon–companion program now use one declared ontology and one constructive finite Hamiltonian chain:

```text
finite source fuel
    <-> high-frequency photon
    <-> lower-frequency photon + neutral chi + converter recoil
    <-> lower-frequency photon + bound chi + capture recoil
    <-> receiver excitation + bound chi + capture recoil
    -> one constrained frame
    -> material-clock shift, massive-body acceleration, and light bending.
```

The companion is `chi = varphi`, the neutral relative-phase excitation. A deposit is a bound state of that same sector. The photon remains the compact-U(1) transverse quantum. The graviton remains the transverse frame excitation. The frame connection remains constrained and auxiliary.

This closes the ambiguity and common-ledger requirements of issue #7. It does **not** show that the branch scales to a healthy 3+1-dimensional continuum or to nature. In particular, it does not establish deconfined QED, chiral matter, nonlinear gravity, a common relativistic cone, a many-body galactic deposit, an astrophysical conversion rate, or agreement with observations.

## Canonical field dictionary

| Object | Selected meaning | Disallowed identification |
|---|---|---|
| Photon | Transverse quantum of the compact electromagnetic link; the high and low bins are frequencies of the same field | Not `chi`, a scalar deposit, or the graviton |
| Companion `chi=varphi` | Neutral scalar relative-phase quantum `theta_m-theta_g` | Not a photon, tensor graviton, frame-connection mode, or fitted halo |
| Deposit | Bound one- or many-`chi` state | Not a second energy reservoir counted in addition to its `chi` energy |
| Matter | Finite source, converter, capture/recoil, and receiver defects | The finite proxies are not yet Standard Model fermions |
| Recoil/phonon | Explicit matter degree of freedom carrying capture energy and momentum | Not an invisible sink |
| Frame `h` | One shared constrained geometry sourced by the deposit and ordinary matter | Not separate motion and lensing fields |
| Graviton | Transverse excitation of the frame | Not the captured companion reservoir |
| Frame connection | Constrained auxiliary comparison variable | Not a low-energy particle |
| `alpha` and `Gamma_cap` | Coarse-grained coefficients to derive from the microscopic theory | Not freely fitted constants per object |

## Finite common Hamiltonian

The executable model uses

\[
H=H_{\rm source}+H_\gamma+H_\chi+H_{\rm matter}+H_{\rm recoil}
 +H_{\rm emission}+H_{\rm conversion}+H_{\rm capture}+H_{\rm receiver}
 +H_{\rm frame}^{\rm red}.
\]

The transition terms occur with their Hermitian reverses:

\[
H_{\rm emission}=g_s(|H\rangle\langle S|+\mathrm{h.c.}),
\]

\[
H_{\rm conversion}=g_c(|\gamma_L,\chi, B\rangle\langle\gamma_H,A|+\mathrm{h.c.}),
\]

\[
H_{\rm capture}=g_b(|\gamma_L,\chi_{\rm bound},r\rangle
\langle\gamma_L,\chi_{\rm free}|+\mathrm{h.c.}),
\]

\[
H_{\rm receiver}=g_r(|R,\chi_{\rm bound},r\rangle
\langle\gamma_L,\chi_{\rm bound},r|+\mathrm{h.c.}).
\]

The bound state is the lowest eigenmode of a seven-site `chi` tight-binding Hamiltonian with a local matter well. Its energy lies `1.2027895096` below the free band and its inverse participation ratio is `0.8906794017`. The reduced constrained-frame term is included in the bound sectors:

\[
H_{\rm frame}[h;\rho]=\frac12 h^TKh-\kappa h^T\rho,
\qquad
H_{\rm frame}^{\rm red}=-\frac{\kappa^2}{2}\rho^TK^{-1}\rho.
\]

A static well alone cannot irreversibly capture a traveling excitation in a closed unitary system. The capture sector therefore includes explicit matter recoil, and release is the reverse matrix element rather than a separately invented process.

## Executed results

### Matter-assisted local conversion vertex

A translation-invariant seven-site calculation uses the sectors

```text
|A; photon(k), matter(p)>
|B; photon(k'), chi(q), matter(p')>.
```

The full Hilbert space has 392 states and the conserved crystal-momentum block has 56. The photon survives at lower energy. For the declared on-shell channel:

```text
photon energy loss       1.0820883461
chi energy               0.7012081585
matter recoil            0.2014204980
matter internal gap      0.1794596896
shell mismatch           0
```

The target conversion probability reaches `0.9743285451`; energy drift is `2.58e-14`; the translation commutator, charge commutator, selected-`Jz` commutator, and forward/reverse probability difference are zero at the reported precision.

### Source-to-receiver calculation

The common five-sector Hamiltonian evolves a finite source through emission, conversion, capture, and reception. At the declared best event:

```text
receiver + bound chi + recoil probability     0.9742646744
maximum bound chi + recoil probability         0.9778685459
energy-expectation drift                        9.77e-15
component-ledger error                          8.88e-16
forward/reverse probability difference          0
closed-boundary energy flux                     0
field-momentum / reciprocal-force residual      4.45e-5
```

The component ledger includes source fuel, photon energy, free and bound `chi`, converter/receiver matter gaps, recoil, signed interaction energies, and the reduced frame energy. The large transition probability is a resonant finite architecture stress test, not an astrophysical conversion rate.

### Shared motion, lensing, and clock response

The dynamically populated bound-mode density is the sole new source passed to the frame solve. The same solution gives:

```text
frame equation residual              2.78e-17
massive-body acceleration            -2.331869e-4
light deflection                     -3.904997e-3
same-field identity residual          0
material clock-gap shift              3.125778e-3
independent motion multipliers        0
independent lensing multipliers       0
```

This is a synthetic weak-field dictionary. It is not a 3D ray-bundle calculation or an observational lensing result.

### Explicit unresolved result

The finite dispersions give

```text
maximum chi group speed / maximum photon group speed = 0.3571428571.
```

Equal-speed propagation is therefore **not** established. Common-cone recovery must be derived under issue #8 or this branch fails.

## Issue #7 acceptance crosswalk

| Required by issue #7 | Evidence | Decision |
|---|---|---|
| One unambiguous companion identity | `chi=varphi`, with the photon, graviton, frame connection, and deposit separately defined | Pass |
| One common emitter/photon/companion/matter/frame Hamiltonian | Source-to-receiver finite Hamiltonian with exactly reduced frame term | Pass at finite architecture level |
| Energy, momentum, source-fuel, and boundary-flow ledger | Frozen component ledger, reciprocal force integral, finite source energy, zero closed-boundary flux | Pass |
| Reverse transition and recoil | Every coupling is Hermitian; matter recoil is explicit | Pass |
| Same stored field for motion and lensing | One bound density and one frame solve; zero independent multipliers | Pass as synthetic dictionary |
| No inserted halo/capture lifetime | Deposit is a computed bound eigenmode populated dynamically; no lifetime is claimed | Pass |
| Dependency matrix across the active program | Twelve-goal and R01–R32 maps below | Pass |
| End-to-end emitter-to-receiver calculation | Executed source-to-receiver trajectory | Pass |
| Disposition of incompatible prior branches | Explicit table below | Pass |

## Active twelve-goal ledger

No item in `research_plan/solution-goal-ledger.md` is marked fully achieved by this finite bridge.

| # | Goal | Bridge supplied | Still required |
|---:|---|---|---|
| 1 | Complete equations and stability | Complete finite Hermitian chain and constrained frame | 3+1D nonlinear equations, long-time stability, causal boundaries (#5, #8) |
| 2 | Remove directional artifacts | Exact lattice translation symmetry in the conversion vertex | Rotational covariance and regulator removal (#4, #8) |
| 3 | Isolated 3D evolution | Closed finite evolution and zero boundary flux | Isolated/open nonlinear 3D evolution (#5, #8) |
| 4 | Physical source budget | Finite source fuel and exact energy ledger | Real source histories, heat, capacity, and duration (#4, #9) |
| 5 | Durable swirl | No false swirl claim; scalar deposit is kept distinct | Any circulation/memory mechanism and angular-momentum transport (#5, #8) |
| 6 | Outer force and stars | One sourced-frame force dictionary | Extended profiles, stable orbits, vertical/noncircular motion (#5, #9) |
| 7 | Full lensing | Same frame controls matter and light | 3D images, shear, magnification, delays, and data (#5, #8, #9) |
| 8 | Light/local effects | Photon survives at lower energy; clock shift and reverse channel exist | Linewidth, polarization, coherence, clocks, equivalence and Solar-System tests (#4, #6, #8, #9) |
| 9 | Observation provenance | No observations used to select the finite parameters | End-to-end provenance for selected datasets (#9) |
| 10 | Shared held-out fit | No incompatible fitted halo introduced | One frozen model on held-out galaxies, clusters, and mergers (#9) |
| 11 | Contextual benchmarks | Claim boundary and parameter count are explicit | Matched comparisons with uncertainty and complexity (#9) |
| 12 | Reproducibility and attribution | Executable script, frozen JSON, controls, CI | Independent reproduction of the final continuum/data pipeline |

## R01–R32 dependency map

| Requirement | Finite bridge status / owner |
|---|---|
| R01 Operational definitions | Photon, `chi`, deposit, matter, frame, graviton, and auxiliary connection fixed here |
| R02 Microscopic conversion | Finite reversible matter-assisted vertex passes; continuum rate remains #4/#6/#8 |
| R03 Conservation | Finite energy, momentum-reaction, source-fuel, and boundary ledgers pass |
| R04 Quantum/classical consistency | Finite Hermiticity/unitarity only; interacting continuum remains #8 |
| R05 Graviton identity | Kept distinct from scalar `chi`; tensor completion remains #2/#5/#8 |
| R06 Source production | Finite source proxy only; stellar/AGN histories remain #9 |
| R07 Redshift | Lower photon energy exists; no distance law or line prediction yet (#8/#9) |
| R08 Time dilation | Not derived (#4/#8/#9) |
| R09 Brightness and photon counts | Photon survives in the selected channel; full transport absent (#6/#8/#9) |
| R10 Angular distances/surface brightness | Open (#8/#9) |
| R11 Spectral/image fidelity | Closed-system coherence only; linewidth, blur, polarization open (#6/#8/#9) |
| R12 Atomic clocks/cavities | Shared frame shift defined; microscopic clock calculation open (#4/#8/#9) |
| R13 Local gravity | Shared source dictionary exists; precision local tests open (#5/#8/#9) |
| R14 Capture/retention | Finite reversible bound-state seed passes; many-body capacity/lifetime open (#5/#8) |
| R15 Halo profiles | No halo inserted; derived extended profiles open (#5/#9) |
| R16 Galaxy rotation/scaling | Open (#5/#9) |
| R17 Other galaxy dynamics | Open (#5/#9) |
| R18 Lensing | Same-field rule fixed; 3D strong/weak predictions open (#5/#8/#9) |
| R19 Clusters | Open (#5/#9) |
| R20 Merging systems | Open (#5/#9) |
| R21 Environmental dependence | Matter-assisted dependence is allowed but not derived macroscopically (#4/#9) |
| R22 Stellar/thermal physics | Recoil is explicit; realistic heating/cooling open (#4/#9) |
| R23 Gravitational waves | Tensor sector separate; nonlinear production/propagation open (#2/#5/#8) |
| R24 Compact objects | Open (#5/#8/#9) |
| R25 Microwave-background spectrum | Open (#8/#9) |
| R26 Microwave-background angular structure | Open (#8/#9) |
| R27 Structure/BAO | Open (#5/#8/#9) |
| R28 Abundances/chronology | Open (#9) |
| R29 Static background/topology | Open (#5/#8/#9) |
| R30 Global radiation/entropy | Finite ledger only; global fuel/entropy open (#8/#9) |
| R31 Statistical identifiability | No per-object parameters used here; full ledger/fits open (#9) |
| R32 Distinguishing predictions | None frozen yet (#9) |

## Disposition of prior and parallel branches

| Branch/object | Disposition |
|---|---|
| RC-1 reversible two-mode oscillator | Retained as an effective precursor, not a second companion identity |
| Direct conversion and bound companions | Selected branch; its empirical `alpha` and `Gamma_cap` remain quantities to derive |
| Cumulative-time/additional `q`-field cause | Archived alternative; cannot be combined without an explicit common action and revalidation |
| Photon as companion | Rejected for this branch; the photon survives and is distinct |
| Tensor graviton as deposit | Rejected; the tensor mode remains the gravity mediator |
| Frame-connection relative mode as companion | Rejected; it is removed/gapped by the second-class lock |
| `X/Y`, rotor, local-transfer, and common-cone fields | Effective candidates only until an action-level reduction identifies them with `chi` or the frame |
| Imposed halo or independent deposit profile | Disallowed |
| Independent motion/lensing multipliers | Disallowed |
| Added screening or clock factors to protect a fit | Inactive unless derived from the common theory |
| Fixed volume as a cosmological-constant solution | Not established; remains an issue-#5 linear global constraint candidate |

## Falsification and next ownership

The finite architecture must be rejected or revised if it cannot:

1. embed in the finite deconfined electromagnetic phase (#6);
2. replace proxy defects with gauge-covariant microscopic matter (#4);
3. derive its coefficients with the electromagnetic and frame sectors from one move set (#3);
4. recover a common causal cone without imposing it (#8);
5. support many-body bound states with finite capacity, acceptable lifetime/release, and nonlinear self-gravity (#5/#8);
6. satisfy local clock, spectral, coherence, equivalence, and gravity constraints (#8/#9);
7. freeze at least one cross-sector prediction before broad astrophysical fitting (#9).

The issue-#7 architecture is solved. The physical theory and empirical program are not.
