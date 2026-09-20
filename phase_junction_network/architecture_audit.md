# Phase Junction Network: Big-Picture Architecture Audit

**Audit date:** 2026-09-20  
**Updated after issue #2, issue #4, and issue #7 finite-gate closures:** 2026-09-20  
**Scope:** determine which elements of a candidate fundamental theory are actually established, which remain regulator prototypes, and what work now controls success or failure.

## Executive conclusion

The project is no longer dominated by a single successful gravity calculation. It now has three substantial finite architecture results:

1. **Issue #2 — linear finite gravity regulator:** an exact finite scalar/vector quotient, written in the dressed frame Weyl algebra, leaves two positive linearly dispersing tensor branches without a transverse-traceless projector.
2. **Issue #4 — finite charged-matter prototype:** odd-strand endpoint defects, exact Gauss-covariant hopping, a primitive anomaly-free chiral charge set, and a protected domain-wall gap hierarchy pass their declared free-regulator gates.
3. **Issue #7 — finite photon–companion integration architecture:** one companion identity, reversible matter-assisted conversion, recoil, binding, reverse transitions, energy/momentum ledgers, and one shared frame for motion and lensing are explicit.

These are meaningful closures. They do **not** make the framework a complete theory. The decisive open risks are now:

- whether the actual finite electromagnetic model has a deconfined two-photon Coulomb phase with dynamical defects;
- whether all sector coefficients come from one microscopic move set rather than independent choices;
- whether nonlinear gravity preserves the constraints and universally includes field self-energy;
- whether the matter mirror can be removed or symmetrically gapped in a healthy interacting theory;
- whether photons, companions, tensor modes, and matter recover one causal cone;
- whether different finite regulators flow to one unitary or reflection-positive continuum theory;
- whether the companion has a viable many-body capacity, lifetime, release, and self-gravity mechanism;
- whether the theory makes a frozen cross-sector prediction that survives data.

Linear gravity should now be treated as a regression-tested input, not the dominant research program. The clock normalization, gravity scale, and compact completion used in issue #2 remain regulator choices until issues #3 and #8 derive or universalize them.

## Relationship to the wider repository

The active program is governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, `research_plan/universe-contract.md`, and the solution-goal ledger. It requires fixed published distances, an operationally nonexpanding universe, complete energy accounting, no inserted dark halo, and joint predictions for redshift, timing, brightness, motion, lensing, gravitational waves, and background observations.

The current microscopic dictionary is:

| Object | Selected role |
|---|---|
| Photon | Transverse excitation of the compact U(1) link connection |
| Companion `chi=varphi` | Neutral matter–geometry relative-phase excitation |
| Deposit | Bound many-body state of the same companion sector |
| Graviton | Transverse dressed-frame excitation |
| Frame connection | Constrained auxiliary comparison variable with a gapped relative band |
| Matter | Finite charged or neutral endpoint defect |

Issue #7 now supplies a finite reversible bridge between photons, companion excitations, matter recoil, a bound mode, and a shared frame. It does not supply a continuum astrophysical conversion rate, a galaxy-scale deposit, a common propagation cone, or empirical success. The finite event probabilities are architecture stress tests, not observed rates.

## Master gate matrix

Status meanings:

- **Established internally:** multiple committed checks establish the stated result under explicit assumptions.
- **Finite prototype pass:** an explicit finite construction passes its declared algebraic or free-regulator gates while physical completion remains open.
- **Partial:** a target structure exists but the physical phase or complete mechanism is not demonstrated.
- **Open:** no construction closes the requirement.

| Pillar | Minimum requirement | Current status | Established evidence | Decisive missing work |
|---|---|---|---|---|
| Ontology and companion integration | One field dictionary and one reversible source-to-storage-to-gravity Hamiltonian | **Finite architecture pass** | Issue #7: one companion identity, Hermitian reverse channels, recoil, finite bound mode, energy/momentum ledgers, shared frame | Replace proxies with full finite matter/EM sectors; derive physical rates, many-body storage, lifetime, release, and common cone |
| Finite electromagnetism | Deconfined 3+1D Coulomb phase, two transverse photons, dynamical charge, `1/r`, Ward identities | **Partial** | Exact finite Gauss symmetry; spin-1 is the smallest tested link with nonconstant electric energy; issue-#4 defects are Gauss covariant | Phase diagram and deconfinement with defects, transverse spectrum, static force, scaling, Ward identities; issue #6 |
| Linear gravity and finite regulator | Four constraints, exactly two positive linear tensor modes, finite positive Hamiltonian/transfer regulator, no TT insertion | **Established internally** | Continuum and exact first-order identities; local/full-real-space reductions; finite Weyl lock; issue-#2 quotient and compact spectra at `p=5,7,11` | Preserve as regression input; derive coefficients in #3, nonlinear completion in #5, universality in #8 |
| Nonlinear gravity | Closed nonlinear constraints, universal self-coupling, strong-field stability, no extra scalar/ghost | **Open** | Linear first-/second-class structure; finite linear regulator; fixed-volume linear control | Nonlinear frame/connection action, self-energy source, closure, strong-field spectrum; issue #5 |
| Matter | Finite charged fermionic defects with chirality, anomaly control, protected gaps, clocks, and healthy continuum | **Finite prototype pass** | Issue #4: odd-strand endpoint, exact hopping, `(-11,-5,-1,-1,9,9)`, one cone per wall, protected gaps, common bare frame derivative | Observed or declared reduced spectrum, mirror completion, interactions, bound-state clocks, radiative stability; issues #3, #6, #8 |
| Shared parameters | One microscopic move set deriving `U_A,K_A,U_g,K_g`, matter localization, and companion couplings with fewer inputs than outputs | **Open** | Current formulas expose required relations; issue-#2 `lambda_g`, issue-#4 `r0,eta`, and issue-#7 couplings are explicit inputs | Derive coefficients and at least one held-out cross-sector relation; issue #3 |
| Common causal geometry | One cone and one frame coupling for photons, companions, tensor modes, and matter | **Partial** | Shared frame variables; matter species use one bare frame derivative; issue #7 uses one frame for motion and lensing | Finite EM phase and interacting cone recovery; finite bridge currently gives `v_chi/v_gamma=0.3571428571`; issues #6 and #8 |
| Quantum continuum consistency | Unitary/reflection-positive interacting limit, positive residues, Ward identities, anomaly control, regulator universality | **Open** | Reduced free tensor transfer matrix; positive compact tensor regulator; finite matter topology/anomaly checks | Interactions, mirror completion, renormalization, Lorentz recovery, regulator comparison, no-go audit; issue #8 |
| Vacuum and background | Stable vacuum, derived volume term, thermodynamic/statistical state, admissible background | **Partial at one linear global mode** | Fixed-volume plus trace-momentum pair removes the homogeneous conformal pair at linear order | Microscopic origin, local vacuum energy, nonlinear closure, background solutions; issue #5 |
| Empirical distinctiveness | Fewer calibration inputs than outputs and at least one frozen cross-sector prediction | **Open** | Parameter roles are becoming explicit; flexible fitting remains prohibited | Parameter ledger, observable dictionary, held-out prediction, staged local-to-cosmological test; issue #9 |
| Reproducibility | Frozen outputs, controls, CI, finite-size tests, retained failures | **Established internally** | Executable issue-#2, #4, and #7 gates; manifests; negative branches; GitHub Actions | Keep every new claim executable and synchronized |

## Issue #2 resolution: finite dressed-frame gravity

### Construction

The finite connection lock enforces the relative variables through a positive Weyl oscillator while preserving the exact dressed frame algebra

\[
\overline Z_i=Z_{h_i},
\qquad
\overline X_i=X_{h_i}\prod_aX_{C_a}^{A_{ai}}.
\]

The scalar and vector finite stabilizers are reduced over `GF(p)`. The code computes—not assumes—the quotient spaces

\[
\ker G/\operatorname{row}C,
\qquad
\ker C/\operatorname{row}G,
\]

then symplectically dualizes them. For every tested odd prime `p=5,7,11`, the quotient contains exactly two canonical logical Weyl pairs.

For each physical polarization, the positive compact regulator is

\[
H_{p,L}=\frac{\lambda_g}{2}\sum_z
\left[(2-X_z-X_z^\dagger)+
(2-Z_{z+1}Z_z^\dagger-Z_zZ_{z+1}^\dagger)\right],
\]

and the tensor Hamiltonian is two identical copies. Equivalently, the parent local evolution is a constrained transfer step with local group averaging over the scalar and vector stabilizers.

This route does not contradict the retained derivative-order obstruction. A frame-only Hamiltonian made from squares of the lowest exact local invariants still gives `omega~k^3` and remains rejected. Locality belongs to the parent constrained transfer step; the physical logarithm acts on the quotient Hilbert space.

### Finite results

- finite primes: `5,7,11`;
- tested lattice sizes: `p=5: L=4–8`, `p=7: L=4–7`, `p=11: L=3–6`;
- exactly two logical Weyl pairs per nonzero propagation block;
- exact two-polarization first-gap degeneracy;
- gap powers between `1.02690` and `1.07049` in `Delta proportional to |k_hat|^s`;
- maximum residual of the linear-plus-cubic fit: `1.663e-3`;
- minimum connection-lock/tensor-gap ratio: `2.491`;
- maximum lock/frame commutator: `1.67e-16`;
- ground-band leakage norm: `5.16e-16`;
- `p=11` tensor compact-cut bond fraction: `1.89e-4`;
- complete three-dimensional relative frequency error: `2.15e-16`;
- polarization split: `2.84e-16`;
- cubic-lattice anisotropy recovery: approximately `L^-2.0166`.

### Claim boundary

Issue #2 closes the **linear finite constrained regulator**. It does not establish:

- nonlinear quantum gravity;
- a derived Newton constant or shared gravity/photon impedance;
- a common photon/companion/matter/tensor cone;
- regulator universality or a unique ultraviolet completion;
- a vacuum-energy mechanism;
- a new empirical prediction.

Those dependencies are owned by issues #3, #5, #6, #8, and #9.

## Issue #4 and issue #7 in the combined architecture

The matter prototype now supplies a real finite charged participant rather than a continuum placeholder. The companion bridge supplies a reversible interaction architecture rather than a variable-name mapping. The next task is not to repeat either standalone demonstration. It is to embed all three finite pieces—the issue-#2 frame regulator, issue-#4 endpoint, and issue-#7 conversion architecture—inside the actual issue-#6 electromagnetic phase and derive their couplings jointly under issue #3.

The combined model fails if the matter endpoint, conversion vertex, photon phase, or frame regulator each requires an independently selected cone or response multiplier.

## Main strategic finding

The project is still not in a general tuning stage. Three finite architecture gates pass, but the remaining uncertainties are structural:

- deconfinement may fail once dynamical matter is included;
- the issue-#2 compact completion may not share an interacting universality class with alternative regulators;
- nonlinear gravity may regenerate forbidden modes;
- the matter mirror may not admit a healthy completion;
- the finite companion speed mismatch may survive the scaling limit;
- shared microscopic coefficients may not close;
- many-body companion storage may lack adequate capacity, lifetime, or energy supply;
- no frozen prediction yet distinguishes the theory from QED plus general relativity with selected coefficients.

Further decimal precision on any already-passed finite benchmark does not address these risks.

## Correct work order

### Priority 0: finite electromagnetic phase with dynamical matter

Issue #6 must demonstrate or reject a deconfined spin-1 Coulomb phase containing the actual issue-#4 endpoints. Required outputs include a phase diagram, two transverse branches, expanding-range `1/r`, charge propagation, finite-representation scaling, and Ward identities.

### Priority 0: shared microscopic coefficients

Issue #3 must derive, from one move set,

\[
U_A,\ K_A,\ U_g,\ K_g,\ r_0,\ \eta,
\]

plus the issue-#7 conversion, binding, and recoil coefficients. A common limiting speed constrains products but does not by itself fix the impedance ratio. The issue-#2 clock point and `lambda_g` cannot be promoted from regulator choices to predictions without this derivation.

### Priority 1: nonlinear gravity and many-body binding

Issue #5 owns nonlinear constraint closure, gravitational self-energy, strong-field stability, the volume/vacuum term, and the gravitational source of many-body companion deposits.

### Priority 1: interacting continuum and common cone

Issue #8 owns mirror completion, unitarity/reflection positivity, Ward identities, anomaly accounting, radiative stability, Lorentz recovery, common-cone recovery, regulator universality, and the relevant emergent-spin-2 no-go assumptions.

### Priority 2: frozen predictions and data

Issue #9 should maintain the parameter and observable ledger now. Empirical unblinding waits until the common action and calibration roles are frozen. Local propagation, free fall, light bending/delay, clock response, and gravitational-wave tests precede galaxy or cosmological flexibility.

## Dependency graph

```text
   issue #2 finite frame regulator   issue #4 matter prototype   issue #7 bridge
                \                         |                         /
                 +------------------------+------------------------+
                                          |
                              finite EM phase #6
                                          |
                             shared coefficients #3
                                          |
                    nonlinear gravity / many-body binding #5
                                          |
                interacting continuum / common cone / mirror #8
                                          |
                         frozen predictions and data #9
```

Exploratory work may proceed in parallel, but downstream claims cannot close while upstream coefficients or phases remain arbitrary.

## Stop rules

1. Do not reopen issue #2 merely for smaller residuals; reopen only for a failed regression or a concrete dependency from #3, #5, #6, or #8.
2. Do not call the issue-#2 finite linear regulator nonlinear quantum gravity or a unique ultraviolet theory.
3. Do not identify the issue-#4 charge set with Standard Model particles without deriving the mapping and additional gauge structure.
4. Do not hide or silently discard the remote matter mirror wall.
5. Do not infer a deconfined QED phase from exact Gauss-covariant hopping alone.
6. Do not treat issue-#7 finite transition probabilities as astrophysical rates or one-particle binding as a halo.
7. Do not impose equal photon and companion speeds after the finite benchmark found a mismatch; derive common-cone recovery or reject the branch.
8. Do not retune photon, companion, matter, and gravity sectors independently to obtain the same speed or response.
9. Do not call the linear fixed-volume control a cosmological-constant solution.
10. Do not fit galaxy or cosmological data before local predictions and parameter roles are frozen.
11. Do not count recovery of QED and general relativity after freely selecting all coefficients as a distinct prediction.

## Decision

The finite linear gravity gate is complete. The matter endpoint and companion bridge also pass their declared finite architecture scopes. The program should now be evaluated by whether these pieces can coexist in one deconfined, nonlinear, interacting, universal theory with fewer microscopic inputs than low-energy outputs.

The next claim threshold is:

> One finite local model contains a deconfined two-polarization photon sector, the two-helicity dressed-frame sector, a healthy protected charged-matter sector, and the declared companion sector; it derives shared coefficients, preserves constraints nonlinearly, reaches one common relativistic continuum, and produces at least one frozen cross-sector prediction.

Nothing currently committed meets that full threshold. Closing issue #2 removes one major uncertainty and makes the remaining integration and continuum risks more sharply testable.
