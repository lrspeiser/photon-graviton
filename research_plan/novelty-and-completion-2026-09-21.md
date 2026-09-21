# Novelty decision and completed next steps — 21 September 2026

## Decision

**The current evidence does not justify saying this is an original, completed photon-graviton theory.** SM-1's static equation is an instance of the established AQUAL framework with a chosen screening function. Its unscreened square-root relation has an exact published antecedent. The existing repository audit had already established these facts; the preceding SM-1 summary should have foregrounded them rather than making implementation progress sound like a new foundational theory.

Keep SM-1 as a frozen phenomenological benchmark. Do not spend the next phase selecting interpolation functions until a galaxy score improves. The research contribution must be a specified microscopic mechanism, its derived macroscopic response, and a discriminating prediction. Recovering a known effective equation would be acceptable, but choosing that equation first and then calling it derived is not.

This decision does not refute companion fields, source-generated memory or the wider alternative-gravity programme. It distinguishes known equations, newly executed project-specific calculations, and unproved physical claims.

## 1. Equation-level provenance

### Exact known formula

Famaey and Binney (2005), equation 5, print the Bekenstein toy interpolation

    mu_B(x) = [sqrt(1+4x)-1]/[sqrt(1+4x)+1].

PM-2A uses

    mu_PM(x) = 4x/[1+sqrt(1+4x)]^2.

These are algebraically identical. With y=g_N/a0, their spherical inverse is

    g/a0 = y+sqrt(y).

Thus neither g=g_N+sqrt(a0*g_N), rationalizing mu, nor the resulting asymptotic v^4=GMa0 is a new discovery of this project. The published function's equivalence does not imply that our full dynamical model is TeVeS; its nonspherical fields, other sectors and relativistic matter couplings would have to match too. The source paper itself gives regime and geometry qualifications.

### Exact known static framework

Bekenstein and Milgrom (1984), equations 2b-3, derive

    div[mu(|grad Phi|/a0) grad Phi] = 4*pi*G*rho

from an aquadratic gradient functional. SM-1's mu is a particular choice in this framework. Deriving a compatible energy by integrating mu and checking its Hamiltonian derivatives are useful consistency/implementation steps, not new principles.

The exact SM-1 denominator 1+(g_N/a0)^2 is a phenomenological choice whose exhaustive priority has not been established. Even a previously unused denominator would be a new interpolation within a known framework, not by itself a new fundamental theory.

### Known architecture, not a unique meaning of memory

Eliminating dynamical variables to obtain a memory kernel is standard projection mathematics. Time-nonlocal response, coupled source/field models, and gravitational models with history dependence have prior art. A local Markovian system for (U,U_t) is not a new fundamental non-Markovian law merely because its field persists after a source changes. A physical prediction for memory must specify which state was eliminated, its initial state and its coupling.

The earlier in-repository audit remains authoritative for its scope: [prior-art/audit.md](prior-art/audit.md). Primary literature inspected for this update includes Bekenstein and Milgrom (1984), Famaey and Binney (2005), and Zwanzig (1973). This is a bounded scientific literature/equivalence check, not exhaustive proof of global priority or a legal novelty opinion.

## 2. Completed SM-2: actual source-based galaxy equations

**Executed, not just proposed.** [Protocol](../research_work/results/screened-memory/protocol-sm2.md), [runner](../research_work/results/screened-memory/run_sm2.py), and [results summary](../research_work/results/screened-memory/results-sm2-summary.json).

Scientific baseline: abd026f. Protocol commit: 284f3bb. Executed code: ad7bde6. GitHub Actions run 35658900397 completed successfully. Success means the execution completed, not that all accuracy gates or the physical theory passed.

All 149 archived galaxies and 3,152 positive-radius rows are retained. Each inferred density is constructed once from the full source inputs. Four equations share its exact grid source: Newtonian, unscreened Bekenstein-root AQUAL, SM-1 screened AQUAL and simple-mu AQUAL. The root and screened scale stays 6.54e-11 m/s^2; the simple-mu control retains its previously selected 8.563e-11. No parameters were fitted here. These are exposed data, not a new blind test.

There are **1,341 field solves**: four equations at each of two grids and a screened outer-boundary control per galaxy. All solve convergence gates pass; the largest residual is 3.09e-10. All source-mass gates pass; the largest relative mass discrepancy is 2.24e-5. All predicted radial forces are positive at the scoring rows.

**Accuracy is not identical to convergence.** The fixed coarse/fine force gate is 2 percent and the screened outer-boundary gate is 0.5 percent. Eleven galaxies fail refinement and eleven fail the boundary test, with two overlapping. Therefore **129 galaxies are commonly certified and 20 remain numerically unresolved**. Both the full provisional sample and the common certified subset are reported; no unfavorable galaxy was silently removed.

### Mean per-galaxy velocity RMSE, km/s

Lower is better. All columns compare the same sources, rows and weighting within their sample.

| Equation | All 149, provisional | Common certified 129 |
|---|---:|---:|
| Same-density Newtonian | 47.48 | 46.30 |
| Published unscreened Bekenstein-root AQUAL | 17.28 | 16.18 |
| **SM-1 screened AQUAL** | **19.33** | **18.16** |
| Simple-mu AQUAL comparison | 17.26 | 16.26 |

On the full sample, screened minus root is +2.058 km/s. The paired 10,000-draw galaxy-bootstrap interval is [+0.573,+3.560] km/s. On the certified subset it is +1.981 km/s, interval [+0.557,+3.472]. These intervals describe sampling variation among the exposed galaxies under the frozen source assumptions; they do not include distance, inclination, mass-to-light, density-reconstruction or other astrophysical systematics, and are not a general theory-exclusion confidence.

Full-sample split means for root / screened / simple are:

- training 89: 16.878 / 19.420 / 17.376 km/s;
- validation 29: 22.230 / 22.406 / 20.965 km/s;
- archived test 31: 13.782 / 16.211 / 13.450 km/s.

The actual nonspherical field equation does not rescue this particular screening choice under the stated source model. It remains substantially better than baryons-only Newtonian gravity in this benchmark, but worse overall than the two known modified-gravity controls. That is not evidence that all screening or companion mechanisms fail.

The density reconstruction is itself an assumption: its Newtonian velocity differs from the tabulated SPARC baryonic velocity by a mean per-galaxy RMS of 4.511 km/s. The corresponding Newtonian observational RMSE changes from 45.571 to 47.475 km/s. Do not attribute every change from SM-1's old algebraic table to nonlinear field geometry or screening; source reconstruction also changed. Within SM-2, this source is matched across all four equations.

## 3. Completed CM-1: memory derived from the existing microscopic Hamiltonian

[Declaration](../phase_junction_network/microscopic/protocol_companion_memory.md), [executable verifier](../phase_junction_network/microscopic/check_companion_memory.py), [results](../phase_junction_network/microscopic/companion_memory_results_summary.json). Executed commit 58b91f7; GitHub Actions run 35659502279 completed successfully.

The unchanged five-state construction contains source fuel, a high-energy photon, traveling companion, bound companion plus recoil, and a receiver plus bound companion/recoil. This remains a finite synthetic model with chosen couplings, not a spatially resolved universe. CM-1 inserts no memory time or damping. It extracts the original Hamiltonian construction and checks its source blob hash before proceeding.

After removing the common diagonal energy 5, retain the two bound-containing states p and eliminate the three source/photon/traveling states q. Partition the original Hamiltonian as A=H_PP, B=H_QQ, V=H_PQ. Exact elimination gives

    i p_dot(t) = A p(t) + V exp(-i B t) q(0)
                 - i integral_0^t K(t-s) p(s) ds,
    K(t) = V exp(-i B t) V_dagger.

The initial source term is essential: dropping it removes the original fuel/initial-state information. This is standard exact projection, not an invented new formalism.

For the unchanged source/conversion/capture couplings a=.19, b=.17, c=.15, the only nonzero kernel element is

    K_00(t) = c^2 [a^2+b^2 cos(Omega*t)]/(a^2+b^2),
    Omega = sqrt(a^2+b^2) = 0.25495097567963926.

Numerically:

    K_00(t) = 0.012496153846153846
              + 0.010003846153846155*cos(0.25495097567963926*t).

Its period is 24.6446803760 in the original model's time units. It has a nondecaying constant plus a coherent oscillation. **There is no derived exponential decay time.** This kernel is nonlocal in the reduced description, while the complete five-state unitary system remains local in time.

### Independent verification paths

The closed kernel agrees with matrix exponentiation to 1.74e-17. A separate Schur-complement resolvent reproduces the projected full resolvent. Independent Gauss-Legendre time convolution reconstructs eliminated amplitudes to 6.52e-14 and the retained equation to 9.33e-15. Omitting the initial-source term produces a 0.128 residual, so that control demonstrably fails. Full norm, energy and component-ledger checks remain below 1e-14.

### The physical distinction exposed by the calculation

Starting with source fuel, the old bound-sector peak is reproduced: **97.7869 percent at time 87.55**. But it falls below **9.8337 percent by time 98** and reaches a sampled minimum of 2.255e-7 in bound probability later in the unchanged evolution. The source return probability reaches **99.9991 percent** near time 750.15.

The exact infinite-time mean bound probability is nevertheless **40.8785 percent**. Thus it would also be wrong to claim that the system has no bound occupancy or no memory. The correct statement is that this closed finite construction gives coherent exchange and recurrence, not a demonstrated irreversible capture process or a persistent deposit with a derived physical lifetime.

All times are dimensionless model times. No conversion to years, galaxy lifetime, merger offset or astrophysical rate has been derived. The long-time extrema are sampled on [0,2200] with spacing .05, not certified global continuous-time extrema.

This is a project-specific calculation that moves the microscopic-to-effective connection forward: the response is now calculated from the existing Hamiltonian rather than hypothesized. It does not identify this companion with SM-1's AQUAL scalar and does not derive a0, a galaxy density profile or a gravitational force law.

## 4. Independent archive readback

The downloaded artifacts underwent a further 1,682-check audit. It independently reconstructs all SM-2 velocities, per-galaxy errors, refinement measures, classification flags and score tables from the saved arrays, and checks 291 CM-1 state/time points by a separate matrix exponential. No audit check failed; the latter probability discrepancy is 3.79e-13. This is not another run of every PDE solve or the entire historical research suite.

The full SM-2 artifact is 10667515859, SHA-256 90fba0a209400acfbdaddda88e6e8044f321a40f30fac5ff22e66d2dc5845bc9. The full CM-1 artifact is 10666995125, SHA-256 f0e3b6ea78122a209a1b2286174b4089903538f1343c8b35d8b3b56fc05d296c. Both were downloaded into the originating conversation, preserving full arrays beyond dependence on Actions retention. Compact results and code are committed. Actions retention is finite; the full arrays should not be assumed permanently hosted by GitHub.

## 5. Completion path: three physical deliverables, not more interpolation choices

The following work is still open; it is not described as completed by SM-2 or CM-1.

### A. One spatial microscopic conversion/capture/recoil model

Replace the abstract recoil/capture sectors with explicit propagating degrees of freedom and compute conversion, capture, escape and release from one Hamiltonian. Use the SAME companion spectrum for traveling and bound states and include source fuel, reverse reactions, local energy/momentum transport and finite capacity. Do not create permanent capture by deleting the reverse matrix element or by adding an unaccounted sink.

Success requires controlled spatial/volume/time refinement, a measured outward energy/momentum ledger, and either a physically justified persistent bound fraction or a derived response kernel over the required regime. Closed-box recurrences, coherent oscillations, radiative escape and a thermodynamic-bath approximation must be distinguished. A useful control is source removal or exhaustion: the model must predict what is retained and for how long without a separately inserted lifetime. Failure of one spatial realization rejects that realization under its assumptions, not every memory mechanism.

CM-1 now gives an exact finite baseline and a working projection procedure for this step. It cannot supply the missing spatial interaction or physical spectrum by renaming its three eliminated levels a bath.

### B. Derive the force response from that microscopic state

Calculate the companion/deposit contribution to the same sourced frame, then obtain the response to ordinary-matter and radiative histories. The static force, memory kernel, characteristic acceleration and propagation coefficients must be outputs of the stated Hamiltonian and state, not independent selected functions. A known effective limit is permitted; the contribution would be its microscopic derivation and additional constrained predictions.

The existing SM-1 equation has no established derivation from the finite bridge. Its ad hoc time-dependent scalar extension has squared principal speeds c_*^2*mu transverse to a nonzero background gradient and c_*^2*(mu+z*mu') along it. These follow by differentiating the static flux but do not fix c_* or damping. At g_N/a0=1, mu=2/3 and mu+z*mu'=4/3; a single chosen c_* therefore does not produce one direction-independent scalar speed. This does not constrain a separate tensor sector by itself. It shows why static ellipticity cannot stand in for a derived common light/gravity cone or an allowed scalar-radiation rate.

### C. One fixed model across the observations that distinguish it

Use one matter/light coupling and one source prescription for galaxy rotation, vertical dynamics, lensing, Solar-system external-field response and binary/scalar radiation. For a memory model, add a predeclared source-history comparison or merger/wake observable. The 149 SPARC galaxies here are already exposed; a new confirmation claim needs independent data or genuinely withheld observables.

The Solar calculation must solve this candidate's actual AQUAL external-field problem; a QUMOND integral or an isolated monopole suppression is not its substitute. Park et al. (2026), revised 28 July, report a Cassini quadrupole constraint Q2=(1.6+/-1.8)e-27 s^-2. Neither SM-1 nor SM-2 has calculated its prediction for that quantity. Neither run supplies a relativistic lensing/binary completion.

## Working research rule

Known mathematics is allowed and should be credited. Reproducing its consequences is a baseline. A novelty claim requires a specified difference in mechanism or physical predictions, not a different name for the same equation. A completion claim requires the physical links between sectors, not merely green implementation tests in separate toys.

**Current outcome:** a documented originality boundary; a completed 149-source field benchmark that disfavors the chosen screening relative to the controls under fixed assumptions; and an exact microscopic memory response showing precisely why a high transfer peak is not yet a durable deposit. The next construction is the spatial mechanism and its derived force response, not another arbitrary acceleration formula.

## Primary literature

- Bekenstein, J. and Milgrom, M. (1984), Does the missing mass problem signal the breakdown of Newtonian gravity?, Astrophysical Journal 286, 7-14; equations 2b-3.
- Famaey, B. and Binney, J. (2005), Modified Newtonian Dynamics in the Milky Way, MNRAS 363, 603-608, arXiv:astro-ph/0506723; equation 5.
- Feshbach, H. (1962), A unified theory of nuclear reactions. II, Annals of Physics 19, 287-313; effective projection formalism.
- Zwanzig, R. (1973), Nonlinear generalized Langevin equations, Journal of Statistical Physics 9, 215-220, DOI 10.1007/BF01008729.
- Park, R. S., Hees, A., Famaey, B., Desmond, H. and Durakovic, A. (2026), Improved constraints on modified Newtonian gravity from Cassini radio tracking data, arXiv:2602.17884v2.
