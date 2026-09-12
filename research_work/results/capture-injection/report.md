# Does capture inject bound material with the required motion?

Under one explicit local-merger rule, capture produces inward motion and zero tangential angular momentum. Weak absorption can inject material that is initially bound, but does not automatically build the circular-orbit distribution from the previous retention calculation. More strongly attenuated examples mostly inject unbound material at the tested weak-gravity scales. This identifies a missing formation mechanism; it does not rule out capture with a momentum recipient, cooling, scattering or a different stored field state.

## Assumptions and equation provenance

Use the existing spherical p6,C1,R100 deposit snapshots, with kappa0=.1,10,1000 and their gravity unchanged. Evaluate the instantaneous incoming source at each snapshot; do not evolve the density. The deposited profiles are conditional inputs, not self-consistent outputs of this new formation rule. Traveling companions are assigned momentum E/c, as photons have. All captured directions in a small region are assumed to merge into one cold massive cohort, with no external impulse and no outgoing energy. This merger is a new candidate postulate for the project, not an established graviton interaction or a derived microscopic transition.

Special-relativistic energy/momentum composition and spherical radiation continuity are known mathematics. Their application below is a conditional derivation. No uniqueness is claimed for the formulas or the resulting source model.

Let J be angle-integrated incident intensity at a point and F_r its radial energy flux. Scalar opacity gives absorbed energy q=kappa J and absorbed radial momentum rate kappa F_r/c. If E and P_r denote a newly captured cohort's total local energy and momentum, then:

`beta=|v_r|/c=c |P_r|/E=|F_r|/J`,

`M_rest c^2=E sqrt(1-beta^2)`,

`K_bulk=E [1-sqrt(1-beta^2)]`.

These satisfy E=M_rest c^2+K_bulk and P_r=gamma M_rest v_r. Treating all absorbed E/c^2 as stationary rest mass would silently discard this bulk kinetic energy and momentum. A cold cohort is a coarse-grained assumption; the calculation does not establish a fixed particle species, particle-number law or mass spectrum.

For stationary spherical absorption with no internal emission:

`F_r=-(1/r^2) integral_0^r q(u) u^2 du`.

Thus F_r is inward. Spherical symmetry cancels local mean tangential momentum under this merger rule, so each resulting cold cohort has j=0. This cannot create the nonzero individual circular angular momenta of the previous supported snapshot. Zero *net* angular momentum would not forbid oppositely oriented circular orbits; the restriction here comes from merging all local directions into a single cold flow. A different interaction rule could preserve velocity dispersion or transfer angular momentum.

## What is tested

The [protocol](protocol.md) was written before computation. Existing512/1024-shell snapshots supply opacity and gravity. Absorption along chords gives q; spherical continuity gives F_r. Independent pointwise angular ray integrals at r=.1,.3,1,3,10 check the inferred beta, with256/512 angle nodes. This tests the energy-to-momentum conversion calculation through two different integrations.

Initial binding is assessed using the weak-field escape condition `beta^2 < 2 Psi/(c/v0)^2`, where v0=sqrt(G M_b/a). The ratios c/v0=300,1000,3000 are sensitivity choices, not measured or fitted model parameters. The largest Psi/c^2 is below1.3e-5. Escape classification is a weak-field diagnostic: accepted particles are slow compared with c, while a high-beta rejection does not justify evolving such particles with nonrelativistic dynamics. No scattering, later energy loss or changing potential is included.

The reported fractions weight the instantaneous absorbed-energy source over the finite R100 domain, not all previously deposited mass or all future arrivals. The finite grid can miss a minute bound region at the origin, where exact symmetry gives beta=0. Consequently a tabulated zero means no selected shell center passed the threshold, not proof of exactly zero bound material. The first shell's total source fraction is retained as a bound on that particular unresolved central contribution; it is not a global numerical error bar or a statement about material outside the domain.

## Formation implications

Initial binding and support are different requirements. Bound cold radial cohorts can fall through the well and change its density; they do not remain at capture positions. Their subsequent trajectories, phase mixing and self-gravity must be evolved. Angular-momentum redistribution, interaction with an existing well constituent, emission of energy, or a retained field stress are explicit alternative physical ingredients, not automatic consequences of the current equations.

The local rest/kinetic ledger does not settle the global formation energy budget. Binding-energy changes and radiation/gravity transport still need accounting. In particular, bulk kinetic energy above or below that required by the previous circular snapshot is not by itself proof of a deficit: infall can change kinetic and potential energies. The earlier deposited mass profile must be recalculated if this source rule is adopted.

This advances the formation and momentum requirements of goals4–5. The photon/time source, redshift and arrival-time mechanism, observational comparisons and withheld evaluation remain open. A particle source derived here is a candidate closure, not a replacement claim that gravity has been explained.

## Reproduction

Run `python research_work/results/capture-injection/run.py`, with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 for practical runtime. It reads the preserved retention profiles; it does not alter them or access observations. The absorption, gravity and direct angular methods are explicit in the code. `export.py` produces the comparison table and final checkpoint. Numerical gates, angular comparisons, source fractions and local ledgers are all retained in results.json.

## Final checkpoint

| kappa0 | Initially bound source energy, c/v0=300 | c/v0=1000 | c/v0=3000 | Source energy converted to bulk kinetic energy |
|---|---:|---:|---:|---:|
| 0.1 | 100% | 100% | 71.5385% | 4.88093e-06% |
| 10 | 1.41264% | 0.15194% | 0.0149494% | 0.0395362% |
| 1000 | 0.00308574% | 0% | 0% | 4.38417% |

All declared gates pass. The largest change in bound-energy fraction across radial resolutions is 0.5394 percentage points. The largest absolute beta disagreement between continuity and direct angular integration is 0.0011051; angular refinement changes beta by at most 1.7926e-05. Local rest-plus-kinetic energy and momentum identities pass to floating-point tolerance.

At c/v0=1000, the middle-opacity example puts only about0.152 percent of its instantaneous captured-energy source into initially bound cold cohorts under this rule. In the weakest-opacity example all sampled source regions pass that initial binding check. Both still inject radial flow rather than circular support.

Zero entries retain the finite-grid qualification above. No source luminosity or physical age is inferred. [comparison.csv](comparison.csv) includes rest-mass fractions and the unresolved central-shell source bounds.
