# Pair creation is possible in the conditional model, but the redshift law is unfinished

The [derivation](derivation.md) adds a neutral scalar-pair production channel to the bound-store research branch. It distinguishes conserved charge from total occupation, includes recoil of the companion store and retains the reverse reaction. It does not identify ordinary gravitons or derive the seed store, bound-mode wavefunctions or absolute interaction rate.

## What was calculated

- **180 exact kinematic cases:** a photon leaves a lower-energy photon and increases a companion store's internal energy, with store recoil included. Maximum relative energy residual is 3.56e-16. Forward transfer is allowed kinematically without mandatory photon deflection; an amplitude must still determine its probability.
- **12 closed quantum-mode runs:** initially empty particle/antiparticle modes can be populated from a finite supply of high-energy photons. The finite basis is fixed by invariants of the retained mode Hamiltonian; projection from the full field operator is an approximation. State evolution agrees with independent matrix exponentials within 2.18e-13; maximum absolute energy-ledger error is 9.89e-16. Charge, photon count and momentum are constant in each retained basis state; the full proposed operator also has photon-number-changing channels.
- **Reverse transfer:** for a single initially available photon on resonance, pair probability is sin^2(g t). The pair appears fully at g t=pi/2 and disappears at g t=pi, returning energy to the photon. A permanently occupied store does not follow from this closed model. Escaping radiation requires a spatial calculation.
- **Finite-source/recoil effects:** for twelve initial photons, maximum sampled pair occupation ranges from about 3.20 to 5.23 across the stipulated masses and couplings. These are finite-history expectation values, not efficiencies for an astronomical source or stable deposits.

The parent interaction also permits pair annihilation into two photons, omitted from the selected mode Hamiltonian. Its rate must be calculated from the same coupling before claiming long storage; merely letting the incident light leave does not settle that question. This is an analytic channel check, not a lifetime measured by these runs.

The seed store's initial mass is explicitly included. It is a companion store so that recoil remains within that sector; its unexplained origin cannot be counted as photon-created gravity. Pair rest energy, negative binding energy, recoil and shared interaction energy are separately accounted. Counting all photon loss as new constituent rest mass would be incorrect.

## The selected transition has a color problem

For a fixed bound-pair energy gap delta and initial receiving mass M, a forward event removes Delta=delta+delta^2/(2M). In the illustrative delta=0.2, M=100 case:

| Initial photon energy | Lost energy in one forward event | Conditional redshift |
| ---: | ---: | ---: |
| 0.5 | 0.2002 | 0.667779 |
| 1 | 0.2002 | 0.250313 |
| 2 | 0.2002 | 0.111235 |
| 4 | 0.2002 | 0.052687 |

These are dimensionless diagnostic energies and conditional single-event shifts. They are not a population prediction or an observational fit. Equal fractional loss would require a compensating event rate or a different receiving spectrum.

The point-electric-dipole, heavy-store limit of the proposed operator instead gives a pair-creation rate proportional to E(E-delta)^3, hence fractional loss rate proportional to delta(E-delta)^3. Relative to E=1, fractional loss rates at E=0.5, 1, 2 and 4 are approximately 0.0527, 1, 11.3906 and 107.1719. The absolute coefficient is not known. This limit therefore does not meet the common-redshift requirement. Finite spatial overlap, a continuum of modes, inverse occupation factors and different interactions need their own derivation; this result does not exclude them.

## What remains before this could serve our universe

This step removes a structural obstacle to population production in a **specified optional scalar extension**: conserving charge need not freeze total occupation. It does not derive ordinary-graviton production, a stable growing halo, a weak-gravity coupling, or the proposed time law. In particular, discrete frequency conversion has not been shown to stretch event arrival intervals.

Next relate illuminated production to radiative pair loss, derive spatial mode overlap and the emitted spectrum, then test whether the same interaction yields an evolving propagation factor with common frequency/event stretching, small spectral and angular distortion and an accounted source history. Reuse the earlier vacuum and matter-assisted constraints. Seed formation, irreversible retention, supported halo structure, available energy, matter clocks and joint motion/lensing remain open. The color constraint applies to using these events directly as the redshift law; pair production as a source of a changing coherent background remains to be tested from the common interaction. The full 20-task/32-area objective stays active.

[Protocol](protocol.json), [reproduction script](check.py) and [saved machine-readable results](bound-pair-production-results.json) accompany this report. Numerical checks verify the conditional calculations; they do not turn their hypotheses into established physics.
