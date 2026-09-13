# Capture with recoil: allowed kinematics, unresolved interaction

13 September 2026. Bounded continuation of the orbital-support analysis. No new
observations, fitting or change to exact-third retention.

## Outcome in plain language

An incoming companion does not necessarily have to throw away most of its energy
to stop. A sufficiently massive local receiver can take its momentum with very
little recoil energy. This can either increase the receiver's internal energy or,
if an appropriate interaction exists, produce a separate stored particle while
the receiver recoils. Conservation alone does not rule out either branch.

This extends the earlier twelve-case capture ledger by comparing recoil budgets
with all ten fitted reservoir support profiles and distinguishing an independent
new particle from simple absorption into an existing object. It does not supply
the interaction, its probability, the new particle's mass spectrum or the required
tangential orbits. A receiver's rest energy is not consumed in the ideal recoil
process; the requirement below is an effective locally participating mass, not a
new energy source. The whole galaxy cannot be assumed to respond instantaneously.

## Derivation and provenance

Use a local inertial frame, initially stationary receiver mass M, incident null
packet energy E, and q=E/(M c squared). The conserved incoming energy and momentum
are M c squared + E and E/c. Standard special-relativistic invariant kinematics
applies; it is not a novel companion law. Reference: [Particle Data Group,
Kinematics, sections 49.1-49.2](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf).

**A. Absorption into a composite receiver.** No outgoing radiation is required
kinematically if the receiver can access the resulting internal state:

    M_final = M sqrt(1+2q)
    Delta_rest_energy/E = 2/[sqrt(1+2q)+1]
    K_recoil/E = q/[1+q+sqrt(1+2q)]
    v_final/c = q/(1+q)

The two fractions sum to one. The stored energy is an increase in the composite's
invariant mass, not automatically a new independently orbiting particle. A
discrete absorber needs an allowed transition or other accessible internal state.

**B. Separate particle at rest, unchanged receiver invariant mass.** Let the
receiver carry all incoming momentum while a newly produced particle has zero
momentum in this frame. Then

    K_recoil/E = q/[sqrt(1+q squared)+1]
    m_new c squared/E = 1 - K_recoil/E

The outgoing receiver energy is sqrt(M squared c fourth + E squared). Adding the
new particle energy gives exactly M c squared + E; momentum remains E/c. These
equations determine the required new mass for this special final state. A fixed
particle mass cannot be freely retuned for every incident energy: an actual
reaction must satisfy its mass thresholds, momenta and other conserved charges.
The particle produced at rest does not yet have the tangential velocity required
by our circular support construction. Falling radially is not circularization.

**C. Two incident null packets without a receiver.** Define
delta=abs(E1-E2)/(E1+E2), with mutual direction angle theta. Their composite speed
before any emission obeys

    beta_CM squared = delta squared
        + (1-delta squared)(1+cos(theta))/2

An equal-energy opposite-direction pair has zero momentum and nonzero invariant
mass. For independent isotropic directions, the geometric fraction with
beta_CM<=b is max[0,(b squared-delta squared)/(1-delta squared)]. At equal energy
this is b squared. Actual collisions weight angles by relative flux and the
unknown interaction cross section. A constant cross section with the usual
massless relative-flux weight 1-cos(theta) instead gives 2b squared-b fourth for
equal energies. Neither number is the actual capture efficiency. A pair with
delta>b cannot meet this speed ceiling without another momentum/energy channel.

All three derivations are applications of established energy-momentum algebra.
Their use as companion capture mechanisms is hypothetical. The shared one-third
law does not follow from these equations.

## Quantitative comparison with the fitted support budgets

For a distinguishing diagnostic, require recoil/input energy <= epsilon, where
epsilon is the refined profile's previously calculated total orbital K/(M_d c
squared). This is a chosen comparison budget, not a measured cooling limit or
necessary condition for capture. Inverting the equations gives

    A: M c squared/E >= (1-epsilon) squared/(2 epsilon)
    B: M c squared/E >= (1-epsilon squared)/(2 epsilon)

These are almost identical for the small epsilon values here. For branch B:

| Profile | Required receiver rest energy / packet energy | Equal-energy geometric slow-pair fraction |
|---|---:|---:|
| J0037 free | 657000 | 0.00000205 |
| J1112 free | 385000 | 0.00000360 |
| J1204 free | 465000 | 0.00000286 |
| J1402 free | 311000 | 0.00000450 |
| J1621 free | 7170 | 0.000201 |
| J1630 free | 425000 | 0.00000288 |
| J1621 scale 0.3 | 677000 | 0.00000171 |
| J1621 scale 10 | 165000 | 0.00000872 |
| J1621 scale 30 | 37100 | 0.0000391 |
| J1621 scale 100 | 7170 | 0.000201 |

For illustration only, a 1 eV packet and ratio 657000 require participating
receiver rest energy at least about 0.657 MeV under the chosen budget. The packet
energy is not inferred from our data. Grouping many packets into one coherent
event increases E and can change this estimate substantially. The pair column
uses each profile's maximum circular speed as a permissive speed ceiling, not
local escape speed or an orbit-assignment criterion. It therefore does not prove
binding anywhere. Fractions are especially conditional for hollow envelopes.

## Consequences for the theory and next discriminant

Receiver-assisted capture is kinematically more permissive than relying only on
fortuitously opposed, closely energy-matched pairs. We cannot rank their rates
without interaction strengths and local packet/receiver distributions. Efficient
ordinary-matter absorption would attach added energy to ordinary matter and does
not by itself produce our independently specified extended density. Separate
stored particles or collective field states need an explicit production and
redistribution law; a pre-existing companion receiver also raises a seeding issue.

Recoil can remain kinetic rather than becoming radiation immediately. Consequently
these calculations do not predict unavoidable heat luminosity, a cooling spectrum
or a black-hole emission signal. Dissipation, settling and repeated directional
kicks require their own evolution and energy accounting. Causal receiver size and
angular-momentum transfer are the next distinguishing requirements. A meaningful
next experiment would specify one local receiver interaction and calculate its
resulting bound-orbit distribution, rather than assigning circular orbits by hand.

## Verification and reproduction

Run `python research_work/results/companion-extensions/capture-recoil-threshold.py`.
It records the refined support-input hash, ten profile comparisons and eight
incident/receiver energy ratios from 1e-12 to 1000. Checks cover energy partition,
the composite invariant, analytic budget inverses and equal/unequal pair limits.
These verify the conditional algebra, not the interaction or physical formation.
The previous goal turn was progress: the repository backup was committed and
release-asset hashes verified. This turn advances the physical capture assessment;
the comprehensive research goal remains unfinished.
