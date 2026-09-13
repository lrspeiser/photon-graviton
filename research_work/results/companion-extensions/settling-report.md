# Conservative capture and shared settling

13 September 2026. Follow-up to the phase-selector pilot. The original one-third capture law and input inventory remain unchanged.

**A single shared parameter set improves the inner Milky Way fit for both ordinary-matter baselines. Its outer predictions are unchanged. Settling can satisfy a conditional energy ledger by releasing energy into an outgoing channel, but neither the capture interaction nor the settling rate has been derived.**

## What changed in the calculation

Previously each matter baseline chose its own phase mass and contraction. We now choose one pair by minimizing the mean of their inner-bin mean-square errors. Both descriptions refer to the same Galaxy; they are alternative assumptions, not independent observations. The outer 18 bins remain excluded from selection, but were previously examined.

The shared result is m=17.78 eV/c² and s=0.8032. That means the selected mobile material moves inward by about 19.7% of its original radius. The ideal-phase prescription selects 13.9% or 15.7% of the inventory under the two baselines. These are effective model parameters, not measured particle properties.

| Matter baseline | Inner RMS before → after | Outer RMS before → after |
|---|---:|---:|
| I | 3.68 → **2.74 km/s** | 8.36 → **8.36 km/s** |
| II | 11.83 → **6.70 km/s** | 9.78 → **9.78 km/s** |

The same spherical-gas approximation is used on both sides. The result improves on the matched controls but is weaker than the separately fitted pilot. It does not demonstrate transfer to another galaxy. Because the redistributed material remains within the outer radii, Newtonian spherical gravity there is unchanged.

## Capture: a wave can transfer energy into a massive reservoir

Start locally with a stationary bound reservoir of mass M0. An incoming packet or group of packets carries total energy E and vector momentum P. If everything is absorbed, special-relativistic conservation gives

\[
M_f^2c^4=(M_0c^2+E)^2-c^2|\mathbf P|^2,
\qquad
E=(M_f-M_0)c^2+K_{\rm recoil}.
\]

**Provenance: known energy–momentum conservation, not a new formula.** Our proposed application is absorption of traveling companions by a bound reservoir. This establishes allowed bookkeeping, not the probability that absorption occurs or the constitution of the reservoir.

For a single light-speed packet, |P|=E/c. Some incident energy becomes rest energy and some becomes recoil. For E much smaller than M0c², recoil is approximately E²/(2M0c²), so most energy can enter the reservoir. If equal packets arrive from opposite directions, their momenta cancel and all their combined energy can enter rest energy in an ideal fully absorbing event. Momentum balance makes your idea of arrivals from many directions physically relevant, although an isotropic average is not exact cancellation for every capture.

A freely traveling, exactly collinear group of massless waves has zero invariant mass on its own. It cannot simply stop and become a massive stationary object without another participant, a change of directions or outgoing momentum. A pre-existing reservoir provides that participant. This preserves a light-speed traveling channel while allowing a different bound state; it does not require assigning rest mass to the traveling wave.

The [capture ledger script](capture-ledger.py) evaluates directional imbalance c|P|/E=0, 0.5 and 1 at four incident-to-reservoir energy ratios. Stored rest energy plus recoil equals input energy to numerical precision. Local flat-frame kinematics are used; a global curved-spacetime capture history is not calculated.

Crucially, the fitted 17.78 eV phase mass is still unexplained. Incrementally energizing a large reservoir is not the same thing as proving it creates bosons of that mass or meets an ideal-gas condensation criterion. Energy accounting removes a bookkeeping obstacle; it does not derive the phase model.

## Settling: account for the energy released by deeper binding

We calculated companion self-energy, energy in the fixed ordinary-matter potential and the kinetic energy required by a global virial condition:

\[
U_d=\frac12\int\rho_d\Phi_d\,dV,\qquad
U_b=\int\rho_d\Phi_b\,dV,
\]
\[
2K=\int\rho_d\,\mathbf r\!\cdot\nabla(\Phi_b+\Phi_d)\,dV,
\qquad E_{\rm mech}=K+U_b+U_d.
\]

**Provenance: standard Newtonian potential energy and scalar virial theorem.** Applying them to our proposed spherical deposit distributions is conditional. The external-potential energy has no one-half factor; the companion self-energy does. In a general external potential, the force virial is not simply the potential energy. See the [stellar-dynamics derivation](https://galaxiesbook.org/chapters/I-04.-Equilibria-of-Collisionless-Stellar-Systems_2-The-virial-theorem.html).

The energy that must leave during fully cooled settling is

\[
Q=E_{\rm mech,initial}-E_{\rm mech,final},\qquad
E_{\rm mech,initial}=E_{\rm mech,final}+Q.
\]

**Provenance: conservation identity; choosing complete removal of Q is our proposed cooling closure.** This energy is released by deeper binding. It is not newly supplied photon energy and must not be counted a second time as an increase of deposited rest mass.

| Baseline | Released energy Q | Q / deposited rest energy |
|---|---:|---:|
| I | 4.527×10^50 J | 1.210×10^-8 |
| II | 4.522×10^50 J | 1.225×10^-8 |

These totals refer to the entire numerically represented deposited profile and its prior fitted inventory, not an independently verified radiation supply. No duration is assumed, so these are energies, not luminosities. Assigning a duration tau would give average power Q/tau, but this study does not choose tau.

All 1,353 mass/contraction choices passed the nonnegative-release condition under both baselines, including zero-change limits. Thus the energy-sign requirement does **not** select the successful parameters or explain the contraction amount. It is a consistency gate, not a derived force or rate law.

## How it fits the original concept

The proposed sequence is: traveling companion energy is absorbed with recoil accounted for; a fraction of the bound reservoir can change state and redistribute; deeper binding releases energy into outgoing companions or radiation. The remaining bound inventory stays deposited. Release of binding energy is different from spontaneous decay of the original deposit.

For this branch, the energy channel carrying Q must actually leave the bound subsystem. It can remain in the universe as traveling energy; globally it is not discarded. If it is reabsorbed locally, it heats or energizes the reservoir, and the fully cooled solution must be recalculated. Ordinary photons and outgoing companions are alternative channels with different observational consequences; neither is demonstrated here.

At Newtonian order, constituent mass is conserved during redistribution while mechanical energy changes and escapes. At relativistic order the total system mass-energy changes by Q/c²; one cannot hold every gravitational mass term exactly fixed and also claim a complete relativistic energy solution. The small ~10^-8 ratio makes this a small correction to the present Newtonian calculation, not permission to omit it from a future relativistic theory.

The scalar virial condition only assigns the necessary *total* support energy. It does not show that each radius has the right pressure or orbit distribution, that the phase temperature is self-consistent, or that the configuration is stable. The ordinary-matter potential is held fixed, so recoil/backreaction of the entire galaxy is outside this calculation. These limitations are why the result is a conditional conservative settling model, not a first-principles derivation of the full phenomenon.

## Verification and next decision

[Protocol](settling-protocol.md) preceded the sweep. [Code and results](settling.py) retain all shared-grid cases and individual predictions in `settling-results.json`. Independent checks cover the self-energy of a thin shell, zero-contraction release, mass conservation inherited from the phase routine, and closure of the mechanical-energy ledger. Doubling shell resolution and increasing angular quadrature changes selected speeds by under 0.0022 km/s and released energies by under 0.0061%. Runs complete with warnings treated as errors.

We now have a possible conservation-consistent route from moving energy to a bound reservoir, and a shared placement rule that retains an inner-fit gain. What is missing is a local law that determines the capture/cooling rate and final radius instead of fitting s. That is the next physical task: specify an outgoing-energy channel and solve pressure/orbital support with the updated density and temperature. Without this, we cannot predict how long settling takes or whether it stops at the fitted radius. No new lensing, redshift-timing, solar-system or other-galaxy success is claimed.
