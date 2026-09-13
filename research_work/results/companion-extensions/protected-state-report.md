# Protected-state storage with thermal return and an energy account

13 September 2026. Alternative finite-site model, not a completed companion interaction.

## Outcome

A lower-energy state can become highly occupied while coupling only indirectly to the incident light. Thermal return then exposes some stored energy to loss through the illuminated state. The calculation supplies a conditional protection mechanism with a cost: energy released while settling into the lower state cannot also be counted as retained gravity.

This is a different state structure from the preceding unlimited bosonic pair mode. Each modeled site permits one excitation. That finite capacity is postulated, not derived from exclusion statistics, the one-third law or galaxy data. The number and physical nature of sites remain unknown.

## Levels and shared thermal rates

**Hypothesis:** a site has empty state 0, illuminated state a with energy E_a, and protected state b with E_b=E_a-epsilon. Photons excite 0->a and can reverse a->0. The protected state has no direct optical transition in this idealized model. Its indirect path is b->a followed by radiation or photon-assisted removal.

**Known thermal-bath relations**, adopted for a weakly coupled bosonic bath:

\[
n_B=(e^{\epsilon/(k_BT)}-1)^{-1},\quad
k_\downarrow=\kappa(n_B+1),\quad k_\uparrow=\kappa n_B,
\quad \frac{k_\uparrow}{k_\downarrow}=e^{-\epsilon/(k_BT)}.
\]

These rates enforce the known thermal detailed-balance ratio, rather than freely choosing a one-way sink. The bath and its coupling to companions are hypotheses; statistical mechanics itself does not establish their existence. See [Tong's quantum-gas treatment](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S3.html) and [kinetic-theory discussion of detailed balance](https://www.damtp.cam.ac.uk/user/tong/kintheory/).

With optical excitation A, inverse rate B and spontaneous a->0 rate gamma_a, the population equations are

\[
\dot P_0=-AP_0+(B+\gamma_a)P_a,
\]
\[
\dot P_a=AP_0-(B+\gamma_a+k_\downarrow)P_a+k_\uparrow P_b,
\qquad \dot P_b=k_\downarrow P_a-k_\uparrow P_b.
\]

These are known rate-equation mathematics applied to the new level proposal. A direct b->0 decay is set to zero only as an optimistic diagnostic. Its absence requires a physical selection rule or derived small matrix element; it cannot be asserted for actual deposited companions. Other heating and collisional channels are also absent.

## Energy tracking

Normalize E_a=1. In general units, net photon input has rate E_a(AP_0-BP_a), net bath energy gain is epsilon(k_down*P_a-k_up*P_b), and escaping a-state radiation has rate E_a*gamma_a*P_a. Starting empty, these integrated accounts obey

\[
E_{\rm photon,net}=E_aP_a+E_bP_b+Q_{\rm bath}+E_{\rm emitted}.
\]

The maintained thermal bath can absorb or supply energy; its signed account is retained during illumination and darkness. A fixed temperature is a reservoir assumption, not an untracked free cooling resource. A real finite bath must evolve when this exchange changes its energy appreciably.

At illumination equilibrium P_b/P_a=exp(epsilon/k_BT). In darkness set A=B=0 but keep spontaneous bright-state loss and both thermal rates. The late-time decay exponent of the populated subspace is

\[
\lambda_{\rm slow}=\frac{S-\sqrt{S^2-4\gamma_a k_\uparrow}}2,
\qquad S=\gamma_a+k_\downarrow+k_\uparrow.
\]

The code uses its rationalized form to avoid cancellation. For k_down much greater than gamma_a and k_up much smaller than k_down, lambda_slow is approximately gamma_a*exp(-epsilon/k_BT). This is a thermal-protection limit, not a measured lifetime. ln(2)/lambda_slow is the half-time of the slow mode, not necessarily the half-time of an arbitrary initial state during its initial transient.

## Tested examples

Set A=1, B=1.01^3 from the preceding ideal reciprocal coefficient at a small gap/photon ratio, and gamma_a=0.001. These are illustrative rate ratios; they do not derive the optical and spontaneous rates from one coupling. Test epsilon/E_a=0.1,0.5,0.9, k_BT/E_a=0.01,0.1,1 and kappa/A=0.1,1,10. Illumination runs for 20/A and then darkness for 1000/A.

The following rows fix k_BT/E_a=0.1 and kappa/A=1. All energies are per initial empty site in units E_a, and all times are in units 1/A.

| Released fraction epsilon/E_a | Protected occupation after illumination | Stored energy after illumination | Net bath energy gain | Stored energy after dark interval | Slow-mode half-time |
|---|---:|---:|---:|---:|---:|
| 0.1 | 0.5723 | 0.7256 | 0.0572 | 0.5546 | 2,578 |
| 0.5 | 0.9860 | 0.4998 | 0.4930 | 0.4964 | 103,667 |
| 0.9 | 0.9991 | 0.1003 | 0.8992 | 0.1000 | 5,622,938 |

The remaining energy is accounted for in emitted radiation and the net optical input. Deeper relaxation gives stronger isolation in this comparison but retains less energy per occupied site. A lower bath temperature is another way to suppress return, with its own physical requirement. None of these values selects a galaxy formation age.

The small creation gaps considered in the recent redshift diagnostics make the assumed bath temperature significant. For k_BT/E_a=0.1, E_a=1e-8 eV corresponds to T approximately 1.16e-5 K; E_a=2e-10 eV corresponds to approximately 2.32e-7 K. Those are illustrative conversions, not measured companion temperatures. A dark bath need not share the electromagnetic radiation temperature, but its decoupling, preparation and heating must be explained. The fitted 17.8 eV phase mass is not silently substituted for this excitation gap.

## Verification and limits

All 27 cases use matrix exponentiation with integrated energy counters. Independent ODE integrations agree within 1e-8 in normalized variables. Probability remains normalized, and illumination plus dark energy ledgers close within 1e-9 E_a. Analytic stationary populations satisfy the rate equations. No time in years, absolute reservoir capacity, mode density, force law, lensing profile or observational fit is inferred.

These calculations show that a lower protected level with an explicit thermal reverse path is mathematically consistent. They do not establish the level spectrum, eliminate all decay channels or repair the earlier ideal-condensate support failure. Under sustained illumination a finite site eventually saturates; this alone cannot produce indefinite accumulated gravity. New available sites, transport, or an extended many-body state would require additional derivation.

## Decision

Retain protected-state relaxation as a candidate connecting absorption to storage, with its released energy and bath explicitly counted. The next physical requirement is a Hamiltonian or selection rule that produces weak direct coupling of b while allowing formation, and a bath whose energy budget supports the assumed temperature. The existing galaxy settling and saturation formulas cannot be identified with these microscopic rates without that connection. Redshift/event timing, photon supply, universal retention and joint lensing/motions remain open.

Reproduce with `python research_work/results/companion-extensions/protected-state.py`. [Source](protected-state.py), [complete results](protected-state-results.json).
