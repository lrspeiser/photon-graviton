# Can an internal oscillator supply small energy transfers?

## Why this calculation matters

Our last comparison kept photons straighter but still transferred about a quarter of their energy per event. Uneven large losses spread sharp colors into broad distributions. This calculation asks whether a spring-like internal material state can favor small energy packets and account for the energy it temporarily stores. Without such a cause, choosing tiny transfers in a formula would not explain why the light behaves that way.

The oscillator is an optional candidate, not an adopted ingredient. Small packets alone are insufficient: different colors must acquire compatible fractional shifts, and the same model must eventually explain directions, timing, capture and gravity.

## 1. Add a dynamical energy receiver

Use the preceding point-target, heavy-material approximation with a canonical scalar companion chi. Introduce an internal coordinate q on the target worldline, in units hbar=c=1:

```
S_q = integral d tau [ (q_dot^2-Omega_0^2 q^2)/2
                      + g q chi(X) - (lambda/2) q E_mu E^mu ]
E_mu = F_mu_nu u^nu,      u.u=1, signature (+---)
```

The bulk scalar and Maxwell actions and dynamical target center of mass are retained from the preceding candidate. In the target rest frame the last term is +(lambda/2)q |electric_field|^2. It replaces the constant chi-dependent polarizability by a material state that responds and then radiates. It neither prescribes changing clocks nor assumes an expanding universe.

Oscillator–scalar-field couplings are a standard comparison framework; see the action in [Iso, Yamamoto and Zhang, On the cancellation mechanism of radiation from the Unruh detector](https://academic.oup.com/ptep/article/2013/6/063B01/1626128). That reference is not evidence for the fictional mechanism. We use a stationary material target here, not its acceleration interpretation. The light coupling and redshift calculation below are specified candidate work.

## 2. Damping has an energy destination

For a classical excited oscillator coupled to outgoing massless scalar waves in three spatial dimensions, the far field is chi approximately g q(t-r)/(4pi r). Integrating its outward energy flux gives:

```
P_companion = g^2 q_dot^2/(4pi)
gamma_chi = g^2/(4pi)
q_ddot + gamma_chi q_dot + Omega_R^2 q = lambda J(t)
J = |electric_field|^2/2
H_q = (q_dot^2+Omega_R^2 q^2)/2
dH_q/dt = lambda J q_dot - P_companion
```

Omega_R is the positive renormalized frequency after treating the near-field frequency shift with a finite-size or ultraviolet regulator. The constant-damping equation is the low-frequency radiative approximation, not proof of arbitrary-frequency point-particle validity. Incoming scalar fluctuations and quantum noise are not represented by this classical driven equation. An unexcited quantum oscillator is not asserted to radiate energy for free.

The runner checks both an initially excited oscillator and an externally driven pulse. In the pulse test, input work is integrated explicitly; stored energy plus radiated energy equals that work. This checks the response subsystem, not a complete self-consistent astronomical photon source.

The q–light vertex also permits two-photon emission from an excited q state. Its rate scales as lambda^2 Omega_R^6 in this heavy point-target approximation, whereas the scalar decay rate scales as g^2. The numerical response below assumes the weak-light-coupling regime where scalar damping dominates. That hierarchy must be checked when selecting physical parameters. Other damping cannot be assigned wholly to companions: ordinary radiation, heat or any other receiver needs its own energy account and branching fraction. No absolute coupling or target abundance has been established.

## 3. A causal frequency response

With Fourier convention exp(-i omega t), the induced scalar/photon coupling is, to leading order in the light coupling:

```
beta_eff(omega) = lambda g / [Omega_R^2-omega^2-i gamma omega]
gamma approximately gamma_chi in the stated regime
```

Its poles give a decaying retarded response when gamma and Omega_R^2 are positive. The numerical tests check the analytic impulse response, energy balance, stable poles and the ideal response's static dispersion sum rule. These are limited response checks, not a proof of the stability or cutoff independence of the full field theory.

The response is largest when the outgoing companion energy omega is near Omega_R. Thus an oscillator frequency much smaller than a photon energy can favor a small transferred fraction. The frequency is a physical scale to explain, not a value inferred or chosen here to claim observational success.

## 4. Predict the transferred-energy distribution

Substitute this response into the preceding **point-target** scalar/polarizability amplitude. Its leading companion-energy distribution is proportional to:

```
d sigma/d omega proportional to E (E-omega)^3 omega |beta_eff(omega)|^2
0 < omega < E
```

This expression includes the same photon derivative factors and three-body phase space as the earlier calculation. At a narrow resonance with E much greater than Omega_R, most events transfer energy near Omega_R. At gamma/Omega_R=0.01, the numerical results are:

| E/Omega_R | Mean fraction transferred per event |
|---:|---:|
| 10 | 0.100028 |
| 100 | 0.0100256 |
| 1,000 | 0.00100308 |
| 1,000,000 | 0.00000100318 |

The response has a tail; the runner calculates its second energy moment as well as the mean. A small mean alone does not establish sufficiently narrow astronomical line profiles after many events. All values are conditional examples, not measured material properties or a new redshift fit.

## 5. Why frequency selection alone does not repair the color dependence

For this point-target interaction, any fixed forward-channel response depending only on the outgoing companion frequency yields:

```
alpha(E) = C integral_0^E (E-omega)^3 W(omega) d omega
W(omega) = omega^2 |beta_eff(omega)|^2 >= 0
C > 0 includes the target density and normalization
```

For a finite differentiable integral with nonzero loss, differentiate with respect to E:

```
d alpha/dE = 3C integral_0^E (E-omega)^2 W(omega) d omega
d ln(alpha)/d ln(E) >= 3
```

The inequality follows because E >= E-omega inside the positive integral. Equivalently, for E_2>E_1, alpha(E_2) >= (E_2/E_1)^3 alpha(E_1). Additional frequency bands with positive weights cannot cancel this dependence. The bound is also numerically checked with three non-Lorentzian positive trial weights.

In our oscillator examples, doubling E from 100 Omega_R to 200 Omega_R increases the initial fractional loss rate by about 8.12–8.15. It approaches the factor eight at larger photon energy relative to the response scale. This is not the desired roughly same fractional loss across colors.

The result is broader than failure of one oscillator parameter choice, but its assumptions matter: the same point-target vertex, positive forward-channel weights, an energy-independent target population, standard dispersion and a response dependent only on transferred frequency. It does not cover a different vertex, an incoming-photon-frequency-dependent response, inverse channels, or a jointly derived extended/multiple-target interaction. It does not exclude the overall companion hypothesis.

## 6. Why the previous forward-scattering result cannot simply be attached

The spatially extended target and this point oscillator each have a response that needs its own derivation. A collection of oscillators radiating the same scalar field can interact through that field. Its response is generally a matrix:

```
D_ret(omega) = [D_0^-1(omega) - Sigma_ret(omega)]^-1
```

Off-diagonal radiation/recoil couplings can alter resonances, damping and which collective states radiate. Multiplying the earlier Gaussian profile by the single-oscillator response without checking these effects would borrow favorable pieces from models that have not been shown compatible. No such combined successful model is claimed here. The field fluctuations and any additional damping recipients also belong in that calculation.

## 7. What we learned and what comes next

An internal state can favor small energy packets while its stored and radiated energy is tracked. In the simple point-target interaction, it cannot also remove the strong color dependence merely by changing its outgoing-frequency response. It also retains the point target's broad photon-angle pattern. Small packets are therefore a partial ingredient, not the completed redshift mechanism.

The next useful calculation must change something the bound actually depends on: derive a joint collective response, an interaction with a different photon-frequency dependence, or the ordinary spin-two emission channel. Testing more scalar frequency filters alone would not resolve the demonstrated bound. No new law is selected by this report. Event timing, supported deposits and joint motion/lensing remain open in the original goal.

See [response, energy and spectrum results](oscillator-response-results.json) and [runner](check_oscillator_response.py). No new astronomical data were fitted. This is preliminary T03/T04/T05 work, not completion of the original 20-task, 32-area program.
