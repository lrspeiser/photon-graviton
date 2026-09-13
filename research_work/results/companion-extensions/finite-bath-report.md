# Protected storage with a finite, warming bath

13 September 2026. Energy-accounted extension of the protected-site candidate.

## Result

The fixed-temperature assumption is consequential. A small bath warms as it receives relaxation energy, increasing the rate at which protected excitations return to the illuminated state. A larger bath or an explicit cooling outlet can preserve protection in the tested examples. Neither is a free component: its capacity, energy and physical realization must be supplied.

This continues the same hypothetical three-level, single-excitation site. It does not derive a galactic reservoir, a bath particle spectrum, an absolute lifetime, or a temperature observed in our universe.

## Added energy equation and provenance

Keep E_a=1, E_b=1-epsilon, and define theta=k_BT/E_a. **Hypotheses:** bath energy changes linearly with theta, U_bath/E_a=C(theta-theta_0), and a linear outlet exchanges heat with an environment at theta_0. C is dimensionless, equal to physical heat capacity per site divided by k_B. It is held constant only for this diagnostic, not claimed to be the actual low-temperature heat capacity of a bosonic bath.

**Known energy accounting applied to those assumptions:**

\[
C\dot\theta=\epsilon(k_\downarrow P_a-k_\uparrow P_b)
-C\ell(\theta-\theta_0).
\]

The prior thermal rates now use this evolving temperature, including k_up/k_down=exp(-epsilon/theta). Integrate the signed outlet power C*ell*(theta-theta_0) as Q_out; negative values would mean energy entering from the ambient reservoir, not disappearing cooling costs. These thermal relations are standard statistical mechanics; the companion bath and outlet are the proposed physics. See [Tong's thermodynamics notes](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S4.html).

Starting empty at theta_0, the complete normalized ledger is

\[
E_{\rm photon,net}=P_a+(1-\epsilon)P_b
+C(\theta-\theta_0)+E_{\rm emitted}+Q_{\rm out}.
\]

The bath energy also obeys the independent identity C(theta-theta_0)+Q_out=epsilon*P_b. A protected excitation is counted at its lower stored energy, and its released energy is counted in the bath or outlet, never in both places.

## Closed-bath consequence

With no outlet, theta=theta_0+epsilon*P_b/C. The former arbitrarily cold fixed bath is therefore replaced by a temperature linked to occupancy. Under steady illumination with A=1, the stationary probability solves

\[
P_b=\{1+(1+B+\gamma_a)
\exp[-\epsilon/(\theta_0+\epsilon P_b/C)]\}^{-1}.
\]

This is a consequence of the stipulated finite bath and rate equations, not a new fundamental law. The right-hand side decreases as occupancy heats the bath, providing a finite self-consistent solution. As C tends toward zero at fixed positive epsilon, its hot-bath limit approaches 1/(2+B+gamma_a), approximately one-third for the illustrative B near one. That is a three-state occupancy limit, unrelated to the retention-law exponent.

## Numerical examples

Use initial theta_0=0.01, illumination A=1, reverse coefficient B=1.01^3, bright-state spontaneous decay gamma_a=0.001, and relaxation scale kappa=1. As before, direct protected-state decay is omitted as an optimistic assumption. Illumination lasts 20/A, followed by 1000/A in darkness. Temperature and rates continue evolving after switch-off.

Thirty-six cases span epsilon=0.1,0.5,0.9; C=0.1,1,10,100; and outlet rate ell=0,0.01,1. The following rows show epsilon=0.5 without an outlet. Energy is normalized to E_a; temperature is theta, and times have no assigned astronomical units.

| C | Protected probability after illumination | Bath temperature | Thermal up/down rate ratio | Stored energy after dark interval |
|---|---:|---:|---:|---:|
| 0.1 | 0.3887 | 1.9536 | 0.7742 | 0.3181 |
| 1 | 0.6749 | 0.3474 | 0.2371 | 0.4066 |
| 10 | 0.9989 | 0.05995 | 2.39e-4 | 0.4997 |
| 100 | 0.9994 | 0.01500 | 3.32e-15 | 0.4998 |

These finite-time results need not equal the stationary solution. A high protected probability is not equivalent to retaining all incident energy: approximately half the excitation energy in these rows goes into the bath during relaxation. Increasing heat capacity preserves protection by limiting warming; adding the outlet exports explicitly recorded energy instead.

## A capacity requirement instead of a free cold bath

For a closed bath to have occupancy P_b at thermal backflow ratio no greater than r_T, the temperature must obey theta<=epsilon/ln(1/r_T). Thus a necessary capacity condition is

\[
C\geq\frac{\epsilon P_b}{\epsilon/\ln(1/r_T)-\theta_0},
\]

provided the denominator is positive. For epsilon=0.5, P_b=0.9, r_T=0.001 and theta_0=0.01, this requires approximately C>=7.21. If the initial bath is already too warm, increasing its capacity cannot meet that temperature target without additional cooling. This is an inverse design condition, not a measurement or a proof that the target occupation is reached dynamically.

## Verification and implications

All cases conserve probability, the full energy ledger and the independent bath/occupancy identity to the stated numerical tolerances. Six representative cases were rerun with tighter ODE tolerances; every recorded normalized quantity changed by less than 1e-6. The code also records stationary roots for all 12 closed-bath combinations. No finite-time result is promoted to a stability proof or a prediction in years.

The candidate remains possible at the level of rate equations, but a physical model must supply the bath's density of states, temperature-dependent heat capacity, volume and transport. A very cold bath can have a much smaller heat capacity than the constant-C assumption; the resulting mass/volume and cooling costs must be compared with the gravitational energy inventory. The outlet must specify where its energy goes and whether it is observable radiation or escaping companions.

The next decisive calculation is that physical bath budget. Choosing a large C without constructing its degrees of freedom would only hide the storage-energy problem. Photon supply, mode selection rules, direct protected-state decay, event timing, galaxy support and joint lensing/motions remain unresolved. The exact-third reference remains unchanged.

Reproduce with `python research_work/results/companion-extensions/finite-bath.py`. [Source](finite-bath.py), [all outcomes](finite-bath-results.json).
