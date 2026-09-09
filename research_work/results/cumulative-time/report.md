# Cumulative time stretching with companion-energy transfer

## The model requirements used here

The user proposes a new property of time/propagation in the alternate universe: light accumulates stretching during its journey through affected space, particularly voids. This is not a proposal to obtain the effect from ordinary static gravitational time dilation. At conversion, all photon energy lost through the stretching goes into the companion/graviton sector. The companions can travel beside the light or on their own paths. These are agreed model requirements, not questions about whether a receiver should be chosen.

This calculation supplies a conditional mathematical bridge from that proposal to observable predictions. It does not claim to derive a fundamental time field or to identify companions as ordinary spin-2 gravitons. The new closure tested here is that **the entire signal** stretches: its wave cycles and its event features share the same stretch factor. Photon-energy loss alone did not imply that closure in the earlier conversion-only calculation.

## 1. Derivation from the declared postulate

Let S(s) be cumulative stretch after path length s. A small additional interval ds adds fractional stretch alpha(s) ds. Successive intervals multiply, so

    S(s+ds) = S(s)[1+alpha(s)ds] + terms of order ds^2,
    d ln S / ds = alpha(s),
    S(R) = exp[integral_0^R alpha(s) ds].

Here alpha is an effective stretch per unit path length. Its environmental dependence is not specified by naming time or voids. If stretching is instead parameterized per unit travel time by kappa, then alpha=kappa*dt/ds; the propagation rule must supply dt/ds. A constant fractional rate is an initial benchmark, not a first-principles result.

For an event over which S is constant, write the signal's arrival map in calibrated source/receiver time units as

    t_out = T + S t_in.

T fixes the arrival origin; differences obey Delta t_out=S Delta t_in. Phase transported by this map obeys phi_out(t_out)=phi_in((t_out-T)/S). Differentiating phase gives nu_out=nu_in/S. Thus the same map produces both redshift and event-duration stretching:

    1+z = S;   duration ratio = S;   photon energy ratio = 1/S.

The last result uses the same measured Planck relation at emission and reception. Counting a clock-unit change as transferable energy would not implement the user's required conversion; actual photon energy measured with the stipulated standards must decrease.

Two segments compose as

    S_total=S_2*S_1;   T_total=T_2+S_2*T_1.

The delay must transform as well as the duration. A physical implementation needs one consistent phase/time evolution; independently resetting the time origin for each flash is not a derived mechanism. A time-independent linear delay kernel cannot implement S other than one for arbitrary signals, because it preserves an imposed translation of the source signal. A new time law, evolving state or other nonstationary/nonlinear mechanism remains available; ordinary gravitational clock comparisons are not being substituted for that law.

## 2. A constructive propagation example through a complete void

To test that the proposed behavior can survive exit from a void, use an optional effective propagation-clock factor n(x,t), in reference c=1 units:

    omega(x,t,k) = k/n(x,t),
    n(x,t) = 1 + alpha*g(x)*t,
    g(x) = sin(pi*x/L)^2,  0<=x<=L.

Both endpoints have n=1 at all times. The clock standards used to compare emission and reception are fixed and equal. This construction is not ordinary static time dilation; it is an explicitly evolving photon propagation law. It changes effective photon propagation speed relative to the reference clock, so it is an optional realization rather than a claim that the user has chosen a varying-speed theory. A complete theory must establish its relation to matter clocks and measured distances.

Hamilton's equations give

    dx/dt=1/n;   dk/dt=k*(partial_x n)/n^2.

The code integrates these equations with x as independent variable. It does not impose a separate photon-frequency-loss equation. Differentiating the Hamiltonian along the ray yields

    d ln omega/dx = -partial_t n = -alpha*g(x).

The arrival equation is dt/dx=n. Varying launch time gives J=partial t(x)/partial t_launch with

    dJ/dx=alpha*g(x)*J;   J(0)=1.

Consequently omega_out/omega_in=1/J and J(L)=exp(alpha*L/2). The redshift and arrival-time stretching survive passage back to n=1. The spatial gradient also changes wave momentum; a full companion/matter action must account for that momentum, not only energy.

Fifty-four cases vary slab length, rate, launch epoch and initial frequency. Two separately launched flashes are followed in each case. Frequency change, their finite arrival separation and the required companion-energy gain agree with the derived relations to better than 3e-12 relative error. Zero rate gives no effect. These checks establish consistency of this prescribed-ray example. They do not show that n follows from self-consistent field equations, that its backreaction can be neglected, or that real astronomical propagation has this law.

The earlier [matched-wave action](../matched-wave/derivation.md) provides an electromagnetic route to omega=k/n and energy exchange with a dynamical field. It also shows why prescribing n without backreaction is incomplete and why matter-clock coupling matters. Identifying that receiving field with the companion sector and deriving traveling excitations, the void profile and capture remains work to do. The agreed energy receiver itself is not being reconsidered.

## 3. The companion energy and brightness predictions

For an initial photon packet energy E_0 and conserved photon count,

    E_gamma=E_0/S,
    Delta E_companion=E_0(1-1/S).

The receiving term is fixed by the user's rule. The ray calculation integrates dE_companion/dx=alpha*g(x)*omega, and checks that it equals the photon loss. This is a required receiving ledger attached to the photon Hamiltonian, not yet a closed Hamiltonian for both sectors.

With unchanged geometry, no extra attenuation/focusing, and the same clocks and detector calibration, photon arrivals per unit time decrease by 1/S. Received bolometric power therefore decreases by 1/S^2. Energy integrated over the complete event decreases only by 1/S, because the event lasts S times longer.

| Redshift z | Duration factor | Photon energy remaining | Energy transferred to companions | Received power relative to unshifted event |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 1.1 | 90.91% | 9.09% | 82.64% |
| 1 | 2 | 50% | 50% | 25% |
| 2 | 3 | 33.33% | 66.67% | 11.11% |
| 5 | 6 | 16.67% | 83.33% | 2.78% |

For source spectral luminosity L_nu and Euclidean distance R, this closure predicts

    F_nu_observed(nu) = L_nu_emitted(S*nu)/(4*pi*R^2*S).

Integrating over observed frequency gives the bolometric S^-2 result. At corresponding event phases the emitted time argument is (t_out-T)/S. This is an isolated-source relation, not a derivation of a thermal cosmic background distribution. Deterministic affine stretching preserves relative feature spacing and fractional line width; stochastic broadening or chromatic effects require an additional physical model.

The previous conversion-only branch predicted power suppression S^-1 with no duration stretch. Its historical results remain unchanged. The whole-signal branch predicts an additional arrival-rate suppression. It has **the same total companion-energy supply** at a given photon redshift, so solving the timing relation does not itself solve an energy shortfall for gravity.

## 4. Capture and distance illustrations

As an optional bookkeeping limit, use a common evolution parameter q, no companion travel loss and permanent deposits:

    dE_gamma/dq=-a E_gamma,
    dE_companion/dq=a E_gamma-k E_companion,
    dE_deposit/dq=k E_companion.

The total is constant. Twelve ODE/analytic comparisons verify positivity and conservation, including equal conversion/capture rates and zero capture. These equations do not require companions to stay on the photon path: the common parameter is a bookkeeping coordinate, not a derived spatial transport map. Own-path motion, focusing, capture probability and supported storage must be calculated separately. No-loss travel and permanent deposits are tested limiting cases, not imposed universal rules.

Using the earlier fitted effective rate solely as an illustration gives S approximately 1.0793, 1.4646 and 2.1450 at 1, 5 and 10 billion light-years of uniformly affected path. These are extrapolations far beyond the reused 10-93 Mpc sample. That fit did not identify void-specific alpha or a fundamental time mechanism. Actual void weighting would require an environment map and an alpha(environment) law.

## 5. What is established and next work

Starting from whole-signal cumulative stretching and the agreed energy-transfer rule, one can write mutually consistent redshift, timing, brightness and companion-energy predictions. A prescribed propagation-clock example realizes the frequency/arrival relation through a complete void with matching endpoint standards. Gaussian pulse integrations independently verify photon-count, energy and duration factors.

This is a derivation **from postulates**, not a derivation of the postulates from established first principles. Next derive or test a companion/time-field coupling that generates the cumulative law with backreaction and momentum conservation, then its propagation and capture. The timing estimator's unresolved numerical calibration and all 32 observational requirements remain in scope. No observation is newly claimed to agree.

```sh
python research_work/results/cumulative-time/check.py
```

Requires NumPy and SciPy. The standard offline diagnostic runner includes this check; mathematical checks and snapshot integrity are not scientific validation.
