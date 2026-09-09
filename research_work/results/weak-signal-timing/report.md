# Weak signals in an evolving companion/time field

## Why this test follows the backreaction result

The previous interacting field calculation showed that a photon packet can alter the field as it transfers energy and momentum. Changing that packet's energy changes the resulting propagation. This does not establish that arbitrarily faint signals need source-specific stretching rates: they may travel through a field primarily set by other radiation and its existing state.

Here an infinitesimal test signal probes the same dynamical field without altering it. Its trajectory and frequency are calculated independently of the field-generating packet. The companion remains the agreed energy receiver. This test concerns the weak-signal limit of that interaction, not a new receiver or ordinary gravitational time dilation.

## 1. The carrier and timing relation

For the effective smoothed field n(t,x), the test ray has Hamiltonian omega=k/n. Its equations are

    dx/dt=1/n;   dk/dt=k*n_x/n^2.

The arrival trajectory can also be written dt/dx=n. Differentiating with respect to the launch time gives the local arrival Jacobian J:

    J=partial t_received/partial t_emitted,
    dJ/dx=n_t*J;   J_at_source=1.

Along the Hamiltonian ray,

    d ln(omega)/dx=-n_t,
    omega_emitted/omega_received=J.

Thus frequency shift and the stretching of sufficiently close event features agree even in a spatially varying, evolving background. This result follows from the candidate's nondispersive propagation law, not from an imposed photon-energy-loss formula. It does not require a separate parameter for each weak signal's brightness. The test-limit photon amplitude is absent; finite signals still require the full feedback equations.

The two sides are integrated separately: a time-parameterized Hamiltonian ray finds reception by an event crossing, and a distance-parameterized variational equation finds J and arrival time. Initialization k=omega_emitted*n_emitted makes their frequency comparison use identical reference-clock standards.

## 2. Calculation and numerical checks

The background retains the previous closed spectral Hamiltonian, with 32 modes, box length 8, field inertia 1, wave speed 0.5, packet smoothing 0.2 and initial rolling rate 0.05. The single packet driving it has initial energy 0, 0.001, 0.01 or 0.1. All units are dimensionless. A single driving packet is not a realistic isotropic ambient radiation field.

Background evolution is calculated through time 6. Its packet energy at time 3 is cross-checked against the archived backreaction calculation. Test signals are launched at five times from 0 through 1 and propagate reference distances 1 and 2. Forty frequency/timing comparisons pass; eight additional signals at three times the carrier frequency verify achromatic behavior of this selected dispersion law.

Maximum disagreement between the frequency and timing factors is below 2e-14 relative, and independently calculated arrival times differ by less than 3e-15. With no driving radiation, n=1+0.05*t and J=exp(0.05*distance) exactly; the numerical control reproduces this. This is an initial rolling-field solution, not a derived universal history.

The weak signals are not added at finite energy to the background budget. Conservation with finite radiation and field feedback belongs to the preceding Hamiltonian test; claiming it again for finite probes that have been omitted from the field equations would be incorrect.

## 3. Local stretching need not preserve a whole long event

If J changes over the emission interval, frequency and local feature spacing still agree at each point, but the entire light curve is not one uniformly stretched copy. For path length 2:

| Energy of the packet driving the field | Smallest local stretch during the event | Largest local stretch | Fractional range |
| --- | ---: | ---: | ---: |
| 0 | 1.10517 | 1.10517 | effectively zero |
| 0.001 | 1.10634 | 1.10730 | 0.086% |
| 0.01 | 1.11693 | 1.12654 | 0.860% |
| 0.1 | 1.22540 | 1.33035 | 8.564% |

These variations describe this finite toy event in this selected field, not measured distortions of supernovae. They show why the current candidate must predict the environment and its evolution before a universal affine signal law can be claimed. The earlier source-loading result should therefore be read in two parts: weak signals can share a common propagation law in a fixed field, while radiation loading changes that field and can distort long events.

## 4. What physical clocks would read remains an explicit requirement

The ray calculation uses reference time. If the same field also changes an atomic clock's frequency per reference time by q(n), then both the measured redshift factor and measured local duration factor become

    S_measured = J * q(n_received)/q(n_emitted).

This is a detector-standard comparison, not a substitution of ordinary gravitational time dilation for the user's new time mechanism. A complete interaction must specify what light and the instruments measuring it both do.

For diagnostic scalings q=n^-p, the result file records p=0,1,2. In the example with driving energy 0.01, path 2 and launch 0, J=1.11693 and n_received=1.10690, with n_emitted=1. The measured factors are 1.11693, 1.00907 and 0.91162 respectively. These are different candidate matter responses; none is adopted by this test. The p=2 control recalls the previously derived atomic-frequency issue for a particular fixed-charge/mass matched electromagnetic action. It is not a rejection of all cumulative-time mechanisms.

The mathematical frequency/timing equality survives these consistent endpoint conversions, but the sign and magnitude of the observable shift can change. A localized void solution with physically established source and receiver standards could behave differently from this homogeneous rolling background plus spatial disturbance.

## 5. Next physical requirements

This pass establishes a common local carrier/event-time mapping for weak signals in the selected dynamical field. It also quantifies finite-event distortion when that field changes. It does not establish finite-probe source independence, a realistic radiation population, void localization, actual atomic coupling, the astronomical rate, capture or extra gravity.

Next derive and test an environmental field configuration and common source/detector matter response, then evaluate weak signals under a distributed radiation population and finite-probe feedback. Use one field solution to predict redshift, event shape and companion energy flow. Preserve the original timing estimator failures and all other observation requirements.

```sh
python research_work/results/weak-signal-timing/check.py
```

Requires NumPy and SciPy and the archived companion-backreaction result for a numerical cross-check. The offline suite includes this diagnostic; its checks are not observational validation.
