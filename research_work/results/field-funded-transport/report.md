# Letting temporal-field energy pay for radiation exchange

## Outcome

This assignment fixes the previous empty-companion problem in the tested spatial domain. All energy reservoirs stay nonnegative, the reference-energy ledger balances, and temporal energy travels out or is captured instead of remaining as stationary memory. Observable frequency changes agree with independently measured event stretching.

It does **not** yet explain the desired astronomical redshift. Signals beginning inside the supplied seed region redshift, while signals from the domain boundary blueshift or show almost no change. The separate radiation ledger shows a small net transfer **into radiation**, not net extraction from it. Deposits in this experiment are funded by the supplied seed budget overall. These distinctions prevent a numerical conservation pass from becoming a false redshift/gravity claim.

## Proposed equations and provenance

**New project constitutive choice, not claimed unique:** radiation exchanges energy directly with the energy M of the temporal field. T is an additional traveling companion reservoir; R is shared electromagnetic/gravitational-wave reference radiation energy; D is deposited reference energy. Set n=1+M/K, v=1/n, b=R/(Kn), and Q=(T-M)/tau:

\[
\begin{aligned}
R_t+(vR)_x&=-bM_t,\\
T_t+(vT)_x&=-Q-\Gamma T,\\
M_t+(\sigma vM)_x&=bM_t+Q,\\
D_t&=\Gamma T.
\end{aligned}
\]

K is an energy-density scale, tau a response time, Gamma a capture rate and sigma the temporal energy flux speed relative to local c0. The calculation uses nondimensional units c0=1, K=2, tau=1 and sigma=0.5. Thus photons, GWs and T share local c0; **M has a separately postulated slower energy flux**. This is not an ordinary-graviton identification or a demonstrated tensor field theory.

M's energy is explicit and transported. Clock readings still obey the proposed d tau_clock=dt/n rule with fixed spatial rods. A microscopic origin of time is not required to test these operational postulates. Physical momentum and a complete measured-energy/stress prescription remain missing.

**Conditional consequences, using known conservation and characteristic mathematics:**

\[
\partial_t(R+T+M+D)+\partial_x[v(R+T)+\sigma vM]=0,
\]

\[
dM_t+\frac{\sigma}{n^2}M_x=Q,\qquad d=1-b.
\]

At T=0 with T_x=0, T_t=-Q=M/tau>=0 regardless of the M gradient. The previous counterexample therefore gives +0.1 rather than a negative companion derivative. This repairs that specific donor assignment.

The repair restores a possible singularity at d=0. The coordinate principal matrix for (R,T,M) is

\[
\begin{pmatrix}
v&0&-R/(Kn^2)-b a_M\\
0&v&-T/(Kn^2)\\
0&0&a_M
\end{pmatrix},\qquad a_M=\frac{\sigma}{n^2d}.
\]

The temporal-state characteristic travels at sigma/(nd) relative to local c0. We require d>0 and nd>sigma along these trajectories, separating that characteristic from the repeated radiation/T speed. The smallest recorded margins are d=0.87791 and nd-sigma=0.38371. Domain checks are made at every solver right-hand-side evaluation; reported minima use accepted steps. This is not a proof that arbitrary initial data stay in that domain.

## Spatial experiment

The domain is x=0..10, evolved to reference time 24. Initial radiation is R=0.22, M=D=0, and the specified seed is T=0.1 exp[-((x-2)/0.4)^2]. Incoming radiation turns off smoothly around time 2. There is no incoming T or M. Capture uses Gamma=1/[1+exp(-(x-8)/0.25)], or Gamma=0 for comparison. A third case has no seed at all.

These are controlled synthetic conditions, not a model of the real cosmic radiation history or a fitted galaxy. Probes are negligible test signals; they do not add appreciable energy, including the late probe after background inflow shuts off.

**Known ray/clock identities applied to the postulated field:**

\[
\frac{dt}{dx}=n,\qquad\frac{d\ln J}{dx}=n_t,
\quad 1+z_{\rm measured}=J\frac{n_{\rm emit}}{n_{\rm observe}}.
\]

J is the derivative of arrival reference time with respect to emission reference time. An independent check integrates source and receiver proper-clock intervals for signals emitted at te +/-0.001. It agrees with the measured spectral factor within 1.28e-5 across all archived runs. Shared EM/GW delay is zero by the identical propagation postulate; this does not predict an astrophysical source emission lag.

Fine-grid results, with 640 cells for seeded runs:

| Probe (start, end, emission time) | No capture: measured z | Capture: measured z |
|---|---:|---:|
| (0, 10, 1) | -0.0284511 | -0.0201074 |
| (2, 9, 0.1), inside seed region | +0.1101740 | +0.1094322 |
| (0, 10, 5), late probe | -0.000001983 | -0.000001658 |

All zero-seed probes have z=0. There is no spontaneous initiation in that control. The tiny late shifts are effectively near zero at the declared absolute-z refinement tolerance; their signs should not be treated as a precision observational prediction. Positive seed-overlapping redshift is a conditional result, not evidence for a universal distance-redshift law.

## Where did the deposited energy come from?

**Additional conditional energy identity:** integrate signed radiation-to-M work

\[
W=\int dt\,dx\,bM_t.
\]

Radiation then satisfies R_final+R_out+W=R_initial+R_in. Positive W would mean net radiation energy supplied to the field; negative W means the reverse.

| Fine-grid budget quantity | No capture | Capture |
|---|---:|---:|
| Initial seed energy | 0.070898154 | 0.070898154 |
| Net radiation-to-M work W | **-0.001553807** | **-0.001567008** |
| Final deposit energy | 0 | 0.048698995 |
| Final traveling T+M | less than 2e-47 | less than 2e-47 |

There is no net radiation funding of the combined deposit/export budget here. The initial seed supplies deposits, exported nonradiation energy and the small net radiation gain. Individual exchanges can occur in both directions; the integrated result cannot establish a net photon-powered gravity reservoir. All quantities in this table are reference energies in the declared nondimensional model, not inferred galactic masses.

## Reproduction and limits

[run.py](run.py), [protocol.md](protocol.md) and [results.json](results.json) retain all results. Run both commands:

```text
python research_work/results/field-funded-transport/run.py
python research_work/results/field-funded-transport/run.py --refine
```

The first command computes the 80/160/320-cell runs; the second adds the two 640-cell seeded refinements. The existing coarse archive predates the separate radiation-work ledger; regeneration adds that diagnostic to those runs without changing the physical equations. Fine runs include both ledgers.

The 160-to-320 seed-overlapping z differences exceeded the original 0.003 gate, so refinement was required. Those failures remain in the archive. At 320-to-640, the maximum difference is 0.001644 and both seeded cases pass that gate. No negative energy was recorded. Total and separate radiation budgets agree to at most 8.0e-15. This tests numerical consistency over the specified trajectory, not global nonlinear stability or physical energy/momentum conservation.

## Next physical question

The remaining obstacle is now sharper: what supplies and evolves the temporal field so that ordinary source-to-observer paths retain the desired redshift, while radiation is a net energy donor? Repeated seed placement tuned to each source would not answer it. Derive and test a common initiation/background history, including the radiation exchange sign and clock endpoints, before calibrating the distance relation. The spatial force/lensing law, real-data joint tests and genuinely unseen evaluation remain required. No observed data, fitted alpha or holdout status changed in this experiment.
