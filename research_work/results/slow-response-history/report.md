# Slow response does not sustain stretching in the fixed-through-flow test

## Outcome

Increasing the response time from one to ten or one hundred does not preserve accumulated redshift in this tested branch. Traveling energy is replenished while the temporal field approaches a steady profile. Deposits continue growing, but the model does not yet allow their gravity to change the temporal field or the capture geometry.

An independent stationary calculation reproduces the late numerical shifts. This distinguishes a physical limit of the specified equations from a time-integration artifact. It is not a proof about all source histories, an isotropic universe, or a model with deposit feedback.

## Single changed postulate

Retain the [radiation-powered initiation model](../radiation-powered-initiation/report.md), including whole-photon conversion I=epsilon R, but set

\[
Q=(T-M)/\tau,\qquad \tau=10\ \mathrm{or}\ 100
\]

instead of tau=1. **This is a parameter sensitivity of a project postulate, not a new first-principles law or a fitted astronomical timescale.** K=2, sigma=0.5, epsilon=0.01, capture strength ten, the two capture regions, constant radiation inflow and zero initial T/M seed are unchanged. The domain length is ten and time runs to forty in the same nondimensional units. The slower temporal-sector energy flux remains explicit.

No source, seed or capture parameter is adjusted for individual emission times. The tau=1 results are reused exposed controls. Six new spatial simulations use 80,160 and320 cells; each tests three signal emission times.

## What the signals do

| Response time | Early measured z, emission0.1 | Late measured z, emission20 | Late photon survival |
|---|---:|---:|---:|
| 1, previous control | +0.00107293 | -0.00085772 | 0.904642 |
| 10 | -0.00034467 | -0.00253196 | 0.904720 |
| 100 | -0.00006897 | -0.00054529 | 0.904818 |

All values above use 320 cells. The late coordinate event-stretch factors are approximately one. The remaining measured shifts are primarily endpoint clock differences, and all are blueshifts here. Changing the response time does not eliminate the roughly 9.5% photon removal imposed by the unchanged initiation channel.

The maximum change of n between reference times30 and40 is 5.94e-13 for tau=10 and6.49e-10 for tau=100 on the fine grid. Small residual derivatives depend on numerical resolution; they are not evidence for an astrophysically persistent signal. The causal explanation in this finite channel is replacement of outgoing energy under fixed inflow, rather than unlimited local accumulation. A long exchange time alone does not make energy remain at one location.

## Independent steady-flow calculation

**Conditional derivation using ordinary flux algebra:** for stationary R,T,M define

\[
F_R=R/n,\quad F_T=T/n,\quad F_M=\sigma M/n.
\]

Since n=1+M/K,

\[
n=\frac{1}{1-F_M/(\sigma K)},\quad R=nF_R,\quad T=nF_T,\quad M=K(n-1).
\]

The stationary equations reduce to spatial ordinary differential equations:

\[
F_R'=-\epsilon R,\qquad
F_T'=-Q-\Gamma T+\epsilon R,\qquad
F_M'=Q.
\]

Boundary fluxes are (0.22,0,0), fixing n at the inlet to one. The flux lost from the traveling sectors equals the deposit growth rate, the integral of Gamma T over x. D itself is **not stationary**: it keeps accumulating. Its effect on gravity has not been included in these equations.

Under the declared observer clock rule, a stationary field has coordinate Jacobian J=1 and measured factor S=n_in/n_out. Thus

\[
z_{\rm steady}=1/n_{\rm out}-1=-F_{M,\rm out}/(\sigma K).
\]

The sign here follows from the specified inlet condition and positive outgoing temporal energy. It is not a universal sign rule for every spatial history. The usual conditional clock/ray identities supply this prediction; no novelty is claimed for the calculus.

| Response time | Independently calculated steady z | Difference from fine evolved z | Continuing deposit growth rate |
|---|---:|---:|---:|
| 1 | -0.000857755 | 3.44e-8 | 0.01974918 |
| 10 | -0.002532324 | 3.68e-7 | 0.01817362 |
| 100 | -0.000545359 | 6.86e-8 | 0.02018763 |

The stationary solver also independently predicts travel time and photon survival. All its flux ledgers close within 4.45e-16. Agreement improves with the spatial refinement; all fine-grid z differences pass the separately declared 1e-5 stationary-comparison gate. This verifies the observed late-state interpretation without asserting global attraction for arbitrary nonlinear initial data.

## Verification and reproduction

[Protocol](protocol.md), [time evolution](run.py), [evolution results](results.json), [stationary solver](stationary.py) and [stationary results](stationary.json) are archived. Run:

```text
python research_work/results/slow-response-history/run.py
python research_work/results/slow-response-history/stationary.py
```

The stationary comparison reads the existing tau=1 initiation archive as well as the new results. Reservoirs remain nonnegative; characteristic-domain checks pass. All original energy, proper-clock/event and 0.003 absolute-z refinement gates pass. These are synthetic reference-energy and kinematic checks, not observations, physical momentum conservation or a full gravity model. The astronomical alpha calibration and holdout status are unchanged.

## Consequence for the next physical step

Further slowing the same local response is not supported as the fix by this test. Two physically distinct extensions remain: a common evolving radiation/environment history, and feedback of the accumulated deposits on gravity and clocks. The latter is already required by goal5 and cannot remain absent in a claimed complete model.

The next derivation should make the deposit-to-gravity/clock response explicit, including its sign and its energy cost. Determine whether it sustains the desired light stretching or produces an opposing clock effect. A new spatial potential must also predict stellar acceleration and light bending; adding an arbitrary positive redshift term would bypass that requirement. No microscopic account of why time exists is needed to state and test these operational rules. All nine goals remain incomplete.
