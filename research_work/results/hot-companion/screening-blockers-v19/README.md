# One rate mechanism for screening and delayed release

## Status

This package contains **new, executed reduced-model calculations** responding to the question:

> Can a physically specified pair of rates produce the screening shape and release length together, without choosing them independently?

Two concrete stochastic constructions generate the requested constant-field response. The more detailed construction ties both factors to one blocker lifetime, but predicts a different response to a changing field. **Neither construction independently predicts the numerical gravitational constants, proves a microscopic origin of gravity, or calculates a gravitational force.** No astronomical data were fitted, no Cassini or wide-binary solver was run, and the adopted repository, law and blog were not changed.

The reference targets are the user's proposed screening `exp(-g/g_d)` and release `1-exp(-r/L)`, with `g_d=2.03e-10 m/s^2`, `u=169 km/s`, `L=0.15 pc`. These are input calibration values, not discoveries of this calculation. The source analysis introduces scalar relaxation in its lines 257–304; the current notebook lists these constants in lines 145–150 and identifies screening as assumed in lines 168–183. Here `g` means the scalar field controlling screening, provisionally `|g_N|` as in that notebook.

## 1. What any literal two-state model can determine

For open fraction R and nonnegative rates a and b:

```
dR/dt = a(g)*(1-R) - b(g)*R
R_equilibrium = a/(a+b)
tau = 1/(a+b)
```

Requiring the exact equilibrium `exp(-g/g_d)` and a constant time `tau0` fixes:

```
a = exp(-g/g_d)/tau0
b = (1-exp(-g/g_d))/tau0
```

Writing those equations is inverse design, not a derivation. Multiplying both rates by a common positive function leaves equilibrium unchanged and changes the kinetics. Equilibrium information alone therefore cannot determine the release time. Thermodynamic constraints on a rate ratio do not in general remove this kinetic freedom.

### A microscopic refresh construction

At independent Poisson refresh events of rate nu, draw two independent unit Gaussian internal quadratures X and Y. Let their energy be

```
E/E_star = (X*X+Y*Y)/2.
```

Let the coupling be open after the refresh exactly when E exceeds a field-dependent threshold `epsilon(g)=zeta*g`. Retain that state until the next refresh. Since a two-dimensional Gaussian's quadratic energy has an exponential distribution:

```
P(E > epsilon) = exp(-zeta*g/E_star)
a = nu*exp(-zeta*g/E_star)
b = nu*(1-exp(-zeta*g/E_star))
g_d = E_star/zeta
tau = 1/nu
L = u/nu
```

Thus the literal two-state pair has a possible specified origin. However, the equilibrium energy scale and refresh clock are independent properties until another microscopic calculation relates them. The construction assumes full independent refresh, a linear threshold, and a state-independent refresh clock. A generic two-level thermal system does not automatically have these properties. Refreshing internal energy entails exchanges with a reservoir not dynamically simulated here.

`reset_gate_stochastic.csv` tests the process by sampling Gaussian quadratures at refresh events. It does not assign `exp(-g/g_d)` as the gate-opening probability. Five seeds of 200,000 realizations each were used per field value; every sampled open fraction was within 2.33 binomial sampling standard errors of the prediction.

## 2. A more structured candidate: temporary blocking excitations

Postulate a small coupling element carried with the companion. It can contain an integer number n of metastable excitations that prevent it from coupling. These are called blockers for clarity, not as a claim that their physical identity has been established.

There are **two elementary transition laws**:

```
n -> n+1     at rate lambda(g) = eta*g
n -> n-1     at rate gamma*n
```

The creation events are independent; each excitation disappears independently with lifetime `1/gamma`. The coupling element is active only when n=0. Blocker multiplicity is retained. This is not, in general, a memoryless two-state projection onto only open/closed.

### A possible local capture realization

Suppose blockers are created by capturing independent packets from the ordered component of the companion. For packet energy E_b, capture cross-section Sigma and incident ordered energy flux F_ord:

```
lambda = Sigma*F_ord/E_b.
```

In the ordered-stream limit already used by the notebook, `F_ord = ell*g/(4*pi*G)`. Then:

```
eta = Sigma*ell/(4*pi*G*E_b).
```

A blocker holds the captured energy until it decays into a specified bath or channel. A launch blocker is paid for by the source. This gives an interpretation for eta and gamma rather than assigning a desired equilibrium probability.

This is **a proposed realization, not a derivation from the existing wave code**. It assumes the relevant independent packet flux exists and is the ordered component controlling screening. A real coupling that captures total wave power, including heat, could screen hot radiation differently. Coupling selectivity, saturation, packet correlations, energy exchange, recoil and feedback all require explicit calculation. Independent blocker lifetimes and field-independent capture cross-section are substantive assumptions.

## 3. Derivation of the two factors

### Equilibrium screening

The stationary master equation implies

```
lambda*p_n = gamma*(n+1)*p_(n+1)
p_n = exp(-x)*x^n/n!,       x=lambda/gamma
```

Therefore:

```
R_equilibrium = P(n=0) = exp(-x) = exp(-g/g_d)
g_d = gamma/eta.
```

No exponential field dependence appears in the creation or removal rate: it emerges from the probability of zero independent blockers.

### Release after emission

Assume that, at launch, the element inherits the **stationary ambient blocker population at the local field** plus **one additional identical blocker** created by emission. In a constant field, the ambient population remains stationary. The launch blocker survives with probability `exp(-gamma*t)`. Hence:

```
R(t,g) = exp(-g/g_d)*(1-exp(-gamma*t))
L = u/gamma
```

for advection at constant speed u. The same gamma controls screening occupancy and launch activation; no separately chosen release lifetime is used.

The stationary preparation plus one added blocker is important. It is not established by the two transition laws alone. One can test it in a source model; it must not be silently assigned because it reproduces a target.

### Initial-condition controls

Starting with exactly one blocker and NO equilibrated ambient population instead gives:

```
R(t,g) = (1-exp(-gamma*t))*exp[-x*(1-exp(-gamma*t))].
```

For x=10 and t=0.1/gamma this is 0.0367435, compared with 0.00000432037 for stationary ambient plus one launch blocker. The latter is the desired product. The difference is physical preparation, not numerical error.

Starting with m independent launch blockers on top of ambient gives:

```
R(t,g) = exp(-x)*(1-exp(-gamma*t))^m.
```

Thus the observed initial activation shape could, in principle, constrain the launch process rather than only a lifetime.

## 4. Tests actually run

### Full master-equation check

`run_checks.py` constructs the complete birth/death generator with 81 and 121 states, independently evolves the probability vector with a matrix exponential, and compares the zero-occupancy probability with the derived product.

Parameters: x = 0, 0.1, 0.3, 1, 3, 5, 10; times gamma*t = 0.1, 0.5, 1, 3, 5; both truncations.

- Maximum absolute error in released fraction: **1.22e-14**.
- Maximum probability-normalization error: **1.79e-14**.
- Largest population at the artificial upper boundary: **4.65e-43**.

These errors check the stochastic mathematics and implementation, not the physical truth of the postulates.

### Independent event/count simulation

Five seeds (11, 23, 37, 53, 71), each with 200,000 independent elements: **1,000,000 realizations per field/time condition**.

The ambient process starts empty, then evolves for 25 blocker lifetimes before a launch blocker is added. It is not initialized with the desired open fraction. Existing blockers survive as independent exponential lifetimes. Births are independently Poisson. Their random birth times and survival are aggregated exactly over intervals, which avoids a timestep approximation. At no point is `exp(-x)` used to set whether a blocker element opens.

At t=1/gamma:

| x=g/g_d | Predicted released fraction | Simulated fraction | Sampling standard error |
|---:|---:|---:|---:|
| 0 | 0.632121 | 0.632369 | 0.000482 |
| 0.3 | 0.468286 | 0.468628 | 0.000499 |
| 1 | 0.232544 | 0.232437 | 0.000422 |
| 3 | 0.0314714 | 0.0312100 | 0.0001746 |
| 5 | 0.00425919 | 0.00426800 | 0.0000651 |
| 10 | 0.0000286982 | 0.0000270000 | 0.00000536 |

All tested values were within **2.17 sampling standard errors**. Rare openings in the highest field have larger relative sampling error; this is not a sub-percent experimental determination there. Measurements at several times on the same realization are correlated.

Blocker counts satisfy `n_final-n_initial=births-deaths` exactly. Attaching an energy per blocker gives a bookkeeping identity for this reduced process. It does **not** independently close the matter-plus-medium energy or momentum budget.

## 5. How the numerical scales become related

From `g_d=gamma/eta` and `L=u/gamma`:

```
g_d*L = u/eta = u*v_b,       v_b := 1/eta.
```

The target constants imply:

| Quantity | Inferred value |
|---|---:|
| gamma | 3.651278e-11 s^-1 |
| 1/gamma | 867.8629 years |
| eta | 0.1798659 s/m |
| v_b=1/eta | 5.559697 m/s |
| v_b/u | 3.289762e-5 |

`v_b` is a velocity-dimension rate coefficient, not a separately discovered propagation speed.

For the illustrative packet-capture realization, using `ell=a*u/2`, the same inputs require:

```
E_b/Sigma = a*g_d*L/(8*pi*G) = 35288.47 J/m^2.
```

This is an independently testable microscopic combination **only after** a microscopic model predicts the packet energy and capture cross-section. Inferring it from g_d and L is not an independent prediction.

The single-scale guess `v_b=u` would imply `L=u^2/g_d=4559.60 pc`, not 0.15 pc: a factor 30,397 difference. It is that guess, not the general blocker model, that fails this calibration. Numerical agreement with the present theory needs an additional small dimensionless ratio, whose physical origin remains open.

**Parameter counting:** two empirical quantities have been re-expressed through two microscopic quantities. There is functional unification and a useful relationship, but no reduction in fitted freedom yet.

## 6. A decisive new prediction: history response

With field g(t) along the companion trajectory, an initially Poisson ambient population stays Poisson with mean m(t):

```
D_t m = eta*g - gamma*m
D_t s = -gamma*s
R = (1-s)*exp(-m)
```

Here s is the probability the additional launch blocker still survives. For an old companion (s=0):

```
D_t R = -gamma*R*[log(R)+g/g_d].
```

So **the logarithm of the released fraction**, not the released fraction itself, relaxes linearly. This is not the scalar relaxation law in the previous proposed model.

For a mature element initially equilibrated at x0, suddenly moved to x1:

```
R_blocker(t) = exp[-x1-(x0-x1)*exp(-gamma*t)]
R_reset(t) = exp(-x1) + [exp(-x0)-exp(-x1)]*exp(-gamma*t).
```

They have the same equilibrium curve and weak-field activation timescale but different finite transients. For x0=10, x1=0:

| t/tau | Individual blockers | Refreshed gate |
|---:|---:|---:|
| 1 | 0.0252534 | 0.632137 |
| 3 | 0.607824 | 0.950215 |
| 5 | 0.934840 | 0.993262 |

Time to 50% release is 2.66910 tau for blockers versus 0.693102 tau for refresh. The independent master-equation calculation also reproduces these step responses.

### Small-signal measurements do not distinguish the two models

For small perturbations around equilibrium g0:

```
delta R(omega)/delta g(omega)
  = -[exp(-g0/g_d)/g_d]/(1+i*omega*tau)
```

in both constructions. Direct sinusoidal ODE checks at omega/gamma=0.1, 1 and 10 confirm the linear response to relative errors around 1.3e-7 or smaller for an input of amplitude 0.001*g_d. These errors are dominated by finite input amplitude, not an asserted universal numerical precision.

Finite field changes, state statistics or fluctuation spectra are needed to distinguish the mechanisms. This sharpens the earlier suggestion to measure a transfer function: a linear transfer function alone is insufficient here.

Stationary normalized gate autocorrelation for the blocker model is

```
C(t)/C(0) = [exp(x*exp(-gamma*t))-1]/[exp(x)-1],
```

whereas a reset gate has `exp(-gamma*t)`. These correlations were checked against the master equation. Observable force noise additionally depends on the number of independent coupling elements and their actual force coupling.

### Spatial path control

`spatial_history.csv` compares the two kinetic models and an endpoint product in a prescribed field `g/g_d=(0.3L/r)^2`, for different starting radii. The blocker model can retain very strong suppression far beyond one L if it begins with a large ambient blocker population. E.g. at start r=0.1L (x0=9), it gives R=0.01032 at r=L, versus 0.40728 for scalar relaxation. At r=10L they nearly converge (0.998265 vs 0.998741).

This is a normalized path experiment, not a solar calculation. It shows why a lifetime length must not be identified with the actual release radius after an arbitrary strong-field history. The equilibrium and launch preparation cannot be separated from an eventual Cassini calculation.

## 7. Robustness and unresolved physics

1. **Creation law:** eta*g must follow from a physical flux and interaction or be declared an effective postulate. A rate linear in all companion intensity could also depend on the heat term, unlike the adopted g_N-only screening.
2. **Initial preparation:** inheritance of ambient blockers and creation of exactly one launch blocker need to be derived from emission, not selected after inspecting the desired transient.
3. **Independent blockers:** finite capacity, correlated formation, cooperative decay or field-dependent lifetimes change the result. For K independent sites with creation rate per vacant site eta*g/K and removal gamma per occupied site, the empty probability is `(1+x/K)^(-K)`, not exactly exp(-x). It approaches the Poisson result as K increases. The analytical capacity sensitivity is saved in `finite_capacity_control.csv`.
4. **Meaning of the release variable:** the calculation gives an active probability. Its identification with a force multiplier must follow from the local wave/receiver interaction. Suppressing source power gives sqrt(R) in an amplitude-force model; suppressing coherent source amplitude gives a different combination rule. No power-to-force conversion was simulated.
5. **Conservation and passivity:** incoming captured energy, the launch cost, decay output, force feedback and recoil must be explicit in the actual medium model. Count conservation alone is not proof of these.
6. **Advection:** L=u/gamma assumes the gating state travels at u. A gate bound to stationary matter does not have that same spatial activation equation.
7. **Local meaning of g:** the candidate uses the screening field of the current phenomenology. Its coupling must be defined in the material/medium frame rather than by an arbitrary observer's acceleration.
8. **Finite-number fluctuations:** a small number of coupling elements would have telegraph-like or occupancy noise. A smooth gravitational multiplier needs an averaging argument and noise bounds.
9. **No astronomy verdict:** no claim is made that either candidate improves Cassini, the Milky Way, dwarfs, binaries or lensing. Each is a new law outside constant conditions and requires the actual field solver.

## 8. Recommended next experiment

Use the existing local matter/medium model as the proposed reservoir and ask whether it supports metastable blocking configurations. Measure **formation rate versus the field** and **single-excitation decay lifetime** independently. Do not fit g_d and L separately in that experiment.

Then test: constant field, strong-to-weak step, weak-to-strong step, a sequence of steps, and different launch preparations. Compare measured gate occupancy and the **actual force** with the two candidate predictions. Keep the local energy and momentum accounts explicit. The force should not be set by multiplying an existing answer by R.

Only after one candidate survives should its derived state equations be introduced as a separately named alternative into the full astrophysical solver, first with the current calibration fixed.

## Reproduction

Tested: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0. Exact library versions are also recorded in `results/summary.json`. Output timings can differ between runs; numerical outputs should reproduce within floating-point and RNG-library conventions.

```
python -m pip install -r requirements.txt
python run_checks.py
python make_plots.py
```

Default five seeds x 200,000 realizations. For a smaller sanity check:

```
python run_checks.py --n-per-seed 20000 --output-dir quick-results
```

The exact interval sampler does not approximate a finite timestep: it separately counts Poisson births and each individual's survival. The continuous probability generator and matrix exponential are independent checks. `results/stochastic_by_seed.csv` retains seed-resolved outputs. `make_plots.py` reads existing full results rather than running fits.

## Background literature, not evidence for gravity

- Julian Lee, *The origin of the Poisson distribution in stochastic dynamics of gene expression*, Physica A 630 (2023) 129201; arXiv:2206.01938. Mathematical background for independent production/removal, time-dependent Poisson populations and composition with surviving initial populations. https://arxiv.org/abs/2206.01938
- Dominic C. Rose, Hugo Touchette, Igor Lesanovsky, Juan P. Garrahan, *Spectral properties of simple classical and quantum reset processes*, Physical Review E 98, 022129 (2018); arXiv:1806.01298. Background for state-independent Poisson resetting. https://arxiv.org/abs/1806.01298

The Poisson and reset mathematics are established, not claimed as new. Their interpretation as a release/screening mechanism is a hypothesis investigated here. No claim of priority or physical validation follows from the toy-model agreement.
