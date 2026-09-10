# Internal storage can evade the restricted transport obstruction, at a stated cost

**An extra internal degree of freedom admits an explicit energy-conserving, shared-speed construction.** A conventional quadratic reservoir does not stay synchronized with the reference light ray; a deliberately linear reservoir does. This changes the next theoretical question from whether such equations can be written to whether the special internal state can be physically created, coupled to light, and captured without spoiling the other requirements.

This is a synthetic transport calculation, not an astronomical fit. No galaxy observations, distances, training parameters or holdout outcomes were changed or opened. It does not establish a strong observational case for the full theory.

## Hypothesis and formula provenance

The earlier [shared-speed restriction](../shared-speed-obstruction/report.md) explicitly excluded evolving internal states. Here we change that assumption. For one forward-moving packet, introduce position/momentum `(x,p)` and internal coordinate/momentum `(s,P)`:

`H = p/(1+a*s) + R(P)`, with `p>0`, `a>=0`, `1+a*s>0` and reference `c0=1`.

**Status: hypothetical effective packet Hamiltonian.** It is not a graviton action, a derivation of time, or a novelty claim. The established mathematical method of adding a time coordinate and conjugate momentum to make a time-dependent Hamiltonian autonomous is described in [Sussman and Wisdom, section 5.5](https://mitp-content-server.mit.edu/books/content/sectbyfn/books_pres_0/9579/sicm_edition_2.zip/chapter005.html). Interpreting the additional variable as a physical reservoir traveling inside a companion is a further hypothesis, not a consequence of that rewriting.

**Known Hamilton equations applied to this proposal:**

`x_dot=1/(1+a*s)`, `p_dot=0`, `s_dot=R'(P)`, `P_dot=p*a/(1+a*s)^2`.

The reference photon has prescribed speed `v_gamma(t)=1/(1+a*t)`. It is a comparison ray, not a dynamically coupled energetic photon field. All parameters and the interval `0<=t<=6` are dimensionless diagnostic choices, without astronomical calibration. The internal variable `s` is not an ordinary atomic clock or an established proper time.

**Conditional energy identity derived using known calculus:**

`d[p/(1+a*s)]/dt = -p*a*R'(P)/(1+a*s)^2`,

`dR/dt = +p*a*R'(P)/(1+a*s)^2`, hence `dH/dt=0`.

Energy moves from the propagation term into the packet's internal reservoir. There is no spatially separate void deposit in this point-packet model. Whether the entire energy actually travels this way in a local field theory still needs an energy-current/stress calculation.

## Conventional quadratic reservoir: conserved energy, unequal speed

Choose `R(P)=P^2/(2M)` and initialize `s=0`, `P=M`, so `s_dot=1` initially. The reservoir already carries `M/2` energy. For positive `a,p`, it gains energy and its rate increases: `s_dot>1`, `s>t`, and therefore `v_c<v_gamma`.

An independent reduction follows from conservation:

`s_dot=sqrt[1+2(p/M)*(1-1/(1+a*s))]`.

Integrating its reciprocal independently reproduces the full Hamiltonian trajectory. At `a=.05`, `t=6`, `p=1`, the reference index reaches 1.3:

| M | Pre-existing reservoir energy / initial propagation energy | Final companion speed / reference photon speed | Retained total excess energy |
|---:|---:|---:|---:|
| 1 | 0.5 | 0.9721314073 | 100% |
| 10 | 5 | 0.9971161691 | 100% |
| 100 | 50 | 0.9997106070 | 100% |
| 1000 | 500 | 0.9999710506 | 100% |

Retained excess means `H-R(P_initial)`, divided by initial propagation energy. Subtracting the seed is an accounting comparison; it does not make that seed physically absent, nongravitating or free to create. Increasing M reduces this branch's speed error but raises its seed energy. No measured speed bound or physical resource limit is inferred from this illustrative table.

## Ideal linear reservoir: an exact mathematical escape

Choose `R(P)=E0+P` for `P>=0`, with `P_initial=0` and `E0>0`. Then `s_dot=1` irrespective of the reservoir's acquired energy. For the selected positive-transfer evolution:

`s=t`, `P=p*[1-1/(1+a*t)]`, `v_c=v_gamma`,

`H=p+E0`, `x=ln(1+a*t)/a` (or `x=t` at `a=0`).

**Status: exact conditional solution of the stipulated Hamiltonian, using the known extended-coordinate construction.** It is not a uniquely new redshift formula. At the same `a=.05,t=6`, the propagation term retains 76.9231% of its initial energy and the internal reservoir gains 23.0769%. Their combined excess retains 100%, while the speed ratio is exactly one. The numerical calculation verifies this for initial p=.25,1,4, and a=0,.05,.1. There is no energy-specific retuning.

It would be incorrect to dismiss this construction solely because a linear function extended to all negative P is unbounded below. An explicitly positive, twice continuously differentiable reservoir completion is:

`R(P)=E0*[1+(sqrt(pi)/2)*erf(P/E0)]` for `P<0`,

`R(P)=E0+P` for `P>=0`.

Its global lower bound is `E0*(1-sqrt(pi)/2)>0`. Its derivative on the negative side is `exp[-(P/E0)^2]`; at zero its value, first derivative and second derivative match the linear side. The positive-transfer trajectory remains at P>=0. This is a project construction from the known error-function integral, not a claimed fundamental energy spectrum. The total H is positive on the declared forward-packet, positive-index domain. This is not a proof of stability or positivity of an eventual interacting field theory outside that domain.

## What the exact option demands

For a>0, shared speed along the synchronized path requires s=t. Differentiating `s_dot=R'(P)=1` gives `R''(P)*P_dot=0`. Since P_dot>0, **the reservoir must have exactly constant slope over the visited interval in this separable family**. Constant slope is an assumption to explain physically, not a parameter determined by observed redshift.

For a positive-domain perturbation `R=E0+P+epsilon*P^2/(2E0)`, keeping the same initial synchronization and E0=p=1, energy remains conserved but final speed ratios become 0.9997106070 for epsilon=.01 and 0.9971161691 for epsilon=.1. These are sensitivity choices, not observational exclusions. Another interaction Hamiltonian could change this conclusion; the condition is specific to the family tested here.

The construction has no specified reservoir capacity ceiling. In this isolated episode its gain is bounded by the initial propagation energy p. Multiple loading interactions, inverse processes and preparation costs require an actual interaction model. The earlier [coherent receiver audit](../coherent-receiver-audit/report.md) supplies a different, finite quantum receiving-state example; its optical-coherence results cannot be imported into this classical packet model without a derivation.

## Consequence for the full research goal

The prior no-internal-state restriction remains valid within its assumptions. The quadratic reservoir is an unfavorable realization, but **energy retention plus shared speed is no longer an unconstructed possibility once we permit the ideal internal state**. Equally, introducing that state is not yet explaining a graviton or photon redshift.

Before this can improve the theory's observational standing, a joint mechanism must answer:

1. **Creation and synchronization:** how photon conversion prepares the reservoir, pays its seed energy and sets s correctly for packets created at different times and places. A clock reset to zero at each birth need not track a common environmental clock.
2. **Observable light:** how the same interaction predicts frequency, packet-arrival histories, brightness and local clock/ruler measurements. Conserving this packet's H does not solve the earlier photon clock-cancellation or supernova-duration problems.
3. **Transport and gravity:** what energy and momentum current carries both terms, whether this can be made into a local field model, and how capture deposits the energy with the required gravitational response. Conserved canonical p in this homogeneous diagnostic is not a complete emission/capture momentum balance.

These are the next discriminating requirements. A further sweep of M alone would not address them. No claim is made that the index 1.3 is an observed redshift or that the model reproduces galaxy lensing/rotation. The three-dimensional stellar likelihood and held-out galaxy tests remain outstanding. The total universe photon-supply budget remains deferred; the explicitly present seed must nevertheless stay in any later energy ledger.

## Reproduction and validation

Run `python research_work/results/internal-companion-reservoir/run.py`. It executes 23 base cases, an energy/momentum-scale control, a tighter integration, independent scalar quadrature for the quadratic branch, exact trajectory controls for the ideal branch, and a derivative check of the positive reservoir completion. Conservation residuals are normalized to initial propagation energy, not the potentially much larger seed. `results.json` includes the source hash and all numerical results. These tests verify the stated mathematics and its numerical solution, not physical existence or observational adequacy.
