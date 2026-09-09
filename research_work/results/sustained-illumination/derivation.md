# Sustained illumination: a rolling field with nearly steady signal stretching

## The physical question

Can the proposed cumulative time/companion mechanism stretch successive signals by nearly the same amount under sustained illumination, even though its field keeps changing? The previous finite train smoothed fluctuations but retained a substantial drift. This calculation derives a possible long-time behavior from the same Hamiltonian and compares it with longer, finite, energy-accounted source histories.

The distinction matters for an explosion: a common stretch preserves its timeline up to an overall duration factor, while a stretch that changes during the explosion distorts its brightening and fading. Matching photon redshift alone does not establish the event's appearance.

The following is a conditional derivation within an optional one-dimensional scalar-field candidate. It does not derive the user's low-gravity environmental cause from fundamental physics, identify ordinary gravitons, establish matter clocks, or provide capture and gravity. The source and region scales remain inputs. A computational photon packet is a radiation bundle, not a new particle species.

## 1. Start from the tested interaction

Use reference light speed c=1, field inertia K>0, companion-wave speed v>0, and a static nonnegative environmental profile g. The normalized smoothing kernel W has unit integral. The canonical photon positions and momenta are X_i and P_i. Define

```
n(X,t) = 1 + integral W(x-X) g(x) phi(x,t) dx,
H = integral K/2 [phi_t^2 + v^2 phi_x^2] dx + sum_i P_i/n(X_i,t).
```

Hamilton's equations, already tested in the finite-train calculation, give

```
Xdot_i = 1/n_i,
Pdot_i = P_i n_X(X_i,t)/n_i^2,
phi_tt - v^2 phi_xx = sum_i P_i W(x-X_i) g(x)/(K n_i^2).
```

The photon energy E_i=P_i/n_i therefore obeys dE_i/dt=-E_i n_t/n_i. Its lost energy is exactly the work done on this field. No independent energy-conversion efficiency or extra gravitational energy multiplier is introduced. Static environmental supports can exchange momentum; their dynamical material realization is still absent.

## 2. A conditional rolling branch

Suppose sustained illumination eventually produces a field whose velocity is spatially uniform over the interaction support:

```
phi(x,t) = u t + B(x),   u>0,
g_bar(X) = integral W(x-X) g(x) dx,
n_t(X,t) = u g_bar(X).
```

This is an ansatz to test, not a stability proof. It permits n itself to grow. Along a ray dt/dX=n, the derivative J of arrival time with respect to emission time obeys dJ/dX=n_t J. With common reference source and detector standards,

```
G(X) = integral from upstream to X of g_bar(Y) dY,
J(X) = exp[u G(X)],
E(X)/E_in = 1/J(X).
```

Thus a constant u gives a constant frequency and local duration stretch at a fixed downstream location, despite the continually growing field. The travel time of successively emitted signals still increases: the arrival map is affine with slope J, not a constant travel delay. Finite source and detector positions approximate the asymptotic endpoints when the kernel tails there are negligible.

Let Q be constant incoming photon power and let q(X) be the local power. Both photon energy and photon arrival rate decrease by J, giving

```
q(X) = Q exp[-2u G(X)].
```

This relation includes the arrival-time effect; using only Q/J would count energy per photon but omit the changed arrival rate.

## 3. Derive the rolling rate, rather than fitting it independently

In a continuous right-moving beam, photon number per unit length equals local number flux times n. Multiplying by P/n^2=E/n makes the continuum source proportional to local power:

```
S_phi(x) = g(x)/K integral W(x-X) q(X) dX.
```

For a steady spatial source of the one-dimensional massless wave equation with outgoing waves on both sides, the retarded solution has phi_t=(1/(2v)) integral S_phi dx once the source's causal influence spans the region. This follows directly from integrating the retarded kernel (1/(2v)) times the integral of S_phi over the past light cone; differentiating the time integral leaves the entire spatial source after it is causally included. Equivalently, B''=-S_phi/v^2 and the outgoing slopes are B'=-u/v on the right and +u/v on the left.

Because the kernel is normalized, G_total=integral g_bar dX=integral g dx. Interchanging the integrals gives

```
integral S_phi dx
  = (1/K) integral g_bar(X) Q exp[-2uG(X)] dX
  = Q [1-exp(-2uG_total)]/(2Ku).
```

Consequently the rolling branch must satisfy

```
4 K v u^2 = Q [1-exp(-2uG_total)],
S_total = exp(uG_total).
```

For Q>0 and G_total>0, the zero root introduced by multiplication by u is not a valid sourced rolling solution. Divide by u and use the finite limit at u=0. There is one positive solution: [1-exp(-2G_total u)]/u decreases from 2G_total to zero, while 4Kvu increases. This is uniqueness within the ansatz, not uniqueness or stability of all solutions of the original evolution equations.

In the weak regime uG_total much less than one,

```
u approximately Q G_total/(2Kv),
z approximately Q G_total^2/(2Kv).
```

The source power and field parameters therefore determine the redshift together. A universal astronomical coefficient has not been derived. The source structure also matters: for the three cosine-squared profiles with half-width 0.5, G_total=1.5. This branch does not justify multiplying independently reset single-region factors.

## 4. Where the apparently missing light energy goes

The outgoing field waves carry energy at rate Kvu^2 on each side, for a total field-energy growth rate 2Kvu^2 in an expanding causal region. In a fixed control region, the same term appears as outward wave-energy flux after the local rolling state is established. Thus

```
field energy growth rate = 2Kvu^2 = (Q-Q_out)/2.
```

The other half of the power difference is not an unaccounted loss. Photons slow in reference coordinates as n grows, so more photon energy remains in transit. The beam energy per unit length is n q. On the rolling branch q is time independent and

```
d/dt integral n q dX
  = integral u g_bar Q exp[-2uG] dX
  = (Q-Q_out)/2 = 2Kvu^2.
```

Hence the full control ledger is

```
incoming photon power - outgoing photon power
  = companion-field production + growth of photon energy still in transit.
```

Every unit actually lost by an individual photon still enters the companion field. The total brightness deficit at a detector is a different quantity: it includes photons whose arrival is delayed. Assigning the entire deficit to new companion energy would overcount the supply in this branch. An unlimited steady source would require unlimited fuel; the numerical tests use finite initial photon energy and make no such supply assumption.

## 5. Tests and interpretation

See [report.md](report.md) for the frozen finite-source comparisons, numerical checks and measured agreement with this branch. The independent frequency and event-time integrations use the same declared reference standards. Actual atomic clock behavior, three-dimensional illumination, environmental dynamics, companion capture, stable support, and a joint rotation/lensing law remain required.

Even a successful late-time comparison is not evidence that the universe has reached this state. A rolling solution has growing n and increasing photon inventory; it is not a static universe-wide propagation law. Turning the finite source off, changing its luminosity, changing g or adding an absorbing receiver changes the equations or state and requires a new joint test.
