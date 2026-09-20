# HH-1: temporary photon-companion coupling

Declared before calculation, following the owner's request for companions that
attach to a photon temporarily and detach after pulling it toward a swirl.
This is an effective coupling-state hypothesis, not a demonstrated bound state
of a real spin-2 graviton and a photon.

Use MS-1's counted matter source and evolving spiral stream without changing
its transport parameters. Treat test photons as probes. Attachment-energy,
dispersion, recoil dynamics and microscopic cross sections remain unresolved;
constant light speed/energy are imposed in this diagnostic, not derived from
a physical binding reaction. Record the missing momentum impulse explicitly.
No dark matter, expanding universe, distance fitting or observational fit.

## Mean state
p is attached probability; m is the mean attached heading vector, weighted by
occupation. Starting p=0,m=0:
a=k_on*u_g*((1+n dot t_g)/2)^2,
p_dot=a*(1-p)-k_off*p,
m_dot=a*(1-p)*t_g-k_off*m,
theta_dot=kappa*cross(n,m), x_dot=n.
A captured heading is carried with the ray until release. Subsequent captures
can follow new local directions. m is not silently reset to a nonexistent
remote stream after the photon leaves it. |m|<=p<=1 is a checked invariant.
This is a mean-field closure; nonlinear trajectory means need not equal a
mean-state trajectory.

k_on=1,10,100; k_off=0.1,1,10; kappa=1,4,16; pitch=0.5,2,4;
positive swirl sign, fixed forward attachment gate. 81 parameter cases.
Ten impact parameters as MS-1:810 mean-state rays, dt0.02 to x=20 or time80.
Record coupling occupancy, direction memory, observer angles and positions,
delay, accumulated absolute turning and the required stream recoil.
Use the same paired inward/sideways decomposition as MS-1.
Report all cases, not only favorable inward bends.

## Discrete events, separately
Three fixed representatives (pitch,k_on,k_off,kappa):
(2,10,1,4), (4,100,0.1,16), (0.5,1,10,1).
At impacts -2,-1,+1,+2, use128 seeded realizations each:1536 histories.
Free photons attach with probability1-exp(-a*dt); on attachment store the
local stream heading. Bound photons detach with probability1-exp(-k_off*dt).
A bound photon turns toward its carried heading; a detached one flies straight
until another attachment. Only one transition per step is allowed. Events are
applied at the step midpoint; finite-step event errors are checked by repeating
at dt0.01 and0.005 with independent seeds, not hidden by matched random draws.
Archive means, standard deviations, nonarrivals, attachment counts and all
trajectories. Angular scatter divided by mean bending is descriptive, not an
observational acceptance threshold or a p-value.

## Checks
- Uniform-stream occupation p=a/(a+b)*(1-exp(-(a+b)t)) and m=p*t_g.
- Zero capture produces straight light; zero guide coupling does not bend.
- Mean-state positivity and |m|<=p throughout every recorded trajectory.
- Compare six fixed mean-state parameter cases at dt0.01: the three stochastic
  representatives plus (0.5,100,0.1,16),(4,1,10,1),(2,10,0.1,4).
  Observer angle/position tolerance0.01 absolute; record failures.
- For stochastic refinement compare all12 ensemble means using their combined
  standard errors and absolute differences. Do not interpret statistical
  agreement as a proof of negligible blur or exact microscopic dynamics.
- Photon momentum impulse must be debited to the stream; a recorded debit is
  not an evolved closed momentum reservoir.
- Rates here are explicitly frequency independent. Achromaticity by this
  assumption does not establish a frequency-independent microscopic binding
  law, and a real attachment mechanism could change color, speed or timing.

Source hashes, declared seeds and all failed outcomes must be saved and checked
in. Prior alignment, relaxation kinetics and random capture/release mathematics
are credited as established constructions. No unprecedented-theory claim.
