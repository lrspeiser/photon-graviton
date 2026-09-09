# Spherical propagation changes the sustained-source result

## Why test geometry?

The previous one-dimensional calculation found a rolling companion field and nearly steady redshift under sustained illumination. A line confines how the field spreads. A star emits into three-dimensional space. This calculation retains the photon/field interaction while allowing spherical spreading, to test whether the rolling result survives.

This is a spherically symmetric sector of a three-dimensional scalar model, not a simulation of arbitrary angular illumination or ordinary spin-2 gravitons. The production profile remains specified rather than derived from weak gravity. Fixed spatial geometry omits the field's gravitational backreaction; that omission is particularly relevant to the user's proposed graviton self-binding.

## 1. Radial Hamiltonian and photon energy transfer

Set reference light speed c=1. For a spherical scalar field phi(r,t), define w=r phi and let each outgoing radiation shell have total canonical radial momentum P_i and radius X_i. The shell represents photons distributed over all directions, not a single spherical particle. Its energy is E_i=P_i/n_i.

```
n(X,t) = 1 + integral_0^infinity W_sigma(r-X) g(r) w(r,t)/r dr,
H_field = (4 pi K)/2 integral_0^infinity [w_t^2 + v^2 w_r^2] dr,
H_total = E_emitter_fuel + H_field + sum_i P_i/n_i.
```

The physical field energy before transformation is (4 pi K)/2 times the integral of r^2[phi_t^2+v^2 phi_r^2]. Integration by parts gives the w expression when the boundary term w^2/r vanishes at zero and infinity. The finite numerical domain uses w=0 at its boundaries, and the outer boundary is beyond relevant wave propagation during each run. Regularity at the origin is w(0)=0; inward waves pass through the physical center according to this regularity condition, rather than encountering an artificial absorbing hole.

The equations derived from this Hamiltonian are

```
Xdot_i = 1/n_i,
Pdot_i = P_i (partial_X n_i)/n_i^2,
w_tt - v^2 w_rr
  = sum_i [P_i/(4 pi K n_i^2)] W_sigma(r-X_i) g(r)/r.
```

Equivalently, the phi equation contains the three-dimensional radial Laplacian phi_rr+2 phi_r/r and a source proportional to 1/r^2. That area factor is part of the geometry and energy normalization, not an adjustable attenuation coefficient.

Along every shell,

```
dE_i/dt = -E_i (partial_t n_i)/n_i,
dH_field/dt = sum_i E_i (partial_t n_i)/n_i.
```

The companion field receives exactly the photon energy lost. No new conversion efficiency is selected. During emission, a finite amount epsilon is removed from the declared emitter fuel and assigned to a new shell at the source radius. Its canonical momentum is epsilon*n_source, so its actual energy is epsilon even if the field has changed at the emitter. The emission schedule is an input, not a derived atomic-clock or stellar-fusion law.

The fixed radial profile can exchange local radial momentum. Perfect spherical symmetry has zero total vector momentum, but that symmetry is not a derivation of the omitted material supports or a general recoil solution.

## 2. Why a localized steady source can become static in three dimensions

For the massless wave equation phi_tt-v^2 Laplacian(phi)=S, the three-dimensional retarded response is

```
phi(x,t) = integral S(y,t-|x-y|/v)/(4 pi v^2 |x-y|) d^3y.
```

The standard retarded Green function and its support on the wave cone are derived in the University of Manchester's [time-dependent wave-equation notes, equations 4.20-4.22](https://oer.physics.manchester.ac.uk/AM/Notes/jsmath/Notesse12.html). The factor v^-2 here follows from using partial_t^2-v^2 Laplacian as the operator rather than v^-2 partial_t^2-Laplacian. Applying that kernel to our candidate source is the present derivation.

If a spatially localized source becomes time independent, and no continuing homogeneous incoming field is imposed, the retarded integral at a fixed location eventually becomes time independent too. This differs from the one-dimensional kernel that integrated a steady source into a rolling field. This argument is conditional on the source becoming steady; the coupled simulations test the approach instead of prescribing a static field.

In the candidate's steady continuous radial photon flow, partial_t n=0 implies unchanged photon energy and unchanged arrival intervals. Local luminosity is then Q. Shell number per unit radius equals the number flux times n, so the continuum source is

```
S(r) = Q g(r) C(r)/(4 pi K r^2),
C(r) = integral_source^infinity W_sigma(r-X) dX.
```

The Gaussian integral C is retained rather than replaced by one near the emitter. The solution regular at the center and vanishing at infinity is

```
phi_static(r) = Q/(4 pi K v^2)
             * integral_0^infinity g(s) C(s)/max(r,s) ds.
```

This formula predicts the static field amplitude without a fit. Outside compact support it is proportional to 1/r. Substitution into the optical factor also predicts a nonzero excess travel time integral(n_static-1) dr.

## 3. A fixed delay is different from continuous stretching

For the same source and detector standards, the ray arrival map satisfies

```
dt/dr = n(r,t),
dJ/dr = (partial_t n) J,
frequency stretch = local event-duration stretch = J.
```

Once n is static, J=1 through the region even if n differs substantially from one. Each signal receives the same extra travel delay. Light sent ten days apart then still arrives ten days apart, shifted later by a common delay. During startup, partial_t n is nonzero, and redshift can accumulate and survive exit.

This is not ordinary gravitational endpoint time dilation and is not cancellation of an earlier redshift. A photon retains the shift it acquired while the propagation field was changing. The problem is that newly emitted photons crossing an already settled region acquire little further shift from this candidate.

The static spatial field still stores energy, but it is neither a demonstrated population of captured gravitons nor a derived gravitational halo. A static store does not by itself supply the continually changing optical factor required here.

## 4. Scope of the finding

The numerical comparison in [report.md](report.md) checks finite fuel, photon/field backreaction, two spatial grids, two shell spacings, a luminosity control and a larger radial extent. The localized rolling behavior is not preserved: the runs approach the static prediction and later redshifts become much smaller than startup redshifts.

A cosmic distribution of many sources is not equivalent to a compact source around one center. A nonlocalized source, finite causal history, evolving illumination, different boundary state, nonlinear self-interaction or dynamical spacetime can change the result. Each requires its own conserved equations and predictions; none is ruled out or adopted here. In particular, the user's proposed companion self-binding is absent from the free scalar propagation model and must not be claimed tested by this calculation.

The next step should examine a specified companion self-interaction or self-gravity mechanism and a distributed finite source history, keeping binding, continued time stretching and the energy budget distinct. Positive attraction alone does not establish capture, long-lived support, a suitable halo shape or an ongoing redshift law.
