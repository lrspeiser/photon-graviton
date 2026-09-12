# Isotropic depleted supply: definitions and derivation

The opacity law `kappa=kappa0 [W/(1+W)]^p` is a project postulate. Exponential absorption, spherical Poisson gravity and the following energy accounting are known mathematics applied to that postulate. No uniqueness is claimed. Companions are modeled as straight traveling rays; deposition adds mass-equivalent density in place. This is a transport/capture diagnostic, not a derived microscopic interaction.

Let I be the isotropic incident specific intensity integrated over companion frequencies. Its units are energy/(area time solid angle). For a sphere of radius R, inward flux through each surface element is pi I. Hence total incoming power is `4 pi R^2 * pi I`.

Write `B=4 pi I T/c^2` for exposure T in mass-equivalent units, and normalized exposure s=t/T. Then incoming mass-equivalent energy per unit s is `B pi R^2`. In the unattenuated limit `d rho_D/ds=B kappa`, so the earlier combined parameter is `C=B kappa0`. Holding C fixed while changing kappa0 changes both opacity and incident supply. It is a test of previously unidentified factors, not a test in which the radiation field stays fixed.

For spherical opacity, each incident direction has the same chord absorption profile. Its projected area integral `dA=2 pi b db=pi d(b^2)` represents one direction; multiplication by the 4 pi angular integral is already included in B. A chord at impact parameter b has a half-length in shell j:

`L_j=sqrt(max(r_outer,j^2-b^2,0))-sqrt(max(r_inner,j^2-b^2,0))`.

Define `tau_j=kappa_j L_j`, `T_half=sum_j tau_j`, and sums `T_lower,j` and `T_higher,j` over shells inside and outside j. The fraction of incident chord energy absorbed in shell j on the inward and outward legs is:

`A_j=(1-exp(-tau_j)) [exp(-T_higher,j)+exp(-T_half-T_lower,j)]`.

The transmitted fraction is `exp(-2 T_half)`. Summing all shell absorption gives `sum_j A_j=1-exp(-2 T_half)`. This telescoping identity is independent of the opacity's dependence on deposited gravity at that instant.

The shell density evolves as `d rho_j/ds = B integral A_j dA / V_j`. Gravity and opacity are recomputed as density evolves. Therefore:

`d M_D/ds + d E_transmitted,eq/ds = B pi R^2`,

and for zero initial deposits over unit exposure, `0<=M_D<=B pi R^2`.

This upper bound scales as R^2, so a balanced supply ledger on a finite domain alone does not prove an isolated finite-mass limit. In particular, an absorbing outer layer can intercept most of the area-dependent input.

The solver assumes the incident ray distribution equilibrates rapidly relative to deposition. It has no explicit radiation-in-flight time derivative, gravity propagation, momentum exchange, mechanical support or binding-energy correction. The identity is an absorption ledger in mass-equivalent units; it is not a full relativistic energy conservation proof. A physical source model must still supply I and T, account for energy spent in state changes/support, and check that the quasistatic transport approximation is justified. Even a balanced absorbed-energy budget does not show that photon redshift supplies the required incoming companion intensity.

## Conditional momentum consequence: incoming capture does not provide outward support

Suppose the traveling companions carry momentum E/c along their rays, as photons do, and absorption transfers that momentum locally. This is an explicit additional companion postulate; the following follows from ordinary energy/momentum accounting, not a new gravity law.

In dimensional units let q(r)>=0 be absorbed energy per volume per time and F_r the net radial energy flux. With no internal companion emission and stationary transport, spherical energy continuity gives:

`d(r^2 F_r)/dr=-r^2 q`, hence `F_r(r)=-(1/r^2) integral_0^r q(u) u^2 du <=0`.

For isotropic scalar opacity, absorbed radial momentum per volume per time is `f_rad=kappa F_r/c <=0`. Incoming absorption therefore pushes inward in this spherical experiment. It cannot supply the outward support needed to hold a deposited particle distribution fixed. This is not a numerical momentum-conservation test; it identifies a missing term and its sign under the stated companion-momentum assumption.

If the stored state is approximated as a static isotropic fluid with density rho_D and inward gravitational acceleration magnitude g, its necessary force balance is `dP/dr=-rho_D g+f_rad`. A zero-pressure outer boundary would require `P(r)=integral_r^R [rho_D(u) g(u)-f_rad(u)] du`. An equation of state, energy source for that stress and stability still need specification; writing this pressure integral alone does not supply them. A bound-particle distribution or an anisotropic field stress would require its corresponding dynamical equations instead. No such support is implemented in the current capture solver.
