# FM-1: finite-speed, matter-funded scalar/vector mediator
Declared 20 September 2026 before implementation and calculation.

## Next step and boundaries
Replace SV-1's instantaneous finite-packet kernel with local wave fields, derive a source term from matter, and derive massive/massless test-particle response from the same positive Hamiltonian. Execute a planar (two spatial dimensions) consistency experiment, NOT an astronomical fit or full theory. Preserve all earlier empirical failures. No dark matter, expansion, prescribed halo, distance changes or optical normalization fit.

## Hamiltonian (dimensionless c=1)
Fields F=(phi,A_x,A_y), momenta Pi, wave speed c_g=1, field frequency omega=0.2:
H_field = integral [|Pi|^2+c_g^2 |grad F|^2+omega^2 |F|^2]/2 d^2x.
For ordinary source particles and probes:
E_i=sqrt(m_i^2+|p_i|^2); u_i=p_i/E_i;
sigma_i=g phi_bar(q_i)+eta A_bar(q_i) dot u_i;
H_i=E_i exp(sigma_i); H_total=H_field+sum H_i.
Bar denotes a normalized compact source/interpolation profile, same for sampling and deposition. In the continuum it is a translated finite-radius profile; this is an extended-particle effective model, not point-local microphysics. It permits instantaneous sampling within a source radius a, so only propagation between nonoverlapping source supports may be called wave-mediated. A point-local continuum completion remains separate.
The full energy is nonnegative, but positive energy alone does not establish kinetic convexity or Lorentz invariance. Particle group speeds can differ from c; measure them.

Derived equations:
qdot_i=exp(sigma_i)[u_i+eta(A_bar-u_i(A_bar dot u_i))]
pdot_i=-H_i[g grad phi_bar+eta (grad A_bar)^T u_i]
Fdot=Pi
Pidot=c_g^2 Laplacian F-omega^2 F
       -sum_i H_i (g, eta*u_ix, eta*u_iy) f_a(x-q_i).

Thus stationary matter sources phi, moving matter also sources A. Initially empty fields are energized at the expense of particle Hamiltonian energy; no extra carrier reservoir is inserted. This is a postulated matter coupling, not a demonstrated photon-conversion process. Field-dependent rest energies/clock response require future local observational constraints.

Weak-field static elimination of the quadratic field gives
delta H_eff = -1/2 integral J(x) G_omega(x-y) J(y) dx dy,
J=(g E, eta p_x, eta p_y).
The vector cross term is proportional to -eta^2 G_omega p_i dot p_j, the qualitative SV-1 coupling with a Green function replacing its Gaussian/particle-count normalization. This is a leading-order static reduction, not equality at finite coupling or finite frequency. Reversing eta for every particle changes A's sign but should leave trajectories unchanged.

For massless probes H=|p| exp[g phi_bar+eta A_bar dot n], homogeneous of degree one in p: in a prescribed field its trajectory is independent of momentum magnitude. Test finite probe backreaction separately. This derivation is not evidence of image sharpness or a measured lensing law.

## Discretization and sources
Periodic square [-8,8)^2, primary N=48 cells per side, central nearest-neighbor Laplacian. Fixed physical profile radius a=0.75 in each coordinate; profile product of (1-(d/a)^2)^3 for |d|<a, normalized by its discrete sum. Differentiate this normalization in particle forces. Deposit field source with cell-area factors. Fixed radius is a declared smoothing scale, not reduced with grid spacing.

Six ordinary particles of mass 1 at radius 0.8, equally spaced. Cases: resting, rotating with p=0.2 e_phi, reverse rotating. No external potential or fixed support. Scalar g in {0,0.08}; vector eta in {-0.08,0,0.08}; probes in {none, photon, massive}: 54 primary runs. Fields and Pi start at zero.
Photon: m=0, p=(1e-4,0), q=(-2.5,0.7).
Massive probe: m=1e-4, initial physical free velocity 0.4 along x, q=(-1.5,0.7). Its momentum is m*v/sqrt(1-v^2).
These are diagnostic paths at different speeds/locations, not matched astrophysical star and photon observations.
RK4 to t=4, dt=0.02. Field source supports and wave travel must stay clear of periodic wraparound; explicitly check and report any limitation.

## Declared checks and follow-ups
Before simulation: finite-difference full discrete Hamiltonian gradients for field, Pi, q and p (relative <=1e-5); zero coupling; zero-source field; source deposition integral; discrete integration by parts; exponential/source energy positivity.
Principal-part propagation: independent 1D compact pulse, exact massless d'Alembert reference at t=2, grids 128/256/512 on [-8,8), RK4 with dt<=0.2 dx. Relative max error <=0.03; ahead-of-front amplitude outside support+ct+2dx <=0.001, retaining failures. This validates a numerical wave operator; it does not remove extended-source nonlocality.
All primary cases: scaled energy drift <=1e-4. Record field gain versus matter energy loss, maximum physical speed, field amplitude, circulation/curl, probe direction and final position, total momentum and angular momentum residuals. Continuum momentum is sum p minus integral Pi grad F. Angular momentum includes orbital field terms and vector-field spin A_x Pi_y-A_y Pi_x. The square grid/profile break exact rotational symmetry; report residuals and convergence rather than claiming exact discrete momentum conservation.
Refine 18 cases (g=.08, eta=0 or .08, every source/probe) at N=64, dt=.01. Also run N=48, dt=.01 for the same 18 to separate time from spatial refinement. Numerical targets: max position discrepancy <=0.01 model units; field-energy relative change <=0.03 (normalize by max(fine field energy,1e-10)). Numerical failures remain failures.
Double photon momentum to 2e-4 for all 6 g=.08, eta=0/.08 source choices. Difference in final physical direction <=0.001 radians; this is a finite-backreaction color diagnostic, not a broadband observational pass.
Record global-eta sign reversal and rotation reversal identities. No hypothesis/settings selected from a favorable result.

## Attribution and interpretation
Wave equations, Hamiltonian mechanics, scalar/vector couplings, compact interpolation, weak-field Green-function elimination and homogeneous ray Hamiltonians are established mathematics. Cite primary/expert sources; no historical novelty claim and no relabeling of an established gravity theory as our invention. The particular exponential coupling is a declared project candidate.
Success means a tested finite-speed mediator with explicit source accounting and derived probe dynamics under the stated assumptions. It does not mean a point-local relativistic completion, exact reproduction of SV-1, a self-generated stable swirl, successful microscopic photon attachment, or improved galaxy/cluster scores.
