# SE-O: candidate observer normalization for clocks and photon energy

Declared 20 September 2026. Current bundle arrival times and Hamiltonian photon
energies are coordinate quantities. Define an observer convention before using
them as measured delays or frequency shifts. This is an explicit extension of
the candidate, not a demonstrated physical clock mechanism.

For the regularized massive Hamiltonian H(q,p,m), define a clock phase by the
massive test-body action. Its rate per unit rest parameter is

    N(q,p,m) = partial H / partial m = sum W alpha m / E.

Euler homogeneity H=p dot H_p + m H_m implies along its Hamiltonian velocity
v=H_p that (H-p dot v)/m=N. Adopt d tau=N dt as the observer clock convention.
For a prescribed observer velocity v_obs, find its canonical momentum by solving
H_p=v_obs, at m=1. A photon with coordinate energy H_gamma and momentum k is
assigned observed energy (H_gamma-k dot v_obs)/N. This is the local exchange
normalization obtained by differentiating the observer's H with respect to
momentum and rest parameter. The convention's actual material-clock and
electromagnetic realization remains to be established.

Test 48 seeded fixtures, seed20260927, cycling n24/32, kernel radius .9/1.2 and
observer v=(0,0,0) or (.1,-.03,.02), on the same smooth analytic phi,A fields as
SE-P. Use random positions in [-2,2]^3 and random nonzero photon momenta.
Solve the observer momentum using scipy root, tolerance1e-11. Require velocity
residual<1e-10 and N>0. Compare N with centered mass derivative of the directly
assembled Hamiltonian at m=1 +/-1e-5, error<1e-8; Euler identity error<1e-12.
Compute observed photon energies and verify their exact scaling for twice the
photon momentum within1e-12. This does not test physical dispersion/polarization.

Add two uniform controls: U=0,beta=0,v_obs=0 gives N=1 and E_obs=|k|;
uniform nonzero U,A with static observer has
N=alpha sqrt(1-|beta|^2/C^2), C=exp(2U), derived by solving the same H_p=0.
Require error<1e-10. Do not use an external gravity theory as target.

Do not retrofit observed delays onto current bundle archives: the time history
of the field at the detector is not archived. Future observer-qualified bundle
runs must integrate detector N over coordinate time and evaluate emission and
detection photon energy using each observer's state. Report raw coordinates
alongside clock-normalized values and test the convention experimentally later.
Hamiltonian homogeneity and local energy exchange are established mathematics;
the choice of action-based clock convention belongs to this candidate.
