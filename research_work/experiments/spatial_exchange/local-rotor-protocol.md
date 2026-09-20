# SE-LR1: local rotor interaction sector

Declared 20 September 2026, before calculation. SE-Rot4 demonstrated an
instantaneous cross-source response from a single rotor sampled across a
finite kernel. Replace that interaction sector by pointwise vector fields
Q(x), P(x). This is not yet a replacement for the full matter/gravity model.
It introduces no dark matter component or expanding background. Its tests
use internal mathematical consistency, not agreement with older gravity laws.

Define B=curl A, K=P-epsilon B cross Q, and the Hamiltonian density

    H = (Pi^2 + |grad A|^2 + omega_A^2 A^2)/2
      + (K^2 + c_Q^2 |grad Q|^2 + Omega^2 Q^2)/2.

The canonical equations are

    A_t = Pi
    Pi_t = Laplacian A - omega_A^2 A + epsilon curl(Q cross K)
    Q_t = K
    P_t = c_Q^2 Laplacian Q - Omega^2 Q - epsilon B cross K.

These equations describe local field exchange. Q is provisionally an internal
excitation, not a demonstrated property of ordinary matter or an identified
graviton. Matter production, fuel, the scalar sector, light coupling and
observer definitions are not supplied here. No pre-existing halo is assumed.

For a unit propagation direction n and frozen background Q0, define
T=(n dot Q0) I - n Q0^T. In variables (v_A,v_Q,d_A,d_Q), where v are time
derivatives and d directional spatial derivatives, the principal matrix is

    M = [[0, -epsilon T^T, I, 0],
         [-epsilon T, 0, 0, c_Q^2 I],
         [I, 0, 0, 0], [0, I, 0, 0]].

S=diag(I,I,I,c_Q^2 I) is a positive symmetrizer for c_Q>0. For each singular
value sigma of T, the squared speeds are roots of

    w^2 - (1+c_Q^2+epsilon^2 sigma^2) w + c_Q^2 = 0.

Use the product relation for the smaller root to avoid cancellation. Frozen
real speeds and this symmetrizer address the principal interaction sector;
they do not prove global nonlinear stability or causality of the full earlier
finite-source model. Speeds need not all equal the unit vacuum A-wave speed.

Preregistered fixtures (seed 20261002): Q magnitudes .03,.1,.3; c_Q=.25,.5,1;
epsilon=0,.4,2,10; eight random direction pairs each, total288. Require
symmetrizer residual <1e-12, numerical eigenvalue imaginary parts <1e-10,
analytic/numerical speed difference <1e-10, positive squared speeds, and
agreement after a common proper rotation <1e-10. Report the full speed range
and faster-than-unit branches without classifying agreement with an older
gravity theory as a gate.

Independently check all four canonical equations against centered directional
energy derivatives on periodic n8 and n12 grids, L8, omega_A=.2, Omega=1,
c_Q=.5; epsilon=0,.4,2,10; six random states each:48 cases. Forward-difference
gradient energy, nearest-neighbor Laplacian and centered curl give an exactly
consistent discrete Hamiltonian. Use perturbations 1e-6 and5e-7; require both
scaled errors <2e-7 and their difference <2e-7. Field amplitudes .03 and unit
Euclidean random test directions. No time evolution is claimed.

This is a candidate coupling choice using established Hamiltonian calculus,
curl identities, singular values and symmetric-hyperbolic analysis, not a
claim that these mathematical techniques or local rotor fields are novel.
The related prior-art discussion in rotor-report.md remains applicable.
Preserve first results and do not alter gates after seeing them.
