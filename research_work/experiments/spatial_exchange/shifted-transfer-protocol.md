# SE-ST1: compatibility of local transfer with matter/light propagation

Declared 20 September 2026 before calculation. The LR3 transfer subsystem
currently has no coupled photon or matter evolution. Test a proposed shared
coefficient extension, not a claim that full integration already exists.

Let C>0,a=C^2,0<c_Q<=1, and bound the drift |beta|<c_Q C. Retain
K=P-lambda A-epsilon(curl A) cross Q. A candidate local energy density is

    H = a(Pi_A^2+K^2)/2 + (|grad A|^2+c_Q^2 |grad Q|^2)/2
        - beta_j [Pi_A dot partial_j A + K dot partial_j Q] + V,

where V is nonnegative. The advection current uses K, not canonical P, in
the Q term. This is a deliberate model change, not an algebraic relabeling
of the old canonical-P shift. Completing squares shows positive energy for
the stated beta bound. A scalar phi sector can use the same a and beta with
its own positive kinetic/gradient block; its nonlinear source derivatives
must be included when implemented.

For a frozen smooth background and direction n, define b=beta dot n,
T=(n dot Q)I-n Q^T. With u=(Pi_A,K,partial_n A,partial_n Q), the proposed
principal matrix u_t=M partial_n u is

    M = [[-bI,-a epsilon T^T,I,epsilon b T^T],
         [-a epsilon T,-bI,epsilon b T,c_Q^2 I],
         [aI,0,-bI,0],[0,aI,0,-bI]].

The energy Hessian is

    S = [[aI,0,-bI,0],[0,aI,0,-bI],
         [-bI,0,I,0],[0,-bI,0,c_Q^2 I]].

S should be positive and SM symmetric. Physical directional propagation
speeds are minus eigenvalues of M. Because T n=0, longitudinal A modes
include speeds b-C and b+C, independent of epsilon; other branches can be
faster. This is coverage of the particle reference cone, not equality of all
field speeds and not a complete nonlinear stability proof.

Pair provisionally with h_m=sqrt(C)*sqrt(m^2+C|p|^2)+beta dot p.
For m0 the velocity is beta+C p/|p|; for m>0 it lies within that velocity
ball. Derive reciprocal field sources from this same h when implementing
matter, rather than adding a force separately. A source-production law and
physical fuel budget for Q/P remain missing. Spatially local matter would
require local fields/distribution dynamics or a controlled point-source limit;
the earlier finite shared-rotor averaging is not a causality proof.

Fixtures: seed20261006, C=.5,1,2; c_Q=.25,.5,1; epsilon0,2,10;
|Q|=.03,.3; |beta|/(c_Q C)=0,.5,.9; four random orientation sets each,
total648. Require min eigenvalue(S)>0, scaled symmetry error<1e-12,
max imaginary eigenvalue<1e-10, longitudinal-mode residual<1e-10,
and extremal field speeds cover b+-C within1e-10. Check particle velocities
for m0/1 with random nonzero p against centered momentum derivatives at
steps1e-6 and5e-7, scaled error<1e-7 and relative-drift speed<=C+1e-12.
No observed curve or older gravity formula is used as an acceptance target.

Positive Hamiltonian calculus, advection, eigenvalue tests and the particle
dispersion ansatz are established mathematical structures. Their use does not
establish historical novelty of a gravity theory. This test does not complete
the coupled equations, real source production, light-quality predictions,
full principal constraints or global evolution.
