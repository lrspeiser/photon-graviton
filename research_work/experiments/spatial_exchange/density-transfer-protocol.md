# SE-DT1: scalar-gradient-driven local transfer

Declared 20 September 2026 before calculation. JE1 generates very little
K/Q energy for slow sources. Test a structural change to the joint Hamiltonian:

    K = P - lambda A - zeta grad phi - epsilon(curl A) cross Q.

Keep all other JT1 energy terms and particle interactions unchanged. Use
centered grad phi in this kinetic shift. The exact reciprocal scalar reaction
is an additional -zeta div(V_Q) in pi_phi_t, where V_Q=Q_t=aK-beta dot grad Q.
All other derivatives are those of JT1 evaluated at P-zeta grad phi. This
substitution is a way to differentiate the new Hamiltonian; it is not a claim
that the shifted variables constitute a canonical transformation.

Energy stays a positive quadratic block in Pi,K,grad phi,grad A,grad Q for
the same bounded beta. For a frozen direction n, group W=(phi,A). The derivative
coupling is D=[zeta n, epsilon T], a3x4 matrix with
T=(n dot Q)I-n Q^T. The ST1 matrices generalize to blocks of dimensions4,3,4,3:

    M=[[-bI4,-aD^T,I4,bD^T],[-aD,-bI3,bD,c_Q^2 I3],
       [aI4,0,-bI4,0],[0,aI3,0,-bI3]],
    S=[[aI4,0,-bI4,0],[0,aI3,0,-bI3],
       [-bI4,0,I4,0],[0,-bI3,0,c_Q^2 I3]].

The longitudinal A direction remains in ker D, so reference photon speeds
b+-C should still be covered. This does not prove global nonlinear stability.

An empty-field initial kick has phi_tt=-g rho_eff and A_tt=-d j, where
d=kappa c_Q eta, j=sum W p/dV and rho_eff=sum W(E+p^2/E)/dV,
E=sqrt(m^2+p^2). Consequently

    Q_ttt(0)=lambda d j + zeta g grad rho_eff.

The new term can seed the companion sector at zero matter velocity. This is
a proposed change in response, not evidence that it produces attractive or
enhanced lensing; the reciprocal term could instead screen the field.

Three declared checks:
1.48 full joint energy derivatives: n8/12,L8,zeta0/.5/2,eight fixtures each,
  epsilon2,lambda1,seed20261008, same random-state and particle setup as JT1.
  Steps1e-6,5e-7, both errors and their difference<2e-7 scaled. For zeta0,
  require full RHS and energy reduce to JT1 within1e-12.
2.192 frozen principal cases: C=.5/2,c_Q=.5,epsilon0/2,zeta0/.5/2,
  |Q|=.03,|beta|/(c_Q C)=0/.9,eight random orientations each. Require positive
  S,scaled symmetry<1e-12,imaginary spectrum<1e-10, photon longitudinal-mode
  residual and extremal coverage<1e-10.
3.18 initial-kick cases: n8/12,zeta0/.5/2,source momentum0/.0007/.2. Six
  mass1 particles on the radius.7 ring; no photon for this source-isolation
  diagnostic. Starting with empty fields, independently reconstruct rho_eff
  and j using full-grid cardinal weights. Compare Q after one RK4 step to
  Q_ttt(0)*dt^3/6 at dt1e-4 and5e-5; normalized residual of6Q/dt^3 against
  predicted third derivative must be<1e-4, denominator max(1,max(abs(pred))).

Preserve all first results. No long production evolution, stronger attraction,
fuel sufficiency or observational success follows from these tests. Nonzero
Q at rest would distinguish this seed from the velocity-suppressed one; it
does not establish useful swirl. No dark matter, expansion or old-gravity
target. Kinetic shifts, gradient couplings and Hamiltonian calculus are
established structures and are not claimed as historically new.
