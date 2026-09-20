# SE-JT1: scalar, local transfer, massive particles and photons

Declared 20 September 2026 before calculation. Implement the proposed ST1
Hamiltonian, not just its frozen matrices. Coordinates are phi,A,Q and their
canonical momenta pi,Pi,P, plus particle positions x_i and momenta p_i.

Set g=.08,eta=.08,kappa=.5,c_Q=.5,omega=.2,Omega1;
C=exp(2g phi), a=C^2, alpha=sqrt(C),
beta=C kappa c_Q eta A/sqrt(1+eta^2 |A|^2),
K=P-lambda A-epsilon(curl A) cross Q. Thus |beta|<c_Q C.

    H_field = integral [a(pi^2+Pi^2+K^2)/2
       +(grad phi^2+|grad A|^2+c_Q^2 |grad Q|^2)/2
       - beta_j(pi partial_j phi+Pi dot partial_j A+K dot partial_j Q)
       +(omega^2(phi^2+A^2)+Omega^2 Q^2)/2].
    H_particles = sum_i <alpha sqrt(m_i^2+C p_i^2)+beta dot p_i>_W.

W is the previously tested positive cardinal cubic point-probe interpolation,
with width shrinking with grid spacing and exact position derivative. Include
its reciprocal field deposition divided by cell volume. It is a numerical
regularization, not a proof of finite-grid rotation symmetry or locality.
Massive particles and a nonzero-momentum massless particle all backreact;
none is an externally prescribed force or a passive probe.

Use the exact Hamiltonian derivatives of C,beta,K and W. With
V_Q=aK-beta dot grad Q, the vector momentum equation includes
+lambda V_Q+epsilon curl(Q cross V_Q), and P_t includes
-epsilon(curl A) cross V_Q-div(beta K). The scalar momentum equation includes
-2g a(pi^2+Pi^2+K^2)+2g beta dot current, besides its wave/mass/particle terms.
The derivative of beta with respect to A supplies its reciprocal current term.
Do not omit these coefficient/source derivatives.

Discrete gradient density uses the symmetric forward/back allocation,
(|Dplus f|^2+|Dminus f|^2)/4, with c_Q^2 weight for Q; its total equals the
earlier nearest-neighbor energy. Use centered curl and centered currents.
The symmetric allocation bounds centered-gradient cross terms pointwise.

Checks:48 random full-state energy derivatives, n8/12,L8,epsilon0/2,
lambda0/1,six fixtures each,seed20261007. Field amplitudes.03; three particles
of masses1,.5,0 at random positions within[-.6,.6]^3, momenta random scale.2
for massive particles and.01 for the photon. Reject only exactly zero photon
momentum. Use unit Euclidean perturbations and steps1e-6,5e-7. Require both
scaled errors and their difference<2e-7. Include field cell-volume factors
and unweighted particle canonical pairing. Require positive total energy and
pointwise photon/matter interpolation cone bounds within1e-12.

Twelve zero-g/zero-eta controls: n8,all four epsilon/lambda pairs,three random
fixtures each. Require exact reduction (error<1e-12) to LR3 plus a free scalar
and freely moving particles. Report field-sector and particle-sector errors.

No evolution or empty-field production success is claimed by these derivative
tests. A later run must start all fields empty, use moving ordinary matter,
and track its energy loss, generated transfer field, light deflection and
angular momentum together. Bare particle masses are fixed; this does not yet
identify physical source fuel or guarantee galactic persistence. Full nonlinear
stability and the implemented joint principal system remain to be tested.
No dark matter or expansion. Hamiltonian calculus, interpolation and the
dispersion ansatz are established tools, not a historical novelty claim.
