# SE-Rot: finite angular reservoir coupled to field curl

Declared 20 September 2026 before numerical evaluation. SE-C found the scalar
wave mechanism barely alters circulation. Test a reciprocal angular coupling
as a new candidate, not as an established graviton mechanism.

Replace a scalar internal oscillator by canonical vectors Q,P and define
Bbar=sum W(q) curl A, K=P-epsilon Bbar cross Q,

    M = m + (|K|^2 + Omega^2 |Q|^2)/2.

Use this M inside the existing Hparticle=sum W[alpha sqrt(M^2+z p^2)+beta.p].
With R=partial Hparticle/partial M, the internal equations are

    Qdot = R K,
    Pdot = -R [Omega^2 Q + epsilon Bbar cross K].

The derivative M_B=-epsilon (Q cross K). Because discrete centered curl is
self-adjoint in the volume inner product, the additional vector momentum
source is +curl[W epsilon R (Q cross K)/dV]. The positional force includes
+epsilon R (Q cross K) dot sum (partial W/partial q) curl A, in addition to
the usual -sum (partial W/partial q) h. Retain all scalar/vector coefficient
sources from the existing particle energy with its new M.

Canonical internal angular momentum is Q cross P, not Q cross K. Add it to
the matter-plus-field Noether ledger. In empty fields, Q=sqrt(Eint)/Omega e_x,
P=+/-sqrt(Eint) e_y gives internal energy Eint and signed angular momentum
+/-Eint/Omega. The reservoir is finite; lowering Omega to add angular momentum
would change the physical source assumption and must be declared, not free fuel.

First-stage checks only: 96 local M derivative and simultaneous-rotation identity
fixtures, seed20260928, epsilon=0,.4,2,10,Omega=1. Q,P normal scale.03,
B normal scale.1. Central derivative step1e-6,error<1e-8; torque identity
|Q cross M_Q+P cross M_P+B cross M_B|<1e-12; M>=m=1.
Then 48 n12,L8,radius1.2 grid interaction fixtures with random A scale.02,
Q,P scale.03,q in[-.5,.5]^3: vary A,q,Q,P jointly and compare M derivative
to the curl-deposited and positional/internal gradients, step1e-6,error<1e-8.
Check curl adjoint identity to1e-12. Preserve all fixtures and source hashes.

This does not yet implement a full evolution. Positive M is not stability or
causality: the new derivative-dependent, regularized coupling needs its own
full-system Hamiltonian, propagation and angular-budget checks before 3D runs.

Attribution: squared shifted momentum, spin/field torque, rotor exchange and
Noether accounting have substantial prior art. Primary abstracts reviewed:
Banerjee and Fehske, reciprocal spin-rotor dynamics (2026),
https://arxiv.org/abs/2604.23768v2 ; Turcati et al., spin-transfer torque with
non-minimal electromagnetic coupling (2020), https://arxiv.org/abs/1907.12196v2 .
These are related mechanisms, not this gravity model, and no originality claim
is made for rotor coupling or angular exchange. Our specific constitutive
choice and proposed gravity interpretation remain hypotheses. Old gravity
formula agreement is not a test gate. All twelve original goals remain active.
