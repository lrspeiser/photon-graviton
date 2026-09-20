# SE-LR3: separate circulation from angular-momentum transfer

Declared 20 September 2026 before execution. SE-LR2 produces circulation with
negligible A-sector angular momentum. A candidate explanation is that its
leading source curl(Q cross P) is an axisymmetric toroidal source when the
free internal spin density is radial and points along z. For an axisymmetric
vector pattern, partial_phi A=e_z cross A, so

    J_A,z = integral Pi dot (e_z cross A - partial_phi A) = 0.

Thus circulation of an axisymmetric vector pattern need not transport spin.
This is a continuum identity, not proof that all finite-coupling LR2 fields
are exactly axisymmetric. Higher-order feedback can break that restriction.

Test a local extension K=P-lambda A-epsilon(curl A) cross Q. Keep every other
term in the positive SE-LR1 Hamiltonian unchanged. The canonical equations are

    A_t=Pi
    Pi_t=Delta A-omega_A^2 A+lambda K+epsilon curl(Q cross K)
    Q_t=K
    P_t=c_Q^2 Delta Q-Omega^2 Q-epsilon(curl A) cross K.

The direct lambda term can excite a rotating vector pattern instead of only
the curl of the spin density. It has no spatial derivatives and leaves the
principal matrix/symmetrizer derived in SE-LR1 unchanged. This does not prove
global stability. Canonical angular momentum still uses P, not K or P-lambda A.
The energy remains a sum of positive squares. Lambda carries the units needed
to match P/A; numerical values below are in model units.

The momentum-shift construction and canonical vector angular momentum are
established mathematics. We do not claim they or this coupling are historically
new. This is a proposed ingredient in our model, not a copied gravity law used
as a success target. Ordinary-matter generation, physical fuel, light and the
scalar sector are still absent. No dark matter or expanding background.

First check48 random full Hamiltonian directional derivatives: n8/12, L8,
epsilon0/2, lambda0/1, six fixtures each, seed20261003, amplitudes.03 and unit
Euclidean perturbation directions. Steps1e-6 and5e-7; scaled errors and their
difference must be<2e-7. Use unchanged local discrete energy/gradient rules.

Then run five cases with the identical compact LR2 rotating excitation, T1,
L8,c_Q.5,Omega1,omega_A.2: (epsilon,lambda)=(2,0),(0,1),(2,1) at n24dt.02;
combined time control n24dt.01; combined space control n32dt.01. Archive
initial/final fields and every-step total/sector energy, total angular and
linear momentum, edge energy and separate A-sector angular momentum.

Retain LR2 individual gates: energy drift<1e-5, total relative angular drift
<1%, edge-energy fraction<1e-5, and128/256 circulation quadrature error<1e-5
at radii1,1.5,2. Require lambda0 final fields reproduce archived LR2 rotating
run to1e-12. For the combined case time comparison require relative r1.5
circulation and J_A,z difference<.001; space comparison require both<.05,
denominators max(abs(finer),1e-10). If a quantity is effectively zero, report
it explicitly rather than interpret relative accuracy as signal detection.

Report actual transfer fractions J_A,z/J_initial,z and E_A/E_initial. No
minimum transfer is stipulated as a numerical pass criterion. Positive transfer
would establish only this finite-excitation toy mechanism, not persistent
swirl or useful gravity. Keep failures and thresholds unchanged.
