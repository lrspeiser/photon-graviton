# SE-LF1: local energy flux identity

Declared 20 September 2026, before calculation. For the unchanged SE-LR3
Hamiltonian, set B=curl A, K=P-lambda A-epsilon B cross Q, S=Q cross K.
The local energy flux candidate is

    F_i = -Pi dot partial_i A - c_Q^2 K dot partial_i Q
          + epsilon [Pi cross (Q cross K)]_i.

Derivation: partial H/partial(partial_i A)=partial_i A+epsilon e_i cross S
and partial H/partial(partial_i Q)=c_Q^2 partial_i Q. The canonical equations
give H_t=div(Pi dot H_gradA + K dot H_gradQ), hence H_t+div F=0. The direct
lambda term contributes through K, not as an additional gradient flux.
This is canonical Hamiltonian calculus, not a novel conservation principle.

Check24 smooth periodic Fourier fixtures: n24/32,L8,epsilon0/2,lambda0/1,
three seeds per combination generated with RNG seed20261004; c_Q.5,Omega1,
omega_A.2. Use eight prescribed integer wavevectors of components at most2
and amplitudes.003 per sine/cosine coefficient for each component of A,Pi,Q,P.
The resulting polynomial products are resolved without spectral aliasing at
these grids. Compute RHS using Fourier derivatives, independently compute
local energy density at state plus/minus delta*RHS, and compare this time
derivative to -div F for delta1e-6 and5e-7. Require maximum absolute residual
divided by max(1,max(abs(div F))) <1e-7 at both steps. Also require absolute
integrated H_t<1e-7 in the periodic box.

This checks the continuum flux formula with resolved Fourier fields, not an
exact discrete local balance for the existing mixed finite-difference energy.
No production run or source is changed. A compatible discrete surface-flux
ledger is still required before flux claims from those runs. No matter/light
completion, physical source fuel or observed gravity is demonstrated here.
