# SE-DT1: companion generation no longer requires initial matter motion

20 September 2026. Protocolcec47f9; implementation37d98ce. All48 full
Hamiltonian derivative fixtures,192 frozen principal fixtures and18 initial
source-kick cases pass. Largest scaled energy-derivative error is8.26e-9;
largest kick residual is1.68e-7. First results and source hashes are preserved
in density-transfer-v1/results.json.

The changed kinetic momentum is

    K=P-lambda A-zeta grad phi-epsilon(curl A) cross Q.

The scalar equation receives the reciprocal term -zeta div(Q_t). It is
implemented by differentiating the original Hamiltonian at P-zeta grad phi
and adding this exact chain-rule reaction, not by appending an independent
force. At zeta0 both energy and RHS reduce to the earlier JT1 equations.
The positive-energy bound on beta is unchanged. The enlarged frozen principal
matrix passes positive-symmetrizer, real-speed and reference-photon coverage
checks. These are not a global nonlinear stability proof.

## The initial-generation difference

Starting with empty fields, six source particles give

    phi_tt(0)=-g rho_eff, A_tt(0)=-d j,
    Q_ttt(0)=lambda d j+zeta g grad rho_eff,
    d=kappa c_Q eta.

At rest, j vanishes but a nonuniform mass profile has a density gradient.
The previous coupling therefore has a zero leading companion seed while
the new term need not. The isolated-source diagnostic excludes the photon
so a photon current cannot masquerade as matter-at-rest generation.

| Source at rest | zeta | Maximum predicted component of Q_ttt, model units |
|---|---:|---:|
|n8|0|0|
|n8|0.5|0.0185744|
|n8|2|0.0742975|
|n12|0|0|
|n12|0.5|0.0448166|
|n12|2|0.1792665|

One-step RK4 measurements at1e-4 and5e-5 reproduce these seeds, as well as
the moving-source cases. The n8/n12 values are not a convergence claim: the
particle interpolation profile changes with the grid, and its density gradient
is more sharply resolved. This is a short-time derivative check, not a full
finite-duration production run. At rest the initial response can be radial;
it does not by itself demonstrate swirl or angular-momentum generation.

The change addresses the velocity-suppressed initial seed exposed by JE1.
It does not show stronger gravity. The reciprocal scalar reaction could
screen or reshape the field instead of enhancing attraction. The next test
must evolve matched slow-source cases, including coupling-sign controls,
track the complete energy budget, and compare actual particle/light motion
under refinement. No initial companion reservoir should be inserted.

Gradient/kinetic couplings and Hamiltonian chain rules are established
mathematical structures, not a claimed original gravity theory. No dark
matter, expansion or agreement with older gravity formulas is used. Physical
fuel, long-term persistence, joint stability and galaxy/cluster observations
remain unresolved.

Concurrent point-bundle campaign: five of seven cases now pass411 independent
archive checks. Its control n48/n64 bend difference is3.77243%, above the
declared1% limit, while its transport-matrix difference0.00013745 is below
0.001. Thus that available spatial comparison fails; the full campaign remains
unfinished with direct-companion n64 and timestep runs pending.
