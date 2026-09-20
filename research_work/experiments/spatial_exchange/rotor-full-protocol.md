# SE-Rot2: embed the rotor in the full discrete Hamiltonian

Declared before implementation, 20 September 2026. Retain the six-field SE-2
Hamiltonian and replace each scalar internal emitter by SE-Rot's vector Q,P.
There is no simultaneous direct X/Y emitter in this variant. Keep all existing
field coefficient derivatives and particle sources, using the new internal M.
Add curl-deposited source, positional reaction and canonical rotor spin exactly
as derived in rotor-protocol.md. Preserve previous numerical sources unchanged.

State layout: six F and Pi fields; six q,p,Q,P vectors; sponge energy ledger.
In empty fields, circular internal states Q=sqrt(Eint) e_x,P=spin sqrt(Eint) e_y
for spin=+/-1 have energy Eint and spin=+/-Eint along z (Omega=1). A spin=0
control uses Q=sqrt(2Eint)e_x,P=0, same energy. Rotate internal and orbital
coordinates together. Orbital speed becomes an explicit parameter, default.2.

First full-system controls: 48 fixtures seed20260929, n12,L10,radius1.2,
rotor coupling0,.4,2,10,chi0/200, both internal-spin signs, compact random
F/Pi scale.03/.02 and perturb q,Q,P by scale.03. Compare complete energy
directional derivative with complete RHS after removing sponge, step2e-5,
scaled error<2e-6. Require M>0 and matter probe-cone excess<1e-10.

For the twelve zero-coupling fixtures, set Q,P along x, compare field/q/p/Qx/Px
RHS and total energy with the pinned SE-1 scalar oscillator at emission=0,
error<1e-12; orthogonal rotor derivatives must vanish to1e-12. This is a
reduction check of the same Hamiltonian, not a target from old gravity.

Do not launch physical claims or a new production evolution from these checks.
The derivative-dependent regularized coupling requires source-region linear
response/propagation tests, timestep bounds and angular-momentum convergence.
Matter retaining its original averaged speed bound does not prove that the
modified field has the same propagation behavior. No full-system stability
claim follows from gradient agreement or positive energy.

Attribution and finite angular-reservoir assumptions remain those explicitly
recorded in rotor-protocol.md. The galaxy/cluster and all twelve original
requirements remain unfulfilled and active.
