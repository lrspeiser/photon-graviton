# Capture and restoration need a temporal-state energy budget

We added an explicit energy U(n) for the accumulated temporal state, rather than treating that state as free. This is a postulate-based constitutive model: it does not explain what time ultimately is and does not need to do so. It changes the feedback and tests whether capture can restore n while leaving positive deposited energy.

## Derivation

In the existing homogeneous reference-time wave law, radiation energy obeys R_dot=-(n_dot/n)R, so R*n is constant. This is known integration conditional on our chosen Hamiltonian. It is reference radiation energy: a local observer with d tau=dt/n assigns energy proportional to n*R for a fixed set of quanta. Reference blueshifting is therefore not automatically an observable blueshift. A full spatial ray/clock comparison is still necessary.

Postulate either U(n)=K(n-1) or U(n)=K(n-1)^2/2 for n>=1. Neither shape or K is an established temporal-energy law or a claimed unique discovery. The point is to specify a cost and follow its consequences.

For a local passive reset without incoming energy, external work or another reservoir,

    D(n) = U(n_i)-U(n) - [R_i*n_i/n-R_i]
    D_dot = (-n_dot) [U'(n)-R(n)/n].

These are conservation identities conditional on that postulate. Restoring n releases stored state energy, but part may be needed by the radiation sector in the chosen reference ledger. We cannot assign the entire released amount to the deposit as well. If D_dot<0, this proposed capture-only operation needs another paying sector; it is not a permitted negative-energy deposit.

For linear U, a complete reset's endpoint budget is (n_i-1)(K-R_i). A nonnegative transfer into deposits throughout the prescribed monotone reset additionally requires K>=R_i*n_i. Positive final budget alone is weaker than a physically allowed positive capture rate at every moment.

For quadratic U, U'(n) approaches zero near n=1, while R/n remains positive. Thus a passive full reset with radiation still in the same closed cell cannot keep the deposit rate nonnegative all the way under these assumptions. This does not rule out export, energy supplied by matter, or different constitutive postulates.

## Numerical checks

Eighteen reset combinations use n_i=1.1, K=0.05,0.5,2 and R_i=0.1,1,10. Three pass nonnegative sampled deposit rates: linear U with (K,R_i)=(0.5,0.1),(2,0.1),(2,1). Their full-reset budgets are 0.04,0.19,0.10 in arbitrary energy units. All other choices fail the passive-rate test. Forced failing paths remain diagnostic results, not adopted physical evolutions. Integrated budgets agree with the exact formula within 3.56e-13.

Adding memory cost also changes the previous growth rule. With n_dot=T, linear U and capture Gamma, the receiving field must obey

    T_dot = T [R/n-K-Gamma],

rather than counting all radiation loss as receiving energy while also creating U for free. Six growth runs account for radiation, receiving energy, memory and deposits, conserving their reference total to 7.55e-15. Initial radiation is 1 and seed 0.01. At K=0.5, the initial growth coefficient is +0.5 without capture and -0.5 with Gamma=1. At K=2 it is already -1 without capture. The price of creating the state therefore changes whether feedback grows; these are illustrative parameters, not a fitted theory.

## Closed-cycle consequence

A further exact check follows the no-capture growth runs and considers restoration to n=1 while keeping the same radiation in the cell. R*n conservation restores the original reference radiation energy of 1. The net endpoint deposit budget is only initial seed minus remaining receiving energy. It is not a net conversion of retained radiation into gravity.

For K=0.5 the reset budget is approximately 0.00923, equal to seed 0.01 minus remaining receiving energy 0.000768. For K=0.05 the endpoint budget is negative at the tested time because substantial energy remains in the receiving sector. This is not a perpetual energy source; every sector has to be tracked. Endpoint cycle budgets do not separately certify positive instantaneous capture rates.

## Implication for the selected hypothesis

The user proposes an open process: radiation travels onward and wells later capture temporal-field/companion energy. That is not the same as the closed cell just tested. Spatial separation or different loading/capture epochs could change the budget, so this calculation does not reject the overall idea. It does establish that restoring clock state is not an automatic consequence of removing receiving energy and cannot be added for free.

The next calculation should carry the new U(n) cost into open spatial transport and derive a capture/restoration response from available energy, including possible energy transfer to ordinary matter. Test the already problematic endpoint clocks and both messenger maps; do not reset observers by hand or interpret reference-energy loss as measured redshift. Spatial field forces and momentum remain outside the present scalar budget and must be added before a full-conservation claim.

Reproduce with `python research_work/results/temporal-memory-energy/run.py`. The reset tests, growth cases and three exact closed-cycle checks are in results.json. No redshift coefficient, observation, source-energy supply conclusion or holdout status changed. This advances the operational postulates and energy/feedback goals; all nine goals remain incomplete.
