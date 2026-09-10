# A local energy current for the ideal internal reservoir

**The existing ideal packet Hamiltonian admits a locally conserved energy current that carries acquired reservoir energy with the packet.** It need not leave that energy stationary in the void. This is a classical kinetic extension of the stipulated packet construction, not a derived electromagnetic/graviton field theory or an observational validation.

The numerical comparison also demonstrates why balanced global energy accounts are insufficient: a traveling reservoir and a stationary reservoir can receive the same total energy yet leave different spatial distributions behind.

## Conditional derivation using known Hamiltonian and kinetic mathematics

Retain the hypothetical effective Hamiltonian in reference units c0=1:

`H(x,p,s,P) = p/(1+a*s) + R(P)`.

On the ideal positive-transfer branch, `R(P)=E0+P`, `p>0`, `P>=0`, and `1+a*s>0`. E0 is a pre-existing reservoir seed energy; it is included rather than counted as free. The earlier positive completion outside this branch is not altered here. This Hamiltonian is a postulate using established extended-state mechanics, not a novel fundamental law.

Known Hamilton equations give `xdot=v=1/(1+a*s)`, `pdot=0`, `sdot=R'(P)` and `Pdot=a*p/(1+a*s)^2`.

For a nonnegative ensemble density f(x,p,s,P,t), the corresponding known Liouville equation is:

`partial_t f + partial_x(v*f) + partial_s(R'(P)*f) + partial_P[a*p*f/(1+a*s)^2] = 0`.

Integrate over the internal and momentum variables, with no unaccounted flux through their boundaries. On the declared positive-transfer trajectories, these support conditions can be imposed. There is no creation or capture term in this transport calculation.

`u_total = integral H*f dp ds dP`
`J_total = integral v*H*f dp ds dP`
`partial_t u_total + partial_x J_total = 0`.

The identity follows because H is constant along each Hamiltonian trajectory. Splitting propagation energy e=p/(1+a*s) from reservoir energy R gives:

`partial_t u_prop + partial_x J_prop = -Q`
`partial_t u_res + partial_x J_res = +Q`
`Q = integral [a*p*R'(P)/(1+a*s)^2]*f dp ds dP`.

Each current contains the same packet velocity v multiplying that component's energy. The transfer is local and opposite in the two equations. The component named "photon" in the diagnostic script is this stipulated propagation-energy term; no Maxwell-field identity has been established for it.

## Synchronized ideal branch and executed test

Choose all initial internal coordinates s=0, internal momenta P=0, p=1 and E0=0.1. Then s=t, n=1+a*t, and all packets move through `X(t)=ln(1+a*t)/a`. For an initial number-density profile N0:

`N(x,t)=N0[x-X(t)]`
`u_prop=N/n`, `u_seed=E0*N`, `u_gain=(1-1/n)*N`
`J_total=v*(1+E0)*N`.

These are exact conditional solutions, not a new redshift law. The initial profile is compact, has unit integral and occupies 1<x<2. The diagnostic uses a=0.1, final t=8 and domain 0<x<12. All quantities are dimensionless; no astronomical rate, age or dataset is calibrated.

A conservative finite-volume calculation evolves the propagation term, seed, moving gain and a stationary-gain control. In the control, transferred energy stays at the transfer location instead of having a reservoir flux. Its seed still moves, making both ledgers include the same initial energy. This control changes the transport model; it is not another solution of the same packet Hamiltonian.

| Quantity at final time | Moving gain | Stationary-gain control |
|---|---:|---:|
| Initial total energy, including seed | 1.100000 | 1.100000 |
| Acquired reservoir energy, entire domain | 0.444444 | 0.444444 |
| Acquired energy remaining left of x=5 | 6.6e-39 | 0.295304 |
| Total energy transported through x=5 | 1.100000 | 0.804696 |

The exact moving packet ends at 6.87787<x<7.87787, entirely to the right of the monitored region. Its final propagation energy is 0.555556, acquired reservoir energy is 0.444444, and seed energy is 0.1. All three components travel; the seed has not disappeared from the energy budget.

The stationary control leaves about 0.2953 units behind x=5 even though it gains the same total reservoir energy. Independent integration over the initial packet gives 0.2951968 for that stationary remainder. Thus global conservation alone cannot choose the user's intended no-void-deposit behavior; the local current is essential.

## Numerical accuracy

| Cells | Number-profile L1 error | Stationary remainder error | Maximum local energy-balance error |
|---:|---:|---:|---:|
| 480 | 0.64125 | 0.000427633 | 1.33e-15 |
| 960 | 0.42335 | 0.000214016 | 2.89e-15 |
| 1920 | 0.25481 | 0.000107058 | 2e-15 |

The upwind scheme is deliberately simple and diffusive: the finest number-profile L1 error is about 0.255, so it is not a precision pulse-shape calculation. Profile and stationary-remainder errors decrease on refinement. Exact component fractions, integrated global balance and the balance in the region x<5 are independently checked. The no-trailing-deposit claim for the moving branch follows from its analytic compact solution, not from interpreting a tiny numerical residual as exact zero.

## What is still missing

This fills a local-energy-current requirement for a classical ensemble of the postulated packets. It does not make a local relativistic field theory. The canonical momentum density `pi=integral p*f` has its own conserved flux, but identifying it with total physical relativistic momentum is an additional step. For this monoenergetic branch `J_total/pi=(1+E0)/n`, which is not generally unity in c0=1 units. An ordinary symmetric flat-spacetime stress-energy interpretation cannot simply be assumed; any medium, preferred frame or changed geometry needs explicit accounting.

We still need a physical electromagnetic channel that changes frequency as the stipulated propagation term changes, preparation and synchronization of newly created packet states, a physical source for E0, ordinary endpoint-clock behavior, and a capture/gravity coupling. The existing supernova timing requirement is not re-tested here: resetting each independently born internal clock remains different from preparing a shared phase.

The spatial mechanism tested here is also distinct from the earlier optional baryon-attached shielding model. The latter remains a poor source-shape fit and has unresolved momentum and support requirements. Transporting energy successfully does not establish how it becomes the required deposited gravitational source.

No holdouts, empirical coefficients, source-energy budget or observational likelihoods are changed. The full theory remains incomplete. The next mechanism step should identify creation/EM coupling and its stress-energy, rather than counting this kinetic construction as a completed graviton theory.

## Reproduction

Run `run.py`, then `report.py`. All three resolutions, analytic controls, global/control-volume balances and source hashes are retained in `results.json`. The previous internal-reservoir code is unchanged.
