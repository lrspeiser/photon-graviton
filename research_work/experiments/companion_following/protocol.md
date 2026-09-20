# GF-1: companion following and collective turning
Declared before this campaign's calculations. Repository baseline main7762282;
reviewed whirlpool branch claude/tests-clusters-lensing-xrhttm at2be9b33.
This is a separate experiment, not an amendment to RW-1 or RUT-10.
The owner asks for hundreds of formulas for moving companions that follow the
turn of a companion ahead, with the response extending outward.

## Scope and claims
Produce a searchable catalogue of600 structurally specified force/memory laws
and18 conservative kinetic-coupling parameter cases. They are618 candidate
formulas, not618 independently invented theories. Test actual trajectories.
No dark matter, expanding background, fitted astronomical distance, inserted
galactic halo or claim of an observed cluster solution. Finite initial packet
energies are counted. Formation from photon conversion is not derived in this
isolated interaction screen, and no supplied population is called a derived
astrophysical reservoir. "Graviton" here means the owner's hypothetical companion;
spin2, microscopic dispersion and gravitational universality remain unresolved.

The existing branch is reviewed as committed evidence, not silently merged.
RW-1 has a protocol/library/driver but no committed result archive at the reviewed
tip. CL-2/NL-1/NK-1/TF-1 reports are not new calculations of this campaign.

## Definitions shared by the600 candidates
Planar packet coordinates x_i, velocities v_i, unit inertial weights.
d_ij=x_j-x_i, r_ij=|d_ij|, e_ij=d_ij/r_ij, n_i=v_i/|v_i|.
J rotates a vector90 degrees counterclockwise; cross is the planar signed cross.
Self pairs have weight0. At exactly zero speed the transverse projector is0.
No speed floor, clipping, friction pump or uncounted stabilizer is used.
A common central attraction/repulsion is the established generalized Morse form:
U=(1/16) sum_{i<j} [2 exp(-r_ij/0.4)-exp(-r_ij/2)].
Its analytic negative gradient is included in every force. The finite short-range
repulsion is a declared regularization, not a point-collapse singularity.
The common coefficient1/16 stays fixed when packet count changes.

Six dimensionless range functions, s=r/2.5:
K0=exp(-s^2); K1=exp(-s); K2=(1+s^2)^(-1);
K3=(1+s^2)^(-2); K4=max(1-s,0)^4(1+4s);
K5=s^2 exp(1-s^2).
Five angular gates, c=n_i dot e_ij, h=n_i dot n_j:
G0=1; G1=(1+c)/2; G2=max(c,0); G3=max(c,0)^2;
G4=max(c,0)(1+h)/2.
W_ij=K_a G_b, S_ij=(W_ij+W_ji)/2; no row normalization.
The positive-part operations define the declared gates/kernels, not output repair.

Four memories supply u_i (a velocity signal):
M0: u_i=v_i.
M1: tau z1dot=v_i-z1; u_i=z1.
M2: tau z1dot=v_i-z1; tau z2dot=z1-z2; u_i=z2.
M3: tau^2 z2ddot+2 zeta tau z2dot+z2=v_i,
    u_i=z2, with zeta=0.35.
All memory states begin at the initial velocity, derivatives0; tau=0.5.
Memory acts on the source's velocity, not on a newly prescribed field.
M3 is an oscillatory causal-in-time filter, not a finite-speed spatial mediator.
All600 have instantaneous spatial interactions; this limitation stays explicit.

Five force operators, acceleration a_i=-grad_i U+0.2 Phi_i:
F0 velocity following: Phi_i=sum_j W_ij(u_j-v_i).
F1 transverse following: Phi_i=P_i sum_j W_ij(u_j-v_i),
   P_i=I-n_i n_i^T.
F2 curvature following: Phi_i=P_i sum_j W_ij(u_j-v_i+tau omega_j J u_j),
   tau omegadot_i=cross(v_i,a_i)/|v_i|^2-omega_i.
   omega_i initially0 in the chain and the known initial circular angular rate
   in the ring. At zero velocity its measured angular rate is0.
F3 reciprocal gyroscopic exchange:
   Phi_i=sum_j S_ij cross(u_i-u_j,e_ij) J(v_i-v_j).
   Pair forces are opposite and do no pair work; angular momentum must still
   be tested rather than declared conserved.
F4 attraction to a trailing target:
   q_ij=d_ij-tau u_j;
   Phi_i=sum_j W_ij q_ij/(1+|q_ij|^2)^(3/2).
   This is an approximate moving target, not an evolved retarded field.
The grammar gives5 x6 x5 x4=600 explicit equations. Every combination is saved.

## Conservative alternative:18 cases
Six K above at eta=0.5,2,8, with no directional gate or auxiliary memory.
A benchmark Lagrangian, using established position-dependent kinetic geometry:
L=1/2 sum_i |v_i|^2 + eta/32 sum_{i<j} K(r_ij/2.5)|v_i-v_j|^2-U.
Let k=eta/16 and let d_i_j=x_i-x_j in the following derivation.
M=I+k Laplacian(K), p_i=v_i+k sum_j K(v_i-v_j).
The Euler-Lagrange system is
M a = -grad U + k sum_j K'_r [
          0.5 e_i_j |v_i-v_j|^2
          -(e_i_j dot(v_i-v_j))(v_i-v_j)].
H=1/2 sum |v_i|^2 + k/2 sum_{i<j} K |v_i-v_j|^2+U.
Total momentum=sum p_i=sum v_i; angular momentum=sum cross(x_i,p_i).
M is positive definite for eta>=0 and K>=0. Verify the derivation by independent
finite differences of L, numerical energy/momentum/angular momentum, and a
two-body constant-K analytic acceleration-sharing control.
This is reciprocal acceleration entrainment. It does not postulate a preferred
leader and does not supply finite-speed relativistic propagation.
Negative eta is an explicit positivity control, not accepted if M has a
nonpositive eigenvalue. No ghost is used to manufacture a successful orbit.

## Declaration of fixtures and interpretation
A. Turn chain:12 packets at(-i,0), v=(1,0), evolved through time32.
The first packet's prescribed speed1 direction is
theta0=0.4 sin^2(pi(t-3)/8) for3<t<11, and0 otherwise.
Its acceleration and external reaction are recorded. This is an input-response
experiment, not spontaneous angular momentum creation. Compare follower angular
signals, lagged correlation, response versus distance, pulse amplification and
closest approach. The no-turn straight-chain symmetry is an analytic control,
and a numerical no-turn control is run for representatives.
A response is scored when the follower exceeds0.02 radians; report earliest
arrival, strongest angle, lagged correlation and last-four-packet mean response.
A behavioral chain pass requires last-four gain>=0.10 relative to the0.4 input,
best lagged correlation>=0.7, no follower gain>3, and min pair distance>0.1.
The signal need not survive arbitrarily far:24-packet extension is a separate test.

B. Prepared whirlpool:16 equally spaced packets at radius3, perturbed radially
by0.02 cos(2theta). Initial tangential speed is the circular equilibrium speed
of the unperturbed central potential; the conservative family uses its own
analytic circular balance including kinetic coupling. If no positive balance
exists, record it as unavailable, not as a silently changed fixture.
Evolve through time64; save initial period and number of periods covered.
No central mass, pinning field, applied torque or external drive.
Score circulation, radius spread, min separation, energy, momentum and angular
momentum drift. Behavioral ring pass: circulation>=0.8, mean radius0.7..1.3 of
initial, relative radial spread<0.25, min separation>0.1.
This tests persistence of prepared rotation, not whirlpool formation from rest.

Numerics: explicit fourth-order Runge-Kutta; default dt0.04, exact end time;
record every0.2. The600 candidates run both fixtures. Save every score and
exceptions/nonfinite states; never discard a failed formula.
Numerical laws with the derived conservative energy must have maximum relative
energy drift<1e-3, momentum drift<1e-5 (normalized initial sum speeds), angular
momentum drift<1e-3 (normalized initial canonical angular scale) in the closed
ring. F1/F2 must separately verify zero instantaneous guide work. F3 must verify
zero net guide force and guide work. F0/F4 work is recorded; integrating a work
residual is not an independent physical energy closure.
Failed laws can remain useful as phenomenological controls, not full successes.

## Refinement and robustness, selected by a fixed rule
Select up to12 highest behavioral candidates from the600, sorted first by the
number of chain/ring behavioral gates passed, then by minimum of chain gain
(capped at1) and circulation, then catalogue ID. This is exposed-fixture
selection, not a holdout. Keep every result, including ties resolved by ID.
For those12: dt0.02 both fixtures; turn amplitudes0.2 and0.6; reversed turn;
24-packet chain; coupling sign reversed. These are follow-ups, never retuned
original outcomes. Report score changes and small-signal odd symmetry.
For the conservative18: run derivative controls and both fixtures; refine the
top3 by the same ordering, and extend their ring integration to time256.
Timestep convergence: dominant response and radius metrics change<5% of max(1,
their magnitudes); conservation thresholds remain unchanged.
Count finite runtime as its actual duration, never a claim of permanent stability.

## Prior work and benchmark checks
Credit velocity alignment (Cucker-Smale), fixed-speed alignment (Vicsek), turning
inertia/spin (Cavagna et al.), and attraction/repulsion swarms (D'Orsogna et al.).
These precedents are mathematical analogies; they do not establish companion
or graviton physics. Do not rename published theories as ours.
Before the search: reproduce a symmetric two-agent alignment exponential and
the linear periodic turning-wave dispersion from a positive spin/spring chain.
These are narrow equation benchmarks, not reproduction of full biological data.
Verify the common potential gradient and the conservative Euler-Lagrange forces
independently. Compare transformations (translation, rotation, reflection),
zero coupling and zero input. No new observational score is invented.

## Evidence and check-ins
Commit protocol and catalogue specification before execution. Keep immutable
evidence directories and exact code/input hashes, runtimes and dependency versions.
Push main at declaration, after screening, and after verification/report.
Fetch before pushes; never force-push. Preserve all RUT, RW, CL, NL, NK and TF
evidence and the running feature branch. Inspect the branch again at final
report time; if RW results appear, report them separately.
No full historical suite with excluded expansion/halo comparison jobs is run.
