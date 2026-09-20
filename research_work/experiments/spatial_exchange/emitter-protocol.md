# SE-2: emission into a fixed combination of the two channels

Declared 20 September 2026 before implementing or calculating this variant.
SE-1S found no enhancement of integrated wave-sector scalar source in its
completed examples. Test whether its X-only emitter is limiting the response.

Retain SE-1's field Hamiltonian and initial empty fields. Replace Xbar_i in the
internal mass by Bbar_i = cos(theta) Xbar_i + sin(theta) Ybar_i, using one
constant theta for every source. Internal mass is

    M_i = m_i + [P_i^2 + Omega^2 (Q_i - emission Bbar_i)^2]/2.

Writing d_i=Q_i-emission Bbar_i and R_i=sum W alpha M_i/E_i, the field
momentum sources are R_i Omega^2 emission d_i W/dV times cos(theta),sin(theta)
for X,Y respectively. The additional positional force uses grad_q Bbar_i.
Qdot=R_i P_i and Pdot=-R_i Omega^2 d_i. The emitter direction is constant,
so there are no omitted field derivatives of theta. All other equations,
physical parameters and source smoothing remain SE-1's. No external field
reservoir is inserted. Theta is a new candidate parameter, not observed physics.

First verify 48 full Hamiltonian directional derivatives with seed 20260924,
theta cycling 0,pi/4,pi/2, chi cycling -50,0,50,200, both mixing signs and
emission=0,.4; use SE-1's compact random states and error gate 2e-6. Verify
theta=0 against the pinned SE-1 RHS and energy on all relevant controls to
1e-12 absolute. Retain frozen-cone tests; do not call these full stability.

Spatial runs, all with n32,L16,radius .9,T4,dt.02,initial excitation .001 per
source, same ring momentum .2: theta pi/4 and pi/2 crossed with chi 0,50,200
(six runs), plus theta pi/2 chi200 with zero excitation (one null run).
Run theta pi/2 chi200 at dt.01, at n40 dt.01, and rotated pi/3 (three more).
SE-1 theta0 archives provide the baseline, not a newly fitted target.

Use SE-1's gates: relative ledger drift<1e-5, cone error<1e-10, edge field
amplitude<1e-5, positive pre-boundary clearance. Zero excitation must produce
zero X/Y to 1e-12. Time comparison: source positions within .001 and each
X propagation, Y propagation and mixing energy within 1%, denominator floor
1e-12. Space and rotation: positions .01 and energies 5%. Record all failures.
Archive initial/final states and .1-time traces. Independently recompute the
full energies and wave-sector scalar source before interpreting enhancement.

No pass on numerical gates counts as a galaxy/cluster solution. Large source
response still requires the correct radial field, stable orbits, full light
propagation, a physical emitter and sustainable energy budget. Established
Hamiltonian and numerical mathematics is credited; no physical novelty claim.
