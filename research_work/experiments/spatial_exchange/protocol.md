# SE-1: matter-funded spatial radiation/companion exchange

Declared 20 September 2026 before calculations. Owner's latest requirements: derive the invented interaction, judge it by its own consistency and observations, allow coherent coupling/release onto different paths, and do not demand agreement with an older gravity formula. SR-2 continues unchanged. SE-1 connects RC-1's local conversion concept to the existing common-cone source dynamics rather than inserting a desired orbit or lens profile.

## Degrees of freedom and a complete candidate Hamiltonian

Six real fields F=(phi,Ax,Ay,Az,X,Y) and conjugate momenta Pi occupy a 3D grid. X is a radiation-like scalar channel; Y is a companion-like channel. Their identification with Maxwell photons and spin-2 gravitons is **not derived**. Matter has positions q_i, momenta p_i and canonical internal oscillators (Q_i,P_i). All six fields and field momenta start at zero. Matter's internal oscillators supply a finite, counted emission reservoir.

Set U=g phi, alpha=exp(U), z=C=exp(2U), a=exp(4U), beta=C kappa eta A/sqrt(1+eta^2 A^2). Use g=.08,kappa=.5,eta=.08. Define k(U)=k0 tanh(U/U0), U0=.001, m2(U)=mu^2 exp(2 chi U), mu=.75. The field energy is

Hf = integral [a sum Pi^2/2 + sum |grad F|^2/2 - beta dot sum Pi grad F
              + omega^2 (phi^2+|A|^2)/2 + m2(U)(Y-k(U)X)^2/2] dV,

with omega=.2. The potential is nonnegative. Its local eigenchannels include a massless direction and a massive direction. At U=0, k=0 so the external X channel is massless and Y has mass coefficient mu; conversion occurs where the self-generated field changes the coupling. A reduced companion mass near a source is a possible retention mechanism, not proof of a bound state or long residence time. All field principal characteristics share beta +/- C; mixing is a lower-derivative term. No prescribed external gravitational potential, expansion or independent dark component is added.

Let W_i be the same normalized compact spherical source profile used in SR-1, Xbar_i=sum_cells W_i X. Define the positive effective particle rest energy

M_i = m_i + [P_i^2+Omega^2(Q_i-lambda Xbar_i)^2]/2,

Hp = sum_i sum_cells W_i [alpha sqrt(M_i^2+z |p_i|^2)+beta dot p_i].

Here m_i=1,Omega=1 and lambda=.4 unless switched off. The squared difference is a specific, declared emission/absorption interaction. It makes the oscillator and radiation source reciprocal and includes its interaction energy inside M_i. This is effective matter modeling, not derived atomic physics.

## Required derivatives

With E_i(x)=sqrt(M_i^2+z p_i^2), R_i=sum W_i alpha M_i/E_i and d_i=Q_i-lambda Xbar_i:

qdot_i=sum W_i[alpha z p_i/E_i+beta];

pdot_i=-sum (partial_q W_i)[alpha E_i+beta.p_i] + R_i Omega^2 lambda d_i sum (partial_q W_i) X;

Qdot_i=R_i P_i; Pdot_i=-R_i Omega^2 d_i;

Pidot_X receives +sum_i R_i Omega^2 lambda d_i W_i/dV.

All original scalar/vector source terms are retained using these dynamic M_i. Field kinetic/shift coefficient derivatives include **all six fields**, not only phi,A. For r=Y-kX and k_U=dk/dU, the mixing-potential derivatives are

partial_phi Vmix=g m2 [chi r^2-k_U X r], partial_X Vmix=-m2 k r, partial_Y Vmix=m2 r.

The first term in particle force differentiates W at fixed M, and the second differentiates Xbar inside M; neither may be omitted. Use the original positive nearest-neighbor gradient energy with centered shift coupling, and its exact discrete adjoint. The finite spherical source regulator remains nonlocal within its declared radius.

## Tests and run declaration

First pin code, then execute 48 seeded full-state directional Hamiltonian derivative tests (seed20260923), spanning chi={-50,0,50,200}, k0={-.5,.5}, lambda={0,.4}, n=12,L=10,radius=1.2. Random compact fields/momenta plus nonzero internal Q,P exercise every reciprocal derivative. Require scaled finite-difference error <2e-6. Require positive M_i and averaged massive-particle cone bounds. Test the constant-coefficient six-field principal symbol and positive potential eigenvalues analytically/numerically; this does not establish nonlinear stability.

Eleven declared 3D runs, in order, all T=4,n=32,L=16,dt=.02,radius=.9 unless changed. Six matter sources form the existing radius .7 rotating ring, initial tangential momentum/mass=.2. Initial internal energy epsilon=.001 per particle, Q_i=sqrt(2 epsilon)/Omega and P_i=0, all equal phase. This prepared internal excitation is counted, not claimed to be a reconstructed stellar emission history.

1. no-emission: lambda=0,chi=0,k0=.5.
2. no-excitation: epsilon=0,lambda=.4,chi=200,k0=.5.
3. emission-only: lambda=.4,chi=0,k0=0.
4. exchange-0: lambda=.4,chi=0,k0=.5.
5. exchange-50: lambda=.4,chi=50,k0=.5.
6. exchange-200: lambda=.4,chi=200,k0=.5.
7. exchange-negative: lambda=.4,chi=-50,k0=.5.
8. sign-mirror: lambda=.4,chi=200,k0=-.5; must reproduce case6 under Y -> -Y.
9. time refinement: case6,dt=.01.
10. spatial refinement: case6,n=40,dt=.01,L=16.
11. rotation: case6 rotated pi/3 about (1,2,3).

Record initial/final raw states and .1-time-unit metrics. Record source-H change, field energy, internal rest-energy component, separate X/Y propagation energy, mixing-potential energy, outlet, momentum and angular momentum. Do not assign the full mixing energy to Y or label its amplitude as a conversion energy fraction. Record squared-amplitude spatial centers and RMS radii for X and Y only when nonzero; differing distributions do not alone prove particle attachment/detachment. No observed orbit or lens profile enters.

Numerical gates: maximum sampled relative H+Qout drift <1e-5; averaged massive-particle cone excess <1e-10; finite state; sampled edge amplitude <1e-5; positive conservative travel clearance to damping entrance using largest sampled characteristic speed and largest source extent. No late-time absorber or strict discrete-causality claim. Null controls require X,Y amplitudes <1e-12; emission-only requires Y<1e-12. Sign mirror requires relative field-energy difference <1e-8 and matter final-state difference <1e-7. Time refinement requires particle position difference <1e-3 and X/Y/mixing energies to agree within1%; spatial and rotation comparisons require particle position difference <.01 and those energies within5%, with a floor1e-12 for vanishing energy. Record all failures; no tuning to make them pass.

Independently reconstruct total energy from raw states and verify source hashes, gate decisions and comparison metrics. Archive the first runs. A positive result establishes only a source-funded spatial exchange mechanism in this candidate. It cannot establish enough gravity, stable bound modes, realistic radiation/spin, long-lived swirl or joint galaxy/cluster agreement.

## Attribution and next physical tests

Canonical Hamiltonian evolution, coupled positive wave channels, internal oscillators and mode conversion are established mathematics; retain CWC-1/RC-1 attribution and prior Gertsenshtein context without importing that physical conversion law. The chosen squared emission coupling and U-dependent mass/mixing are explicit effective assumptions. No historical uniqueness is claimed. If they pass, test signal propagation, binding/residence and release, spectral/image quality, force on separate stellar/photon probes and generation at realistic source speeds/energy allowances before any shared observation fit. Older gravity formulas are not acceptance gates; conservation requirements apply to this declared Hamiltonian itself.
