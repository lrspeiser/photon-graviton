# CWC-1 equations and their status
The canonical variational method, Maxwell constitutive framework, scalar
gradient/mass energy and symmetric splitting are established mathematics.
This document derives their consequences for the declared *candidate*.
It does not attribute those structures to this project.

Set a=exp(-g phi). The continuum equations from the protocol's H are

A_t=a D,
D_t=div(a grad A),
phi_t=Pi,
Pi_t=v_phi^2 laplacian(phi)-m_phi^2 phi-lambda phi^3
     + g a (D^2+|grad A|^2)/2
     + sum b g I_j Omega_j delta(x-X_j),
Xdot_j=P_j/M_j,
Pdot_j=b g I_j Omega_j grad phi(X_j),
theta_dot_j=Omega_j=exp(-b g phi(X_j)),
Idot_j=0.

Pairing canonical derivatives gives dH/dt=0 for periodic boundaries. This is
an energy theorem of this Hamiltonian, not proof that it is physical gravity.
The continuum translation charge is integral(-D grad A-Pi grad phi)+sum P_j.
Its lattice analogue need not be exactly conserved and is tested separately.
Material clocks can donate energy and produce phi even without EM when b>0;
calling all receiving energy photon conversion would then be incorrect.

For homogeneous phi, electromagnetic oscillator actions are constant and
H_EM=A0 exp(-g phi). The scalar equation is
phi_tt = g A0 exp(-g phi)+b g C0 exp(-b g phi)
          -m^2 phi-lambda phi^3.
A ray travels dx/dt=a. For fixed path distance and two neighboring launch
times, coordinate arrival Jacobian J=a_emit/a_receive. Normalizing emission
and reception with the stipulated material clock gives
S_time=J Omega_receive/Omega_emit
      =exp[(1-b)g(phi_receive-phi_emit)]=S_spectral.
This exact identity of the model is inherited characteristic/clock algebra.
For an increasing index: b0 permits redshift, b1 cancels it, b2 reverses it.
The field may oscillate, so the numerical histories must determine the sign.

With fixed coordinate rods the modeled local optical speed in clock units is
a/Omega=exp[(b-1)g phi]. Exact equality to the reference for every phi needs
b1; that same homogeneous choice cancels measured redshift. These are
conditional limitations of the family, not a theorem against all time models.

At fixed phi and gradient the extra acceleration is
b g (I_j/M_j) exp(-b g phi) grad phi.
Two internal-action/mass ratios need not accelerate equally. b0 gives no
extra material force. This does not establish a universal gravity coupling.
Treating I as an independently fitted galaxy mass would violate our source
and attribution rules. Rest-mass/atomic binding/covariant stress completion
is missing; the model's oscillator energy is an explicitly effective sector.

The lattice places phi,A,D,Pi on nodes and uses positive edge-averaged a for
gradient(A)^2. Variations of this *same discrete H* determine both EM forces
and scalar backreaction. Material interpolation and deposition use the same
multilinear weights, with their actual spatial derivatives for recoil.
Symmetric V/2, K_EM/2, K_phi+K_matter, K_EM/2, V/2 subflows are individually
exact. Hamiltonian derivative, reversal and deliberately missing-reciprocity
checks verify the implementation independently.

No Maxwell-to-spin2 particle amplitude, general-relativistic field equation,
photon-number operator, physical atomic Hamiltonian or cosmic formation
history is derived here. Those are explicit promotion gates.
