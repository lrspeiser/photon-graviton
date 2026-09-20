# CC-1: shared local propagation geometry

Declared 20 September 2026 before numerical implementation. Parent evidence: FM-1 at 1fab2fc. This is the next equation/3D-source audit, not the completion of the twelve-item research goal.

## Proposed change

Use the same lapse alpha and drift beta in the particle Hamiltonian and field wave operator. A coordinate speed above one is not itself a causal violation when the local mediator cone changes with it. Test equality of the local photon and mediator characteristic cones instead of clipping photon velocities after integration.

Three-dimensional fields F=(phi,A_x,A_y,A_z), canonical Pi, and Euclidean spatial gradients:

alpha=exp(g phi); b=kappa eta A/sqrt(1+eta^2 |A|^2); beta=alpha b, with 0<=kappa<1.

K=(sum Pi_a^2+sum |grad F_a|^2)/2; J=sum Pi_a grad F_a.
H_field=integral [alpha K-beta dot J+omega^2 |F|^2/2] d^3x.
H_particle=sum [alpha(q_i) sqrt(m_i^2+p_i^2)+beta(q_i) dot p_i].

The positive, constant-coefficient mass potential is outside alpha. This is a declared choice, not the covariant scalar action or Einstein gravity. Lapse/shift optical Hamiltonians are established mathematics. The dynamical constitutive choice for alpha,beta is the project candidate; historical originality is not claimed.

Canonical equations, with comma-a denoting derivative with respect to F_a:

qdot=alpha p/E+beta; pdot=-E grad alpha-(grad beta)^T p.
Fdot_a=alpha Pi_a-beta dot grad F_a.
Pidot_a=div(alpha grad F_a-beta Pi_a)-omega^2 F_a
         -alpha_,a K+beta_,a dot J
         -sum_i [E_i alpha_,a+p_i dot beta_,a] delta(x-q_i).

All field self-source terms must be retained in any future solver. These expressions describe the point-source continuum formalism; point-source energy regularization and finite source implementation remain open. The spherical kernel below is a tested regulator, not a proof of point locality.

## Analytic checks to document

1. |b|<kappa implies positive field kinetic/gradient energy and H_i>=alpha(E-kappa|p|)>0 for nonzero energy. Particle Hessian is alpha(I/E-pp^T/E^3): positive for m>0, semidefinite for photons. No global nonlinear stability conclusion follows.
2. Frozen principal symbol gives frequencies beta dot k +/- alpha |k|. Photon velocity is beta+alpha n; massive velocity is inside this local cone. The field Hamiltonian needs the MINUS beta dot J sign for this match; retain the opposite sign as a negative control.
3. About source-free zero-field vacuum, all four modes have omega_mode^2=|k|^2+omega^2. For nonzero homogeneous frozen coefficients the dispersion is (frequency-beta dot k)^2=alpha^2 |k|^2+alpha omega^2. Such nonzero backgrounds need not be solutions: this is a coefficient/principal audit, not their nonlinear stability.
4. For an energy-homogeneous massless Hamiltonian H=|p|h(n,F), if its vacuum velocity saturates a fixed unit speed bound and fields vary smoothly with both signs, differentiating the bound implies the first-order field derivative of h vanishes at vacuum. State all assumptions. This explains why smooth fixed-speed clipping tends to remove the desired linear light coupling; it is not a no-go theorem for all modified gravity.
5. Weak-field source is (g E,kappa eta p), preserving a leading-order current-current term on static elimination. Beyond that order the field self-coupling is essential.

## Declared numerical audit

Seed 20260920. Use 600 three-dimensional local samples: g and eta each drawn from {-0.3,-0.08,0,0.08,0.3}; kappa from {0,0.25,0.5,0.9}; phi uniform [-4,4], each A component uniform [-4,4]; m cycles 0,0.01,1; p and wavevector normal random (nonzero). These are equation stress tests, not fitted parameter variants. omega=0.2.

Check analytic particle velocity, field-source gradient and momentum Hessian against centered finite differences: scaled error <=2e-5. Check photon cone residual <=1e-12, massive cone excess <=1e-12, positive field kinetic block and particle energy, and frozen Fourier eigenvalue error <=1e-10. Check rotation covariance of alpha,beta,H,qdot using seeded orthogonal proper rotations, <=1e-11. Rescale photon momentum by 2 and require unchanged direction/velocity <=1e-12. Record negative-control opposite-sign cone mismatches without counting them as model successes. Do not silently omit extreme or failing samples.

Spherical 3D profile W proportional to max(1-|x-q|^2/a^2,0)^3, normalized on the grid. Differentiate the discrete normalization. Check analytic continuous normalization via radial quadrature and first/second moments; grid sum and gradient sum; adjoint deposition identity; finite-difference sampled-field force. Use domains [-4,4)^3, N=32,48,64, a=0.6,0.9,1.2; 12 fixed rotations about a non-axis-aligned unit axis. q=R(0.3,0.2,0.1), sample function exp(-|x-R(0.7,-0.2,0.4)|^2/2). Rotation spread relative to the zero-rotation sample <=0.02 at N=64; record all grid errors. Exact continuous profile rotation covariance <=1e-12; force scaled error <=2e-5; normalized weight and derivative sums <=1e-12. No smoothing radius is chosen from its best result. Changing radius changes the regulated model; report that sensitivity instead of demanding it vanish.

Before execution, commit the implementation. Preserve first results and source hashes in a new evidence-v1 directory. Document all failures and keep the twelve-goal ledger incomplete. No dark matter, expansion, astronomical data fit, or claim of full 3D nonlinear evolution occurs in CC-1.

## Attribution

Gibbons, Herdeiro, Warnick and Werner, Stationary Metrics and Optical Zermelo-Randers-Finsler Geometry: https://arxiv.org/abs/0811.2877. Hamiltonian geometry and stationary optical drift are precedents, not project inventions. Our dynamical scalar/vector wave action is not derived from that paper and does not inherit Einstein constraints or observational validation.

For later local comparisons, the primary GW170817/GRB analysis constrains propagation differences under its emission-time and source assumptions: https://dcc-lho.ligo.org/LIGO-P1700308/public. CC-1 does not import a cosmological distance model or claim to reproduce the event.
