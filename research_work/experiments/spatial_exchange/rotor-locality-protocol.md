# SE-Rot4: finite-kernel locality check

Declared 20 September 2026. SE-Rot3's positive frozen Hessian does not prove a
finite propagation front. The term D* T D samples curl A throughout a finite
kernel before redepositing a source throughout that same kernel. Test whether
separated grid cells couple directly in the instantaneous RHS.

Use one normalized spherical kernel at q=0,radius1.2,L8, n24 and48. Internal
Q=sqrt(.001)e_x, T=|Q|^2I-QQ^T; epsilon=0,2,10. Put a unit perturbation in A_z
at x=+2/3,y=z=0. Measure the A_z stiffness response at x=-2/3,y=z=0, separation
4/3 model units, with all other perturbation cells zero. The original nearest-
neighbor Laplacian has zero matrix element for this separated pair at both
resolutions. The response is a Hessian matrix element per unit nodal amplitude;
it need not have the same magnitude across resolutions.

Compute the new coupling by actual centered curls and by the independent
factorized expression +epsilon^2 |Q|^2 (D_x W at input)(D_x W at output)/dV.
Require agreement to1e-12 absolute plus1e-10 relative. Require zero coupling
when epsilon0 and proportionality epsilon^2 for the nonzero cases to1e-10.
Use both directions to check reciprocity to1e-12.

Report nonzero remote response as a locality limitation, not a numerical test
failure to be hidden. Normalize additionally by cell volume to inspect whether
the kernel response persists as h decreases. A nonzero entry means these
fixed-width regularized equations do not have a strict pointwise local
continuum domain of dependence inside a source: distant cells can influence
accelerations at the same time. This does not itself give a measured signal
velocity or invalidate every effective finite-size model.

Do not use this to reject new gravity for disagreeing with older gravity.
It tests the claimed scope of our own finite-front interpretation. A causal
local completion would need internal source degrees of freedom distributed
through matter with their own propagation equations, rather than one instantly
shared rotor over the whole kernel. That completion remains undeveloped.

Pre-execution algebra correction: the factorized sign is positive. Summing
-W D_x delta(A_z) gives +D_x W at the input; the output curl contributes
+D_x W. The two spatial derivatives have opposite signs at the chosen points,
so the evaluated off-diagonal matrix element is negative. This corrects the
written formula before any numerical execution; the implemented model is unchanged.
