# SE-C: distinguish scalar enhancement from circulation

Declared 20 September 2026 before evaluation. For every completed SE-1/SE-2
final state, measure circulation of the actual vector shift beta, not scalar
source strength. Also decompose angular momentum between matter and fields.

At radii 1,2,3,4 in the source's rotated equatorial plane, compute the line
integral integral beta dot dl. Use trilinear interpolation of beta at 128,256
and 512 equally spaced angles. Test absolute difference between the last two
<1e-9+0.005 abs(the 512 result); report failures without hiding small values.
This is a quadrature check, not a grid convergence claim. No static force law
is inferred from this circulation. Compare equal-radius values across cases.

For F=(phi,A,X,Y), J=sum Pi grad_centered F, field angular momentum is
-integral x cross J plus integral A cross Pi_A. Separate the orbital X/Y
part from other fields and spin, and add sum_i q_i cross p_i for matter.
Require reconstructed total to match the archived final angular vector within
1e-10. Report projections along the initial source axis and the X/Y share of
initial total angular momentum (zero denominators reported undefined).

Report all cases including no-excitation and sign mirror. Use exact documented
rotation matrix, no fitted axis. Validate interpolation/integration against
analytic beta=(-b y,b x,0), b=.001, whose loop integral is 2 pi b r^2, error
<1e-12. This identity is vector calculus, not a gravity benchmark. Save source
and state hashes, all radii and quadrature results, and explicit gate statuses.

Interpret only the short generated state at T4. A larger scalar source is not
necessarily larger circulation. A transient loop integral does not prove a
stable vortex, force enhancement, continued outward transport or physical
graviton attraction. Full spatial/time convergence and observations remain.
