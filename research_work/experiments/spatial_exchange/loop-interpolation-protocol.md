# SE-LI1: isolate circulation interpolation error

Declared 20 September 2026 before execution. LR4 fails its2% rotation gate
using trilinear loop sampling. Keep that result. Test interpolation separately
without changing any evolved fields, source or coupling.

Compare trilinear interpolation with tensor-product four-point cubic Lagrange
interpolation on nodes(-1,0,1,2). These are established numerical methods.
For fractional coordinate t, weights are

    [-t(t-1)(t-2)/6, (t+1)(t-1)(t-2)/2,
     -(t+1)t(t-2)/2, (t+1)t(t-1)/6].

Negative cubic weights are acceptable here for measuring a saved vector field;
this does not replace positive interpolation in the matter/light Hamiltonian.

Controls use A(x)=.02 exp(-r^2/1.4^2) (axis cross x), L8, n24/32/48/64,
loop radii1,1.5,2 and axes e_z or the LR4 rotated e_z. Rotate the loop plane
together. Exact circulation is2*pi*.02*radius^2*exp(-radius^2/1.4^2).
Use256/512 loop samples. Require at n64 relative error<2% for linear and
<0.1% for cubic, at512 samples, for every radius/orientation. Report errors
at all grids and256/512 differences; do not infer spatial accuracy just from
angular quadrature. Test cubic weight sum and polynomial reproduction for
powers0..3 at33 fractions in[0,1], absolute error<1e-12.

Then measure the four saved LR4 final A fields with both methods, both loop
counts and all three radii. Recompute the n48/n64 spatial difference and the
unrotated/rotated n48 difference at radius1.5. Report the previously declared
5% spatial and2% rotation thresholds for both methods but do not relabel the
original campaign. This is a diagnostic follow-up, not an independent new
evolution or opportunity to select a preferred result after execution.

Archive source and input hashes. Even a successful higher-order measurement
does not certify field evolution, all radii, another observable or physical
rotational invariance. It can only identify how much the measured discrepancy
depends on this interpolation choice. A new independent refinement/rotation
campaign would still be needed for a replacement accuracy claim.
