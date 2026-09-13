# External collection and optional depth-based capture

Status: hypothetical transport calculation, not an observational fit. No universe age or size is fixed. All six research objectives remain open.

## Collecting area

For an isotropic incoming companion bath with energy density u_c, moving at c,

    deposited power = c u_c sigma_eff
                    = pi R^2 c u_c f_capture.

This uses known isotropic transport geometry: the inward flux per unit surface is c u_c/4, and the sphere's full surface is 4 pi R^2. Do not multiply the full surface by c u_c. For fixed incoming density and capture fraction, doubling R quadruples power. An effective cross section already includes capture efficiency:

    sigma_eff = 2 pi integral_0^infinity b [1-exp(-tau(b))] db,
    tau(b) = integral_along_ray kappa ds.

For uniform external stellar luminosity density j_star, straight propagation, isotropic emitters, and source distance D much larger than the receiving well, a shell contributes

    dP = sigma_eff j_star C(D) dD.

The shell's 4 pi D^2 increase in sources cancels inverse-square dilution. Under the proposed lossless conversion rule with constant coefficient alpha and no intervening capture, C(D)=1-exp(-alpha D). This fraction approaches one; it cannot exceed the emitted energy. With intervening capture, use the surviving companion fraction instead. The finite-history and competing-capture reports implement that correction. An arbitrarily extended, eternally luminous source distribution without competing absorption is not a finite supply prediction.

## Optional depth law

To test the user's emphasis on well depth, compare kappa=chi*(-Phi)^4, with Phi(infinity)=0 in an isolated well. Here -Phi=v_escape^2/2 is specific escape energy. This is an optional phenomenological postulate, not a first-principles result or a claim of unique mathematics. The fourth power supplies a convergent outer collecting area in this example; observations have not selected it. A physical potential reference in a general many-well universe remains necessary.

For a fixed Plummer potential Phi=-GM/sqrt(r^2+a^2), define x=r/a and A=chi*(GM)^4/a^3. Then a*kappa=A/(1+x^2)^2. The isotropic attenuated deposition profile obeys

    q(0)/(c u_infinity/a) = A exp(-pi A/4),
    q(x)/q(0) = 1 + [-2+pi A/4+A^2/6] x^2 + O(x^4).

The curvature changes sign at A=1.8332747453568032. Weak capture permits a central local maximum. Strong capture intercepts incoming energy before it reaches the center, yielding a central local minimum even though opacity is highest there. At the threshold the quadratic term alone does not classify the center. Three-dimensional curvature does not by itself establish the sign of projected lensing shear, or a viable equilibrium distribution of stored material.

| A | Effective area / a^2 | Half-deposition radius / a |
|---:|---:|---:|
| 0.1 | 0.968146 | 2.30862 |
| 1 | 8.383788 | 2.67757 |
| 1.833275 | 13.921098 | 2.97147 |
| 3 | 20.517390 | 3.31865 |
| 10 | 49.645273 | 4.63936 |

Between two identical wells, accelerations cancel at the midpoint but depths add. This candidate has nonzero midpoint opacity, eight times the sum of the separate fourth-power opacities there. That is a local opacity ratio, not an eightfold increase in total captured power; attenuation for the pair has not been solved.

## Verification and limits

run.py independently compares the ray cross section with the volume integral of absorbed power (relative discrepancy below 3.2e-8) and checks central curvature by finite differences (absolute coefficient discrepancy below 1.1e-6). results.json preserves all cases. These are mathematical consistency tests using synthetic wells. No actual cluster lensing or dynamics is fitted, chi is uncalibrated, and storage support, capture microphysics, backreaction, and a complete gravitational energy ledger remain unresolved. This candidate does not replace the acceleration-squared candidate as an established choice.
