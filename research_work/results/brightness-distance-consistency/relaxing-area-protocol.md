# Relaxing optical response: calibration and transfer

This exploratory candidate is declared before its fit. Previous comparisons have exposed both the supernova partitions and the six lens systems; this is not a blind discovery test.

The retained rational response B=1+f/(1+q*f) improves supernova brightness but leaves a high-redshift residual and worsens the six-system lens aggregate. Test a single alternative shape with the same number of parameters:

    f = z/(1+z)
    B(f) = 1 + f*exp(-q*f), q >= 0
    D_G = D*sqrt(B)
    D_A = D_G/(1+z)

B is a proposed beam-area response, not a derived time or companion interaction. Exponential relaxation is known mathematics; no novelty of the functional form is claimed. Its rationale is a finite optical response which can turn over rather than remain monotonic. The rate q is empirical. B(0)=1 and B'(0)=1 preserve the necessary regular observer condition of the retained metric-geodesic branch. Its observer focusing is alpha^2*(3*q+13/4), as for the rational branch at the same q. This is not a global metric/source solution.

Keep alpha and event-stretch exponent b=1 fixed. Fit q in [0,100] using the same 466 rows at 0.1<=zHD<0.3 and the same 77 calibrator rows and covariance. Freeze q before evaluating the 494 farther rows. Use the same scoring and calibration propagation as the rational candidate. No additional parameter or high-redshift refitting.

Propagate this frozen H(f)=-(1-f)*ln(1-f)*sqrt(B) through the existing Jacobi lens-distance calculation. Keep SPARC capture parameters fixed; refit only the previously declared ordinary stellar mass and orbital nuisance parameters to inner stellar measurements. Outer stellar bins and lens angles remain predictions. Preserve earlier outputs. Compare both farther brightness score and lens fractional RMS with the rational response; an improvement in only one does not establish a joint improvement. The earlier morphology, source geometry, spherical gravity, energy and deposit-support limitations still apply.
