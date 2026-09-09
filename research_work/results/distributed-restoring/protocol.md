# Restoring-field robustness test

Freeze before calculation, 9 September 2026. Extend the exact homogeneous periodic finite-radiation sector with a positive quadratic scalar potential. This is an established potential form used as an optional candidate, not a new force claimed as unique. No fit to observed redshifts.

Use c0=K=g=volume=1, U0=0.003, n(0)=1, n_dot(0)=0. Compare m=0,0.01,0.03,0.1. For positive m integrate to 8 pi/m; the m=0 comparison runs to the same maximum time as m=0.01. No continuing radiation injection or external energy supply. The original radiation reservoir exchanges reversibly with the field.

Derived here within the proposed model; originality unverified: n_ddot=U0/n^2-m^2(n-1). Established energy accounting applied to the model: H=n_dot^2/2+m^2(n-1)^2/2+U0/n. Check total-energy relative error below 1e-7, convergence below 1e-6, and positive n. DOP853 tolerances 1e-9 and 1e-11 with absolute tolerance 1e-13 and 1e-15; limit steps to at most 0.5 to resolve turning events.

For m>0 derive the upper turning point and compare to the numerical first maximum. Detect the first return to n=1 and compare it to twice the first turning time. Verify the rate changes sign. Report late constant-rate failure if present; do not treat a temporary positive-rate interval as a permanent result. Distinguish homogeneous potential from environmental screening confined to dense regions. A failure here does not exclude zero restoring force in voids or a different driven mechanism.
