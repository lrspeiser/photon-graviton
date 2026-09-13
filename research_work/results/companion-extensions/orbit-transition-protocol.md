# Free orbital transition radius in the all-motion diagnostic

13 September 2026. Written before the fit.

Retain the all-motion experiment's six systems, two alternative population proxies, full covariance, fixed exact-third companion model, lens-calibrated stellar masses, gradient profile and bounds. Change only the orbit transition radius from Re to ra, fitted over 0.1<=ra/Re<=10. Optimize log(ra/Re). The stellar mass-to-light gradient still transitions at Re. No new capture parameter or lensing multiplier is introduced.

The known orbital family becomes beta(r)=beta0+(beta_infinity-beta0) r^2/(r^2+ra^2). Update its Jeans integrating factor and necessary gamma>=2 beta cap consistently with ra. Keep the central asymptotic limit beta0<=0.375; validate the final condition on the same doubled 8193-point grid. This is a necessary condition in the specified separable augmented-density class, not proof of positive distribution functions or stability.

Fit all 40 motion bins, including outer bins, with catalogue lens angles consumed as stellar-mass calibration. Thus every score is descriptive, not independent validation. Add one nuisance parameter per galaxy, giving four fitted stellar parameters rather than three. Use the prior 22 starting configurations at ra/Re=1 plus the prior all-motion optimum at ra/Re=0.1, 1/3, 1, 3 and 10. Require reproduction of the fixed-radius result at ra=Re and a final objective no worse than that feasible result. Report all systems, boundaries, inner/outer decomposition, transitions and failed optimizer starts; do not call multistart optimization proof of a global optimum.

A useful improvement must be assessed alongside remaining residuals and the additional flexibility. No branch is promoted to the companion reference on this exposed-data fit. If the radius reaches a limit, record that dependency rather than claim a measured scale.
