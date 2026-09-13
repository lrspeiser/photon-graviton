# All-motion diagnostic with calibrated lensing and constrained orbits

13 September 2026. Recorded before execution.

Fit the existing three-parameter stellar gradient/radial-anisotropy family to every released stellar-motion bin, using the full covariance. Keep exact-third companions, galaxy-trained amplitude, optical geometry, PSF, light profiles, gradient and orbital bounds, central cusp limit and sampled gamma>=2 beta constraint unchanged. At each gradient the catalogue lens angle still calibrates stellar mass.

This consumes the formerly reserved outer bin. It is a representational diagnostic on exposed data, not a prediction or a validation score. Compare total motion chi-squared at the previous inner-only constrained optimum against the new all-bin optimum, and report the inner and covariance-conditioned outer contributions separately. The standard block-covariance identity requires total chi-squared=inner chi-squared+conditional outer residual squared. Verify it numerically.

Use the previous 21 starts plus the inner-only constrained optimum. Require recovery of the earlier nested-model inner scores, positive masses and moments, lens equation closure, and central/refined slope margins. The all-bin objective must not exceed its value at the inner-only optimum. Report every galaxy, parameter bounds, optimizer outcomes, and residuals. No statistical significance or absolute rejection is inferred from local optimization alone. Compare model limitations rather than searching for a passing data subset.
