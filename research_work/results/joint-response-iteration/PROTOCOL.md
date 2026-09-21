# JR-1: iterative real-data companion-response construction

21 September 2026. User authorizes changes to equations and parameters, including exchange-time assumptions. This is exploratory model development, not a preregistered discovery or a completed microscopic derivation. Baseline c22d188949874ea0ef05b7e782d341973d421385. Input export 5cc2205ce71f68c694a4702c4a42a67bf6e02f4d is reproducibility infrastructure only.

## Theory first

Matter/radiation may build an extended companion state whose same gravitational response steers stars and bends light. This experiment estimates a spatial effective closure needed by that state. It does not equate a chosen halo-shaped ansatz with a microscopically produced deposit. No object receives a freely fitted companion amplitude, and no independent photon/lensing gain is used.

## Observations and roles

Use the archived 149 SPARC galaxies and all positive-radius raw rows, with the preserved 89 train / 29 validation / 31 comparison split. The raw baryonic velocity components give the disk's ordinary radial field; source normalization uses photometric stellar integration and catalog gas mass, not R*v_observed^2/G or force-equivalent enclosed mass. SPARC is principally gas-traced rotation, not 3152 independent stellar measurements.

The more stringent joint sample is the six training-role SLACS systems with released usable radial kinematics, covariance, detailed light profiles and conditional population masses. J0330-0020 is excluded by the existing release flag; J1538+5817 lacks detailed components. Fit global laws using J0037-0942, J1112+0826, J1204+0358 and J1402+6321. Reserve J1621+3931 and J1630+4520 from fitting in THIS experiment. All six were examined historically, so these two are out-of-fit transfer checks, not genuinely new blind data. Do not use their scores to select revisions.

Use the archived static-Euclidean conditional lens distances and energy_loss_and_event_stretch stellar normalization, Chabrier baseline. These are explicitly model-dependent conversions, not direct geometric distances or a freshly fitted stellar-population posterior. The released SIE Einstein angles are image-model summaries, not raw lens images. Their full uncertainties/covariance are unavailable here. Fit resolved V_rms with its actual covariance and the fixed published light profiles and PSF.

## Initial spatial closure

For source mass S, photometric effective radius Re and C=[G*S/Re^2]/[G*1e10 Msun/(1 kpc)^2], define

    A = Aref*(S/1e10 Msun)^p*C^dA
    rc = c*Re*C^dc
    rt = t*Re
    gchi(r) = A/r * (r/rc)^q/[1+(r/rc)^q] / sqrt[1+(r/rt)^2].

Here A has units (km/s)^2. The source mass is ordinary stellar plus gas mass; it is a production proxy, not a measured photon-energy conversion rate. The initial root uses p=.5, c=1, q=2, t=20, dA=dc=0 and the archived rounded acceleration scale for Aref. The finite companion mass implied by this force is A*rt/G. Its enclosed mass increases monotonically for q>0. A spherical companion envelope plus the catalog disk baryonic field is an explicit geometry choice, not an AQUAL solve or a microscopic result.

Temporal and spatial potentials are both sourced by the same total force (equal-potential hypothesis). Calculate annular Jeans second moments and lens deflection using that same force. No lensing multiplier or separately fitted halo per galaxy.

## Revisions

Compare the unmodified starting guess with successively fitted universal amplitude/core, source-mass/compactness response, core compactness and shape/outer scale, and shared stellar-normalization/orbital nuisance parameters. A further variant uses the steady production/single-release/pair-release occupation F(s) proportional to sqrt(1+4*s/s0)-1, normalized at s=1, with s=(S/Sref)^(2p). This connects a rate-balance shape but does not derive its physical rates. A radial-anisotropy extension may replace constant beta with beta0+(betainf-beta0)*r^2/(r^2+ra^2), universal ra/Re.

Every added parameter and bound is recorded by the executable before optimization; retain all candidates and optimizer starts. Candidate selection uses the training objective and SPARC validation only. Release SPARC comparison and the two out-of-fit lenses after selection; do not revise based on those scores in this experiment.

## Objective and reporting

Give SPARC rotation, resolved lens stellar kinematics, and lensing separate equal-weight mean standardized-square blocks so thousands of rotation rows cannot silently overwhelm six lens systems. Within each block use equal-object weighting. SPARC uses quoted velocity errors with no fitted error floor; resolved V_rms uses the supplied full covariance. Lensing uses a declared 5-percent working fractional deflection scale solely as an optimization tradeoff, NOT as a measured uncertainty or proof of statistical acceptance. Report physical-unit and fractional errors, raw kinematic chi-square, per-object residuals and exact predicted Einstein angles separately. This balanced development objective is not a normalized joint likelihood or a p-value.

## Numerical and physical limits

Use a single coupled mass-to-light scale only in explicitly labeled nuisance revisions. Do not rescale observed velocities or alter their errors. Double radial/angular/deprojection quadrature at the selected parameters, and compare lens deflections with independent adaptive integration. Report failed numerical comparisons. Source-mass and coefficient positivity are numerical/conditional consistency, not full dynamical stability.

Static observables do not identify a clock rate for the exchange. A common scaling of production and release rates changes time-dependent relaxation but not their stationary balance; include and report this exact degeneracy rather than claiming that a timeless exchange was measured. An instantaneous-equilibrium effective limit is allowed; nonlocal propagation, observable delays, causal ordering and energy flux remain separate open physics.

Even a good fit does not close microscopic production, energy supply, source selection/formation, lens-image geometry, Solar-system internal response, merger dynamics or a common gravitational-wave theory. Preserve those distinctions. The supplied portfolio permits multiple sectors and nonlocal/memory laws; failure of this spatial closure does not reject that portfolio.
