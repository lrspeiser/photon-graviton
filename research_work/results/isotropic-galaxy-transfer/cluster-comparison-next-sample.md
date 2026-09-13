# Beyond the six Coma bins

## What the present result resolves

The initial inner-three/outer-three split gave an unstable impression of NFW's performance. The extended diagnostic reproduces the published full-profile NFW score approximately and shows comparable leave-one-bin-out errors for NFW and transparent companions. Strong interception does not retain its original apparent advantage. A point mass remains only a deliberately inadequate control; it is not in the physical candidate shortlist.

The existing comparison is about projected shapes. It is not yet a calculation of cluster gravity from incoming radiation. For any constant target-level retention eta, a fitted normalization B makes B*eta indistinguishable from any other positive eta after refitting B. Therefore the one-third exponent cannot be established by this experiment.

## Concrete next comparisons

1. **An ensemble of resolved cluster shear profiles.** [LoCuSS, Okabe and Smith 2016](https://arxiv.org/abs/1507.04493) presents measurements of 50 clusters and reports good agreement of its stacked profile with NFW. This is an important countercheck to any suggestion that the Coma outer-only ranking generalizes. Use the individual angular reduced-shear observations, their source calibration and covariance; do not use NFW-derived masses as observations. The publication was inspected; a machine-readable profile package was not acquired during this turn.
2. **Lensing plus independent observables.** [CLASH, Umetsu et al. 2014](https://arxiv.org/abs/1404.1375) combines shear and magnification for 20 clusters. This offers more than one lensing observable, but its reconstructed mass profiles carry geometry and gravity assumptions. Refit observables with declared geometry before using them to test our universe. No CLASH numerical fit has been performed here.
3. **Gas and stars in the same targets.** For each target, include gas and stellar density profiles in all three models. This prevents the companions or fitted halo from silently replacing ordinary matter. MOND must use this extended ordinary-matter distribution, not the compact point-baryon proxy used in our exploratory Coma shape check.

## Equations to use, with status

Known spherical Newtonian mass bookkeeping:

    M_b(<r) = 4*pi*integral_0^r rho_b(s)*s^2 ds
    g_b(r) = G*M_b(<r)/r^2

Known NFW comparison, with cluster gas and stars kept separate:

    rho_h(r) = rho_s/[x*(1+x)^2], x=r/r_s
    g_total(r) = G*(M_b(<r)+M_h(<r))/r^2

Known simple MOND interpolation under a spherical, isolated approximation:

    g_total = g_b/2 + sqrt(g_b^2/4 + a0*g_b)

Its lensing completion must be specified; this acceleration equation alone is insufficient. Fix the previously fitted galaxy a0 for a transfer test. If the cluster requires additional matter, report that residual rather than silently adding it.

Our proposed companion interpretation, using standard transport and energy bookkeeping:

    d u_dep / dt = eta * P_absorbed_per_volume - P_released_per_volume
    rho_dep = u_dep/c^2                       [assumed ordinary gravitational coupling]
    g_total = G*(M_b(<r)+M_dep(<r))/r^2       [spherical, nonrelativistic limit]
    eta = X^(1/3)/(1+X^(1/3))                [known Hill function; proposed use]

Do not insert a newly fitted cluster X without declaring a new candidate: the galaxy input used disk luminosity and disk scale, which do not automatically define a cluster quantity. The deposited component also needs a stress/support prescription before the density can be assumed stationary. No universe age or size is fixed to manufacture accumulation.

Known weak-lensing notation:

    gamma_t = (mean_Sigma_inside_R - Sigma_at_R)/Sigma_crit
    kappa = Sigma/Sigma_crit
    reduced_shear = gamma_t/(1-kappa)

These equations assume the usual lensing response to mass; alternative field stress requires its own lensing response. Compare the measured quantity (often reduced shear) instead of treating all profiles as weak gamma. Adopt distances explicitly under the user's stipulated-distance convention, or compute them with the candidate optical law. Neither published halo masses nor cosmological critical-density normalizations are required as input facts.

## Freeze rules

Specify candidate source/capture laws and target split before new fitting. An initial cluster fit can calibrate shared quantities, followed by predictions on different clusters. Report separately a freely fitted halo per cluster, a shared halo mapping, fixed galaxy MOND, and companion transfer with fixed normalization. They use different amounts of target information. Do not rank them as equal-information predictions without stating those differences.

The current turn completed the Coma diagnostic and literature comparison, not these ensemble fits. All six goals remain open.
