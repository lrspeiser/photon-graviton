# A likelihood component for uncertain explosion durations

## Why this calculation is needed

The original timing estimator forced each event into one shape and returned one duration. Its injected-event checks failed, and the subsequent diagnosis found that fitted-duration errors could dominate acceptance selection. We now begin a different estimator that retains uncertainty instead of treating every fitted duration as exact.

For example, if observations end before an explosion has faded, the data may permit several durations. The revised method should carry that ambiguity into the population result. Rejecting the event because one fitted crossing falls outside the observations can select events according to duration and redshift; blindly extrapolating the crossing can produce a precise-looking wrong answer. Neither problem is fixed just by accepting more events.

This stage derives and independently checks a normalized probability calculation for uncertain brightness and background. It does **not** yet implement or validate the complete duration estimator. The original failed gate is unchanged.

## Shape and duration

Let W be the full width at half maximum, tau the observer-time peak, and r the fraction of that width before the peak. Define u=t-tau and

    s(t) = exp[-ln(2) (|u|/(r W))^p_rise]          for u<0,
           exp[-ln(2) (|u|/((1-r) W))^p_fall]     for u>=0.

The peak is one; the half-maximum points are exactly tau-rW and tau+(1-r)W, whose separation is W. Independent powers allow different rise and decline shapes. This empirical family is an optional statistical approximation, not a physical explosion law. Some shoulder structures and multiple peaks remain outside it and must be included in calibration trials. Flexibility alone is not proof of accurate recovery.

All times and W are in observer days. The individual-event function takes no redshift and imposes no division by (1+z).

## Integrating uncertain brightness and background

For a fixed proposed shape, measured fluxes obey the statistical model

    y = X theta + noise,  X=[s,1], theta=[A,B].

The diagonal noise covariance is D=diag(error_i^2). Before restricting amplitude to be positive, take theta~Normal(mu,P), with mu=(1,0) and P=diag(0.5^2,0.1^2). These are declared artificial-flux priors, not facts about actual supernovae. Marginalizing these two nuisance parameters gives

    C = D + X P X^T,
    L_unrestricted = Normal(y; X mu, C).

Positive brightness is imposed by conditioning the amplitude prior on A>0. By Bayes' rule the exact normalized likelihood is

    L_positive = L_unrestricted * Pr(A>0 | y) / Pr(A>0).

The unrestricted Gaussian posterior has

    V = (P^-1 + X^T D^-1 X)^-1,
    m = mu + V X^T D^-1 (y-X mu).

Thus the probability ratio is Phi(m_A/sqrt(V_AA))/Phi(mu_A/sqrt(P_AA)), where Phi is the standard normal cumulative probability. The implementation uses log probabilities, a 2-by-2 Cholesky solve and the matrix determinant identity. The determinant and prior normalization cannot be dropped: they are part of the probability, particularly when different trial durations or shapes accommodate different ranges of brightness and baseline.

The independent check forms the full observation covariance C instead, and separately integrates over positive A numerically after integrating B. It also checks that the likelihood integrates to one over a single measured flux. These checks verify the formula and implementation under the stated Gaussian model. They do not verify that the noise, prior or shape model describes real explosions.

## What still has to be built and tested

The population model must integrate each event likelihood over its uncertain W, peak and shape. A candidate width distribution is

    log W ~ Normal(a + b log(1+z), sigma_intrinsic^2).

Here b is inferred rather than fixed. A broad event likelihood can then provide weak information, without inventing an exact duration. This treatment is still sensitive to shape and population priors; an unobserved tail cannot be reconstructed from nothing. Numerical integration limits, prior sensitivity, independent-shape injections, uncertainty coverage and population-bias criteria must be recorded before evaluating the revised estimator.

The component protocol specifies what is fixed now and which estimator settings are still unspecified. It must not be presented as a complete preregistration of a population inference. No real flux is read. Real filter throughput, source evolution, survey selection and calibration remain separate requirements. In particular, intrinsic duration evolution and propagation may enter through the same sum in b; successful measurement alone will not identify the physical cause.

## Reproduction

Requires Python, NumPy and SciPy (no network or DES cache for this component):

```sh
python research_work/results/timing-likelihood/check.py
```

Generated output goes to `research_work/generated/timing-likelihood`. This optional component check does not add a default physics-suite job or establish a passed timing feasibility gate.
