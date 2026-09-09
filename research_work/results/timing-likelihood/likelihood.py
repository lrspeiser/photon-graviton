"""Normalized flux likelihood with uncertain positive amplitude and baseline.

This is a building block, not a calibrated duration or population estimator.
"""
import numpy as np
from scipy.special import log_ndtr


def log_flux_likelihood(flux, error, shape, amplitude_mean=1.,
                        amplitude_sigma=.5, baseline_sigma=.1):
    """Integrate A>0 and B for y=A*s+B+noise.

    A has a normal prior truncated at zero; B has a zero-mean normal
    prior. Errors are known independent Gaussian standard deviations.
    Priors are part of this optional statistical model, not measured facts.
    """
    y, err, s = (np.asarray(x, dtype=float) for x in (flux, error, shape))
    if y.ndim != 1 or y.shape != err.shape or y.shape != s.shape:
        raise ValueError('Matching one-dimensional arrays are required')
    if not all(np.all(np.isfinite(x)) for x in (y, err, s)):
        raise ValueError('Nonfinite inputs')
    if np.any(err <= 0) or not np.isfinite(amplitude_mean):
        raise ValueError('Positive errors and finite prior mean required')
    if not (np.isfinite(amplitude_sigma) and amplitude_sigma > 0
            and np.isfinite(baseline_sigma) and baseline_sigma > 0):
        raise ValueError('Positive finite prior standard deviations required')
    if not len(y):
        return 0.
    # Integrate the two Gaussian nuisance parameters using a 2x2 system.
    X = np.column_stack((s, np.ones_like(s)))
    mu = np.array([amplitude_mean, 0.])
    prior_variance = np.array([amplitude_sigma**2, baseline_sigma**2])
    residual = y-X@mu
    weights = 1/err**2
    precision = np.diag(1/prior_variance) + X.T@(weights[:, None]*X)
    chol = np.linalg.cholesky(precision)
    posterior_covariance = np.linalg.solve(chol.T, np.linalg.solve(chol, np.eye(2)))
    rhs = X.T@(weights*residual)
    posterior_mean = mu+posterior_covariance@rhs
    quadratic = residual@(weights*residual)-rhs@posterior_covariance@rhs
    logdet = (np.sum(np.log(err**2))+np.sum(np.log(prior_variance))
              +2*np.sum(np.log(np.diag(chol))))
    log_untruncated = -.5*(len(y)*np.log(2*np.pi)+logdet+quadratic)
    log_positive_posterior = log_ndtr(posterior_mean[0]/np.sqrt(posterior_covariance[0, 0]))
    log_positive_prior = log_ndtr(amplitude_mean/amplitude_sigma)
    return float(log_untruncated+log_positive_posterior-log_positive_prior)


def width_shape(time, width, peak, rise_fraction, rise_power, fall_power):
    """Positive unimodal pulse with exact full width at half maximum W.

    Neither redshift nor a predetermined time-stretch exponent is an input.
    """
    if not (width > 0 and 0 < rise_fraction < 1 and rise_power > 0 and fall_power > 0):
        raise ValueError('Invalid pulse parameters')
    t = np.asarray(time, dtype=float)-peak
    halfwidth = np.where(t < 0, rise_fraction*width, (1-rise_fraction)*width)
    power = np.where(t < 0, rise_power, fall_power)
    return np.exp(-np.log(2.)*(np.abs(t)/halfwidth)**power)
