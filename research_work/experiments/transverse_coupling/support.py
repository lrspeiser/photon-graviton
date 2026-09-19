"""The steering law's mean support factor f(beta) = 1 - <v_r^2/v^2> for an anisotropic Gaussian velocity
distribution with sigma_t^2/sigma_r^2 = 1 - beta, by the Laplace identity 1/S = int exp(-sS) ds
(established mathematics); kept free of other imports so the lens stage can use it."""
import numpy as np
from scipy.integrate import quad


def support_factor(beta):
    q = 1 - beta
    val = quad(lambda s: (1 + 2*s)**-1.5/(1 + 2*s*q), 0, np.inf, epsabs=0, epsrel=1e-12, limit=400)[0]
    return 1 - val
