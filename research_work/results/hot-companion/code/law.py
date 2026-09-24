"""Hot-companion gravity: the law and its ingredients.

Mechanism (proposed here; originality unverified -- see ../README.md for prior art):

 1. All ordinary matter feeds a companion field.
 2. Ordered (cold) matter feeds it COHERENTLY. Contributions from different places keep a
    fixed relative phase and add as vectors, so they cancel exactly the way Newton's pull
    does. The coherent intensity at a point equals |g_N| there (Gauss's law: established).
 3. Randomly moving (hot) matter feeds it INCOHERENTLY, in proportion to its random kinetic
    energy measured against the companion's own speed u: weight k = 3 sigma^2 / u^2.
    Incoherent contributions add as intensities -- no cancellation:
        S_hot(x) = G * integral k(x') rho(x') / |x - x'|^2 dV'
 4. The companion pulls with its AMPLITUDE, the square root of its total intensity times a
    constant a, along the net (coherent) pull.
 5. Strong fields hold the companion ATTACHED (inactive); it is released where the
    ordinary pull is weak, like a thermally activated escape. The released fraction is
        f = exp( -|g_N| / g_d ),     g_d = lambda * a
    This switch-off is required by the Solar System, where any surviving companion pull
    would be measured, and it also improves the galaxy fits.
 6. Total pull on matter and light (no slip assumed):
        g = g_N + f * sqrt( a * ( |g_N| + S_hot ) )

Three constants: a and lambda (fitted once, on SPARC rotation curves) and u (fitted once,
on X-COP cluster masses). Everything else is predicted.

Units: kpc, km/s, Msun; accelerations in (km/s)^2/kpc unless marked _SI.
"""
from __future__ import annotations
import numpy as np

G = 4.30091727003628e-6                     # kpc (km/s)^2 / Msun
KMS2_PER_KPC = 1e6 / 3.0856775814913673e19  # (km/s)^2/kpc -> m/s^2
C_KMS = 299792.458


def w_inside(x):
    """Shell-averaged 1/d^2 for a shell of radius s = x R inside the field point (x<1),
    in units of 1/R^2: artanh(x)/x. Established geometry."""
    x = np.clip(np.asarray(x, float), 1e-12, 1 - 1e-12)
    return np.arctanh(x) / x


def w_outside(x):
    """Same for a shell outside (x>1): arcoth(x)/x -> 1/x^2 far away."""
    x = np.maximum(np.asarray(x, float), 1 + 1e-12)
    return 0.5 * np.log((x + 1) / (x - 1)) / x


def shell_weights(R, s):
    """Matrix W[i,j]: 1/d^2 averaged over shell j (radius s_j) at radius R_i, times R_i^2."""
    R = np.atleast_1d(R)[:, None]; s = np.asarray(s)[None, :]
    return np.where(s < R, w_inside(s / R), w_outside(s / R))


def scalar_sum(R, s, dm, k):
    """S_hot at radii R for spherical shells of mass dm at radii s with heat weight k."""
    R = np.atleast_1d(R)
    return G * (shell_weights(R, s) @ (np.asarray(k) * np.asarray(dm))) / R ** 2


HEAT_P = 2.0   # round 14: the heat weight's exponent, k = 3 (sigma / u)^p; the law has 2. Set by the regression
               # suite from the law's 'heat_exponent' (regression/common.py) before any test runs.


def heat_weight(sigma_kms, u_kms):
    """k = 3 sigma^2 / u^2: random kinetic energy against the companion speed. With HEAT_P = p != 2,
    k = 3 (sigma / u)^p (round 14)."""
    if HEAT_P == 2.0:
        return 3.0 * np.asarray(sigma_kms) ** 2 / u_kms ** 2
    return 3.0 * (np.abs(np.asarray(sigma_kms)) / u_kms) ** HEAT_P


def k_from_sig2(sig2, u_kms, pref=3.0):
    """The heat weight from a mean-square speed: pref sig2 / u^2 (pref = 3 for a 1D dispersion sig2; 1 for a
    3D mean square; 3 - 2 beta for an anisotropic radial one). With HEAT_P = p != 2, the same effective
    1D dispersion sigma_eff^2 = pref sig2 / 3 enters as k = 3 (sigma_eff / u)^p (round 14). At p = 2 the
    arithmetic is the scripts' own, so results are unchanged to the last digit."""
    if HEAT_P == 2.0:
        return pref * sig2 / u_kms ** 2
    return 3.0 * (np.maximum(pref * np.asarray(sig2) / 3.0, 0.0) / u_kms ** 2) ** (HEAT_P / 2.0)


def released(gN, a, lam):
    """Fraction of the companion released from attachment: exp(-|g_N| / (lambda a))."""
    return np.exp(-np.abs(np.asarray(gN)) / (lam * a))


def total(gN, S_hot, a, lam=np.inf):
    """g = g_N + f sqrt(a (|g_N| + S_hot)); lam=inf switches the release factor off."""
    gN = np.asarray(gN)
    f = 1.0 if not np.isfinite(lam) else released(gN, a, lam)
    return gN + f * np.sqrt(a * (np.abs(gN) + S_hot))


def point_mass(a, u, lam=np.inf):
    """Isolated point mass with internal dispersion sigma: S_hot = k g_N exactly."""
    def f(r, M, sigma_kms):
        gN = G * M / np.asarray(r) ** 2
        return total(gN, heat_weight(sigma_kms, u) * gN, a, lam)
    return f


# ---- context laws (NOT ours; used only as labelled benchmarks) ----------------------
def mond_simple(gN, a):
    """Milgrom's MOND with the 'simple' interpolation function. Context only."""
    return gN / 2 + np.sqrt(gN ** 2 / 4 + gN * a)


def nfw_velocity2(r, v200, c, h=0.7):
    """NFW halo circular velocity^2 (Navarro, Frenk & White 1996). Dark-matter context."""
    r200 = v200 / (10 * h * 0.1)      # kpc for v200 in km/s with H0 = 100h km/s/Mpc
    x = np.asarray(r) / r200
    m = lambda y: np.log(1 + y) - y / (1 + y)
    return v200 ** 2 * m(c * x) / (x * m(c))
