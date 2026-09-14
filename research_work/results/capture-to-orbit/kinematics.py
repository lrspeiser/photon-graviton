"""Exact two-body threshold kinematics for receiver-assisted production (RB-1).

c=1 and the product rest energy m=1 throughout. Receiver mass M is in units of
m and the incident null companion energy E is measured in the bath frame. The
kinematics and the Wigner s-wave threshold law are established physics; using
them for companion capture is the project hypothesis.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

C_KMS = 299792.458


def boost(p, beta):
    """Four-momenta p[..., 4]=(E, px, py, pz) seen from a frame moving with 3-velocity beta."""
    p, beta = np.asarray(p, float), np.asarray(beta, float)
    b2 = np.sum(beta*beta, axis=-1, keepdims=True)
    gam = 1/np.sqrt(1 - b2)
    bp = np.sum(beta*p[..., 1:], axis=-1, keepdims=True)
    coef = np.where(b2 > 0, (gam - 1)*bp/np.where(b2 > 0, b2, 1.), 0.) - gam*p[..., :1]
    return np.concatenate([gam*(p[..., :1] - bp), p[..., 1:] + coef*beta], axis=-1)


def isotropic(n, rng):
    mu = rng.uniform(-1, 1, n)
    ph = rng.uniform(0, 2*np.pi, n)
    s = np.sqrt(1 - mu*mu)
    return np.stack([s*np.cos(ph), s*np.sin(ph), mu], axis=-1)


def threshold(M):
    """Companion energy that reaches threshold on a receiver at rest: m + m^2/(2M)."""
    return 1 + 1/(2*M)


def react(E, n_c, u_rec, M, n_star=None):
    """Exact c + R -> R + X with isotropic s-wave emission in the CM frame.

    Returns bath-frame four-momenta of the incident companion k, the initial
    receiver P, the product X and the final receiver R, a kinematic mask and
    the s-wave weight sigma ~ p_f/p_i (constant |M|^2). The caller applies the
    flux factor (1 - n_c.u_rec). With n_star=None, X and R are averages over the
    isotropic emission direction; boosts are linear, so this is exact.
    """
    E = np.asarray(E, float)
    g = 1/np.sqrt(1 - np.sum(u_rec*u_rec, axis=1))
    k = np.concatenate([E[:, None], E[:, None]*n_c], axis=1)
    P = M*np.concatenate([g[:, None], g[:, None]*u_rec], axis=1)
    Ep = E*g*(1 - np.sum(n_c*u_rec, axis=1))          # receiver-frame companion energy
    W2 = M*M + 2*M*Ep
    W = np.sqrt(W2)
    above = 2*M*(Ep - 1) - 1                           # W^2 - (M+1)^2 without cancellation
    ok = above > 0
    pf = np.sqrt(np.maximum(above, 0)*(2*M*(Ep + 1) - 1))/(2*W)
    pi = M*Ep/W
    b = (k + P)[:, 1:]/(k + P)[:, :1]                  # CM velocity
    dirs = np.zeros((len(E), 3)) if n_star is None else n_star*pf[:, None]
    X = boost(np.concatenate([np.sqrt(1 + pf*pf)[:, None], dirs], axis=1), -b)
    R = boost(np.concatenate([np.sqrt(M*M + pf*pf)[:, None], -dirs], axis=1), -b)
    return dict(k=k, P=P, X=X, R=R, ok=ok, weight=np.where(ok, pf/pi, 0.), E_receiver=Ep)


def sigma_swave(x):
    """s-wave p_f/p_i for a heavy receiver at receiver-frame excess x=E'/m-1."""
    x = np.asarray(x, float)
    return np.sqrt(np.maximum(x*(2 + x), 0))/(1 + x)


# Declared reference spectrum: flat in bath-frame energy from threshold to twice threshold.
REFERENCE_RATE = np.sqrt(3) - np.pi/3                              # int_0^1 sigma dx
REFERENCE_POWER = (2*np.sqrt(3) - np.log(2 + np.sqrt(3)))/2         # int_0^1 sigma (1+x) dx
REFERENCE_CHI = np.arccosh(2.)/REFERENCE_POWER                       # <dln sigma/dln E>_P


def drag_coefficient(chi):
    """First-order receiver momentum loss per absorbed energy, in units of v/c^2.

    Absorption of an isotropic bath with sigma(E') gives -(beta/3)(1+chi) per
    absorbed energy, where chi is the power-weighted <dln sigma/dln E>; a product
    emitted isotropically in the receiver frame carries a further beta. chi=0 is
    the familiar gray (Compton-drag) value 4/3.
    """
    return 4/3 + chi/3


def bound_energy_fraction(v_esc_kms):
    """Retained rest energy per absorbed energy for the declared reference spectrum.

    Near threshold the s-wave rate per unit velocity-space volume is 1/(4 pi) in
    c=1 units, so the bound rate is v_esc^3/3; absorbed power per unit spectral
    density is REFERENCE_POWER. Nonrelativistic in v_esc.
    """
    v = np.asarray(v_esc_kms, float)/C_KMS
    return v**3/3/REFERENCE_POWER


def line_bound_fraction(delta, beta, v_esc):
    """Nonrelativistic bound fraction of reacting events for a bath-frame line.

    The receiver-frame excess is uniform on [delta-beta, delta+beta]; with the
    s-wave rate the products fill the velocity shell sqrt(2 max(0,delta-beta)) <
    u < sqrt(2(delta+beta)) uniformly. Valid when the shell contains the whole
    escape ball shifted by the receiver speed (delta<=beta and
    2(delta+beta) >= (v_esc+beta)^2); returns nan otherwise.
    """
    u_hi = np.sqrt(2*np.maximum(delta + beta, 0))
    ok = (delta <= beta) & (u_hi >= v_esc + beta)
    return np.where(ok, v_esc**3/np.where(u_hi > 0, u_hi, 1.)**3, np.nan)


def drag_quadrature(beta, E_range, M=1e4, n_energy=4000, n_mu=400):
    """Exact-kinematics receiver momentum loss per absorbed bath-frame energy / beta.

    Flat bath spectrum on E_range, isotropic directions, flux factor and s-wave
    weight included; emission averaged analytically. Returns the coefficient and
    the maximum four-momentum imbalance relative to M.
    """
    x, xw = leggauss(n_energy)
    E = (E_range[1] - E_range[0])*(x + 1)/2 + E_range[0]
    ew = xw*(E_range[1] - E_range[0])/2
    mu, mw = leggauss(n_mu)
    EE, MU = np.meshgrid(E, mu, indexing='ij')
    W = np.outer(ew, mw/2).ravel()
    n_c = np.stack([np.sqrt(1 - MU.ravel()**2), np.zeros(MU.size), MU.ravel()], axis=1)
    u = np.zeros((MU.size, 3))
    u[:, 2] = beta
    ev = react(EE.ravel(), n_c, u, M)
    w = W*ev['weight']*(1 - beta*MU.ravel())
    dpz = ev['R'][:, 3] - ev['P'][:, 3]
    kappa = -np.sum(w*dpz)/(np.sum(w*EE.ravel())*beta)
    imbalance = np.max(np.abs((ev['k'] + ev['P'] - ev['X'] - ev['R'])[ev['ok']]))/M
    return float(kappa), float(imbalance)
