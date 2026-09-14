"""Unified incident companion field for the RB-1 consistency revision.

c=1, product rest energy m=1, heavy receivers (threshold E'=1 in the receiver
frame). A field is a callable I(E, n): number intensity per unit bath-frame
energy and solid angle, averaged over the sphere to 1 for an isotropic unit
bath. Bound production, absorbed power, receiver force and receiver energy
change are all computed from the same I(E, n). See protocol-addendum.md.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss, legval
import kinematics as kin

SPECTRA = dict(threshold_cut=(1., 2.), extends_below=(.5, 2.))   # flat dN/dE, units of the rest threshold


def flat(spec):
    lo, hi = SPECTRA[spec]
    return lambda E: ((E >= lo) & (E <= hi)).astype(float)


def separable(spec, A):
    """I(E, n) = S(E) A(n_z), with the local outward radial direction along z."""
    S = flat(spec)
    return lambda E, n: S(E)*A(n[..., 2])


def boosted_isotropic(spec, w):
    """Isotropic flat bath in a frame moving with velocity w: I_N(E) = D^2 I'_N(E/D)."""
    S, w = flat(spec), np.asarray(w, float)
    g = 1/np.sqrt(1 - w@w)

    def field(E, n):
        D = 1/(g*(1 - n@w))
        return D*D*S(E/D)
    return field


def reference_attenuation(x, T, mu):
    """exp(-tau) of the reference model (paper equation 6) at x=r/a, propagation cosine mu."""
    x = np.asarray(x, float)[..., None]
    t = x*mu
    B2 = 1 + x*x*(1 - mu*mu)
    B = np.sqrt(B2)
    tau = T*(t/(2*B2*(B2 + t*t)) + (np.arctan(t/B) + np.pi/2)/(2*B**3))
    return np.exp(-np.maximum(tau, 0))


def legendre_field(coeffs):
    return lambda mu: legval(mu, coeffs)


def moments(A, n_mu=192):
    """Energy density, radial flux and pressure components of an axisymmetric A(mu), sphere-averaged."""
    mu, w = leggauss(n_mu)
    vals = A(mu)
    w = w/2
    u = vals@w
    Pp = (vals*(1 - mu*mu)/2)@w
    return dict(u=u, F_r=(vals*mu)@w, P_rr=(vals*mu*mu)@w, P_perp=Pp, p_perp=Pp/u, zeta=((vals*mu)@w)/u)


def spectrum_constants():
    """Reacting-support integrals, identical for both flat controls (sigma=0 below threshold)."""
    return dict(rate=kin.REFERENCE_RATE, power=kin.REFERENCE_POWER, chi=kin.REFERENCE_CHI)


def first_order_force(mom, beta_vec):
    """Force per unit (Sigma_P u): F_A/u - (1+chi) P.beta/u - beta, with z the radial direction."""
    chi = kin.REFERENCE_CHI
    b = np.asarray(beta_vec, float)
    Pb = np.array([mom['P_perp']*b[0], mom['P_perp']*b[1], mom['P_rr']*b[2]])
    return (np.array([0., 0., mom['F_r']]) - (1 + chi)*Pb)/mom['u'] - b


def directions(n_mu=48, n_phi=48):
    mu, w = leggauss(n_mu)
    phi = (np.arange(n_phi) + .5)*2*np.pi/n_phi
    M, PH = np.meshgrid(mu, phi, indexing='ij')
    s = np.sqrt(1 - M*M)
    n = np.stack([s*np.cos(PH), s*np.sin(PH), M], axis=-1).reshape(-1, 3)
    return n, np.repeat(w/2, n_phi)/n_phi


def bound_rate(field, beta_vec, v_esc, n_u=160, n_mu=64, n_phi=64):
    """Rate of galaxy-frame bound products per unit time, same normalization as absorbed power.

    Receiver-frame product speed u has E' = 1/sqrt(1-u^2) and sigma = u. Per unit
    E' the rate from direction n is I(E_b, n)/gamma with E_b = E'/(gamma(1-n.beta)).
    Bound means |beta + u n*| < v_esc for isotropic emission n* (nonrelativistic
    addition, adequate at the speeds used).
    """
    b = np.asarray(beta_vec, float)
    beta = np.linalg.norm(b)
    g = 1/np.sqrt(1 - beta*beta)
    n, wn = directions(n_mu, n_phi)
    edges = sorted({0., abs(v_esc - beta), v_esc + beta})
    total = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        x, w = leggauss(n_u)
        u = (hi - lo)*(x + 1)/2 + lo
        wu = w*(hi - lo)/2
        Ep = 1/np.sqrt(1 - u*u)
        dxdu = u/(1 - u*u)**1.5
        Eb = Ep[:, None]/(g*(1 - n@b))[None, :]
        D = (field(Eb, np.broadcast_to(n, Eb.shape + (3,)))@wn)/g
        if beta > 0:
            Pb = np.clip((v_esc**2 - (u - beta)**2)/(4*u*beta), 0, 1)
            Pb[u <= v_esc - beta] = 1.
        else:
            Pb = (u < v_esc).astype(float)
        total += np.sum(wu*u*dxdu*D*Pb)
    return total


def energy_nodes(lo, hi, beta, n):
    """Gauss nodes on a flat support, split across the Doppler window around threshold."""
    g = 1/np.sqrt(1 - beta*beta)
    cuts = sorted({lo, hi, *[c for c in (1/(g*(1 + beta)), 1/(g*(1 - beta))) if lo < c < hi]})
    E, W = [], []
    for a, b in zip(cuts[:-1], cuts[1:]):
        x, w = leggauss(n)
        E.append((b - a)*(x + 1)/2 + a)
        W.append(w*(b - a)/2)
    return np.concatenate(E), np.concatenate(W)


def absorbed_power(field, beta_vec, spec_range, n_E=400, n_mu=64, n_phi=64):
    """Bath-frame absorbed energy rate with exact receiver-frame energies (heavy receiver)."""
    b = np.asarray(beta_vec, float)
    beta = np.linalg.norm(b)
    g = 1/np.sqrt(1 - beta*beta)
    n, wn = directions(n_mu, n_phi)
    E, wE = energy_nodes(*spec_range, beta, n_E)
    flux = 1 - n@b
    total = 0.
    for Ei, wi in zip(E, wE):
        Ep = Ei*g*flux
        total += wi*Ei*np.sum(wn*flux*kin.sigma_swave(Ep - 1)*field(np.full(len(n), Ei), n))
    return total


def exact_channels_boosted(spec, w, beta_vec, M=1e7, n_E=300, n_mu=48, n_phi=48):
    """exact_channels for a flat bath isotropic in a frame moving with w.

    Integrates over bath-rest-frame energies E0, whose support edges are fixed:
    E = D E0 with D = 1/(gamma_w (1 - n.w)), and I_N dE = D^3 S0(E0) dE0.
    """
    b, w = np.asarray(beta_vec, float), np.asarray(w, float)
    n, wn = directions(n_mu, n_phi)
    D = 1/((1/np.sqrt(1 - w@w))*(1 - n@w))
    # Threshold E'=1 sits at E0 = 1/(D gamma (1-n.beta)); split nodes across its range.
    kink = 1/(D*(1/np.sqrt(1 - b@b))*(1 - n@b))
    lo, hi = SPECTRA[spec]
    cuts = sorted({lo, hi, *[c for c in (kink.min(), kink.max()) if lo < c < hi]})
    E0, W0 = [], []
    for a, c in zip(cuts[:-1], cuts[1:]):
        x, q = leggauss(n_E)
        E0.append((c - a)*(x + 1)/2 + a)
        W0.append(q*(c - a)/2)
    u = np.broadcast_to(b, (len(n), 3))
    K, X, dR, rate = np.zeros(4), np.zeros(4), np.zeros(4), 0.
    for Ei, wi in zip(np.concatenate(E0), np.concatenate(W0)):
        ev = kin.react(D*Ei, n, u, M)
        wt = ev['weight']*(1 - n@b)*D**3*wn*wi
        K += wt@ev['k']
        X += wt@ev['X']
        dR += wt@(ev['R'] - ev['P'])
        rate += wt.sum()
    return dict(rate=rate, K=K, X=X, dR=dR, P_abs=K[0],
                residual=float(np.max(np.abs(K - X - dR))/K[0]))


def exact_channels(field, beta_vec, spec_range, M=1e7, n_E=300, n_mu=48, n_phi=48):
    """Exact reaction-weighted sums: incident K, product X and receiver change dR four-momenta."""
    b = np.asarray(beta_vec, float)
    n, wn = directions(n_mu, n_phi)
    E, wE = energy_nodes(*spec_range, np.linalg.norm(b), n_E)
    u = np.broadcast_to(b, (len(n), 3))
    K, X, dR, rate = np.zeros(4), np.zeros(4), np.zeros(4), 0.
    for Ei, wi in zip(E, wE):
        ev = kin.react(np.full(len(n), Ei), n, u, M)
        w = ev['weight']*(1 - n@b)*field(np.full(len(n), Ei), n)*wn*wi
        K += w@ev['k']
        X += w@ev['X']
        dR += w@(ev['R'] - ev['P'])
        rate += w.sum()
    return dict(rate=rate, K=K, X=X, dR=dR, P_abs=K[0],
                residual=float(np.max(np.abs(K - X - dR))/K[0]))
