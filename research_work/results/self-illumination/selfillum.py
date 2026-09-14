"""Self-illumination pilot: companion fields derived from declared rotating emitters.

Emitters emit isotropically in their own frames with a flat rest-frame number
spectrum. Steady, optically thin transport gives each emitter element a lab
beam with number flux per unit lab energy  w D s(E/D)/(4 pi d^2),  where
D = 1/(gamma_e (1 - n.beta_e)) and n points from emitter to receiver. The
RB-1 reaction is applied beam by beam with exact kinematics. c=1, m=1.
"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

CAPTURE = Path(__file__).resolve().parents[1]/'capture-to-orbit'
sys.path.insert(0, str(CAPTURE))
import incident as inc  # noqa: E402
import kinematics as kin  # noqa: E402


def ring(R, beta, n=2048):
    th = (np.arange(n) + .5)*2*np.pi/n
    pos = np.stack([R*np.cos(th), R*np.sin(th), 0*th], 1)
    vel = beta*np.stack([-np.sin(th), np.cos(th), 0*th], 1)
    return pos, vel, np.full(n, 1/n)


def disk(Rd, beta, Rmax, h, nR=96, nth=256, nz=8):
    """Exponential disk, flat rotation curve, uniform slab of half-thickness h (softens the near field)."""
    x, wx = leggauss(nR)
    R, wR = Rmax*(x + 1)/2, wx*Rmax/2
    th = (np.arange(nth) + .5)*2*np.pi/nth
    z, wz = leggauss(nz)
    RR, TT, ZZ = np.meshgrid(R, th, h*z, indexing='ij')
    wt = (np.exp(-R/Rd)*R*wR)[:, None, None]*np.full(nth, 2*np.pi/nth)[None, :, None]*(wz/2)[None, None, :]
    pos = np.stack([RR*np.cos(TT), RR*np.sin(TT), ZZ], -1).reshape(-1, 3)
    vel = beta*np.stack([-np.sin(TT), np.cos(TT), 0*TT], -1).reshape(-1, 3)
    return pos, vel, (wt/wt.sum()).ravel()


def steady_flow_kappa(w_over_beta):
    """First-order drag at the centre of a fixed shell through which emitters flow with velocity w.

    kappa = (4/3 + chi/3) - (w/beta)(1 + chi/3): the isotropic value at w=0 and exactly 1/3 when the
    emitters co-move with the receiver. A steady emitting structure is not a co-moving bath; a bath
    isotropic in the receiver frame would give zero.
    """
    chi = kin.REFERENCE_CHI
    return 4/3 + chi/3 - w_over_beta*(1 + chi/3)


def shell(n_mu=32, n_phi=64, velocity=(0., 0., 0.), radius=1.):
    """Fixed spherical shell of emitting positions around the origin; emitters flow through it with one
    velocity (steady flow). Used for the isotropic and steady-co-moving controls."""
    n, w = inc.directions(n_mu, n_phi)
    return radius*n, np.broadcast_to(np.asarray(velocity, float), n.shape).copy(), w


def beams(pos, vel, w, x_r):
    dv = np.asarray(x_r, float) - pos
    d = np.linalg.norm(dv, axis=1)
    n = dv/d[:, None]
    D = 1/((1/np.sqrt(1 - np.sum(vel*vel, 1)))*(1 - np.sum(n*vel, 1)))
    return n, D, w/(4*np.pi*d*d)


def channels(n, D, fw, beta_r, spec, M=1e7, n_E=160):
    """Exact reaction-weighted sums over beams: incident K, products X, receiver change dR."""
    b = np.asarray(beta_r, float)
    gr = 1/np.sqrt(1 - b@b)
    Q = D*gr*(1 - n@b)                     # emitter-frame to receiver-frame energy factor
    lo, hi = inc.SPECTRA[spec]
    kink = 1/Q
    cuts = sorted({lo, hi, *[c for c in (kink.min(), kink.max()) if lo < c < hi]})
    K, X, dR = np.zeros(4), np.zeros(4), np.zeros(4)
    u = np.broadcast_to(b, n.shape)
    for a, c in zip(cuts[:-1], cuts[1:]):
        x, q = leggauss(n_E)
        for Ee, wq in zip((c - a)*(x + 1)/2 + a, q*(c - a)/2):
            ev = kin.react(D*Ee, n, u, M)
            wt = ev['weight']*(1 - n@b)*fw*D*D*wq
            K += wt@ev['k']
            X += wt@ev['X']
            dR += wt@(ev['R'] - ev['P'])
    return dict(K=K, X=X, dR=dR, P_abs=K[0], residual=float(np.max(np.abs(K - X - dR))/K[0]))


def bound_rate(n, D, fw, beta_r, v_esc, spec, n_u=96):
    """Galaxy-frame bound products per unit time, same normalization as channels()."""
    b = np.asarray(beta_r, float)
    beta = np.linalg.norm(b)
    gr = 1/np.sqrt(1 - beta*beta)
    Q = D*gr*(1 - n@b)
    lo, hi = inc.SPECTRA[spec]
    total = 0.
    for a, c in zip(*[sorted({0., abs(v_esc - beta), v_esc + beta})[i:] for i in (0, 1)]):
        x, q = leggauss(n_u)
        uu, wu = (c - a)*(x + 1)/2 + a, q*(c - a)/2
        Ep = 1/np.sqrt(1 - uu*uu)
        Ee = Ep[:, None]/Q[None, :]
        rate = ((((Ee >= lo) & (Ee <= hi))*(fw*D)[None, :]).sum(1))/gr
        if beta > 0:
            Pb = np.clip((v_esc**2 - (uu - beta)**2)/(4*uu*beta), 0, 1)
            Pb[uu <= v_esc - beta] = 1.
        else:
            Pb = (uu < v_esc).astype(float)
        total += np.sum(wu*uu*uu/(1 - uu*uu)**1.5*rate*Pb)
    return total


def emitter_loss(vel, w, spec):
    """Lab power and momentum-loss rate of isotropically emitting elements: P = gamma w <E_e>, dp = beta P."""
    lo, hi = inc.SPECTRA[spec]
    g = 1/np.sqrt(1 - np.sum(vel*vel, 1))
    P = g*w*(lo + hi)/2
    return P, vel*P[:, None]


def emission_check(beta=.01, spec='extends_below', n_mu=64, n_phi=64):
    """Integrate one element's beams over a sphere: number w, power gamma w <E>, momentum beta gamma w <E>."""
    n, wn = inc.directions(n_mu, n_phi)
    vel = np.array([beta, 0., 0.])
    D = 1/((1/np.sqrt(1 - beta*beta))*(1 - n@vel))
    lo, hi = inc.SPECTRA[spec]
    mean = (lo + hi)/2
    number = np.sum(wn*D*D)
    power = np.sum(wn*D**3)*mean
    momentum = np.sum(wn[:, None]*(D**3)[:, None]*n, 0)*mean
    g = 1/np.sqrt(1 - beta*beta)
    return dict(number=number, power_over_expected=power/(g*mean), momentum_over_expected=momentum[0]/(g*beta*mean))
