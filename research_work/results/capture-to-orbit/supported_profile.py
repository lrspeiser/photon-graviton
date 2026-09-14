"""Supported-profile interface for the capture-to-orbit test.

Spherical monopole potentials, orbit apsides and time fractions, steady
phase-mixed population profiles and explicitly labeled support closures.
No fitting and no file I/O. Units: kpc, km/s, Msun; time unit kpc/(km/s).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

G = 4.30091727003628e-6            # kpc (km/s)^2 / Msun
GYR = 0.9777922216807891           # Gyr per kpc/(km/s)


class SphericalPotential:
    """Newtonian monopole of a nondecreasing enclosed mass M(<r) on a log grid.

    Inside the first grid radius the density is uniform; beyond the last radius
    all mass is enclosed, so the potential is exactly Keplerian there. Use a
    fine grid: M is interpolated linearly in ln r.
    """

    def __init__(self, r, mass):
        r = np.asarray(r, float)
        mass = np.maximum.accumulate(np.clip(np.asarray(mass, float), 0, None))
        assert r[0] > 0 and np.all(np.diff(r) > 0) and mass[-1] > 0
        self.lr, self.m = np.log(r), mass
        self.r0, self.r1, self.mtot = r[0], r[-1], mass[-1]
        integrand = G*mass/r                      # dPhi/dln r
        steps = (integrand[1:] + integrand[:-1])/2*np.diff(self.lr)
        outside = np.concatenate([np.cumsum(steps[::-1])[::-1], [0.]])
        self.phi_grid = -G*self.mtot/self.r1 - outside

    def mass(self, r):
        r = np.asarray(r, float)
        inner = np.interp(np.log(np.maximum(r, 1e-300)), self.lr, self.m)
        return np.where(r < self.r0, self.m[0]*(r/self.r0)**3, np.where(r > self.r1, self.mtot, inner))

    def phi(self, r):
        r = np.asarray(r, float)
        inner = np.interp(np.log(np.maximum(r, 1e-300)), self.lr, self.phi_grid)
        core = self.phi_grid[0] - G*self.m[0]*(self.r0**2 - r*r)/(2*self.r0**3)
        return np.where(r < self.r0, core, np.where(r > self.r1, -G*self.mtot/np.maximum(r, 1e-300), inner))

    def vesc(self, r):
        return np.sqrt(-2*self.phi(r))

    def vcirc(self, r):
        r = np.asarray(r, float)
        return np.sqrt(G*self.mass(r)/r)


def apsides(pot, E, L, r_start, iterations=60):
    """Pericenter and apocenter of bound orbits (E<0) passing through r_start."""
    E, L, r_start = (np.asarray(v, float) for v in np.broadcast_arrays(E, L, r_start))
    assert np.all(E < 0)

    def F(lr):
        r = np.exp(lr)
        return 2*(E - pot.phi(r)) - (L/r)**2

    lo, hi = np.log(r_start) - 45., np.log(r_start)
    for _ in range(iterations):              # F(lo)<0 at the angular barrier
        mid = (lo + hi)/2
        neg = F(mid) < 0
        lo, hi = np.where(neg, mid, lo), np.where(neg, hi, mid)
    rp = np.where(L > 0, np.exp(hi), 0.)
    # Phi >= -G Mtot/r everywhere, so F<0 beyond G Mtot/|E|.
    lo = np.log(r_start)
    hi = np.maximum(np.log(1.01*G*pot.mtot/np.abs(E)), lo + 1e-12)
    for _ in range(iterations):
        mid = (lo + hi)/2
        pos = F(mid) >= 0
        lo, hi = np.where(pos, mid, lo), np.where(pos, hi, mid)
    return rp, np.exp(lo)


def circular_radius(pot, E, iterations=60):
    """Radius of the circular orbit with energy E."""
    E = np.asarray(E, float)
    lo = np.full(E.shape, np.log(pot.r0) - 20.)
    hi = np.maximum(np.log(1.01*G*pot.mtot/np.abs(E)), lo + 1.)
    for _ in range(iterations):
        mid = (lo + hi)/2
        r = np.exp(mid)
        low = pot.phi(r) + G*pot.mass(r)/(2*r) < E
        lo, hi = np.where(low, mid, lo), np.where(low, hi, mid)
    return np.exp((lo + hi)/2)


def _tables(pot, E, L, rp, ra, n_theta):
    """Radial phase tables using r = rm - dr cos(theta), which removes the
    turning-point singularity of dt = dr/|v_r|."""
    rm, dr = (ra + rp)/2, (ra - rp)/2
    t = (np.arange(n_theta) + .5)*np.pi/n_theta
    rr = rm[:, None] - dr[:, None]*np.cos(t)
    F = 2*(E[:, None] - pot.phi(rr)) - (L[:, None]/rr)**2
    h = dr[:, None]*np.sin(t)/np.sqrt(np.maximum(F, 1e-300))
    circular = (dr < 1e-9*rm) | ~np.all(np.isfinite(h), axis=1)
    h[circular] = 1.     # harmonic limit: uniform in theta
    half = h.sum(axis=1)*np.pi/n_theta
    return rr, h, half, F, circular


def population_profile(pot, r_inj, E, L, weight, r_eval, n_theta=96, chunk=4096, bin_edges=None):
    """Steady phase-mixed profile of injected orbits in a fixed spherical potential.

    Each injected orbit contributes its time-averaged enclosed-mass fraction.
    Returns enclosed mass on r_eval, orbit diagnostics and, when bin_edges is
    given, time-weighted radial/tangential velocity moments in radial bins.
    """
    r_eval = np.asarray(r_eval, float)
    rp, ra = apsides(pot, E, L, r_inj)
    mass = np.zeros(len(r_eval))
    period = np.empty(len(E))
    nb = 0 if bin_edges is None else len(bin_edges) - 1
    bm, bvr, bvt = np.zeros(nb), np.zeros(nb), np.zeros(nb)
    for s in range(0, len(E), chunk):
        sl = slice(s, s + chunk)
        rr, h, half, F, circular = _tables(pot, E[sl], L[sl], rp[sl], ra[sl], n_theta)
        C = np.concatenate([np.zeros((len(h), 1)), np.cumsum(h, axis=1)], axis=1)/(h.sum(axis=1, keepdims=True))
        rm, dr = (ra[sl] + rp[sl])/2, (ra[sl] - rp[sl])/2
        arg = np.clip((rm[:, None] - r_eval[None, :])/np.where(dr > 0, dr, 1.)[:, None], -1, 1)
        pos = np.arccos(arg)*n_theta/np.pi
        i = np.clip(np.floor(pos).astype(int), 0, n_theta - 1)
        f = pos - i
        rows = np.arange(len(rm))[:, None]
        P = C[rows, i]*(1 - f) + C[rows, i + 1]*f
        P[circular] = (r_eval[None, :] >= rm[circular, None])
        mass += weight[sl] @ P
        period[sl] = np.where(circular, np.nan, 2*half)
        if nb:
            w = (weight[sl]/h.sum(axis=1))[:, None]*h
            vt2 = (L[sl, None]/rr)**2
            vr2 = np.maximum(F, 0)
            vr2[circular] = 0
            k = np.digitize(rr.ravel(), bin_edges) - 1
            ok = (k >= 0) & (k < nb)
            bm += np.bincount(k[ok], w.ravel()[ok], nb)
            bvr += np.bincount(k[ok], (w*vr2).ravel()[ok], nb)
            bvt += np.bincount(k[ok], (w*vt2).ravel()[ok], nb)
    Lc = np.empty(len(E))
    for s in range(0, len(E), chunk):
        rc = circular_radius(pot, E[s:s + chunk])
        Lc[s:s + chunk] = rc*pot.vcirc(rc)
    out = dict(r=r_eval, mass=mass, rp=rp, ra=ra, period_gyr=period*GYR, circularity=np.clip(L/Lc, 0, 1))
    if nb:
        with np.errstate(invalid='ignore', divide='ignore'):
            sr2, st2 = bvr/bm, bvt/bm
            out.update(bin_edges=np.asarray(bin_edges), bin_mass=bm, sigma_r2=sr2, sigma_t2=st2,
                       beta=1 - st2/(2*sr2))
    return out


def uniform_ball_injections(pot, r_sites, site_weight, n_speed=24, n_angle=12, s_split=.95, e_min=1e-6):
    """Bound products uniform in the galaxy-frame velocity ball |v|<v_esc at each site.

    Site weight is the local reaction rate; the bound rate per site scales with
    the ball volume, v_esc^3. Speeds s=v/v_esc use two Gauss panels: s in
    [0, s_split] for the bulk and ln(1-s^2) in [ln e_min, ln(1-s_split^2)] for
    the weakly bound orbits that set the outer tail. The omitted shell
    1-s^2<e_min carries a fraction of about 1.5 e_min. Radial direction cosine
    w in [0,1] covers both signs of v_r, which give the same (E, L).
    """
    nb = max(2, int(round(2*n_speed/3)))
    x, xw = leggauss(nb)
    s1 = s_split*(x + 1)/2
    w1 = 3*s1**2*xw*s_split/2                          # int 3 s^2 ds
    y, yw = leggauss(n_speed - nb)
    lo, hi = np.log(e_min), np.log(1 - s_split**2)
    e = np.exp((hi - lo)*(y + 1)/2 + lo)
    s2 = np.sqrt(1 - e)
    w2 = 1.5*s2*e*yw*(hi - lo)/2                       # 3 s^2 ds = (3/2) s e dln e
    s, sw = np.r_[s1, s2], np.r_[w1, w2]
    n_speed = len(s)
    y, yw = leggauss(n_angle)
    w, ww = (y + 1)/2, yw/2
    r0 = np.repeat(np.asarray(r_sites, float), n_speed*n_angle)
    ve = pot.vesc(r0)
    S = np.tile(np.repeat(s, n_angle), len(r_sites))
    W = np.tile(np.tile(w, n_speed), len(r_sites))
    weight = (np.repeat(np.asarray(site_weight, float)*pot.vesc(np.asarray(r_sites, float))**3/3, n_speed*n_angle)
              * np.tile(np.outer(sw, ww).ravel(), len(r_sites)))
    E = pot.phi(r0)*(1 - S*S)
    L = r0*ve*S*np.sqrt(1 - W*W)
    return dict(r=r0, E=E, L=L, weight=weight, speed_fraction=S)


def summarize(r, mass, radii):
    """Mass radii, fractions beyond given radii and the outer logarithmic slope."""
    frac = mass/mass[-1]
    half, ninety = (float(np.interp(q, frac, r)) for q in (.5, .9))
    lr, lm = np.log(r), np.log(np.maximum(np.gradient(mass, r)/(4*np.pi*r*r), 1e-300))
    outer = (frac > .99) & (frac < .9999)
    slope = float(np.polyfit(lr[outer], lm[outer], 1)[0]) if outer.sum() > 3 else float('nan')
    return dict(r50=half, r90=ninety, outer_density_slope=slope,
                fraction_outside={str(v): float(1 - np.interp(v, r, frac)) for v in radii})


# Explicitly labeled support closures -------------------------------------------------

def plummer_K(M, a):
    """Polytropic constant of the isolated n=5 Plummer solution, P=K rho^(6/5)."""
    return G/6*(4*np.pi/3)**.2*M**.8*a**-.4


def plummer_scale_for_K(M, K):
    """Scale radius at fixed K: a = [G (4 pi/3)^(1/5) M^(4/5)/(6 K)]^(5/2), so a ~ M^2."""
    return (G/6*(4*np.pi/3)**.2*M**.8/K)**2.5


def plummer_state(M, a, r):
    """Density, potential and pressure of the isolated self-gravitating Plummer sphere."""
    rho = 3*M/(4*np.pi*a**3)*(1 + (r/a)**2)**-2.5
    phi = -G*M/np.sqrt(r*r + a*a)
    return dict(rho=rho, phi=phi, pressure=rho*G*M/(6*np.sqrt(r*r + a*a)))
