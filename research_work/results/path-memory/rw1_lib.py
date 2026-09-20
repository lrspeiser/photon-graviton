"""RW-1 library: the whirlpool beyond the rut (protocol-rw1.md).
Provenance: the Plummer-cored power-law reach f(s) = s (s^2 + w^2)^(-(p+1)/2) r_*^(p-2) is the owner's proposal (p = 2 is
Plummer's softened Newtonian force and p = 1 the force of the two-dimensional logarithmic potential, both established);
the ring average is 128-node Gauss-Legendre quadrature (established); the shell average's closed form and binomial series
are elementary integrals derived here; the fits use stage 2's equal-galaxy speed loss, Lawson-Hanson NNLS and stage 2's
lens rows with the three-dimensional bend (RPG-1's quadrature)."""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import poch, ellipk, ellipe
from scipy.optimize import minimize_scalar, nnls
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.linalg import solve_triangular
import cl2_sources as CS
import cl2s2_lib as L2
import steady_field as SF

G, C_KMS = SF.G, SF.C_KMS
R_STAR = 1.                                   # kpc: lambda is the pull relative to Newton's at 1 kpc from a point mass
P_GRID = (.5, .75, 1., 1.25, 1.5, 1.75, 2.)
W_GAL = (.25, .5, 1., 2., 4., 8., 16.)
W_CL = (1., 2., 4., 8., 16., 32., 64., 128., 256.)
KERNELS_GAL = [(p, w) for p in P_GRID for w in W_GAL]
KERNELS_CL = [(p, w) for p in P_GRID for w in W_CL]
N128, N256 = leggauss(128), leggauss(256)
SWITCH, TERMS = .2, 40


# ------------------------------------------------------------------ the kernels
def whirl_ring(R, Rp, p, w, nodes=N128):
    """Inward pull at R (in the plane) from a unit-mass coplanar ring at Rp: the azimuthal average of
    (R - Rp cos phi) (s^2 + w^2)^(-(p+1)/2), s^2 = R^2 + Rp^2 - 2 R Rp cos phi. Shape (len(R), len(Rp))."""
    x, wt = nodes
    c = np.cos(np.pi*(x + 1)/2)[None, None, :]
    wh = (wt/2)[None, None, :]
    R = np.atleast_1d(np.asarray(R, float))[:, None, None]
    Rp = np.asarray(Rp, float)[None, :, None]
    s2 = R*R + Rp*Rp - 2*R*Rp*c
    return ((R - Rp*c)*wh*(s2 + w*w)**(-(p + 1)/2)).sum(-1)


def whirl_shell(r, rp, p, w, switch=SWITCH, terms=TERMS):
    """Inward pull at r from a unit-mass spherical shell at rp: (1/2) int_{-1}^{1} (r - rp mu) (C - B mu)^(-a) dmu with
    C = r^2 + rp^2 + w^2, B = 2 r rp, a = (p+1)/2. Closed form (two power integrals, differences through expm1 and
    log1p) where x = B/C > switch; the binomial series in x elsewhere. Shape (len(r), len(rp))."""
    r = np.atleast_1d(np.asarray(r, float))[:, None]
    rp = np.asarray(rp, float)[None, :]
    a = (p + 1)/2
    shape = np.broadcast(r, rp).shape
    R, P = np.broadcast_to(r, shape), np.broadcast_to(rp, shape)
    C = R*R + P*P + w*w
    B = 2*R*P
    x = B/C
    out = np.empty(shape)
    hi = x > switch
    if hi.any():
        rr, pp, Bh, Ch = R[hi], P[hi], B[hi], C[hi]
        u0 = (rr - pp)**2 + w*w
        q = np.log1p(2*Bh/u0)
        D1 = q if abs(a - 1) < 1e-12 else u0**(1 - a)*np.expm1((1 - a)*q)/(1 - a)
        D2 = u0**(2 - a)*np.expm1((2 - a)*q)/(2 - a)
        out[hi] = ((rr*Bh - pp*Ch)*D1 + pp*D2)/(2*Bh*Bh)
    lo = ~hi
    if lo.any():
        rr, pp, xx, CC = R[lo], P[lo], x[lo], C[lo]
        s = np.zeros_like(rr)
        xk = np.ones_like(rr)
        for k in range(terms):
            coef = poch(a, k)/math.factorial(k)
            s = s + coef*xk*(rr/(k + 1) if k % 2 == 0 else -pp/(k + 2))
            xk = xk*xx
        out[lo] = CC**(-a)*s
    return out


def shell_apply(r, nodes_r, dM, p, w, chunk=1000):
    """G r_*^(p-2) sum_shells dM K(r, r'): the inward pull per unit lambda of spherical shells, chunked over r."""
    r = np.atleast_1d(np.asarray(r, float))
    out = np.empty(len(r))
    for i in range(0, len(r), chunk):
        out[i:i + chunk] = whirl_shell(r[i:i + chunk], nodes_r, p, w)@dM
    return G*R_STAR**(p - 2)*out


def shell_p1_exact(r, rp):
    """Shell average of the 1/s pull (p = 1, w = 0): 1/(2r) + (r^2 - rp^2)/(4 r^2 rp) ln((r + rp)/|r - rp|)."""
    r = np.atleast_1d(np.asarray(r, float))[:, None]
    rp = np.asarray(rp, float)[None, :]
    return 1/(2*r) + (r*r - rp*rp)/(4*r*r*rp)*np.log((r + rp)/np.abs(r - rp))


def ring_newton_exact(R, Rp):
    """In-plane inward pull of a unit ring at Rp with the 1/s^2 kernel (p = 2, w = 0): the derivative of the ring
    potential -(2/pi) K(k)/(R + Rp), k^2 = 4 R Rp/(R + Rp)^2."""
    R = np.atleast_1d(np.asarray(R, float))[:, None]
    Rp = np.asarray(Rp, float)[None, :]
    k2 = 4*R*Rp/(R + Rp)**2
    K, E = ellipk(k2), ellipe(k2)
    dk2 = 4*Rp*(Rp - R)/(R + Rp)**3
    dK = (E/(k2*(1 - k2)) - K/k2)/2*dk2
    return -(2/np.pi)*(dK/(R + Rp) - K/(R + Rp)**2)


def kernel_gates():
    """K1 to K6 of the protocol on synthetic rings and shells; returns the magnitudes and the pass flags."""
    out = {}
    Rp = np.array([1.])
    ratios = np.array([.3, .7, .9, 1.1, 1.5, 3., 30.])
    k = whirl_ring(ratios, Rp, 1., 1e-6)[:, 0]
    out['K1_ring_p1_gauss'] = float(np.max(np.abs(k - np.where(ratios > 1, 1/ratios, 0.))))
    k = whirl_ring(ratios, Rp, 2., 1e-6)[:, 0]
    ex = ring_newton_exact(ratios, Rp)[:, 0]
    out['K2_ring_p2_elliptic'] = float(np.max(np.abs(k - ex)/np.abs(ex)))
    k = whirl_shell(ratios, Rp, 2., 1e-6)[:, 0]
    out['K3_shell_p2_theorem'] = float(np.max(np.abs(k - np.where(ratios > 1, 1/ratios**2, 0.))))
    k = whirl_shell(ratios, Rp, 1., 1e-6)[:, 0]
    ex = shell_p1_exact(ratios, Rp)[:, 0]
    out['K4_shell_p1_closed'] = float(np.max(np.abs(k - ex)/np.abs(ex)))
    worst = 0.
    for p in P_GRID:
        a = (p + 1)/2
        for w in (0., .01, .3, 3., 300.):
            for r in (1e-9, 1e-6, 1e-3, .05, .3, .5, .8, .99, 1., 1.01, 1.2, 2., 10., 1e3, 1e6):
                if w == 0. and r == 1.:
                    continue
                f = lambda mu: (r - mu)*(r*r + 1 + w*w - 2*r*mu)**(-a)/2
                q = quad(f, -1, 1, limit=800, epsabs=0, epsrel=2e-14)[0]
                kk = whirl_shell(np.array([r]), Rp, p, w)[0, 0]
                scale = max(abs(q), (r*r + 1 + w*w)**(-a)*min(r, 1.))
                worst = max(worst, abs(kk - q)/scale)
    out['K5_shell_vs_quad'] = float(worst)
    worst = 0.
    r = np.geomspace(1e-3, 1e3, 4000)
    for p in P_GRID:
        for w in (.01, .3, 3.):
            x = 2*r/(r*r + 1 + w*w)
            m = (x > .15) & (x < .3)
            k1 = whirl_shell(r[m], Rp, p, w, switch=.05)[:, 0]
            k2 = whirl_shell(r[m], Rp, p, w, switch=.5)[:, 0]
            worst = max(worst, float(np.max(np.abs(k1 - k2)/np.abs(k2))))
    out['K5_branch_agreement'] = worst
    worst = 0.
    for p in P_GRID:
        ks = whirl_shell(np.array([100.]), Rp, p, 1.)[0, 0]
        kr = whirl_ring(np.array([100.]), Rp, p, 1.)[0, 0]
        worst = max(worst, abs(ks*100.**p - 1), abs(kr*100.**p - 1))
    out['K6_far_field'] = float(worst)
    out['passes'] = dict(K1=out['K1_ring_p1_gauss'] < 1e-9, K2=out['K2_ring_p2_elliptic'] < 1e-8, K3=out['K3_shell_p2_theorem'] < 1e-8,
                         K4=out['K4_shell_p1_closed'] < 1e-9, K5=out['K5_shell_vs_quad'] < 1e-6 and out['K5_branch_agreement'] < 1e-8,
                         K6=out['K6_far_field'] < 1e-3)
    return out


# ------------------------------------------------------------------ galaxy columns
def galaxy_source(gal, n_grid=3000, r_max=300.):
    """Ring masses of the razor-thin stellar and gas disks on the geometric grid, the bulge's shell masses, the
    half-mass radius of the whole and the disk scale length."""
    comp = CS.BAR.sparc_components(gal)
    Rg = np.geomspace(1e-3, r_max, n_grid)
    sig = comp['sigma_star'](Rg) + (comp['sigma_gas'](Rg) if comp['sigma_gas'] is not None else 0.)
    Rm = np.sqrt(Rg[1:]*Rg[:-1])
    mass = np.pi*(Rg[1:]**2 - Rg[:-1]**2)*.5*(sig[1:] + sig[:-1])
    cum = np.concatenate([[0.], np.cumsum(mass)])
    dMb = None
    if comp['m_bulge'] is not None:
        Mb = comp['m_bulge'](Rg)
        dMb = np.diff(Mb)
        cum = cum + Mb
    R_half = float(np.interp(.5*cum[-1], cum, Rg))
    return dict(Rm=Rm, mass=mass, dMb=dMb, R_half=R_half, M_grid=float(cum[-1]), rd=float(comp['rd']))


def whirl_columns_disk(R, Rm, mass, kernels, nodes=N128, chunk=24):
    """Inward pull per unit lambda at the radii R from coplanar rings of masses `mass` at Rm, one column per kernel."""
    x, wt = nodes
    c = np.cos(np.pi*(x + 1)/2)[None, None, :]
    wh = (wt/2)[None, None, :]
    R = np.asarray(R, float)
    out = np.zeros((len(R), len(kernels)))
    Rp = np.asarray(Rm, float)[None, :, None]
    for i in range(0, len(R), chunk):
        RR = R[i:i + chunk][:, None, None]
        s2 = RR*RR + Rp*Rp - 2*RR*Rp*c
        num = (RR - Rp*c)*wh
        for k, (p, w) in enumerate(kernels):
            K = (num*(s2 + w*w)**(-(p + 1)/2)).sum(-1)
            out[i:i + chunk, k] = G*R_STAR**(p - 2)*(K@mass)
    return out


def galaxy_whirl_columns(gal, R, kernels=KERNELS_GAL, n_grid=3000, r_max=300., nodes=N128, point_source=False):
    """The whirlpool pull per unit lambda at the observed radii: the disks through the ring average, the bulge through
    the shell average; or, for the N2 control, the whole mass on the grid placed at the centre."""
    src = galaxy_source(gal, n_grid, r_max)
    R = np.asarray(R, float)
    if point_source:
        cols = np.array([G*R_STAR**(p - 2)*src['M_grid']*R*(R*R + w*w)**(-(p + 1)/2) for (p, w) in kernels]).T
        return cols, src
    cols = whirl_columns_disk(R, src['Rm'], src['mass'], kernels, nodes)
    if src['dMb'] is not None:
        cols = cols + np.array([shell_apply(R, src['Rm'], src['dMb'], p, w) for (p, w) in kernels]).T
    return cols, src


# ------------------------------------------------------------------ cluster columns
def cluster_whirl_columns(cl, kernels=KERNELS_CL, n_eval=1200, full=False):
    """Pressure columns (the hydrostatic integral of n_e g, stage 2's construction) of the whirlpool pull per unit lambda
    sourced by the cluster's gas and stars: g on a 1200-point evaluation grid, cubic spline in ln r to the working grid,
    or the full-grid evaluation for the K8 gate. Returns (ops, g on the grid)."""
    r, src = cl['r'], cl['source']
    re = r if full else np.geomspace(r[0], r[-1], n_eval)
    ops = np.zeros((len(cl['rp']), len(kernels)))
    gs = np.zeros((len(r), len(kernels)))
    for k, (p, w) in enumerate(kernels):
        ge = shell_apply(re, src.nodes, src.dM, p, w)
        g = ge if full else CubicSpline(np.log(re), ge)(np.log(r))
        gs[:, k] = g
        ops[:, k] = cl['pressure_of'](g)
    return ops, gs


def cluster_fit(col_w, nu_w, y_w):
    """One cluster, one kernel: min |col lambda + nu pout - y|^2 over lambda, pout >= 0 by NNLS; KKT residual relative
    to |M^T y|."""
    M = np.column_stack([col_w, nu_w[:, 0]])
    x, res = nnls(M, y_w)
    g = M.T@(M@x - y_w)
    act = x > 0
    kkt = max(float(np.max(np.abs(g[act]))) if act.any() else 0., max(0., -float(np.min(g[~act]))) if (~act).any() else 0.)
    return dict(lam=float(x[0]), pout=float(x[1]), chi2=float(res*res), kkt=kkt/max(float(np.linalg.norm(M.T@y_w)), 1e-300))


def cluster_fit_universal(cols_w, nus_w, ys_w):
    """One lambda for every cluster with its own boundary pressure: NNLS on the stacked whitened rows."""
    N = sum(len(y) for y in ys_w)
    M = np.zeros((N, 1 + len(cols_w)))
    off = 0
    for k, (c, nu) in enumerate(zip(cols_w, nus_w)):
        n = len(c)
        M[off:off + n, 0] = c
        M[off:off + n, 1 + k] = nu[:, 0]
        off += n
    y = np.concatenate(ys_w)
    x, res = nnls(M, y)
    g = M.T@(M@x - y)
    act = x > 0
    kkt = max(float(np.max(np.abs(g[act]))) if act.any() else 0., max(0., -float(np.min(g[~act]))) if (~act).any() else 0.)
    per = []
    off = 0
    for k, c in enumerate(cols_w):
        n = len(c)
        per.append(float(np.sum((M[off:off + n]@x - y[off:off + n])**2)))
        off += n
    return dict(lam=float(x[0]), pouts=x[1:].tolist(), chi2=float(res*res), per_cluster_chi2=per, kkt=kkt/max(float(np.linalg.norm(M.T@y)), 1e-300))


# ------------------------------------------------------------------ the certified one-parameter fit
def fit_lambda(W, gN, R, v, wts):
    """min over lambda >= 0 of sum wts (sqrt(R (gN + lambda W)) - v)^2 on its domain: bisection on the analytic
    derivative (solver A) checked by a bounded Brent minimisation of the value (solver B); F'' > 0 analytically."""
    W, gN, R, v, wts = (np.asarray(a, float) for a in (W, gN, R, v, wts))
    neg = W < 0
    lam_max = float(np.min(gN[neg]/(-W[neg]))) if neg.any() else np.inf
    top = lam_max*(1 - 1e-9) if np.isfinite(lam_max) else None

    def F(l):
        u = R*(gN + l*W)
        return np.inf if np.any(u <= 0) else float(np.sum(wts*(np.sqrt(u) - v)**2))

    def dF(l):
        u = R*(gN + l*W)
        return float(np.sum(wts*(1 - v/np.sqrt(u))*R*W))

    d0 = dF(0.)
    if d0 >= 0:
        lam, how = 0., 'bound'
    else:
        hi = top
        if hi is None:
            hi = 1.
            while dF(hi) < 0:
                hi *= 4.
        lo = 0.
        for _ in range(300):
            mid = .5*(lo + hi)
            if dF(mid) < 0:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 1e-14*hi:
                break
        lam, how = .5*(lo + hi), 'interior'
    FA = F(lam)
    upper = top if top is not None else max(4*lam, 1.)
    res = minimize_scalar(F, bounds=(0., upper), method='bounded', options=dict(xatol=1e-14*max(upper, 1e-30), maxiter=3000))
    FB = float(res.fun)
    u = R*(gN + lam*W)
    return dict(lam=float(lam), F=FA, rmse=float(np.sqrt(FA)), F_brent=FB, lam_brent=float(res.x), agree=abs(FA - FB)/max(FA, 1e-300),
                how=how, dF=dF(lam), dF0=d0, d2F=float(np.sum(wts*v*R*R*W*W/(2*u**1.5))), lam_max=lam_max)


def regime_masks(R, R_half, n_edge=5):
    R = np.asarray(R, float)
    inner = R < R_half
    edge = np.zeros(len(R), bool)
    edge[-min(n_edge, len(R)):] = True
    return dict(inner=inner, outer=~inner, edge=edge, whole=np.ones(len(R), bool))


# ------------------------------------------------------------------ the lenses
class LensWhirl(L2.LensRows):
    """A lens with whirlpool rows: for each kernel the inward pull per unit lambda per 1e11 Msun of the deprojected stars
    on the model grid (a 1000-point evaluation with the cubic spline in ln r, or the full grid), and its bend at b_E by the
    three-dimensional route with the r^-p tail beyond the grid."""
    def __init__(self, name, kernels=KERNELS_GAL, scenario=L2.LENS_GEOMETRY, n_eval=1000, full=False, cached=None):
        super().__init__(name, scenario, widths=np.array([1.]))
        L = self.L
        r = L.model.r
        self.kernels = list(kernels)
        if cached is not None:
            self.basis, self.basis_bend = cached
            return
        src = SF.SphericalSource(r, 1e11*L.frac)
        re = r if full else np.geomspace(r[0], r[-1], n_eval)
        rows = []
        for (p, w) in kernels:
            ge = shell_apply(re, src.nodes, src.dM, p, w)
            rows.append(ge if full else CubicSpline(np.log(re), ge)(np.log(r)))
        self.basis = np.array(rows)
        self.basis_bend = np.array([SF.deflection_from_g(self._tail(r, row, p), L.bE) for row, (p, w) in zip(self.basis, self.kernels)])

    @staticmethod
    def _tail(r, row, p):
        def g(x):
            return np.interp(x, r, row, left=row[0]) if x <= r[-1] else row[-1]*(r[-1]/x)**p
        return g

    def coefficients(self, beta):
        key = ('whirl', round(float(beta), 12))
        if key not in self.cache:
            self.L.model.forces = np.vstack([self.L.starforce, self.basis])
            self.cache[key] = self.L.model.coefficients(beta)
        return self.cache[key]

    def rows(self, imf, beta):
        m = self.L.pop[imf]/1e11
        co = self.coefficients(beta)
        A = (m*co[1:]).T
        y = self.L.y**2 - m*co[0]
        return dict(Ak=solve_triangular(self.Lc, A, lower=True), yk=solve_triangular(self.Lc, y, lower=True),
                    ae=m*self.basis_bend, ye=self.L.need - m*self.L.starbend, se=L2.EINSTEIN_LIMIT*self.L.need, m=m, co=co)

    def einstein_lambda(self, imf, k):
        """The strength of kernel k that meets the Einstein radius exactly (zero if the stars alone overshoot)."""
        m = self.L.pop[imf]/1e11
        return max(float((self.L.need - m*self.L.starbend)/(m*self.basis_bend[k])), 0.)


# ------------------------------------------------------------------ deflection curves and Coma
def deflection_curve(g_grid_r, g_grid, p_tail, bs):
    """alpha(b) by the three-dimensional route for a tabulated inward pull with an r^-p_tail beyond the grid (p_tail
    None: zero beyond the grid, the truncated Newtonian part)."""
    r, g = np.asarray(g_grid_r, float), np.asarray(g_grid, float)

    def fun(x):
        if x <= r[-1]:
            return float(np.interp(x, r, g, left=g[0]))
        return float(g[-1]*(r[-1]/x)**p_tail) if p_tail is not None else 0.
    return np.array([SF.deflection_from_g(fun, b) for b in bs])


def project_effective(r, g, R_out, n_R=800):
    """Sigma(R) and the projected mass M_p(<R) of the effective density of a spherical inward pull g(r) on the grid r
    (rho_eff = d(r^2 g/G)/dr/(4 pi r^2), any sign), on a geometric R grid to R_out."""
    r, g = np.asarray(r, float), np.asarray(g, float)
    Meff = r*r*g/G
    rho = np.gradient(Meff, r)/(4*np.pi*r*r)
    Rg = np.geomspace(r[0], R_out, n_R)
    Sigma = SF.project_density(r, rho, Rg)
    Mp = np.concatenate([[0.], np.cumsum(.5*(2*np.pi*Rg[1:]*Sigma[1:] + 2*np.pi*Rg[:-1]*Sigma[:-1])*np.diff(Rg))]) + np.pi*Rg[0]**2*Sigma[0]
    return Rg, Sigma, Mp


# ------------------------------------------------------------------ the Milky Way columns (the cache builder's construction)
def milky_way_columns(variant, R, kernels=KERNELS_GAL, n_grid=3000, r_max=300.):
    """Whirlpool pull per unit lambda at the radii R for baseline I or II: razor-thin disks (the archived projected
    surface density) through the ring average, the Plummer bulge through the shell average; returns (cols, R_half, M)."""
    sigma, m_bulge = CS.milky_way_source(variant)
    Rg = np.geomspace(1e-3, r_max, n_grid)
    sig = sigma(Rg)
    Rm = np.sqrt(Rg[1:]*Rg[:-1])
    mass = np.pi*(Rg[1:]**2 - Rg[:-1]**2)*.5*(sig[1:] + sig[:-1])
    cols = whirl_columns_disk(R, Rm, mass, kernels)
    cum = np.concatenate([[0.], np.cumsum(mass)])
    if m_bulge is not None:
        Mb = m_bulge(Rg)
        cols = cols + np.array([shell_apply(R, Rm, np.diff(Mb), p, w) for (p, w) in kernels]).T
        cum = cum + Mb
    return cols, float(np.interp(.5*cum[-1], cum, Rg)), float(cum[-1])


def regress_loglog(M, lam, seed=7, n_boot=1000):
    """log10 lam against log10 M by least squares over the systems with lam > 0; bootstrap 16th and 84th percentiles."""
    M, lam = np.asarray(M, float), np.asarray(lam, float)
    keep = lam > 0
    n = int(keep.sum())
    if n < 3:
        return dict(n=n, zero=int((~keep).sum()), slope=None, intercept=None)
    x, y = np.log10(M[keep]), np.log10(lam[keep])
    a, b = np.polyfit(x, y, 1)
    resid = y - (a*x + b)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.ptp(x[idx]) > 0:
            boots.append(np.polyfit(x[idx], y[idx], 1))
    boots = np.array(boots)
    return dict(n=n, zero=int((~keep).sum()), slope=float(a), intercept=float(b), scatter_dex=float(np.sqrt(np.mean(resid**2))),
                slope_16_84=[float(np.percentile(boots[:, 0], 16)), float(np.percentile(boots[:, 0], 84))],
                intercept_16_84=[float(np.percentile(boots[:, 1], 16)), float(np.percentile(boots[:, 1], 84))],
                log10_M_range=[float(x.min()), float(x.max())])
