"""CL-2 stage 2 library: the reconstructed Newtonian disk force, the exact equal-galaxy speed loss with two
solvers and a certificate, the release's SZ correlation matrices, and the corrected boundary fit.
Provenance: ring potential and force by complete elliptic integrals (established); the uniform-sheet
subtraction is derived here from established mathematics; the convex speed loss is the review's identity."""
import io
import tarfile
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, nnls
from scipy.special import ellipk, ellipe, i0, i1, k0, k1
from scipy.linalg import cholesky, solve_triangular
import cl2_sources as CS
import cl2_response as CR

G = CS.G
WIDTHS = CR.WIDTHS[CR.WIDTHS >= CR.GALAXY_MIN_WIDTH]          # 21 widths, 0.158 kpc to 10 Mpc, every block


# ------------------------------------------------------------------ S2-3: the razor-thin disk's Newtonian force
def ring_force(m, a, R):
    """Inward in-plane force of a thin ring of mass m and radius a at radius R (R != a): the R-derivative of the
    ring potential Phi = -(2 G m/pi) K(p)/(R+a), p = 4 R a/(R+a)^2 (established)."""
    R, a = np.asarray(R, float), np.asarray(a, float)
    q = ((a - R)/(R + a))**2                      # 1 - p, computed without cancellation
    p = 1 - q
    K, E = ellipk(p), ellipe(p)
    with np.errstate(divide='ignore', invalid='ignore'):
        dK = (E - q*K)/(2*p*q)
        dp = 4*a*(a - R)/(R + a)**3
        out = -(2*G*m/np.pi)*(dK*dp/(R + a) - K/(R + a)**2)
    return np.where(q > 0, out, 0.)              # the ring's own radius is a set of measure zero in the integrals


def disk_force(sigma, R, kinks=(), Rmax=None, rmin=1e-6):
    """g(R) = int 2 pi s [Sigma(s) - Sigma(R)] F(R,s) ds - Sigma(R) int_{Rmax}^inf 2 pi s F(R,s) ds, using the
    identity that an infinite uniform sheet exerts no in-plane force. kinks: radii where Sigma has slope breaks."""
    R = np.atleast_1d(np.asarray(R, float))
    out = np.zeros_like(R)
    kinks = np.asarray(sorted(set(float(k) for k in kinks)), float)
    for i, x in enumerate(R):
        Rm = Rmax if Rmax is not None else max(60.*x, 400.)
        S0 = float(sigma(np.array([x]))[0])
        reg = lambda s: 2*np.pi*s*(float(sigma(np.array([s]))[0]) - S0)*ring_force(1., s, x)
        lo_pts = [k for k in kinks if rmin < k < x]
        hi_pts = [k for k in kinks if x < k < Rm]
        I1 = quad(reg, rmin, x, points=lo_pts or None, limit=800)[0] + quad(reg, x, Rm, points=hi_pts or None, limit=800)[0]
        tail = lambda s: 2*np.pi*s*ring_force(1., s, x)
        I2 = -(quad(tail, Rm, 20*Rm, limit=800)[0] + quad(tail, 20*Rm, np.inf, limit=800)[0])
        out[i] = I1 + S0*I2
    return out


def plain_ring_sum(sigma, R, n=8000, Rmax=400.):
    """The control: rings between grid nodes summed directly; it does not converge near the target radius."""
    Rg = np.geomspace(1e-3, Rmax, n)
    Sig = sigma(Rg)
    Rm = np.sqrt(Rg[1:]*Rg[:-1])
    mass = np.pi*(Rg[1:]**2 - Rg[:-1]**2)*.5*(Sig[1:] + Sig[:-1])
    return np.array([np.sum(ring_force(mass, Rm, x)) for x in np.atleast_1d(R)])


def freeman_force(Sig0, Rd, R):
    y = np.asarray(R, float)/(2*Rd)
    return 4*np.pi*G*Sig0*Rd*y*y*(i0(y)*k0(y) - i1(y)*k1(y))/np.asarray(R, float)


def disk_force_gates():
    Sig0, Rd = 5e8, 2.
    R = np.array([.3, .5, 1., 2., 4., 8., 16., 30.])
    sig = lambda s: Sig0*np.exp(-np.asarray(s, float)/Rd)
    ref = freeman_force(Sig0, Rd, R)
    g = disk_force(sig, R)
    ctrl = plain_ring_sum(sig, R)
    # a compact ring set (uniform annulus 0.9-1.1 kpc) at ten times its radius against the point mass
    ann = lambda s: np.where((np.asarray(s, float) > .9) & (np.asarray(s, float) < 1.1), 1e9, 0.)
    mass = np.pi*(1.1**2 - .9**2)*1e9
    pm = disk_force(ann, np.array([10.]), kinks=(.9, 1.1))[0]
    # the bulge: a Plummer sphere's enclosed mass gives G M(<r)/r^2 exactly by construction of the block builder
    return dict(freeman_relative=float(np.max(np.abs(g/ref - 1))), control_relative=float(np.max(np.abs(ctrl/ref - 1))),
                point_mass_relative=float(abs(pm/(G*mass/100.) - 1)), radii_kpc=R.tolist(), freeman=ref.tolist(), rings=g.tolist())


def reconstructed_newton(gal, upsilon_disk=CS.I.UPSILON_DISK):
    """The Newtonian force of the same components that source the written field: razor-thin stellar and gas
    disks by the ring route, the spherical bulge by G M(<r)/r^2; plus the integrated baryonic mass."""
    comp = CS.BAR.sparc_components(gal)
    R, vobs, ev, gN_tab = CS._galaxy_arrays(gal, upsilon_disk)
    rows = gal['rotmod']
    use = (rows[:, 0] > 0) & np.isfinite(rows[:, 6]) & (rows[:, 6] > 0)
    kinks = rows[use, 0]
    sigma = lambda s: comp['sigma_star'](s) + (comp['sigma_gas'](s) if comp['sigma_gas'] is not None else 0.)
    g = disk_force(sigma, R, kinks=kinks)
    if comp['m_bulge'] is not None:
        g = g + G*comp['m_bulge'](R)/R**2
    return dict(R=R, vobs=vobs, ev=ev, gN_reconstructed=g, gN_tabulated=gN_tab, mass=float(CS.BAR.sparc_total_mass(comp)))


# ------------------------------------------------------------------ S2-4: the exact equal-galaxy speed loss
class SpeedLoss:
    """F(L) = mean over galaxies of the mean over points of (sqrt(u) - v)^2, u = R (gN + A L); convex on u >= 0.
    Internally the amplitudes are scaled column by column (x_j = L_j c_j with c_j the column's largest force per
    unit amplitude) so that both solvers work on O(1) variables; every public value is in the original units."""
    def __init__(self, A, gN, R, vobs, sizes, weights=None):
        self.A, self.gN, self.R, self.v = np.asarray(A, float), np.asarray(gN, float), np.asarray(R, float), np.asarray(vobs, float)
        self.sizes = list(sizes)
        self.w = np.concatenate([np.full(n, 1./(n*len(self.sizes))) for n in self.sizes]) if weights is None else np.asarray(weights, float)
        self.c = np.maximum(np.max(np.abs(self.A), axis=0), 1e-300)
        self.As = self.A/self.c[None, :]

    # ---- the objective in original units
    def u(self, L):
        return self.R*(self.gN + self.A@L)

    def value(self, L):
        u = self.u(L)
        if np.any(u < 0):
            return np.inf
        return float(np.sum(self.w*(np.sqrt(u) - self.v)**2))

    def gradient(self, L):
        u = self.u(L)
        su = np.sqrt(np.maximum(u, 1e-300))
        return self.A.T@(self.w*(1 - self.v/su)*self.R)

    def hessian(self, L):
        u = np.maximum(self.u(L), 1e-300)
        cc = self.w*self.v/(2*u**1.5)*self.R**2
        return (self.A*cc[:, None]).T@self.A

    def rmse(self, L):
        return float(np.sqrt(self.value(L)))

    def per_galaxy_rmse(self, L):
        res = (np.sqrt(np.maximum(self.u(L), 0)) - self.v)**2
        out, off = [], 0
        for n in self.sizes:
            out.append(float(np.sqrt(np.mean(res[off:off + n]))))
            off += n
        return out

    # ---- scaled versions for the solvers
    def _fx(self, x):
        return self.value(x/self.c)

    def _gx(self, x):
        return self.gradient(x/self.c)/self.c

    def _hx(self, x):
        return self.hessian(x/self.c)/np.outer(self.c, self.c)

    # solver A: L-BFGS-B on the scaled variables, the objective infinite outside the domain
    def solve_lbfgs(self, L0, maxiter=50000):
        def f(x):
            val = self._fx(x)
            if not np.isfinite(val):
                return 1e300, np.zeros_like(x)
            return val, self._gx(x)
        x0 = np.maximum(np.asarray(L0, float)*self.c, 0)
        res = minimize(f, x0, jac=True, method='L-BFGS-B', bounds=[(0, None)]*len(x0),
                       options=dict(maxiter=maxiter, ftol=1e-16, gtol=1e-14, maxcor=50, maxls=100))
        return res.x/self.c, res

    # solver B: a trust-region Newton method with the analytic Hessian (SciPy's trust-constr), an independent code path
    def solve_trust_newton(self, L0=None, fallback=None, maxiter=5000):
        from scipy.optimize import Bounds
        n = self.A.shape[1]
        x = np.zeros(n) if L0 is None else np.maximum(np.asarray(L0, float)*self.c, 0)
        start = 'origin' if L0 is None else 'given'
        if not np.isfinite(self._fx(x)) or not np.all(self.u(x/self.c) > 0):
            for label, cand in (('fallback', None if fallback is None else np.maximum(np.asarray(fallback, float)*self.c, 0)), ('small_positive', np.full(n, 1e-9))):
                if cand is not None and np.isfinite(self._fx(cand)) and np.all(self.u(cand/self.c) > 0):
                    x, start = cand, label
                    break
            else:
                raise ValueError('trust-region Newton: no feasible start')
        big = 1e300

        def f(z):
            v = self._fx(z)
            return v if np.isfinite(v) else big
        res = minimize(f, x, jac=self._gx, hess=self._hx, method='trust-constr', bounds=Bounds(np.zeros(n), np.full(n, np.inf)),
                       options=dict(gtol=1e-13, xtol=1e-16, maxiter=maxiter, verbose=0))
        return res.x/self.c, dict(iterations=int(res.nit), final=self.value(res.x/self.c), start=start, message=str(res.message))

    solve_projected_newton = solve_trust_newton
    solve_active_set = solve_trust_newton

    def projected_gradient(self, L):
        g = self._gx(np.asarray(L, float)*self.c)
        x = np.asarray(L, float)*self.c
        return float(np.linalg.norm(np.where(x > 0, g, np.minimum(g, 0))))

    def kkt(self, L, scale_point):
        """Optimality in the scaled variables: gradient components relative to the gradient's size at the start point."""
        x = np.asarray(L, float)*self.c
        g = self._gx(x)
        scale = float(np.max(np.abs(self._gx(np.maximum(np.asarray(scale_point, float)*self.c, 0)))))
        act = x > 0
        return dict(max_abs_gradient_active=float(np.max(np.abs(g[act]))/scale) if act.any() else 0.,
                    min_gradient_inactive=float(np.min(g[~act])/scale) if (~act).any() else 0., scale=scale)

    def curvature_check(self, L, n_dirs=200, seed=5, eps=1e-6):
        """Numerical Hessian-vector products in the scaled variables along random feasible directions (components at
        the bound are not moved below it), the smallest divided by the largest so the reading is scale-free."""
        rng = np.random.default_rng(seed)
        x = np.asarray(L, float)*self.c
        worst, best = np.inf, -np.inf
        for _ in range(n_dirs):
            d = rng.normal(size=len(x))
            d[(x <= 0) & (d < 0)] = 0.
            if not np.any(d):
                continue
            d /= np.linalg.norm(d)
            h = eps*max(1., np.linalg.norm(x))
            neg = d < 0
            if neg.any():
                h = min(h, .5*float(np.min(x[neg]/(-d[neg]))))
            if h <= 0:
                continue
            xp, xm = x + h*d, x - h*d
            if np.any(xm < 0):
                xm = x
                curv = np.dot(d, self._gx(xp) - self._gx(xm))/h
            else:
                curv = np.dot(d, self._gx(xp) - self._gx(xm))/(2*h)
            worst, best = min(worst, curv), max(best, curv)
        return float(worst/max(best, 1e-300))


def nnls_start(A, y, s, weights):
    """Stage 1's acceleration problem: nonnegative least squares on (y - A L)/s with per-point weights."""
    Aw = A*np.sqrt(weights)[:, None]/s[:, None]
    yw = y*np.sqrt(weights)/s
    x, _ = nnls(Aw, yw, maxiter=100000)
    return x


# ------------------------------------------------------------------ S2-6: the SZ correlation matrices and the corrected boundary fit
def sz_correlations(tar_path, names):
    """The release's per-cluster pressure correlation matrix (COVMAT block, fourth HDU of the Y-PROF-COVMAT files)
    and the block's radii, in kpc."""
    from astropy.io import fits
    out = {}
    with tarfile.open(tar_path) as t:
        members = {m.name: m for m in t.getmembers()}
        for name in names:
            f = [k for k in members if k.startswith(name + '/') and 'COVMAT' in k][0]
            h = fits.open(io.BytesIO(t.extractfile(members[f]).read()))
            d = h[4].data[0]
            C = np.array(d['COVMAT'], float)
            sd = np.sqrt(np.diag(C))
            out[name] = dict(corr=C/np.outer(sd, sd), radius_kpc=np.array(d['RW'], float), file=f)
    return out


def sz_submatrix(cl, corr_record, tol=.02):
    """The correlation sub-matrix of the SZ points the block retains (those inside the cluster's outer radius),
    matched to the release's block by radius within tol."""
    kinds = np.asarray(cl['kinds'])
    sz = np.flatnonzero(kinds == 'sz')
    corr, radii = corr_record['corr'], corr_record.get('radius_kpc')
    if radii is None or len(corr) == len(sz):
        return sz, corr[:len(sz), :len(sz)] if radii is None else corr
    idx = []
    for r in cl['rp'][sz]:
        j = int(np.argmin(np.abs(radii - r)))
        if abs(radii[j] - r)/r > tol:
            raise ValueError('SZ point at %.1f kpc has no release radius within %.0f%%' % (r, 100*tol))
        idx.append(j)
    return sz, corr[np.ix_(idx, idx)]


def whitened_cluster(cl, corr_record, ops):
    """Rows of one cluster with the SZ covariance corr o (err err^T) on the retained SZ points and independent X-ray
    errors: returns (A_w, nuisance_w, y_w, Lc) such that chi2 = |A_w L + nu_w P_out - y_w|^2."""
    err = cl['eP']
    n = len(err)
    Cov = np.diag(err**2)
    sz, corr = sz_submatrix(cl, corr_record)
    Cov[np.ix_(sz, sz)] = corr*np.outer(err[sz], err[sz])
    Lc = cholesky(Cov, lower=True)
    A_w = solve_triangular(Lc, ops, lower=True)
    nu_w = solve_triangular(Lc, cl['thermal'][:, None], lower=True)
    y_w = solve_triangular(Lc, cl['Pobs'] - cl['P_N'], lower=True)
    return A_w, nu_w, y_w, Lc


def boundary_fit(P_model, thermal, Pobs, Lc):
    """chi2 of a fixed model (thermal fraction already applied inside P_model) with a nonnegative boundary
    pressure under the covariance Lc Lc^T: the corrected fit_boundary."""
    r = solve_triangular(Lc, Pobs - P_model, lower=True)
    t = solve_triangular(Lc, thermal, lower=True)
    pout = max(float(np.dot(t, r)/np.dot(t, t)), 0.)
    return float(np.sum((r - pout*t)**2)), pout


# ------------------------------------------------------------------ S2-5: the lens rows (the well coupling, TF-1's construction restated here to avoid a module-name clash)
from scipy.optimize import minimize_scalar   # noqa: E402

LENS_GEOMETRY = 'PF1_static_euclidean'
CHI2_LIMIT, EINSTEIN_LIMIT, EINSTEIN_WEIGHT = 128., .03, 1e3


class LensRows:
    """One lens through CL-2's LensSystem under the static geometry: whitened kinematics rows and the Einstein row,
    linear in the amplitudes at fixed anisotropy; the well coupling (light factor two, stars the full force)."""
    def __init__(self, name, scenario=LENS_GEOMETRY, widths=WIDTHS):
        self.name = name
        self.L = CS.LensSystem(name, scenario, widths=widths)
        D = np.diag(2*self.L.y)
        self.Lc = cholesky(D@self.L.cov@D, lower=True)
        self.cache = {}
        self.n_bins = len(self.L.y)

    def coefficients(self, beta):
        key = round(float(beta), 12)
        if key not in self.cache:
            self.L.model.forces = np.vstack([self.L.starforce, self.L.memforce])
            self.cache[key] = self.L.model.coefficients(beta)
        return self.cache[key]

    def rows(self, imf, beta):
        m = self.L.pop[imf]/1e11
        co = self.coefficients(beta)
        A = (m*co[1:]).T
        y = self.L.y**2 - m*co[0]
        return dict(Ak=solve_triangular(self.Lc, A, lower=True), yk=solve_triangular(self.Lc, y, lower=True),
                    ae=m*self.L.membend, ye=self.L.need - m*self.L.starbend, se=EINSTEIN_LIMIT*self.L.need, m=m, co=co)

    def evaluate(self, imf, beta, amps):
        r = self.rows(imf, beta)
        v2 = r['m']*(r['co'][0] + np.asarray(amps)@r['co'][1:])
        return dict(chi2=self.L.chi2_v(v2), einstein_residual=float((r['ae']@amps - r['ye'])/self.L.need), beta=float(beta))

    def fit_beta(self, imf, amps):
        ev = lambda b: self.evaluate(imf, b, amps)['chi2']
        lo, hi = self.L.beta_bounds
        grid = np.linspace(lo, hi, 11)
        vals = [ev(b) for b in grid]
        i = int(np.argmin(vals))
        res = minimize_scalar(ev, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 10)]), method='bounded', options={'xatol': 1e-6})
        return (float(res.fun), float(res.x)) if res.fun < vals[i] else (float(vals[i]), float(grid[i]))
