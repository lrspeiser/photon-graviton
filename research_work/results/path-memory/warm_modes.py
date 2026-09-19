"""RUT-1 stage 8: the linear modes of a warm annulus coupled to the two-stage memory response (protocol-rut8.md).

A perturbation delta C_m(r) e^{i m theta + s t} of the memory field is a potential perturbation -delta C for
the bodies. In action-angle variables of the unperturbed potential Phi0 = -GM/r - C0(r) the linearised
collisionless response, harmonic by harmonic in the radial angle, is

    delta f_l(J) = F_l(J) psi_l(J) / (nu_l(J) - i s),   nu_l = l Omega_r + m Omega_theta,
    F_l = nu_l df0/dE + m df0/dL|_E,

and the memory field answers the density locally in the Laplace variable,

    delta C_m(r) = alpha H(s) 2 pi int r' dr' k_m(r, r') delta Sigma_m(r'),
    H(s) = tau_keep / [(1 + s tau_keep)(1 + s tau_form)],   k_m(r, r') = exp[-(r-r')^2/2w^2] I_m^e(r r'/w^2).

Expanding delta C_m in the kernel's own harmonics k_m(r_b, .) centred on radial nodes, and collocating there,

    T(s) c = 0,   T(s) = P + alpha H(s) M(s),   P_ab = k_m(r_a, r_b),
    M_ab(s) = (2 pi)^2 sum_l int dJ_r dL  G_al(J) G_bl(J) F_l(J) / (nu_l(J) - i s),

with G_al the orbit transform of k_m(r_a, .). It is solved in the orthonormal basis of P's well-represented
subspace, which has the same roots and a determinant of order one. T is analytic for Re s > 0 and its roots
there are the growing modes. A DAMPED mode is invisible to this search: no root above a floor is not stability.

The contour primitives (the rectangle quadrature's Gauss-Legendre panels, the refined winding number, the
bounded Newton polish) are stage 7's, imported unchanged. The moment routine is NOT: `spectrum.located()`
takes the determinant on the contour with the 2x2 formula whatever the size of T, which stage 7 never
exercised and which returns the winding number of the wrong function here. `located()` below is this
stage's own, with the full determinant.
"""
import numpy as np
from scipy.fft import dct
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import ive

import spectrum as SP

GM = 1.
T0 = 2*np.pi


class Disk:
    """The unperturbed state, from a population's archived state: potential, DF parameters, response."""

    def __init__(self, state, tau_form=3*T0, field=True):
        s = state
        self.L0, self.dL, self.dE, self.w = s['L0'], s['dL'], s['dE'], s['w']
        self.alpha, self.tau_keep, self.tau_form, self.reach = s['alpha'], s['tau_keep'], float(tau_form), s['reach']
        self.amplitude, self.mass = s['amplitude'], s['mass']
        self.r = np.asarray(s['r'], float)
        self.C = np.asarray(s['C'], float)*(1. if field else 0.)     # field=False is L1's negative control
        self.sigma = np.asarray(s['sigma'], float)
        self.Cs = CubicSpline(self.r, self.C)
        self.dCs = self.Cs.derivative()
        self.d2Cs = self.Cs.derivative(2)

    def phi(self, r):
        return -GM/r - self.Cs(r)

    def dphi(self, r):
        return GM/r**2 - self.dCs(r)

    def d2phi(self, r):
        return -2*GM/r**3 - self.d2Cs(r)

    def H(self, s):
        return self.tau_keep/((1 + s*self.tau_keep)*(1 + s*self.tau_form))

    def circular(self, L):
        """The global minimum of the effective potential: r_c, E_circ, Omega_c = dE_circ/dL, and kappa."""
        eff = self.phi(self.r) + .5*(L/self.r)**2
        i = int(np.clip(np.argmin(eff), 1, len(self.r) - 2))
        rc = brentq(lambda r: self.dphi(r) - L*L/r**3, self.r[i - 1], self.r[i + 1], xtol=1e-14)
        kappa2 = float(self.d2phi(rc) + 3*self.dphi(rc)/rc)
        return rc, float(self.phi(rc) + .5*(L/rc)**2), L/rc**2, np.sqrt(kappa2) if kappa2 > 0 else np.nan


class Orbits:
    """Every populated orbit on a Gauss-Legendre grid in (L, u), E = E_circ(L) + u^2/2, with its frequencies,
    its radial angle and its azimuthal lag along the orbit, to spectral accuracy in an eccentric-type angle."""

    def __init__(self, disk, n_L, n_u, n_eta):
        d = self.disk = disk
        xl, wl = np.polynomial.legendre.leggauss(n_L)
        L, wL = d.L0 + d.reach*d.dL*xl, d.reach*d.dL*wl
        u_max = np.sqrt(2*d.reach*d.dE)
        xu, wu = np.polynomial.legendre.leggauss(n_u)
        u, wU = .5*u_max*(xu + 1), .5*u_max*wu
        fine = np.linspace(d.r[0], d.r[-1], 4000)
        phi_fine = d.phi(fine)
        rows, circ = [], []
        self.narrow = self.two_wells = 0
        for Li, wLi in zip(L, wL):
            rc, Ec, Oc, kap = d.circular(Li)
            circ.append((Li, rc, Ec, Oc, kap))
            for uj, wUj in zip(u, wU):
                E = Ec + .5*uj*uj
                pos = 2*(E - phi_fine) - (Li/fine)**2 > 0
                if pos[0] or pos[-1]:
                    raise RuntimeError('a populated orbit reaches the edge of the radial domain')
                edges = np.nonzero(np.diff(pos.astype(int)))[0]
                Qf = lambda r: 2*(E - d.phi(r)) - (Li/r)**2
                spans = [(fine[edges[k]], fine[edges[k] + 1], fine[edges[k + 1]], fine[edges[k + 1] + 1])
                         for k in range(0, len(edges), 2)]
                if not any(sp[0] <= rc <= sp[3] for sp in spans):
                    step = 1e-7                  # nearly circular, narrower than the search grid
                    while Qf(rc - step) > 0 or Qf(rc + step) > 0:
                        step *= 2.
                    spans.append((rc - step, rc, rc, rc + step))
                    self.narrow += 1
                self.two_wells += len(spans) > 1
                for lo0, lo1, hi0, hi1 in spans:
                    rows.append((Li, uj, E, Ec, Oc, brentq(Qf, lo0, lo1, xtol=1e-15), brentq(Qf, hi0, hi1, xtol=1e-15),
                                 wLi*wUj*uj))                       # dE = u du
        self.L, self.u, self.E, self.Ec, self.Oc, self.rp, self.ra, self.wt = np.array(rows).T
        self.circular_rows = np.array(circ)
        self.n = len(rows)
        eta = (np.arange(n_eta) + .5)*np.pi/n_eta
        mid, half = .5*(self.ra + self.rp), .5*(self.ra - self.rp)
        self.rk = mid[:, None] - half[:, None]*np.cos(eta)[None, :]
        Q = 2*(self.E[:, None] - d.phi(self.rk)) - (self.L[:, None]/self.rk)**2
        q = Q/((self.rk - self.rp[:, None])*(self.ra[:, None] - self.rk))
        if not np.all(q > 0):
            raise RuntimeError('the radial speed is not real along a populated orbit')
        dt = 1./np.sqrt(q)                                       # dt / d eta
        dth = (self.L[:, None]/self.rk**2)*dt                    # d theta / d eta
        ca, cb = dct(dt, type=2, axis=1)/n_eta, dct(dth, type=2, axis=1)/n_eta
        a0, b0 = .5*ca[:, 0], .5*cb[:, 0]
        self.Om_r, self.Om_th = 1./a0, b0/a0
        j = np.arange(1, n_eta)
        S = np.sin(j[None, :]*eta[:, None])/j[None, :]
        t_osc, th_osc = ca[:, 1:] @ S.T, cb[:, 1:] @ S.T
        self.wr = eta[None, :] + self.Om_r[:, None]*t_osc        # radial angle along the orbit
        self.lag = th_osc - self.Om_th[:, None]*t_osc            # theta - w_theta
        self.dwr = self.Om_r[:, None]*dt                         # d w_r / d eta
        x, y = (self.L - d.L0)/d.dL, .5*self.u**2/d.dE
        self.f0 = d.amplitude*np.exp(-.5*x*x - .5*y*y)
        self.dfdE = -(y/d.dE)*self.f0
        self.dfdL = (-x/d.dL + (y/d.dE)*self.Oc)*self.f0
        self.measure = (2*np.pi)**2*self.wt/self.Om_r            # (2 pi)^2 dJ_r dL

    def transform(self, g_of_r, m, harmonics):
        """G_l(J) = (1/pi) int_0^pi dw_r g(r) cos(l w_r - m lag): shape (orbits, harmonics)."""
        g = g_of_r(self.rk)*self.dwr
        return np.stack([np.mean(g*np.cos(l*self.wr - m*self.lag), axis=1) for l in harmonics], axis=1)

    def moment(self, g_of_r):
        """int g Sigma_0 d^2x through the orbits."""
        return float(np.sum(self.measure*self.f0*np.mean(g_of_r(self.rk)*self.dwr, axis=1)))


class Modes:
    """T(s) for one azimuthal number on one population."""

    def __init__(self, disk, orbits, m, l_max, spacing, margin, clip, rank_cutoff, harmonics=None, prune=1e-15):
        self.d, self.o, self.m = disk, orbits, int(m)
        self.l = np.arange(-l_max, l_max + 1) if harmonics is None else np.asarray(harmonics)
        lo = max(clip[0], float(orbits.rp.min()) - margin)
        hi = min(clip[1], float(orbits.ra.max()) + margin)
        self.nodes = lo + spacing*np.arange(int(np.floor((hi - lo)/spacing + 1e-9)) + 1)
        self.P = self.kernel(self.nodes[:, None], self.nodes[None, :])
        # Overlapping Gaussians are nearly dependent: cond(P) is 2e5 at a node spacing of w/2 and 5e8 at 3w/8.
        # The problem is posed instead in the orthonormal basis of P's well-represented subspace, Q^T P Q = 1,
        # which has the same roots, a determinant of order one, and a size set by the kernel's bandwidth and
        # not by how finely the nodes are laid down.
        lam, U = np.linalg.eigh(self.P)
        keep_basis = lam > rank_cutoff*lam.max()
        self.Q = U[:, keep_basis]/np.sqrt(lam[keep_basis])[None, :]
        self.condition = float(lam.max()/lam.min()) if lam.min() > 0 else float('inf')   # of P itself
        self.basis_size = int(keep_basis.sum())
        G = np.stack([orbits.transform(lambda r, ra=ra: self.kernel(ra, r), self.m, self.l) for ra in self.nodes],
                     axis=2).reshape(orbits.n*len(self.l), len(self.nodes))
        nu = (self.l[None, :]*orbits.Om_r[:, None] + self.m*orbits.Om_th[:, None]).ravel()
        F = (nu.reshape(orbits.n, -1)*orbits.dfdE[:, None] + self.m*orbits.dfdL[:, None]).ravel()
        FW = F*np.repeat(orbits.measure, len(self.l))
        size = np.abs(FW)*np.sum(G*G, axis=1)
        keep = size > prune*size.max()                           # rows that cannot matter at double precision
        self.rows_total, self.rows_kept = int(len(keep)), int(keep.sum())
        self.G_nodes = np.ascontiguousarray(G[keep])             # in the node basis: what L2 compares with
        self.G, self.nu, self.FW = np.ascontiguousarray(self.G_nodes @ self.Q), nu[keep], FW[keep]

    def kernel(self, ra, r):
        w = self.d.w
        return np.exp(-(r - ra)**2/(2*w*w))*ive(self.m, r*ra/(w*w))

    def M(self, s):
        """The response matrix in the NODE basis, M_ab(s) of the module docstring."""
        return (self.G_nodes*(self.FW/(self.nu - 1j*s))[:, None]).T @ self.G_nodes

    def T(self, s):
        """Q^T [P + alpha H(s) M(s)] Q = 1 + alpha H(s) Q^T M(s) Q."""
        reduced = (self.G*(self.FW/(self.nu - 1j*s))[:, None]).T @ self.G
        return np.eye(self.basis_size) + self.d.alpha*self.d.H(s)*reduced

    def T_many(self, S):
        return np.stack([self.T(s) for s in np.atleast_1d(np.asarray(S, complex))])


# ---------------------------------------------------------------- the contour
def contour(rect, left_panel, first_panel, growth, right_panel, order=16):
    """Counter-clockwise rectangle. The left edge runs beside the imaginary axis, where the discretised
    response has its poles, so its panels are short and uniform; the horizontal edges are graded, short at
    the left corner and growing geometrically away from it."""
    re_lo, re_hi, im_lo, im_hi = rect
    xg, wg = np.polynomial.legendre.leggauss(order)

    def graded(length):
        cuts, step = [0.], first_panel
        while cuts[-1] + step < length:
            cuts.append(cuts[-1] + step)
            step *= growth
        return np.array(cuts + [length])

    def uniform(length, panel):
        return np.linspace(0., length, max(int(np.ceil(length/panel - 1e-9)), 1) + 1)

    z, dz = [], []

    def edge(a, b, cuts):
        span = abs(b - a)
        unit = (b - a)/span
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            mid, half = .5*(lo + hi), .5*(hi - lo)
            z.append(a + unit*(mid + half*xg))
            dz.append(unit*half*wg)

    h = graded(re_hi - re_lo)
    edge(complex(re_lo, im_lo), complex(re_hi, im_lo), h)                                   # bottom, left to right
    edge(complex(re_hi, im_lo), complex(re_hi, im_hi), uniform(im_hi - im_lo, right_panel))  # right, upward
    edge(complex(re_hi, im_hi), complex(re_lo, im_hi), (re_hi - re_lo) - h[::-1])            # top, right to left
    edge(complex(re_lo, im_hi), complex(re_lo, im_lo), uniform(im_hi - im_lo, left_panel))   # left, downward
    return np.concatenate(z), np.concatenate(dz)


def located(T_many, rect, quad, probes=4, moments=8, seed=0, noise=1e-11):
    """Beyn's moments on the contour, the FULL determinant there, and the winding number from the same
    evaluations. Rank is chosen at the largest gap of the singular spectrum; the caller decides whether that
    gap is large enough to mean anything."""
    z, dz = contour(rect, **quad)
    T = T_many(z)
    n = T.shape[-1]
    centre = complex(.5*(rect[0] + rect[1]), .5*(rect[2] + rect[3]))
    rho = .5*abs(complex(rect[1] - rect[0], rect[3] - rect[2]))
    rng = np.random.default_rng(seed)
    ell = min(probes, n)
    V = rng.standard_normal((n, ell)) + 1j*rng.standard_normal((n, ell))
    X = np.linalg.solve(T, np.broadcast_to(V, (len(z), n, ell)))
    zeta = (z - centre)/rho
    K = moments//2
    A, pw = np.empty((moments, n, ell), complex), np.ones(len(z), complex)
    for p in range(moments):
        A[p] = np.einsum('s,sab->ab', dz*pw, X)/(2j*np.pi*rho)
        pw = pw*zeta
    B0 = np.block([[A[i + j] for j in range(K)] for i in range(K)])
    B1 = np.block([[A[i + j + 1] for j in range(K)] for i in range(K)])
    U, sv, Wh = np.linalg.svd(B0)
    # as in stage 7: singular values under the quadrature's own noise floor are not eigenvalues, so a gap
    # among them is not a rank. The winding number, which does not depend on this, is the safeguard.
    floor = noise*float(np.sum(np.abs(dz)*np.linalg.norm(X, axis=(1, 2)))/(2*np.pi*rho))
    live = int(np.sum(sv > floor))
    if live == 0:
        k, gap = 1, 0.
    else:
        ratios = sv[:live]/np.maximum(np.append(sv[1:live], floor), 1e-300)
        k = int(np.argmax(ratios)) + 1
        gap = float(ratios[k - 1])
    red = U[:, :k].conj().T @ B1 @ Wh[:k].conj().T @ np.diag(1./sv[:k])
    raw = centre + rho*np.linalg.eigvals(red)
    dets = np.linalg.det(T)
    return dict(raw=raw, rank=k, gap=gap, capacity=int(K*ell), singular_values=sv, contour=z, det_on_contour=dets,
                nodes=int(len(z)), median_det=float(np.median(np.abs(dets))))


def solve_rectangle(modes, rect, L, quad, determinant=None):
    """Locate, polish, then test; and count independently. `determinant` exists for L5's negative control."""
    det_many = (lambda zz: np.linalg.det(modes.T_many(zz))) if determinant is None else \
               (lambda zz: determinant(modes.T_many(zz)))
    loc = located(modes.T_many, rect, quad)
    values = loc['det_on_contour'] if determinant is None else determinant(modes.T_many(loc['contour']))
    wind = SP.winding(det_many, loc['contour'], values, L['winding_phase_step'])
    clear = loc['gap'] >= L['singular_gap']
    roots = []
    if clear:
        det1 = lambda s: np.linalg.det(modes.T(s))
        inside = lambda zz: rect[0] < zz.real < rect[1] and rect[2] < zz.imag < rect[3]
        for s in loc['raw']:
            to_edge = min(abs(s.real - rect[0]), abs(s.real - rect[1]), abs(s.imag - rect[2]), abs(s.imag - rect[3]))
            sharp = SP.polish(det1, s, max(min(10*L['polish_movement'], to_edge), 1e-9))
            if inside(sharp) and not any(abs(sharp - q['s']) < 1e-7 for q in roots):
                roots.append(dict(s=complex(sharp), raw=complex(s), moved=float(abs(sharp - s)),
                                  relative_residual=float(abs(det1(sharp))/loc['median_det'])))
    return dict(rect=[float(x) for x in rect], roots=roots, located=len(roots), rank=loc['rank'] if clear else 0,
                gap=loc['gap'], clear_rank=bool(clear), winding=wind, counts_agree=bool(len(roots) == wind['count']),
                saturated=bool(clear and loc['rank'] >= loc['capacity']), nodes=loc['nodes'],
                singular_values=[float(x) for x in loc['singular_values'][:8]])


def strip_rectangles(m, L, shifted=False):
    """The declared strip for one m, cut into rectangles of the declared height; the second partition is
    offset by half a rectangle, so the two share no cut."""
    lo = -(L['strip_im_per_m']*m + L['strip_im_offset'])
    hi = L['strip_im_retrograde']
    h = L['rectangle_height']
    start = lo - (.5*h if shifted else 0.)
    count = int(np.ceil((hi - start)/h - 1e-9))
    return [(L['strip_re'][0], L['strip_re'][1], start + k*h, start + (k + 1)*h) for k in range(count)]
