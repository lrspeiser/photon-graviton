"""RUT-1 stage 9: the corrected linear response of a SHARPLY TRUNCATED warm annulus.

Stage 8's `warm_modes.py` is imported unchanged and is not edited: its archive is pinned to its bytes. It
differentiates the population's Gaussian inside the declared bounds and omits what the bounds themselves
contribute (found by the owner in the review of 3a80fec). The population that is built, drawn and simulated is

    f0(E, L) = g(E, L) Theta(E_max(L) - E) Theta(L - L_lo) Theta(L_hi - L),
    g = A exp[-(L-L0)^2/2dL^2] exp[-(E-E_c(L))^2/2dE^2],  E_max = E_c(L) + reach dE,  L_lo,hi = L0 -+ reach dL,

so that, with dE_max/dL = dE_c/dL = Omega_c(L),

    df0/dE      = g_E Theta... - g delta(E - E_max) ...,
    df0/dL |_E  = g_L Theta... + g Omega_c delta(E - E_max) ... + g [delta(L - L_lo) - delta(L - L_hi)] ...,

and F_l = nu_l df0/dE + m df0/dL|_E gains three families of boundary rows,

    energy edge   (m Omega_c - nu_l) g(E_max(L), L)      on the orbits (L, E_max(L)),     measure (2 pi)^2 dL / Omega_r,
    lower L edge  + m g(E, L_lo)                         on the orbits (L_lo, E),         measure (2 pi)^2 dE / Omega_r,
    upper L edge  - m g(E, L_hi)                         on the orbits (L_hi, E),         measure (2 pi)^2 dE / Omega_r.

The lower energy bound E = E_c(L) is not a truncation: it is J_r = 0, which no perturbation carries a body
across, and g_E vanishes there.

`OrbitSet` is stage 8's orbit construction for an arbitrary tensor of (L, u) nodes; on stage 8's grid it
reproduces `warm_modes.Orbits` exactly, which the checks verify. `Response` is stage 8's `Modes` with the
boundary rows, a kernel-node rule that is an argument and not a constant, and the left eigenvector needed for
an eigenvalue sensitivity.
"""
import numpy as np
from scipy.fft import dct
from scipy.optimize import brentq
from scipy.special import ive

import warm_modes as WM

GM = WM.GM
T0 = WM.T0


def gauss_legendre(lo, hi, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return .5*(hi - lo)*x + .5*(hi + lo), .5*(hi - lo)*w


def stage8_grids(disk, n_L, n_u):
    """The (L, u) quadrature of stage 8, in its own arithmetic: L0 + reach dL x, and u_max (x + 1)/2."""
    xl, wl = np.polynomial.legendre.leggauss(n_L)
    xu, wu = np.polynomial.legendre.leggauss(n_u)
    u_max = np.sqrt(2*disk.reach*disk.dE)
    return (disk.L0 + disk.reach*disk.dL*xl, disk.reach*disk.dL*wl), (.5*u_max*(xu + 1), .5*u_max*wu)


class OrbitSet:
    """Orbits on a tensor of (L, u) nodes, E = E_c(L) + u^2/2, with their frequencies, radial angle and
    azimuthal lag: `warm_modes.Orbits` with the nodes as arguments. `energy_measure` multiplies the weight
    by u (dE = u du); an energy-edge family, which is not integrated over energy, turns it off."""

    def __init__(self, disk, L_grid, u_grid, n_eta, energy_measure=True):
        d = self.disk = disk
        (L, wL), (u, wU) = L_grid, u_grid
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
                                 wLi*wUj*uj if energy_measure else wLi*wUj))
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
        self.x, self.y = (self.L - d.L0)/d.dL, .5*self.u**2/d.dE
        self.g = d.amplitude*np.exp(-.5*self.x*self.x - .5*self.y*self.y)   # the Gaussian, bounds aside
        self.measure = (2*np.pi)**2*self.wt/self.Om_r

    def transform(self, g_of_r, m, harmonics):
        """G_l(J) = (1/pi) int_0^pi dw_r g(r) cos(l w_r - m lag): shape (orbits, harmonics)."""
        g = g_of_r(self.rk)*self.dwr
        return np.stack([np.mean(g*np.cos(l*self.wr - m*self.lag), axis=1) for l in harmonics], axis=1)

    def moment(self, g_of_r, f0=None):
        return float(np.sum(self.measure*(self.g if f0 is None else f0)*np.mean(g_of_r(self.rk)*self.dwr, axis=1)))


# ---------------------------------------------------------------- what each family of rows contributes to F_l
def interior(o, nu, m):
    """nu g_E + m g_L|_E inside the bounds: stage 8's F_l, in stage 8's arithmetic."""
    d = o.disk
    dfdE = -(o.y/d.dE)*o.g
    dfdL = (-o.x/d.dL + (o.y/d.dE)*o.Oc)*o.g
    return nu*dfdE[:, None] + m*dfdL[:, None]


def edge_energy(o, nu, m):
    return (m*o.Oc[:, None] - nu)*o.g[:, None]


def edge_L_lower(o, nu, m):
    return m*o.g[:, None]*np.ones_like(nu)


def edge_L_upper(o, nu, m):
    return -m*o.g[:, None]*np.ones_like(nu)


def sharp_families(disk, n_L, n_u, n_eta, edges=True):
    """The interior on stage 8's quadrature and, with `edges`, the three boundary families on the same nodes."""
    Lg, ug = stage8_grids(disk, n_L, n_u)
    fam = [(OrbitSet(disk, Lg, ug, n_eta), interior)]
    if edges:
        one = np.ones(1)
        u_max = np.sqrt(2*disk.reach*disk.dE)
        fam.append((OrbitSet(disk, Lg, (u_max*one, one), n_eta, energy_measure=False), edge_energy))
        fam.append((OrbitSet(disk, ((disk.L0 - disk.reach*disk.dL)*one, one), ug, n_eta), edge_L_lower))
        fam.append((OrbitSet(disk, ((disk.L0 + disk.reach*disk.dL)*one, one), ug, n_eta), edge_L_upper))
    return fam


def kernel_nodes(extent, spacing, margin, clip, shift=0.):
    """Kernel centres `spacing` apart from `margin` inside the populated orbits to `margin` outside them;
    `shift`, in units of the spacing, moves every centre and is a convergence variant."""
    lo = max(clip[0], extent[0] - margin) + shift*spacing
    hi = min(clip[1], extent[1] + margin)
    return lo + spacing*np.arange(int(np.floor((hi - lo)/spacing + 1e-9)) + 1)


class Response:
    """T(s) for one azimuthal number: families of (orbits, F) rows on declared kernel nodes."""

    def __init__(self, disk, families, m, l_max, nodes, rank_cutoff, harmonics=None, prune=1e-15):
        self.d, self.m = disk, int(m)
        self.families = families
        self.o = families[0][0]
        self.l = np.arange(-l_max, l_max + 1) if harmonics is None else np.asarray(harmonics)
        self.nodes = np.asarray(nodes, float)
        self.P = self.kernel(self.nodes[:, None], self.nodes[None, :])
        lam, U = np.linalg.eigh(self.P)
        keep_basis = lam > rank_cutoff*lam.max()
        self.Q = U[:, keep_basis]/np.sqrt(lam[keep_basis])[None, :]
        self.condition = float(lam.max()/lam.min()) if lam.min() > 0 else float('inf')
        self.basis_size = int(keep_basis.sum())
        Gs, nus, FWs, ls, self.family_rows = [], [], [], [], []
        for orbits, F_of in families:
            G = np.stack([orbits.transform(lambda r, ra=ra: self.kernel(ra, r), self.m, self.l) for ra in self.nodes],
                         axis=2).reshape(orbits.n*len(self.l), len(self.nodes))
            nu = self.l[None, :]*orbits.Om_r[:, None] + self.m*orbits.Om_th[:, None]
            FW = (F_of(orbits, nu, self.m)*orbits.measure[:, None]).ravel()
            Gs.append(G)
            nus.append(nu.ravel())
            FWs.append(FW)
            ls.append(np.tile(self.l, orbits.n))
            self.family_rows.append(len(FW))
        G, nu, FW = np.concatenate(Gs), np.concatenate(nus), np.concatenate(FWs)
        self.family_of_row = np.repeat(np.arange(len(families)), self.family_rows)
        size = np.abs(FW)*np.sum(G*G, axis=1)
        keep = size > prune*size[:self.family_rows[0]].max()     # the threshold is the interior's, as in stage 8
        self.rows_total, self.rows_kept = int(len(keep)), int(keep.sum())
        self.keep = keep
        self.family_of_row, self.l_of_row = self.family_of_row[keep], np.concatenate(ls)[keep]
        self.G_nodes = np.ascontiguousarray(G[keep])
        self.G, self.nu, self.FW = np.ascontiguousarray(self.G_nodes @ self.Q), nu[keep], FW[keep]

    def kernel(self, ra, r):
        w = self.d.w
        return np.exp(-(r - ra)**2/(2*w*w))*ive(self.m, r*ra/(w*w))

    def M(self, s, family=None):
        """The response matrix in the NODE basis; `family` restricts it to one family of rows."""
        sel = slice(None) if family is None else self.family_of_row == family
        G = self.G_nodes[sel]
        return (G*(self.FW[sel]/(self.nu[sel] - 1j*s))[:, None]).T @ G

    def reduced(self, s, family=None, derivative=False, harmonic=None):
        """Q^T M(s) Q, or its derivative in s; `family` and `harmonic` restrict the rows."""
        sel = np.ones(len(self.nu), bool)
        if family is not None:
            sel &= self.family_of_row == family
        if harmonic is not None:
            sel &= self.l_of_row == harmonic
        G, den = self.G[sel], self.nu[sel] - 1j*s
        return (G*(self.FW[sel]*(1j/den**2 if derivative else 1./den))[:, None]).T @ G

    def T(self, s):
        return np.eye(self.basis_size) + self.d.alpha*self.d.H(s)*self.reduced(s)

    def T_many(self, S):
        return np.stack([self.T(s) for s in np.atleast_1d(np.asarray(S, complex))])

    def dT_ds(self, s):
        d = self.d
        H = d.H(s)
        dH = -H*(d.tau_keep/(1 + s*d.tau_keep) + d.tau_form/(1 + s*d.tau_form))
        return d.alpha*(dH*self.reduced(s) + H*self.reduced(s, derivative=True))

    def eigenvectors(self, s):
        """Right and left null vectors of T(s) at a root, from the smallest singular triple."""
        U, sv, Vh = np.linalg.svd(self.T(s))
        return Vh[-1].conj(), U[:, -1], float(sv[-1]/sv[0])        # T v = 0,  u^H T = 0

    def shift(self, s, delta_T):
        """First-order movement of the root at s under a change delta_T of T: -(u^H dT v)/(u^H T' v)."""
        v, u, _ = self.eigenvectors(s)
        return complex(-(u.conj() @ delta_T @ v)/(u.conj() @ self.dT_ds(s) @ v))

    def field_coefficients(self, s):
        """The mode's field in the node basis: delta C_m(r) = sum_b c_b k_m(r_b, r), c = Q v."""
        v, _, ratio = self.eigenvectors(s)
        return self.Q @ v, ratio

    def off_node_residual(self, s, points):
        """The mode equation at radii that are NOT collocation nodes: the field the response regenerates,
        -alpha H(s) [M c](r), against the field imposed, sum_b c_b k_m(r, r_b), over the largest imposed value.
        Collocation makes this vanish on the nodes whatever the basis; between them it measures the basis."""
        v, _, _ = self.eigenvectors(s)
        c = self.Q @ v
        points = np.atleast_1d(np.asarray(points, float))
        Gt = np.concatenate([np.stack([orbits.transform(lambda r, rt=rt: self.kernel(rt, r), self.m, self.l)
                                       for rt in points], axis=2).reshape(orbits.n*len(self.l), len(points))
                             for orbits, _ in self.families])[self.keep]
        weight = self.FW/(self.nu - 1j*s)
        regenerated = -self.d.alpha*self.d.H(s)*(Gt.T @ (weight*(self.G_nodes @ c)))
        imposed = self.kernel(points[:, None], self.nodes[None, :]) @ c
        return np.abs(regenerated - imposed)/np.max(np.abs(imposed))
