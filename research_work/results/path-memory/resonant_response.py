"""RUT-1 stage 11: the resonant integral of the linear response, evaluated analytically cell by cell.

Stage 9's `warm_response.py` integrates

    M_ab(s) = (2 pi)^2 sum_l int dL du  W(L, u) G_al G_bl F_l / (nu_l(L, u) - i s)

with the Gauss-Legendre rule of stage 8. That rule replaces the continuum of orbital frequencies by a finite
set of them, so the discretised response has a POLE at s = -i nu for every node, all of them sitting on the
imaginary axis. As long as the growth rate Re s stays well above the spacing of those poles the sum is a good
quadrature; when Re s falls below it, it is not a quadrature of anything. That is stage 8's floor of 0.006, an
e-folding time of 26.5 reference periods, and the reason its search says nothing about slower modes.

This module removes the floor rather than lowering it. On each cell of a refined (L, u) grid the frequency and
the numerator are interpolated LINEARLY from the vertices, and the resulting integral is done in closed form.
For a triangle whose vertex frequencies are n1 <= n2 <= n3 the integral of each barycentric coordinate against
1/(nu_lin - z) needs only

    L0(x) = log(1+x)/x,   L1(x) = [x - log(1+x)]/x^2,   L2(x) = [log(1+x) - x + x^2/2]/x^3,

evaluated by their series where the closed forms cancel. The result is the exact integral of a piecewise-linear
model of the true integrand, so its error is second order in the cell size UNIFORMLY in Re s: the pole is
integrated, not sampled. Two refinements are carried and Richardson-extrapolated, and the difference between
the pair and a finer pair is the declared error estimate.

Only the harmonics that can resonate in the rectangle being searched are treated this way; the rest keep stage
9's rows, whose poles are then a declared margin away from the contour. Stage 9's `Response` supplies the
orbits, the kernel basis and the transform, unchanged and by import: this module never rebuilds them.

z = i s, so Im z = Re s. For Re s > 0 the integral above IS the Laplace transform, with no continuation, and
that is the only region this module is used in. A damped mode needs the continuation across the cut on
Re s = 0 and is not attempted here.
"""
import numpy as np

import warm_response as WR

_SMALL = .01
_TERMS = 8             # (0.01)^8 = 1e-16; at |x| = 0.01 the closed form of L2 still keeps 11 digits


def logs(x, want0=True):
    """L0, L1, L2 together: one Horner pass for the series, one log1p for the rest. `want0` drops L0, which
    the two-dimensional rule never uses."""
    x = np.asarray(x, complex)
    small = np.abs(x) < _SMALL
    mx = -np.where(small, x, 0.)
    a0 = np.zeros_like(mx) if want0 else None
    a1 = np.zeros_like(mx)
    a2 = np.zeros_like(mx)
    for k in range(_TERMS, -1, -1):
        if want0:
            a0 = a0*mx + 1./(k + 1)
        a1 = a1*mx + 1./(k + 2)
        a2 = a2*mx + 1./(k + 3)
    xl = np.where(small, 1., x)
    lg = np.log1p(xl)
    b1 = (xl - lg)/(xl*xl)
    b2 = (lg - xl + .5*xl*xl)/(xl**3)
    return (np.where(small, a0, lg/xl) if want0 else None), np.where(small, a1, b1), np.where(small, a2, b2)


class Triangles:
    """A tensor grid of vertices cut into triangles, each with one frequency per vertex. `weights(z)` returns
    the vector omega_v(z) = sum over triangles of int phi_v dA / (nu_lin - z), phi_v the barycentric hat."""

    def __init__(self, xa, xb, nu, pattern='alternate'):
        na, nb = len(xa), len(xb)
        idx = np.arange(na*nb).reshape(na, nb)
        v00, v10, v01, v11 = idx[:-1, :-1], idx[1:, :-1], idx[:-1, 1:], idx[1:, 1:]
        area = .5*np.outer(np.diff(xa), np.diff(xb))
        if pattern == 'alternate':                      # the diagonal alternates, so no direction is favoured
            flip = (np.add.outer(np.arange(na - 1), np.arange(nb - 1)) % 2).astype(bool)
        elif pattern == 'fixed':
            flip = np.zeros((na - 1, nb - 1), bool)
        elif pattern == 'other':
            flip = np.ones((na - 1, nb - 1), bool)
        else:
            raise ValueError('unknown triangulation %r' % (pattern,))
        t1 = np.where(flip[..., None], np.stack([v00, v10, v01], -1), np.stack([v00, v10, v11], -1))
        t2 = np.where(flip[..., None], np.stack([v10, v11, v01], -1), np.stack([v00, v11, v01], -1))
        tri = np.concatenate([t1.reshape(-1, 3), t2.reshape(-1, 3)])
        A = np.concatenate([area.ravel(), area.ravel()])
        nuv = np.asarray(nu, float).ravel()[tri]
        order = np.argsort(nuv, axis=1)                 # the closed form is written for n1 <= n2 <= n3
        self.tri = np.ascontiguousarray(np.take_along_axis(tri, order, axis=1))
        nus = np.take_along_axis(nuv, order, axis=1)
        self.n1, self.n3 = np.ascontiguousarray(nus[:, 0]), np.ascontiguousarray(nus[:, 2])
        d21, d31, d32 = nus[:, 1] - nus[:, 0], nus[:, 2] - nus[:, 0], nus[:, 2] - nus[:, 1]
        self.flat = d31 <= 0.                           # a cell the frequency is constant on
        safe = np.where(self.flat, 1., d31)
        self.A, self.r21, self.r32 = A, d21/safe, d32/safe
        self.nvert, self.ntri = na*nb, len(tri)

    def weights(self, z):
        D1, D3 = self.n1 - z, self.n3 - z
        _, a1, a2 = logs(self.r21*(self.n3 - self.n1)/D1, want0=False)
        _, b1, b2 = logs(-self.r32*(self.n3 - self.n1)/D3, want0=False)
        first, second = self.A*self.r21/D1, self.A*self.r32/D3
        W1 = 2*first*(a1 - .5*a2*(self.r21 + 1.)) + second*self.r32*b2
        W2 = first*a2 + second*b2
        W3 = first*self.r21*a2 + 2*second*(b1 - .5*b2*(self.r32 + 1.))
        if self.flat.any():
            const = self.A/(3.*D1)
            W1, W2, W3 = (np.where(self.flat, const, W) for W in (W1, W2, W3))
        out = np.zeros(self.nvert, complex)
        for k, Wk in enumerate((W1, W2, W3)):
            out.real += np.bincount(self.tri[:, k], weights=Wk.real, minlength=self.nvert)
            out.imag += np.bincount(self.tri[:, k], weights=Wk.imag, minlength=self.nvert)
        return out


class Segments:
    """The same rule in one dimension, for the boundary families, which are integrated along one coordinate."""

    def __init__(self, x, nu):
        self.h = np.diff(np.asarray(x, float))
        self.nu = np.asarray(nu, float).ravel()
        self.dn = np.diff(self.nu)
        self.nvert = len(self.nu)
        self.ntri = len(self.h)

    def weights(self, z):
        D = self.nu[:-1] - z
        l0, l1, _ = logs(self.dn/D)
        out = np.zeros(self.nvert, complex)
        out[:-1] += self.h*(l0 - l1)/D
        out[1:] += self.h*l1/D
        return out


def nu_range(resp, l):
    sel = resp.l_of_row == l
    return float(resp.nu[sel].min()), float(resp.nu[sel].max())


def resonant_harmonics(resp, im_lo, im_hi, margin):
    """The harmonics whose frequencies come within `margin` of a resonance -Im s of the rectangle. Every other
    harmonic keeps its Gauss-Legendre rows, whose poles are then at least that margin from the contour."""
    lo, hi = -im_hi - margin, -im_lo + margin
    out = []
    for l in resp.l:
        a, b = nu_range(resp, int(l))
        if b >= lo and a <= hi:
            out.append(int(l))
    return out


def bary_matrix(n, x_new):
    """Barycentric Lagrange interpolation from the n Gauss-Legendre nodes of [-1, 1] to x_new: the orbit
    quantities are known there to spectral accuracy, so moving them to the refined grid costs nothing."""
    nodes, w = np.polynomial.legendre.leggauss(n)
    lam = (-1.)**np.arange(n)*np.sqrt((1 - nodes**2)*w)
    diff = x_new[:, None] - nodes[None, :]
    hit = np.abs(diff) < 1e-14
    diff = np.where(hit, 1., diff)
    M = lam[None, :]/diff
    M = M/M.sum(axis=1, keepdims=True)
    rows = hit.any(axis=1)
    M[rows] = hit[rows].astype(float)
    return M


class _Fine:
    """What the family's F function reads off a set of orbits: the same attributes `OrbitSet` exposes."""
    __slots__ = ('disk', 'Oc', 'x', 'y', 'g')

    def __init__(self, disk, Oc, x, y):
        self.disk, self.Oc, self.x, self.y = disk, Oc, x, y
        self.g = disk.amplitude*np.exp(-.5*x*x - .5*y*y)


class Hybrid:
    """One refinement. `resp` is a `warm_response.Response` built with prune <= 0, so that its rows are the
    full (orbit, harmonic) tensor of each family and its transform can be interpolated instead of redone."""

    def __init__(self, resp, n_L, n_u, refine, pattern='alternate'):
        self.resp, self.d, self.m = resp, resp.d, resp.m
        self.n_L, self.n_u, self.refine, self.pattern = int(n_L), int(n_u), int(refine), pattern
        if not resp.keep.all():
            raise RuntimeError('the analytic rule needs an unpruned response')
        self.nl = len(resp.l)
        self.l_index = {int(l): i for i, l in enumerate(resp.l)}
        NL, NU = self.refine*self.n_L, self.refine*self.n_u
        A_full, B_full = bary_matrix(self.n_L, np.linspace(-1, 1, NL + 1)), bary_matrix(self.n_u, np.linspace(-1, 1, NU + 1))
        d = self.d
        u_max = np.sqrt(2*d.reach*d.dE)
        Lf = d.L0 + d.reach*d.dL*np.linspace(-1, 1, NL + 1)
        uf = .5*u_max*(np.linspace(-1, 1, NU + 1) + 1)
        self.fam, start = [], 0
        for k, (orbits, F_of) in enumerate(resp.families):
            if orbits.two_wells:
                raise RuntimeError('two wells: the tensor layout of the orbit rows is lost')
            nL = self.n_L if k in (0, 1) else 1
            nU = self.n_u if k in (0, 2, 3) else 1
            if orbits.n != nL*nU:
                raise RuntimeError('family %d does not hold the expected %d x %d orbits' % (k, nL, nU))
            rows = orbits.n*self.nl
            self.fam.append(dict(index=k, orbits=orbits, F_of=F_of, nL=nL, nU=nU, pieces={},
                                 G=resp.G[start:start + rows].reshape(nL, nU, self.nl, resp.basis_size),
                                 A=A_full if nL > 1 else np.ones((1, 1)), B=B_full if nU > 1 else np.ones((1, 1)),
                                 Lg=Lf if nL > 1 else np.array([orbits.L[0]]),
                                 ug=uf if nU > 1 else np.array([orbits.u[0]]), energy_measure=k != 1))
            start += rows
        if start != len(resp.nu):
            raise RuntimeError('the response rows do not account for the families')
        self.gl_mask = np.ones(len(resp.nu), bool)
        self.done = []
        self._refresh()

    def _refresh(self):
        self.G_gl = np.ascontiguousarray(self.resp.G[self.gl_mask])
        self.nu_gl = self.resp.nu[self.gl_mask]
        self.FW_gl = self.resp.FW[self.gl_mask]

    def resonant_for(self, im_lo, im_hi, margin):
        return resonant_harmonics(self.resp, im_lo, im_hi, margin)

    def pole_margin(self, im_lo, im_hi):
        """How far the nearest remaining Gauss-Legendre pole is from the rectangle's imaginary range."""
        keep = self.nu_gl
        return float(np.min(np.maximum(np.maximum(-im_hi - keep, keep + im_lo), 0.))) if len(keep) else float('inf')

    def cells(self):
        return sum(p['rule'].ntri for f in self.fam for p in f['pieces'].values())

    def add(self, harmonics):
        for l in harmonics:
            l = int(l)
            if l in self.done:
                continue
            j = self.l_index[l]
            for f in self.fam:
                A, B, Lg, ug, o = f['A'], f['B'], f['Lg'], f['ug'], f['orbits']
                shape = (f['nL'], f['nU'])
                interp = lambda v: A @ v.reshape(shape) @ B.T
                Om_r, Om_th, Oc = interp(o.Om_r), interp(o.Om_th), interp(o.Oc)
                nu = l*Om_r + self.m*Om_th
                LL, UU = np.meshgrid(Lg, ug, indexing='ij')
                d = self.d
                fine = _Fine(d, Oc.ravel(), ((LL - d.L0)/d.dL).ravel(), (.5*UU*UU/d.dE).ravel())
                c = f['F_of'](fine, nu.reshape(-1, 1), self.m)[:, 0] \
                    * ((2*np.pi)**2*(UU if f['energy_measure'] else np.ones_like(UU))/Om_r).ravel()
                Gf = np.tensordot(A, np.tensordot(B, f['G'][:, :, j, :], axes=([1], [1])), axes=([1], [1]))
                rule = Triangles(Lg, ug, nu, self.pattern) if f['nL'] > 1 and f['nU'] > 1 \
                    else Segments(Lg if f['nL'] > 1 else ug, nu)
                f['pieces'][l] = dict(G=np.ascontiguousarray(Gf.reshape(-1, self.resp.basis_size)), c=c, rule=rule)
            self.gl_mask &= self.resp.l_of_row != l
            self.done.append(l)
        self._refresh()

    def reduced(self, s):
        z = 1j*s
        out = (self.G_gl*(self.FW_gl/(self.nu_gl - z))[:, None]).T @ self.G_gl
        for f in self.fam:
            for p in f['pieces'].values():
                wt = p['c']*p['rule'].weights(z)
                G = p['G']
                out = out + (G*wt.real[:, None]).T @ G + 1j*((G*wt.imag[:, None]).T @ G)
        return out

    def T(self, s):
        return np.eye(self.resp.basis_size) + self.d.alpha*self.d.H(s)*self.reduced(s)


class Pair:
    """Two refinements of the same construction, Richardson-extrapolated: the rule is second order, so
    (4 b - a)/3 for the pair (1, 2) and the same combination for any pair with a doubled grid."""

    def __init__(self, resp, n_L, n_u, refines, pattern='alternate'):
        self.refines = (int(refines[0]), int(refines[1]))
        if self.refines[1] != 2*self.refines[0]:
            raise ValueError('the Richardson combination is written for a doubled grid')
        self.h = [Hybrid(resp, n_L, n_u, k, pattern) for k in self.refines]
        self.resp, self.d, self.m = resp, resp.d, resp.m

    def add(self, harmonics):
        for h in self.h:
            h.add(harmonics)

    def resonant_for(self, im_lo, im_hi, margin):
        return self.h[0].resonant_for(im_lo, im_hi, margin)

    def pole_margin(self, im_lo, im_hi):
        return self.h[0].pole_margin(im_lo, im_hi)

    def cells(self):
        return [h.cells() for h in self.h]

    def reduced(self, s, raw=False):
        a, b = (h.reduced(s) for h in self.h)
        rich = b + (b - a)/3.
        return (a, b, rich) if raw else rich

    def T(self, s):
        return np.eye(self.resp.basis_size) + self.d.alpha*self.d.H(s)*self.reduced(s)

    def T_many(self, S):
        return np.stack([self.T(s) for s in np.atleast_1d(np.asarray(S, complex))])

    def spread(self, s):
        """The relative difference between the pair and its coarser member: the error the rule declares."""
        a, b, rich = self.reduced(s, raw=True)
        return float(np.linalg.norm(rich - b)/np.linalg.norm(rich))


def unpruned(disk, cfg, m):
    """Stage 9's representation, with every row kept so the transform can be reused on a refined grid."""
    fam = WR.sharp_families(disk, cfg['n_L'], cfg['n_u'], cfg['n_eta'], edges=cfg['edges'])
    inner = fam[0][0]
    nodes = WR.kernel_nodes((float(inner.rp.min()), float(inner.ra.max())), cfg['node_spacing'],
                            cfg['node_margin'], cfg['node_clip'], cfg['node_shift'])
    return WR.Response(disk, fam, m, cfg['l_max'], nodes, cfg['kernel_rank_cutoff'], prune=-1.)


# ---------------------------------------------------------------- the reference integral of gate A1
def reference_integral(z, a, b, c, dL, dE, reach, n_in=600, n_out=600):
    """int dL du  exp(-x^2/2 - y^2/2) u / (a + b (L - L0) + c u^2 - z) over the same sharply bounded support as
    a population, x = (L - L0)/dL, y = u^2/(2 dE), by subtracting the pole in L, where the frequency is linear,
    and integrating what is left. It has nothing to do with the model: it is an integral whose value is known
    to fourteen digits at any growth rate, so a rule that fails at a small one is seen to fail."""
    xl, wl = np.polynomial.legendre.leggauss(n_in)
    X, wX = reach*dL*xl, reach*dL*wl
    xu, wu = np.polynomial.legendre.leggauss(n_out)
    u_max = np.sqrt(2*reach*dE)
    U, wU = .5*u_max*(xu + 1), .5*u_max*wu
    f = np.exp(-.5*(X/dL)**2)
    total = 0.
    for u, w in zip(U, wU):
        Xs = (z - a - c*u*u)/b
        fs = np.exp(-.5*(Xs/dL)**2)
        smooth = np.sum(wX*(f - fs)/(X - Xs))
        edge = fs*(np.log(reach*dL - Xs) - np.log(-reach*dL - Xs))
        total += w*np.exp(-.5*(u*u/(2*dE))**2)*u*(smooth + edge)/b
    return complex(total)


def reference_by_rule(z, a, b, c, dL, dE, reach, n_L, n_u, refine, pattern='alternate', gauss=False):
    """The same integral by the rule under test: `gauss` gives stage 9's quadrature instead."""
    u_max = np.sqrt(2*reach*dE)
    if gauss:
        xl, wl = np.polynomial.legendre.leggauss(n_L)
        xu, wu = np.polynomial.legendre.leggauss(n_u)
        X, wX = reach*dL*xl, reach*dL*wl
        U, wU = .5*u_max*(xu + 1), .5*u_max*wu
        XX, UU = np.meshgrid(X, U, indexing='ij')
        N = np.exp(-.5*(XX/dL)**2 - .5*(UU*UU/(2*dE))**2)*UU
        return complex(np.sum(np.outer(wX, wU)*N/(a + b*XX + c*UU*UU - z)))
    X = np.linspace(-reach*dL, reach*dL, refine*n_L + 1)
    U = np.linspace(0., u_max, refine*n_u + 1)
    XX, UU = np.meshgrid(X, U, indexing='ij')
    N = np.exp(-.5*(XX/dL)**2 - .5*(UU*UU/(2*dE))**2)*UU
    return complex(np.sum(Triangles(X, U, a + b*XX + c*UU*UU, pattern).weights(z)*N.ravel()))
