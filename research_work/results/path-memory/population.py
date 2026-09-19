"""RUT-1 stage 7, part S: a population defined with its support, and drawn faithfully (protocol-rut7.md).

Stage 6's annuli had four defects, none of which its own checks could see. The sampler multiplied its node
weights by 2 pi R, although in (R, L, v_r) the measure is f dR dL dv_r dtheta with no factor of R, so it
drew a different population from the one the solver built (found by the owner). "10% support" was read at
the argmax grid node, beside the field's gradient zero, where it jumps by half the target from one node to
the next, so the solve had several roots and the archived alpha moves 10-35% with the grid. E_circ(L) was
obtained by sorting and interpolating r -> L_c(r), which is not monotonic in these potentials. And the
distribution's support was never declared: five-width integration windows stood in for it, and the most
extended orbit they allow leaves the simulator's box.

Here the distribution is

    f(E, L) = A exp[-(L-L0)^2/2dL^2] exp[-(E-E_circ(L))^2/2dE^2]   for |L-L0| <= reach dL and
                                                                    0 <= E - E_circ(L) <= reach dE,
            = 0                                                     otherwise,

with E_circ(L) the GLOBAL MINIMUM over radius of Phi(r) + L^2/2r^2, A fixed by a declared total mass, and
the mass-weighted support as the population's descriptor. The historical equilibrium.py is imported
unchanged and subclassed; nothing in it is edited.

BODIES ARE DRAWN WITH THE OWNER'S CORRECTED NODE SAMPLER (research_work/annulus_sampling, v2), the sampler
of record. `draw_continuous` is an independent cross-check that places no body on a quadrature node.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import i0e

import equilibrium as EQ

GM = EQ.GM
T0 = EQ.T0


def repository_root():
    """The directory holding research_work/, found from the historical module this one subclasses."""
    return Path(EQ.__file__).resolve().parents[3]


def owner_sampler():
    """The owner's package lives outside the results tree and is imported from the repository root."""
    root = str(repository_root())
    if root not in sys.path:
        sys.path.insert(0, root)
    import research_work.annulus_sampling as AS
    return AS


class Annulus(EQ.WarmAnnulus):
    """The warm annulus with its support declared, its mass fixed and its descriptor well defined."""

    def __init__(self, *args, mass=1., reach=4., n_L=160, n_vr=96, **kw):
        super().__init__(*args, **kw)
        self.mass = None if mass is None else float(mass)
        self.reach, self.n_L, self.n_vr = float(reach), int(n_L), int(n_vr)
        self.amplitude = 1.

    # ---------------------------------------------------------------- a single-valued E_circ(L)
    def circular_energy(self, L, r, C, sigma=None):
        """The global minimum over radius of Phi(r) + L^2/2r^2, refined parabolically: single-valued for
        every L, a lower bound on the energy of every orbit of that L, and a function of L alone."""
        L = np.atleast_1d(np.asarray(L, float))
        eff = (-GM/r - C)[None, :] + .5*(L[:, None]/r[None, :])**2
        i = np.clip(eff.argmin(1), 1, len(r) - 2)
        k = np.arange(len(L))
        y0, y1, y2 = eff[k, i - 1], eff[k, i], eff[k, i + 1]
        curv = y0 - 2*y1 + y2
        shift = np.clip(np.where(curv > 0, .5*(y0 - y2)/np.where(curv > 0, curv, 1.), 0.), -1., 1.)
        return y1 - .25*(y0 - y2)*shift, r[i] + shift*(r[1] - r[0])

    # the declared support and quadrature are routed through every inherited integral
    def _grids(self, C, n_L=None, n_vr=None, reach=None, sigma=None):
        return super()._grids(C, self.n_L if n_L is None else n_L, self.n_vr if n_vr is None else n_vr,
                              self.reach if reach is None else reach, sigma)

    def surface_density(self, C, n_L=None, n_vr=None, reach=None, sigma=None):
        return super().surface_density(C, self.n_L if n_L is None else n_L,
                                       self.n_vr if n_vr is None else n_vr,
                                       self.reach if reach is None else reach, sigma)

    def _map(self, sigma):
        C = self.write_field(sigma)
        raw = self.surface_density(C, sigma=sigma)
        if self.mass is not None:
            total = float(np.trapezoid(2*np.pi*self.r*raw, self.r))
            self.amplitude = self.mass/total if total > 0 else 1.
            raw = raw*self.amplitude
        return raw, C

    def consistency_residual(self):
        """Both sides of the coupled equations on the solver's own discretization: a DISCRETE fixed point.
        Discretization accuracy is a different question, answered by refinement (gate V4)."""
        sigma, C = self._map(self.sigma)
        return dict(density_relative=float(np.max(np.abs(sigma - self.sigma))/max(np.max(np.abs(sigma)), 1e-300)),
                    field_relative=float(np.max(np.abs(C - self.C))/max(np.max(np.abs(C)), 1e-300)))

    def write_field_at(self, r):
        """C0 at arbitrary radii, including beyond the solver's grid: the same ring integral."""
        r = np.atleast_1d(np.asarray(r, float))
        R = self.r
        weight = 2*np.pi*R*self.sigma*self.alpha*self.tau_keep
        kern = np.exp(-(r[:, None] - R[None, :])**2/(2*self.w**2))*i0e(r[:, None]*R[None, :]/self.w**2)
        return np.trapezoid(kern*weight[None, :], R, axis=1)

    # ---------------------------------------------------------------- the descriptor
    def support_profile(self):
        """Inward memory acceleration over Newtonian: the profile, and its mass-weighted moments. The value
        at the density peak is reported and is NOT used: it sits beside the field's gradient zero."""
        sup = -self.write_field_gradient(self.sigma)/(GM/self.r**2)
        wgt = 2*np.pi*self.r*self.sigma
        wgt = wgt/np.trapezoid(wgt, self.r)
        mean = float(np.trapezoid(wgt*sup, self.r))
        cdf = np.cumsum(wgt*np.gradient(self.r))
        lo, hi = (float(self.r[min(np.searchsorted(cdf, q), len(self.r) - 1)]) for q in (.05, .95))
        return dict(mean=mean, std=float(np.sqrt(np.trapezoid(wgt*(sup - mean)**2, self.r))),
                    at_density_peak_node=float(sup[int(np.argmax(self.sigma))]),
                    at_5_percent_mass=float(np.interp(lo, self.r, sup)),
                    at_95_percent_mass=float(np.interp(hi, self.r, sup)),
                    maximum=float(sup.max()), minimum=float(sup.min()), field_peak=float(self.C.max()))

    def solve_mean_support(self, target, alpha_lo=1e-3, alpha_hi=.2, xtol=1e-13):
        """alpha for a declared MASS-WEIGHTED mean support, with a bracketed root finder."""
        evals = []

        def miss(alpha):
            self.alpha = float(alpha)
            self.fixed_point()
            s = self.support_profile()['mean']
            evals.append((float(alpha), float(s)))
            return s - target

        lo, hi = miss(alpha_lo), miss(alpha_hi)
        if lo*hi > 0:
            raise RuntimeError(f'mean support {target} is not bracketed by alpha in [{alpha_lo}, {alpha_hi}]')
        alpha = brentq(miss, alpha_lo, alpha_hi, xtol=xtol, rtol=1e-14)
        self.support_error = abs(miss(alpha))
        self.alpha_evaluations = evals
        return self

    def at_alpha(self, alpha):
        self.alpha = float(alpha)
        self.fixed_point()
        self.support_error = 0.
        return self

    # ---------------------------------------------------------------- the declared support
    def orbit_extent(self, reach=None, n_L=161):
        """The most extended orbit the support allows, and whether E_circ's radius runs smoothly in L."""
        reach = self.reach if reach is None else reach
        L = np.linspace(self.L0 - reach*self.dL, self.L0 + reach*self.dL, n_L)
        Ec, rc = self.circular_energy(L, self.r, self.C, self.sigma)
        phi = self.phi_total(self.r, self.C)
        peri, apo = np.inf, 0.
        for ell, ec in zip(L, Ec):
            allowed = phi + .5*(ell/self.r)**2 <= ec + reach*self.dE
            if allowed.any():
                peri, apo = min(peri, float(self.r[allowed][0])), max(apo, float(self.r[allowed][-1]))
        step = np.diff(rc)
        return dict(reach=float(reach), pericentre_min=float(peri), apocentre_max=float(apo),
                    energy_max=float(np.max(Ec + reach*self.dE)), all_bound=bool(np.max(Ec + reach*self.dE) < 0.),
                    solver_domain=[float(self.r[0]), float(self.r[-1])],
                    circular_radius_increasing=bool(np.all(step > 0)),
                    circular_radius_largest_jump=float(np.max(np.abs(step))),
                    circular_radius_range=[float(rc[0]), float(rc[-1])])

    # ---------------------------------------------------------------- quadrature moments, correct measure
    def moments(self, nodes=None):
        """From the owner's node distribution: every moment a drawn sample is compared with."""
        nodes = self.nodes() if nodes is None else nodes
        p = nodes.probabilities
        R, L, E, vr = nodes.node_radius, nodes.angular_momentum, nodes.energy, nodes.radial_speed
        vt = L/R
        m = lambda q: float(np.sum(p*q))
        out = dict(R_mean=m(R), L_mean=m(L), E_mean=m(E), vt_mean=m(vt))
        out.update(R_var=m((R - out['R_mean'])**2), L_var=m((L - out['L_mean'])**2),
                   E_var=m((E - out['E_mean'])**2), vr_var=m(vr*vr), vt_var=m((vt - out['vt_mean'])**2),
                   R_fourth=m((R - out['R_mean'])**4), vr_fourth=m(vr**4))
        return out

    # ---------------------------------------------------------------- drawing bodies
    def nodes(self):
        return owner_sampler().build_distribution(self, n_L=self.n_L, n_vr=self.n_vr, reach=self.reach)

    def verify_nodes(self, nodes=None):
        return owner_sampler().verify_distribution(self, nodes, n_L=self.n_L, n_vr=self.n_vr, reach=self.reach)

    def draw(self, n, seed=0, nodes=None):
        """The sampler of record: the owner's corrected node sampler."""
        return (self.nodes() if nodes is None else nodes).draw(n, seed)

    def draw_continuous(self, n, seed=0, batch=400000, safety=1.5):
        """The cross-check: rejection in (R, L, u) with v_r = a + u(b - a). The proposal is uniform in R, L
        and u, so the target density is f(E, L)(b - a): no factor of R, and (b - a) is the Jacobian of u."""
        rng = np.random.default_rng(seed)
        phi = CubicSpline(self.r, self.phi_total(self.r, self.C))
        Lg = np.linspace(self.L0 - self.reach*self.dL, self.L0 + self.reach*self.dL, 400)
        Ec = CubicSpline(Lg, self.circular_energy(Lg, self.r, self.C, self.sigma)[0])

        def target(R, L, u):
            E_min = phi(R) + .5*(L/R)**2
            ec = Ec(L)
            hi = 2*(ec + self.reach*self.dE - E_min)
            lo = np.maximum(2*(ec - self.reach*self.dE - E_min), 0.)
            ok = hi > lo
            a, b = np.sqrt(np.where(ok, lo, 0.)), np.sqrt(np.where(ok, hi, 1.))
            vr = a + u*(b - a)
            f = np.exp(-(L - self.L0)**2/(2*self.dL**2) - (E_min + .5*vr*vr - ec)**2/(2*self.dE**2))
            return np.where(ok, f*(b - a), 0.), vr

        scout = target(rng.uniform(self.r[0], self.r[-1], batch), rng.uniform(Lg[0], Lg[-1], batch),
                       rng.uniform(0., 1., batch))[0]
        bound = safety*float(scout.max())
        Rs, Ls, Vs, got, proposals = [], [], [], 0, 0
        while got < n:
            R = rng.uniform(self.r[0], self.r[-1], batch)
            L = rng.uniform(Lg[0], Lg[-1], batch)
            wgt, vr = target(R, L, rng.uniform(0., 1., batch))
            if wgt.max() > bound:
                raise RuntimeError('rejection bound violated: the sample would be biased')
            keep = rng.uniform(0., bound, batch) < wgt
            Rs.append(R[keep])
            Ls.append(L[keep])
            Vs.append(vr[keep])
            got += int(keep.sum())
            proposals += batch
        R, L, vr = (np.concatenate(q)[:n] for q in (Rs, Ls, Vs))
        vr = vr*rng.choice([-1., 1.], n)
        th = rng.uniform(0., 2*np.pi, n)
        c, s = np.cos(th), np.sin(th)
        vt = L/R
        self.acceptance = n/proposals
        return np.stack([R*c, R*s], axis=1), np.stack([vr*c - vt*s, vr*s + vt*c], axis=1)

    # ---------------------------------------------------------------- moving between processes
    def state(self):
        return dict(L0=self.L0, dL=self.dL, dE=self.dE, w=self.w, tau_keep=self.tau_keep, alpha=self.alpha,
                    mass=self.mass, reach=self.reach, n_L=self.n_L, n_vr=self.n_vr, amplitude=self.amplitude,
                    r=self.r.tolist(), sigma=self.sigma.tolist(), C=self.C.tolist())

    @classmethod
    def from_state(cls, s):
        r = np.asarray(s['r'], float)
        a = cls(L0=s['L0'], dL=s['dL'], dE=s['dE'], w=s['w'], tau_keep=s['tau_keep'], alpha=s['alpha'],
                mass=s['mass'], reach=s['reach'], n_L=s['n_L'], n_vr=s['n_vr'],
                r_lo=float(r[0]), r_hi=float(r[-1]), n_r=len(r))
        a.sigma, a.C, a.amplitude = np.asarray(s['sigma'], float), np.asarray(s['C'], float), s['amplitude']
        return a


def reweighted_by_radius(x, v, n, seed=0):
    """The stage 6 defect, reproduced on purpose for the negative controls: resample a correct draw with
    probability proportional to R, which is exactly what multiplying the node weights by 2 pi R did."""
    rng = np.random.default_rng(seed)
    r = np.linalg.norm(x, axis=1)
    pick = rng.choice(len(r), size=n, replace=False, p=r/r.sum())
    return x[pick], v[pick]


# ---------------------------------------------------------------- what a population is doing
def decompose(x, v, harmonics=2):
    """Separate organised motion from dispersion.

    sqrt(mean(v_r^2)) is not a dispersion: a breathing annulus and a two-lobed stream both raise it without
    any body's random motion changing, which is how stage 6 came to call a growing m = 2 mode "heating".
    Removed here by least squares are the mean radial flow, a breathing term linear in radius, and the low
    azimuthal harmonics; what is left is the residual dispersion, with the variance a p-parameter fit
    absorbs from noise subtracted.
    """
    x = np.asarray(x, float)
    v = np.asarray(v, float)
    r = np.linalg.norm(x, axis=1)
    th = np.arctan2(x[:, 1], x[:, 0])
    vr = np.sum(v*x, axis=1)/r
    vt = (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])/r
    n = len(r)
    cols = [np.ones(n), r - r.mean()]
    for k in range(1, harmonics + 1):
        cols += [np.cos(k*th), np.sin(k*th)]
    X = np.stack(cols, axis=1)
    p = X.shape[1]
    out = dict(mean_radius=float(r.mean()), radial_spread=float(r.std()),
               L_mean=float(np.mean(r*vt)), L_spread=float(np.std(r*vt)), vr_rms=float(np.sqrt(np.mean(vr*vr))))
    for name, y in (('radial', vr), ('azimuthal', vt)):
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X@beta
        out[f'{name}_mean_flow'] = float(beta[0])
        out[f'{name}_breathing'] = float(beta[1])
        out[f'{name}_m1_streaming'] = float(np.hypot(beta[2], beta[3]))
        out[f'{name}_m2_streaming'] = float(np.hypot(beta[4], beta[5])) if harmonics >= 2 else 0.
        out[f'residual_{name}_dispersion'] = float(np.sqrt(np.sum(resid**2)/max(n - p, 1)))
    # the shape of the annulus itself: displacement harmonics, and the bare density harmonics
    dr = r - r.mean()
    out['m1_displacement'] = float(2*np.abs(np.mean(dr*np.exp(-1j*th))))
    out['m2_displacement'] = float(2*np.abs(np.mean(dr*np.exp(-2j*th))))
    out['density_harmonics'] = [float(np.abs(np.mean(np.exp(-1j*m*th)))) for m in range(1, 5)]
    return out
