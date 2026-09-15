"""RC-2, part 1: the open-region donor-and-companion model (see protocol.md and its Amendment 1).

Stage 2B-F1's field channel in stage 2A's engine, with no incident bath and no zone of influence. Companions are born
throughout an open region of radius R_comp (the engine's R_b) and every one is followed, bound or unbound; a companion
that crosses R_comp leaves, and nothing enters. What gravitates is the change from the initial state: the baryons,
every companion, and the donor's depletion, -q t per unit volume, because conversion removes donor mass wherever a
companion is born. The static bath and its subtraction of the incident density are gone, so the donor is counted once.
Each tracer carries its birth radius, birth time and birth angular momentum, and the time it first came within 25 kpc.
"""
import math
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-source'))
import field as FD  # noqa: E402  (puts companion-formation on the path)
import mc  # noqa: E402

NET_RUNAWAY = FD.RUNAWAY_COST       # stop once the net contrast (companions minus depletion) passes this many mass units
N_GRID_REF, N_SHELL_REF = 1024, 40  # the engine's grid nodes and birth shells for a region of the reference size
TAGS = ('birth_r_kpc', 'birth_t_Gyr', 'birth_j_kpc_kms', 'first_inside_25_Gyr')
ARRIVAL_R = 25.                     # kpc: the radius whose first crossing is recorded


def grid_for(R_comp, R_ref, r_half):
    """Grid nodes and birth shells for a region of radius R_comp at the log spacing the engine uses for a region of
    radius R_ref, so that every region resolves the same radii alike (the engine's inner edge is 10^-3 r_half)."""
    k = math.log(R_comp/(1e-3*r_half))/math.log(R_ref/(1e-3*r_half))
    return max(int(round(N_GRID_REF*k)), 64), max(int(round(N_SHELL_REF*k)), 4)


def wquant(values, weights, qs):
    """Weighted quantiles (the value below which the fraction q of the weight lies)."""
    o = np.argsort(values)
    v, c = values[o], np.cumsum(weights[o])
    c = c/c[-1]
    return {f'{q:g}': float(np.interp(q, c, v)) for q in qs}


class DonorModel(FD.FieldModel):
    """2B-F1's FieldModel without collisions or bath gravity, with the donor's depletion in the potential.

    track_unbound keeps every birth as a tracer; without it, births unbound at birth leave at once, as in 2B-F1. quiet
    selects the quiet representation (Amendment 1, B2) instead of Poisson draws from the pools. t_off (internal time
    units) ends production, and the depletion with it. mass_unit (default: the baryons inside the region) scales the
    ratios and the runaway stop, for runs without baryons. score_radii adds each snapshot's enclosed companion mass and
    depletion there. The births inside each aperture, and their Poisson variance, are counted as they are drawn, so that
    each aperture's net contrast splits exactly into (births inside - depletion inside) and the net inflow."""

    def __init__(self, *args, apertures=(), depletion=True, mass_ref_kpc=None, track_unbound=True, t_off=None,
                 mass_unit=None, score_radii=None, quiet=False, **kw):
        kw['channels'] = ()
        kw['bath_gravity'] = False
        self.apertures = tuple(float(a) for a in apertures)
        self.depletion = bool(depletion)
        self.mass_ref = None if mass_ref_kpc is None else float(mass_ref_kpc)
        self.track_unbound = bool(track_unbound)
        self.quiet = bool(quiet)
        self.t_off = None if t_off is None else float(t_off)
        self.score_radii = None if score_radii is None else np.asarray(score_radii, float)
        self.born_inside = np.zeros(len(self.apertures))
        self.born_m2_inside = np.zeros(len(self.apertures))
        self.carry = None
        super().__init__(*args, **kw)
        self.mass_unit = float(self.M_b[-1]) if mass_unit is None else float(mass_unit)

    def depleted_mass(self, r, t=None):
        """The donor mass removed inside r by time t (every birth): q min(t, t_off) (4/3) pi r^3."""
        t = self._clock if t is None else t
        if self.t_off is not None:
            t = min(t, self.t_off)
        r = np.asarray(r, float)
        return self.q*t*4/3*math.pi*r**3 if self.depletion else np.zeros_like(r)

    def total_mass(self):
        return super().total_mass() - self.depleted_mass(self.r)

    def make_pools(self, *a, **kw):
        out = super().make_pools(*a, **kw)
        if self.track_unbound:
            # every draw is kept, so the pools do not depend on the potential. NaN turns off the engine's regeneration on
            # a 1% change of the potential at r_half, which divides by that potential (zero in a run without baryons),
            # and leaves the regular regeneration every `regen` steps
            self.pool_phi = float('nan')
        return out

    def make_field_pools(self, n_per=1000, n_min=100, n_cap=20000):
        """With track_unbound, every draw is kept, bound or not (f = 1), and the fraction bound at birth is only
        recorded. Without it, 2B-F1's pools, which keep the bound draws alone."""
        if not self.track_unbound:
            return super().make_field_pools(n_per, n_min, n_cap)
        edges = np.exp(np.linspace(self.lr[0], self.lr[-1], self.n_shell + 1))
        pools, frac = [], np.zeros(self.n_shell)
        for j in range(self.n_shell):
            a3, b3 = edges[j]**3, edges[j + 1]**3
            r = np.cbrt(a3 + self.rng.uniform(0., 1., n_per)*(b3 - a3))
            mu = self.rng.uniform(-1., 1., n_per)         # cosine of the angle to the radial direction
            psi = self.rng.uniform(0., 2*math.pi, n_per)
            vt = self.v_d*np.sqrt(1. - mu*mu)
            vloc = np.stack([self.v_d*mu, vt*np.cos(psi), vt*np.sin(psi)], axis=1)
            x, v = np.empty((n_per, 3)), np.empty((n_per, 3))
            mc.place(r.copy(), vloc.copy(), x, v)
            E = mc.energies(x, v, self.lr0, self.dl, self.n, self.phi, self.g, self.r_lo)
            frac[j] = float(mc.confined_many(r, E, (r*vt)**2, self.lr0, self.dl, self.n, self.phi, self.r_lo).mean())
            pools.append(dict(V=4/3*math.pi*(b3 - a3), f=1., bound=frac[j], r=r, vloc=vloc, E_unconfined=0.))
        self.field_pools, self.field_frac = pools, frac

    def set_field_masses(self, T, n_field, weight_range=10.):
        """With mass_ref_kpc set, tracer masses follow one function of radius in every region: 2B-F1's rule (a shell's
        tracer mass proportional to the square root of its production, within weight_range of the heaviest) for a region
        of radius mass_ref_kpc in which every birth is kept and births total about n_field over the span. Beyond that
        radius the square-root rule continues, so a larger region adds births outside without thinning those inside.
        Without it, 2B-F1's rule applies."""
        if self.mass_ref is None:
            return super().set_field_masses(T, n_field, weight_range)
        dls = (math.log(self.mass_ref) - self.lr[0])/N_SHELL_REF
        e = np.exp(self.lr[0] + dls*np.arange(N_SHELL_REF + 1))
        root = np.sqrt(self.q*4/3*math.pi*(e[1:]**3 - e[:-1]**3))
        c = T*float(root.sum())/n_field
        edges = np.exp(np.linspace(self.lr[0], self.lr[-1], self.n_shell + 1))
        rc3 = (edges[1:]*edges[:-1])**1.5                  # the cube of each shell's geometric centre
        # a shell of the reference width centred on r holds (4/3) pi r^3 2 sinh(1.5 dls)
        self.field_m = np.maximum(c*np.sqrt(self.q*4/3*math.pi*rc3*2*math.sinh(1.5*dls)), c*float(root.max())/weight_range)

    def _quiet_births(self, Delta, L):
        """The quiet representation. Each shell receives its converted mass exactly, as whole pairs of tracers with the
        remainder carried forward. Radii are stratified in volume within the shell; each pair shares a tangential speed
        and has opposite radial velocities, and the pairs' direction cosines are stratified."""
        if self.carry is None:
            self.carry = np.zeros(len(self.field_pools))
        edges = np.exp(np.linspace(self.lr[0], self.lr[-1], self.n_shell + 1))
        xs, vs, ms = [], [], []
        for j, p in enumerate(self.field_pools):
            made = self.q*p['V']*Delta
            L['field_exported_mass'] += made*(1 - p['f'])      # zero when every birth is kept
            self.carry[j] += made*p['f']
            mj = float(self.field_m[j])
            h = int(self.carry[j]/(2*mj))                        # whole pairs
            if not h:
                continue
            self.carry[j] -= 2*h*mj
            a3, b3 = edges[j]**3, edges[j + 1]**3
            r = np.cbrt(a3 + (np.arange(h) + self.rng.random())/h*(b3 - a3))
            mu = ((np.arange(h) + self.rng.random())/h)[self.rng.permutation(h)]   # |cos| stratified, not tied to r
            psi = self.rng.uniform(0., 2*math.pi, h)
            vt = self.v_d*np.sqrt(1. - mu*mu)
            half = np.stack([self.v_d*mu, vt*np.cos(psi), vt*np.sin(psi)], axis=1)
            vloc = np.concatenate([half, half*np.array([-1., 1., 1.])])
            x, v = np.empty((2*h, 3)), np.empty((2*h, 3))
            mc.place(np.concatenate([r, r]), vloc.copy(), x, v)
            xs.append(x); vs.append(v); ms.append(np.full(2*h, mj))
        if not ms:
            return np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0), 0.
        x, v, m = np.concatenate(xs), np.concatenate(vs), np.concatenate(ms)
        E = float(m @ mc.energies(x, v, self.lr0, self.dl, self.n, self.phi, self.g, self.r_lo))
        L['field_births'] += float(len(m)); L['field_born_mass'] += float(m.sum()); L['E_field_births'] += E
        return x, v, m, E

    def _source_births(self, Delta, L):
        if self.t_off is not None and self._clock >= self.t_off*(1 - 1e-12):
            self._clock += Delta                           # production has stopped; the clock runs on
            return None
        if self.quiet and self.q > 0 and self.field_pools is not None:
            for key in ('field_births', 'field_born_mass', 'E_field_births', 'field_exported_mass', 'E_field_exported',
                        'seedless_predicted', 'seedless_birth_variance'):
                L.setdefault(key, 0.)
            self._clock += Delta
            out = self._quiet_births(Delta, L)
        else:
            out = super()._source_births(Delta, L)
        if out is None:
            return None
        for key in ('field_made_mass', 'field_expected_births', 'field_expected_bound', 'field_birth_m2'):
            L.setdefault(key, 0.)
        V = np.array([p['V'] for p in self.field_pools])
        f = np.array([p['f'] for p in self.field_pools])
        b = np.array([p.get('bound', p['f']) for p in self.field_pools])
        L['field_made_mass'] += self.q*Delta*float(V.sum())            # the donor mass converted this step
        L['field_expected_births'] += self.q*Delta*float(V @ f)        # the births' expectation
        L['field_expected_bound'] += self.q*Delta*float(V @ b)         # of which bound at birth
        x, v, m = out[0], out[1], out[2]
        if not len(m):
            return out
        L['field_birth_m2'] += float(m @ m)                             # their Poisson variance
        rr = np.linalg.norm(x, axis=1)
        for i, a in enumerate(self.apertures):
            s = m[rr < a]
            self.born_inside[i] += float(s.sum())
            self.born_m2_inside[i] += float(s @ s)
        tags = np.column_stack([rr, np.full(len(m), self._clock/mc.PER_GYR), np.linalg.norm(np.cross(x, v), axis=1),
                                np.full(len(m), np.nan)])
        return (*out, tags)

    def _advance_source(self, t_next, Delta):
        # no collisions and no bath gravity, so the incident density plays no part; the early stop looks at the net
        # contrast, since in a large region the companion mass is mostly cancelled by the donor's depletion
        if self.q > 0:
            self.rho_inf = self.q*(t_next + .5*Delta)
        if self.tags is not None and len(self.x):          # the first time within ARRIVAL_R
            new = np.nonzero(~np.isfinite(self.tags[:, 3]))[0]
            if len(new):
                xn = self.x[new]
                self.tags[new[np.einsum('ij,ij->i', xn, xn) < ARRIVAL_R**2], 3] = t_next/mc.PER_GYR
        net = float(self.m.sum()) - float(self.depleted_mass(self.R_b, t_next))
        if net > NET_RUNAWAY*self.mass_unit:
            raise mc.RunawayError(f'net contrast above {NET_RUNAWAY:g} mass units at t={t_next/mc.PER_GYR:.3f} Gyr')

    def run(self, *a, **kw):
        self.born_inside[:] = 0.
        self.born_m2_inside[:] = 0.
        self.carry = None
        self.tags = np.zeros((0, len(TAGS))) if not len(self.m) else None
        return super().run(*a, **kw)

    def snapshot(self, t, L):
        out = super().snapshot(t, L)
        out.update(apertures=self.aperture_table(), depletion_in_region=float(self.depleted_mass(self.R_b)),
                   escaped_mass=float(L['escaped_mass']), birth_history=self.birth_history())
        if self.score_radii is not None:
            rr = np.linalg.norm(self.x, axis=1) if len(self.x) else np.zeros(0)
            o = np.argsort(rr)
            cm = np.r_[0., np.cumsum(self.m[o])]
            out.update(companions_at_R=cm[np.searchsorted(rr[o], self.score_radii, side='right')].tolist(),
                       depletion_at_R=self.depleted_mass(self.score_radii).tolist())
        return out

    def aperture_table(self):
        """At each aperture: the companions inside (the positive inventory), the donor's depletion inside, their
        difference (the non-baryonic mass contrast), the baryons inside and the circular speed of all three, the births
        inside so far with their Poisson standard error and their departure from the depletion in those standard errors
        (control C2), and the net inflow of companions across the aperture (companions inside minus births inside)."""
        rr = np.linalg.norm(self.x, axis=1) if len(self.x) else np.zeros(0)
        rows = []
        for i, a in enumerate(self.apertures):
            c = float(self.m[rr < a].sum())
            dpl = float(self.depleted_mass(a))
            mb = float(np.interp(math.log(a), self.lr, self.M_b))
            born, se = float(self.born_inside[i]), math.sqrt(float(self.born_m2_inside[i]))
            rows.append(dict(r_kpc=a, companions=c, depletion=dpl, net=c - dpl, baryons=mb,
                             v_circ=math.sqrt(mc.G*max(mb + c - dpl, 0.)/a), born_inside=born, born_inside_se=se,
                             births_minus_depletion_z=(born - dpl)/se if se > 0 else None, inflow=c - born,
                             tracers=int((rr < a).sum())))
        return rows

    def birth_history(self, radii=(10., 25., 100., 200.)):
        """For the companions inside each radius: mass-weighted quantiles (10%, 50%, 90%, 99%) of birth radius, birth
        time and birth angular momentum, histograms of birth radius (logarithmic bins) and birth time (1 Gyr bins), and
        the mass that has been within 25 kpc."""
        if self.tags is None or not len(self.m):
            return None
        rr = np.linalg.norm(self.x, axis=1)
        r_edges = np.geomspace(1e-2, self.R_b, 29)
        t_edges = np.arange(0., math.ceil(self._clock/mc.PER_GYR - 1e-9) + 1.)
        out = {}
        for R in radii:
            s = rr < R
            m = self.m[s]
            if not len(m):
                out[f'{R:g}'] = dict(mass=0., tracers=0)
                continue
            T = self.tags[s]
            row = dict(mass=float(m.sum()), tracers=int(s.sum()),
                       arrived_within_25_fraction=float(m[np.isfinite(T[:, 3])].sum()/m.sum()))
            for col, name in enumerate(TAGS[:3]):
                row[name] = wquant(T[:, col], m, (.1, .5, .9, .99))
            row['birth_r_hist'] = np.histogram(T[:, 0], bins=r_edges, weights=m)[0].tolist()
            row['birth_t_hist'] = np.histogram(T[:, 1], bins=t_edges, weights=m)[0].tolist()
            out[f'{R:g}'] = row
        return dict(radii=out, r_edges_kpc=r_edges.tolist(), t_edges_Gyr=t_edges.tolist())


def depletion_potential_check(model, t):
    """D2: the grid potential of the depletion alone against the analytic uniform sphere, 2/3 pi G q t (R_b^2 - r^2),
    zero at the region's edge. Returns the largest error over the grid relative to the central depth."""
    M = -model.depleted_mass(model.r, t)
    g = mc.G*M/model.r**2
    phi = cumulative_trapezoid((g*model.r)[::-1], model.lr[::-1], initial=0)[::-1]
    exact = 2/3*math.pi*mc.G*model.q*min(t, model.t_off or t)*(model.R_b**2 - model.r**2)
    return float(np.max(np.abs(phi - exact)))/float(np.max(np.abs(exact)))
