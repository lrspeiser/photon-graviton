"""CC-2 stage 2B-F1: a homogeneous coherent field decaying into companion pairs (see protocol.md).

The decay kinematics, the field as an added source channel of stage 2A's engine (mc.Model), and the quadratures the
validation compares with. The field is at rest in each system's frame and decays at a constant rate (its lifetime is
much longer than the span). Companions are born at a rate q per unit volume everywhere, each with the decay speed v_d
in an isotropic direction. Born companions that are confined become tracers; the rest join the far-field bath, whose
density at R_b grows as q t. Units are mc's: kpc, km/s, Msun, and kpc/(km/s) for time.
"""
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-formation'))
import mc  # noqa: E402

C_KMS = 299792.458
RUNAWAY_COST = 1000.      # stop a run whose confined mass passes this many baryon masses (the cost gate is 10, lenient 20)
OPAQUE_TAU = 3.           # stop a collisional run whose bath optical depth to r_half passes this (the regime ends at 0.3)


class OpaqueBath(RuntimeError):
    """The bath became opaque inside r_half: outside the model's regime, and the run is stopped."""


# --- decay kinematics: chi -> C + C at rest, with m_chi = 2 m_C (1 + eps)
def decay_speed(eps):
    """The daughters' common speed, c sqrt(1 - (1+eps)^-2), written without cancellation."""
    return C_KMS*math.sqrt(eps*(2 + eps))/(1 + eps)


def eps_for_speed(v_kms):
    """The mass defect giving daughters of speed v: 1 + eps = gamma, gamma - 1 written without cancellation."""
    b2 = (v_kms/C_KMS)**2
    s = math.sqrt(1 - b2)
    return b2/(s*(1 + s))


def per_field_energy(eps):
    """Per unit of field rest energy: the companions' rest-mass and kinetic fractions, and the pair's momentum."""
    return dict(rest_mass_fraction=1/(1 + eps), kinetic_fraction=eps/(1 + eps), net_momentum=0.)


def sample_decays(n, eps, rng):
    """n decays at rest, in units m_C = c = 1: each daughter has energy 1 + eps and momentum sqrt(eps (2 + eps)), back
    to back in an isotropic direction. Returns the energies and momenta of both daughters."""
    p = math.sqrt(eps*(2 + eps))
    mu = rng.uniform(-1., 1., n)
    phi = rng.uniform(0., 2*math.pi, n)
    s = np.sqrt(1 - mu*mu)
    p1 = p*np.stack([s*np.cos(phi), s*np.sin(phi), mu], axis=1)
    E = np.full(n, 1 + eps)
    return E, p1, E.copy(), -p1


# --- the field as a source channel of stage 2A's engine
class FieldModel(mc.Model):
    """Stage 2A's engine with companions born throughout the zone by the field (rate q per unit volume, speed v_d,
    isotropic) and an incident bath at speed v_d whose density grows as q t. The constructor's rho_inf should be the
    density at the middle of the first step, q Delta/2; u_k and f_k should describe the single speed v_d."""

    def __init__(self, *args, q=0., v_d=1., **kw):
        self.q, self.v_d = float(q), float(v_d)
        self.field_pools = self.field_m = self.field_frac = None
        self._clock = 0.                                   # time since the start of the run, for F4's bookkeeping
        self.pool_factor = 1.                              # scales the pools' draws (the revision's numerical control)
        super().__init__(*args, **kw)

    def make_field_pools(self, n_per=1000, n_min=100, n_cap=20000):
        """One pool of births per logarithmic shell: positions uniform in volume, speed v_d, isotropic directions, kept
        when confined in the current potential. After the first pool (n_per draws), a shell's number of draws is set
        from the previous pool's confined fraction, aiming at n_min confined samples, as in stage 2A's pools. pool_factor
        multiplies n_per and n_min."""
        n_per, n_min = int(round(n_per*self.pool_factor)), int(round(n_min*self.pool_factor))
        edges = np.exp(np.linspace(self.lr[0], self.lr[-1], self.n_shell + 1))
        if self.field_frac is None:
            sizes = np.full(self.n_shell, n_per, np.int64)
        else:
            sizes = np.clip(np.ceil(n_min/np.maximum(self.field_frac, n_min/n_cap)), n_per, n_cap).astype(np.int64)
        pools, frac = [], np.zeros(self.n_shell)
        for j in range(self.n_shell):
            a3, b3 = edges[j]**3, edges[j + 1]**3
            k = int(sizes[j])
            r = np.cbrt(a3 + self.rng.uniform(0., 1., k)*(b3 - a3))
            mu = self.rng.uniform(-1., 1., k)             # cosine of the angle to the radial direction
            psi = self.rng.uniform(0., 2*math.pi, k)
            vt = self.v_d*np.sqrt(1. - mu*mu)
            vloc = np.stack([self.v_d*mu, vt*np.cos(psi), vt*np.sin(psi)], axis=1)
            x, v = np.empty((k, 3)), np.empty((k, 3))
            mc.place(r.copy(), vloc.copy(), x, v)
            E = mc.energies(x, v, self.lr0, self.dl, self.n, self.phi, self.g, self.r_lo)
            conf = mc.confined_many(r, E, (r*vt)**2, self.lr0, self.dl, self.n, self.phi, self.r_lo)
            frac[j] = float(conf.mean())
            Eu = E[~conf]
            pools.append(dict(V=4/3*math.pi*(b3 - a3), f=frac[j], r=r[conf], vloc=vloc[conf],
                              E_unconfined=float(Eu.mean()) if len(Eu) else 0.))
        self.field_pools, self.field_frac = pools, frac

    def make_pools(self, *a, **kw):
        out = super().make_pools(*a, **kw)
        if self.q > 0:
            self.make_field_pools()
        return out

    def set_field_masses(self, T, n_field, weight_range=10.):
        """Tracer mass per shell for the field's births, proportional to the square root of the shell's confined
        production rate and within weight_range of the heaviest (stage 2A's rule). A shell with no confined birth yet
        takes the nearest producing shell's mass."""
        rates = np.array([self.q*p['V']*p['f'] for p in self.field_pools])
        root = np.sqrt(np.maximum(rates, 0.))
        if not root.max() > 0:
            self.field_m = np.ones(len(root))
            return
        c = T*root.sum()/n_field
        top = c*float(root.max())
        on = np.nonzero(root > 0)[0]
        near = on[np.abs(np.arange(len(root))[:, None] - on[None, :]).argmin(axis=1)]
        self.field_m = np.clip(c*root[near], top/weight_range, top)

    def _source_births(self, Delta, L):
        if self.q <= 0 or self.field_pools is None:
            return None
        for key in ('field_births', 'field_born_mass', 'E_field_births', 'field_exported_mass', 'E_field_exported',
                    'seedless_predicted', 'seedless_birth_variance'):
            L.setdefault(key, 0.)
        if 'ii' in self.channels:
            self._book_seedless(Delta, L)
        self._clock += Delta
        xs, vs, ms = [], [], []
        for j, p in enumerate(self.field_pools):
            made = self.q*p['V']*Delta                      # companion mass the field makes in this shell this step
            L['field_exported_mass'] += made*(1 - p['f'])
            L['E_field_exported'] += made*(1 - p['f'])*p['E_unconfined']
            if not len(p['r']):
                continue
            k = self.rng.poisson(made*p['f']/self.field_m[j])
            if not k:
                continue
            idx = self.rng.integers(0, len(p['r']), k)
            x, v = np.empty((k, 3)), np.empty((k, 3))
            mc.place(p['r'][idx].copy(), p['vloc'][idx].copy(), x, v)
            xs.append(x); vs.append(v); ms.append(np.full(k, self.field_m[j]))
        if not ms:
            return np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0), 0.
        x, v, m = np.concatenate(xs), np.concatenate(vs), np.concatenate(ms)
        E = float(m @ mc.energies(x, v, self.lr0, self.dl, self.n, self.phi, self.g, self.r_lo))
        L['field_births'] += float(len(m)); L['field_born_mass'] += float(m.sum()); L['E_field_births'] += E
        return x, v, m, E

    def _book_seedless(self, Delta, L):
        """F4's bookkeeping, which draws no random numbers. Stage 2A's engine redraws its seedless pools every ten steps,
        so F4's prediction integrates the production of the pools each step used, weighted by the exact integral of
        q² t² over the step (the engine uses the step's midpoint density), together with the variance of the births
        drawn from them: tracer masses differ between shells, and one event can confine both partners. Called after the
        step's seedless births, whose pools, density and tracer masses are still in place."""
        t = self._clock
        w = self.q*self.q*((t + Delta)**3 - t**3)/3
        for j, p in enumerate(self.pools):
            if p['p'] is None or not np.isfinite(self.shell_m[j]):
                continue
            g, nc = p['ev'][:, 1], p['ev'][:, 2]
            unit = .5*self.sm*p['W']
            L['seedless_predicted'] += unit*p['g_n']*w
            L['seedless_birth_variance'] += unit*self.rho_inf**2*Delta*self.shell_m[j]*float(np.sum(g*nc*nc))/p['draws']

    def _advance_source(self, t_next, Delta):
        if self.q > 0:
            self.rho_inf = self.q*(t_next + .5*Delta)          # the incident density at the middle of the next step
        # early stops, far beyond the declared gates and regime: a born-bound population whose own gravity binds ever
        # more births runs away, and an opaque collisional bath is outside the model and slow to integrate
        Mb = float(self.M_b[-1])
        if float(self.m.sum()) > RUNAWAY_COST*Mb:
            raise mc.RunawayError(f'non-baryonic mass above {RUNAWAY_COST:g} times the baryons at t={t_next/mc.PER_GYR:.3f} Gyr')
        if self.sm > 0 and float(np.interp(math.log(self.r_half), self.lr, self.tau)) > OPAQUE_TAU:
            raise OpaqueBath(f'optical depth to r_half above {OPAQUE_TAU:g} at t={t_next/mc.PER_GYR:.3f} Gyr')

    def run(self, T_gyr, Delta_max_gyr, snaps_gyr, n_target=16000, n_field=16000, **kw):
        """As mc.Model.run. Tracer masses come from the end-of-span incident density q T for the seedless channel
        (whose births then total about n_target/3, since the rate grows as t²) and from the current potential for the
        field's births (about n_field over the span)."""
        self._clock = 0.
        if self.q > 0:
            T = T_gyr*mc.PER_GYR
            rho0 = self.rho_inf
            self.rho_inf = self.q*T
            self.bath_tables()
            mc.Model.make_pools(self)
            if self.shell_m is None:
                self.set_shell_masses(T, n_target)
            self.rho_inf = rho0
            self.bath_tables()
            self.make_field_pools()
            self.set_field_masses(T, n_field)
        return super().run(T_gyr, Delta_max_gyr, snaps_gyr, n_target=n_target, **kw)


# --- quadratures the validation compares with
def confined_fraction(model, r, v_d, n_mu=4001):
    """The fraction of isotropic directions at speed v_d that are confined at radius r in the model's current potential,
    by quadrature over mu, the cosine of the angle to the radial direction (the test depends on E and J² only)."""
    mu = np.linspace(0., 1., n_mu)
    rr = np.full(n_mu, float(r))
    vloc = np.stack([v_d*mu, v_d*np.sqrt(1 - mu*mu), np.zeros(n_mu)], axis=1)
    x, v = np.empty((n_mu, 3)), np.empty((n_mu, 3))
    mc.place(rr.copy(), vloc.copy(), x, v)
    E = mc.energies(x, v, model.lr0, model.dl, model.n, model.phi, model.g, model.r_lo)
    conf = mc.confined_many(rr, E, (rr*vloc[:, 1])**2, model.lr0, model.dl, model.n, model.phi, model.r_lo)
    return float(np.trapezoid(conf.astype(float), mu))


def born_bound_mass(model, q, T, v_d, n_r=600, n_sub=64):
    """F2's prediction: q T times the integral of the confined fraction over the zone's volume. For hot decays the
    fraction falls from 1 to 0 within a sliver of radius (angular momentum barely moves the escape energy), and a fixed
    grid misplaces that edge by up to half an interval: with 600 points the first canonical run's quadrature was 2.6%
    high at 300 km/s. So every interval whose ends differ, or lie strictly between 0 and 1, is split n_sub ways, and
    wholly confined intervals are integrated exactly. n_sub = 0 gives the plain 600-point trapezoid."""
    r = np.geomspace(model.r_lo, model.R_b, n_r)
    P = np.array([confined_fraction(model, ri, v_d) for ri in r])
    if not n_sub:
        return q*T*float(np.trapezoid(P*4*math.pi*r**3, np.log(r)))
    lr, tot = np.log(r), 0.
    for i in range(n_r - 1):
        if P[i] == P[i + 1] and P[i] in (0., 1.):
            tot += P[i]*(r[i + 1]**3 - r[i]**3)/3
        else:
            x = np.linspace(lr[i], lr[i + 1], n_sub + 1)
            Ps = np.array([confined_fraction(model, math.exp(xi), v_d) for xi in x])
            tot += float(np.trapezoid(Ps*np.exp(3*x), x))
    return 4*math.pi*q*T*tot


def cold_density(model, q, T, r_eval, n_u=64, n_s=4000):
    """F3's prediction: the density at T of companions born at rest (radial orbits through the centre) at rate q per
    unit volume since t = 0, in the model's current potential. An orbit born at r0 passes radius r at the fall time t_f
    and at Ph - t_f, repeating every Ph = 2 t_f(0); every passage within the span adds dr/|v| of time at r."""
    phi = lambda rr: model.phi_grid_at(np.asarray(rr, float))
    u, wu = np.polynomial.legendre.leggauss(n_u)
    u, wu = .5*(u + 1), .5*wu

    def t_fall(r, r0):
        """Time from rest at r0 to radius r, with r' = r0 - (r0 - r) u^2 removing the turning-point singularity."""
        r, r0 = np.broadcast_arrays(np.asarray(r, float), np.asarray(r0, float))
        d = (r0 - r)[..., None]
        rp = r0[..., None] - d*u**2
        dphi = np.maximum(phi(r0)[..., None] - phi(rp), 1e-300)
        return np.sum(wu*2*d*u/np.sqrt(2*dphi), axis=-1)
    out = np.zeros(len(r_eval))
    s = np.linspace(0., 1., n_s + 1)[1:]
    for i, r in enumerate(np.asarray(r_eval, float)):
        r0 = r + (model.R_b - r)*s**2
        tf = t_fall(np.full_like(r0, r), r0)
        Ph = 2*t_fall(np.zeros_like(r0), r0)
        N = np.where(T >= tf, np.floor((T - tf)/Ph) + 1, 0.) + np.floor((T + tf)/Ph)
        v = np.sqrt(np.maximum(2*(phi(r0) - phi(r)), 1e-300))
        f = r0**2*N/v*2*(model.R_b - r)*s                    # dr0 = 2 (R_b - r) s ds
        out[i] = q/r**2*float(np.trapezoid(np.r_[0., f], np.r_[0., s]))
    return out


def growing_bath_seedless(model, q, T):
    """Seedless production after T with the incident density growing as q t, in the model's current (frozen, transparent)
    potential, from one set of pools: their production at unit density times q² T³/3. F4 first compared the run with
    this, from pools drawn after the run; one set's sampling error (about 4% here) is not in the births' standard error,
    so F4 now uses the pools the run itself used (FieldModel._book_seedless) and keeps this value for the record."""
    rho = model.rho_inf
    model.rho_inf = 1.
    rate = model.production_rate()
    model.rho_inf = rho
    return rate*q*q*T**3/3
