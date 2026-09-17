"""RUT-1 stage 3: bodies and their written field evolved together, from an empty field.

The memory field's GRADIENT is carried on a fixed Cartesian grid. From

    dPhi/dt = -Phi/tau - q sum_i exp[-|x - X_i|^2/(2 w^2)]

it follows that

    dg/dt = -g/tau + q sum_i (x - X_i)/w^2 exp[-|x - X_i|^2/(2 w^2)],   g = grad Phi,   a_mem = -g,

so each component has an analytic source and the acceleration is read by interpolating g rather than by
differencing a stored potential -- which would reintroduce the readout error PM-2A stage B measured.
Motion is planar and the sources lie in the plane, where the three-dimensional Gaussian factorizes
exactly, so the two-dimensional grid is not an approximation.

The decay is integrated exactly over a step and the source deposited at the midpoint,

    g(t+h) = exp(-h/tau) g(t) + tau (1 - exp(-h/tau)) S(x_mid),

whose weight tends to h as h -> 0, so halving the timestep does not double the writing rate.
"""
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.special import i0e

GM = 1.


class MemoryField:
    """grad Phi on a square grid, with local Gaussian deposition and exact decay."""

    def __init__(self, half_width, spacing, w, tau, reach=4.):
        self.w, self.tau, self.reach = float(w), float(tau), float(reach)
        n = int(np.ceil(2*half_width/spacing)) + 1
        self.x0, self.d, self.n = -half_width, 2*half_width/(n - 1), n
        self.axis = self.x0 + self.d*np.arange(n)
        self.g = np.zeros((2, n, n))                      # [component, ix, iy]
        self.pad = int(np.ceil(reach*w/self.d))

    def deposit(self, positions, rates, h):
        """One exact-decay step with the sources evaluated at `positions`."""
        self.g *= np.exp(-h/self.tau)
        weight = self.tau*(1 - np.exp(-h/self.tau))       # -> h as h -> 0
        for (px, py), q in zip(positions, rates):
            i0 = int(round((px - self.x0)/self.d))
            j0 = int(round((py - self.x0)/self.d))
            lo_i, hi_i = max(i0 - self.pad, 0), min(i0 + self.pad + 1, self.n)
            lo_j, hi_j = max(j0 - self.pad, 0), min(j0 + self.pad + 1, self.n)
            if lo_i >= hi_i or lo_j >= hi_j:
                continue
            dx = self.axis[lo_i:hi_i, None] - px
            dy = self.axis[None, lo_j:hi_j] - py
            e = np.exp(-(dx*dx + dy*dy)/(2*self.w*self.w))/(self.w*self.w)
            self.g[0, lo_i:hi_i, lo_j:hi_j] += weight*q*dx*e
            self.g[1, lo_i:hi_i, lo_j:hi_j] += weight*q*dy*e

    def sample(self, positions):
        """g at arbitrary points, by cubic-spline interpolation of the grid (order 3, not linear)."""
        p = np.atleast_2d(np.asarray(positions, float))
        c = ((p - self.x0)/self.d).T                      # [axis, point]
        return np.stack([map_coordinates(self.g[k], c, order=3, mode='nearest') for k in (0, 1)], axis=1)

    def prime_with_ring(self, R, A):
        """Initialise with the mature circular-history field of amplitude A = q tau, analytically:
        Phi(r) = -A exp[-(r-R)^2/(2w^2)] I0e(r R/w^2), so g = (dPhi/dr) r_hat."""
        X, Y = np.meshgrid(self.axis, self.axis, indexing='ij')
        r = np.maximum(np.hypot(X, Y), 1e-12)
        h = 1e-5*R

        def phi(rr):
            return -A*np.exp(-(rr - R)**2/(2*self.w**2))*i0e(rr*R/(self.w**2))
        dphi = (phi(r + h) - phi(r - h))/(2*h)
        self.g[0], self.g[1] = dphi*X/r, dphi*Y/r


def _accel(pos, field, memory=True):
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    a = -GM*pos/np.maximum(r, 1e-12)**3
    if memory:
        a = a - field.sample(pos)                          # a_mem = -grad Phi
    return a


def run(positions, velocities, rates, w, tau, t_end, h, half_width=2.5, spacing=None,
        memory=True, prime_ring_R=None, samples=400):
    """Evolve bodies and field together with velocity-Verlet, the field advanced across each step.

    No prescribed trajectory, no velocity resets, no target speed, and no writing rate adjusted to reach
    a chosen orbit. The body's own footprint is not excluded from its own force.
    """
    spacing = spacing if spacing is not None else w/5
    field = MemoryField(half_width, spacing, w, tau)
    if prime_ring_R is not None:
        field.prime_with_ring(prime_ring_R, A=float(np.sum(rates))*tau)
    x = np.array(positions, float)
    v = np.array(velocities, float)
    rates = np.asarray(rates, float)
    n_steps = int(round(t_end/h))
    every = max(1, n_steps//samples)
    rec = dict(t=[], r=[], L=[], support=[], newtonian=[])
    a = _accel(x, field, memory)
    escaped = False
    for k in range(n_steps):
        v_half = v + .5*h*a
        x_mid = x + .5*h*v_half
        x = x + h*v_half
        if np.max(np.abs(x)) > .9*half_width:
            escaped = True
            break
        if memory:
            field.deposit(x_mid, rates, h)
        a = _accel(x, field, memory)
        v = v_half + .5*h*a
        if k % every == 0 or k == n_steps - 1:
            r = np.linalg.norm(x, axis=1)
            g = field.sample(x) if memory else np.zeros_like(x)
            rec['t'].append((k + 1)*h)
            rec['r'].append(r.copy())
            rec['L'].append(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])
            rec['support'].append(np.sum(-(-g)*x/r[:, None], axis=1))   # inward memory acceleration
            rec['newtonian'].append(GM/r**2)
    out = {k: np.array(val) for k, val in rec.items()}
    out['escaped'] = escaped
    out['steps'] = k + 1
    out['field'] = field
    return out


def ring_start(n, R, q_total, phase_jitter=0., speed_jitter=0., seed=0, extra_inward=0.):
    """n bodies on a circle of radius R sharing q_total, at the circular speed for GM plus any extra
    inward acceleration already present. extra_inward = 0 launches at the Newtonian circular speed, which
    is the honest choice from an EMPTY field; a body released into an already-built track must instead be
    launched in equilibrium with it, or the control measures a body dropped into a deeper well."""
    rng = np.random.default_rng(seed)
    th = 2*np.pi*np.arange(n)/n + phase_jitter*rng.standard_normal(n)
    v_c = np.sqrt(GM/R + extra_inward*R)*(1 + speed_jitter*rng.standard_normal(n))
    pos = np.stack([R*np.cos(th), R*np.sin(th)], axis=1)
    vel = np.stack([-v_c*np.sin(th), v_c*np.cos(th)], axis=1)
    return pos, vel, np.full(n, q_total/n)


def summarize(out, R0, period):
    """Support, stability, settling and formation cost, reported separately and not conflated."""
    t, r, L = out['t'], out['r'], out['L']
    if len(t) < 4:
        return dict(terminated=True, escaped=bool(out['escaped']))
    late = t >= t[-1] - 5*period
    frac = out['support']/out['newtonian']
    early = t <= min(t[-1], 5*period)
    amp = lambda m: float(np.mean(r[m].max(axis=0) - r[m].min(axis=0)))
    return dict(
        escaped=bool(out['escaped']), orbits=float(t[-1]/period),
        support_fraction_mean_late=float(np.mean(frac[late])),
        support_fraction_final=float(np.mean(frac[-1])),
        radius_initial=float(R0), radius_final_mean=float(np.mean(r[-1])),
        radius_change_fraction=float(np.mean(r[-1])/R0 - 1),
        radius_spread_final=float(r[-1].max() - r[-1].min()),
        radius_min_over_run=float(r.min()), radius_max_over_run=float(r.max()),
        radius_change_note=('measured at TERMINATION, which for an escaped run is the moment a body left '
                            'the domain, not a settled state: read it with radius_min_over_run and '
                            'radius_max_over_run' if bool(out['escaped']) else 'measured at the end of a '
                            'completed run'),
        angular_momentum_final_over_initial=float(np.mean(L[-1])/np.mean(L[0])),
        oscillation_amplitude_early=amp(early), oscillation_amplitude_late=amp(late),
        settling_ratio=float(amp(late)/amp(early)) if amp(early) > 0 else float('nan'),
        bounded=bool(not out['escaped'] and np.all(np.isfinite(r[-1]))))
