"""RUT-1 stage 3: bodies and their written field evolved together, from an empty field.

Two constitutive models share this machinery, differing only in how a disturbance becomes force-producing.
With S the writing pattern, E an excitation field and C the matured, force-producing depth:

    one-stage:   dC/dt = S - C/tau_keep
    two-stage:   tau_form dE/dt = S - E,   dC/dt = E - C/tau_keep

and in both Phi_mem = -C with a_mem = +grad C. The one-stage model is exactly tau_form = 0, and the steady
state is the same in both -- constant forcing gives C = tau_keep S -- so the comparison changes the
maturation process and not the strength of the equilibrium field.

Each field is carried on a fixed Cartesian grid as its VALUE and its two GRADIENT components, because
differencing a stored potential would reintroduce the readout error PM-2A stage B measured, while the
value is needed to state a body's specific orbital energy -- which is not conserved while the field evolves,
so it supports a consistency check rather than an escape test. Every source term
is analytic. Motion is planar and the sources lie in the plane, where the 3-D Gaussian factorizes exactly.

The timestep is adaptive, set by the shortest local dynamical time, because a fixed step badly mishandles
close passages: at pericentre 0.063 with h = 0.01 a bound Kepler orbit acquires a 4.4% energy error while
its angular momentum stays good to 1e-14. The exponential deposit weight tau(1 - exp(-h/tau)) is a
function of the step actually taken, so variable stepping needs no other change and halving the step still
does not double the writing rate.
"""
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.special import i0e

GM = 1.
COMPLETED, LEFT_DOMAIN, UNRESOLVED_CENTRE = 'completed', 'left_domain', 'entered_unresolved_centre'


class MemoryField:
    """One- or two-stage memory field: value and gradient of C, plus of E when tau_form > 0."""

    def __init__(self, half_width, spacing, w, tau_keep, tau_form=0., reach=5.):
        self.w, self.tau_keep, self.tau_form = float(w), float(tau_keep), float(tau_form)
        n = int(np.ceil(2*half_width/spacing)) + 1
        self.x0, self.d, self.n = -half_width, 2*half_width/(n - 1), n
        self.axis = self.x0 + self.d*np.arange(n)
        self.C = np.zeros((3, n, n))                      # value, d/dx, d/dy
        self.E = np.zeros((3, n, n)) if tau_form > 0 else None
        self.pad = int(np.ceil(reach*w/self.d))

    def _source(self, positions, rates):
        """S and its gradient, deposited only where the Gaussian is not negligible."""
        out = np.zeros((3, self.n, self.n))
        for (px, py), q in zip(positions, rates):
            i0 = int(round((px - self.x0)/self.d))
            j0 = int(round((py - self.x0)/self.d))
            a, b = max(i0 - self.pad, 0), min(i0 + self.pad + 1, self.n)
            c, e = max(j0 - self.pad, 0), min(j0 + self.pad + 1, self.n)
            if a >= b or c >= e:
                continue
            dx = self.axis[a:b, None] - px
            dy = self.axis[None, c:e] - py
            g = q*np.exp(-(dx*dx + dy*dy)/(2*self.w*self.w))
            out[0, a:b, c:e] += g
            out[1, a:b, c:e] += -g*dx/(self.w*self.w)     # grad of exp(-|x-X|^2/2w^2)
            out[2, a:b, c:e] += -g*dy/(self.w*self.w)
        return out

    def advance(self, positions, rates, h):
        """One step, integrated exactly for a source held constant across it."""
        S = self._source(positions, rates)
        rk = np.exp(-h/self.tau_keep)
        if self.E is None:
            self.C = rk*self.C + self.tau_keep*(1 - rk)*S
            return
        rf = np.exp(-h/self.tau_form)
        E0 = self.E
        self.E = rf*E0 + (1 - rf)*S
        # dC/dt = E(t) - C/tau_keep with E(t) = S + (E0 - S) exp(-t/tau_form)
        tk, tf = self.tau_keep, self.tau_form
        # C(h) = rk C0 + tk(1 - rk) S + (E0 - S) tk tf/(tf - tk) (rf - rk), by integrating
        # dC/dt = S + (E0 - S) exp(-t/tau_form) - C/tau_keep exactly across the step
        mix = (tk*tf/(tf - tk))*(rf - rk) if abs(tk - tf) > 1e-12*tk else h*rk
        self.C = rk*self.C + tk*(1 - rk)*S + mix*(E0 - S)

    def sample_gradient(self, positions):
        p = np.atleast_2d(np.asarray(positions, float))
        c = ((p - self.x0)/self.d).T
        return np.stack([map_coordinates(self.C[k], c, order=3, mode='nearest') for k in (1, 2)], axis=1)

    def sample_value(self, positions):
        p = np.atleast_2d(np.asarray(positions, float))
        return map_coordinates(self.C[0], ((p - self.x0)/self.d).T, order=3, mode='nearest')

    def prime_with_ring(self, R, A):
        """The mature circular-history field: Phi = -A exp[-(r-R)^2/2w^2] I0e(rR/w^2), so C = -Phi."""
        X, Y = np.meshgrid(self.axis, self.axis, indexing='ij')
        r = np.maximum(np.hypot(X, Y), 1e-12)
        hh = 1e-5*R
        phi = lambda rr: -A*np.exp(-(rr - R)**2/(2*self.w**2))*i0e(rr*R/(self.w**2))
        dphi = (phi(r + hh) - phi(r - hh))/(2*hh)
        self.C[0] = -phi(r)
        self.C[1], self.C[2] = -dphi*X/r, -dphi*Y/r
        if self.E is not None:
            self.E[:] = self.C/self.tau_keep              # the excitation in equilibrium with that C


def _accel(pos, field, memory=True):
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    a = -GM*pos/np.maximum(r, 1e-12)**3
    return a + field.sample_gradient(pos) if memory else a


def run(positions, velocities, rates, w, tau_keep, t_end, h_max, tau_form=0., half_width=2.5,
        spacing=None, memory=True, prime_ring_R=None, write=True, freeze_field=False,
        r_min=.25, eta=.005, samples=400):
    """Evolve bodies and field together. No prescribed trajectory, no velocity resets, no target speed.

    `write=False` lets an existing field decay without new deposition and `freeze_field=True` holds it
    fixed; with both defaults on, the three together are the split primed controls.
    """
    spacing = spacing if spacing is not None else w/5
    field = MemoryField(half_width, spacing, w, tau_keep, tau_form)
    if prime_ring_R is not None:
        field.prime_with_ring(prime_ring_R, A=float(np.sum(rates))*tau_keep)
    x = np.array(positions, float)
    v = np.array(velocities, float)
    rates = np.asarray(rates, float)
    zero = np.zeros_like(rates)
    rec = dict(t=[], r=[], L=[], support=[], newtonian=[], vr=[], specific_orbital_energy=[])
    a = _accel(x, field, memory)
    t, status, next_sample = 0., COMPLETED, 0.
    while t < t_end:
        r = np.linalg.norm(x, axis=1)
        if r.min() < r_min:
            status = UNRESOLVED_CENTRE
            break
        h = min(h_max, eta*float(np.min(np.sqrt(r**3/GM))), t_end - t)
        v_half = v + .5*h*a
        x_mid = x + .5*h*v_half
        x = x + h*v_half
        if np.max(np.abs(x)) > .9*half_width:
            status = LEFT_DOMAIN
            break
        if memory and not freeze_field:
            field.advance(x_mid, rates if write else zero, h)
        a = _accel(x, field, memory)
        v = v_half + .5*h*a
        t += h
        if t >= next_sample or t >= t_end:
            next_sample = t + t_end/samples
            r = np.linalg.norm(x, axis=1)
            g = field.sample_gradient(x) if memory else np.zeros_like(x)
            phi_mem = -field.sample_value(x) if memory else np.zeros(len(x))
            rec['t'].append(t)
            rec['r'].append(r.copy())
            rec['L'].append(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])
            rec['vr'].append(np.sum(v*x, axis=1)/r)
            rec['support'].append(-np.sum(g*x/r[:, None], axis=1))   # inward memory acceleration
            rec['newtonian'].append(GM/r**2)
            rec['specific_orbital_energy'].append(.5*np.sum(v*v, axis=1) - GM/r + phi_mem)
    out = {k: np.array(val) for k, val in rec.items()}
    out.update(status=status, field=field, t_final=t)
    return out


def ring_start(n, R, q_total, phase_jitter=0., speed_jitter=0., seed=0, extra_inward=0.):
    """n bodies on a circle of radius R sharing q_total, at the circular speed for GM plus any extra
    inward acceleration already present. extra_inward = 0 is the honest launch from an EMPTY field; a body
    released into a built track must be launched in equilibrium with it, or the control measures a body
    dropped into a deeper well. Phase and speed jitter are separate so their effects can be told apart."""
    rng = np.random.default_rng(seed)
    th = 2*np.pi*np.arange(n)/n + phase_jitter*rng.standard_normal(n)
    v_c = np.sqrt(GM/R + extra_inward*R)*(1 + speed_jitter*rng.standard_normal(n))
    pos = np.stack([R*np.cos(th), R*np.sin(th)], axis=1)
    vel = np.stack([-v_c*np.sin(th), v_c*np.cos(th)], axis=1)
    return pos, vel, np.full(n, q_total/n)


def summarize(out, R0, period):
    """Support, stability, settling, heating and formation cost, reported separately."""
    t, r, L = out['t'], out['r'], out['L']
    if len(t) < 4:
        return dict(status=out['status'], completed=False, orbits=float(out['t_final']/period),
                    terminated_early=True)
    late = t >= t[-1] - 5*period
    early = t <= min(t[-1], 5*period)
    amp = lambda m: float(np.mean(r[m].max(axis=0) - r[m].min(axis=0)))
    dL = L[-1]/L[0]
    return dict(
        status=out['status'], completed=bool(out['status'] == COMPLETED),
        orbits=float(t[-1]/period),
        support_fraction_mean_late=float(np.mean(out['support'][late]/out['newtonian'][late])),
        radius_initial=float(R0), radius_final_mean=float(np.mean(r[-1])),
        radius_change_fraction=float(np.mean(r[-1])/R0 - 1),
        radius_min_over_run=float(r.min()), radius_max_over_run=float(r.max()),
        radius_spread_final=float(r[-1].max() - r[-1].min()),
        angular_momentum_final_over_initial=float(np.mean(L[-1])/np.mean(L[0])),
        per_body_L_ratio_min=float(dL.min()), per_body_L_ratio_max=float(dL.max()),
        per_body_L_ratio_spread=float(dL.max() - dL.min()),
        radial_velocity_rms_late=float(np.sqrt(np.mean(out['vr'][late]**2))),
        radial_velocity_rms_early=float(np.sqrt(np.mean(out['vr'][early]**2))),
        oscillation_amplitude_early=amp(early), oscillation_amplitude_late=amp(late),
        settling_ratio=float(amp(late)/amp(early)) if amp(early) > 0 else float('nan'),
        max_specific_orbital_energy_final=float(out['specific_orbital_energy'][-1].max()),
        any_body_positive_specific_orbital_energy=bool(np.any(out['specific_orbital_energy'][-1] > 0)),
        note='the specific orbital energy |v|^2/2 - GM/r - C is NOT conserved in an evolving field: along '
             'the declared equations its rate of change is -dC/dt at the body, so a positive value at one '
             'instant is not a permanent-escape test and a negative one is not a boundness guarantee. '
             'Leaving the domain and entering the unresolved centre are different statements again')
