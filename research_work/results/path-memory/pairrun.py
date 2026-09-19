"""RUT-1 stage 7, part X: the same bodies evolved with the field live and with it frozen (protocol-rut7.md).

Stage 6's H6 mixed two questions -- was the intended equilibrium represented correctly, and does the live
matter-field system stay put -- compared against a run with no field, which is not the constructed
population at all, and reported only azimuthal averages, which are second order in any mode amplitude and
so could not see the m = 2 mode growing in both of its runs.

Here the frozen run holds the constructed field fixed: it is the stationary reference, and any change in it
is sampling or integration, not physics. The live run lets the same bodies write. Every record separates
organised motion from dispersion, carries the support the bodies feel and the support on a fixed ring, and
saves the COMPLEX azimuthal coefficients of the source S, the excitation E, the force-producing field C
and the inward pull P on that ring -- phases, not only powers -- so a growing pattern is visible as what it
is, with its rate and its pattern speed.

A run can be stopped with its full state -- bodies, E and C -- and resumed bit for bit; restarting from
the bodies alone would discard the very history under study.
"""
import time

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import map_coordinates

import formation as FM
import population as PO

GM = FM.GM
T0 = PO.T0
HALF_WIDTH = 3.5            # the simulator box; the declared support must fit inside 0.9 of it (gate V3)
MODES = 9                   # complex azimuthal coefficients m = 0 .. 8
RING_POINTS = 128


def prime(field, annulus):
    """The constructed axisymmetric field on the Cartesian grid: value and ANALYTIC radial gradient,
    continued smoothly beyond the solver's radial domain rather than cut to zero there."""
    X, Y = np.meshgrid(field.axis, field.axis, indexing='ij')
    rr = np.maximum(np.hypot(X, Y), 1e-12)
    r_in = np.linspace(min(float(rr.min()), .5*annulus.r[0]), annulus.r[0], 48)[:-1]
    r_out = np.linspace(annulus.r[-1], float(rr.max()) + 1e-9, 64)[1:]
    r_all = np.concatenate([r_in, annulus.r, r_out])
    C_all = np.concatenate([annulus.write_field_at(r_in), annulus.C, annulus.write_field_at(r_out)])
    dC_all = np.concatenate([annulus.write_field_gradient(annulus.sigma, r_in),
                             annulus.write_field_gradient(annulus.sigma),
                             annulus.write_field_gradient(annulus.sigma, r_out)])
    field.C[0] = CubicSpline(r_all, C_all)(rr)             # linear interpolation here costs 4e-4 of the force
    der = CubicSpline(r_all, dC_all)(rr)
    field.C[1], field.C[2] = der*X/rr, der*Y/rr
    if field.E is not None:
        field.E[:] = field.C/field.tau_keep
    return field


def ring_coefficients(field, positions, rates, radius, w):
    """c_m = (1/M) sum_k f(phi_k) exp(-i m phi_k) for S, E, C and the inward pull P in units of GM/R^2.
    A pattern f = cos[m(phi - Omega_p t)] gives c_m proportional to exp(-i m Omega_p t)."""
    ph = 2*np.pi*np.arange(RING_POINTS)/RING_POINTS
    unit = np.stack([np.cos(ph), np.sin(ph)], axis=1)
    pts = radius*unit
    c = ((pts - field.x0)/field.d).T
    sample = lambda arr: map_coordinates(arr, c, order=3, mode='nearest')
    d2 = ((pts[:, None, :] - positions[None, :, :])**2).sum(-1)
    series = dict(S=(np.exp(-d2/(2*w*w))*rates[None, :]).sum(1), C=sample(field.C[0]),
                  P=-(sample(field.C[1])*unit[:, 0] + sample(field.C[2])*unit[:, 1])/(GM/radius**2))
    if field.E is not None:
        series['E'] = sample(field.E[0])
    return {k: [[float(z.real), float(z.imag)] for z in np.fft.fft(y)[:MODES]/RING_POINTS]
            for k, y in series.items()}


class Run:
    """Velocity-Verlet in the constructed field, written to (live) or held fixed (frozen)."""

    def __init__(self, annulus, x, v, rates, live=True, tau_form=3*T0, h_max=.01, eta=.01,
                 record_every=.25*T0, ring_radius=None, light=False, primed=True):
        self.w, self.live, self.h_max, self.eta = annulus.w, bool(live), h_max, eta
        self.record_every, self.light = record_every, bool(light)
        self.field = FM.MemoryField(HALF_WIDTH, annulus.w/5, annulus.w, annulus.tau_keep, tau_form)
        if primed:
            prime(self.field, annulus)
        self.x, self.v = np.array(x, float), np.array(v, float)
        self.rates = np.asarray(rates, float)
        self.ring_radius = float(ring_radius) if ring_radius is not None else float(np.mean(np.linalg.norm(self.x, axis=1)))
        self.t, self.k, self.status, self.rows, self.steps = 0., 0, 'completed', [], 0

    def _accel(self, q):
        return (-GM*q/np.maximum(np.linalg.norm(q, axis=1, keepdims=True), 1e-12)**3
                + self.field.sample_gradient(q))

    def _record(self):
        x, v = self.x, self.v
        r = np.linalg.norm(x, axis=1)
        if self.light:
            vr = np.sum(v*x, axis=1)/r
            vt = (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])/r
            row = dict(mean_radius=float(r.mean()), radial_spread=float(r.std()),
                       vr_rms=float(np.sqrt(np.mean(vr*vr))), vt_mean=float(vt.mean()),
                       vr_mean=float(vr.mean()))
        else:
            row = PO.decompose(x, v)
            g = self.field.sample_gradient(x)
            row['support_on_bodies'] = float(np.mean(-np.sum(g*x, axis=1)/r/(GM/r**2)))
            row['coefficients'] = ring_coefficients(self.field, x, self.rates, self.ring_radius, self.w)
            row['support_on_ring'] = row['coefficients']['P'][0][0]
        row['t_periods'] = float(self.t/T0)
        self.rows.append(row)

    def advance(self, t_end, wall_limit=None):
        """Integrate to the absolute time t_end, landing exactly on every record time on the way."""
        # record times are COUNTED, k * record_every, never accumulated: a resumed run must land on exactly
        # the same instants as an uninterrupted one, and a horizon must not be missed by a rounding error.
        if not self.rows:
            self._record()
        target = int(round(t_end/self.record_every))
        started = time.time()
        a = self._accel(self.x)
        while self.k < target and self.status == 'completed':
            next_rec = (self.k + 1)*self.record_every
            r = np.linalg.norm(self.x, axis=1)
            if r.min() < .25:
                self.status = 'entered_unresolved_centre'
                break
            h = min(self.h_max, self.eta*float(np.min(np.sqrt(r**3/GM))), next_rec - self.t)
            v_half = self.v + .5*h*a
            x_mid = self.x + .5*h*v_half
            self.x = self.x + h*v_half
            if np.max(np.abs(self.x)) > .9*HALF_WIDTH:
                self.status = 'left_domain'
                break
            if self.live:
                self.field.advance(x_mid, self.rates, h)
            a = self._accel(self.x)
            self.v = v_half + .5*h*a
            if not (np.isfinite(self.x).all() and np.isfinite(self.v).all()):
                self.status = 'non_finite'               # fail fast: never grind on through NaNs
                break
            if wall_limit is not None and self.steps % 256 == 0 and time.time() - started > wall_limit:
                self.status = 'wall_budget_exceeded'
                break
            self.t += h
            self.steps += 1
            if self.t >= next_rec - 1e-12:
                self.t = next_rec
                self.k += 1
                self._record()
        return self

    # ---------------------------------------------------------------- the full state
    def snapshot(self):
        return dict(x=self.x.copy(), v=self.v.copy(), C=self.field.C.copy(),
                    E=None if self.field.E is None else self.field.E.copy(),
                    t=self.t, k=self.k, status=self.status, steps=self.steps,
                    rows=[dict(r) for r in self.rows])

    def restore(self, snap, fields=True):
        """`fields=False` restarts from the bodies alone, leaving whatever field this run was primed with:
        the mistake a full restart state exists to prevent, kept here as V7's negative control."""
        self.x, self.v = snap['x'].copy(), snap['v'].copy()
        if fields:
            self.field.C = snap['C'].copy()
            if snap['E'] is not None:
                self.field.E = snap['E'].copy()
        self.t, self.k, self.status = snap['t'], snap['k'], snap['status']
        self.steps, self.rows = snap['steps'], [dict(r) for r in snap['rows']]
        return self

    def result(self):
        return dict(status=self.status, t_final_periods=float(self.t/T0), live=self.live, steps=int(self.steps),
                    ring_radius=self.ring_radius, rows=self.rows)


def run(annulus, x, v, rates, horizon_periods, **kw):
    return Run(annulus, x, v, rates, **kw).advance(horizon_periods*T0).result()
