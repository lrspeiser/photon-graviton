"""CC-1 bound systems in the co-scaling completion (see protocol.md). Covers a classical Coulomb atom in the
field frame (N1) and in the conformal frame (N2), and a one-dimensional quantum clock in the field frame
(N3). Units c = hbar = m = k = 1.

Matter couples to g_m = -dt^2 + n(t)^2 dx^2 with fixed masses and charges. Maxwell's equations are
conformally invariant, so the field of a charge at rest is the flat Coulomb field in (eta, x), where
eta = int dt/n. Measured locally, that field is 1/(proper distance)^2. The atom's Lagrangian is then
    field frame (t, x):        L  = (1/2) n^2 |dx/dt|^2 + 1/(n |x|)
    conformal frame (eta, x):  L~ = (n/2) |dx/deta|^2 + 1/|x|     (static metric, mass proportional to n)
These are one action written with two time coordinates (L dt = L~ deta). Each check measures its own
frame's quantities; nothing is converted between frames by assumption.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import circulant, eigh

RTOL, ATOL = 1e-12, 1e-14


class Profile:
    """n = 1 before t = 0, then n = 1.5 - 0.5 cos(pi t/T) during the ramp, then n = 2 after t = T.
    n and ndot are continuous. eta(t) = int_0^t dt'/n and its inverse are in closed form."""

    def __init__(self, T):
        self.T, self.w, self.eta_T = float(T), math.pi/T, T/math.sqrt(2)

    def n(self, t):
        if t <= 0:
            return 1.
        if t >= self.T:
            return 2.
        return 1.5 - .5*math.cos(self.w*t)

    def dn(self, t):
        return .5*self.w*math.sin(self.w*t) if 0 < t < self.T else 0.

    def ddn(self, t):
        return .5*self.w**2*math.cos(self.w*t) if 0 < t < self.T else 0.

    def eta(self, t):
        if t <= 0:
            return t
        if t >= self.T:
            return self.eta_T + (t - self.T)/2
        h = self.w*t/2   # int dt/(1 + sin^2(wt/2)) = (sqrt2/w) arctan(sqrt2 tan(wt/2))
        return math.sqrt(2)/self.w*math.atan2(math.sqrt(2)*math.sin(h), math.cos(h))

    def t_of_eta(self, e):
        if e <= 0:
            return e
        if e >= self.eta_T:
            return self.T + 2*(e - self.eta_T)
        u = self.w*e/math.sqrt(2)
        return 2/self.w*math.atan2(math.sin(u), math.sqrt(2)*math.cos(u))

    def n_of_eta(self, e):
        return self.n(self.t_of_eta(e))


def integrate(rhs, y0, nodes, events):
    """Integrate across the profile's kinks one segment at a time. Returns event times and states per event."""
    y, times, states = np.asarray(y0, float), [[] for _ in events], [[] for _ in events]
    for a, b in zip(nodes[:-1], nodes[1:]):
        sol = solve_ivp(rhs, (a, b), y, method='DOP853', rtol=RTOL, atol=ATOL, events=events)
        if sol.status != 0 or not np.all(np.isfinite(sol.y[:, -1])):
            raise FloatingPointError(f'integration failed on [{a}, {b}]: {sol.message}')
        for i in range(len(events)):
            keep = sol.t_events[i] > a + 1e-9   # an event exactly at a segment start is reported twice
            times[i].extend(sol.t_events[i][keep])
            states[i].extend(sol.y_events[i][keep])
        y = sol.y[:, -1]
    return [np.array(t) for t in times], [np.array(s) for s in states]


def _orbit_summary(cross_t, pre_end, post_start):
    """Periods between successive ascending crossings, and the orbits wholly inside each pad."""
    P = np.diff(cross_t)
    lo, hi = cross_t[:-1], cross_t[1:]
    return P, hi <= pre_end, lo >= post_start


def classical_field_frame(prof, pad_orbits=20, e=.3):
    """N1. Canonical momentum p = n^2 dx/dt. The orbit starts at pericentre with semi-major axis 1 (period 2 pi)."""
    def rhs(t, y):
        n = prof.n(t)
        r3 = (y[0]*y[0] + y[1]*y[1])**1.5
        return [y[2]/n**2, y[3]/n**2, -y[0]/(n*r3), -y[1]/(n*r3)]

    def cross(t, y):
        return y[1]
    cross.direction = 1

    def radial(t, y):   # proportional to d(proper r^2)/dt: ndot |x|^2 + x.p/n
        return prof.dn(t)*(y[0]**2 + y[1]**2) + (y[0]*y[2] + y[1]*y[3])/prof.n(t)

    pad = 2*math.pi*pad_orbits
    y0 = [1 - e, 0., 0., math.sqrt((1 + e)/(1 - e))]
    (tc, tr), (_, yr) = integrate(rhs, y0, [-pad, 0., prof.T, prof.T + pad], [cross, radial])
    P, pre, post = _orbit_summary(np.r_[-pad, tc], 0., prof.T)
    omega = 2*math.pi/P
    # proper radius at the radial extrema; pair successive minimum and maximum into semi-major axes
    n_r = np.array([prof.n(t) for t in tr])
    r = n_r*np.hypot(yr[:, 0], yr[:, 1])
    a = .5*(r[:-1] + r[1:])[::2]
    a_coord = (.5*(r[:-1]/n_r[:-1] + r[1:]/n_r[1:]))[::2]
    w_pre, w_post = omega[pre].mean(), omega[post].mean()
    return dict(
        orbits=int(len(P)), T=prof.T, eccentricity=e,
        omega_pre=float(w_pre), omega_post=float(w_post), clock_factor=float(w_post/w_pre),
        max_frequency_deviation=float(np.max(np.abs(omega/w_pre - 1))),
        max_proper_size_deviation=float(np.max(np.abs(a/a[0] - 1))),
        coordinate_size_ratio_post_to_pre=float(a_coord[-1]/a_coord[0]),
        crossing_times=np.r_[-pad, tc],
        passed=bool(np.max(np.abs(omega/w_pre - 1)) < 1e-6 and np.max(np.abs(a/a[0] - 1)) < 1e-6))


def classical_conformal_frame(prof, pad_orbits=20, e=.3):
    """N2. Time eta, canonical momentum P = n dx/deta: a static Coulomb law with mass n(eta)."""
    def rhs(s, y):
        n = prof.n_of_eta(s)
        r3 = (y[0]*y[0] + y[1]*y[1])**1.5
        return [y[2]/n, y[3]/n, -y[0]/r3, -y[1]/r3]

    def cross(s, y):
        return y[1]
    cross.direction = 1

    pad_pre, pad_post = 2*math.pi*pad_orbits, math.pi*pad_orbits   # the eta-period halves at n = 2
    y0 = [1 - e, 0., 0., math.sqrt((1 + e)/(1 - e))]
    (se,), _ = integrate(rhs, y0, [-pad_pre, 0., prof.eta_T, prof.eta_T + pad_post], [cross])
    s = np.r_[-pad_pre, se]
    P_eta, pre, post = _orbit_summary(s, 0., prof.eta_T)
    t = np.array([prof.t_of_eta(v) for v in s])
    n_bar = np.diff(t)/P_eta                        # orbit average of n over eta
    w_eta = 2*math.pi/P_eta
    ratio = w_eta/n_bar                             # orbit-averaged eta-frequency divided by n
    return dict(
        orbits=int(len(P_eta)),
        omega_eta_pre=float(w_eta[pre].mean()), omega_eta_post=float(w_eta[post].mean()),
        clock_factor=float(w_eta[post].mean()/w_eta[pre].mean()),
        max_deviation_of_omega_eta_over_n=float(np.max(np.abs(ratio/ratio[pre].mean() - 1))),
        crossing_times_as_t=t,
        passed=bool(np.max(np.abs(ratio/ratio[pre].mean() - 1)) < 1e-6))


def soft_coulomb(r):
    return -1/np.sqrt(r*r + 1)


def _upsample(f, up):
    """Band-limited interpolation of a periodic grid function onto a grid up times finer."""
    N = len(f)
    c = np.fft.fft(f)
    C = np.zeros(N*up, complex)
    C[:N//2], C[-(N//2):] = c[:N//2], c[N//2:]
    return np.fft.ifft(C)*up


def _window_fit(t, theta, freqs):
    """Least-squares slope of a phase with sinusoidal ripples at known frequencies removed."""
    s = t - t.mean()
    cols = [np.ones_like(s), s] + [f(s) for w in freqs for f in (lambda u, w=w: np.cos(w*u), lambda u, w=w: np.sin(w*u))]
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), theta, rcond=None)
    return coef[1]


def quantum_clock(prof, N=2048, L=160., dt=.05, pad=2000., ramp_windows=20, sample_every=40, log=print):
    """N3. i phi_t = -phi_xx/(2 n^2) + U(n x) phi on a comoving grid, Strang split. Starts in
    (|0> + |1> + |2>)/sqrt3 of the n = 1 soft-Coulomb atom. The phases of the projections onto the
    instantaneous eigenstates sqrt(n) chi_j(n x) give the beat frequencies w10 and w20 in t."""
    dx = L/N
    x = (np.arange(N) - N//2)*dx
    k = 2*np.pi*np.fft.fftfreq(N, dx)
    H = circulant(np.real(np.fft.ifft(k*k/2))) + np.diag(soft_coulomb(x))
    E, V = eigh(.5*(H + H.T), subset_by_index=[0, 2])
    chi = V/math.sqrt(dx)
    chi *= np.sign(chi[np.argmax(np.abs(chi), axis=0), range(3)])
    x2 = (chi**2*(x*x)[:, None]).sum(0)*dx
    up = 16
    xf = x[0] + np.arange(N*up)*dx/up
    splines = [CubicSpline(xf, _upsample(chi[:, j], up).real, extrapolate=False) for j in range(3)]

    def states(n):
        return np.nan_to_num(np.array([math.sqrt(n)*s(n*x) for s in splines]))

    phi = chi.sum(1).astype(complex)/math.sqrt(3)
    t, steps = -pad, int(round((prof.T + 2*pad)/dt))
    h_prev = np.exp(-.5j*dt*soft_coulomb(prof.n(t)*x))
    rec_t, rec_n, rec_c = [], [], []
    for s in range(steps + 1):
        if s % sample_every == 0:
            n = prof.n(t)
            rec_t.append(t)
            rec_n.append(n)
            rec_c.append(states(n) @ phi*dx)
        if s == steps:
            break
        n_mid, n_new = prof.n(t + dt/2), prof.n(t + dt)
        h_new = np.exp(-.5j*dt*soft_coulomb(n_new*x))
        phi = h_new*np.fft.ifft(np.exp(-.5j*dt*k*k/n_mid**2)*np.fft.fft(h_prev*phi))
        h_prev, t = h_new, t + dt
        if s % 100000 == 0:
            norm = float(np.sum(np.abs(phi)**2)*dx)
            if not math.isfinite(norm) or abs(norm - 1) > 1e-9:
                raise FloatingPointError(f'quantum clock lost normalization at t={t:.1f}: {norm}')
            log(f'  N3 step {s}/{steps}, t={t:.0f}, n={prof.n(t):.4f}')
    rt, c = np.array(rec_t), np.array(rec_c)
    th10 = np.unwrap(np.angle(c[:, 1]*np.conj(c[:, 0])))
    th20 = np.unwrap(np.angle(c[:, 2]*np.conj(c[:, 0])))
    edges = np.r_[-pad, np.linspace(0, prof.T, ramp_windows + 1), prof.T + pad]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (rt >= lo) & (rt <= hi)
        w10, w20 = -_window_fit(rt[m], th10[m], ()), -_window_fit(rt[m], th20[m], ())
        ripple = (w10, w20, w20 - w10)
        w10, w20 = -_window_fit(rt[m], th10[m], ripple), -_window_fit(rt[m], th20[m], ripple)
        rows.append(dict(t_mid=float(.5*(lo + hi)), n_mid=float(prof.n(.5*(lo + hi))), w10=float(w10), w20=float(w20)))
    w10 = np.array([r['w10'] for r in rows])
    w20 = np.array([r['w20'] for r in rows])
    dev = dict(w10=float(np.max(np.abs(w10/w10[0] - 1))), w20=float(np.max(np.abs(w20/w20[0] - 1))),
               ratio=float(np.max(np.abs((w20/w10)/(w20[0]/w10[0]) - 1))))
    return dict(
        grid=dict(N=N, L=L, dt=dt, pad=pad, T=prof.T),
        eigenvalues=E.tolist(), x2_of_states=x2.tolist(),
        w10_pre=float(w10[0]), w20_pre=float(w20[0]), w10_post=float(w10[-1]), w20_post=float(w20[-1]),
        splitting_offset=dict(w10=float(w10[0]/(E[1] - E[0]) - 1), w20=float(w20[0]/(E[2] - E[0]) - 1)),
        clock_factor=float(w10[-1]/w10[0]), clock_ratio_pre=float(w20[0]/w10[0]),
        max_deviation=dev, windows=rows,
        passed=bool(max(dev.values()) < 1e-6))
