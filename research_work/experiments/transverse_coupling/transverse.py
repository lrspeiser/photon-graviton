"""TF-1: the transverse laws and the exact-property gates G1, G2, G5, G6, G7, G8, G9.

Provenance (see provenance.md): the Lorentz structure of the vector channel is established; the
steering law is proposed here; the Lagrangian family, the ray-speed argument and the velocity-moment
quadrature are derived here from established mathematics; CWC-1's archived field is reused by hash.
"""
import sys
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp, quad, simpson
from scipy.interpolate import RectBivariateSpline
from common import CWC, sha, read

sys.path.insert(0, str(CWC))
import two_dimensional as CWC2D   # noqa: E402  (CWC-1's frozen-snapshot index rays, reused)

IMPACTS = [-2., -1., -.5, .5, 1., 2.]
ARCHIVED_FIELD = CWC/'evidence/spatial2d-v1/primary-final-field.npz'
ARCHIVED_2D = CWC/'evidence/spatial2d-v1/results.json'
G_INDEX = .2                      # CWC-1's primary coupling; the steering coefficient is set equal to it


# ------------------------------------------------------------------ the archived generated field
def load_field():
    z = np.load(ARCHIVED_FIELD)
    return dict(axis=np.array(z['axis']), phi=np.array(z['phi']), sha256=sha(ARCHIVED_FIELD),
                max_abs_phi=float(np.max(np.abs(z['phi']))))


def field_spline(axis, phi):
    return RectBivariateSpline(axis, axis, phi, kx=3, ky=3, s=0)


# ------------------------------------------------------------------ rays in two dimensions
def steering_rays_2d(axis, phi, g=G_INDEX, impacts=IMPACTS, refined=False, zero=False):
    """The steering law for light: direction angle theta with dn/dt = g P_perp(n) grad phi, unit speed by
    construction; x is the independent variable as in CWC-1's ray notes."""
    sp = field_spline(axis, phi*0 if zero else phi)
    rows = []
    for b in impacts:
        def rhs(x, z):
            y, th = z
            fx = float(sp.ev(x, y, dx=1))
            fy = float(sp.ev(x, y, dy=1))
            return [np.tan(th), g*(-np.sin(th)*fx + np.cos(th)*fy)/np.cos(th)]
        sol = solve_ivp(rhs, [-8, 8], [b, 0.], method='DOP853', rtol=2e-11 if refined else 2e-9,
                        atol=2e-13 if refined else 2e-11, max_step=.025 if refined else .1, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        x = np.linspace(-8, 8, 641)
        y, th = sol.sol(x)
        born = g*simpson(sp.ev(x, np.full_like(x, b), dy=1), x=x)
        rows.append(dict(impact=b, angle=float(th[-1]), born_angle=float(born), exit_y=float(y[-1]),
                         x=x[::8].tolist(), y=y[::8].tolist()))
    return rows


def index_rays_2d(axis, phi, g=G_INDEX, refined=False, zero=False):
    """CWC-1's own Hamiltonian index rays on the same field (reused code), for the live comparison."""
    return CWC2D.rays(SimpleNamespace(axis=axis, phi=phi, g=g), refined=refined, zero=zero)


def compare_rays(steering, index):
    a = np.array([r['angle'] for r in steering])
    b = np.array([r['angle'] for r in index])
    scale = max(float(np.max(np.abs(b))), 1e-30)
    return dict(max_bend_steering=float(np.max(np.abs(a))), max_bend_index=scale,
                relative_difference=float(np.max(np.abs(a - b))/scale),
                all_inward=bool(np.all(a*np.array(IMPACTS) < 0)),
                mirror_relative=float(np.max(np.abs(a + a[::-1]))/max(float(np.max(np.abs(a))), 1e-30)),
                signed_inward=(-a*np.sign(IMPACTS)).tolist())


# ------------------------------------------------------------------ G8: the solenoidal model field
A0_WEAK = .005                    # keeps every deflection in the small-angle regime the focusing metric assumes


def solenoidal_field(x, A0=A0_WEAK, w=1.):
    """B = curl(A_phi phi-hat) with A_phi = A0 R exp(-(R^2+z^2)/2w^2): a fixture, not a source model."""
    x = np.atleast_2d(np.asarray(x, float))
    X, Y, Z = x[:, 0], x[:, 1], x[:, 2]
    R = np.hypot(X, Y)
    e = np.exp(-(R*R + Z*Z)/(2*w*w))
    BR = A0*R*Z/w**2*e
    Bz = A0*e*(2 - R*R/w**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        cx = np.where(R > 0, X/R, 0.)
        cy = np.where(R > 0, Y/R, 0.)
    return np.stack([BR*cx, BR*cy, Bz], axis=-1)


def divergence_check(n=200, seed=3):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-2, 2, (n, 3))
    h = 1e-5
    div = np.zeros(n)
    for d in range(3):
        e = np.zeros(3)
        e[d] = h
        div += (solenoidal_field(pts + e)[:, d] - solenoidal_field(pts - e)[:, d])/(2*h)
    scale = np.max(np.linalg.norm(solenoidal_field(pts), axis=1))
    return float(np.max(np.abs(div))/scale)


def toward_axis_force(x, A0=A0_WEAK, w=1.):
    """A force of the same magnitude as B, pointing at the z-axis: the well and steering comparison."""
    x = np.atleast_2d(np.asarray(x, float))
    mag = np.linalg.norm(solenoidal_field(x, A0, w), axis=-1)
    R = np.maximum(np.hypot(x[:, 0], x[:, 1]), 1e-300)
    rho = np.stack([x[:, 0]/R, x[:, 1]/R, np.zeros_like(R)], -1)
    return -mag[:, None]*rho


def perpendicular_frame(d):
    d = np.asarray(d, float)/np.linalg.norm(d)
    trial = np.array([0., 1., 0.]) if abs(d[1]) < .9 else np.array([1., 0., 0.])
    e1 = trial - np.dot(trial, d)*d
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    return d, e1, e2


def rays_3d(law, direction, offsets, g=1.):
    """Rays of unit speed through the model field under 'vector' (dn/dt = n x B) or 'transverse'
    (dn/dt = P_perp(n) F, the transverse part of the toward-axis force; the same for the steering law and
    for the light of a well up to its factor). Returns deflections and the focusing fraction."""
    d = np.asarray(direction, float)/np.linalg.norm(direction)
    rows = []
    for off in offsets:
        off = np.asarray(off, float)
        x0 = -8*d + off

        def rhs(t, z):
            x, n = z[:3], z[3:]
            if law == 'vector':
                dn = np.cross(n, solenoidal_field(x)[0])
            else:
                F = g*toward_axis_force(x)[0]
                dn = F - n*np.dot(n, F)
            return np.concatenate([n, dn])
        sol = solve_ivp(rhs, [0, 16], np.concatenate([x0, d]), method='DOP853', rtol=1e-11, atol=1e-13)
        if not sol.success:
            raise RuntimeError(sol.message)
        n_end = sol.y[3:, -1]
        delta = n_end - d
        p = off - np.dot(off, d)*d
        toward = -np.array([p[0], p[1], 0.])
        toward = toward - np.dot(toward, d)*d          # a ray can only deflect transversely to its own motion
        toward /= np.linalg.norm(toward)
        mag = float(np.linalg.norm(delta))
        rows.append(dict(offset=off.tolist(), deflection=delta.tolist(), magnitude=mag,
                         cosine=float(np.dot(delta, toward)/max(mag, 1e-300)),
                         speed_drift=float(np.max(np.abs(np.linalg.norm(sol.y[3:], axis=0) - 1)))))
    mags = np.array([r['magnitude'] for r in rows])
    cos = np.array([r['cosine'] for r in rows])
    return dict(direction=d.tolist(), rows=rows,
                focusing_fraction=float(np.sum(mags*cos)/max(float(np.sum(mags)), 1e-300)),
                max_speed_drift=float(max(r['speed_drift'] for r in rows)),
                max_deflection=float(np.max(mags)))


def focusing_geometry():
    views = dict(along_axis=[0., 0., 1.], in_plane=[1., 0., 0.], oblique=[1., 0., 1.])
    out = {}
    for name, d in views.items():
        d, e1, e2 = perpendicular_frame(d)
        offsets = [b*e1 for b in IMPACTS] + ([b*e2 for b in IMPACTS] if name == 'oblique' else [])
        out[name] = {law: rays_3d(law, d, offsets) for law in ('vector', 'transverse')}
    C_KMS = 299792.458

    def global_fraction(law):
        rows = [r for v in out.values() for r in v[law]['rows']]
        mags = np.array([r['magnitude'] for r in rows])
        cos = np.array([r['cosine'] for r in rows])
        return float(np.sum(mags*cos)/max(float(np.sum(mags)), 1e-300))
    return dict(views=out, divergence_relative=divergence_check(), field_amplitude=A0_WEAK,
                max_deflection=float(max(v[l]['max_deflection'] for v in out.values() for l in v)),
                focusing_fraction=dict(vector=global_fraction('vector'),
                                       transverse=float(np.min([v['transverse']['focusing_fraction'] for v in out.values()])),
                                       vector_max_abs=abs(global_fraction('vector')),
                                       per_view={k: {l: v[l]['focusing_fraction'] for l in v} for k, v in out.items()}),
                max_speed_drift=float(max(v[l]['max_speed_drift'] for v in out.values() for l in v)),
                light_to_matter_ratio_vector={f'{v:.0f}_km_s': C_KMS/v for v in (200., 300.)})


# ------------------------------------------------------------------ G6: the Lagrangian family
def lagrangian_family(alpha, n_states=20, seed=1):
    """L = |v|^2/2 + psi(x)|v|^alpha: the Euler-Lagrange acceleration, checked against finite differences of
    L along the motion, and its parallel and perpendicular parts against the formulas derived in the protocol."""
    rng = np.random.default_rng(seed)
    k = np.array([.3, -.2, .5])

    def psi(x):
        return .25 + .2*np.exp(-np.dot(x, x)/2) + .05*np.dot(k, x)

    def grad_psi(x):
        return -.2*x*np.exp(-np.dot(x, x)/2) + .05*k

    def L(x, v):
        return .5*np.dot(v, v) + psi(x)*np.linalg.norm(v)**alpha

    def d4(f, eps):
        """fourth-order central difference: truncation eps^4, round-off 1e-16/eps"""
        return (8*(f(eps) - f(-eps)) - (f(2*eps) - f(-2*eps)))/(12*eps)

    def dL_dv(x, v, eps=1e-3):
        out = np.zeros(3)
        for i in range(3):
            e = np.zeros(3)
            e[i] = 1.
            out[i] = d4(lambda h: L(x, v + h*e), eps)
        return out

    def dL_dx(x, v, eps=1e-3):
        out = np.zeros(3)
        for i in range(3):
            e = np.zeros(3)
            e[i] = 1.
            out[i] = d4(lambda h: L(x + h*e, v), eps)
        return out

    worst_el, worst_par, worst_perp, speed_rate = 0., 0., 0., 0.
    for _ in range(n_states):
        x = rng.normal(size=3)
        v = rng.normal(size=3)
        v *= rng.uniform(.5, 2.)/np.linalg.norm(v)
        s = np.linalg.norm(v)
        vh = v/s
        p = psi(x)
        gp = grad_psi(x)
        M = np.eye(3) + alpha*p*s**(alpha - 2)*(np.eye(3) + (alpha - 2)*np.outer(vh, vh))
        rhs = s**alpha*gp - alpha*np.dot(gp, v)*s**(alpha - 2)*v
        a = np.linalg.solve(M, rhs)
        a_par = np.dot(a, vh)
        a_perp = a - a_par*vh
        f_par = (1 - alpha)*s**alpha*np.dot(gp, vh)/(1 + alpha*(alpha - 1)*p*s**(alpha - 2))
        f_perp = s**alpha*(gp - np.dot(gp, vh)*vh)/(1 + alpha*p*s**(alpha - 2))
        dt = 1e-3
        dP = d4(lambda h: dL_dv(x + v*h, v + a*h), dt)      # time derivative along the motion, fourth order
        el = np.linalg.norm(dP - dL_dx(x, v))/max(np.linalg.norm(dL_dx(x, v)), 1e-12)
        na = max(np.linalg.norm(a), 1e-12)
        worst_el = max(worst_el, el)
        worst_par = max(worst_par, abs(a_par - f_par)/na)
        worst_perp = max(worst_perp, np.linalg.norm(a_perp - f_perp)/na)
        speed_rate = max(speed_rate, abs(a_par)/na)
    return dict(alpha=alpha, euler_lagrange_residual=worst_el, parallel_formula_error=worst_par,
                perpendicular_formula_error=worst_perp, speed_change_fraction=speed_rate)


# ------------------------------------------------------------------ G7: velocity averages
from support import support_factor   # noqa: E402


def velocity_averages(N=1_000_000, seed=7):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(N, 3))
    B = np.array([0., 0., 1.])
    iso = np.mean(np.cross(v, B), axis=0)
    rot = np.mean(np.cross(v*.1 + np.array([0., 1., 0.]), B), axis=0)   # co-rotating at speed 1, dispersion 0.1
    betas = [-2., -1., 0., .45]
    mc, qd = [], []
    for beta in betas:
        q = 1 - beta
        x = rng.normal(size=N)
        y = rng.normal(size=N)*np.sqrt(q)
        z = rng.normal(size=N)*np.sqrt(q)
        mc.append(1 - float(np.mean(x*x/(x*x + y*y + z*z))))
        qd.append(support_factor(beta))
    return dict(N=N, isotropic_mean_abs=float(np.linalg.norm(iso)), standard_error=float(1/np.sqrt(N)),
                rotating_mean=rot.tolist(), rotating_error=float(abs(rot[0] - 1.) + abs(rot[1]) + abs(rot[2])),
                betas=betas, support_quadrature=qd, support_monte_carlo=mc,
                support_max_difference=float(np.max(np.abs(np.array(qd) - np.array(mc)))),
                isotropic_exact_error=float(abs(support_factor(0.) - 2/3)))


# ------------------------------------------------------------------ G9: orbits in the spherical steering field
def g_field(r, g0=1.):
    return g0*r/(r + 1)**2


def acceleration(x, v, law, g0=1.):
    r = np.linalg.norm(x, axis=-1)
    rhat = x/np.maximum(r, 1e-300)[:, None]
    F = -g_field(r, g0)[:, None]*rhat
    if law == 'well':
        return F
    s = np.linalg.norm(v, axis=-1)
    vh = v/np.maximum(s, 1e-300)[:, None]
    F = F - np.sum(F*vh, axis=1)[:, None]*vh
    return np.where(s[:, None] > 0, F, 0.)


def rk4(x, v, dt, steps, law, record_every=None):
    rec = []
    for k in range(steps):
        k1x, k1v = v, acceleration(x, v, law)
        k2x, k2v = v + .5*dt*k1v, acceleration(x + .5*dt*k1x, v + .5*dt*k1v, law)
        k3x, k3v = v + .5*dt*k2v, acceleration(x + .5*dt*k2x, v + .5*dt*k2v, law)
        k4x, k4v = v + dt*k3v, acceleration(x + dt*k3x, v + dt*k3v, law)
        x = x + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
        v = v + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        if record_every and (k + 1) % record_every == 0:
            rec.append((x.copy(), v.copy()))
    return x, v, rec


def circular_and_radial(g0=1.):
    out = {}
    v0 = np.sqrt(g_field(1., g0))          # circular speed at r = 1
    T = 2*np.pi/v0
    for law in ('steering', 'well'):
        def rhs(t, z):
            return np.concatenate([z[3:], acceleration(z[None, :3], z[None, 3:], law, g0)[0]])
        sol = solve_ivp(rhs, [0, 10*T], [1., 0., 0., 0., v0, 0.], method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
        t = np.linspace(0, 10*T, 20001)
        z = sol.sol(t)
        r = np.linalg.norm(z[:3], axis=0)
        ang = np.unwrap(np.arctan2(z[1], z[0]))
        period = 10*T*2*np.pi/ang[-1]        # ten turns in 10 T if the period is T
        out[law] = dict(radius_drift=float(np.max(np.abs(r - 1))), period_over_analytic=float(period/T),
                        speed_drift=float(np.max(np.abs(np.linalg.norm(z[3:], axis=0) - v0))))
    # a radial orbit through the centre under steering: starts at r = 3 moving inward at speed 0.3
    def rhs_s(t, z):
        return np.concatenate([z[3:], acceleration(z[None, :3], z[None, 3:], 'steering', g0)[0]])
    sol = solve_ivp(rhs_s, [0, 20], [3., 0., 0., -.3, 0., 0.], method='DOP853', rtol=1e-12, atol=1e-14)
    v_end = sol.y[3:, -1]
    out['radial_steering'] = dict(direction_change=float(np.linalg.norm(v_end/np.linalg.norm(v_end) - np.array([-1., 0., 0.]))),
                                  speed_change=float(abs(np.linalg.norm(v_end) - .3)), final_x=float(sol.y[0, -1]))
    return out


def ensemble(N=2000, dynamical_times=20, dt=.005, seed=11, g0=1.):
    rng = np.random.default_rng(seed)
    u = rng.uniform(size=(N, 3))
    r = 2*u[:, 0]**(1/3)
    cth = 2*u[:, 1] - 1
    sth = np.sqrt(1 - cth*cth)
    ph = 2*np.pi*u[:, 2]
    x0 = np.stack([r*sth*np.cos(ph), r*sth*np.sin(ph), r*cth], 1)
    d = rng.normal(size=(N, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    v0 = d*(.7*np.sqrt(g_field(r, g0)*r))[:, None]
    T = 2*np.pi/np.sqrt(g_field(1., g0))
    steps = int(round(dynamical_times*T/dt))
    out = {}
    for law in ('steering', 'well'):
        x, v, _ = rk4(x0.copy(), v0.copy(), dt, steps, law)
        rr = np.linalg.norm(x, axis=1)
        bound = rr < 10
        if law == 'steering':
            drift = np.max(np.abs(np.linalg.norm(v, axis=1)/np.linalg.norm(v0, axis=1) - 1))
        else:
            pot = lambda s: g0*(np.log(s + 1) + 1/(s + 1))
            E0 = .5*np.sum(v0*v0, 1) + pot(np.linalg.norm(x0, axis=1))
            E1 = .5*np.sum(v*v, 1) + pot(rr)
            drift = np.max(np.abs(E1 - E0)/np.abs(E0))

        def anisotropy(xx, vv):
            rh = xx/np.linalg.norm(xx, axis=1)[:, None]
            vr = np.sum(vv*rh, 1)
            vt2 = np.sum(vv*vv, 1) - vr*vr
            return float(1 - np.mean(vt2)/(2*np.mean(vr*vr)))
        out[law] = dict(bound_fraction=float(np.mean(bound)), conservation_drift_per_orbit=float(drift/dynamical_times),
                        anisotropy_initial=anisotropy(x0, v0), anisotropy_final_bound=anisotropy(x[bound], v[bound]),
                        median_radius_initial=float(np.median(r)), median_radius_final_bound=float(np.median(rr[bound])),
                        max_radius_final=float(np.max(rr)))
    return dict(N=N, dynamical_times=dynamical_times, dt=dt, laws=out)


# ------------------------------------------------------------------ G1 and G2 on the archived field
def probes_on_field(axis, phi, g=G_INDEX, b_clock=2):
    """Two probes at one position and velocity with internal action over mass differing by two: the steering
    and vector accelerations are independent of the action; CWC-1's material force is the control."""
    sp = field_spline(axis, phi)
    x, y = 0., 1.
    grad = np.array([float(sp.ev(x, y, dx=1)), float(sp.ev(x, y, dy=1))])
    ph = float(sp.ev(x, y))
    v = np.array([.3, .4])
    vh = v/np.linalg.norm(v)
    steer = g*(grad - np.dot(grad, vh)*vh)
    Bz = .5
    vec = np.array([v[1]*Bz, -v[0]*Bz])
    probes = [dict(M=1., I=.0005), dict(M=1., I=.001)]
    out = []
    for p in probes:
        omega = np.exp(-b_clock*g*ph)
        cwc = b_clock*g*p['I']*omega*grad/p['M']
        out.append(dict(**p, steering=steer.tolist(), vector=vec.tolist(), cwc1_material=cwc.tolist()))
    diff = lambda key: float(np.max(np.abs(np.array(out[0][key]) - np.array(out[1][key]))))
    return dict(probes=out, steering_difference=diff('steering'), vector_difference=diff('vector'),
                cwc1_ratio=float(np.linalg.norm(out[1]['cwc1_material'])/np.linalg.norm(out[0]['cwc1_material'])))


def index_speed_control(axis, phi, g=G_INDEX):
    """CWC-1's index rays travel at a = exp(-g phi): the largest deviation from unit speed in the field."""
    return float(np.max(np.abs(np.exp(-g*phi) - 1)))
