"""RUT-1 stage 12 (protocol-rut12.md): the units of work of formation from an empty field.

Everything that PRODUCES a raw result lives here. rut12.py assembles, evaluates the gates and reads. Every
threshold is read from the machine-readable block of the protocol; this file holds none. Stage 4R's reciprocal
field, stage 7's population builder and the owner's sampler are used BY IMPORT, unchanged; the angular-momentum
ledger is this stage's own module, `reciprocal_budget.py`.
"""
import json
import os
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut11_tasks as R11       # noqa: E402  (imports stages 10, 9, 8 and 7, one thread per worker)
from rut11_tasks import PO, T0, _find, _safe, _seed, _sha, np   # noqa: E402
import reciprocal as RC         # noqa: E402
import reciprocal_budget as RB  # noqa: E402

CODE = ('rut12.py', 'rut12_tasks.py', 'reciprocal_budget.py') + tuple(R11.CODE)
PROTOCOL = 'protocol-rut12.md'
TAU_FORM = 3*T0


def protocol_path():
    override = os.environ.get('RUT12_PROTOCOL')      # smoke tests on a toy protocol only; recorded in the archive
    return Path(override) if override else _find(PROTOCOL)


def thresholds():
    text = protocol_path().read_text(encoding='utf-8')
    return json.loads(re.search(r'## Declared thresholds.*?```json\s*(\{.*?\})\s*```', text, re.S).group(1))


TH = thresholds()
F = TH['F']


# ---------------------------------------------------------------- the two populations of the experiment
def target_state(name):
    return R11.population(name)['state']


def bare_annulus(st):
    """The same f(E, L) parameters at alpha = 0: an exact equilibrium of the bare point mass, which knows
    nothing about the field it is about to write."""
    return PO.Annulus(L0=st['L0'], dL=st['dL'], dE=st['dE'], w=st['w'], tau_keep=st['tau_keep'],
                      mass=st['mass'], reach=st['reach'], r_lo=st['r'][0], r_hi=st['r'][-1],
                      n_r=len(st['r'])).at_alpha(0.)


def start_state(name, start):
    st = target_state(name)
    return (bare_annulus(st) if start == 'bare' else PO.Annulus.from_state(st)), st


# ---------------------------------------------------------------- structure, measured the same way everywhere
def circular_radii(r_grid, C, L):
    eff = (-1./r_grid - C)[None, :] + .5*(np.asarray(L)[:, None]/r_grid[None, :])**2
    i = np.clip(eff.argmin(1), 1, len(r_grid) - 2)
    k = np.arange(len(L))
    y0, y1, y2 = eff[k, i - 1], eff[k, i], eff[k, i + 1]
    curv = y0 - 2*y1 + y2
    shift = np.where(curv > 0, .5*(y0 - y2)/np.where(curv > 0, curv, 1.), 0.)
    return r_grid[i] + np.clip(shift, -1, 1)*(r_grid[1] - r_grid[0])


def structure(x, v, m, r_grid, C):
    """The quantities the prediction is about: angular momentum, the circular radii it implies, and the
    epicyclic amplitude left over. C is the field the bodies are actually in, on the radial grid."""
    L = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
    wm = m/np.sum(m)
    r = np.linalg.norm(x, axis=1)
    phi = -1./r_grid - C
    E = .5*np.sum(v*v, axis=1) + np.interp(r, r_grid, phi)
    rc = circular_radii(r_grid, C, L)
    Ec = np.interp(rc, r_grid, phi) + .5*(L/rc)**2
    d1 = np.gradient(phi, r_grid)
    kap2 = np.interp(rc, r_grid, np.gradient(d1, r_grid)) + 3*np.interp(rc, r_grid, d1)/rc
    kap = np.sqrt(np.maximum(kap2, 1e-12))
    Er = np.maximum(E - Ec, 0.)
    amp = np.sqrt(2*Er)/kap
    mean = lambda q: float(np.sum(wm*q))
    return dict(mean_L=mean(L), sigma_L=float(np.sqrt(mean((L - mean(L))**2))),
                mean_rc=mean(rc), sigma_rc=float(np.sqrt(mean((rc - mean(rc))**2))),
                rms_epicycle=float(np.sqrt(mean(amp*amp))), mean_kappa=mean(kap),
                mean_radius=mean(r), sigma_radius=float(np.sqrt(mean((r - mean(r))**2))),
                mean_Er=mean(Er), bound=float(np.sum(wm*(E < 0))))


def radial_field(field, r_grid):
    """C(r) of the spectral field, azimuthally averaged, on a declared sampling grid and interpolated onto the
    population's own radial grid. The readout builds a (points, modes, modes) phase tensor, so it is chunked:
    the whole grid at once would ask for a gigabyte."""
    rs = np.linspace(r_grid[0], r_grid[-1], F['sample_radii'])
    ang = np.linspace(0, 2*np.pi, F['ring_points'], endpoint=False)
    pts = np.stack([np.outer(rs, np.cos(ang)).ravel(), np.outer(rs, np.sin(ang)).ravel()], axis=1)
    out = np.concatenate([field.C_and_grad(pts[i:i + F['chunk']])[0]
                          for i in range(0, len(pts), F['chunk'])])
    return np.interp(r_grid, rs, out.reshape(len(rs), len(ang)).mean(axis=1))


def prime_axisymmetric(field, r_grid, sigma):
    """The smooth field the population writes, not a snapshot of point sources.

    The steady state of the h equation for a static source is h = tau_keep alpha W*rho, so with an axisymmetric
    surface density sigma(r) the coefficients follow from its Hankel transform,
    rho_hat(k) = 2 pi int sigma(r) J_0(k r) r dr. Priming instead from the instantaneous positions puts every
    body at the peak of its own Gaussian; the field is then as lumpy as the sample, it relaxes violently, and
    what breaks is the integrator rather than anything physical."""
    from scipy.special import j0
    kk = np.sqrt(field.kx**2 + field.ky**2)
    r = np.asarray(r_grid, float)
    integrand = np.asarray(sigma, float)*r
    rho_hat = 2*np.pi*np.trapezoid(j0(np.outer(kk.ravel(), r))*integrand[None, :], r, axis=1).reshape(kk.shape)
    field.h = field.tau_keep*(field.alpha/field.area)*field.What*rho_hat
    field.hd = np.zeros_like(field.h)


def harmonics(x, m, count=4):
    th = np.arctan2(x[:, 1], x[:, 0])
    tot = float(np.sum(m))
    return [float(abs(np.sum(m*np.exp(-1j*k*th)))/tot) for k in range(1, count + 1)]


# ---------------------------------------------------------------- the prediction, computed before any run
def task_prediction(name, n, seed):
    """What adiabatic invariance says must form: the start's own L carried over, so the start's bodies read in
    the TARGET's field give the circular radii, and the epicyclic amplitude shrinks by sqrt(kappa_b/kappa_f).
    Nothing here uses a formation run."""
    t0 = time.time()
    st = target_state(name)
    a0, at = bare_annulus(st), PO.Annulus.from_state(st)
    r_grid, C = np.asarray(st['r'], float), np.asarray(st['C'], float)
    m = np.full(n, float(st['mass'])/n)
    xb, vb = a0.draw(n, seed=seed)
    xt, vt = at.draw(n, seed=seed)
    start_own = structure(xb, vb, m, r_grid, np.zeros_like(C))
    start_in_target = structure(xb, vb, m, r_grid, C)
    tgt = structure(xt, vt, m, r_grid, C)
    scale = float(np.sqrt(start_own['mean_kappa']/start_in_target['mean_kappa']))
    # and, while the target is in hand: does the spectral field at the declared resolution reproduce the field
    # the population was built with? The two representations must agree before either is used as a control.
    G = TH['G']
    f = RC.ReciprocalField(G['box'], G['modes'], st['w'], st['tau_keep'], TAU_FORM, float(st['alpha']))
    prime_axisymmetric(f, r_grid, np.asarray(at.sigma, float))
    C_spec = radial_field(f, r_grid)
    inside = (r_grid >= F['compare_r'][0]) & (r_grid <= F['compare_r'][1])
    scale_C = float(np.max(np.abs(C)))
    field_check = dict(max_relative=float(np.max(np.abs(C_spec[inside] - C[inside]))/scale_C),
                       rms_relative=float(np.sqrt(np.mean((C_spec[inside] - C[inside])**2))/scale_C),
                       spectral_max=float(np.max(C_spec)), archived_max=scale_C)
    return dict(name=name, bodies=int(n), seed=int(seed), spectral_field_against_the_population=field_check,
                start_in_its_own_field=start_own,
                start_read_in_the_target_field=start_in_target, target=tgt, kappa_scale=scale,
                predicted=dict(mean_L=start_own['mean_L'], sigma_L=start_own['sigma_L'],
                               mean_rc=start_in_target['mean_rc'], sigma_rc=start_in_target['sigma_rc'],
                               rms_epicycle=start_own['rms_epicycle']*scale),
                hotter_than_target_by=float(start_own['rms_epicycle']*scale/tgt['rms_epicycle']),
                seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- part A: the budget itself
def task_budget(step, box, modes, offset, orbits, name=None):
    """The reference run for the ledgers: the declared population, a short horizon, both budgets carried."""
    t0 = time.time()
    A = TH['A']
    name = name or A['population']
    # the reference IS the formation configuration at a short horizon: the equilibrium of the bare point mass
    # and an empty field. Priming a mature field from a few point sources instead would put every body at the
    # bottom of its own Gaussian well, and what then fails is the integrator, not the ledger.
    src, st = start_state(name, A['start'])
    m = np.full(A['bodies'], float(st['mass'])/A['bodies'])
    x, v = src.draw(A['bodies'], seed=A['seed'])
    field = RC.ReciprocalField(box, modes, st['w'], st['tau_keep'], TAU_FORM, float(st['alpha']))
    rows = advance(field, x, v, m, orbits*T0, step, offset=offset, samples=A['samples'])
    last = rows[-1]
    return dict(step=float(step), box=float(box), modes=int(modes), offset=float(offset), orbits=float(orbits),
                name=name, rows=rows, energy_balance=last['energy_balance'], angular_balance=last['angular_balance'],
                L_field=last['L_field'], L_matter=last['L_matter'], L_total_0=last['L_total_0'],
                dissipated=last['dissipated'], torque_integral=last['torque_integral'],
                seconds=round(time.time() - t0, 1))


def advance(field, x, v, m, t_end, step, offset=0., samples=20, r_grid=None, alive=True, record=None):
    """Velocity-Verlet for the matter, each field mode advanced exactly for the midpoint source, with both
    ledgers carried by the trapezoid rule. `alive` false freezes the writing: the alpha = 0 control."""
    x, v = np.array(x, float), np.array(v, float)

    def accel(xx):
        r = np.linalg.norm(xx, axis=1, keepdims=True)
        return -xx/r**3 + field.C_and_grad(xx)[1]

    def H(xx, vv):
        C = field.C_and_grad(xx)[0]
        return float(np.sum(m*(.5*np.sum(vv*vv, axis=1) - 1/np.linalg.norm(xx, axis=1) - C))) + field.field_energy()

    ledger = RB.Ledger(field, x, v, m, H, offset)
    a = accel(x)
    t, rows = 0., []
    while t < t_end - 1e-12:
        s = min(step, t_end - t)
        v_half = v + .5*s*a
        x_mid = x + .5*s*v_half
        x = x + s*v_half
        if alive:
            field.advance(x_mid, m, s)
        a = accel(x)
        v = v_half + .5*s*a
        ledger.step(s)
        t += s
        if len(rows) < samples and t >= (len(rows) + 1)*t_end/samples - 1e-12:
            row = dict(t_periods=t/T0, **ledger.read(x, v))
            if record is not None:
                row.update(record(x, v, field))
            rows.append(row)
    return rows if record is None else (rows, x, v)


# ---------------------------------------------------------------- part B: formation
def task_formation(name, start, seed, variant='base'):
    """From an empty field (or the declared variant) to the horizon, carrying both ledgers and the structure."""
    t0 = time.time()
    G = TH['G']
    src, st = start_state(name, start)
    n = int(G['bodies_by_variant'].get(variant, G['bodies']))
    step = float(G['step'])*(0.5 if variant == 'half_step' else 1.)
    modes = int(G['modes_by_variant'].get(variant, G['modes']))
    box = float(G['box_by_variant'].get(variant, G['box']))
    horizon = float(G['horizon_by_variant'].get(variant, G['horizon']))
    alive = variant != 'no_writing'
    m = np.full(n, float(st['mass'])/n)
    x0, v0 = src.draw(n, seed=seed)
    field = RC.ReciprocalField(box, modes, st['w'], st['tau_keep'], TAU_FORM, float(st['alpha']))
    r_grid = np.asarray(st['r'], float)
    if start == 'mature':
        prime_axisymmetric(field, r_grid, np.asarray(src.sigma, float))

    def record(x, v, fld):
        C = radial_field(fld, r_grid)
        s = structure(x, v, m, r_grid, C)
        return dict(structure=s, harmonics=harmonics(x, m),
                    support=float(np.sum(m*np.interp(np.linalg.norm(x, axis=1), r_grid, C))
                                  / np.sum(m/np.linalg.norm(x, axis=1))),
                    max_C=float(np.max(C)))

    rows, x, v = advance(field, x0, v0, m, horizon*T0, step, samples=G['samples'], alive=alive, record=record)
    C_final = radial_field(field, r_grid)
    return dict(name=name, start=start, seed=int(seed), variant=variant, bodies=n, step=step, modes=modes,
                box=box, horizon=horizon, rows=rows, initial_sha256=_sha(x0, v0),
                final_structure=structure(x, v, m, r_grid, C_final),
                final_C=C_final.tolist(), final_harmonics=harmonics(x, m),
                initial_structure=rows[0]['structure'] if rows else None,
                alpha=float(st['alpha']), seconds=round(time.time() - t0, 1))


def run_task(kind, key, args):
    t0 = time.time()
    fn = dict(prediction=task_prediction, budget=task_budget, formation=task_formation)[kind]
    return kind, key, _safe(fn(*args)), round(time.time() - t0, 1)
