"""RUT-1 stage 10 (protocol-rut10.md): the units of work of the seeded, full-state warm-mode experiment.

A sample with exact four-fold symmetry carries no m = 1, 2 or 3 shot noise, and the simulator's Cartesian field
keeps that symmetry, so an m = 2 mode can only grow from what is deliberately put there. The seed is the
PREDICTED EIGENMODE in the full state: the bodies carry delta f because they have been moved by the mode's own
field, imposed at its predicted complex rate, from eight e-foldings earlier; the force-producing field carries
delta C; and the excitation carries delta E = (s + 1/tau_keep) delta C, which is what a mode growing as e^{st}
requires and what the equilibrium relation E = C/tau_keep would get wrong.

Stage 9's corrected response calculation, stage 7's simulator and population builder and the owner's sampler
are used BY IMPORT, unchanged. Every threshold is read from the protocol.
"""
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut9_tasks as R9        # noqa: E402  (imports stages 8 and 7, which set one thread per worker)
from rut9_tasks import PO, R7, SP, TAU_FORM, WM, _find, _safe, _seed, np   # noqa: E402
from scipy.interpolate import CubicSpline   # noqa: E402
import formation as FM         # noqa: E402
import pairrun as PR           # noqa: E402

T0 = PO.T0
CODE = ('rut10.py', 'rut10_tasks.py') + tuple(R9.CODE)
PROTOCOL = 'protocol-rut10.md'


def protocol_path():
    override = os.environ.get('RUT10_PROTOCOL')             # smoke tests on a toy protocol only; recorded in the archive
    return Path(override) if override else _find(PROTOCOL)


def thresholds():
    text = protocol_path().read_text(encoding='utf-8')
    return json.loads(re.search(r'## Declared thresholds.*?```json\s*(\{.*?\})\s*```', text, re.S).group(1))


TH = thresholds()


def _sha(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a, dtype=float).tobytes())
    return h.hexdigest()


# ---------------------------------------------------------------- whole arrays, saved exactly and reproducibly
def dump_arrays(path, arrays):
    """Little-endian float64, in sorted key order, behind a one-line JSON index; gzip with no timestamp, so the
    same arrays always give the same bytes. JSON text costs three times as much for a field on a grid."""
    import gzip
    index, blobs, offset = {}, [], 0
    for key in sorted(arrays):
        a = np.ascontiguousarray(arrays[key], dtype='<f8')
        index[key] = dict(shape=list(a.shape), offset=offset)
        blobs.append(a.tobytes())
        offset += a.nbytes
    head = json.dumps(index, sort_keys=True).encode('utf-8') + b'\n'
    Path(path).write_bytes(gzip.compress(head + b''.join(blobs), mtime=0))


def load_arrays(path):
    import gzip
    raw = gzip.decompress(Path(path).read_bytes())
    cut = raw.index(b'\n')
    index, body = json.loads(raw[:cut].decode('utf-8')), raw[cut + 1:]
    out = {}
    for key, meta in index.items():
        n = int(np.prod(meta['shape']))
        out[key] = np.frombuffer(body, dtype='<f8', count=n, offset=meta['offset']).reshape(meta['shape']).copy()
    return out


# ---------------------------------------------------------------- the predicted mode
_MODE = {}


def predicted_mode(pop):
    """Stage 9's corrected calculation, by import: the m = 2 root, the mode's field on a fine radial grid and
    its node coefficients, scaled to the declared peak."""
    if pop['name'] not in _MODE:
        resp = R9.response_for(pop['state'], pop['name'], 2, R9.representation())
        got = WM.solve_rectangle(resp, tuple(R9.TH['N']['family_rect']), R9.REP, R9.QUAD)
        s = max((r['s'] for r in got['roots']), key=lambda z: z.real)
        rgrid, psi, c, v, u, ratio = R9.mode_field(resp, s, TH['P']['peak_field'])
        _MODE.clear()
        _MODE[pop['name']] = dict(resp=resp, s=s, rgrid=rgrid, psi=psi, c=c, counts_agree=got['counts_agree'],
                                  located=got['located'])
    return _MODE[pop['name']]


def amplitude_to_eps(pop, mode, target):
    """delta C = eps Re[psi(r) e^{2 i phi}] has ring coefficient C_2 = eps psi(R)/2: the eps that puts |C_2|/C_0
    at `target` on the ring at the population's own mean radius."""
    R = pop['description']['mean_radius']
    disk = mode['resp'].d
    psi_R = complex(np.interp(R, mode['rgrid'], mode['psi'].real), np.interp(R, mode['rgrid'], mode['psi'].imag))
    return float(target*float(disk.Cs(R))/(abs(psi_R)/2.)), float(R)


# ---------------------------------------------------------------- the quiet sample
def quiet_sample(a, quarter, seed):
    """`quarter` bodies from the sampler of record and their images under rotation by 90, 180 and 270 degrees,
    which are exact in floating point: (x, y) -> (-y, x)."""
    x, v = a.draw(quarter, seed=seed)
    rot = lambda q: np.stack([-q[:, 1], q[:, 0]], axis=1)
    xs, vs = [x], [v]
    for _ in range(3):
        xs.append(rot(xs[-1]))
        vs.append(rot(vs[-1]))
    return np.concatenate(xs), np.concatenate(vs)


def sample_for(a, pop, quarter, k, quiet=True):
    seed = _seed('stage10', pop['name'], int(quarter), int(k))
    if quiet:
        return quiet_sample(a, quarter, seed), seed
    return a.draw(4*quarter, seed=seed), seed               # the control: an ordinary draw of the same size


def source_harmonics(x, w, radius, rates, count=6):
    """|S_m|/S_0 on the ring, m = 1..count: the shot noise the bodies would write."""
    ph = 2*np.pi*np.arange(PR.RING_POINTS)/PR.RING_POINTS
    pts = radius*np.stack([np.cos(ph), np.sin(ph)], axis=1)
    d2 = ((pts[:, None, :] - x[None, :, :])**2).sum(-1)
    S = (np.exp(-d2/(2*w*w))*rates[None, :]).sum(1)
    c = np.fft.fft(S)/PR.RING_POINTS
    return [float(abs(c[m])/abs(c[0])) for m in range(1, count + 1)]


def task_Q(pop, quarter, k):
    """The quiet sample IS the population (its moments against the owner's node distribution, in standard
    errors of the INDEPENDENT bodies) and carries no m = 1, 2, 3 source; an ordinary draw is the control."""
    a = PO.Annulus.from_state(pop['state'])
    nodes = a.nodes()
    ref = a.moments(nodes)
    R = pop['description']['mean_radius']
    out = {}
    for label, quiet in (('quiet', True), ('ordinary', False)):
        (x, v), seed = sample_for(a, pop, quarter, k, quiet)
        n = len(x)
        r = np.linalg.norm(x, axis=1)
        vr = np.sum(v*x, axis=1)/r
        L = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
        independent = quarter if quiet else n
        z = dict(R_mean=(r.mean() - ref['R_mean'])/np.sqrt(ref['R_var']/independent),
                 L_mean=(L.mean() - ref['L_mean'])/np.sqrt(ref['L_var']/independent),
                 vr_var=(np.mean(vr*vr) - ref['vr_var'])/np.sqrt(max(ref['vr_fourth'] - ref['vr_var']**2, 1e-300)/independent))
        out[label] = dict(bodies=int(n), seed=seed, z_scores={k2: float(v2) for k2, v2 in z.items()},
                          source_harmonics=source_harmonics(x, a.w, R, np.full(n, a.alpha*a.mass/n)),
                          sha256=_sha(x, v))
    return dict(name=pop['name'], quarter=int(quarter), realization=int(k), **out)


# ---------------------------------------------------------------- the full-state seed
def seed_field(field, mode, eps, static_relation=False):
    """delta C and its gradient on the simulator's grid, and delta E = (s + 1/tau_keep) delta C. With
    `static_relation` the excitation is given the EQUILIBRIUM relation delta E = delta C/tau_keep instead: the
    mistake this construction exists to avoid, kept as a negative control."""
    m, s, resp = 2, mode['s'], mode['resp']
    X, Y = np.meshgrid(field.axis, field.axis, indexing='ij')
    rr = np.maximum(np.hypot(X, Y), 1e-9)
    th = np.arctan2(Y, X)
    rg = np.linspace(1e-6, float(rr.max()) + 1e-6, 6000)
    psi = sum(cb*resp.kernel(rb, rg) for cb, rb in zip(mode['c'], resp.nodes))
    pr, pi = CubicSpline(rg, psi.real), CubicSpline(rg, psi.imag)
    P = pr(rr) + 1j*pi(rr)
    dP = pr.derivative()(rr) + 1j*pi.derivative()(rr)
    e = np.exp(1j*m*th)
    comps = [P*e, (dP*np.cos(th) - 1j*m*P/rr*np.sin(th))*e, (dP*np.sin(th) + 1j*m*P/rr*np.cos(th))*e]
    factor = (1./field.tau_keep) if static_relation else (s + 1./field.tau_keep)
    for k in range(3):
        field.C[k] += eps*comps[k].real
        field.E[k] += eps*(factor*comps[k]).real


def prepared_bodies(pop, mode, x, v, eps):
    """The bodies moved by the mode's own field, imposed at its predicted complex rate from `efoldings` earlier."""
    P = TH['P']
    n = len(x)
    q, p, sign, steps = R9._integrate(mode['resp'].d, mode['psi'], mode['rgrid'], 2, mode['s'], abs(eps), x, v,
                                      P['efoldings'], P['step'])
    half = slice(0, n) if eps > 0 else slice(n, 2*n)
    return q[half], p[half], int(steps)


def build_run(pop, quarter, k, target, variant):
    """One run, ready to advance. `variant`: 'base'; 'half_step'; 'fine_grid'; 'coarse_grid'; 'unprepared_bodies' and
    'static_excitation', the two negative controls of the preparation; 'ordinary_draw', the same seed on a sample
    that is NOT four-fold symmetric. target = 0 is the unseeded run."""
    G = TH['G']
    a = PO.Annulus.from_state(pop['state'])
    (x, v), seed = sample_for(a, pop, quarter, k, quiet=(variant != 'ordinary_draw'))
    n = len(x)
    info = dict(name=pop['name'], bodies=int(n), quarter=int(quarter), realization=int(k), target=float(target),
                variant=variant, sample_seed=seed, sample_sha256=_sha(x, v))
    mode, eps, R = None, 0., float(pop['description']['mean_radius'])
    if target != 0.:
        mode = predicted_mode(pop)
        eps, R = amplitude_to_eps(pop, mode, target)
        if variant != 'unprepared_bodies':
            x, v, steps = prepared_bodies(pop, mode, x, v, eps)
            info['preparation_steps'] = steps
        info.update(eps=eps, root=[float(mode['s'].real), float(mode['s'].imag)],
                    mode_coefficients=[[float(z.real), float(z.imag)] for z in mode['c']],
                    kernel_nodes=[float(z) for z in mode['resp'].nodes])
    run = PR.Run(a, x, v, np.full(n, a.alpha*a.mass/n), live=True, tau_form=TAU_FORM,
                 record_every=G['record_every']*T0, ring_radius=R,
                 h_max=G['h_max']*(.5 if variant == 'half_step' else 1.),
                 eta=G['eta']*(.5 if variant == 'half_step' else 1.))      # the step is min(h_max, eta sqrt(r^3/GM)): halve both
    if variant in G['grids_per_w']:                      # 'fine_grid', and 'coarse_grid', which is G2's negative control
        run.field = FM.MemoryField(PR.HALF_WIDTH, a.w/G['grids_per_w'][variant], a.w, a.tau_keep, TAU_FORM)
        PR.prime(run.field, a)
    if target != 0.:
        seed_field(run.field, mode, eps, static_relation=(variant == 'static_excitation'))
    info.update(ring_radius=R, initial_sha256=_sha(x, v), field_sha256=_sha(run.field.C, run.field.E),
                grid_points=int(run.field.n))
    return run, info, (x, v)


def task_run(pop, quarter, k, target, variant, horizon):
    t0 = time.time()
    run, info, (x, v) = build_run(pop, quarter, k, target, variant)
    out = run.advance(horizon*T0, wall_limit=R7.WALL_LIMIT).result()
    info.update(status=out['status'], t_final_periods=out['t_final_periods'], steps=out['steps'], rows=out['rows'],
                initial_positions=x.tolist(), initial_velocities=v.tolist(), seconds=round(time.time() - t0, 1))
    return info


def task_field_state(pop, target):
    """The complete initial field of a seeded run, which does not depend on the realization: saved whole."""
    run, info, _ = build_run(pop, TH['G']['quarters'][0], 0, target, 'unprepared_bodies')
    return dict(target=float(target), eps=info['eps'], field_sha256=info['field_sha256'],
                C=run.field.C.tolist(), E=run.field.E.tolist(), axis=run.field.axis.tolist())


def task_background_sensitivity(pop):
    """How far the predicted root moves when the equilibrium field C0 is scaled by 1 +- delta with the
    distribution function's parameters held: a first-order number to read a live run's slowly changing
    axisymmetric field against. A uniform scaling is a PROXY for whatever the live field actually does."""
    import copy
    delta = TH['G']['background_delta']
    rows = {}
    for label, factor in (('base', 1.), ('plus', 1. + delta), ('minus', 1. - delta)):
        state = copy.deepcopy(pop['state'])
        state['C'] = (np.asarray(state['C'], float)*factor).tolist()
        resp = R9.response_for(state, f"{pop['name']}:C0x{factor:.6f}", 2, R9.representation())
        got = WM.solve_rectangle(resp, tuple(R9.TH['N']['family_rect']), R9.REP, R9.QUAD)
        s = max((r['s'] for r in got['roots']), key=lambda z: z.real)
        rows[label] = [float(s.real), float(s.imag)]
    ds = complex(*rows['plus']) - complex(*rows['minus'])
    return dict(name=pop['name'], delta=delta, roots=rows,
                d_root_per_unit_relative_change_of_C0=[float(ds.real/(2*delta)), float(ds.imag/(2*delta))])


def run_task(kind, key, args):
    t0 = time.time()
    fn = dict(Q=task_Q, run=task_run, field_state=task_field_state, background_sensitivity=task_background_sensitivity)[kind]
    return kind, key, _safe(fn(*args)), round(time.time() - t0, 1)
