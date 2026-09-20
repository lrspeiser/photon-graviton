"""RUT-1 stage 11 (protocol-rut11.md): the units of work of the quiet region.

Everything that PRODUCES a raw result lives here, as in stages 7 to 10: rut11.py assembles, evaluates the gates
and reads. Every threshold is read from the machine-readable block of the protocol; this file holds none.
Stage 9's response, stage 10's quiet sample, stage 7's simulator and population builder and the owner's sampler
are used BY IMPORT, unchanged.
"""
import json
import os
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut10_tasks as R10      # noqa: E402  (imports stages 9, 8 and 7, which set one thread per worker)
from rut10_tasks import PO, PR, R7, R9, T0, TAU_FORM, WM, _find, _safe, _seed, _sha, np   # noqa: E402
from scipy.interpolate import CubicSpline    # noqa: E402
from scipy.special import ive                # noqa: E402
import formation as FM         # noqa: E402
import resonant_response as RR  # noqa: E402

CODE = ('rut11.py', 'rut11_tasks.py', 'resonant_response.py') + tuple(R10.CODE)
PROTOCOL = 'protocol-rut11.md'


def protocol_path():
    override = os.environ.get('RUT11_PROTOCOL')     # smoke tests on a toy protocol only; recorded in the archive
    return Path(override) if override else _find(PROTOCOL)


def thresholds():
    text = protocol_path().read_text(encoding='utf-8')
    return json.loads(re.search(r'## Declared thresholds.*?```json\s*(\{.*?\})\s*```', text, re.S).group(1))


TH = thresholds()
REP = TH['R']
QUAD = dict(R9.QUAD)
QUAD.update(REP['contour'])


# ---------------------------------------------------------------- populations
def population(name):
    """A population by name: stage 9's archived eight, or a family B member built at the declared width."""
    archived = R9.archived_populations()
    if name in archived:
        p = archived[name]
        return p if 'state' in p else dict(name=name, state=p)
    mo = re.fullmatch(r'B_dE([0-9.]+)', name)
    if not mo:
        raise KeyError(name)
    built = R9.task_family_population(float(mo.group(1)))
    built['name'] = name                      # stage 9's own name rounds to three decimals; keep the asked-for one
    return built


def mean_radius(state):
    a = PO.Annulus.from_state(state)
    return float(np.trapezoid(2*np.pi*a.r*a.sigma*a.r, a.r))


# ---------------------------------------------------------------- the response, and the analytic rule on it
_RESP, _PAIR = {}, {}


def response_for(name, m, cfg=None):
    cfg = cfg or REP['base']
    key = (name, int(m), json.dumps(cfg, sort_keys=True))
    if key not in _RESP:
        _RESP.clear()
        _PAIR.clear()
        state = population(name)['state']
        _RESP[key] = RR.unpruned(WM.Disk(state, TAU_FORM), cfg, int(m))
    return _RESP[key]


def pair_for(name, m, harmonics, refines, cfg=None):
    """The analytic rule on the declared harmonics, cached within a worker: building it is a tensor
    interpolation of a transform stage 9 already did, so the cost is in the orbits, not here."""
    resp = response_for(name, m, cfg)
    key = (id(resp), tuple(sorted(harmonics)), tuple(refines))
    if key not in _PAIR:
        _PAIR.clear()
        p = RR.Pair(resp, (cfg or REP['base'])['n_L'], (cfg or REP['base'])['n_u'], refines, REP['pattern'])
        p.add(sorted(harmonics))
        _PAIR[key] = p
    return _PAIR[key]


def rule_for(name, m, rect, refines=None, cfg=None):
    resp = response_for(name, m, cfg)
    refines = tuple(refines or REP['refines'])
    harmonics = RR.resonant_harmonics(resp, rect[2], rect[3], REP['resonant_margin'])
    return pair_for(name, m, harmonics, refines, cfg), harmonics


def _rect_row(rule, rect, spread_points=9):
    """One rectangle: the located roots, the winding number, and the error the rule declares on the contour."""
    got = WM.solve_rectangle(rule, tuple(rect), R9.REP, QUAD)
    edge = [complex(rect[0], y) for y in np.linspace(rect[2], rect[3], spread_points)] \
        + [complex(x, .5*(rect[2] + rect[3])) for x in np.linspace(rect[0], rect[1], 3)]
    spread = [rule.spread(s) for s in edge] if hasattr(rule, 'spread') else [float('nan')]
    return dict(rect=[float(x) for x in rect], located=got['located'], winding=got['winding']['count'],
                counts_agree=bool(got['counts_agree']), gap=float(got['gap']), contour_nodes=int(got['nodes']),
                roots=[[float(r['s'].real), float(r['s'].imag)] for r in got['roots']],
                residuals=[float(r['relative_residual']) for r in got['roots']],
                max_spread=float(np.max(spread)), spread_at_floor=float(spread[0]))


def strip_rectangles(m, shifted=False):
    """The declared strip for one m in the LOW box: Re s from the new floor to stage 8's."""
    L = TH['L']
    lo = -(L['strip_im_per_m']*m + L['strip_im_offset'])
    hi, h = L['strip_im_retrograde'], L['rectangle_height']
    start = lo - (.5*h if shifted else 0.)
    count = int(np.ceil((hi - start)/h - 1e-9))
    return [(REP['gamma_min'], REP['gamma_stage8'], start + k*h, start + (k + 1)*h) for k in range(count)]


# ---------------------------------------------------------------- part A: the rule itself
def task_A1():
    """An integral with the shape of the response's own and a reference good to fourteen digits: the rule at
    the declared refinements against it, and stage 9's quadrature against it, at every declared growth rate."""
    A = TH['A']['toy']
    out = []
    for gam in A['growth_rates']:
        z = complex(A['omega'], gam)             # z = i s, so Im z is the growth rate
        ref = RR.reference_integral(z, A['a'], A['b'], A['c'], A['dL'], A['dE'], A['reach'])
        vals = {k: RR.reference_by_rule(z, A['a'], A['b'], A['c'], A['dL'], A['dE'], A['reach'], A['n_L'], A['n_u'], k)
                for k in (REP['refines'] + REP['control_refines'])}
        r1 = vals[REP['refines'][1]] + (vals[REP['refines'][1]] - vals[REP['refines'][0]])/3.
        r2 = vals[REP['control_refines'][1]] + (vals[REP['control_refines'][1]] - vals[REP['control_refines'][0]])/3.
        gauss = RR.reference_by_rule(z, A['a'], A['b'], A['c'], A['dL'], A['dE'], A['reach'], A['n_L'], A['n_u'], 1,
                                     gauss=True)
        out.append(dict(growth_rate=float(gam), reference=[ref.real, ref.imag],
                        pair_error=float(abs(r1 - ref)/abs(ref)), refined_pair_error=float(abs(r2 - ref)/abs(ref)),
                        gauss_legendre_error=float(abs(gauss - ref)/abs(ref)),
                        single_errors={str(k): float(abs(v - ref)/abs(ref)) for k, v in vals.items()}))
    return dict(rows=out)


def task_A2(name, m, rect):
    """The same root, by stage 9's quadrature and by the analytic rule, in a rectangle where both are valid;
    and the rule's own convergence there."""
    t0 = time.time()
    resp9 = R9.response_for(population(name)['state'], name, m, R9.representation())
    got9 = WM.solve_rectangle(resp9, tuple(rect), R9.REP, R9.QUAD)
    s9 = max((r['s'] for r in got9['roots']), key=lambda z: z.real, default=None)
    rule, harmonics = rule_for(name, m, rect)
    row = _rect_row(rule, rect)
    s11 = max((complex(*r) for r in row['roots']), key=lambda z: z.real, default=None)
    fine, _ = rule_for(name, m, rect, refines=REP['control_refines'])
    rowf = _rect_row(fine, rect)
    sf = max((complex(*r) for r in rowf['roots']), key=lambda z: z.real, default=None)
    out = dict(name=name, m=int(m), rect=[float(x) for x in rect], harmonics=harmonics,
               stage9_root=None if s9 is None else [s9.real, s9.imag], stage9_located=got9['located'],
               analytic=row, refined=rowf, seconds=round(time.time() - t0, 1))
    if s9 is not None and s11 is not None:
        out['difference_from_stage9'] = float(abs(s11 - s9)/abs(s9))
        out['growth_difference_from_stage9'] = float((s11.real - s9.real)/s9.real)
    if sf is not None and s11 is not None:
        out['difference_under_refinement'] = float(abs(sf - s11)/abs(s11))
    return out


# ---------------------------------------------------------------- part L and part M: the search
def task_box(name, m, index, shifted):
    t0 = time.time()
    rect = strip_rectangles(int(m), bool(shifted))[int(index)]
    rule, harmonics = rule_for(name, m, rect)
    row = _rect_row(rule, rect)
    row.update(name=name, m=int(m), index=int(index), shifted=bool(shifted), harmonics=harmonics,
               pole_margin=float(rule.pole_margin(rect[2], rect[3])), seconds=round(time.time() - t0, 1))
    return row


def task_box_control(name, m, index, shifted):
    """The same rectangle under stage 9's quadrature, whose poles lie inside it: the negative control."""
    t0 = time.time()
    rect = strip_rectangles(int(m), bool(shifted))[int(index)]
    resp = R9.response_for(population(name)['state'], name, m, R9.representation())
    got = WM.solve_rectangle(resp, tuple(rect), R9.REP, R9.QUAD)      # stage 9's rule AND stage 9's own contour
    return dict(name=name, m=int(m), index=int(index), shifted=bool(shifted), rect=[float(x) for x in rect],
                located=got['located'], winding=got['winding']['count'], counts_agree=bool(got['counts_agree']),
                roots=[[float(r['s'].real), float(r['s'].imag)] for r in got['roots']],
                seconds=round(time.time() - t0, 1))


def task_norm(name, m, refined=False):
    """max ||alpha H(s) M(s)|| over the declared strip on a declared grid. Where that stays below one, T = 1 +
    alpha H M cannot be singular and the population has no mode there at any resolution; the bound is only as
    good as the grid, so a refined grid is run as its control."""
    t0 = time.time()
    M = TH['M']
    L = TH['L']
    im_lo, im_hi = -(L['strip_im_per_m']*m + L['strip_im_offset']), L['strip_im_retrograde']
    f = M['refine_factor'] if refined else 1
    gammas = np.concatenate([[REP['gamma_min']],
                             np.geomspace(REP['gamma_stage8'], M['re_max'], f*M['n_re'])])
    omegas = np.linspace(im_lo, im_hi, f*int(np.ceil((im_hi - im_lo)/M['im_spacing'])) + 1)
    best, at, rows = 0., None, 0
    resp = response_for(name, m)
    rule_cache = {}
    for om in omegas:
        harmonics = tuple(RR.resonant_harmonics(resp, om, om, REP['resonant_margin'])) if M['analytic_at_floor'] else ()
        for gam in gammas:
            s = complex(gam, om)
            if gam <= REP['gamma_stage8'] and harmonics:
                if harmonics not in rule_cache:
                    rule_cache.clear()
                    rule_cache[harmonics] = pair_for(name, m, harmonics, REP['refines'])
                red = rule_cache[harmonics].reduced(s)
            else:
                red = resp.reduced(s)
            n = float(np.linalg.norm(resp.d.alpha*resp.d.H(s)*red, 2))
            rows += 1
            if n > best:
                best, at = n, s
    return dict(name=name, m=int(m), refined=bool(refined), points=rows, max_norm=best,
                at=[float(at.real), float(at.imag)], n_re=len(gammas), n_im=len(omegas),
                seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- part T: where the family's threshold is
def ladder_name(dE):
    """A name that round-trips: stage 9's `family_name` rounds to three decimals, which would collapse the
    declared widths of the ladder onto two names and build the wrong populations."""
    return 'B_dE%g' % float(dE)


def task_threshold(dE, m=2):
    """One family B member, built as stage 9 built them, searched from the new floor upward for its fastest
    root: the family's boundary bracketed by resolution instead of extrapolated."""
    t0 = time.time()
    name = ladder_name(dE)
    pop = population(name)
    rect = (REP['gamma_min'], TH['T']['re_max'], TH['T']['im_lo'], TH['T']['im_hi'])
    rule, harmonics = rule_for(name, m, rect)
    row = _rect_row(rule, rect)
    s = max((complex(*r) for r in row['roots']), key=lambda z: z.real, default=None)
    out = dict(name=name, dE=float(dE), m=int(m), harmonics=harmonics, mean_radius=mean_radius(pop['state']),
               alpha=float(pop['state']['alpha']), seconds=round(time.time() - t0, 1), **row)
    if s is not None:
        fine, _ = rule_for(name, m, rect, refines=REP['control_refines'])
        rf = _rect_row(fine, rect)
        sf = max((complex(*r) for r in rf['roots']), key=lambda z: z.real, default=None)
        out['refined_root'] = None if sf is None else [sf.real, sf.imag]
        out['difference_under_refinement'] = None if sf is None else float(abs(sf - s)/abs(s))
        out['e_folding_periods'] = float(1./(s.real*T0))
    return out


# ---------------------------------------------------------------- part D: the controlled disturbance
def disturbance_field(field, R, w, m, eps, relation):
    """delta C = eps Re[k_m(R, r) e^{i m theta}], the memory kernel centred on the population's own mean
    radius: a declared shape with no free parameter but its amplitude. `relation` sets the excitation the
    memory is given: 'settled' is delta E = delta C/tau_keep, the state a slowly written field would be in;
    'none' leaves the excitation alone."""
    X, Y = np.meshgrid(field.axis, field.axis, indexing='ij')
    rr = np.maximum(np.hypot(X, Y), 1e-9)
    th = np.arctan2(Y, X)
    rg = np.linspace(1e-6, float(rr.max()) + 1e-6, 6000)
    sp = CubicSpline(rg, np.exp(-(rg - R)**2/(2*w*w))*ive(m, rg*R/(w*w)))
    P, dP = sp(rr), sp.derivative()(rr)
    e = np.exp(1j*m*th)
    comps = [P*e, (dP*np.cos(th) - 1j*m*P/rr*np.sin(th))*e, (dP*np.sin(th) + 1j*m*P/rr*np.cos(th))*e]
    factor = {'settled': 1./field.tau_keep, 'none': 0.}[relation]
    for k in range(3):
        field.C[k] += eps*comps[k].real
        field.E[k] += eps*(factor*comps[k]).real


def folded_sample(a, quarter, seed, fold):
    """`quarter` bodies from the sampler of record and their images under rotation by 2 pi / fold, which are
    exact in floating point for fold 2 and 4. A fold of 4 removes the sample's own m = 1, 2 and 3; a fold of 2
    removes its odd harmonics, which is what an odd m needs."""
    x, v = a.draw(quarter, seed=seed)
    rot = (lambda q: np.stack([-q[:, 1], q[:, 0]], axis=1)) if fold == 4 else (lambda q: -q)
    xs, vs = [x], [v]
    for _ in range(fold - 1):
        xs.append(rot(xs[-1]))
        vs.append(rot(vs[-1]))
    return np.concatenate(xs), np.concatenate(vs)


def build_disturbance(name, m, quarter, k, target, variant, pop=None):
    D = TH['D']
    pop = pop or population(name)
    a = PO.Annulus.from_state(pop['state'])
    fold = 4 if m % 2 == 0 else 2
    seed = _seed('stage11', name, int(quarter), int(k), int(m))
    if variant == 'ordinary_draw':                 # the control: the same body count, but not folded
        x, v = a.draw(fold*int(quarter), seed=seed)
    else:
        x, v = folded_sample(a, quarter, seed, fold)
    n = len(x)
    R = mean_radius(pop['state'])
    C0 = float(a.C[int(np.argmin(np.abs(a.r - R)))])
    psi_R = float(np.exp(0.)*ive(m, R*R/(a.w*a.w)))
    eps = float(target*C0/(psi_R/2.))
    half = variant == 'half_step'
    run = PR.Run(a, x, v, np.full(n, a.alpha*a.mass/n), live=True, tau_form=TAU_FORM,
                 record_every=D['record_every']*T0, ring_radius=R,
                 h_max=D['h_max']*(.5 if half else 1.), eta=D['eta']*(.5 if half else 1.))
    if variant in D['grids_per_w']:
        run.field = FM.MemoryField(PR.HALF_WIDTH, a.w/D['grids_per_w'][variant], a.w, a.tau_keep, TAU_FORM)
        PR.prime(run.field, a)
    if target != 0.:
        disturbance_field(run.field, R, a.w, m, eps, 'none' if variant == 'no_excitation' else 'settled')
    info = dict(name=name, m=int(m), bodies=int(n), fold=int(fold), quarter=int(quarter), realization=int(k),
                target=float(target), variant=variant, sample_seed=seed, eps=eps, ring_radius=R, C0=C0,
                initial_sha256=_sha(x, v), field_sha256=_sha(run.field.C, run.field.E), grid_points=int(run.field.n),
                source_harmonics=R10.source_harmonics(x, a.w, R, np.full(n, a.alpha*a.mass/n)))
    return run, info, (x, v)


def task_disturb(name, m, quarter, k, target, variant, horizon):
    t0 = time.time()
    run, info, (x, v) = build_disturbance(name, m, quarter, k, target, variant)
    out = run.advance(horizon*T0, wall_limit=R7.WALL_LIMIT).result()
    info.update(status=out['status'], t_final_periods=out['t_final_periods'], steps=out['steps'], rows=out['rows'],
                horizon=float(horizon), initial_positions=x.tolist(), initial_velocities=v.tolist(),
                seconds=round(time.time() - t0, 1))
    return info


def run_task(kind, key, args):
    t0 = time.time()
    fn = dict(A1=task_A1, A2=task_A2, box=task_box, box_control=task_box_control, norm=task_norm,
              threshold=task_threshold, disturb=task_disturb)[kind]
    return kind, key, _safe(fn(*args)), round(time.time() - t0, 1)
