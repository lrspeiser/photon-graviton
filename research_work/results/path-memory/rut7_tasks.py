"""RUT-1 stage 7 (protocol-rut7.md, amendment 1): the units of work of the campaign.

Everything that PRODUCES a raw result lives here -- the spectrum of one rung at one m, a population and its
correctness checks, one frozen or live run -- so that rut7.py, which assembles, evaluates the gates that are
computed from the archived series, and analyses part X, can be corrected without invalidating an hour of
finished runs: the campaign's resume cache is keyed on the hashes of THIS file and the modules it drives.

Every threshold is read from the machine-readable block of the protocol: this file holds no tolerance of
its own. The sampler of record at every source-sampling call site is the owner's corrected node sampler
(research_work/annulus_sampling, annulus-phase-space-v2). The historical equilibrium.py and beyn.py are
imported where a negative control needs the stage 6 behaviour and are never edited.
"""
import hashlib
import json
import math
import os
import re
import sys
import time
from dataclasses import replace
from pathlib import Path

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '1')              # one thread per worker: the pool is the parallelism

import numpy as np                              # noqa: E402
from scipy.interpolate import CubicSpline       # noqa: E402
from scipy.special import i0e                   # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
sys.path.insert(0, str(HERE))
import evidence_io          # noqa: E402
import beyn                 # noqa: E402  (stage 6's solver: only the negative controls call it)
import equilibrium as EQ    # noqa: E402  (historical, byte-identical: only a negative control calls it)
import formation as FM      # noqa: E402
import ring_modes as RM     # noqa: E402
import rut1 as U            # noqa: E402
import rut3 as R3           # noqa: E402
import spectrum as SP       # noqa: E402
import population as PO     # noqa: E402
import pairrun as PR        # noqa: E402

T0 = RM.T0
R0 = 1.
GM = FM.GM
TAU_KEEP, TAU_FORM = 10*T0, 3*T0
W = .2
LABEL = .1
CODE = ('rut7.py', 'rut7_tasks.py', 'spectrum.py', 'population.py', 'pairrun.py', 'equilibrium.py', 'beyn.py', 'ring_modes.py',
        'formation.py', 'rut3.py', 'rut1.py')
OWNER_PACKAGE = ('__init__.py', 'sampling.py', 'adapter.py')


def _find(name):
    """A committed file of this stage: beside this driver, or on the import path while it is being built."""
    for d in [HERE] + [Path(p) for p in sys.path if p]:
        if (d/name).is_file():
            return d/name
    raise FileNotFoundError(name)


def thresholds():
    """The declared thresholds, parsed out of the protocol. The code holds none of its own."""
    text = _find('protocol-rut7.md').read_text(encoding='utf-8')
    return json.loads(re.search(r'## Declared thresholds.*?```json\s*(\{.*?\})\s*```', text, re.S).group(1))


TH = thresholds()


# ---------------------------------------------------------------- small helpers
def _safe(obj, digits=None):
    """JSON-ready: numpy scalars to Python, complex to [re, im], non-finite to null (an archive holding a
    NaN never compares equal to itself), and optionally floats rounded to `digits` significant figures."""
    if isinstance(obj, dict):
        return {str(k): _safe(v, digits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_safe(v, digits) for v in obj]
    if isinstance(obj, np.ndarray):
        return _safe(obj.tolist(), digits)
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (complex, np.complexfloating)):
        return [_safe(float(obj.real), digits), _safe(float(obj.imag), digits)]
    if isinstance(obj, (float, np.floating)):
        x = float(obj)
        if not math.isfinite(x):
            return None
        return float(f'{x:.{digits}g}') if digits else x
    return obj


def _seed(*parts):
    """A declared, reproducible seed for a named draw."""
    return int(hashlib.sha256(':'.join(str(p) for p in parts).encode()).hexdigest()[:8], 16)


def _sha_arrays(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a, dtype=float).tobytes())
    return h.hexdigest()


def _match(a, b):
    """Greedy one-to-one pairing of two root lists: the worst paired distance, and what found no partner."""
    a, b = [complex(z) for z in a], [complex(z) for z in b]
    used, worst, lone = set(), 0., []
    for x in a:
        cand = [(abs(x - y)/max(1., abs(x)), i) for i, y in enumerate(b) if i not in used]
        if not cand:
            lone.append(x)
            continue
        d, i = min(cand)
        used.add(i)
        worst = max(worst, d)
    return dict(worst_relative=float(worst), unpaired_first=[[z.real, z.imag] for z in lone],
                unpaired_second=[[b[i].real, b[i].imag] for i in range(len(b)) if i not in used],
                sizes=[len(a), len(b)])


def _expanded(roots, key='s'):
    return [r[key] for r in roots for _ in range(r['multiplicity'])]


# ================================================================ PART E: the cold-ring spectrum
RUNGS = ('none', 'instantaneous', 'one_stage', 'two_stage')


def ring(model, n=32, w=W, fraction=LABEL):
    """The stage 5 configuration: 32 writers, w/R = 0.2, the 10% label."""
    return RM.Ring(model, n, R3._label_rate(w, TAU_KEEP, fraction)/n, R=R0, w=w,
                   tau_keep=TAU_KEEP, tau_form=TAU_FORM)


def _strip_summary(st):
    keys = ('counts_agree', 'any_saturated', 'smallest_gap', 'worst_residual', 'worst_polish_movement',
            'worst_winding_remainder', 'worst_winding_phase_step', 'smallest_det_on_contour_over_median')
    out = {k: st[k] for k in keys}
    out.update(rectangles=len(st['rectangles']), roots_with_multiplicity=len(_expanded(st['roots'])),
               edge_rule_applied=st['reruns'], first_pass_before_edge_rule=st['first_pass'],
               still_near_an_edge=st['still_near_an_edge'])
    return out


def _strip_detail(st):
    return [dict(rect=r['rect'], rank=r['rank'], gap=r['gap'], saturated=r['saturated'], nodes=r['nodes'],
                 beyn_count=r['count_with_multiplicity'], winding=r['winding'],
                 singular_values=r['singular_values'][:12], noise_floor=r['floor'],
                 min_det_on_contour=r['min_det_on_contour'], median_det_on_contour=r['median_det_on_contour'],
                 discarded_outside=r['discarded_outside']) for r in st['rectangles']]


def task_spectrum(model, m):
    """Every root of one rung at one m in the declared strip, from two partitions that share no cut."""
    E = TH['E']
    rg = ring(model)
    base = SP.strip(rg, m, E)
    shifted = SP.strip(rg, m, E, shifted=True)
    roots = SP.classify(base['roots'], rg.omega, m, E['unstable_floor'])
    try:
        cont = SP.continuation_check(rg, m, roots, E['strip_re'], E['strip_im'])
    except Exception as err:                                   # a cross-check only: never a gate
        cont = dict(error=repr(err))
    return dict(model=model, m=int(m), omega=float(rg.omega), roots=roots,
                base=_strip_summary(base), shifted=_strip_summary(shifted),
                partition_agreement=_match(_expanded(base['roots']), _expanded(shifted['roots'])),
                continuation_cross_check=cont,
                detail=dict(base=_strip_detail(base), shifted=_strip_detail(shifted)))


def _tiny_residue_problem(eps=1e-9):
    """T(s) with roots a and b whose residues in T^-1 differ by `eps`: the moment method's blind spot."""
    a, b, c = .3 + .2j, .6 - .1j, .7
    th = .4
    Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])

    def T(zz):
        zz = np.asarray(zz, complex)
        core = np.zeros((len(zz), 2, 2), complex)
        core[:, 0, 0], core[:, 0, 1], core[:, 1, 1] = zz - a, c, (zz - b)/eps
        return np.einsum('ab,sbc,dc->sad', Q, core, Q)
    return T, (a, b), (.1, 1., -.5, .5)


def gate_E3_control():
    """A root whose residue is 1e-9 of its neighbour's must be FOUND or FLAGGED, never silently dropped."""
    T, (a, b), rect = _tiny_residue_problem()
    got = SP.solve_generic(T, rect, TH['E'])
    found = [any(abs(r['s'] - z) < TH['E']['partition_agreement'] for r in got['roots']) for z in (a, b)]
    return dict(what='a synthetic T(s) with two roots whose residues in T^-1 differ by 1e-9',
                beyn_count=got['count_with_multiplicity'], winding=got['winding']['count'],
                small_residue_root_located=bool(found[1]), counts_agree=bool(got['counts_agree']),
                singular_values=got['singular_values'][:6],
                rejected=bool(found[1] or not got['counts_agree']))


def gate_E1():
    """Spectra known exactly. Accuracy is asked of what the contour LOCATED, before any polish."""
    E, rows = TH['E'], {}
    tol = E['raw_known_spectra']

    # (a) exp(-s) = 1/2 on a 2x2 identity: three roots, each a double eigenvalue -- six with multiplicity
    T = lambda zz: (np.exp(-np.asarray(zz, complex)) - .5)[:, None, None]*np.eye(2)[None]
    got = SP.solve_generic(T, (.1, 2., -8., 8.), E)
    exact = [np.log(2) + 2j*np.pi*k for k in (-1, 0, 1)]
    raws = [z for r in got['roots'] for z in r['raws']]
    err = max([min(abs(z - e)/max(1., abs(e)) for e in exact) for z in raws] or [1.])
    z6, dz6 = beyn.rectangle(.1, 2., -8., 8., 12, 16)
    stage6_raw, _, _ = beyn.eigenvalues(lambda s: (np.exp(-s) - .5)*np.eye(2), z6, dz6, 2, moments=12,
                                        rank_tol=1e-8)
    rows['transcendental'] = dict(
        distinct_roots=len(got['roots']), with_multiplicity=got['count_with_multiplicity'],
        winding=got['winding']['count'], worst_raw_relative_error=float(err), gap=got['gap'],
        passed=bool(len(got['roots']) == 3 and got['count_with_multiplicity'] == 6
                    and got['winding']['count'] == 6 and err < tol),
        negative_control=dict(what='the stage 6 rank rule (1e-8 of the largest singular value, unscaled '
                                   'moments) on the same problem', returned=int(len(stage6_raw)), expected=6,
                              rejected=bool(len(stage6_raw) != 6)))

    # (b) the zero-lag rung: its determinant is a quartic whose roots follow from the companion matrix
    inst = ring('instantaneous')
    worst, located, expected = 0., 0, 0
    for m in range(E['m_max'] + 1):
        st = SP.strip(inst, m, E)
        exact = [z for z in inst._quartic_roots(inst.B(m, 0.))
                 if E['strip_re'][0] < z.real < E['strip_re'][1] and abs(z.imag) < E['strip_im']]
        raws = [z for r in st['roots'] for z in r['raws']]
        located += len(raws)
        expected += len(exact)
        if raws:
            worst = max(worst, max(min(abs(z - e)/max(1., abs(e)) for e in exact) if exact else 1. for z in raws))
    rows['zero_lag_quartic'] = dict(modes=E['m_max'] + 1, located=located, expected=expected,
                                    worst_raw_relative_error=float(worst),
                                    passed=bool(located == expected and worst < tol))

    # (c) the Kepler control, s^2 (s^2 + Omega^2): all four roots LOCATED, none called unstable
    kep = ring('none')
    floor = E['unstable_floor']
    rows_k, ok = {}, True
    for m in (0, 2, 5):
        st = SP.strip(kep, m, E)
        cls = SP.classify(st['roots'], kep.omega, m, floor)
        got_roots = _expanded(st['roots'])
        pair = _match(got_roots, [0j, 0j, 1j*kep.omega, -1j*kep.omega])
        unstable = [c for c in cls if c['kind'] == 'unstable']
        good = (len(got_roots) == 4 and not pair['unpaired_first'] and not pair['unpaired_second']
                and pair['worst_relative'] <= floor and not unstable and st['counts_agree'])
        rows_k[f'm{m}'] = dict(located_with_multiplicity=len(got_roots), worst_distance=pair['worst_relative'],
                               classified_unstable=len(unstable), counts_agree=st['counts_agree'])
        ok = ok and good
    rows['kepler_control'] = dict(rows=rows_k, matched_within='the declared resolution floor, unstable_floor',
                                  passed=bool(ok))
    controls = bool(rows['transcendental']['negative_control']['rejected'])
    return dict(rows=rows, tolerance=tol, negative_control_rejected=controls,
                passed=bool(all(r['passed'] for r in rows.values()) and controls))


def gate_E2():
    """Raw robustness: the UNPOLISHED located roots under doubled quadrature and under a shifted contour."""
    E = TH['E']
    tol = E['raw_robustness']
    rows, worst = {}, 0.
    lo, hi = E['strip_re']
    raw = lambda st: [z for r in st['roots'] for z in r['raws']]
    for model, m in (('one_stage', 2), ('one_stage', 5), ('two_stage', 1), ('two_stage', 2), ('two_stage', 3)):
        rg = ring(model)
        base = raw(SP.strip(rg, m, E))
        doubled = raw(SP.strip(rg, m, E, panel=(.01, .025, .125, .025)))
        moved = raw(SP.strip(rg, m, E, shifted=True, re=(.75*lo, hi + .2)))
        # the shifted contour's left edge is nearer the axis, so it may hold fewer stable roots: compare
        # the roots both contours enclose
        keep = lambda zs: [z for z in zs if z.real > .75*lo + 1e-3]
        a, b = _match(base, doubled), _match(keep(base), keep(moved))
        rows[f'{model}:m{m}'] = dict(roots=len(base), doubled_quadrature=a['worst_relative'],
                                     shifted_contour=b['worst_relative'],
                                     unpaired=len(a['unpaired_first']) + len(a['unpaired_second'])
                                     + len(b['unpaired_first']) + len(b['unpaired_second']))
        worst = max(worst, a['worst_relative'], b['worst_relative'])
    two = ring('two_stage')
    star = max(_expanded(SP.strip(two, 2, E)['roots']), key=lambda z: z.real)
    near = []
    for panels in (12, 24):
        z, dz = beyn.rectangle(1e-3, 1.5, -6., 6., panels, 16)
        raw6, _, _ = beyn.eigenvalues(lambda s: two.M(2, s), z, dz, 2)
        near.append(raw6[int(np.argmin(np.abs(raw6 - star)))])
    control = dict(what='the stage 6 quadrature, 12 panels of 16 nodes on the stage 6 rectangle, against '
                        'its own doubling, for the two-stage m = 2 root', moved=float(abs(near[0] - near[1])),
                   raw_error_at_stage6_resolution=float(abs(near[0] - star)),
                   rejected=bool(abs(near[0] - near[1]) > tol))
    clean = all(r['unpaired'] == 0 for r in rows.values())
    return dict(rows=rows, worst=float(worst), tolerance=tol, negative_control=control,
                negative_control_rejected=control['rejected'],
                passed=bool(worst < tol and clean and control['rejected']))


def gate_E4():
    """A root placed `edge_offset` inside each edge of a sub-rectangle in turn."""
    E = TH['E']
    off = E['edge_offset']
    two = ring('two_stage')
    cuts = SP.partition(two.omega, E['strip_im'])
    star = max(_expanded(SP.strip(two, 2, E)['roots']), key=lambda z: z.real)
    k = max(i for i, c in enumerate(cuts) if c < star.imag)
    lo, hi = E['strip_re']
    placements = dict(top=(lo, hi, cuts[k], star.imag + off), bottom=(lo, hi, star.imag - off, cuts[k + 1]),
                      left=(star.real - off, hi, cuts[k], cuts[k + 1]),
                      right=(lo, star.real + off, cuts[k], cuts[k + 1]))
    rows = {}
    for name, rect in placements.items():
        got = SP.solve_rectangle(two, 2, rect, E)
        hit = [r for r in got['roots'] if abs(r['s'] - star) < E['partition_agreement']]
        rows[name] = dict(found=bool(hit), counts_agree=got['counts_agree'],
                          winding_phase_step=got['winding']['worst_phase_step'],
                          raw_error=float(abs(hit[0]['raw'] - star)) if hit else None)
    s6 = dict(top=(1e-3, 1.5, -6., star.imag + off), bottom=(1e-3, 1.5, star.imag - off, 6.),
              left=(star.real - off, 1.5, -6., 6.), right=(1e-3, star.real + off, -6., 6.))
    lost = {}
    for name, region in s6.items():
        got = beyn.solve_region(lambda s: two.M(2, s), region, residual=lambda s: two.det(2, s),
                                polish=lambda s: two.det(2, s))
        lost[name] = not any(abs(r['s'] - star) < E['partition_agreement'] for r in got['roots'])
    control = dict(what='the stage 6 solver, which filters on the region and only then polishes, with the same '
                        'root the same distance inside each edge of its own rectangle',
                   lost_at=sorted(k for k, v in lost.items() if v), rejected=bool(any(lost.values())))
    return dict(root=[star.real, star.imag], offset=off, rows=rows, negative_control=control,
                negative_control_rejected=control['rejected'],
                passed=bool(all(r['found'] and r['counts_agree'] for r in rows.values()) and control['rejected']))


# ================================================================ PART S: the population
L0, DL = 1., .06
POPULATIONS = (
    dict(name='A_cold', family='A', dE=.010, alpha_from=None),
    dict(name='A_warm', family='A', dE=.030, alpha_from=None),
    dict(name='B_mid', family='B', dE=.020, alpha_from=None),          # calibrates the family's alpha
    dict(name='B_cold', family='B', dE=.010, alpha_from='B_mid'),
    dict(name='B_warm', family='B', dE=.030, alpha_from='B_mid'),
)
GRID = dict(r_lo=.3, r_hi=3., n_r=280, n_L=160, n_vr=96)
REFINEMENTS = dict(
    radial_grid=dict(n_r=559),                                   # the same nodes, and one between each pair
    L_quadrature=dict(n_L=320),
    vr_quadrature=dict(n_vr=192),
    radial_domain=dict(r_lo=.2, r_hi=3.3),                       # n_r is rescaled to keep the spacing
    all_together=dict(n_r=559, n_L=320, n_vr=192, r_lo=.2, r_hi=3.3),
)


def _grid(**change):
    g = dict(GRID, **change)
    if ('r_lo' in change or 'r_hi' in change):
        per_unit = (change.get('n_r', GRID['n_r']) - 1)/(GRID['r_hi'] - GRID['r_lo'])
        g['n_r'] = int(round(per_unit*(g['r_hi'] - g['r_lo']))) + 1
    return g


def build(spec, alpha=None, **change):
    """One population on one discretization: alpha solved for the declared mean support, or imposed."""
    a = PO.Annulus(L0=L0, dL=DL, dE=spec['dE'], w=W, tau_keep=TAU_KEEP, mass=1., reach=TH['S']['reach'],
                   **_grid(**change))
    return a.solve_mean_support(TH['S']['mean_support']) if alpha is None else a.at_alpha(alpha)


def _peak(r, y):
    """The maximum of a sampled profile, from the parabola through its largest node and that node's
    neighbours. The largest NODE VALUE is not the peak: it carries a sampling error of order y'' dr^2/8,
    about 3e-4 of the peak here, which is the mistake this stage records against stage 6's label."""
    i = int(np.clip(np.argmax(y), 1, len(y) - 2))
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    curv = y0 - 2*y1 + y2
    shift = .5*(y0 - y2)/curv if curv < 0 else 0.
    return float(y1 - .25*(y0 - y2)*shift), float(r[i] + shift*(r[1] - r[0]))


def describe(a):
    """What a population is, in the quantities the refinement gate follows."""
    wgt = 2*np.pi*a.r*a.sigma
    mass = float(np.trapezoid(wgt, a.r))
    mean_r = float(np.trapezoid(wgt*a.r, a.r)/mass)
    width = float(np.sqrt(np.trapezoid(wgt*(a.r - mean_r)**2, a.r)/mass))
    sup = a.support_profile()
    disp = a.dispersions(R=mean_r, n_L=a.n_L, n_vr=a.n_vr, reach=a.reach)
    return dict(alpha=float(a.alpha), mass=mass, total_writing_alpha_times_mass=float(a.alpha*mass),
                mean_radius=mean_r, radial_width=width, peak_field=_peak(a.r, a.C)[0],
                peak_field_radius=_peak(a.r, a.C)[1], largest_field_node=float(a.C.max()),
                peak_density=_peak(a.r, a.sigma)[0], support=sup, support_solve_error=float(a.support_error),
                fixed_point_residual=a.consistency_residual(),
                dispersions_at_mean_radius=dict(sigma_r=disp['sigma_r'], sigma_t=disp['sigma_t'],
                                                mean_azimuthal_speed=disp['v_mean']),
                grid=dict(r_lo=float(a.r[0]), r_hi=float(a.r[-1]), n_r=len(a.r), n_L=a.n_L, n_vr=a.n_vr,
                          reach=a.reach))


def task_population(spec, alpha=None):
    a = build(spec, alpha)
    nodes = a.nodes()
    return dict(name=spec['name'], family=spec['family'], dE=spec['dE'], dL=DL, L0=L0,
                alpha_source='solved for the declared mean support' if alpha is None
                else f"held at the value calibrated on {spec['alpha_from']}",
                description=describe(a), moments=a.moments(nodes), node_count=int(len(nodes.node_mass)),
                state=a.state())


def _annulus(pop):
    return PO.Annulus.from_state(pop['state'])


def _phi_spline(a):
    return CubicSpline(a.r, a.phi_total(a.r, a.C))


def _body_quantities(a, x, v, phi=None):
    phi = _phi_spline(a) if phi is None else phi
    r = np.linalg.norm(x, axis=1)
    vr = np.sum(v*x, axis=1)/r
    L = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
    return dict(R=r, L=L, E=.5*np.sum(v*v, axis=1) + phi(r), vr=vr)


# ---------------------------------------------------------------- V1: the measure
def task_V1(pop):
    S = TH['S']
    a = _annulus(pop)
    before = _sha_arrays(a.r, a.sigma, a.C)
    nodes = a.nodes()
    ok = PO.owner_sampler().verify_distribution(a, nodes, n_L=a.n_L, n_vr=a.n_vr, reach=a.reach,
                                                tolerance=S['measure_identity'])
    unmutated = before == _sha_arrays(a.r, a.sigma, a.C)
    biased = replace(nodes, node_mass=nodes.node_mass*nodes.node_radius)      # the stage 6 weights
    bad = PO.owner_sampler().verify_distribution(a, biased, n_L=a.n_L, n_vr=a.n_vr, reach=a.reach,
                                                 tolerance=S['measure_identity'])
    mo = nodes.moments()['radius']
    shift = biased.moments()['radius']['mean'] - mo['mean']
    predicted = mo['variance']/mo['mean']
    control = dict(what='the stage 6 node weights, with their extra factor of R',
                   verification_passed=bool(bad['numerical_verification_passed']),
                   mean_radius_shift=float(shift), variance_over_mean=float(predicted),
                   shift_identity_error=float(abs(shift - predicted)),
                   rejected=bool(not bad['numerical_verification_passed']
                                 and abs(shift - predicted) < S['measure_identity']))
    keep = ('numerical_verification_passed', 'sampler_version', 'radial_marginal_max_absolute_error',
            'mass_relative_error', 'radius_moment_relative_errors', 'node_count')
    return dict(name=pop['name'], owner_verification={k: ok[k] for k in keep}, source_unmutated=bool(unmutated),
                negative_control=control,
                passed=bool(ok['numerical_verification_passed'] and unmutated and control['rejected']))


# ---------------------------------------------------------------- V2: drawn moments, two samplers
def _moment_test(q, exact, n, S):
    """Means within `mean_sigma` standard errors of quadrature; variances within `variance_relative`."""
    out, ok = {}, True
    for key in ('R', 'L', 'E'):
        se = math.sqrt(exact[f'{key}_var']/n)
        z = (float(q[key].mean()) - exact[f'{key}_mean'])/se
        out[f'{key}_mean_z'] = z
        ok = ok and abs(z) <= S['mean_sigma']
    for key in ('R', 'L', 'E'):
        ratio = float(q[key].var())/exact[f'{key}_var']
        out[f'{key}_variance_ratio'] = ratio
        ok = ok and abs(ratio - 1) <= S['variance_relative']
    ratio = float(np.mean(q['vr']**2))/exact['vr_var']
    out['vr_variance_ratio'] = ratio
    return out, bool(ok and abs(ratio - 1) <= S['variance_relative'])


def task_V2(pop):
    S = TH['S']
    a = _annulus(pop)
    n = int(S['draws'])
    exact, phi = pop['moments'], _phi_spline(a)
    samples = dict(owner_node_sampler=a.draw(n, seed=_seed('V2', pop['name'], 'owner')),
                   continuous_rejection=a.draw_continuous(n, seed=_seed('V2', pop['name'], 'continuous')))
    rows, ok, q = {}, True, {}
    for name, (x, v) in samples.items():
        q[name] = _body_quantities(a, x, v, phi)
        rows[name], good = _moment_test(q[name], exact, n, S)
        rows[name]['distinct_radii'] = int(len(np.unique(np.round(q[name]['R'], 12))))
        rows[name]['passed'] = good
        ok = ok and good
    rows['continuous_rejection']['acceptance'] = float(a.acceptance)
    between = {}
    for key in ('R', 'L', 'E'):
        se = math.sqrt(2*exact[f'{key}_var']/n)
        between[f'{key}_mean_z'] = float(q['owner_node_sampler'][key].mean()
                                         - q['continuous_rejection'][key].mean())/se
    agree = all(abs(z) <= S['mean_sigma'] for z in between.values())
    x, v = samples['owner_node_sampler']
    m = n//4                                   # a quarter, so drawing without replacement stays unbiased
    xb, vb = PO.reweighted_by_radius(x, v, m, seed=_seed('V2', pop['name'], 'control'))
    bad, bad_ok = _moment_test(_body_quantities(a, xb, vb, phi), exact, m, S)
    control = dict(what='a correct draw re-weighted by R, which is what stage 6 drew', bodies=m, moments=bad,
                   rejected=bool(not bad_ok))
    return dict(name=pop['name'], draws=n, samplers=rows, between_samplers=between,
                samplers_agree=bool(agree), negative_control=control,
                passed=bool(ok and agree and control['rejected']))


# ---------------------------------------------------------------- V3: the declared support fits
def _support_fits(ext, a):
    S = TH['S']
    box = S['box_fraction']*PR.HALF_WIDTH
    return bool(ext['all_bound'] and ext['pericentre_min'] > a.r[0] and ext['apocentre_max'] < a.r[-1]
                and ext['apocentre_max'] < box and ext['circular_radius_increasing']
                and ext['circular_radius_continuous'])


def _extent(a, reach):
    """The most extended orbit of a support of `reach` widths, and whether E_circ's minimising radius runs
    continuously in L. A continuous function's largest step halves when the spacing in L halves; a switch
    of the global minimum from one well to another keeps its size. The rule used is that the largest step
    on 321 points is under three quarters of the largest step on 161."""
    ext = a.orbit_extent(reach)
    jump = {}
    for n in (161, 321):
        L = np.linspace(a.L0 - reach*a.dL, a.L0 + reach*a.dL, n)
        rc = a.circular_energy(L, a.r, a.C, a.sigma)[1]
        jump[n] = float(np.max(np.abs(np.diff(rc))))
    ext['circular_radius_largest_step_161_points'] = jump[161]
    ext['circular_radius_largest_step_321_points'] = jump[321]
    ext['circular_radius_continuous'] = bool(jump[321] < .75*jump[161])
    ext['simulator_limit'] = float(TH['S']['box_fraction']*PR.HALF_WIDTH)
    return ext


def task_V3(pop):
    a = _annulus(pop)
    ext = _extent(a, a.reach)
    wide = _extent(a, 5.)
    return dict(name=pop['name'], declared_support=ext, fits=_support_fits(ext, a),
                five_width_support=wide, five_width_fits=_support_fits(wide, a))


# ---------------------------------------------------------------- V4: discretization
def task_V4(spec, variant, alpha=None):
    """One population on one refined discretization. `alpha` is imposed for the fixed-physics populations."""
    t0 = time.time()
    a = build(spec, alpha, **REFINEMENTS[variant]) if variant != 'base' else build(spec, alpha)
    d = describe(a)
    return dict(name=spec['name'], variant=variant, description=d, r=a.r.tolist(), sigma=a.sigma.tolist(),
                C=a.C.tolist(), seconds=round(time.time() - t0, 1))


def task_V4_control(dE, n_r):
    """Stage 6's label: the support read at the argmax grid node, solved with its unbracketed secant."""
    t0 = time.time()
    a = EQ.WarmAnnulus(L0=L0, dL=DL, dE=dE, w=W, tau_keep=TAU_KEEP, n_r=n_r)
    a.solve(target_support=TH['S']['mean_support'])
    return dict(dE=dE, n_r=n_r, alpha=float(a.alpha), seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- V5: the drawn source writes the field
def ring_averaged_field(a, radii):
    """Azimuthal average of the steady field written by equal-mass bodies at the given radii."""
    out = np.zeros_like(a.r)
    for chunk in np.array_split(np.asarray(radii, float), max(1, len(radii)//2048)):
        out += (np.exp(-(a.r[:, None] - chunk[None, :])**2/(2*a.w**2))
                * i0e(a.r[:, None]*chunk[None, :]/a.w**2)).sum(1)
    return a.tau_keep*a.alpha*a.mass/len(radii)*out


def field_error_cell(a, n, k, name, nodes, biased=False):
    seed = _seed('V5', name, n, k, 'control' if biased else 'draw')
    if biased:
        x, v = nodes.draw(4*n, seed)
        x, v = PO.reweighted_by_radius(x, v, n, seed=seed)
    else:
        x, v = nodes.draw(n, seed)
    err = ring_averaged_field(a, np.linalg.norm(x, axis=1)) - a.C
    return float(np.sqrt(np.mean(err**2))/a.C.max()), float(np.max(np.abs(err))/a.C.max())


def task_V5(pop, counts=(256, 1024, 4096, 16384)):
    S = TH['S']
    a = _annulus(pop)
    nodes = a.nodes()
    reps = int(S['field_realizations'])
    rows = {}
    for n in counts:
        cells = [field_error_cell(a, n, k, pop['name'], nodes) for k in range(reps)]
        rms = np.array([c[0] for c in cells])
        rows[str(n)] = dict(rms_error_over_peak=float(np.sqrt(np.mean(rms**2))), realizations=reps,
                            scatter=float(rms.std()), max_error_over_peak=float(np.mean([c[1] for c in cells])))
    lo, hi = S['field_exponent']
    ns = np.log(np.array(counts, float))
    rms_of = lambda table: np.array([table[str(n)]['rms_error_over_peak'] for n in counts])
    exponent = float(np.polyfit(ns, np.log(rms_of(rows)), 1)[0])
    last = str(counts[-1])
    good = bool(lo <= exponent <= hi and rows[last]['rms_error_over_peak'] < S['field_error_16384'])
    # the negative control goes through the SAME gate: the same body counts, realizations, fit and bound
    broken = {}
    for n in counts:
        rms = np.array([field_error_cell(a, n, k, pop['name'], nodes, biased=True)[0] for k in range(reps)])
        broken[str(n)] = dict(rms_error_over_peak=float(np.sqrt(np.mean(rms**2))))
    bad_exponent = float(np.polyfit(ns, np.log(rms_of(broken)), 1)[0])
    floor = broken[last]['rms_error_over_peak']
    control = dict(what='the R-re-weighted sample through the same gate', rows=broken,
                   fitted_exponent=bad_exponent, error_at_largest_count=floor,
                   floor_above_the_bound=bool(floor > S['field_error_16384']),
                   exponent_outside_the_range=bool(not lo <= bad_exponent <= hi),
                   rejected=bool(not (lo <= bad_exponent <= hi and floor < S['field_error_16384'])))
    return dict(name=pop['name'], rows=rows, fitted_exponent=exponent, negative_control=control,
                error_norm='rms over the solver radial grid of the difference from C0, over the peak of C0',
                requirement_met=good, passed=bool(good and control['rejected']))


# ---------------------------------------------------------------- V6: the draw is stationary
V6_CONTROLS = ('kicked_outward', 'velocities_scaled', 'reweighted_by_radius')


def stationary_bodies(a, name, k, control=None, n=None, nodes=None):
    """The declared draw for realization k, or one of the three broken versions of it."""
    S = TH['S']
    n = int(S['stationary_bodies']) if n is None else n
    nodes = a.nodes() if nodes is None else nodes
    seed = _seed('V6', name, k)
    if control == 'reweighted_by_radius':
        x, v = nodes.draw(4*n, seed)
        return PO.reweighted_by_radius(x, v, n, seed=seed)
    x, v = nodes.draw(n, seed)
    if control == 'kicked_outward':
        v = v + S['kick']*x/np.linalg.norm(x, axis=1, keepdims=True)
    elif control == 'velocities_scaled':
        v = v*S['velocity_scale']
    return x, v


def stationarity_statistics(rows, exact, n, S):
    """The ensemble at t = 0 against its own average over 5 T0 to the horizon, in sampling standard errors,
    and the rms fluctuation of the instantaneous mean radius about that average."""
    t = np.array([r['t_periods'] for r in rows])
    late = t >= 5.
    var_R, var_vr = exact['R_var'], exact['vr_var']
    se = dict(mean_radius=math.sqrt(var_R/n),
              radial_spread=math.sqrt(max(exact['R_fourth'] - var_R**2, 0.)/(4*var_R*n)),
              vr_rms=math.sqrt(max(exact['vr_fourth'] - var_vr**2, 0.)/(4*var_vr*n)),
              vt_mean=math.sqrt(exact['vt_var']/n))
    z = {}
    for key in se:
        y = np.array([r[key] for r in rows])
        z[key] = float((y[0] - y[late].mean())/se[key])
    mr = np.array([r['mean_radius'] for r in rows])[late]
    fluct = float(np.sqrt(np.mean((mr - mr.mean())**2))/se['mean_radius'])
    ok = all(abs(v) <= S['stationary_sigma'] for v in z.values()) and fluct <= S['stationary_fluctuation_sigma']
    return dict(z_initial_against_time_average=z, mean_radius_fluctuation_in_standard_errors=fluct,
                standard_errors=se, stationary=bool(ok))


def task_V6(pop, k, control=None, horizon=None, n=None):
    S = TH['S']
    t0 = time.time()
    a = _annulus(pop)
    x, v = stationary_bodies(a, pop['name'], k, control, n)
    n = len(x)
    horizon = S['stationary_horizon'] if horizon is None else horizon
    run = PR.Run(a, x, v, np.full(n, a.alpha*a.mass/n), live=False, tau_form=TAU_FORM, light=True,
                 record_every=TH['X']['record_every']*T0)
    out = run.advance(horizon*T0, wall_limit=3600.).result()
    whole = out['status'] == 'completed' and out['t_final_periods'] > 5.      # a replayed prefix has no late window
    stats = stationarity_statistics(out['rows'], pop['moments'], n, S) if whole else None
    return dict(name=pop['name'], realization=int(k), control=control, bodies=n, status=out['status'],
                t_final_periods=out['t_final_periods'], steps=out['steps'], statistics=stats,
                rows=out['rows'], seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- V7: restart
def task_V7(pop, n=64, half=2., full=4.):
    a = _annulus(pop)
    x, v = a.draw(n, seed=_seed('V7', pop['name']))
    rates = np.full(n, a.alpha*a.mass/n)
    kw = dict(live=True, tau_form=TAU_FORM, record_every=TH['X']['record_every']*T0)
    whole = PR.Run(a, x, v, rates, **kw).advance(full*T0)
    first = PR.Run(a, x, v, rates, **kw).advance(half*T0)
    snap = first.snapshot()
    resumed = PR.Run(a, x, v, rates, **kw).restore(snap).advance(full*T0)
    bare = PR.Run(a, x, v, rates, **kw).restore(snap, fields=False).advance(full*T0)
    same = dict(positions=bool(np.array_equal(whole.x, resumed.x)), velocities=bool(np.array_equal(whole.v, resumed.v)),
                C=bool(np.array_equal(whole.field.C, resumed.field.C)),
                E=bool(np.array_equal(whole.field.E, resumed.field.E)),
                records=bool(json.dumps(whole.rows) == json.dumps(resumed.rows)),
                steps=bool(whole.steps == resumed.steps))
    gap = float(np.max(np.abs(whole.x - bare.x)))
    control = dict(what='a restart from the bodies alone, the field re-primed and its history discarded',
                   max_position_difference=gap, rejected=bool(gap > 0.))
    return dict(population=pop['name'], bodies=n, saved_at_periods=half, horizon_periods=full,
                bit_for_bit=same, negative_control=control,
                passed=bool(all(same.values()) and control['rejected']))


# ================================================================ PART R: the refinement stage 5 quoted
R_CASES = (dict(name='w5_h0.01', spacing_divisor=5, step=.01, primed=True),
           dict(name='w10_h0.01', spacing_divisor=10, step=.01, primed=True),
           dict(name='w5_h0.005', spacing_divisor=5, step=.005, primed=True),
           dict(name='control_no_primed_field', spacing_divisor=5, step=.01, primed=False))


def task_R(case, horizon=None):
    """The unseeded two-stage cold ring: the m = 2 mode emerging alone, measured against linear theory."""
    Rt = TH['R']
    t0 = time.time()
    horizon = Rt['horizon'] if horizon is None else horizon
    n = 32
    q = R3._label_rate(W, TAU_KEEP, LABEL)
    extra = U.ring_inward_acceleration(R0, R0, W, q*TAU_KEEP)
    pos, vel, rates = FM.ring_start(n, R0, q, 0., 0., 1, extra_inward=extra)
    out = RM.run_ring(pos, vel, rates, W, TAU_KEEP, horizon*T0, case['step'], tau_form=TAU_FORM,
                      eta=case['step'], spacing=W/case['spacing_divisor'],
                      prime_ring_R=R0 if case['primed'] else None)
    c2 = np.mean((out['r'] - R0)*np.exp(-2j*np.pi*2*np.arange(n)[None, :]/n), axis=1)
    return dict(case=case, status=out['status'], t_final_periods=float(out['t_final']/T0),
                t_periods=(out['t']/T0).tolist(), c2=[[z.real, z.imag] for z in c2],
                seconds=round(time.time() - t0, 1))


# ================================================================ PART X: live against frozen
WALL_LIMIT = 6*3600.        # a run that grinds past this is stopped and says so; never a silent hang


def x_seed(name, n, k):
    return _seed('X', name, n, k)


def task_X(pop, n, k, live, horizon=None):
    """One half of a pair. Both halves draw the same bodies from the same declared seed."""
    X = TH['X']
    t0 = time.time()
    a = _annulus(pop)
    seed = x_seed(pop['name'], n, k)
    x, v = a.draw(n, seed)
    run = PR.Run(a, x, v, np.full(n, a.alpha*a.mass/n), live=live, tau_form=TAU_FORM,
                 record_every=X['record_every']*T0)
    out = run.advance((X['horizon'] if horizon is None else horizon)*T0, wall_limit=WALL_LIMIT).result()
    rows = out['rows']
    if not live:                        # a frozen field's coefficients never change: keep the source only
        for r in rows:
            r['coefficients'] = dict(S=r['coefficients']['S'])
    return dict(name=pop['name'], bodies=int(n), realization=int(k), live=bool(live), seed=seed,
                initial_sha256=_sha_arrays(x, v), status=out['status'],
                t_final_periods=out['t_final_periods'], steps=out['steps'], ring_radius=out['ring_radius'],
                rows=rows, seconds=round(time.time() - t0, 1))


# ================================================================ the campaign
def run_task(kind, key, args):
    """One unit of work in a worker process."""
    t0 = time.time()
    fn = dict(spectrum=task_spectrum, E1=gate_E1, E2=gate_E2, E3_control=gate_E3_control, E4=gate_E4,
              V1=task_V1, V2=task_V2, V3=task_V3, V4=task_V4, V4_control=task_V4_control, V5=task_V5,
              V6=task_V6, V7=task_V7, R=task_R, X=task_X, population=task_population)[kind]
    out = _safe(fn(*args))
    return kind, key, out, round(time.time() - t0, 1)
