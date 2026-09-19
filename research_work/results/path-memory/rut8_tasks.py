"""RUT-1 stage 8 (protocol-rut8.md): the units of work of the campaign.

Everything that PRODUCES a raw result lives here, as in stage 7, so that rut8.py -- assembly, the gates that
are evaluated on finished tasks, the predictions and the readings -- can be corrected without invalidating
finished work: the resume cache is keyed on the hashes of this file and the modules it drives.

Every threshold is read from the machine-readable block of the protocol: this file holds no tolerance of its
own. Populations are built, checked, drawn and evolved with STAGE 7'S CODE BY IMPORT, unchanged; the owner's
corrected node sampler remains the sampler of record.
"""
import gzip
import json
import math
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut7_tasks as R7        # noqa: E402  (sets one thread per worker before numpy is imported)
from rut7_tasks import PO, RM, SP, T0, TAU_FORM, TAU_KEEP, W, _find, _safe, _seed, np  # noqa: E402
from scipy.interpolate import CubicSpline   # noqa: E402
import warm_modes as WM        # noqa: E402

CODE = ('rut8.py', 'rut8_tasks.py', 'warm_modes.py', 'rut7_tasks.py', 'rut7.py', 'spectrum.py', 'population.py',
        'pairrun.py', 'equilibrium.py', 'ring_modes.py', 'formation.py', 'rut3.py', 'rut1.py')
STAGE7 = ('A_cold', 'A_warm', 'B_cold', 'B_mid', 'B_warm')


def thresholds():
    text = _find('protocol-rut8.md').read_text(encoding='utf-8')
    return json.loads(re.search(r'## Declared thresholds.*?```json\s*(\{.*?\})\s*```', text, re.S).group(1))


TH = thresholds()
QUAD = TH['L']['contour']


# ---------------------------------------------------------------- populations
def stage7_populations():
    """The five populations of stage 7, from its committed series: the state the simulations used."""
    with gzip.open(_find('rut7-results.json').parent/'rut7-series'/'populations.json.gz', 'rt', encoding='utf-8') as f:
        return json.load(f)


def alpha_B():
    return stage7_populations()['B_mid']['description']['alpha']


def task_population(name, dE):
    """A member of family B: stage 7's construction, with the family's coupling imposed."""
    return R7.task_population(dict(name=name, family='B', dE=float(dE), alpha_from='B_mid'), alpha_B())


# ---------------------------------------------------------------- T(s), cached within a worker
_ORBITS, _MODES = {}, {}


def modes_for(pop, m, harmonics=None, **override):
    L = dict(TH['L'], **override)
    okey = (pop['name'], pop['state']['alpha'], pop['state']['dE'], L['n_L'], L['n_u'], L['n_eta'])
    if okey not in _ORBITS:
        _ORBITS.clear()
        _MODES.clear()
        _ORBITS[okey] = WM.Orbits(WM.Disk(pop['state'], TAU_FORM), L['n_L'], L['n_u'], L['n_eta'])
    mkey = okey + (m, L['l_max'], L['node_spacing'], L['node_margin'], None if harmonics is None else tuple(harmonics))
    if mkey not in _MODES:
        if len(_MODES) >= 2:
            _MODES.clear()
        orb = _ORBITS[okey]
        _MODES[mkey] = WM.Modes(orb.disk, orb, m, L['l_max'], L['node_spacing'], L['node_margin'], L['node_clip'],
                                L['kernel_rank_cutoff'], harmonics=harmonics)
    return _MODES[mkey]


def _rectangle_row(got):
    rows = [dict(re=r['s'].real, im=r['s'].imag, raw=[r['raw'].real, r['raw'].imag], moved=r['moved'],
                 relative_residual=r['relative_residual']) for r in got['roots']]
    return dict(rect=got['rect'], roots=rows, located=got['located'], rank=got['rank'], gap=got['gap'],
                clear_rank=got['clear_rank'], winding=got['winding'], counts_agree=got['counts_agree'],
                saturated=got['saturated'], nodes=got['nodes'], singular_values=got['singular_values'])


def task_rectangle(pop, m, rect, shifted):
    """One sub-rectangle of one declared strip: located by moments, counted by the argument principle."""
    modes = modes_for(pop, m)
    row = _rectangle_row(WM.solve_rectangle(modes, tuple(rect), TH['L'], QUAD))
    row.update(name=pop['name'], m=int(m), shifted=bool(shifted), kernel_nodes=int(len(modes.nodes)),
               basis_size=modes.basis_size, kernel_condition=modes.condition, rows_kept=modes.rows_kept,
               rows_total=modes.rows_total)
    return row


def task_strip(pop, m, shifted):
    """Every sub-rectangle of one declared strip, in one worker: the orbit library and the transforms are
    built once."""
    return dict(name=pop['name'], m=int(m), shifted=bool(shifted),
                rectangles=[task_rectangle(pop, m, rect, shifted) for rect in WM.strip_rectangles(m, TH['L'], shifted)])


def task_family(dE):
    """One member of family B: built, and its m = 2 root located in the declared rectangle."""
    pop = task_population(f'B_dE{dE:.3f}', dE)
    modes = modes_for(pop, 2)
    row = _rectangle_row(WM.solve_rectangle(modes, tuple(TH['P']['family_rect']), TH['L'], QUAD))
    d = pop['description']
    row.update(dE=float(dE), alpha=d['alpha'], mean_support=d['support']['mean'],
               sigma_r=d['dispersions_at_mean_radius']['sigma_r'], mean_radius=d['mean_radius'])
    return row


# ---------------------------------------------------------------- L1: the orbit library
def gate_L1(pop, field=True):
    L = TH['L']
    try:
        disk = WM.Disk(pop['state'], TAU_FORM, field=field)
        orb = WM.Orbits(disk, L['n_L'], L['n_u'], L['n_eta'])
    except RuntimeError as err:                       # the control's orbits may not even fit the domain
        return dict(name=pop['name'], field=field, error=str(err), passed=False)
    r, sigma = disk.r, disk.sigma
    wgt = 2*np.pi*r*sigma
    mean_r = float(np.trapezoid(wgt*r, r)/np.trapezoid(wgt, r))
    tests = dict(mass=lambda x: np.ones_like(x), R=lambda x: x, R2=lambda x: x*x, inverse_R=lambda x: 1./x,
                 window=lambda x: np.exp(-(x - mean_r)**2/(2*disk.w**2)))
    moments = {}
    for key, g in tests.items():
        direct = float(np.trapezoid(wgt*g(r), r))
        moments[key] = orb.moment(g)/direct - 1.
    # the most nearly circular orbit of every L against the epicyclic limit of the same potential
    first = np.array([np.nonzero(orb.L == Li)[0][np.argmin(orb.u[orb.L == Li])] for Li in orb.circular_rows[:, 0]])
    kappa, omega_c = orb.circular_rows[:, 4], orb.circular_rows[:, 3]
    epicyclic = dict(radial=float(np.nanmax(np.abs(orb.Om_r[first]/kappa - 1))),
                     azimuthal=float(np.max(np.abs(orb.Om_th[first]/omega_c - 1))))
    worst = float(max(abs(v) for v in moments.values()))
    return dict(name=pop['name'], field=field, orbits=int(orb.n), narrow=int(orb.narrow), two_wells=int(orb.two_wells),
                moments_relative=moments, worst_moment=worst, epicyclic=epicyclic,
                frequency_ranges=dict(radial=[float(orb.Om_r.min()), float(orb.Om_r.max())],
                                      azimuthal=[float(orb.Om_th.min()), float(orb.Om_th.max())]),
                passed=bool(worst < L['moments'] and max(epicyclic.values()) < L['epicyclic']))


# ---------------------------------------------------------------- L2: the response matrix, in the time domain
def task_L2(pop, s_pair):
    """Test particles in the frozen potential plus an imposed growing perturbation, +eps and -eps on identical
    bodies so that shot noise cancels at leading order. At t = 0 the kernel-weighted m-th harmonic of the
    difference, over eps, must be -M_ab(s), for the column b whose kernel sits nearest the mean radius."""
    L = TH['L']
    t0 = time.time()
    m, s = 2, complex(*s_pair)
    modes = modes_for(pop, m)
    disk, nodes = modes.d, modes.nodes
    b = int(np.argmin(np.abs(nodes - pop['description']['mean_radius'])))
    a = PO.Annulus.from_state(pop['state'])
    n = int(L['time_domain_bodies'])
    x, v = a.draw(n, seed=_seed('L2', pop['name']))
    rg = np.linspace(disk.r[0], disk.r[-1], 4000)
    psi = CubicSpline(rg, modes.kernel(nodes[b], rg))
    dpsi = psi.derivative()
    gamma, beta, eps, h = s.real, s.imag, L['time_domain_eps'], L['time_domain_step']
    steps = int(round(L['time_domain_efoldings']/gamma/h))
    q = np.concatenate([x, x])
    p = np.concatenate([v, v])
    sign = np.concatenate([np.ones(n), -np.ones(n)])

    def accel(q, t):
        r = np.hypot(q[:, 0], q[:, 1])
        th = np.arctan2(q[:, 1], q[:, 0])
        ur = q/r[:, None]
        amp = sign*eps*math.exp(gamma*t)
        phase = m*th + beta*t
        fr = -disk.dphi(r) + amp*dpsi(r)*np.cos(phase)
        ft = -amp*psi(r)*m*np.sin(phase)/r
        return fr[:, None]*ur + ft[:, None]*np.stack([-ur[:, 1], ur[:, 0]], axis=1)

    t = -steps*h
    acc = accel(q, t)
    for _ in range(steps):
        p_half = p + .5*h*acc
        q = q + h*p_half
        t += h
        acc = accel(q, t)
        p = p_half + .5*h*acc
    if not np.isfinite(q).all():
        raise FloatingPointError('non-finite test-particle state')
    r = np.hypot(q[:, 0], q[:, 1])
    ph = np.exp(-1j*m*np.arctan2(q[:, 1], q[:, 0]))
    proj = np.array([(modes.kernel(ra, r)*ph*sign).sum() for ra in nodes])*a.mass/n
    measured = proj/eps                 # [<.>_+ - <.>_-]/(2 eps), times 2: the real perturbation is half e^{+i m theta}
    predicted = -modes.M(s)[:, b]
    only_l0 = -modes_for(pop, m, harmonics=[0]).M(s)[:, b]
    scale = float(np.max(np.abs(predicted)))
    worst = float(np.max(np.abs(measured - predicted))/scale)
    control = float(np.max(np.abs(measured - only_l0))/scale)
    tol = L['time_domain_relative']
    return dict(name=pop['name'], s=[s.real, s.imag], bodies=n, steps=steps, column_node=float(nodes[b]),
                kernel_nodes=[float(z) for z in nodes], measured=[[z.real, z.imag] for z in measured],
                predicted=[[z.real, z.imag] for z in predicted], worst_relative=worst,
                negative_control=dict(what='the prediction with only the l = 0 harmonic kept',
                                      worst_relative=control, rejected=bool(control > tol)),
                passed=bool(worst < tol and control > tol), seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- L3: the ring limit
def task_L3(config):
    """A narrow cold annulus against the exact-history RING of stage 5, evaluated at the annulus's own mean
    radius and total writing rate with stage 7's certified contour solver."""
    dL, dE, n_r = config
    t0 = time.time()
    a = PO.Annulus(L0=1., dL=dL, dE=dE, w=W, tau_keep=TAU_KEEP, mass=1., reach=R7.TH['S']['reach'],
                   r_lo=.5, r_hi=1.6, n_r=int(n_r))
    a.solve_mean_support(TH['L']['ring_limit_support'])
    wgt = 2*np.pi*a.r*a.sigma
    mean_r = float(np.trapezoid(wgt*a.r, a.r))
    width = float(np.sqrt(np.trapezoid(wgt*(a.r - mean_r)**2, a.r)))
    pop = dict(name=f'ring_limit_dL{dL}', state=a.state())
    got = WM.solve_rectangle(modes_for(pop, 2), tuple(TH['L']['ring_limit_rect']), TH['L'], QUAD)
    annulus = max((r['s'] for r in got['roots']), key=lambda z: z.real, default=None)

    def ring_root(R):
        rg = RM.Ring('two_stage', 64, a.alpha*a.mass/64, R=R, w=W, tau_keep=TAU_KEEP, tau_form=TAU_FORM)
        roots = [r['s'] for r in SP.strip(rg, 2, R7.TH['E'])['roots']]
        s = max(roots, key=lambda z: z.real)
        return complex(s.real, s.imag - 2*rg.omega), float(rg.omega)        # rotating frame -> inertial

    ring, omega = ring_root(mean_r)
    wrong, _ = ring_root(1.)
    rel = lambda target: None if annulus is None else float(abs(annulus.real/target.real - 1))
    return dict(dL=dL, dE=dE, n_r=int(n_r), alpha=float(a.alpha), mean_radius=mean_r, radial_width=width,
                mean_support=float(a.support_profile()['mean']), counts_agree=got['counts_agree'],
                annulus_root=None if annulus is None else [annulus.real, annulus.imag],
                ring_root_inertial=[ring.real, ring.imag], ring_omega=omega, growth_relative_difference=rel(ring),
                ring_at_stage5_radius=[wrong.real, wrong.imag], growth_relative_difference_wrong_radius=rel(wrong),
                seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- L4: convergence of a root
L4_VARIANTS = dict(orbit_quadrature=dict(n_L=2, n_u=2), orbit_angle=dict(n_eta=2), harmonics=dict(l_max=12),
                   kernel_nodes=dict(node_spacing=.075), kernel_range=dict(node_margin=.3))


def task_L4(pop, variant, rect):
    L = TH['L']
    if variant == 'base':
        modes = modes_for(pop, 2)
    elif variant == 'control_l0_only':
        modes = modes_for(pop, 2, harmonics=[0])
    else:
        change = dict(L4_VARIANTS[variant])
        for key in ('n_L', 'n_u', 'n_eta'):
            if key in change:
                change[key] = int(L[key]*change[key])
        modes = modes_for(pop, 2, **change)
    got = WM.solve_rectangle(modes, tuple(rect), L, QUAD)
    root = max((r['s'] for r in got['roots']), key=lambda z: z.real, default=None)
    return dict(name=pop['name'], variant=variant, root=None if root is None else [root.real, root.imag],
                counts_agree=got['counts_agree'], kernel_nodes=int(len(modes.nodes)), basis_size=modes.basis_size,
                kernel_condition=modes.condition)


# ---------------------------------------------------------------- L5's negative control
def gate_L5_control(pop, rect):
    """The defect found in spectrum.located(): the 2x2 determinant formula applied to an n x n T(s)."""
    modes = modes_for(pop, 2)
    right = WM.solve_rectangle(modes, tuple(rect), TH['L'], QUAD)
    wrong = WM.solve_rectangle(modes, tuple(rect), TH['L'], QUAD, determinant=SP.det_many)
    return dict(name=pop['name'], rect=list(rect), located=right['located'],
                winding_with_full_determinant=right['winding']['count'],
                winding_with_2x2_formula=wrong['winding']['count'],
                what='the 2x2 determinant formula of spectrum.det_many applied to the n x n problem',
                rejected=bool(wrong['winding']['count'] != right['located'] and right['counts_agree']))


# ---------------------------------------------------------------- one unit of work in a worker
def run_task(kind, key, args):
    t0 = time.time()
    fn = dict(population=task_population, rectangle=task_rectangle, strip=task_strip, family=task_family, L1=gate_L1,
              L2=task_L2,
              L3=task_L3, L4=task_L4, L5_control=gate_L5_control,
              V1=R7.task_V1, V2=R7.task_V2, V3=R7.task_V3, V5=R7.task_V5, V6=R7.task_V6, X=R7.task_X)[kind]
    return kind, key, _safe(fn(*args)), round(time.time() - t0, 1)
