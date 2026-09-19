"""RUT-1 stage 9 (protocol-rut9.md): the units of work of the corrected response calculation.

Everything that PRODUCES a raw result lives here, as in stages 7 and 8: rut9.py assembles, evaluates the gates
and reads, and can be corrected without invalidating finished work, because the resume cache is keyed on the
hashes of this file and the modules it drives.

Every threshold is read from the machine-readable block of the protocol: this file holds no tolerance of its
own. Stage 8's library, stage 7's population builder and the owner's sampler are used BY IMPORT, unchanged.
"""
import json
import math
import os
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut8_tasks as R8        # noqa: E402  (imports stage 7's tasks, which set one thread per worker)
from rut8_tasks import PO, R7, RM, SP, TAU_FORM, TAU_KEEP, W, WM, _find, _safe, _seed, np  # noqa: E402
from scipy.interpolate import CubicSpline   # noqa: E402
import warm_response as WR     # noqa: E402

CODE = ('rut9.py', 'rut9_tasks.py', 'warm_response.py') + tuple(R8.CODE)
PROTOCOL = 'protocol-rut9.md'


def protocol_path():
    """The protocol this process reads. RUT9_PROTOCOL exists for smoke tests on a toy protocol: stage 8's
    smoke run silently picked up the real protocol because the repository copy is found first. The path and
    its hash are recorded in every archive, and the suite job refuses an archive made from any other file."""
    override = os.environ.get('RUT9_PROTOCOL')
    return Path(override) if override else _find(PROTOCOL)


def thresholds():
    text = protocol_path().read_text(encoding='utf-8')
    return json.loads(re.search(r'## Declared thresholds.*?```json\s*(\{.*?\})\s*```', text, re.S).group(1))


TH = thresholds()
REP = TH['R']
QUAD = REP['contour']


# ---------------------------------------------------------------- populations
def archived_populations():
    """Stage 7's five and stage 8's three, from their committed series: the states the simulations used."""
    import gzip
    out = dict(R8.stage7_populations())
    with gzip.open(_find('rut8-predictions.json').parent/'rut8-series'/'new-populations.json.gz', 'rt',
                   encoding='utf-8') as f:
        out.update(json.load(f))
    return out


def family_name(dE):
    return f'B_dE{dE:.3f}'


def task_family_population(dE):
    """A member of family B on stage 7's grid at the family's coupling: stage 8's construction, by import."""
    return R8.task_population(family_name(dE), dE)


def task_refined_population(name, family, dE):
    """The same population on stage 7's `all_together` refinement of its own discretization. Family A is
    DEFINED by its mean support, so its coupling is solved again; family B by the family's coupling."""
    t0 = time.time()
    spec = dict(name=name, family=family, dE=float(dE), alpha_from=None if family == 'A' else 'B_mid')
    a = R7.build(spec, None if family == 'A' else R8.alpha_B(), **R7.REFINEMENTS['all_together'])
    d = R7.describe(a)
    return dict(name=name, family=family, dE=float(dE), state=a.state(), alpha=d['alpha'], mean_radius=d['mean_radius'],
                mean_support=d['support']['mean'], grid=d['grid'], seconds=round(time.time() - t0, 1))


def task_annulus(index, refined):
    """One ring-limit annulus, solved for the declared mean support, with the exact-history ring of stage 5
    at the annulus's own mean radius and total writing rate, and at R = 1 for the negative control."""
    N = TH['N']
    dL, dE, n_r = N['ring_limit'][index]
    t0 = time.time()
    fine = N['refined_annulus']
    kw = dict(n_r=int(fine['n_r_factor']*n_r), n_L=int(fine['n_L']), n_vr=int(fine['n_vr'])) if refined else dict(n_r=int(n_r))
    a = PO.Annulus(L0=1., dL=dL, dE=dE, w=W, tau_keep=TAU_KEEP, mass=1., reach=R7.TH['S']['reach'], r_lo=.5, r_hi=1.6, **kw)
    a.solve_mean_support(N['ring_limit_support'])
    wgt = 2*np.pi*a.r*a.sigma
    mean_r = float(np.trapezoid(wgt*a.r, a.r))
    width = float(np.sqrt(np.trapezoid(wgt*(a.r - mean_r)**2, a.r)))

    def ring_root(R):
        rg = RM.Ring('two_stage', 64, a.alpha*a.mass/64, R=R, w=W, tau_keep=TAU_KEEP, tau_form=TAU_FORM)
        s = max((r['s'] for r in SP.strip(rg, 2, R7.TH['E'])['roots']), key=lambda z: z.real)
        return [float(s.real), float(s.imag - 2*rg.omega)]                 # rotating frame -> inertial

    return dict(index=int(index), refined=bool(refined), dL=dL, dE=dE, n_r=int(kw['n_r']), alpha=float(a.alpha),
                mean_radius=mean_r, radial_width=width, mean_support=float(a.support_profile()['mean']),
                ring_root_inertial=ring_root(mean_r), ring_at_stage5_radius=ring_root(1.), state=a.state(),
                seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- the representation
_FAMILIES = {}


def representation(**change):
    cfg = dict(REP['base'])
    cfg.update(change)
    return cfg


def response_for(state, key, m, cfg, harmonics=None):
    """T(s) for one population state under one representation. The orbit families are cached within a worker
    on everything they depend on; `key` names the state."""
    disk = WM.Disk(state, TAU_FORM)
    okey = (key, state['alpha'], state['dE'], len(state['r']), cfg['n_L'], cfg['n_u'], cfg['n_eta'])
    if okey not in _FAMILIES:
        _FAMILIES.clear()
        _FAMILIES[okey] = WR.sharp_families(disk, cfg['n_L'], cfg['n_u'], cfg['n_eta'], edges=True)
    fam = _FAMILIES[okey]
    inner = fam[0][0]
    nodes = WR.kernel_nodes((float(inner.rp.min()), float(inner.ra.max())), cfg['node_spacing'], cfg['node_margin'],
                            cfg['node_clip'], cfg['node_shift'])
    return WR.Response(disk, fam if cfg['edges'] else fam[:1], m, cfg['l_max'], nodes, cfg['kernel_rank_cutoff'],
                       harmonics=harmonics)


def _fastest(got):
    return max((r['s'] for r in got['roots']), key=lambda z: z.real, default=None)


def _root_row(resp, rect, residual_points=9):
    """The fastest root in a rectangle, counted and located, with the mode equation tested BETWEEN the nodes."""
    got = WM.solve_rectangle(resp, tuple(rect), REP, QUAD)
    s = _fastest(got)
    row = dict(root=None if s is None else [float(s.real), float(s.imag)], located=got['located'],
               winding=got['winding']['count'], counts_agree=got['counts_agree'], kernel_nodes=int(len(resp.nodes)),
               basis_size=resp.basis_size, kernel_condition=resp.condition, rows_kept=resp.rows_kept)
    if s is not None:
        inner = resp.o
        pts = np.linspace(float(inner.rp.min()), float(inner.ra.max()), residual_points)
        row['off_node_residual'] = float(resp.off_node_residual(s, pts).max())
    return row


def task_variant(state, key, variant, rect, refined_state=None):
    """One representation variant of one state: N1 on an annulus, N3 on a population. `variant` is 'base', a
    key of the declared variants, or a declared control rule with an optional ':variant' refinement of it."""
    t0 = time.time()
    rule, _, inner_variant = variant.partition(':')
    if rule in ('stage8_rule', 'coarse_rule'):
        change = dict(REP[rule])
        if inner_variant:
            for k, v in REP['control_refinement'][inner_variant].items():
                change[k] = change[k]*v if k == 'node_spacing' else v
    else:
        change = {} if variant == 'base' else dict(REP['variants'][variant])
    use = state
    if change.pop('refined', False):
        if refined_state is None:
            raise ValueError('the population_grid variant needs the refined state')
        use, key = refined_state, key + ':refined'
    row = _root_row(response_for(use, key, 2, representation(**change)), rect)
    row.update(key=key, variant=variant, seconds=round(time.time() - t0, 1))
    return row


# ---------------------------------------------------------------- part B: the boundary rows
def task_B1(pop, archived_root):
    """Reproduction: without the boundary rows, on stage 8's nodes, this library IS stage 8's; and with them it
    is not (the control)."""
    L8 = R8.TH['L']
    disk = WM.Disk(pop['state'], TAU_FORM)
    fam = WR.sharp_families(disk, L8['n_L'], L8['n_u'], L8['n_eta'], edges=True)
    old_orbits = WM.Orbits(disk, L8['n_L'], L8['n_u'], L8['n_eta'])
    mine = fam[0][0]
    arrays = ('L', 'u', 'E', 'rp', 'ra', 'wt', 'Om_r', 'Om_th', 'wr', 'lag', 'dwr', 'measure')
    identical = bool(all(np.array_equal(getattr(old_orbits, k), getattr(mine, k)) for k in arrays))
    old = WM.Modes(disk, old_orbits, 2, L8['l_max'], L8['node_spacing'], L8['node_margin'], L8['node_clip'],
                   L8['kernel_rank_cutoff'])
    nodes = WR.kernel_nodes((float(mine.rp.min()), float(mine.ra.max())), L8['node_spacing'], L8['node_margin'],
                            L8['node_clip'])
    without = WR.Response(disk, fam[:1], 2, L8['l_max'], nodes, L8['kernel_rank_cutoff'])
    with_rows = WR.Response(disk, fam, 2, L8['l_max'], nodes, L8['kernel_rank_cutoff'])
    s8 = complex(*archived_root)
    det = lambda resp: (lambda z: np.linalg.det(resp.T(z)))
    s0 = SP.polish(det(without), s8, 10*REP['polish_movement'])
    s1 = SP.polish(det(with_rows), s0, 10*REP['polish_movement'])
    parts = {label: with_rows.shift(s0, disk.alpha*disk.H(s0)*with_rows.reduced(s0, family=k))
             for k, label in ((1, 'energy_edge'), (2, 'lower_L_edge'), (3, 'upper_L_edge'))}
    size = {label: float(np.linalg.norm(with_rows.reduced(s0, family=k))/np.linalg.norm(with_rows.reduced(s0, family=0)))
            for k, label in ((1, 'energy_edge'), (2, 'lower_L_edge'), (3, 'upper_L_edge'))}
    return dict(name=pop['name'], orbit_arrays_identical=identical, nodes_identical=bool(np.array_equal(old.nodes, nodes)),
                T_difference=float(np.max(np.abs(old.T(s8) - without.T(s8)))),
                T_difference_with_boundary_rows=float(np.max(np.abs(old.T(s8) - with_rows.T(s8)))),
                archived_root=[s8.real, s8.imag], root_without_rows=[s0.real, s0.imag], root_with_rows=[s1.real, s1.imag],
                root_difference_from_archive=float(abs(s0 - s8)),
                shift=[float((s1 - s0).real), float((s1 - s0).imag)], growth_relative=float(s1.real/s0.real - 1),
                first_order={k: [float(v.real), float(v.imag)] for k, v in parts.items()},
                boundary_over_interior=size, edge_orbits=[f[0].n for f in fam[1:]])


def _smoothstep(z):
    z = np.clip(z, 0., 1.)
    return z*z*z*(10. - 15.*z + 6.*z*z)


def _tapered(disk, eps, n_panel, cfg, which, step):
    """The population with each declared bound replaced by a quintic smoothstep `eps` widths wide, on a
    quadrature with a dedicated panel per taper. Its derivatives are CENTRAL DIFFERENCES OF THE FUNCTION."""
    eE, eL = eps*disk.dE, eps*disk.dL
    L_lo, L_hi = disk.L0 - disk.reach*disk.dL, disk.L0 + disk.reach*disk.dL
    u_max = np.sqrt(2*disk.reach*disk.dE)
    u_a = np.sqrt(2*(disk.reach*disk.dE - eE)) if 'energy' in which else u_max
    cuts = [L_lo] + ([L_lo + eL] if 'lower_L' in which else []) + ([L_hi - eL] if 'upper_L' in which else []) + [L_hi]
    counts = ([n_panel] if 'lower_L' in which else []) + [cfg['n_L']] + ([n_panel] if 'upper_L' in which else [])
    Ls, wLs = zip(*[WR.gauss_legendre(a, b, n) for a, b, n in zip(cuts[:-1], cuts[1:], counts)])
    us, wUs = zip(*([WR.gauss_legendre(0., u_a, cfg['n_u'])] +
                    ([WR.gauss_legendre(u_a, u_max, n_panel)] if 'energy' in which else [])))
    orb = WR.OrbitSet(disk, (np.concatenate(Ls), np.concatenate(wLs)), (np.concatenate(us), np.concatenate(wUs)), cfg['n_eta'])

    def f(E, L, Ec):
        x, y = (L - disk.L0)/disk.dL, (E - Ec)/disk.dE
        out = disk.amplitude*np.exp(-.5*x*x - .5*y*y)
        if 'energy' in which:
            out = out*_smoothstep((Ec + disk.reach*disk.dE - E)/eE)
        if 'lower_L' in which:
            out = out*_smoothstep((L - L_lo)/eL)
        if 'upper_L' in which:
            out = out*_smoothstep((L_hi - L)/eL)
        return out

    hE, hL = step*eE, step*eL
    shifted = {Li: (disk.circular(Li - hL)[1], disk.circular(Li + hL)[1]) for Li in np.unique(orb.L)}
    Ec_m = np.array([shifted[Li][0] for Li in orb.L])
    Ec_p = np.array([shifted[Li][1] for Li in orb.L])
    dfdE = (f(orb.E + hE, orb.L, orb.Ec) - f(orb.E - hE, orb.L, orb.Ec))/(2*hE)
    dfdL = (f(orb.E, orb.L + hL, Ec_p) - f(orb.E, orb.L - hL, Ec_m))/(2*hL)
    return orb, (lambda o, nu, m: nu*dfdE[:, None] + m*dfdL[:, None])


def task_B2(pop, s_pair, which):
    """The boundary rows against the limit of a resolved taper. `which` names the tapered bounds."""
    B = TH['B']
    t0 = time.time()
    cfg = representation()
    s = complex(*s_pair)
    sharp = response_for(pop['state'], pop['name'], 2, cfg)
    disk, nodes = sharp.d, sharp.nodes
    family_of = dict(energy=1, lower_L=2, upper_L=3)
    interior = sharp.reduced(s, family=0)
    target = sum(sharp.reduced(s, family=family_of[w]) for w in which)
    rows, diffs = [], []
    for eps in B['taper_eps']:
        orb, F_of = _tapered(disk, eps, int(B['taper_panel']), cfg, which, B['taper_relative_step'])
        tap = WR.Response(disk, [(orb, F_of)], 2, cfg['l_max'], nodes, cfg['kernel_rank_cutoff'])
        D = tap.reduced(s) - interior
        diffs.append(D)
        rows.append(dict(eps=eps, with_rows=float(np.linalg.norm(D - target)/np.linalg.norm(target)),
                         without_rows=float(np.linalg.norm(D)/np.linalg.norm(target)), orbits=int(orb.n)))
    ratio = B['taper_eps'][0]/B['taper_eps'][1]
    if not (abs(ratio - 2.) < 1e-12 and abs(B['taper_eps'][1]/B['taper_eps'][2] - 2.) < 1e-12):
        raise ValueError('the declared taper widths must halve: the extrapolation assumes it')
    limit = (8*diffs[2] - 6*diffs[1] + diffs[0])/3.                 # Richardson, linear and quadratic terms removed
    norm = float(np.linalg.norm(target))
    return dict(name=pop['name'], which=list(which), s=[s.real, s.imag], rows=rows,
                boundary_over_interior=float(norm/np.linalg.norm(interior)),
                extrapolated_with_rows=float(np.linalg.norm(limit - target)/norm),
                extrapolated_without_rows=float(np.linalg.norm(limit)/norm), seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- part N4 and part F: strips and the family
def task_strip(pop, m, config, shifted):
    """Every sub-rectangle of one declared strip under one declared configuration."""
    t0 = time.time()
    change = {} if config == 'base' else dict(TH['N']['root_free_configs'][config])
    resp = response_for(pop['state'], pop['name'], m, representation(**change))
    rects = []
    for rect in WM.strip_rectangles(m, TH['N'], shifted):
        got = WM.solve_rectangle(resp, tuple(rect), REP, QUAD)
        rects.append(dict(rect=got['rect'], roots=[[float(r['s'].real), float(r['s'].imag)] for r in got['roots']],
                          located=got['located'], winding=got['winding']['count'], counts_agree=got['counts_agree'],
                          gap=got['gap'], relative_residuals=[r['relative_residual'] for r in got['roots']],
                          moved=[r['moved'] for r in got['roots']], saturated=got['saturated'], nodes=got['nodes']))
    return dict(name=pop['name'], m=int(m), config=config, shifted=bool(shifted), rectangles=rects,
                kernel_nodes=int(len(resp.nodes)), basis_size=resp.basis_size, rows_kept=resp.rows_kept,
                seconds=round(time.time() - t0, 1))


def task_family(dE):
    """One member of family B under the corrected representation: its m = 2 root in the declared rectangle,
    with and without the boundary rows."""
    pop = task_family_population(dE)
    rect = TH['N']['family_rect']
    row = _root_row(response_for(pop['state'], pop['name'], 2, representation()), rect)
    bare = _root_row(response_for(pop['state'], pop['name'], 2, representation(edges=False)), rect)
    d = pop['description']
    row.update(dE=float(dE), alpha=d['alpha'], mean_support=d['support']['mean'], mean_radius=d['mean_radius'],
               sigma_r=d['dispersions_at_mean_radius']['sigma_r'], root_without_boundary_rows=bare['root'])
    return row


# ---------------------------------------------------------------- part S: what the root is sensitive to
def task_S(pop, guess):
    S = TH['S']
    resp = response_for(pop['state'], pop['name'], 2, representation())
    disk = resp.d
    s = SP.polish(lambda z: np.linalg.det(resp.T(z)), complex(*guess), 10*REP['polish_movement'])
    aH = disk.alpha*disk.H(s)
    frac = S['relative_change']
    v, u, ratio = resp.eigenvectors(s)
    whole = resp.shift(s, frac*aH*resp.reduced(s))
    resolved = SP.polish(lambda z: np.linalg.det(np.eye(resp.basis_size) + (1 + frac)*disk.alpha*disk.H(z)*resp.reduced(z)),
                         s, 100*REP['polish_movement']) - s
    # the control: the memory transfer function held fixed in dT/ds
    dT = frac*aH*resp.reduced(s)
    wrong = complex(-(u.conj() @ dT @ v)/(u.conj() @ (aH*resp.reduced(s, derivative=True)) @ v))
    c = lambda z: [float(z.real), float(z.imag)]
    harmonics = {str(int(l)): c(resp.shift(s, frac*aH*resp.reduced(s, harmonic=int(l)))) for l in resp.l}
    families = {label: c(resp.shift(s, frac*aH*resp.reduced(s, family=k)))
                for k, label in enumerate(('interior', 'energy_edge', 'lower_L_edge', 'upper_L_edge'))}
    d = disk
    memory = {}
    for label, tk, tf in (('tau_keep', (1 + frac)*d.tau_keep, d.tau_form), ('tau_form', d.tau_keep, (1 + frac)*d.tau_form)):
        H2 = tk/((1 + s*tk)*(1 + s*tf))
        memory[label] = c(resp.shift(s, d.alpha*(H2 - d.H(s))*resp.reduced(s)))
    return dict(name=pop['name'], root=c(s), smallest_singular_ratio=ratio, relative_change=frac,
                left_is_conjugate_of_right=float(1. - abs(np.vdot(u, v.conj()))/(np.linalg.norm(u)*np.linalg.norm(v))),
                whole_response=c(whole), whole_response_resolved=c(resolved),
                first_order_error=float(abs(whole - resolved)/abs(resolved)),
                control_with_H_held_fixed=c(wrong), control_error=float(abs(wrong - resolved)/abs(resolved)),
                growth_per_percent_of_response=float(whole.real/s.real/frac/100.),
                by_harmonic=harmonics, by_family=families, memory_times=memory)


# ---------------------------------------------------------------- part E: the mode's own eigenfunction, in time
def _integrate(disk, psi, rgrid, m, s, eps, x, v, efolds, h):
    """+eps and -eps on identical bodies from t = -efolds/gamma to 0, in the frozen potential plus
    eps Re[psi(r) e^{i m theta + s t}] as a perturbation of the memory field (a force +grad)."""
    pr, pi = CubicSpline(rgrid, psi.real), CubicSpline(rgrid, psi.imag)
    dpr, dpi = pr.derivative(), pi.derivative()
    gamma, beta = s.real, s.imag
    steps = int(round(efolds/gamma/h))
    n = len(x)
    q, p = np.concatenate([x, x]), np.concatenate([v, v])
    sign = np.concatenate([np.ones(n), -np.ones(n)])

    def accel(q, t):
        r = np.hypot(q[:, 0], q[:, 1])
        th = np.arctan2(q[:, 1], q[:, 0])
        ur = q/r[:, None]
        amp = sign*eps*math.exp(gamma*t)
        cs, sn = np.cos(m*th + beta*t), np.sin(m*th + beta*t)
        fr = -disk.dphi(r) + amp*(dpr(r)*cs - dpi(r)*sn)
        ft = -amp*(m/r)*(pr(r)*sn + pi(r)*cs)
        return fr[:, None]*ur + ft[:, None]*np.stack([-ur[:, 1], ur[:, 0]], axis=1)

    t = -steps*h
    a = accel(q, t)
    for _ in range(steps):
        ph = p + .5*h*a
        q = q + h*ph
        t += h
        a = accel(q, t)
        p = ph + .5*h*a
    if not (np.isfinite(q).all() and np.isfinite(p).all()):
        raise FloatingPointError('non-finite test-body state')
    return q, p, sign, steps


def mode_field(resp, s, peak):
    """The mode's field on a fine radial grid, its node coefficients and both eigenvectors, scaled so that the
    largest |delta C_m(r)| is `peak` and the largest node coefficient is real and positive."""
    v, u, ratio = resp.eigenvectors(s)
    c = resp.Q @ v
    rgrid = np.linspace(resp.d.r[0], resp.d.r[-1], 4000)
    psi = sum(cb*resp.kernel(rb, rgrid) for cb, rb in zip(c, resp.nodes))
    scale = peak/np.max(np.abs(psi))*np.exp(-1j*np.angle(c[np.argmax(np.abs(c))]))
    return rgrid, psi*scale, c*scale, v*scale, u, ratio


def task_E(pop, guess, batch, variant):
    """One batch: the predicted eigenfunction imposed at the predicted complex rate on bodies drawn with the
    sampler of record. `variant` is 'base'; 'replay', a small batch the suite job recomputes exactly; or one of
    the refinement group -- 'refine_base', 'half_step', 'double_eps' -- which share ONE draw and one duration,
    so that their differences are deterministic and carry no sampling noise. (A longer switch-on is NOT such a
    refinement: the same bodies reach t = 0 in other orbital phases, which is another sample of the response.)"""
    E = TH['E']
    t0 = time.time()
    resp = response_for(pop['state'], pop['name'], 2, representation())
    disk = resp.d
    s = SP.polish(lambda z: np.linalg.det(resp.T(z)), complex(*guess), 10*REP['polish_movement'])
    rgrid, psi, c, v, u, ratio = mode_field(resp, s, E['peak_field'])
    eps = E['eps']*(2. if variant == 'double_eps' else 1.)
    h = E['step']*(.5 if variant == 'half_step' else 1.)
    efolds = E['efoldings']
    group = 'batch' if variant in ('base', 'replay') else 'refinement'      # the refinements share one draw
    n = int(E['bodies'] if variant == 'base' else E['replay_bodies'] if variant == 'replay' else E['refinement_bodies'])
    a = PO.Annulus.from_state(pop['state'])
    x, vel = a.draw(n, seed=_seed('E', pop['name'], group, int(batch)))
    q, p, sign, steps = _integrate(disk, psi, rgrid, 2, s, eps, x, vel, efolds, h)
    r = np.hypot(q[:, 0], q[:, 1])
    ph = np.exp(-2j*np.arctan2(q[:, 1], q[:, 0]))
    measured = np.array([(resp.kernel(ra, r)*ph*sign).sum() for ra in resp.nodes])*a.mass/n/eps
    predicted = -resp.M(s) @ c
    cx = lambda z: [float(z.real), float(z.imag)]
    return dict(name=pop['name'], batch=int(batch), variant=variant, bodies=n, steps=int(steps), eps=eps, step=h,
                efoldings=efolds,
                root=cx(s), kernel_nodes=[float(z) for z in resp.nodes], measured=[cx(z) for z in measured],
                predicted=[cx(z) for z in predicted], smallest_singular_ratio=ratio, seconds=round(time.time() - t0, 1))


def implied_shift(resp, s, v, u, measured, predicted, scale=1.):
    """delta_s = -(u^H dT v)/(u^H T' v) with dT v = alpha H Q^T [(-measured) - scale M c], M c = -predicted."""
    d = resp.d
    dTv = d.alpha*d.H(s)*(resp.Q.T @ (scale*np.asarray(predicted) - np.asarray(measured)))
    return complex(-(u.conj() @ dTv)/(u.conj() @ resp.dT_ds(s) @ v))


def task_E_reduce(pop, guess, batches):
    """Every batch's measured node response turned into the root shift it implies, in one place, so that the
    projection uses ONE pair of eigenvectors. Also the control: the prediction scaled by the declared factor."""
    E = TH['E']
    resp = response_for(pop['state'], pop['name'], 2, representation())
    s = SP.polish(lambda z: np.linalg.det(resp.T(z)), complex(*guess), 10*REP['polish_movement'])
    _, _, c, v, u, _ = mode_field(resp, s, E['peak_field'])
    predicted = -resp.M(s) @ c
    rows = []
    for b in batches:
        measured = np.array([complex(*z) for z in b['measured']])
        ds = implied_shift(resp, s, v, u, measured, predicted)
        control = implied_shift(resp, s, v, u, measured, predicted, scale=E['control_scale'])
        regenerated = resp.d.alpha*resp.d.H(s)*measured
        imposed = resp.P @ c
        big = np.abs(imposed) > .1*np.max(np.abs(imposed))
        rows.append(dict(batch=b['batch'], variant=b['variant'], delta_s=[float(ds.real), float(ds.imag)],
                         control_delta_s=[float(control.real), float(control.imag)],
                         worst_node=float(np.max(np.abs(measured - predicted))/np.max(np.abs(predicted))),
                         regenerated_over_imposed_mean=[float(np.mean((regenerated/imposed)[big]).real),
                                                        float(np.mean((regenerated/imposed)[big]).imag)]))
    return dict(name=pop['name'], root=[float(s.real), float(s.imag)], rows=rows)


# ---------------------------------------------------------------- one unit of work in a worker
def run_task(kind, key, args):
    t0 = time.time()
    fn = dict(family_population=task_family_population, refined_population=task_refined_population, annulus=task_annulus,
              variant=task_variant, B1=task_B1, B2=task_B2, strip=task_strip, family=task_family, S=task_S, E=task_E,
              E_reduce=task_E_reduce)[kind]
    return kind, key, _safe(fn(*args)), round(time.time() - t0, 1)
