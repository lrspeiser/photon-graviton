"""RC-2, part 1: a time-dependent donor-and-companion calculation in an open region (see protocol.md, Amendment 1).

    python rc2a.py [--canonical] [--output-dir DIR]

The reference case: cold collisionless decay at 3 km/s. The Milky Way's companions, bound and unbound, and the donor's
depletion are evolved together for 10 Gyr in open regions, in three declared stages:
  1. controls (C1a-C1c, C3, C3b, D4; C2 in every run) and the region-size comparison, 150-2,400 kpc, at 2B-F1's best
     sampled rate for the matching combination (collisionless, bath gravity on);
  2. a direct RMSE search at 1,200 kpc, only if C1a passes and C1c demonstrates convergence;
  3. the region-size comparison over 600-2,400 kpc at the selected rate, frozen.
Convergence is judged on the gravitational signal with a precision requirement, and every tracer's birth radius, birth
time and angular momentum are traced to where it ends. RC2A_WORKERS sets the pool, RC2A_BUDGET each run's wall-clock
budget, and RC2A_SMOKE=1 runs a reduced set (in which stages 2 and 3 run regardless, to exercise them). Seeds are
drawn in a fixed order, and results return in the order the tasks were issued.
"""
import json
import math
import os
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'companion-source'))
import f1 as F  # noqa: E402  (also puts field, mc, formation and evidence_io on the path)
import formation as FM  # noqa: E402
import mc  # noqa: E402
import evidence_io  # noqa: E402
import open_region as OR  # noqa: E402

SMOKE = os.environ.get('RC2A_SMOKE') == '1'
V_REF = 3.                                                                # km/s: the reference case
REGIONS = (150., 300.) if SMOKE else (150., 300., 600., 1200., 2400.)    # kpc: the open region's radius R_comp
R_SEARCH = 300. if SMOKE else 1200.
STAGE3_REGIONS = (150., 300.) if SMOKE else (600., 1200., 2400.)
R_C1, R_C1C, R_C3 = REGIONS[-1], R_SEARCH, R_SEARCH                     # the controls' regions
SEEDS = 2 if SMOKE else 3
APERTURES = (50., 100., 150., 200., 300.)                                 # kpc
T_REF = 1. if SMOKE else 10.                                              # Gyr
SNAPS = (.25, .5, .75, 1.) if SMOKE else (2., 5., 8., 10.)
T_OFF = .5 if SMOKE else 5.                                               # Gyr: the production episode of C3 and C3b
N_REF = 2000 if SMOKE else 16000         # field births over the span inside the reference region (2B-F1's count)
N_HI = 4*N_REF                            # C1b's and C1c's higher resolution
BUDGET = float(os.environ.get('RC2A_BUDGET', '900' if SMOKE else '14400'))
WORKERS = int(os.environ.get('RC2A_WORKERS', '4' if SMOKE else '8'))
BRACKET, MAX_EXT = (-4, -2, 0, 2, 4), 2  # rates q0 3^(e/8) and outward extensions per side, as in the 2B-F1 revision
COMBO = 'v_d {:g} km/s | 0 cm2/g | bath gravity on'
FLOOR = dict(rmse_38=1., slope_8_20=.05)  # declared tolerances: km/s, and the 8-20 kpc slope
V_TOL = 2.                                # km/s: circular speeds at the apertures and rotation speeds at the Eilers radii
MASS_FLOOR = .05                          # the companion inventory inside an aperture (reported, not decisive)
Z_DRAWS = 4.                              # D1 and C2: births against their expectation, in Poisson standard errors
C1_LIMIT = .1                             # C1a: without baryons, at most this fraction of the reference's companions
SOURCE = Path(os.environ.get('F1R_SOURCE', HERE.parent/'companion-source'/'f1-results.json'))
T0 = time.time()
_MW = []


def log(msg):
    print(f'[{time.time() - T0:9.1f}s] {msg}', flush=True)


def mw():
    if not _MW:
        _MW.append(F.mw_data())
    return _MW[0]


def horizon(key):
    """The radius from which matter at rest falls onto the system's baryons, taken as a point mass, within T."""
    T = T_REF*mc.PER_GYR
    return (2*mc.G*float(FM.system(key)['M_b'])*(2*T/math.pi)**2)**(1/3)


def n_max_for(R, n):
    """Four times the births the mass rule expects (they grow as R^1.5 beyond the reference region), so that the engine
    never thins the population."""
    return int(4*n*max(1., (R/FM.system('MW')['R_b'])**1.5))


def make(spec, s, depletion=True, track_unbound=True):
    """The donor model for one run: 2B-F1's constructor, the region's grid at the reference spacing, tracer masses fixed
    by the Milky Way's 2B-F1 zone (mass_ref None restores 2B-F1's rule), and the quiet representation unless spec says
    otherwise. spec may remove the baryons (C1) or end production at t_off_gyr (C3, C3b). Ratios are to the baryons in
    the 2B-F1 zone."""
    dmax = F.DELTA_MAX[s['kind']]
    s0 = FM.system(spec['key'])
    n_grid, n_shell = OR.grid_for(s['R_b'], s0['R_b'], s['r_half'])
    t_off = spec.get('t_off_gyr')
    model = OR.DonorModel(s['r'], np.zeros_like(s['M']) if spec.get('no_baryons') else s['M'], s['R_b'], s['r_half'],
                          np.array([spec['v_d']]), np.array([1.]), 0., spec['q']*dmax/2, q=spec['q']/mc.PER_GYR,
                          v_d=spec['v_d'], self_gravity=True, freeze=False, seed=spec['rng'],
                          apertures=tuple(a for a in APERTURES if a < s['R_b']), depletion=depletion,
                          mass_ref_kpc=spec.get('mass_ref', s0['R_b']), track_unbound=track_unbound,
                          quiet=spec.get('quiet', True), t_off=None if t_off is None else t_off*mc.PER_GYR,
                          mass_unit=float(s0['M_b']), score_radii=mw()['R'], n_grid=n_grid, n_shell=n_shell)
    return model, dmax


# ---- scores
def scores(M, D):
    """2B-F1's Milky Way scores for a net non-baryonic enclosed mass M at the 38 Eilers radii."""
    M = np.asarray(M, float)
    v = np.sqrt(np.maximum(D['vb']**2 + mc.G*M/D['R'], 0.))
    rmse = lambda n=None: float(np.sqrt(np.mean((v[:n] - D['y'][:n])**2)))
    sel = D['sel']
    slope = float(np.polyfit(np.log(D['R'][sel]), np.log(M[sel]), 1)[0]) if (M[sel] > 0).all() else None
    return dict(rmse_38=rmse(), rmse_inner20=rmse(20), slope_8_20=slope, predicted_kms=v.tolist())


def mw_score(model, D):
    """2B-F1's Milky Way scores, with the donor's depletion inside each radius subtracted from the companions'."""
    rr = np.linalg.norm(model.x, axis=1) if len(model.x) else np.zeros(0)
    o = np.argsort(rr)
    rs, cm = rr[o], np.r_[0., np.cumsum(model.m[o])]
    dep = model.depleted_mass(D['R'])
    M = cm[np.searchsorted(rs, D['R'], side='right')] - dep
    out = scores(M, D)
    out.update(extra_mass_over_required={f'{r0:g} kpc': float(np.interp(r0, D['R'], M)/np.interp(r0, D['R'], D['M_req']))
                                         for r0 in (6., 10., 15., 20., 24.)},
               M_extra_at_R=M.tolist(), depletion_at_R=dep.tolist())
    return out


def gates(sc, D):
    """G1 and G2 as in 2B-F1. There is no cost gate: the inventory, the contrast and the source's energy are reported
    separately instead."""
    if sc is None or sc.get('rmse_38') is None:
        return None
    g1 = bool(sc['rmse_38'] <= F.G1_RMSE and sc['rmse_inner20'] <= D['baryons_rmse_inner20'])
    g2 = bool(sc['slope_8_20'] is not None and abs(sc['slope_8_20'] - D['slope_req']) <= F.G2_TOL)
    return dict(G1_fit=g1, G2_shape=g2, profile_gate_passed=g1 and g2)


# ---- tasks
def task(spec):
    try:
        return _task(spec)
    except Exception as e:                       # keep the pool alive; the driver reports and stops
        return dict(spec=spec, status=f'error: {type(e).__name__}: {e}', traceback=traceback.format_exc())


def d4(spec, T_gyr=None):
    """D4: with the depletion off, unbound births removed at once, Poisson draws and 2B-F1's mass rule, the donor model
    (its birth tags included) against 2B-F1's field model (collisionless, bath gravity off), from the same seed and rate,
    bit for bit."""
    T = T_REF if T_gyr is None else T_gyr
    s = F.system_for(dict(key=spec['key']))
    runs = {}
    for name in ('field', 'donor'):
        if name == 'field':
            model, dmax = F.make_model(dict(spec, sm=0., gravity=False), s)
        else:
            model, dmax = make(dict(spec, mass_ref=None, quiet=False), s, depletion=False, track_unbound=False)
        res = model.run(T, dmax, [T], n_target=spec['n'], n_field=spec['n'], n_max=4*spec['n'], budget_s=BUDGET)
        runs[name] = (model, res['ledger'])
    (a, La), (b, Lb) = runs['field'], runs['donor']
    same_state = all(np.array_equal(u, w) for u, w in ((a.x, b.x), (a.v, b.v), (a.m, b.m)))
    common = sorted(set(La) & set(Lb))
    same_ledger = all(La[k] == Lb[k] for k in common)
    tags_follow = bool(b.tags is not None and len(b.tags) == len(b.m) and np.isfinite(b.tags[:, :3]).all())
    return dict(tracers=int(len(a.m)), same_state=bool(same_state), same_ledger=bool(same_ledger),
                ledger_entries_compared=len(common), tags_follow_tracers=tags_follow,
                passed=bool(same_state and same_ledger and tags_follow and len(a.m) > 0))


def c2_z(rows):
    return [abs(r['births_minus_depletion_z']) for r in rows if r['births_minus_depletion_z'] is not None]


def _task(spec):
    t0 = time.time()
    if spec['role'] == 'D4':
        return dict(spec=spec, status='completed', **d4(spec), wall_s=time.time() - t0)
    s = F.system_for(dict(key=spec['key'], R_b=spec['R_comp']))
    model, dmax = make(spec, s)
    out = dict(spec=spec)
    try:
        res = model.run(T_REF, dmax, SNAPS, n_target=spec['n'], n_field=spec['n'], n_max=spec['n_max'], budget_s=BUDGET)
    except TimeoutError as e:
        return dict(out, status=f'not completed: {e}', wall_s=time.time() - t0)
    except mc.RunawayError as e:
        return dict(out, status=f'not completed: runaway, {e}', wall_s=time.time() - t0)
    L = res['ledger']
    Mu = model.mass_unit
    born, made, esc = (float(L.get(k, 0.)) for k in ('field_born_mass', 'field_made_mass', 'escaped_mass'))
    expected, m2 = float(L.get('field_expected_births', 0.)), float(L.get('field_birth_m2', 0.))
    dep = float(model.depleted_mass(model.R_b))
    M = float(model.m.sum())
    D = mw()
    ap = model.aperture_table()
    hist = []
    for h in res['history']:
        e = dict(t_Gyr=h['t_Gyr'], tracers=h['tracers'], companions=h['M_confined'], escaped=h['escaped_mass'],
                 depletion=h['depletion_in_region'], apertures=h['apertures'], birth_history=h['birth_history'])
        if not spec.get('no_baryons'):
            M_R = np.asarray(h['companions_at_R']) - np.asarray(h['depletion_at_R'])
            e['profile'] = {k: x for k, x in scores(M_R, D).items() if k != 'predicted_kms'}
        hist.append(e)
    rr = np.linalg.norm(model.x, axis=1) if len(model.x) else np.zeros(0)
    out.update(status='completed', tracers=int(len(model.m)), apertures=ap, history=hist,
               birth_history=model.birth_history(),
               inner=dict(companions_within_25=float(model.m[rr < 25.].sum()), depletion_within_25=float(model.depleted_mass(25.))),
               ledger=dict(made=made, born=born, expected_births=expected,
                           bound_at_birth_fraction=float(L.get('field_expected_bound', 0.))/max(made, 1e-300),
                           inside_at_T_fraction=M/max(born, 1e-300), escaped=esc, thinning_mass=float(L['M_thinning']),
                           thinnings=int(L['thinnings']), E_births=float(L.get('E_field_births', 0.)),
                           E_escaped=float(L['E_orbit_escapes']), E_potential_work=float(L['E_potential_work']),
                           closure_max=float(L['closure_max']), capped=float(L['capped']), steps=int(L['steps'])),
               D1=dict(mass_balance=abs(M + esc - born - float(L['M_thinning']))/max(born, 1e-300),
                       depletion_vs_production=abs(dep - made)/max(made, 1e-300),
                       draws_z=(born - expected)/math.sqrt(max(m2, 1e-300)), energy_closure=float(L['closure_max']),
                       tags_follow_tracers=bool(model.tags is not None and len(model.tags) == len(model.m))),
               C2=max(c2_z(ap) + [z for h in res['history'] for z in c2_z(h['apertures'])], default=0.),
               D2=OR.depletion_potential_check(model, T_REF*mc.PER_GYR),
               region=dict(companions_over_mass_unit=M/Mu, depletion_over_mass_unit=dep/Mu, net_over_mass_unit=(M - dep)/Mu,
                           escaped_over_mass_unit=esc/Mu),
               wall_s=time.time() - t0)
    if not spec.get('no_baryons'):
        out['mw'] = mw_score(model, D)
        out['gates'] = gates(out['mw'], D)
    return out


def short_run():
    """The suite's regression anchor: the Milky Way at v_d = 10 km/s in a 600 kpc region, q = 1.5e3, 0.5 Gyr, quiet."""
    spec = dict(role='anchor', key='MW', v_d=10., q=1.5e3, R_comp=600., rng=4343)
    s = F.system_for(dict(key='MW', R_b=spec['R_comp']))
    model, dmax = make(spec, s)
    res = model.run(.5, dmax, [.5], n_target=2000, n_field=2000, n_max=40000, budget_s=600)
    L = res['ledger']
    ap = model.aperture_table()
    M, born, esc = float(model.m.sum()), float(L['field_born_mass']), float(L['escaped_mass'])
    dep, made = float(model.depleted_mass(model.R_b)), float(L['field_made_mass'])
    return dict(M=M, tracers=int(len(model.m)), born=born, escaped=esc, companions_within_100=ap[1]['companions'],
                net_within_300=ap[4]['net'], closure=float(L['closure_max']),
                mass_balance=abs(M + esc - born - float(L['M_thinning']))/born, depletion_vs_production=abs(dep - made)/made,
                c2_max=max(c2_z(ap), default=0.), tags_follow_tracers=bool(len(model.tags) == len(model.m)))


def run_tasks(specs, label):
    if not specs:
        return []
    log(f'{label}: {len(specs)} runs on {min(WORKERS, len(specs))} workers')
    out = []
    with Pool(min(WORKERS, len(specs))) as pool:
        for i, r in enumerate(pool.imap_unordered(task, specs), 1):
            sp = r['spec']
            extra = ''
            if r['status'] == 'completed' and r.get('mw'):
                w, sl = r['mw'], r['mw']['slope_8_20']
                extra = (f"RMSE {w['rmse_38']:.1f} slope {'-' if sl is None else f'{sl:.2f}'} "
                         f"net {r['region']['net_over_mass_unit']:.3g} tracers {r['tracers']} {r['wall_s']:.0f}s")
            elif r['status'] == 'completed' and sp['role'] == 'D4':
                extra = f"passed={r['passed']}"
            elif r['status'] == 'completed':
                extra = (f"inner {r['inner']['companions_within_25']:.3g} net {r['region']['net_over_mass_unit']:.3g} "
                         f"tracers {r['tracers']} {r['wall_s']:.0f}s")
            log(f"  [{i}/{len(specs)}] {sp['role']} v_d={sp['v_d']:g} R_comp={sp['R_comp']:.0f} q={sp['q']:.4g} "
                f"n={sp.get('n')} {'quiet' if sp.get('quiet', sp['role'] != 'D4') else 'random'}: {r['status'][:90]} {extra}")
            if r['status'].startswith('error'):
                log('ERROR TRACEBACK\n' + r.get('traceback', ''))
                raise RuntimeError(f"task failed: {r['status']}")
            out.append(r)
    # imap_unordered hands results back as they finish; restore the order the tasks were issued in (each spec's seed is
    # unique and increasing), so that nothing downstream, the next round's seeds in particular, depends on timing
    out.sort(key=lambda r: r['spec']['rng'])
    return out


def record(r):
    return dict({k: x for k, x in r.items() if k not in ('spec', 'traceback')},
                **{k: r['spec'][k] for k in ('role', 'v_d', 'q', 'R_comp', 'rng', 'step', 'n', 'quiet', 't_off_gyr',
                                             'no_baryons') if k in r['spec']})


# ---- statistics, signal-level convergence and controls
def stats(vals):
    vals = [v for v in vals if v is not None]
    n = len(vals)
    return dict(mean=float(np.mean(vals)) if n else None, sd=float(np.std(vals, ddof=1)) if n > 1 else None, n=n)


def verdict(a, b, tol):
    """The declared comparison of two seed means (agree: within two standard errors of the difference or the tolerance,
    whichever is larger) and the precision requirement (resolved: two standard errors below the tolerance)."""
    if a['mean'] is None or b['mean'] is None or a['n'] < 2 or b['n'] < 2:
        return dict(agree=False, resolved=False, diff=None, se=None, tol=tol)
    se = math.sqrt(a['sd']**2/a['n'] + b['sd']**2/b['n'])
    d = b['mean'] - a['mean']
    return dict(agree=bool(abs(d) <= max(2*se, tol)), resolved=bool(2*se <= tol), diff=d, se=se, tol=tol)


def compare(A, B):
    """Two sets of completed Milky Way runs (seeds) compared on every declared quantity: the signal (RMSE, slope, the
    rotation speed at each Eilers radius, the circular speed and net contrast inside each aperture both contain) and the
    positive inventory. Demonstrated convergence needs every signal quantity to agree and be resolved; the declared
    (statistical) verdict needs the RMSE, the slope and the inventories to agree."""
    D = mw()
    out = dict(signal={}, inventory={})
    for k, tol in FLOOR.items():
        out['signal'][k] = verdict(stats([r['mw'][k] for r in A]), stats([r['mw'][k] for r in B]), tol)
    worst = None
    for i, R0 in enumerate(D['R']):
        v = verdict(stats([r['mw']['predicted_kms'][i] for r in A]), stats([r['mw']['predicted_kms'][i] for r in B]), V_TOL)
        if worst is None or (v['agree'], v['resolved']) < (worst['agree'], worst['resolved']) or \
                (v['diff'] is not None and worst['diff'] is not None and abs(v['diff']) > abs(worst['diff'])):
            worst = dict(v, r_kpc=float(R0))
        out['signal'][f'v_rot {R0:.2f} kpc'] = v
    out['worst_rotation_speed'] = worst
    common = [a['r_kpc'] for a in A[0]['apertures'] if a['r_kpc'] in {b['r_kpc'] for b in B[0]['apertures']}]
    for a in common:
        ia = next(i for i, x in enumerate(A[0]['apertures']) if x['r_kpc'] == a)
        ib = next(i for i, x in enumerate(B[0]['apertures']) if x['r_kpc'] == a)
        vc_a, vc_b = (stats([r['apertures'][ia]['v_circ'] for r in A]), stats([r['apertures'][ib]['v_circ'] for r in B]))
        out['signal'][f'v_circ {a:g} kpc'] = verdict(vc_a, vc_b, V_TOL)
        vc = .5*(vc_a['mean'] + vc_b['mean'])
        dM = 2*a*vc*V_TOL/mc.G                                          # the mass that moves v_circ by V_TOL
        out['signal'][f'net {a:g} kpc'] = verdict(stats([r['apertures'][ia]['net'] for r in A]),
                                                  stats([r['apertures'][ib]['net'] for r in B]), dM)
        ca, cb = stats([r['apertures'][ia]['companions'] for r in A]), stats([r['apertures'][ib]['companions'] for r in B])
        out['inventory'][f'companions {a:g} kpc'] = verdict(ca, cb, MASS_FLOOR*max(abs(ca['mean'] or 0.), abs(cb['mean'] or 0.)))
    out['demonstrated'] = bool(all(v['agree'] and v['resolved'] for v in out['signal'].values()))
    out['statistical'] = bool(all(out['signal'][k]['agree'] for k in FLOOR) and all(v['agree'] for v in out['inventory'].values()))
    out['unresolved'] = sorted(k for k, v in out['signal'].items() if v['agree'] and not v['resolved'])
    out['disagree'] = sorted(k for k, v in out['signal'].items() if not v['agree'])
    return out


def mean_birth_history(rs):
    """Seed means of the birth-history quantiles and histograms at T."""
    bh = [r['birth_history'] for r in rs if r.get('birth_history')]
    if not bh:
        return None
    out = dict(r_edges_kpc=bh[0]['r_edges_kpc'], t_edges_Gyr=bh[0]['t_edges_Gyr'], radii={})
    for R in bh[0]['radii']:
        rows = [b['radii'][R] for b in bh if b['radii'][R].get('tracers')]
        if not rows:
            continue
        e = dict(mass=float(np.mean([x['mass'] for x in rows])), seeds=len(rows),
                 arrived_within_25_fraction=float(np.mean([x['arrived_within_25_fraction'] for x in rows])))
        for name in OR.TAGS[:3]:
            e[name] = {q: float(np.mean([x[name][q] for x in rows])) for q in rows[0][name]}
        for h in ('birth_r_hist', 'birth_t_hist'):
            e[h] = np.mean([x[h] for x in rows], axis=0).tolist()
        out['radii'][R] = e
    return out


def region_summary(rs):
    out = {k: stats([r['mw'][k] for r in rs]) for k in ('rmse_38', 'rmse_inner20', 'slope_8_20')}
    out['apertures'] = {f"{a['r_kpc']:g}": {k: stats([r['apertures'][i][k] for r in rs])
                                             for k in ('companions', 'depletion', 'net', 'baryons', 'v_circ', 'born_inside',
                                                       'inflow', 'tracers')}
                        for i, a in enumerate(rs[0]['apertures'])}
    out['region'] = {k: stats([r['region'][k] for r in rs]) for k in rs[0]['region']}
    out['inner'] = stats([r['inner']['companions_within_25'] for r in rs])
    out['birth_history'] = mean_birth_history(rs)
    out['wall_s'] = stats([r['wall_s'] for r in rs])
    return out


def convergence(runs, regions):
    """The region-size comparison: every doubling of completed regions compared (compare), R_conv (demonstrated) and
    R_conv (statistical). Each needs an actual comparison with a larger completed region, so the largest region never
    qualifies by itself."""
    by = {R: [r for r in runs if r['spec']['R_comp'] == R and r['status'] == 'completed'] for R in regions}
    summ = {R: region_summary(rs) for R, rs in by.items() if len(rs) == SEEDS}
    rows = []
    for R, R2 in zip(regions[:-1], regions[1:]):
        if R not in summ or R2 not in summ:
            rows.append(dict(R_kpc=R, R2_kpc=R2, missing=True, demonstrated=False, statistical=False))
            continue
        rows.append(dict(R_kpc=R, R2_kpc=R2, **compare(by[R], by[R2])))

    def first(key):
        return next((row['R_kpc'] for i, row in enumerate(rows) if all(x[key] for x in rows[i:])), None)
    return dict(summaries={f'{R:g}': x for R, x in summ.items()}, doublings=rows, R_conv_demonstrated=first('demonstrated'),
                R_conv_statistical=first('statistical'))


def rms(vals):
    """The root-mean-square over seeds: a noise excursion can have either sign, so a seed mean could hide it."""
    vals = [v for v in vals if v is not None]
    return float(np.sqrt(np.mean(np.square(vals)))) if vals else None


def control_c1a(c1, ref):
    """C1a: the quiet source without baryons, against the reference runs in the same region at the same rate; the net
    contrast inside each aperture as a root-mean-square over the seeds."""
    done = [r for r in c1 if r['status'] == 'completed']
    if len(done) < SEEDS or len(ref) < SEEDS:
        return dict(passed=False, reason='runs missing', statuses=[r['status'] for r in c1 + ref])
    rows, ok = [], True
    for i, a in enumerate(done[0]['apertures']):
        nets = [r['apertures'][i]['net'] for r in done]
        comp = float(np.mean([r['apertures'][i]['companions'] for r in ref]))
        ratio = rms(nets)/comp if comp > 0 else None
        ok = ok and ratio is not None and ratio < C1_LIMIT
        rows.append(dict(r_kpc=a['r_kpc'], net_without_baryons=stats(nets), net_rms=rms(nets), companions_with_baryons=comp,
                         ratio=ratio))
    inner = stats([r['inner']['companions_within_25'] for r in done])
    inner_ref = float(np.mean([r['inner']['companions_within_25'] for r in ref]))
    inner_ratio = inner['mean']/inner_ref if inner_ref > 0 else None
    ok = ok and inner_ratio is not None and inner_ratio < C1_LIMIT
    return dict(apertures=rows, companions_within_25=inner, companions_within_25_reference=inner_ref,
                inner_ratio=inner_ratio, limit=C1_LIMIT, passed=bool(ok), growth=growth(done), runs=[record(r) for r in c1])


def growth(runs):
    """The seed-mean net contrast inside each aperture at each snapshot."""
    return [dict(t_Gyr=h['t_Gyr'], net=[stats([r['history'][k]['apertures'][i]['net'] for r in runs])
                                        for i in range(len(h['apertures']))]) for k, h in enumerate(runs[0]['history'])]


def control_c1b(lo, hi):
    """C1b: random sampling without baryons at two resolutions: the net contrast inside each aperture and the excess of
    the companions inside 25 kpc over the uniform expectation, as root-mean-squares over the seeds, and their ratios
    between resolutions. Growth seeded by noise should fall as the square root of the number of births (a ratio near
    0.5 for four times the births)."""
    out = {}
    for name, rs in (('lo', lo), ('hi', hi)):
        done = [r for r in rs if r['status'] == 'completed']
        out[name] = dict(n=rs[0]['spec']['n'] if rs else None, completed=len(done), runs=[record(r) for r in rs])
        if len(done) == len(rs) and done:
            nets = {f"{a['r_kpc']:g}": [r['apertures'][i]['net'] for r in done] for i, a in enumerate(done[0]['apertures'])}
            ex = [r['inner']['companions_within_25'] - r['inner']['depletion_within_25'] for r in done]
            out[name].update(apertures={a: dict(stats(v), rms=rms(v)) for a, v in nets.items()},
                             excess_within_25=dict(stats(ex), rms=rms(ex)), growth=growth(done))
    if 'apertures' in out['lo'] and 'apertures' in out['hi']:
        out['hi_over_lo_rms_net'] = {a: (out['hi']['apertures'][a]['rms']/out['lo']['apertures'][a]['rms']
                                         if out['lo']['apertures'][a]['rms'] else None) for a in out['lo']['apertures']}
        lo25, hi25 = out['lo']['excess_within_25']['rms'], out['hi']['excess_within_25']['rms']
        out['hi_over_lo_excess_within_25'] = hi25/lo25 if lo25 else None
    return out


def control_c3(c3, c3b, cont):
    """C3 (production stopped at T_OFF) and C3b (twice the rate until T_OFF: the same total conversion), against the
    continuous reference runs in the same region, at every snapshot; and the declared descriptive readings."""
    D = mw()
    hists = dict(continuous=cont, stopped=c3, early_matched=c3b)
    out = dict(t_off_Gyr=T_OFF, histories={}, readings=[])
    passes = {}
    for name, rs in hists.items():
        done = [r for r in rs if r['status'] == 'completed']
        if len(done) < SEEDS:
            out['histories'][name] = dict(reason='runs missing', statuses=[r['status'] for r in rs])
            continue
        snaps = {}
        for k, h in enumerate(done[0]['history']):
            prof = [r['history'][k]['profile'] for r in done]
            mean = {q: stats([p[q] for p in prof])['mean'] for q in ('rmse_38', 'rmse_inner20', 'slope_8_20')}
            g = gates(mean, D)
            snaps[f"{h['t_Gyr']:g}"] = dict(mean=mean, gates=g,
                                            net_by_aperture={f"{a['r_kpc']:g}": stats([r['history'][k]['apertures'][i]['net'] for r in done])
                                                             for i, a in enumerate(h['apertures'])})
            passes.setdefault(name, {})[round(h['t_Gyr'], 6)] = bool(g and g['profile_gate_passed'])
        out['histories'][name] = dict(snapshots=snaps, runs=[record(r) for r in rs] if name != 'continuous' else None)
    end = round(T_REF, 6)
    if passes.get('early_matched', {}).get(end):
        out['readings'].append('persistent reservoir: the matched early episode passes G1 and G2 at T')
    if passes.get('continuous', {}).get(end) and not passes.get('stopped', {}).get(end) and not passes.get('early_matched', {}).get(end):
        out['readings'].append('continuing infall: only the continuous history passes at T')
    for name, p in passes.items():
        if any(v for t, v in p.items() if t < end) and not p.get(end):
            out['readings'].append(f'transient match: the {name} history passes at an earlier snapshot but not at T')
    out['passes'] = {name: {f'{t:g}': v for t, v in p.items()} for name, p in passes.items()}
    return out


def best_sampled(c):
    """A 2B-F1 combination's best sampled rate and RMSE, over its completed ladder runs and its rescaling-root
    verification: the centre of the 2B-F1 revision's bracket (revision.py's targets)."""
    runs = [(r['q'], r['rmse_38']) for r in c['ladder'] if r['status'] == 'completed' and r.get('rmse_38') is not None]
    v = c.get('verify') or {}
    if v.get('status') == 'completed' and v.get('mw'):
        runs.append((c['q_star'], v['mw']['rmse_38']))
    return min(runs, key=lambda p: p[1])


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    D = mw()
    q0, rm = best_sampled(json.loads(SOURCE.read_text(encoding='utf-8'))['combos'][COMBO.format(V_REF)])
    log(f'{COMBO.format(V_REF)}: best sampled RMSE {rm:.2f} km/s at q = {q0:.5g} Msun/kpc^3/Gyr')
    cnt = [30000]

    def nxt():
        cnt[0] += 1
        return cnt[0]

    def spec(role, R, q=None, e=0, n=N_REF, **extra):
        return dict(role=role, key='MW', v_d=V_REF, q=q0*3**(e/8) if q is None else q, R_comp=float(R), n=n,
                    n_max=n_max_for(R, n), rng=nxt(), step=e, **extra)

    out = dict(scope='RC-2 part 1: a time-dependent donor-and-companion calculation in an open region (Amendment 1)',
               smoke=SMOKE, T_Gyr=T_REF, reference_decay_speed_kms=V_REF, regions_kpc=REGIONS, search_region_kpc=R_SEARCH,
               stage3_regions_kpc=STAGE3_REGIONS, apertures_kpc=APERTURES, seeds=SEEDS, tracers_reference=N_REF,
               tracers_high=N_HI, rate_source=SOURCE.name, fixed_rate=dict(combination=COMBO.format(V_REF), q=q0),
               tolerances=dict(FLOOR, speeds_kms=V_TOL, inventory_fraction=MASS_FLOOR),
               mass_unit_Msun=float(FM.system('MW')['M_b']),
               horizons_kpc={k: horizon(k) for k in ('MW', 'J1630', 'Coma low', 'Coma high')}, validation={}, controls={})
    # stage 1: D4, the reference over every region at the fixed rate, and the controls
    specs = [dict(role='D4', key='MW', v_d=V_REF, q=q0, R_comp=FM.system('MW')['R_b'], n=N_REF//4, rng=nxt())]
    specs += [spec('reference', R) for R in REGIONS for _ in range(SEEDS)]
    specs += [spec('C1a', R_C1, no_baryons=True) for _ in range(SEEDS)]
    specs += [spec('C1b', R_C1, n=n, no_baryons=True, quiet=False) for n in (N_REF, N_HI) for _ in range(SEEDS)]
    specs += [spec('C1c', R_C1C, n=N_HI) for _ in range(SEEDS)]
    specs += [spec('C3', R_C3, t_off_gyr=T_OFF) for _ in range(SEEDS)]
    specs += [spec('C3b', R_C3, q=2*q0, t_off_gyr=T_OFF) for _ in range(SEEDS)]
    res = run_tasks(specs, 'stage 1 (D4, the reference over every region, and controls C1a-C1c, C3, C3b)')
    role = lambda k: [r for r in res if r['spec']['role'] == k]
    ok_ref = lambda R: [r for r in role('reference') if r['spec']['R_comp'] == R and r['status'] == 'completed']
    out['validation']['D4'] = {k: x for k, x in role('D4')[0].items() if k not in ('spec', 'wall_s')}
    conv = convergence(role('reference'), REGIONS)
    conv['runs'] = [record(r) for r in role('reference')]
    out['stage1'] = dict(q=q0, convergence=conv)
    log(f"stage 1: demonstrated from R_comp = {conv['R_conv_demonstrated']}, statistically from {conv['R_conv_statistical']}; "
        + ', '.join(f"{x['R_kpc']:g}->{x['R2_kpc']:g} {x['demonstrated']}/{x['statistical']}" for x in conv['doublings']))
    out['controls']['C1a'] = control_c1a(role('C1a'), ok_ref(R_C1))
    out['controls']['C1b'] = control_c1b([r for r in role('C1b') if r['spec']['n'] == N_REF],
                                         [r for r in role('C1b') if r['spec']['n'] == N_HI])
    hi = [r for r in role('C1c') if r['status'] == 'completed']
    lo = ok_ref(R_C1C)
    c1c = compare(lo, hi) if len(hi) == SEEDS and len(lo) == SEEDS else dict(demonstrated=False, reason='runs missing')
    out['controls']['C1c'] = dict(c1c, runs=[record(r) for r in role('C1c')])
    out['controls']['C3'] = control_c3(role('C3'), role('C3b'), ok_ref(R_C3))
    log(f"C1a passed={out['controls']['C1a']['passed']} (inner ratio {out['controls']['C1a'].get('inner_ratio')}); "
        f"C1b hi/lo within 25 kpc {out['controls']['C1b'].get('hi_over_lo_excess_within_25')}; "
        f"C1c demonstrated={c1c['demonstrated']}; C3 readings {out['controls']['C3']['readings']}")
    interpretable = bool(out['controls']['C1a']['passed'] and c1c['demonstrated'])
    out['stage2_gate'] = dict(interpretable=interpretable, C1a=out['controls']['C1a']['passed'],
                              C1c_demonstrated=c1c['demonstrated'], forced_by_smoke=bool(SMOKE and not interpretable))
    all_runs = res[:]
    if interpretable or SMOKE:
        # stage 2: the direct search at the search region; the reference runs there are its step 0
        evals = {0: ok_ref(R_SEARCH)}

        def complete():
            return {e: float(np.mean([r['mw']['rmse_38'] for r in rs])) for e, rs in evals.items()
                    if len(rs) >= SEEDS and all(r['status'] == 'completed' for r in rs)}

        def absorb(rs):
            for r in rs:
                evals.setdefault(r['spec']['step'], []).append(r)
            all_runs.extend(rs)
        absorb(run_tasks([spec('search', R_SEARCH, e=e) for e in BRACKET if e != 0 for _ in range(SEEDS)],
                         'stage 2 (the bracket at the search region)'))
        for ext in range(MAX_EXT):
            mm = complete()
            e = min(mm, key=mm.get) if mm else None
            new = e + 2 if e is not None and e == max(evals) else e - 2 if e is not None and e == min(evals) else None
            if new is None:
                break
            absorb(run_tasks([spec('search', R_SEARCH, e=new) for _ in range(SEEDS)], f'stage 2.{ext + 1} (extension)'))
        mm = complete()
        if mm:
            e = min(mm, key=mm.get)
            absorb(run_tasks([spec('search', R_SEARCH, e=k) for k in (e - 1, e + 1) if k not in evals for _ in range(SEEDS)],
                             'stage 2 (refinement)'))
        mm = complete()
        rates = [dict(step=k, q=q0*3**(k/8), mean_rmse_38=mm.get(k), statuses=[r['status'] for r in evals[k]])
                 for k in sorted(evals)]
        if mm:
            e = min(mm, key=mm.get)
            S = region_summary(evals[e])
            mean = {k: S[k]['mean'] for k in ('rmse_38', 'rmse_inner20', 'slope_8_20')}
            q_sel = q0*3**(e/8)
            out['stage2'] = dict(q_sel=q_sel, step=e, mean=mean, gates=gates(mean, D),
                                 seed_gates=[r['gates'] for r in evals[e]], summary=S, rates=rates,
                                 runs=[record(r) for k in sorted(evals) if k != 0 for r in evals[k]])
            log(f"stage 2: q_sel = {q_sel:.5g} (step {e}), mean RMSE {mean['rmse_38']:.2f} km/s, gates {out['stage2']['gates']}")
            # stage 3: the region-size comparison at the selected rate, frozen; the search's runs at q_sel serve for R_SEARCH
            specs = [spec('stage 3', R, q=q_sel, e=e) for R in STAGE3_REGIONS if R != R_SEARCH for _ in range(SEEDS)]
            s3 = run_tasks(specs, 'stage 3 (the region-size comparison at the selected rate)')
            all_runs.extend(s3)
            pool3 = s3 + [dict(r, spec=dict(r['spec'], R_comp=R_SEARCH)) for r in evals[e]]
            conv3 = convergence(pool3, STAGE3_REGIONS)
            conv3['runs'] = [record(r) for r in s3]
            out['stage3'] = dict(q=q_sel, convergence=conv3)
            log(f"stage 3: demonstrated from R_comp = {conv3['R_conv_demonstrated']}, statistically from "
                f"{conv3['R_conv_statistical']}")
        else:
            out['stage2'] = dict(q_sel=None, rates=rates)
    # validation and control C2 over every run
    done = [r for r in all_runs if r['spec']['role'] != 'D4' and r['status'] == 'completed']
    V1 = dict(runs=len(done), mass_balance_max=max((r['D1']['mass_balance'] for r in done), default=math.inf),
              depletion_vs_production_max=max((r['D1']['depletion_vs_production'] for r in done), default=math.inf),
              draws_z_max=max((abs(r['D1']['draws_z']) for r in done), default=math.inf),
              energy_closure_max=max((r['D1']['energy_closure'] for r in done), default=math.inf),
              tags_follow_tracers=all(r['D1']['tags_follow_tracers'] for r in done))
    V1['passed'] = bool(V1['mass_balance_max'] < 1e-9 and V1['depletion_vs_production_max'] < 1e-9
                        and V1['draws_z_max'] < Z_DRAWS and V1['energy_closure_max'] < 1e-9 and V1['tags_follow_tracers'])
    out['validation']['D1'] = V1
    out['validation']['D2'] = dict(max_relative_error=max((r['D2'] for r in done), default=math.inf))
    out['validation']['D2']['passed'] = bool(out['validation']['D2']['max_relative_error'] < 1e-4)
    out['controls']['C2'] = dict(runs=len(done), max_abs_z=max((r['C2'] for r in done), default=math.inf), limit=Z_DRAWS)
    out['controls']['C2']['passed'] = bool(out['controls']['C2']['max_abs_z'] < Z_DRAWS)
    out['energy'] = {}
    for name, q in (('fixed_rate', q0), ('selected_rate', (out.get('stage2') or {}).get('q_sel'))):
        if q:
            es = F.energy_supply(dict(q=q, v_d=V_REF))
            es['parent_stock_over_microwave_background'] = {f'{g:g}': es['field_energy_over_microwave_background']/g
                                                            for g in (.01, .1, 1.)}
            es['parent_stock_over_cosmic_mean'] = {f'{g:g}': es['over_cosmic_mean']*(1 + es['eps'])/g for g in (.01, .1, 1.)}
            out['energy'][name] = es
    out['checks_short_run'] = short_run()
    # D1, D2, D4 and C2 check the code and its bookkeeping; C1 and C3 are physical controls, reported with their verdicts
    out['all_validation_passed'] = bool(all(x['passed'] for x in out['validation'].values()) and out['controls']['C2']['passed'])
    out['runtime_seconds'] = time.time() - T0
    log('validation: ' + ', '.join(f"{k}={x['passed']}" for k, x in out['validation'].items())
        + f", C2={out['controls']['C2']['passed']}")
    text = json.dumps(out, indent=1, default=float)
    if SMOKE:
        od = evidence_io.output_dir(args, 'donor-companion-smoke')
        (od/'rc2a-smoke.json').write_text(text, encoding='utf-8')
        log(f'smoke output {od}')
        return 0 if out['all_validation_passed'] else 1
    status = evidence_io.finish(args, 'donor-companion', text, HERE/'rc2a-results.json', ignore={'/runtime_seconds'},
                                rules=((r'.*/(wall_s|runtime_s)$', None, None),))
    return status if out['all_validation_passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
