"""CC-2 stage 2B-F1: a coherent field decaying into companion pairs, as a source and phase-space test (see protocol.md).

    python f1.py [--canonical] [--output-dir DIR]      (F1_SMOKE=1 runs a reduced set; F1_WORKERS sets the pool size)

For each decay speed, sigma/m and treatment of the bath's gravity, a ladder in the field's production rate q brackets
the normalization that best fits the Milky Way's rotation speeds. A verification run there is scored on the declared
profile gate (fit, shape and cost) and on energy supply. The best-scoring source then runs unchanged in J1630 and Coma,
and the Milky Way is repeated with a 150 kpc zone. Fails fast and logs its progress.
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
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
for _p in (RESULTS/'gravitational-focusing', RESULTS/'companion-supply', RESULTS/'companion-formation',
           RESULTS/'capture-to-orbit', HERE):
    sys.path.insert(0, str(_p))
import cf1  # noqa: E402  (also puts the shared result folders on the path)
import mc  # noqa: E402
import formation as FM  # noqa: E402
import field as FD  # noqa: E402
import inputs as I  # noqa: E402
import evidence_io  # noqa: E402

SMOKE = os.environ.get('F1_SMOKE') == '1'
V_D = (10., 300.) if SMOKE else (3., 10., 30., 100., 300.)         # km/s, the decay speeds (universal)
SIGMA_M = (0., 1.) if SMOKE else (0., .1, 1.)                      # cm^2/g
GRAVITY = (False,) if SMOKE else (True, False)                     # the bath's focused excess gravitates, or not
T_REF = 2. if SMOKE else 10.                                       # Gyr
SNAPS = tuple(t for t in (0., .25, .5, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.) if t <= T_REF)
N_TARGET, N_FIELD, N_MAX = (3000, 3000, 10000) if SMOKE else (16000, 16000, 40000)
DELTA_MAX = {'galaxy': .005, 'lens host': .005, 'cluster': .05}   # Gyr, as stage 2A
BUDGET = float(os.environ.get('F1_BUDGET', '900' if SMOKE else '2400'))
WORKERS = int(os.environ.get('F1_WORKERS', '4' if SMOKE else '16'))
LADDER = (1/27, 1/9, 1/3, 1., 3.)                                  # multiples of the starting estimate q0
STEP, MAX_EXTEND = 3., 3
G1_RMSE, G2_TOL, G3_COST, G3_LENIENT = 20., .3, 10., 20.
KPC_M, MSUN_KG, GYR_S = 3.0856775814913673e19, 1.98847e30, 3.15576e16
A_RAD, C_MS, T_CMB = 7.565723e-16, 299792458., 2.72548
RHO_CMB = A_RAD*T_CMB**4/C_MS**2*KPC_M**3/MSUN_KG                  # the microwave background's energy as mass, Msun/kpc^3
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:9.1f}s] {msg}', flush=True)


# ---- the Milky Way target (CR-2's rule and its interior diagnostic)
def mw_data():
    base = next(r for r in I.milky_way_runs() if r['baryons'] == 'I' and abs(r['rd'] - 2.6) < 1e-9 and abs(r['lf'] - 1) < 1e-9)
    R, y, vb = (np.asarray(base[k], float) for k in ('R', 'y', 'vb'))
    M_req = R*(y**2 - vb**2)/mc.G
    sel = (R >= 8) & (R <= 20)
    return dict(R=R, y=y, vb=vb, M_req=M_req, sel=sel, slope_req=float(np.polyfit(np.log(R[sel]), np.log(M_req[sel]), 1)[0]),
                baryons_rmse_38=float(np.sqrt(np.mean((vb - y)**2))), baryons_rmse_inner20=float(np.sqrt(np.mean((vb[:20] - y[:20])**2))))


def mw_score(model, D):
    """Rotation fit, enclosed-mass slope and cost of one run: all gravitating non-baryonic mass counts."""
    rr = np.linalg.norm(model.x, axis=1) if len(model.x) else np.zeros(0)
    o = np.argsort(rr)
    rs, cm = rr[o], np.r_[0., np.cumsum(model.m[o])]
    Mc = cm[np.searchsorted(rs, D['R'], side='right')]
    Mbx = np.interp(np.log(D['R']), model.lr, model.M_bx) if model.bath_gravity else np.zeros(len(D['R']))
    M = Mc + Mbx
    v = np.sqrt(np.maximum(D['vb']**2 + mc.G*M/D['R'], 0.))
    rmse = lambda n=None: float(np.sqrt(np.mean((v[:n] - D['y'][:n])**2)))
    sel = D['sel']
    slope = float(np.polyfit(np.log(D['R'][sel]), np.log(M[sel]), 1)[0]) if (M[sel] > 0).all() else None
    Mb = float(model.M_b[-1])
    bath = float(model.M_bx[-1]) if model.bath_gravity else 0.
    return dict(rmse_38=rmse(), rmse_inner20=rmse(20), slope_8_20=slope, cost=(float(model.m.sum()) + bath)/Mb,
                confined_over_baryons=float(model.m.sum())/Mb, bath_excess_in_R_b_over_baryons=bath/Mb,
                extra_mass_over_required={f'{r0:g} kpc': float(np.interp(r0, D['R'], M)/np.interp(r0, D['R'], D['M_req']))
                                          for r0 in (6., 10., 15., 20., 24.)},
                predicted_kms=v.tolist(), M_extra_at_R=M.tolist())


def rescale(sc, D):
    """The factor by which a run's extra mass would have to be scaled, its shape held fixed, to minimize the RMSE (the
    way CR-2 fits a profile's mass). Near 1 the run sits at its best normalization."""
    M = np.asarray(sc['M_extra_at_R'])
    if not (M > 0).any():
        return None
    f = lambda ls: float(np.sqrt(np.mean((np.sqrt(np.maximum(D['vb']**2 + math.exp(ls)*mc.G*M/D['R'], 0.)) - D['y'])**2)))
    return float(math.exp(minimize_scalar(f, bounds=(-6., 6.), method='bounded', options={'xatol': 1e-6}).x))


def gates(sc, D):
    if sc is None:
        return None
    g1 = bool(sc['rmse_38'] <= G1_RMSE and sc['rmse_inner20'] <= D['baryons_rmse_inner20'])
    g2 = bool(sc['slope_8_20'] is not None and abs(sc['slope_8_20'] - D['slope_req']) <= G2_TOL)
    g3 = bool(sc['cost'] <= G3_COST)
    return dict(G1_fit=g1, G2_shape=g2, G3_cost=g3, lenient_cost=bool(sc['cost'] <= G3_LENIENT), profile_gate_passed=g1 and g2 and g3)


def energy_supply(spec):
    """What the field must supply: the companion mass made by T, the field's minimum density, its energy-loss rate."""
    qT = spec['q']*T_REF
    eps = FD.eps_for_speed(spec['v_d'])
    q_si = spec['q']*MSUN_KG/KPC_M**3/GYR_S                        # kg m^-3 s^-1 of companion rest mass
    return dict(eps=eps, per_field_energy=FD.per_field_energy(eps), companion_mass_made_by_T_Msun_kpc3=qT,
                over_cosmic_mean=qT/cf1.RHO_MEAN, field_min_density_Msun_kpc3=qT*(1 + eps),
                field_energy_over_microwave_background=qT*(1 + eps)/RHO_CMB,
                field_energy_loss_W_m3=q_si*(1 + eps)*C_MS**2)


# ---- systems
def system_for(spec):
    s = dict(FM.system(spec['key']))
    if spec.get('R_b'):
        s['R_b'] = float(spec['R_b'])
        s['M_b'] = float(np.interp(s['R_b'], s['r'], s['M']))
    return s


def make_model(spec, s, freeze=False, channels=None):
    dmax = DELTA_MAX[s['kind']]*spec.get('dmax_factor', 1.)
    sm = spec['sm']
    if channels is None:
        channels = ('ii', 'ib', 'bb') if sm > 0 else ()
    # spec may double the potential grid (n_grid) or scale the field's pools (pool_factor): the revision's controls
    model = FD.FieldModel(s['r'], s['M'], s['R_b'], s['r_half'], np.array([spec['v_d']]), np.array([1.]), sm,
                          spec['q']*dmax/2, q=spec['q']/mc.PER_GYR, v_d=spec['v_d'], self_gravity=not freeze,
                          freeze=freeze, bath_gravity=spec['gravity'], channels=channels, seed=spec['rng'],
                          **({'n_grid': int(spec['n_grid'])} if spec.get('n_grid') else {}))
    model.pool_factor = float(spec.get('pool_factor', 1.))
    return model, dmax


def summarize_f1(model, res, s, spec, D):
    out = FM.summarize(model, res, s, spec.get('full', False), 0.)
    L = res['ledger']
    for k in ('field_births', 'field_born_mass', 'E_field_births', 'field_exported_mass', 'E_field_exported',
              'birth_mass_doublings', 'collisional_mass_fraction_max'):
        out['ledger'][k] = float(L.get(k, 0.))
    # regime diagnostics for a bath at speed v_d growing as q t (stage 2A's assume a fixed 300 km/s bath)
    qT = spec['q']*T_REF
    sig = spec['v_d']/math.sqrt(3)
    Gr = mc.G*qT
    inflow = math.pi*model.R_b**2*spec['v_d']*spec['q']/mc.PER_GYR*(T_REF*mc.PER_GYR)**2/2
    out['regime'].update(incident_density_at_T_Msun_kpc3=qT, inflow_mass_over_span=inflow,
                         capture_efficiency=float(out['M_T'])/inflow if inflow > 0 else None,
                         jeans_length_kpc=sig*math.sqrt(math.pi/Gr), jeans_length_over_R_b=sig*math.sqrt(math.pi/Gr)/model.R_b,
                         kJ_R_b_over_pi=math.sqrt(4*math.pi*Gr)/sig*model.R_b/math.pi,
                         jeans_growth_time_Gyr=1/math.sqrt(4*math.pi*Gr)/mc.PER_GYR)
    out['energy'] = energy_supply(spec)
    if spec['key'] == 'MW':
        out['mw'] = mw_score(model, D)
        out['gates'] = gates(out['mw'], D)
    return out


# ---- tasks
def task(spec):
    try:
        return _task(spec)
    except Exception as e:                       # keep the pool alive; the driver reports and stops
        return dict(spec=spec, status=f'error: {type(e).__name__}: {e}', traceback=traceback.format_exc())


def _task(spec):
    t0 = time.time()
    if spec['role'] == 'validation':
        return dict(spec=spec, status='completed', **VALIDATIONS[spec['name']](spec), wall_s=time.time() - t0)
    D = mw_data()
    s = system_for(spec)
    model, dmax = make_model(spec, s)
    out = dict(spec=spec, bath_static=bool(model.bath_static), static_fail_reason=model.static_fail_reason)
    if not model.bath_static:
        out['status'] = 'no static bath'
        return out
    nf = spec.get('n_factor', 1.)
    try:
        res = model.run(T_REF, dmax, SNAPS, n_target=int(N_TARGET*nf), n_field=int(N_FIELD*nf), n_max=int(N_MAX*nf),
                        budget_s=BUDGET, regen=int(spec.get('regen', 10)))
    except TimeoutError as e:
        out['status'] = f'not completed: {e}'
        return out
    except mc.RunawayError as e:
        out['status'] = f'not completed: runaway, {e}'
        return out
    except mc.StaticBathLost as e:
        out['status'] = f'not completed: static bath lost, {e}'
        return out
    except FD.OpaqueBath as e:
        out['status'] = f'not completed: opaque, {e}'
        return out
    out['status'] = 'completed'
    out.update(summarize_f1(model, res, s, spec, D))
    out['wall_s'] = time.time() - t0
    return out


# ---- validation tasks (F1-F4)
def v_f1(spec):
    rng = np.random.default_rng(spec['rng'])
    rows = {}
    for v in V_D:
        eps = FD.eps_for_speed(v)
        E1, p1, E2, p2 = FD.sample_decays(spec['n'], eps, rng)
        mu = p1[:, 2]/np.linalg.norm(p1, axis=1)
        rows[f'{v:g}'] = dict(eps=eps, energy=float(np.max(np.abs(E1 + E2 - 2*(1 + eps)))/(2*(1 + eps))),
                              momentum=float(np.max(np.abs(p1 + p2))), mass_shell=float(np.max(np.abs(E1**2 - (p1**2).sum(1) - 1))),
                              speed_relative=float(np.mean(np.linalg.norm(p1, axis=1)/E1))*FD.C_KMS/v - 1,
                              mean_mu_over_se=float(mu.mean()*math.sqrt(3*len(mu))),
                              mean_mu2_minus_third_over_se=float(((mu**2).mean() - 1/3)/(2/math.sqrt(45*len(mu)))))
    ok = all(r['energy'] < 1e-15 and r['momentum'] < 1e-15 and r['mass_shell'] < 1e-14 and abs(r['speed_relative']) < 1e-12
             and abs(r['mean_mu_over_se']) < 3 and abs(r['mean_mu2_minus_third_over_se']) < 3 for r in rows.values())
    return dict(rows=rows, passed=bool(ok))


def se_mass(m):
    return float(np.sqrt(np.sum(m*m)))


def v_f2(spec):
    s = system_for(spec)
    rows = {}
    for v in (100., 300.):
        sp = dict(spec, v_d=v)
        model, dmax = make_model(sp, s, freeze=True, channels=())
        res = model.run(T_REF, dmax, [T_REF], n_field=spec['n'], budget_s=BUDGET)
        M, err = float(model.m.sum()), se_mass(model.m)
        pred = FD.born_bound_mass(model, model.q, T_REF*mc.PER_GYR, v)
        rows[f'{v:g}'] = dict(engine=M, standard_error=err, quadrature=pred, z=(M - pred)/err, closure=res['ledger']['closure_max'])
    return dict(rows=rows, passed=bool(all(abs(r['z']) < 3 for r in rows.values())))


def v_f3(spec):
    """Cold births (v_d = 1 km/s) in the frozen Milky Way potential against the radial-orbit quadrature, in ten 2 kpc
    bins over 5-25 kpc. A statistical check: a bin's tolerance is three standard errors or 5%, whichever is larger, and
    the 5% term controls only where the relative standard error is below 1.67%. spec['n_max'] caps the tracers; without
    it the engine's default of 40,000 applies and thins a larger population, as it did in 2B-F1's run. The thinnings,
    birth-mass doublings and each bin's effective sample size, (sum w)^2 / sum w^2, are reported."""
    s = system_for(spec)
    sp = dict(spec, v_d=1.)
    model, dmax = make_model(sp, s, freeze=True, channels=())
    kw = dict(n_max=int(spec['n_max'])) if spec.get('n_max') else {}
    L = model.run(T_REF, dmax, [T_REF], n_field=spec['n'], budget_s=spec.get('budget', BUDGET), **kw)['ledger']
    rr = np.linalg.norm(model.x, axis=1)
    edges = np.linspace(5., 25., 11)
    rows, worst = [], 0.
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (rr >= a) & (rr < b)
        m = model.m[sel]
        V = 4/3*math.pi*(b**3 - a**3)
        rho_e, err = float(m.sum())/V, se_mass(m)/V
        sub = np.linspace(a, b, 9)
        pr = FD.cold_density(model, model.q, T_REF*mc.PER_GYR, sub)
        rho_p = float(np.trapezoid(pr*sub**2, sub)/np.trapezoid(sub**2, sub))
        tol = max(3*err, .05*rho_p)
        worst = max(worst, abs(rho_e - rho_p)/tol)
        rows.append(dict(r_kpc=(float(a), float(b)), engine=rho_e, standard_error=err, quadrature=rho_p, tracers=int(sel.sum()),
                         n_eff=float(m.sum()**2/(m @ m)) if len(m) else 0.))
    return dict(rows=rows, worst_over_tolerance=worst, passed=bool(worst < 1), n_max=kw.get('n_max', 40000),
                thinnings=int(L['thinnings']), birth_mass_doublings=int(L['birth_mass_doublings']), tracers=int(len(model.m)))


def v_f4(spec):
    s = system_for(spec)
    sp = dict(spec, v_d=300., sm=1e-3)
    model, dmax = make_model(sp, s, freeze=True, channels=('ii',))
    res = model.run(T_REF, dmax, [T_REF], n_target=spec['n'], n_field=2000, budget_s=BUDGET)
    L = res['ledger']
    # the prediction from the pools the run used, and the births' own variance (FieldModel._book_seedless)
    born, pred, n = L['born_mass'], L['seedless_predicted'], L['births']
    err = math.sqrt(max(L['seedless_birth_variance'], 1e-300))
    model.make_pools()                                   # for the record: one set of pools drawn after the run
    one = FD.growing_bath_seedless(model, model.q, T_REF*mc.PER_GYR)
    return dict(engine=born, standard_error=err, prediction=pred, z=(born - pred)/err, births=n,
                prediction_from_one_pool_after_run=one,
                tau_r_half=float(np.interp(math.log(s['r_half']), model.lr, model.tau)), passed=bool(abs((born - pred)/err) < 3))


VALIDATIONS = dict(F1=v_f1, F2=v_f2, F3=v_f3, F4=v_f4)


def short_run():
    """The suite's regression anchor: Milky Way, v_d = 10 km/s, 0.1 cm^2/g, q = 3e3 Msun/kpc^3/Gyr, 0.5 Gyr."""
    spec = dict(key='MW', v_d=10., sm=.1, q=3e3, gravity=False, rng=4243, role='anchor')
    s = system_for(spec)
    model, dmax = make_model(spec, s)
    res = model.run(.5, dmax, [.5], n_target=2000, n_field=2000, budget_s=600)
    L = res['ledger']
    return dict(M=float(model.m.sum()), tracers=int(len(model.m)), field_births=float(L.get('field_births', 0.)),
                births=float(L['births']), captures=float(L['captures']), closure=float(L['closure_max']))


# ---- the driver
def run_tasks(specs, label):
    log(f'{label}: {len(specs)} runs on {WORKERS} workers')
    out = []
    with Pool(min(WORKERS, max(len(specs), 1))) as pool:
        for i, r in enumerate(pool.imap_unordered(task, specs), 1):
            sp = r['spec']
            extra = ''
            if r['status'] == 'completed' and r.get('mw'):
                w = r['mw']
                extra = f"RMSE {w['rmse_38']:.1f} slope {w['slope_8_20'] if w['slope_8_20'] is None else round(w['slope_8_20'], 2)} cost {w['cost']:.3g}"
            elif r['status'] == 'completed' and sp['role'] == 'validation':
                extra = f"passed={r['passed']}"
            log(f"  [{i}/{len(specs)}] {sp['role']} {sp.get('key', '')} v_d={sp.get('v_d', '')} sm={sp.get('sm', '')} "
                f"q={sp.get('q', 0):.3g} gravity={sp.get('gravity', '')}: {r['status'][:90]} {extra}")
            if r['status'].startswith('error'):
                log('ERROR TRACEBACK\n' + r.get('traceback', ''))
                raise RuntimeError(f"task failed: {r['status']}")
            out.append(r)
    # imap_unordered hands results back as they finish; restore the order the tasks were issued in (each spec's seed is
    # unique and increasing), so that nothing downstream, the next round's seeds in particular, depends on timing
    out.sort(key=lambda r: r['spec']['rng'])
    return out


def q_start(key, v_d, D):
    """The starting estimate: the rate whose born-bound mass alone reaches ten times the baryons (the cost scale) in
    the frozen baryons-only potential (F2's quadrature on its plain grid)."""
    s = system_for(dict(key=key))
    model = FD.FieldModel(s['r'], s['M'], s['R_b'], s['r_half'], np.array([v_d]), np.array([1.]), 0., 1., q=0., v_d=v_d,
                          bath_gravity=False, freeze=True, channels=(), seed=1)
    # kpc^3 per unit rate and time. A ladder's start needs only a rough value, so it keeps the plain 600-point grid of
    # the first canonical run (n_sub=0), and a rerun reproduces that run's ladders exactly
    V = FD.born_bound_mass(model, 1., 1., v_d, n_sub=0)
    return 10*s['M_b']/(T_REF*mc.PER_GYR*V)*mc.PER_GYR if V > 0 else None


def next_q(runs, D):
    """The rate at which the extra mass reaches its best normalization: the root of ln s(q) = 0, where s is the factor by
    which a completed run's extra mass (shape held fixed) would have to be scaled to minimize the RMSE. s falls as q
    rises, but not in proportion (self-gravity concentrates more mass as q grows), so the root is interpolated in log q
    between the runs that bracket it, as stage 2A interpolates its benchmarks. Without a bracket, one step from the
    nearest run, capped at a factor STEP and kept below any stopped run at a higher rate (the runaway threshold)."""
    pts = sorted([(r['spec']['q'], rescale(r['mw'], D)) for r in runs if r['status'] == 'completed' and r.get('mw')],
                 key=lambda p: p[0])
    pts = [(q, s) for q, s in pts if s]
    if not pts:
        return None, 'none completed'
    for (qa, sa), (qb, sb) in zip(pts[:-1], pts[1:]):
        if sa >= 1 >= sb:
            if sa == sb:
                return float(math.sqrt(qa*qb)), 'bracketed'
            return float(math.exp(math.log(qa) + math.log(sa)*(math.log(qb) - math.log(qa))/(math.log(sa) - math.log(sb)))), 'bracketed'
    if pts[-1][1] > 1:                        # every completed run needs more mass
        q, s = pts[-1]
        qn, where = q*min(s, STEP), 'high edge'
    else:                                     # every completed run has too much
        q, s = pts[0]
        qn, where = q*max(s, 1/STEP), 'low edge'
    above = [r['spec']['q'] for r in runs if r['status'] != 'completed' and r['spec']['q'] > q]
    if above:
        qn = min(qn, math.sqrt(q*min(above)))
    return float(qn), where


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    D = mw_data()
    out = dict(scope='CC-2 stage 2B-F1: a coherent field decaying into companion pairs, as a source and phase-space test',
               smoke=SMOKE, T_Gyr=T_REF, decay_speeds_kms=V_D, sigma_over_m_cm2g=SIGMA_M,
               mw_target=dict(required_slope_8_20=D['slope_req'], baryons_rmse_38=D['baryons_rmse_38'],
                              baryons_rmse_inner20=D['baryons_rmse_inner20'], gates=dict(G1_rmse=G1_RMSE, G2_tolerance=G2_TOL,
                                                                                      G3_cost=G3_COST, lenient_cost=G3_LENIENT)),
               validation={}, combos={})
    rng = [2000]

    def nxt():
        rng[0] += 1
        return rng[0]
    nv = (1, .2) if SMOKE else (1, 1)
    vspecs = [dict(role='validation', name='F1', n=int(1e6*nv[1]), rng=nxt()),
              dict(role='validation', name='F2', key='MW', sm=0., q=1e3, gravity=False, n=int(32000*nv[1]), rng=nxt()),
              dict(role='validation', name='F3', key='MW', sm=0., q=1e3, gravity=False, n=int(64000*nv[1]), rng=nxt()),
              dict(role='validation', name='F4', key='MW', q=1e4, gravity=False, n=int(64000*nv[1]), rng=nxt())]
    q0 = {v: q_start('MW', v, D) for v in V_D}
    log('starting rates (Msun/kpc^3/Gyr): ' + ', '.join(f'{v:g} km/s: {q:.3g}' for v, q in q0.items()))
    ladders = {}
    specs = list(vspecs)
    for v in V_D:
        for sm in SIGMA_M:
            for grav in GRAVITY:
                ladders[(v, sm, grav)] = []
                specs += [dict(role='ladder', key='MW', v_d=v, sm=sm, gravity=grav, q=q0[v]*f, rng=nxt()) for f in LADDER]
    res = run_tasks(specs, 'round 1')
    for r in res:
        if r['spec']['role'] == 'validation':
            out['validation'][r['spec']['name']] = {k: v for k, v in r.items() if k not in ('spec', 'status', 'wall_s')}
        else:
            sp = r['spec']
            ladders[(sp['v_d'], sp['sm'], sp['gravity'])].append(r)
    log('validation: ' + ', '.join(f"{k}={v['passed']}" for k, v in out['validation'].items()))
    # extend ladders whose best point is at an edge (at most MAX_EXTEND times)
    for ext in range(MAX_EXTEND):
        specs = []
        for combo, runs in ladders.items():
            qn, where = next_q(runs, D)
            if where in ('low edge', 'high edge') and qn and not any(abs(math.log(qn/r['spec']['q'])) < 1e-3 for r in runs):
                specs.append(dict(role='ladder', key='MW', v_d=combo[0], sm=combo[1], gravity=combo[2], q=qn, rng=nxt()))
        if not specs:
            break
        for r in run_tasks(specs, f'round 2.{ext + 1}'):
            sp = r['spec']
            ladders[(sp['v_d'], sp['sm'], sp['gravity'])].append(r)
    # verification at the interpolated rate, then refinement: every verification joins the runs the root is interpolated
    # from, and the verification kept is the one whose own best normalization lies closest to it (stopping once within
    # 10%, at most MAX_EXTEND refinements). A refinement that stops (runaway, opacity) caps the next try below it.
    where_, verify, refined = {}, {}, {}
    pending = {}
    for combo, runs in ladders.items():
        q, where = next_q(runs, D)
        where_[combo] = where
        if q:
            pending[combo] = q
    dist = lambda r: abs(math.log(rescale(r['mw'], D) or 1e-300)) if r['status'] == 'completed' and r.get('mw') else math.inf
    for ref in range(MAX_EXTEND + 1):
        specs = [dict(role='verify', key='MW', v_d=c[0], sm=c[1], gravity=c[2], q=q, full=True, rng=nxt()) for c, q in pending.items()]
        if not specs:
            break
        pending = {}
        for r in run_tasks(specs, 'round 3 (verification)' if ref == 0 else f'round 3.{ref} (refinement)'):
            sp = r['spec']
            combo = (sp['v_d'], sp['sm'], sp['gravity'])
            refined[combo] = ref
            if combo not in verify or dist(r) < dist(verify[combo]):
                if combo in verify:
                    ladders[combo].append(verify[combo])
                verify[combo] = r
            else:
                ladders[combo].append(r)
            if ref < MAX_EXTEND and dist(verify[combo]) > math.log(1.1):
                q, _ = next_q(ladders[combo] + [verify[combo]], D)
                if q and not any(abs(math.log(q/x['spec']['q'])) < 1e-3 for x in ladders[combo] + [verify[combo]]):
                    pending[combo] = q
    # aggregate
    for combo, runs in ladders.items():
        v = verify.get(combo)
        name = f"v_d {combo[0]:g} km/s | {combo[1]:g} cm2/g | bath gravity {'on' if combo[2] else 'omitted'}"
        out['combos'][name] = dict(
            v_d_kms=combo[0], sigma_over_m_cm2g=combo[1], bath_gravity=combo[2], q0=q0[combo[0]], best_point=where_[combo],
            refinements=refined.get(combo, 0),
            rescale_at_verification=(rescale(v['mw'], D) if v and v.get('mw') else None),
            ladder=[dict(q=r['spec']['q'], status=r['status'], **({k: r['mw'][k] for k in ('rmse_38', 'rmse_inner20', 'slope_8_20', 'cost')}
                                                                  if r.get('mw') else {})) for r in sorted(runs, key=lambda r: r['spec']['q'])],
            verify=({k: x for k, x in v.items() if k != 'spec'} if v else None), q_star=(v['spec']['q'] if v else None))
    # the best-scoring source: the full gate, then G1 and G2, then G1, then the lowest RMSE; lowest cost within a tier
    cands = [(name, c) for name, c in out['combos'].items() if c['verify'] and c['verify'].get('status') == 'completed']

    def tier(c):
        g = c['verify']['gates']
        return 0 if g['profile_gate_passed'] else 1 if (g['G1_fit'] and g['G2_shape']) else 2 if g['G1_fit'] else 3
    best = min(cands, key=lambda nc: (tier(nc[1]), nc[1]['verify']['mw']['cost'] if tier(nc[1]) < 3 else nc[1]['verify']['mw']['rmse_38'])) if cands else None
    out['best'] = dict(name=best[0], tier=tier(best[1])) if best else None
    if best:
        c = best[1]
        base = dict(v_d=c['v_d_kms'], sm=c['sigma_over_m_cm2g'], gravity=c['bath_gravity'], q=c['q_star'], full=True)
        specs = [dict(base, role='universality', key=k, rng=nxt()) for k in ('J1630', 'Coma low', 'Coma high')]
        qs150 = [c['q_star']*f*8 for f in LADDER]                  # a 150 kpc zone holds an eighth of the volume
        specs += [dict(base, role='sensitivity ladder', key='MW', R_b=150., q=q, rng=nxt(), full=False) for q in qs150]
        res = run_tasks(specs, 'round 4 (universality and the 150 kpc zone)')
        out['universality'] = {r['spec']['key']: {k: x for k, x in r.items() if k != 'spec'} for r in res if r['spec']['role'] == 'universality'}
        sens = [r for r in res if r['spec']['role'] == 'sensitivity ladder']
        q150, where150 = next_q(sens, D)
        sv = run_tasks([dict(base, role='sensitivity verify', key='MW', R_b=150., q=q150, rng=nxt())], 'round 5 (150 kpc verification)')[0] if q150 else None
        out['sensitivity_150kpc'] = dict(ladder=[dict(q=r['spec']['q'], status=r['status'], **({k: r['mw'][k] for k in ('rmse_38', 'slope_8_20', 'cost')} if r.get('mw') else {}))
                                                 for r in sorted(sens, key=lambda r: r['spec']['q'])],
                                         best_point=where150, verify=({k: x for k, x in sv.items() if k != 'spec'} if sv else None))
    # validation across every run: collisions (V5) and ledgers (V6)
    allruns = [r for runs in ladders.values() for r in runs] + list(verify.values())
    done = [r for r in allruns if r['status'] == 'completed']
    out['validation']['V5'] = dict(max_kinetic=max((r['ledger']['max_dK'] for r in done), default=0.),
                                   max_momentum=max((r['ledger']['max_dP'] for r in done), default=0.))
    out['validation']['V5']['passed'] = bool(out['validation']['V5']['max_kinetic'] < 1e-12 and out['validation']['V5']['max_momentum'] < 1e-12)
    out['validation']['V6'] = dict(energy_closure_max=max((r['ledger']['closure_max'] for r in done), default=0.),
                                   mass_ledger_max=max((r['mass_ledger_relative'] for r in done), default=0.))
    out['validation']['V6']['passed'] = bool(out['validation']['V6']['energy_closure_max'] < 1e-9 and out['validation']['V6']['mass_ledger_max'] < 1e-9)
    out['checks_short_run'] = short_run()
    out['all_validation_passed'] = all(v['passed'] for v in out['validation'].values())
    out['runtime_seconds'] = time.time() - T0
    log('validation: ' + ', '.join(f"{k}={v['passed']}" for k, v in out['validation'].items()))
    text = json.dumps(out, indent=1, default=float)
    if SMOKE:
        d = evidence_io.output_dir(args, 'companion-source-smoke')
        (d/'f1-smoke.json').write_text(text, encoding='utf-8')
        log(f'smoke output {d}')
        return 0
    status = evidence_io.finish(args, 'companion-source', text, HERE/'f1-results.json', ignore={'/runtime_seconds'},
                                rules=((r'.*/(wall_s|runtime_s)$', None, None),))
    return status if out['all_validation_passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
