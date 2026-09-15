"""CC-2 stage 2A: complete-channel formation and supply revision (see protocol.md).

    python formation.py [--canonical] [--output-dir DIR]        (CC2A_SMOKE=1 runs a reduced set)

The incident density needed to form a given bound population within the benchmark time is recomputed with all three
collision classes of CF-1's law acting together, and with mass, energy, orbits, bath depletion and gravity evolving
together (mc.py). The fixed-background control C1 gives the starting estimate. Tracer runs bracket the requirement
and interpolate it for benchmarks B1 and B2, and controls at the requirement test the seed and self-gravity.
The benchmarks, the 10 Gyr span and the incident distributions are labeled trials; the cosmic mean density is a
comparison unit only. Fails fast and logs its progress.
"""
import json
import math
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
for _p in (RESULTS/'gravitational-focusing', RESULTS/'companion-supply', HERE):
    sys.path.insert(0, str(_p))
import cf1  # noqa: E402  (also puts the shared result folders on the path)
import focus as F  # noqa: E402
import cc2a  # noqa: E402
import mc  # noqa: E402
import evidence_io  # noqa: E402

SMOKE = os.environ.get('CC2A_SMOKE') == '1'
SIGMA_M = (1., 1000.) if SMOKE else (.1, 1., 10., 100., 1000.)  # cm^2/g, trial scales
U = 300.                                                        # km/s, CF-1's trial speed
BATHS = ('S',) if SMOKE else ('S', 'M')
T_REF = 2. if SMOKE else 10.                                    # Gyr, labeled benchmark
SNAPS = tuple(s for s in (0., .25, .5, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.) if s <= T_REF)
N_TARGET, N_SEED, N_MAX = (3000, 2000, 10000) if SMOKE else (16000, 10000, 40000)
N_V4, R_V4, NK_V4 = (20000, 20, 60000) if SMOKE else (200000, 100, 400000)
DELTA_MAX = {'galaxy': .005, 'lens host': .005, 'cluster': .05}     # Gyr
BUDGET = float(os.environ.get('CC2A_BUDGET', '900' if SMOKE else '2400'))   # s per tracer run
WORKERS = int(os.environ.get('CC2A_WORKERS', '4' if SMOKE else '20'))
TEST_RADII = {'MW': (5., 10., 25.), 'J1630': (3., 10., 21.), 'Coma low': (300., 1000., 2000., 3000.),
              'Coma high': (300., 1000., 2000., 3000.)}
SLOPE_RANGE = {'MW': (5., 25.), 'J1630': (3., 21.), 'Coma low': (300., 3000.), 'Coma high': (300., 3000.)}
MAX_EXTEND = 3
LEDGER_KEYS = ('births', 'born_mass', 'captures', 'captured_mass', 'ejections', 'ejected_mass', 'collisions_ib',
               'collisions_bb', 'evaporations', 'evaporated_mass', 'orbit_escapes', 'escaped_mass', 'E_births', 'E_bath_in',
               'E_exported', 'E_potential_work', 'E_integration', 'E_thinning', 'E_resample', 'M_thinning', 'M_resample',
               'closure_max', 'max_dK', 'max_dP', 'P_max', 'P_over_1', 'thinnings', 'birth_mass_doublings', 'splits',
               'ib_substeps_max', 'bb_level_max', 'capped_mass_fraction_max', 'collisional_mass_fraction_max',
               'bath_max_change', 'steps')
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:8.1f}s] {msg}', flush=True)


def bath(name):
    return mc.single_speed(U) if name == 'S' else mc.maxwell_components(U/math.sqrt(3))


def key_of(s):
    if s['kind'] == 'galaxy':
        return 'MW'
    if 'J1630' in s['name']:
        return 'J1630'
    return 'Coma low' if 'n_e0=0.0025' in s['name'] else 'Coma high'


def load_systems():
    out = {}
    for s in cf1.systems():
        if not (s['kind'] in ('galaxy', 'cluster') or 'J1630' in s['name']):
            continue
        pot = F.Potential(s['r'], s['M'], s['R_b'], s['name'])
        Mb = float(pot.M[-1])
        rho, sj = cf1.seed_profile(pot)
        k = key_of(s)
        out[k] = dict(key=k, name=s['name'], kind=s['kind'], r=s['r'], M=s['M'], R_b=float(s['R_b']), pot=pot, M_b=Mb,
                      r_half=float(np.interp(.5*Mb, pot.M, pot.r)), seed_rho=rho, seed_sigma=sj)
    return out


_SYS = {}


def system(k):
    if not _SYS:
        _SYS.update(load_systems())
    return _SYS[k]


def model_for(s, bath_name, sm, rho, **kw):
    u, f = bath(bath_name)
    return mc.Model(s['r'], s['M'], s['R_b'], s['r_half'], u, f, sm, rho, **kw)


# ----------------------------------------------------------------------------------------------------------------
# summaries of one tracer run
def projected(radii, m, R):
    """Projected mass within cylinder radius R for tracers at radii with masses m."""
    w = np.where(radii <= R, 1., 1 - np.sqrt(np.maximum(0., 1 - (R/np.maximum(radii, 1e-300))**2)))
    return float(m @ w)


def projected_profile(model, M_grid, R):
    r = model.r
    dM = np.diff(M_grid)
    rm = .5*(r[1:] + r[:-1])
    w = np.where(rm <= R, 1., 1 - np.sqrt(np.maximum(0., 1 - (R/rm)**2)))
    return float(M_grid[0] + dM @ w)


def mean_relative_speed(u, f, n=200000, seed=5):
    rng = np.random.default_rng(seed)
    k1, k2 = rng.choice(len(u), n, p=f), rng.choice(len(u), n, p=f)
    d1, d2 = rng.normal(size=(n, 3)), rng.normal(size=(n, 3))
    d1 /= np.linalg.norm(d1, axis=1)[:, None]
    d2 /= np.linalg.norm(d2, axis=1)[:, None]
    return float(np.mean(np.linalg.norm(u[k1, None]*d1 - u[k2, None]*d2, axis=1)))


def summarize(model, res, s, full, M0):
    L, hist = res['ledger'], res['history']
    rr = np.linalg.norm(model.x, axis=1) if len(model.x) else np.zeros(0)
    m = model.m
    Mb, rh = s['M_b'], s['r_half']
    M_T, M_rh = float(m.sum()), float(m[rr < rh].sum())
    grow = None
    if len(hist) >= 2 and hist[-2]['M_confined'] > 0 and M_T > 0:
        grow = math.log(M_T/hist[-2]['M_confined'])/math.log(hist[-1]['t_Gyr']/hist[-2]['t_Gyr'])
    sigma_eff = U if model.uk.size == 1 else U/math.sqrt(3)
    g_bb = mean_relative_speed(model.uk, model.fk)
    inflow = math.pi*model.R_b**2*model.rho_inf*float(model.fk @ model.uk)*T_REF*mc.PER_GYR
    kJ = math.sqrt(4*math.pi*mc.G*model.rho_inf)/sigma_eff
    # the bath's focused excess from its current tables, whether or not its gravity is included in the run
    Mbx, Mb_grid = model.bath_excess_mass(), model.M_b
    lrh = math.log(rh)
    regime = dict(tau_r_half=float(np.interp(lrh, model.lr, model.tau)), tau_centre=float(model.tau[0]),
                  tau_r_half_max_in_history=max(h['tau_r_half'] for h in hist) if hist else None,
                  bath_self_collision_time_Gyr=(1/(model.sm*model.rho_inf*g_bb)/mc.PER_GYR if model.sm > 0 else None),
                  bath_excess_within_r_half_over_baryons=float(np.interp(lrh, model.lr, Mbx)/np.interp(lrh, model.lr, Mb_grid)),
                  uniform_bath_mass_in_R_b_over_baryons=4/3*math.pi*model.R_b**3*model.rho_inf/Mb,
                  capture_efficiency=(M_T - M0)/inflow,
                  jeans_length_over_R_b=sigma_eff*math.sqrt(math.pi/(mc.G*model.rho_inf))/model.R_b,
                  kJ_R_b_over_pi=kJ*model.R_b/math.pi)
    # an added source's births (CC-2 stage 2B's field) enter the ledger too; stage 2A has none
    mass_ledger = M_T - (M0 + L['born_mass'] + L.get('field_born_mass', 0.) + L['captured_mass'] - L['ejected_mass']
                         - L['evaporated_mass'] - L['escaped_mass'] + L['M_thinning'] + L['M_resample'])
    out = dict(M_T=M_T, M_within_r_half_T=M_rh, B1_ratio=M_T/Mb, B2_ratio=M_rh/(.5*Mb), growth_dlnM_dlnt_at_T=grow,
               tracers=len(m), regime=regime, runtime_s=res['runtime_s'],
               mass_ledger_relative=abs(mass_ledger)/max(M_T, M0, 1e-300),
               ledger={k: float(L[k]) for k in LEDGER_KEYS},
               history=[dict(t_Gyr=h['t_Gyr'], M=h['M_confined'], M_r_half=h['M_within_r_half'], tau_r_half=h['tau_r_half'],
                             bath_excess_r_half=h['bath_excess_within_r_half']) for h in hist])
    if not full:
        return out
    # enclosed and projected mass at the test radii (the acceleration ratio equals the enclosed ratio)
    test = {}
    for q in TEST_RADII[s['key']] + (rh,):
        lq = math.log(q)
        mb_q = float(np.interp(lq, model.lr, Mb_grid))
        mc_q = float(m[rr < q].sum())
        test[f'{q:g}'] = dict(companions=mc_q, baryons=mb_q, bath_excess=float(np.interp(lq, model.lr, Mbx)),
                              companions_over_baryons=mc_q/mb_q, tracers=int((rr < q).sum()),
                              projected_companions=projected(rr, m, q), projected_baryons=projected_profile(model, Mb_grid, q),
                              projected_bath_excess=projected_profile(model, Mbx, q))
    out['test_radii'] = test
    a, b = SLOPE_RANGE[s['key']]
    xs = np.geomspace(a, b, 8)
    ys = np.array([float(m[rr < q].sum()) for q in xs])
    out['log_slope_of_enclosed_mass'] = float(np.polyfit(np.log(xs), np.log(ys), 1)[0]) if (ys > 0).all() else None
    # profiles: shell statistics belong to the shell's geometric centre, enclosed masses to its outer edge r_hi. The raw
    # second moments sigma_r and beta include coherent infall; the _about_mean versions subtract the shell's mean v_r
    edges = np.geomspace(1e-2*rh, model.R_b, 25)
    prof = []
    rhat = model.x/np.maximum(rr, 1e-300)[:, None]
    vr = np.einsum('ij,ij->i', model.v, rhat)
    vt2 = np.einsum('ij,ij->i', model.v, model.v) - vr*vr
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (rr >= lo) & (rr < hi)
        ms = float(m[sel].sum())
        row = dict(r_lo_kpc=float(lo), r_hi_kpc=float(hi), r_center_kpc=float(math.sqrt(lo*hi)),
                   density=ms/(4/3*math.pi*(hi**3 - lo**3)), tracers=int(sel.sum()),
                   enclosed=float(m[rr < hi].sum()), baryons_enclosed=float(np.interp(math.log(hi), model.lr, Mb_grid)),
                   bath_excess_enclosed=float(np.interp(math.log(hi), model.lr, Mbx)))
        if sel.sum() >= 10:
            s_r2 = float(m[sel] @ vr[sel]**2)/ms
            s_t2 = float(m[sel] @ vt2[sel])/ms/2
            v_mean = float(m[sel] @ vr[sel])/ms                   # coherent infall (negative) or outflow
            s_rr = max(s_r2 - v_mean*v_mean, 0.)                    # random radial motion about that mean
            row.update(sigma_r=math.sqrt(s_r2), sigma_t=math.sqrt(s_t2), beta=(1 - s_t2/s_r2) if s_r2 > 0 else None,
                       v_r_mean=v_mean, sigma_r_about_mean=math.sqrt(s_rr),
                       beta_about_mean=(1 - s_t2/s_rr) if s_rr > 0 else None)
        prof.append(row)
    out['profile'] = prof
    # energy distribution and the speed distribution at r_half
    phi0 = abs(float(model.phi[0]))
    E = mc.energies(model.x, model.v, model.lr0, model.dl, model.n, model.phi, model.g, model.r_lo) if len(m) else np.zeros(0)
    hE, eE = np.histogram(E/phi0, bins=np.r_[np.linspace(-1, 0, 21), 1.], weights=m)
    out['energy_distribution'] = dict(edges_over_central_depth=eE.tolist(), mass=hE.tolist(), central_depth=phi0)
    shell = (rr > .8*rh) & (rr < 1.25*rh)
    shell &= model.phi_grid_at(rr) < 0          # an escape speed exists only where the potential is negative
    if shell.sum():
        ve = np.sqrt(-2*model.phi_grid_at(rr[shell]))
        sp = np.linalg.norm(model.v[shell], axis=1)
        hv, ev = np.histogram(sp/ve, bins=np.linspace(0, 1.5, 16), weights=m[shell])
        out['speed_distribution_at_r_half'] = dict(edges_over_v_esc=ev.tolist(), mass=hv.tolist(), tracers=int(shell.sum()),
                                                   fraction_above_0_9_v_esc=float(m[shell][sp/ve > .9].sum()/m[shell].sum()))
    # exported companions, by channel
    ex = res['exits']
    names = {0: 'orbit escape', 1: 'ejection by an incoming companion', 2: 'evaporation', 3: 'seedless partner'}
    bins = np.array([0., 50., 100., 200., 300., 500., 1000., 2000., 5000., 1e9])
    spec = {}
    for c, nm in names.items():
        sel = ex['channel'] == c
        vinf = np.sqrt(2*np.maximum(ex['E'][sel], 0.))
        h, _ = np.histogram(vinf, bins=bins, weights=ex['m'][sel])
        spec[nm] = dict(mass=float(ex['m'][sel].sum()), energy=float(ex['m'][sel] @ ex['E'][sel]), mass_by_speed=h.tolist())
    out['exported'] = dict(speed_edges_kms=bins.tolist(), channels=spec)
    return out


def v1_task(spec):
    """V1: the counted seed with no interaction in a static potential for T."""
    s = system(spec['key'])
    mdl = model_for(s, 'S', 1., 1., bath_gravity=False, freeze=True, seed=spec['rng'], channels=())
    mdl.sample_seed(cf1.F_SEED*s['M_b'], np.interp(mdl.lr, s['pot'].lr, s['seed_sigma']), N_SEED)
    E0 = mc.energies(mdl.x, mdl.v, mdl.lr0, mdl.dl, mdl.n, mdl.phi, mdl.g, mdl.r_lo)
    J0 = np.linalg.norm(np.cross(mdl.x, mdl.v), axis=1)
    M0 = float(mdl.m.sum())
    res = mdl.run(T_REF, DELTA_MAX[s['kind']], SNAPS, budget_s=BUDGET)
    out = dict(escapes=res['ledger']['orbit_escapes'], mass_change=float(mdl.m.sum()) - M0, runtime_s=res['runtime_s'],
               substeps=res['ledger']['substeps'], capped=res['ledger']['capped'])
    if len(mdl.x) == len(E0):
        E1 = mc.energies(mdl.x, mdl.v, mdl.lr0, mdl.dl, mdl.n, mdl.phi, mdl.g, mdl.r_lo)
        J1 = np.linalg.norm(np.cross(mdl.x, mdl.v), axis=1)
        sc = abs(mdl._phi1(s['r_half']))
        dE = np.abs(E1 - E0)/sc
        out.update(energy_error_p99=float(np.percentile(dE, 99)), energy_error_max=float(dE.max()),
                   angular_momentum_error_max=float(np.max(np.abs(J1/J0 - 1))))
    out['passed'] = bool(out['escapes'] == 0 and out['mass_change'] == 0 and out.get('energy_error_p99', 1.) < 1e-3
                         and out.get('angular_momentum_error_max', 1.) < 1e-9)
    return out


def v4_task(spec):
    """V4: the incoming-bound kernel over the counted seed (frozen potential, CF-1's E < 0 criterion) against CF-1's
    kernel recomputed with more radii and samples, and against the archived value."""
    s = system(spec['key'])
    pot, rh = s['pot'], s['r_half']
    lo, hi = max(pot.r[1], 1e-2*rh), min(20*rh, .9*pot.R_b)
    mdl = model_for(s, 'S', 1., 1., bath_gravity=False, freeze=True, seed=spec['rng'], channels=())
    mdl.sample_seed(cf1.F_SEED*s['M_b'], np.interp(mdl.lr, pot.lr, s['seed_sigma']), N_V4)
    rr = np.linalg.norm(mdl.x, axis=1)
    keep = (rr >= lo) & (rr <= hi)          # CF-1's kernel integrates over [lo, hi]
    x0, v0, m0 = mdl.x[keep], mdl.v[keep], mdl.m[keep]
    rr = rr[keep]
    sp = np.linalg.norm(v0, axis=1)
    rb = np.interp(np.log(rr), mdl.lr, mdl.rho_bath)
    wmax = np.sqrt(U*U - 2*mdl.phi_grid_at(rr))
    Delta = 1./float(np.max(mdl.sm*rb*(wmax + sp)))
    est = []
    for _ in range(R_V4):
        x, v, m = x0.copy(), v0.copy(), m0.copy()
        alive = np.ones(len(x), np.bool_)
        sx, sv, smass, er = np.zeros_like(x), np.zeros_like(v), np.zeros(len(x)), np.zeros(len(x))
        lb, ps = np.zeros(12), np.zeros(2)
        mc.bath_bound(x, v, m, alive, np.ones(len(x), np.int64), Delta, mdl.sm, mdl.rho_inf, mdl.lr0, mdl.dl, mdl.n,
                      mdl.dtab, mdl.stab, mdl.uk, mdl.phi, mdl.g, mdl.r_lo, True, sx, sv, smass, er, lb, ps)
        est.append(lb[10]/(mdl.sm*mdl.rho_inf*Delta*len(x0)))
    K_mc, se = float(np.mean(est)), float(np.std(est, ddof=1)/math.sqrt(len(est)))
    radii = np.geomspace(lo, hi, 128)
    K_ref = cf1.kernel(pot, U, radii, s['seed_rho'], s['seed_sigma'], NK_V4)['K_kms']
    arch = json.loads((RESULTS/'gravitational-focusing/cf1-results.json').read_text(encoding='utf-8'))
    K_arch = arch['systems'][s['name']]['trials']['300']['K_kms']
    return dict(K_monte_carlo=K_mc, standard_error=se, K_cf1_recomputed=K_ref, K_archived=K_arch,
                relative_to_recomputed=abs(K_mc/K_ref - 1), relative_to_archived=abs(K_mc/K_arch - 1),
                passed=bool(abs(K_mc/K_ref - 1) < .05))


def static_limit_task(spec, lo=1e-2, hi=1e9, rtol=.02):
    """The largest incident density at which the bath with its own excess gravity has a stable static solution
    (mc.Model.relax_bath), by bisection in log density; with the reason the next density up fails."""
    s = system(spec['key'])

    def probe(rho):
        mdl = model_for(s, spec['bath'], spec['sm'], rho, seed=spec['rng'])
        return mdl.bath_static, mdl.static_fail_reason, mdl.bath_gain, mdl
    ok_lo = probe(lo)[0]
    if not ok_lo:
        return dict(limit=None, reason='no stable static bath even at the lowest density probed', lowest=lo)
    if probe(hi)[0]:
        return dict(limit=None, reason='stable static bath up to the highest density probed', highest=hi)
    a, c, why = math.log(lo), math.log(hi), None
    while c - a > math.log(1 + rtol):
        mid = .5*(a + c)
        ok, reason, gain, _ = probe(math.exp(mid))
        if ok:
            a = mid
        else:
            c, why = mid, reason
    ok, reason, gain, mdl = probe(math.exp(a))
    lrh = math.log(s['r_half'])
    Mb_rh = float(np.interp(lrh, mdl.lr, mdl.M_b))
    diag = dict(circular_speed_at_R_b_over_fastest_incident=math.sqrt(mc.G*float(mdl.M[-1])/mdl.R_b)/float(mdl.uk.max()),
                bath_excess_in_R_b_over_baryons=float(mdl.M_bx[-1]/mdl.M_b[-1]),
                bath_excess_in_r_half_over_baryons=float(np.interp(lrh, mdl.lr, mdl.M_bx))/Mb_rh,
                lowest_enclosed_over_baryons=float(np.min(mdl.M/np.maximum(mdl.M_b, 1e-300))),
                central_depth_over_baryons_only=float(mdl.phi[0]/model_for(s, spec['bath'], spec['sm'], math.exp(a),
                                                                             bath_gravity=False).phi[0]),
                tau_r_half=float(np.interp(lrh, mdl.lr, mdl.tau)), tau_centre=float(mdl.tau[0]))
    return dict(limit=math.exp(a), limit_over_mean=math.exp(a)/cf1.RHO_MEAN, reason_above=why,
                gains_at_limit=gain, at_limit=diag)


def task(spec):
    """One run in a worker process. An unexpected exception is returned as an 'error' status with its traceback, so
    that one failing run is reported with its specification instead of stopping every other run."""
    try:
        return _task(spec)
    except Exception as e:                                  # noqa: BLE001 (reported, never silently dropped)
        import traceback
        return dict(spec=spec, status=f'error: {type(e).__name__}: {e}', traceback=traceback.format_exc())


def _task(spec):
    np.seterr(over='raise', invalid='raise', divide='raise')
    t0 = time.time()
    if spec['role'] == 'static limit':
        return dict(spec=spec, status='completed', **static_limit_task(spec), wall_s=time.time() - t0)
    if spec['role'] == 'V1':
        return dict(spec=spec, status='completed', **v1_task(spec), wall_s=time.time() - t0)
    if spec['role'] == 'V4':
        return dict(spec=spec, status='completed', **v4_task(spec), wall_s=time.time() - t0)
    s = system(spec['key'])
    model = model_for(s, spec['bath'], spec['sm'], spec['rho'], self_gravity=not spec.get('frozen', False),
                      freeze=spec.get('frozen', False), bath_gravity=spec.get('bath_gravity', True), seed=spec['rng'])
    out = dict(spec=spec, bath_static=bool(model.bath_static), relax_iterations=int(model.relax_iterations),
               bath_gain=model.bath_gain, static_fail_rho=model.static_fail_rho, static_fail_reason=model.static_fail_reason)
    if not model.bath_static:
        out['status'] = 'no static bath'
        return out
    M0 = 0.
    if spec.get('seed_control'):
        # seed tracers no lighter than the lightest birth tracer (see mc.Model.set_shell_masses)
        model.make_pools()
        model.set_shell_masses(T_REF*mc.PER_GYR, int(N_TARGET*spec.get('n_factor', 1.)))
        M_seed = cf1.F_SEED*s['M_b']
        m_seed = min(max(M_seed/N_SEED, model.m_floor), model.m_top)
        n_seed = max(50, int(round(M_seed/m_seed)))
        model.sample_seed(M_seed, np.interp(model.lr, s['pot'].lr, s['seed_sigma']), n_seed)
        M0 = float(model.m.sum())
    try:
        res = model.run(T_REF, DELTA_MAX[s['kind']]*spec.get('dmax_factor', 1.), SNAPS,
                        n_target=int(N_TARGET*spec.get('n_factor', 1.)), n_max=int(N_MAX*spec.get('n_factor', 1.)),
                        budget_s=BUDGET)
    except TimeoutError as e:
        out['status'] = f'not completed: {e}'
        return out
    except mc.RunawayError as e:
        out['status'] = f'not completed: runaway, {e}'
        return out
    except mc.StaticBathLost as e:
        out['status'] = f'not completed: static bath lost, {e}'
        return out
    out['status'] = 'completed'
    out.update(summarize(model, res, s, spec.get('full', False), M0))
    return out


def short_run():
    """The suite's regression anchor: a zero-seed Milky Way run (bath S, 1 cm^2/g, 2e4 Msun/kpc^3, 0.5 Gyr)."""
    mdl = model_for(system('MW'), 'S', 1., 2e4, seed=4242)
    res = mdl.run(.5, .005, (.25, .5), n_target=2000, budget_s=600)
    L = res['ledger']
    return dict(M=float(mdl.m.sum()), tracers=len(mdl.m), births=L['births'], captures=L['captures'],
                ejections=L['ejections'], closure=L['closure_max'], max_dK=L['max_dK'])


def run_tasks(specs, label):
    log(f'{label}: {len(specs)} runs on {WORKERS} workers')
    results = []
    with Pool(WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(task, specs)):
            results.append(r)
            sp = r['spec']
            extra = ''
            if sp['role'] == 'static limit':
                extra = f"limit={r.get('limit')} ({r.get('reason_above') or r.get('reason')})  {r.get('wall_s', 0):.0f}s"
            elif r['status'] == 'completed' and 'B1_ratio' in r:
                extra = f"M(T)/M_b={r['B1_ratio']:.3g}  M(<r_half)/M_b(<r_half)={r['B2_ratio']:.3g}  {r['runtime_s']:.0f}s"
            elif r['status'] == 'completed':
                extra = f"passed={r.get('passed')}  {r.get('wall_s', 0):.0f}s"
            log(f"  [{i + 1}/{len(specs)}] {sp['role']} {sp['key']} {sp.get('bath', '')} {sp.get('sm', '')} "
                f"rho={sp.get('rho', 0):.3g} gravity={sp.get('bath_gravity', True)}: {r['status']} {extra}")
            if r['status'].startswith('error'):
                log('  ERROR TRACEBACK\n' + r.get('traceback', ''))
            if time.time() - T0 > 10*3600:
                raise RuntimeError('overall wall-clock budget exhausted')
    # imap_unordered hands results back as they finish; restore the order the tasks were issued in (each spec's seed is
    # unique and increasing), so that nothing downstream, the next round's seeds in particular, depends on timing
    results.sort(key=lambda r: r['spec']['rng'])
    return results


# ----------------------------------------------------------------------------------------------------------------
# controls
def c0_regression(arch):
    """C0: CF-1's growth rule with its archived kernel tables reproduces its archived exposures."""
    rows = {}
    for name, a in arch['systems'].items():
        t = a['trials']['300']
        A = t['exposure_for_retained_equal_to_baryons']
        if not A:
            continue
        u_tab, K_tab = np.array(a['kernel_table']['u_kms']), np.array(a['kernel_table']['K_kms'])
        q = cf1.growth(u_tab, K_tab, 300., [A])[0]['q_final']
        rows[name] = dict(archived_exposure=A, q_final_at_archived=q, relative=abs(q - 1))
    return dict(rows=rows, passed=bool(all(r['relative'] < 1e-5 for r in rows.values())))


def production_coefficient(s, u, f, n_per=4000, seed=99):
    """Seedless production per (sigma/m) rho_inf^2 in the baryonic potential with a transparent, weightless bath."""
    mdl = mc.Model(s['r'], s['M'], s['R_b'], s['r_half'], u, f, 1., 1., bath_gravity=False, seed=seed)
    mdl.make_pools(n_per=n_per)
    return mdl.production_rate()/(mdl.sm*mdl.rho_inf**2)


def c1_required(I_S, K, I_evap, M_seed, target, sm_cm2g):
    """C1: dM/dt = sm rho^2 I_S + sm rho K M - sm I_evap (M/M_seed)^2 from M = 0; the rho giving M(T) = target."""
    T = T_REF*mc.PER_GYR
    sm = sm_cm2g*mc.CM2_PER_G

    def MT(lr):
        rho = math.exp(lr)

        def rhs(t, y):
            M = max(y[0], 0.)
            return [sm*rho*rho*I_S + sm*rho*K*M - sm*I_evap*(M/M_seed)**2]
        sol = solve_ivp(rhs, (0, T), [0.], method='LSODA', rtol=1e-9, atol=1e-9*target)
        return float(sol.y[0][-1])
    f = lambda lr: math.log(max(MT(lr), 1e-300)) - math.log(target)
    lo, hi = math.log(1e-6), math.log(1e14)
    if f(lo) > 0 or f(hi) < 0:
        return None
    return math.exp(brentq(f, lo, hi, xtol=1e-9))


def interpolate(runs, key, target):
    """Log-log interpolation of a run quantity against rho; returns (rho or None, bracketed)."""
    pts = sorted((r['spec']['rho'], r[key]) for r in runs if r['status'] == 'completed' and r[key] > 0)
    if len(pts) < 2:
        return None, False
    lr, lm, lt = np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), math.log(target)
    for i in range(len(pts) - 1):
        if (lm[i] - lt)*(lm[i + 1] - lt) <= 0 and lm[i + 1] != lm[i]:
            return float(math.exp(lr[i] + (lt - lm[i])*(lr[i + 1] - lr[i])/(lm[i + 1] - lm[i]))), True
    return None, False


def ladder_spec(k, b, sm, rho, rng, gravity=True, role='ladder', **kw):
    return dict(key=k, bath=b, sm=sm, rho=rho, role=role, rng=rng, bath_gravity=gravity, **kw)


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    arch = json.loads((RESULTS/'gravitational-focusing/cf1-results.json').read_text(encoding='utf-8'))
    arch2 = json.loads((RESULTS/'companion-supply/cc2a-results.json').read_text(encoding='utf-8'))
    system('MW')
    S = _SYS
    out = dict(scope='CC-2 stage 2A: supply needed with all three collision classes of CF-1\'s law evolving together '
                     'with mass, energy, orbits, bath depletion and gravity. Benchmarks, the 10 Gyr span and the '
                     'incident distributions are labeled trials; the cosmic mean density is a comparison unit.',
               smoke=SMOKE, T_Gyr=T_REF, validation={}, controls={}, combos={})
    out['controls']['C0'] = c0_regression(arch)
    log(f"C0 passed={out['controls']['C0']['passed']}")
    # ---- V7: potential and a sampled Plummer sphere
    mw = S['MW']
    m0 = model_for(mw, 'S', 1., 1., bath_gravity=False)
    rr = np.geomspace(1e-2*mw['r_half'], .9*mw['R_b'], 40)
    phi_c = np.interp(np.log(rr), mw['pot'].lr, mw['pot'].phi)
    pointwise = float(np.max(np.abs(m0.phi_grid_at(rr)/phi_c - 1)))
    central = float(np.max(np.abs(m0.phi_grid_at(rr) - phi_c))/abs(mw['pot'].phi[1]))
    rng = np.random.default_rng(3)
    Np, a = 20000, 5.
    rp = a/np.sqrt(rng.random(Np)**(-2/3) - 1)
    rgrid = np.geomspace(1e-3, 1e4, 4000)
    m_pl = mc.Model(rgrid, 1e10*rgrid**3/(rgrid**2 + a*a)**1.5, 300., a, *mc.single_speed(300.), 1., 1., bath_gravity=False)
    dirs = rng.normal(size=(Np, 3))
    m_pl.add_tracers(rp[:, None]*dirs/np.linalg.norm(dirs, axis=1)[:, None], np.zeros((Np, 3)), np.full(Np, 1./Np))
    ana = np.minimum(m_pl.r, 300.)**3/(np.minimum(m_pl.r, 300.)**2 + a*a)**1.5
    ks = float(np.max(np.abs(m_pl.enclosed_tracer_mass() - ana)))
    out['validation']['V7'] = dict(phi_pointwise_to_0_9_R_b=pointwise, phi_over_central_depth=central,
                                   plummer_max_cumulative_deviation=ks, plummer_limit=1.63/math.sqrt(Np),
                                   passed=bool(pointwise < 1e-4 and central < 1e-4 and ks < 1.63/math.sqrt(Np)))
    log(f"V7 {out['validation']['V7']}")
    # ---- V2: bath density against CF-1's orbit integration, one radius barrier-limited (slowest M component)
    v2 = {}
    for r_ in (.5, 5., 50., 200., 295.):
        v2[f'S, r={r_:g}'] = abs(float(np.interp(math.log(r_), m0.lr, m0.rho_bath))/F.density_by_orbits(mw['pot'], r_, U) - 1)
    # a barrier limits the density only below the circular speed at R_b (Milky Way 38 km/s, Coma about 140 km/s):
    # use the slowest component of case M that has one, searching the systems in order
    uM, _ = bath('M')
    r_lim = None
    for kk in ('MW', 'J1630', 'Coma low', 'Coma high'):
        for u_slow in uM[:8]:
            sk = S[kk]
            mslow = mc.Model(sk['r'], sk['M'], sk['R_b'], sk['r_half'], [float(u_slow)], [1.], 1., 1., bath_gravity=False)
            lim = np.nonzero(mslow.stab[0] < .99)[0]
            if len(lim):
                r_lim = float(mslow.r[lim[len(lim)//2]])
                v2[f'{kk}, {u_slow:.1f} km/s, barrier-limited r={r_lim:.4g}'] = abs(
                    float(np.interp(math.log(r_lim), mslow.lr, mslow.rho_bath))/F.density_by_orbits(sk['pot'], r_lim, float(u_slow)) - 1)
                break
        if r_lim:
            break
    samp = mc.sample_bath_many(200000, 5., m0.lr0, m0.dl, m0.n, m0.dtab, m0.stab, m0.uk, m0.phi, m0.g, m0.r_lo)
    w5 = math.sqrt(U*U - 2*m0._phi1(5.))
    speed_err = float(np.max(np.abs(np.linalg.norm(samp, axis=1)/w5 - 1)))
    inward, se_in = float(np.mean(samp[:, 0] < 0)), math.sqrt(.25/len(samp))
    out['validation']['V2'] = dict(density_relative=v2, barrier_limited_radius=r_lim, speed_max_relative=speed_err,
                                   inward_fraction=inward, inward_standard_error=se_in,
                                   passed=bool(max(v2.values()) < 1e-3 and r_lim is not None and speed_err < 1e-12
                                               and abs(inward - .5) < 3*se_in))
    log(f"V2 passed={out['validation']['V2']['passed']} max={max(v2.values()):.2e} barrier r={r_lim}")
    # ---- V3: seedless production against stage 1's analytic kernel where no barrier applies
    v3 = {}
    for r_ in (.5, 3., mw['r_half'], 20., 60.):
        ve = math.sqrt(mw['pot'].v_esc2(r_))
        an = cc2a.seedless_analytic(U, ve)[0]
        ev = np.empty((400000, 15))
        mc.seedless_events(np.full(400000, r_), m0.lr0, m0.dl, m0.n, m0.dtab, m0.stab, m0.uk, m0.phi, m0.g, m0.r_lo, ev)
        mcb = float(np.mean(ev[:, 1]*ev[:, 14]))
        v3[f'v_esc/u={ve/U:.2f}'] = dict(analytic=an, monte_carlo=mcb, relative=abs(mcb/an - 1),
                                         confined_over_bound=float(np.mean(ev[:, 1]*ev[:, 2]))/mcb)
    out['validation']['V3'] = dict(rows=v3, passed=bool(max(q['relative'] for q in v3.values()) < 2e-2))
    log(f"V3 passed={out['validation']['V3']['passed']}")
    out['checks_short_run'] = short_run()
    log(f"suite anchor {out['checks_short_run']}")
    # ---- C1 per system and incident distribution
    c1 = {}
    for k, s in S.items():
        a = arch['systems'][s['name']]
        A_old, K300 = a['trials']['300']['exposure_for_retained_equal_to_baryons'], a['trials']['300']['K_kms']
        tab_u, tab_K = np.array(a['kernel_table']['u_kms']), np.array(a['kernel_table']['K_kms'])
        ch = arch2['systems'][s['name']]['trials']['300']['channels']['1 cm2/g']
        M_seed = cf1.F_SEED*s['M_b']
        I_evap = ch['evaporation_over_seed_capture']*A_old*K300*M_seed/mc.CM2_PER_G
        c1[k] = dict(A_star_old=A_old, I_evap=I_evap, baths={})
        for b in BATHS:
            I_S = production_coefficient(s, *bath(b))
            K = K300 if b == 'S' else cc2a.maxwell_kernel(tab_u, tab_K, U/math.sqrt(3))[0]
            c1[k]['baths'][b] = dict(I_S=I_S, K_kms=K,
                                     rho_required_B1={f'{sm:g}': c1_required(I_S, K, I_evap, M_seed, s['M_b'], sm) for sm in SIGMA_M},
                                     old_rho={f'{sm:g}': A_old/(sm*mc.CM2_PER_G) for sm in SIGMA_M})
        log(f"C1 {k}: " + '; '.join(f"{b}: " + ' '.join(f"{v:.3g}" for v in c1[k]['baths'][b]['rho_required_B1'].values())
                                    for b in BATHS))
    out['controls']['C1'] = c1
    # ---- the whole process's velocity dependence with frozen coefficients (C1), at 1 cm^2/g
    if not SMOKE:
        scan = {}
        for k, s in S.items():
            a = arch['systems'][s['name']]
            tab_u, tab_K = np.array(a['kernel_table']['u_kms']), np.array(a['kernel_table']['K_kms'])
            M_seed = cf1.F_SEED*s['M_b']
            rows = {}
            for disp in (100., 200., 300., 500., 1000.):
                I_S = production_coefficient(s, *mc.maxwell_components(disp), n_per=2000, seed=98)
                K = cc2a.maxwell_kernel(tab_u, tab_K, disp)[0]
                rows[f'Maxwellian {disp:g}'] = dict(I_S=I_S, K_kms=K, rho_required_B1_at_1=c1_required(I_S, K, c1[k]['I_evap'], M_seed, s['M_b'], 1.))
            K300, K3000 = a['trials']['300']['K_kms'], a['trials']['3000']['K_kms']
            for fr in (0., .0034, .01, .1, .34, .5):
                uf = mc.mixture([(1 - fr, mc.single_speed(300.)), (fr, mc.single_speed(3000.))]) if fr else mc.single_speed(300.)
                I_S = production_coefficient(s, *uf, n_per=2000, seed=97)
                K = (1 - fr)*K300 + fr*K3000
                rows[f'fast fraction {fr:g}'] = dict(I_S=I_S, K_kms=K, rho_required_B1_at_1=c1_required(I_S, K, c1[k]['I_evap'], M_seed, s['M_b'], 1.))
            scan[k] = rows
            log(f'C1 velocity scan {k} done')
        out['controls']['C1_velocity_scan_at_1cm2g'] = scan
    # ---- round 1: V1, V4 and the ladders
    rng_seed = [1000]

    def nxt():
        rng_seed[0] += 1
        return rng_seed[0]
    combos = [(k, b, sm) for k in S for b in BATHS for sm in SIGMA_M]
    # ---- round 0: the stable static-bath limit of every combination (bath gravity on)
    lim_runs = run_tasks([dict(key=k, bath=b, sm=sm, role='static limit', rng=nxt()) for (k, b, sm) in combos], 'round 0')
    limits = {(r['spec']['key'], r['spec']['bath'], r['spec']['sm']): r for r in lim_runs}
    out['static_bath_limit'] = {f'{k} | {b} | {sm:g} cm2/g': {q: v for q, v in r.items() if q != 'spec'}
                                for (k, b, sm), r in limits.items()}
    specs = [dict(key=k, role=r, rng=nxt()) for k in S for r in ('V1', 'V4')]
    ladders = defaultdict(list)                     # (key, bath, sm, gravity) -> runs
    for (k, b, sm) in combos:
        rho1 = c1[k]['baths'][b]['rho_required_B1'][f'{sm:g}']
        lim = limits[(k, b, sm)].get('limit')
        lim = math.inf if lim is None and 'highest' in limits[(k, b, sm)] else (lim or 0.)
        pts = [rho1*fac for fac in (1/3, 1., 3.)]
        # bath gravity on: the ladder points below the static limit, plus one just below the limit
        on = [q for q in pts if q < .95*lim]
        if math.isfinite(lim) and lim > 0 and 3*rho1 > .95*lim:
            on.append(.9*lim)
            if len(on) < 2:
                on.append(.3*lim)
        for q in sorted(set(on)):
            specs.append(ladder_spec(k, b, sm, q, nxt()))
        # bath gravity omitted (labeled sensitivity), always
        for q in pts:
            specs.append(ladder_spec(k, b, sm, q, nxt(), gravity=False))
    v8 = not SMOKE          # resolution: the bath-gravity-omitted ladder, which exists at every density
    if v8:
        rho1 = c1['MW']['baths']['S']['rho_required_B1']['1']
        for fac in (1/3, 1., 3.):
            specs.append(ladder_spec('MW', 'S', 1., rho1*fac, nxt(), gravity=False, role='V8', n_factor=2., dmax_factor=.5))
    results = run_tasks(specs, 'round 1')
    v1 = {r['spec']['key']: r for r in results if r['spec']['role'] == 'V1'}
    v4 = {r['spec']['key']: r for r in results if r['spec']['role'] == 'V4'}
    v8runs = [r for r in results if r['spec']['role'] == 'V8']
    for r in results:
        if r['spec']['role'] == 'ladder':
            sp = r['spec']
            ladders[(sp['key'], sp['bath'], sp['sm'], sp['bath_gravity'])].append(r)
    out['validation']['V1'] = dict(rows={k: {q: v for q, v in r.items() if q != 'spec'} for k, r in v1.items()},
                                   passed=bool(all(r['passed'] for r in v1.values())))
    out['validation']['V4'] = dict(rows={k: {q: v for q, v in r.items() if q != 'spec'} for k, r in v4.items()},
                                   passed=bool(all(r['passed'] for r in v4.values())))
    log(f"V1 passed={out['validation']['V1']['passed']}  V4 passed={out['validation']['V4']['passed']}")
    # ---- rounds 2: extend brackets; a bath-gravity-omitted sensitivity ladder where the bath has no static solution
    targets = lambda k: (('M_T', S[k]['M_b']), ('M_within_r_half_T', .5*S[k]['M_b']))
    for ext in range(MAX_EXTEND + 1):
        new = []
        for (k, b, sm, grav), runs in list(ladders.items()):
            done = {round(math.log(r['spec']['rho']), 9) for r in runs}
            if ext == MAX_EXTEND:
                continue
            comp = [r for r in runs if r['status'] == 'completed']
            if not comp:
                continue
            failed_hi = min([r['spec']['rho'] for r in runs if r['status'] != 'completed'] or [math.inf])
            if grav:
                lim = limits[(k, b, sm)].get('limit')
                if lim:
                    failed_hi = min(failed_hi, .95*lim)
            for key_, target in targets(k):
                rho_i, ok = interpolate(comp, key_, target)
                if ok or not comp:
                    continue
                vals = sorted((r['spec']['rho'], r[key_]) for r in comp)
                pos = [(q, v) for q, v in vals if v > 0]

                def step(p0, p1, toward):
                    """Factor for the next density: a log-log secant step toward the target, overshooting by 1.5
                    so that it brackets, limited to 3-100."""
                    slope = 2.
                    if p0 and p1 and p1[0] != p0[0]:
                        slope = min(4., max(.5, math.log(p1[1]/p0[1])/math.log(p1[0]/p0[0])))
                    return min(100., max(3., (toward**(1/slope))*1.5))
                top = max(comp, key=lambda q: q['spec']['rho'])
                if all(v < target for _, v in vals) and top['regime']['tau_r_half'] > 3:
                    continue            # the bath is already opaque inside r_half: higher densities leave the model
                if all(v < target for _, v in vals):
                    last = pos[-1] if pos else None
                    prev = pos[-2] if len(pos) > 1 else None
                    fac = step(prev, last, target/last[1]) if last else 10.
                    rho_n = vals[-1][0]*fac
                    if rho_n >= failed_hi:
                        if failed_hi/vals[-1][0] > 1.3:
                            rho_n = math.sqrt(vals[-1][0]*failed_hi)
                        else:
                            continue
                elif all(v > target for _, v in vals):
                    first = pos[0]
                    second = pos[1] if len(pos) > 1 else None
                    rho_n = vals[0][0]/step(first, second, first[1]/target)
                else:
                    continue
                if round(math.log(rho_n), 9) not in done and not any(abs(math.log(x['rho']/rho_n)) < 1e-9 for x in new
                                                                     if (x['key'], x['bath'], x['sm'], x['bath_gravity']) == (k, b, sm, grav)):
                    new.append(ladder_spec(k, b, sm, rho_n, nxt(), gravity=grav))
        if not new:
            break
        for r in run_tasks(new, f'round 2.{ext + 1}'):
            sp = r['spec']
            ladders[(sp['key'], sp['bath'], sp['sm'], sp['bath_gravity'])].append(r)
    # ---- round 3: verification at the B1 requirement, controls, and a B2 verification
    specs = []
    req = {}
    for (k, b, sm, grav), runs in ladders.items():
        rB1, okB1 = interpolate(runs, 'M_T', S[k]['M_b'])
        rB2, okB2 = interpolate(runs, 'M_within_r_half_T', .5*S[k]['M_b'])
        req[(k, b, sm, grav)] = (rB1, rB2)
        if rB1:
            specs.append(ladder_spec(k, b, sm, rB1, nxt(), gravity=grav, role='verify B1', full=True))
            specs.append(ladder_spec(k, b, sm, rB1, nxt(), gravity=grav, role='seed control', seed_control=True))
            specs.append(ladder_spec(k, b, sm, rB1, nxt(), gravity=grav, role='frozen control', frozen=True))
        elif grav and limits[(k, b, sm)].get('limit'):
            # B1 not reached below the static limit: the state that the highest stable static bath forms
            specs.append(ladder_spec(k, b, sm, .9*limits[(k, b, sm)]['limit'], nxt(), gravity=True,
                                     role='at the static limit', full=True))
        if rB2:
            # B2 lies where the envelope has grown to many times the baryons, and uniform thinning leaves few
            # tracers inside r_half: these runs carry four times the tracers
            specs.append(ladder_spec(k, b, sm, rB2, nxt(), gravity=grav, role='verify B2', full=True, n_factor=4.))
    if v8:
        rV8, _ = interpolate(v8runs, 'M_T', S['MW']['M_b'])
        if rV8:
            specs.append(ladder_spec('MW', 'S', 1., rV8, nxt(), gravity=False, role='V8 verify', n_factor=2., dmax_factor=.5))
    final = run_tasks(specs, 'round 3')
    by_role = defaultdict(dict)
    for r in final:
        sp = r['spec']
        by_role[(sp['key'], sp.get('bath'), sp.get('sm'), sp.get('bath_gravity', True))][sp['role']] = r
    # ---- round 4: where growth steepens between ladder points, a verification run can miss its benchmark by far;
    # add it to the ladder, interpolate again, and verify again (up to three times) until it lies within 25%
    for roles in by_role.values():          # the controls compare with the verification at their own density
        if 'verify B1' in roles:
            roles['verify B1 at the controls'] = roles['verify B1']
    for ref in range(3):
        specs = []
        for combo, roles in by_role.items():
            if combo[0] not in S or combo not in ladders:
                continue
            k = combo[0]
            for role, key_, target in (('verify B1', 'M_T', S[k]['M_b']), ('verify B2', 'M_within_r_half_T', .5*S[k]['M_b'])):
                v = roles.get(role)
                if not v or v['status'] != 'completed' or .8 <= v[key_]/target <= 1.25:
                    continue
                ladders[combo].append(v)
                rho_n, ok = interpolate(ladders[combo], key_, target)
                if ok and abs(math.log(rho_n/v['spec']['rho'])) > 1e-3:
                    specs.append(ladder_spec(k, combo[1], combo[2], rho_n, nxt(), gravity=combo[3], role=role, full=True,
                                             n_factor=(4. if role == 'verify B2' else 1.)))
        if not specs:
            break
        for r in run_tasks(specs, f'round 4.{ref + 1}'):
            sp = r['spec']
            by_role[(sp['key'], sp['bath'], sp['sm'], sp['bath_gravity'])][sp['role']] = r
    for combo in ladders:
        rB1, _ = interpolate(ladders[combo], 'M_T', S[combo[0]]['M_b'])
        rB2, _ = interpolate(ladders[combo], 'M_within_r_half_T', .5*S[combo[0]]['M_b'])
        req[combo] = (rB1, rB2)
    # ---- aggregation and declared labels
    worst_dK = worst_dP = worst_closure = worst_mass = 0.
    worst_where = None
    for (k, b, sm, grav), runs in ladders.items():
        s = S[k]
        rB1, rB2 = req[(k, b, sm, grav)]
        roles = by_role.get((k, b, sm, grav), {})
        allruns = runs + list(roles.values())
        for r in allruns:
            if r['status'] == 'completed':
                worst_dK = max(worst_dK, r['ledger']['max_dK'])
                worst_dP = max(worst_dP, r['ledger']['max_dP'])
                if r['ledger']['closure_max'] > worst_closure:
                    worst_closure = r['ledger']['closure_max']
                    worst_where = f"{k} | {b} | {sm:g} | gravity {grav} | {r['spec']['role']} | rho {r['spec']['rho']:.4g}"
                worst_mass = max(worst_mass, r['mass_ledger_relative'])
        old = c1[k]['baths'][b]['old_rho'][f'{sm:g}']
        ver = roles.get('verify B1')
        entry = dict(system=s['name'], bath=b, sigma_over_m_cm2g=sm, bath_gravity=grav,
                     C1_rho_B1=c1[k]['baths'][b]['rho_required_B1'][f'{sm:g}'], old_rho=old,
                     ladder=[dict(rho=r['spec']['rho'], role=r['spec']['role'], status=r['status'], M_T=r.get('M_T'),
                                  M_r_half_T=r.get('M_within_r_half_T'),
                                  tau_r_half=(r['regime']['tau_r_half'] if r['status'] == 'completed' else None),
                                  closure=(r['ledger']['closure_max'] if r['status'] == 'completed' else None),
                                  static=r.get('bath_static'), static_fail_reason=r.get('static_fail_reason'))
                             for r in sorted(runs, key=lambda q: q['spec']['rho'])],
                     rho_B1=rB1, rho_B2=rB2,
                     A_B1_per_kpc=(sm*mc.CM2_PER_G*rB1 if rB1 else None),
                     rho_B1_over_mean=(rB1/cf1.RHO_MEAN if rB1 else None), rho_B2_over_mean=(rB2/cf1.RHO_MEAN if rB2 else None),
                     rho_B1_over_old=(rB1/old if rB1 else None), rho_B2_over_old=(rB2/old if rB2 else None),
                     runs={role: {q: v for q, v in r.items() if q != 'spec'} for role, r in roles.items()})
        labels = {}
        if rB1:
            labels['seedless_lowers_supply_substantially'] = bool(rB1 <= .1*old)
        if ver and ver['status'] == 'completed':
            tmax = ver['regime']['tau_r_half_max_in_history']
            labels['transparent'] = bool(tmax < .1)
            labels['outside_model_regime'] = bool(tmax > .3)
            labels['bath_gravity_minor_at_tested_radii'] = bool(abs(ver['regime']['bath_excess_within_r_half_over_baryons']) < .1)
            labels['bath_self_relaxed'] = bool(ver['regime']['bath_self_collision_time_Gyr'] < T_REF)
            # added after the first canonical run: sub-cycled collisions stop at 2^10 sub-steps per step
            labels['collisions_resolved'] = bool(ver['ledger']['capped_mass_fraction_max'] < .01 and ver['ledger']['P_over_1'] == 0)
            sc, fz = roles.get('seed control'), roles.get('frozen control')
            vc = roles.get('verify B1 at the controls') or ver      # the zero-seed run at the controls' density
            if sc and sc['status'] == 'completed' and vc['status'] == 'completed':
                labels['seed_independent'] = bool(abs(sc['M_T'] - vc['M_T'])/vc['M_T'] < .1)
            if fz and fz['status'] == 'completed' and vc['status'] == 'completed':
                labels['self_gravity_matters'] = bool(abs(fz['M_T'] - vc['M_T'])/vc['M_T'] > .1)
        entry['labels'] = labels
        out['combos'][f'{k} | {b} | {sm:g} cm2/g | bath gravity {"on" if grav else "omitted"}'] = entry
    out['validation']['V5'] = dict(max_kinetic=worst_dK, max_momentum=worst_dP, passed=bool(worst_dK < 1e-12 and worst_dP < 1e-12))
    out['validation']['V6'] = dict(energy_closure_max=worst_closure, worst_closure_run=worst_where, mass_ledger_max=worst_mass,
                                   passed=bool(worst_closure < 1e-9 and worst_mass < 1e-9))
    if v8:
        rp = req.get(('MW', 'S', 1., False), (None, None))[0]
        rV8, _ = interpolate(v8runs, 'M_T', S['MW']['M_b'])
        out['validation']['V8'] = dict(rho_B1_primary=rp, rho_B1_doubled_tracers_halved_step=rV8,
                                       relative=(abs(rV8/rp - 1) if (rp and rV8) else None),
                                       passed=bool(rp and rV8 and abs(rV8/rp - 1) < .1))
    out['all_validation_passed'] = all(v['passed'] for v in out['validation'].values())
    out['runtime_seconds'] = time.time() - T0
    text = json.dumps(out, indent=1, default=float)
    log(f"validation: " + ', '.join(f"{k}={v['passed']}" for k, v in out['validation'].items()))
    if SMOKE:
        d = evidence_io.output_dir(args, 'companion-formation-smoke')
        (d/'formation-smoke.json').write_text(text, encoding='utf-8')
        log(f'smoke output {d}')
        return 0
    # wall-clock timings (wall_s, runtime_s) differ on every run and are not results; everything else must reproduce
    status = evidence_io.finish(args, 'companion-formation', text, HERE/'formation-results.json', ignore={'/runtime_seconds'},
                                rules=((r'.*/(wall_s|runtime_s)$', None, None),))
    return status if out['all_validation_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
