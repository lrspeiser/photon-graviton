"""RUT-1 stage 4: the instrumented hundred-period campaign (protocol-rut1.md, stage 4).

    python rut4.py [--output-dir DIR] [--canonical] [--workers N]

The equation is frozen: the two-stage model with tau_keep = 10 T0 and tau_form = 3 T0, the one-stage
baseline (tau_form = 0), the Gaussian footprint, the 10% writing label, and 2% phase and speed jitter,
exactly as in rut3.py. What this stage adds is time -- a hundred initial reference periods, past the
two-stage field's constant-source maturity of 99.04% at 50 T0 -- and the measurements needed to say WHY
the remaining radial motion exists: is it initial adjustment, coherent oscillation, or a slow instability?

Every run is declared below before execution, and the refinement runs cover all three two-stage
configurations over the full horizon, so no case is selected for refinement by its outcome.

This is the science driver and takes tens of minutes on many cores; rut4_checks.py is the suite job. It
writes rut4-results.json plus one compressed series per run under rut4-series/.
"""
import argparse
import gzip
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))

T0 = 2*np.pi
R0 = 1.
TAU_KEEP, TAU_FORM = 10*T0, 3*T0
JITTER, HORIZON, LABEL = .02, 100., .1
CONFIGS = (('16_w0.1_s1', 16, .1, 1), ('32_w0.2_s1', 32, .2, 1), ('32_w0.2_s2', 32, .2, 2))
WINDOWS = dict(early=(0., 20.), mid=(20., 50.), late=(50., 100.))
QUANTITIES = ('mean_radius', 'radial_dispersion', 'azimuthal_dispersion', 'mean_radial_flow', 'support_mean',
              'L_mean', 'L_spread', 'radius_spread')
CHECK_TIMES = (5., 10., 20., 50., 100.)
BLOCK = 5.
# declared thresholds (protocol-rut1.md, stage 4)
MIGRATION_MIN, L_DRIFT_MIN = .01, .01
DISPERSION_REL_MIN, DISPERSION_ABS_MIN = .10, .001
EXCURSION_REL_MIN, EXCURSION_ABS_MIN = .10, .005
RESOLVE_SIGMA, RESOLVE_NUMERICAL = 3., 3.


def _spec(name, cfg, model, prime=False, rotate=0., h_max=.01, eta=.01, spacing_factor=1., band=True):
    label, n, wr, seed = cfg
    return dict(name=name, config=label, writers=n, width_ratio=wr, seed=seed, model=model,
                memory=model != 'no_memory', prime=prime, rotate=rotate, h_max=h_max, eta=eta,
                spacing_factor=spacing_factor, band_diagnostics=band and model != 'no_memory')


def declared_runs():
    runs = []
    for cfg in CONFIGS:
        for model in ('two_stage', 'one_stage', 'no_memory'):
            runs.append(_spec(f'{model}_{cfg[0]}', cfg, model))
    for cfg in CONFIGS[:2]:
        for model in ('two_stage', 'one_stage'):
            runs.append(_spec(f'primed_{model}_{cfg[0]}', cfg, model, prime=True))
    runs.append(_spec('two_stage_32_w0.2_s1_rotated', CONFIGS[1], 'two_stage', rotate=.3))
    for cfg in CONFIGS:
        runs.append(_spec(f'two_stage_{cfg[0]}_half_step', cfg, 'two_stage', h_max=.005, eta=.005, band=False))
        runs.append(_spec(f'two_stage_{cfg[0]}_half_spacing', cfg, 'two_stage', spacing_factor=.5, band=False))
    return runs


def execute(spec, horizon=HORIZON, checkpoints=(5, 10, 20, 50, 100)):
    import formation as FM
    import longrun as LR
    import rut1 as U
    import rut3 as R3
    w = spec['width_ratio']*R0
    q = R3._label_rate(w, TAU_KEEP, LABEL)
    extra = U.ring_inward_acceleration(R0, R0, w, q*TAU_KEEP) if spec['prime'] else 0.
    pos, vel, rates = FM.ring_start(spec['writers'], R0, q, JITTER, JITTER, spec['seed'], extra_inward=extra)
    out = LR.run_instrumented(
        pos, vel, rates, w, TAU_KEEP, horizon*T0, spec['h_max'],
        tau_form=TAU_FORM if spec['model'] == 'two_stage' else 0., eta=spec['eta'],
        spacing=w/5*spec['spacing_factor'], memory=spec['memory'], prime_ring_R=R0 if spec['prime'] else None,
        band_diagnostics=spec['band_diagnostics'], rotate=spec['rotate'],
        checkpoints=tuple(c for c in checkpoints if c <= horizon))
    out['spec'] = spec
    out['writing_rate_total'] = q
    return out


def _worker(spec, series_dir):
    t0 = time.time()
    out = execute(spec)
    out['wall_seconds'] = round(time.time() - t0, 1)
    path = Path(series_dir)/f"{spec['name']}.json.gz"
    # mtime=0 so the compressed bytes, and therefore their hash, depend only on the content
    path.write_bytes(gzip.compress(json.dumps(out, default=float).encode('utf-8'), mtime=0))
    return spec['name'], str(path), out['status'], out['t_final_periods'], out['wall_seconds']


# ---------------------------------------------------------------- analysis
def _trend(t, y, lo, hi, block=None):
    """Window mean and linear trend on block means, which tames the autocorrelation of dense samples."""
    block = BLOCK if block is None else block
    edges = np.arange(lo, hi + 1e-9, block)
    c, m = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (t > a + 1e-9) & (t <= b + 1e-9)
        if sel.sum() >= 3:
            c.append(.5*(a + b))
            m.append(float(np.mean(y[sel])))
    if len(c) < 3:
        return None
    c, m = np.array(c), np.array(m)
    x = c - c.mean()
    slope = float(np.sum(x*(m - m.mean()))/np.sum(x*x))
    resid = m - (m.mean() + slope*x)
    se = float(np.sqrt(np.sum(resid**2)/(len(c) - 2)/np.sum(x*x)))
    return dict(mean=float(m.mean()), slope_per_T0=slope, slope_se=se, change=slope*(hi - lo),
                change_se=se*(hi - lo), blocks=len(c), covers_to=float(c.max() + block/2))


def _series_arrays(series):
    recs = series['records']
    t = np.array([r['t_periods'] for r in recs])
    q = {k: np.array([r[k] for r in recs], float) for k in QUANTITIES if k != 'L_spread'}
    q['L_spread'] = np.array([r['L_max'] - r['L_min'] for r in recs])
    return t, q, recs


def _at(t, y, when):
    if t[-1] < when - 1e-6:
        return None
    return float(y[int(np.argmin(np.abs(t - when)))])


def _band_diagnosis(series, recs, t):
    if not series['spec']['band_diagnostics']:
        return None
    band_recs = [r for r in recs if 'C_band_power' in r]
    if not band_recs:
        return None
    tb = np.array([r['t_periods'] for r in band_recs])
    lo_, hi_ = WINDOWS['late']
    late = [r for r in band_recs if r['t_periods'] >= lo_ - 1e-9]
    out = dict()
    if len(late) >= 2:
        first, last = late[0], late[-1]
        out['late_window_net_torque_by_band'] = {k: last['J_band_sum_over_bodies'][k] -
                                                 first['J_band_sum_over_bodies'][k]
                                                 for k in first['J_band_sum_over_bodies']}
        out['late_window_net_work_by_band'] = {k: last['W_band_sum_over_bodies'][k] -
                                               first['W_band_sum_over_bodies'][k]
                                               for k in first['W_band_sum_over_bodies']}
        out['late_window_net_torque_total'] = last['J_total_sum'] - first['J_total_sum']
        out['late_window_net_work_total'] = last['W_total_sum'] - first['W_total_sum']
        mean_r = float(np.mean([r['mean_radius'] for r in recs if r['t_periods'] >= lo_]))
        rings = (.8, .9, 1., 1.1)
        ring = int(np.argmin(np.abs(np.array(rings) - mean_r)))
        out['ring_used'] = rings[ring]
        tl = np.array([r['t_periods'] for r in late])
        growth, suppression = {}, {}
        for band in late[0]['C_band_power']:
            cpow = np.array([r['C_band_power'][band][ring] for r in late])
            spow = np.array([r['S_band_power'][band][ring] for r in late])
            lp = np.log10(np.maximum(cpow, 1e-300))
            x = tl - tl.mean()
            growth[band] = float(np.sum(x*(lp - lp.mean()))/np.sum(x*x)) if np.sum(x*x) > 0 else 0.
            suppression[band] = float(np.mean(cpow)/max(np.mean(spow), 1e-300))
        out['C_band_log10_growth_per_T0_late'] = growth
        out['C_over_S_band_power_late'] = suppression
    cps = series['checkpoints']
    ka, kb = f'{lo_:g}', f'{hi_:g}'
    if ka in cps and kb in cps:
        per_body = {}
        for band in cps[kb]['J_band']:
            d = np.array(cps[kb]['J_band'][band]) - np.array(cps[ka]['J_band'][band])
            per_body[band] = dict(min=float(d.min()), median=float(np.median(d)), max=float(d.max()))
        out['late_window_per_body_torque_by_band'] = per_body
        dL = np.array(cps[kb]['L_ratio']) - np.array(cps[ka]['L_ratio'])
        out['late_window_per_body_L_change'] = dict(min=float(dL.min()), median=float(np.median(dL)),
                                                    max=float(dL.max()))
    return out


def analyze_run(series):
    t, q, recs = _series_arrays(series)
    trends = {wname: {k: _trend(t, q[k], lo, hi) for k in QUANTITIES} for wname, (lo, hi) in WINDOWS.items()}
    late_mask = t >= WINDOWS['late'][0] - 1e-9
    breathing = None
    if late_mask.sum() > 10:
        breathing = dict(mean_radial_flow_std=float(np.std(q['mean_radial_flow'][late_mask])),
                         radial_dispersion_mean=float(np.mean(q['radial_dispersion'][late_mask])))
    return dict(
        status=series['status'], t_final_periods=series['t_final_periods'],
        wall_seconds=series.get('wall_seconds'),
        eps_residual_rel_max=float(max(r['eps_residual_rel'] for r in recs)),
        torque_residual_max=float(max(r['torque_residual_max'] for r in recs)),
        band_sum_check_worst=series.get('band_sum_check_worst'),
        specific_orbital_energy_max=float(max(r['specific_orbital_energy_max'] for r in recs)),
        windows=trends, late_breathing=breathing,
        at_times={f'{c:g}': {k: _at(t, q[k], c) for k in ('mean_radius', 'radial_dispersion', 'support_mean',
                                                          'L_mean', 'mean_radial_flow')} for c in CHECK_TIMES},
        mean_radius_initial=float(q['mean_radius'][0]),
        mean_radius_at_late_start=_at(t, q['mean_radius'], WINDOWS['late'][0]),
        band_diagnosis=_band_diagnosis(series, recs, t))


def _numerical_uncertainty(analyses):
    unc = {}
    for label, *_ in CONFIGS:
        base = analyses.get(f'two_stage_{label}')
        refs = [analyses.get(f'two_stage_{label}_half_step'), analyses.get(f'two_stage_{label}_half_spacing')]
        if base is None or any(r is None for r in refs):
            continue
        unc[label] = {}
        for wname in WINDOWS:
            unc[label][wname] = {}
            for k in QUANTITIES:
                b = base['windows'][wname][k]
                rs = [r['windows'][wname][k] for r in refs]
                if b is None or any(r is None for r in rs):
                    unc[label][wname][k] = None
                    continue
                unc[label][wname][k] = dict(mean=float(max(abs(b['mean'] - r['mean']) for r in rs)),
                                            change=float(max(abs(b['change'] - r['change']) for r in rs)))
    return unc


def _resolved(trend, unc):
    if trend is None:
        return False
    if abs(trend['change']) <= RESOLVE_SIGMA*trend['change_se']:
        return False
    return unc is None or abs(trend['change']) > RESOLVE_NUMERICAL*unc['change']


def classify(a, unc_cfg, reference_late_radius=None):
    late = a['windows']['late']
    if late['mean_radius'] is None or late['mean_radius']['covers_to'] < WINDOWS['late'][1] - BLOCK + 1e-9:
        return dict(verdict='late window not reached', reached=False)
    u = (lambda k: unc_cfg['late'][k]) if unc_cfg else (lambda k: None)
    migration = _resolved(late['mean_radius'], u('mean_radius')) and abs(late['mean_radius']['change']) > MIGRATION_MIN
    l_drift = _resolved(late['L_mean'], u('L_mean')) and abs(late['L_mean']['change']) > L_DRIFT_MIN
    disp = late['radial_dispersion']
    dispersion_growth = (_resolved(disp, u('radial_dispersion')) and disp['change'] > 0
                         and disp['change'] > DISPERSION_REL_MIN*disp['mean'] and disp['change'] > DISPERSION_ABS_MIN)
    exc = late['radius_spread']
    excursion_growth = (_resolved(exc, u('radius_spread')) and exc['change'] > 0
                        and exc['change'] > EXCURSION_REL_MIN*exc['mean'] and exc['change'] > EXCURSION_ABS_MIN)
    br = a['late_breathing']
    breathing = bool(br and br['mean_radial_flow_std'] > br['radial_dispersion_mean'])
    ref = reference_late_radius
    # measured against the matched no-memory control's late mean: a jittered ring's instantaneous mean radius
    # breathes at the orbital frequency, so comparing two instants would measure epicyclic phase instead
    adjustment = bool(ref is not None and abs(late['mean_radius']['mean'] - ref) > MIGRATION_MIN
                      and not migration)
    sup = late['support_mean']
    positive_support = bool(sup['mean'] > 0 and (unc_cfg is None or sup['mean'] > 3*u('support_mean')['mean']))
    if migration or l_drift or dispersion_growth or excursion_growth:
        verdict = 'candidate heating or instability'
    elif breathing:
        verdict = 'bounded coherent oscillation'
    elif adjustment:
        verdict = 'formation adjustment, then stationary'
    else:
        verdict = 'stationary'
    return dict(reached=True, verdict=verdict, secular_migration=bool(migration), angular_momentum_drift=bool(l_drift),
                dispersion_growth=bool(dispersion_growth), excursion_growth=bool(excursion_growth),
                coherent_breathing=breathing, formation_adjustment=adjustment, positive_support=positive_support,
                no_memory_reference_late_radius=ref,
                numerical_uncertainty_applied=unc_cfg is not None)


def analyze(series_by_name):
    analyses = {name: analyze_run(s) for name, s in series_by_name.items()}
    unc = _numerical_uncertainty(analyses)
    classes = {}
    for name, a in analyses.items():
        spec = series_by_name[name]['spec']
        if not spec['memory'] and 'no_memory' not in name:
            continue
        if name.endswith(('_half_step', '_half_spacing')):
            continue
        unc_cfg = unc.get(spec['config']) if (spec['model'] == 'two_stage' and not spec['prime']
                                               and not spec['rotate']) else None
        control = analyses.get(f"no_memory_{spec['config']}")
        ref = None
        if control is not None and spec['model'] != 'no_memory' and control['windows']['late']['mean_radius']:
            ref = control['windows']['late']['mean_radius']['mean']
        classes[name] = classify(a, unc_cfg, ref)
    common = {}
    for label, *_ in CONFIGS:
        common[label] = {m: analyses[f'{m}_{label}']['at_times'] for m in ('two_stage', 'one_stage', 'no_memory')
                         if f'{m}_{label}' in analyses}
    primed = {}
    for label, *_ in CONFIGS[:2]:
        p2, p1, e2 = (analyses.get(f'primed_two_stage_{label}'), analyses.get(f'primed_one_stage_{label}'),
                      analyses.get(f'two_stage_{label}'))
        if p2 and p1:
            primed[label] = dict(two_stage=p2['at_times'], one_stage=p1['at_times'],
                                 two_stage_class=classes.get(f'primed_two_stage_{label}'),
                                 one_stage_class=classes.get(f'primed_one_stage_{label}'),
                                 empty_field_two_stage_late={k: (e2['windows']['late'][k] or {}).get('mean')
                                                             for k in QUANTITIES} if e2 else None,
                                 primed_two_stage_late={k: (p2['windows']['late'][k] or {}).get('mean')
                                                        for k in QUANTITIES})
    rotation = None
    base, rot, s2 = (analyses.get('two_stage_32_w0.2_s1'), analyses.get('two_stage_32_w0.2_s1_rotated'),
                     analyses.get('two_stage_32_w0.2_s2'))
    if base and rot and s2:
        rows, ok = {}, True
        for k in ('mean_radius', 'radial_dispersion', 'support_mean', 'L_mean'):
            b, r, s = (x['windows']['late'][k] for x in (base, rot, s2))
            if b is None or r is None or s is None:
                rows[k] = None
                continue
            u_ = (unc.get('32_w0.2_s1', {}).get('late', {}).get(k) or {}).get('mean', 0.)
            tol = max(3*u_, abs(b['mean'] - s['mean']))
            diff = abs(b['mean'] - r['mean'])
            rows[k] = dict(unrotated=b['mean'], rotated=r['mean'], difference=diff, tolerance=tol,
                           within=bool(diff <= tol))
            ok = ok and diff <= tol
        rotation = dict(rows=rows, passed=bool(ok),
                        note='a grid-imprinted mode would make the rotated and unrotated runs differ by more '
                             'than the seed-to-seed spread and three times the numerical uncertainty')
    return analyses, unc, classes, common, primed, rotation


def gates(series_by_name, analyses, unc, rotation):
    eps_worst = max(a['eps_residual_rel_max'] for n, a in analyses.items() if series_by_name[n]['spec']['memory'])
    eps_ratio = {}
    for label, *_ in CONFIGS:
        b, h = analyses.get(f'two_stage_{label}'), analyses.get(f'two_stage_{label}_half_step')
        if b and h:
            eps_ratio[label] = float(b['eps_residual_rel_max']/max(h['eps_residual_rel_max'], 1e-300))
    torque_worst = max(a['torque_residual_max'] for a in analyses.values())
    band_rows = {}
    for name, s in series_by_name.items():
        if not s['spec']['band_diagnostics'] or not s['checkpoints']:
            continue
        cp = s['checkpoints'][max(s['checkpoints'], key=float)]
        J = np.array(cp['J_total'])
        Jb = sum(np.array(v) for v in cp['J_band'].values())
        W = np.array(cp['W_total'])
        Wb = sum(np.array(v) for v in cp['W_band'].values())
        band_rows[name] = dict(
            torque_relative=float(np.sum(np.abs(Jb - J))/(np.sum(np.abs(J)) + 1e-6)),
            work_relative=float(np.sum(np.abs(Wb - W))/(np.sum(np.abs(W)) + 1e-6)))
    band_ok = all(r['torque_relative'] <= .02 and r['work_relative'] <= .05 for r in band_rows.values())
    return dict(
        eps_consistency=dict(worst_relative=eps_worst, tolerance=1e-3,
                             refinement_ratio_base_over_half_step=eps_ratio,
                             passed=bool(eps_worst < 1e-3 and all(v > 2 for v in eps_ratio.values()))),
        torque_identity=dict(worst=torque_worst, tolerance=1e-10, passed=bool(torque_worst < 1e-10)),
        band_attribution=dict(rows=band_rows, tolerances=dict(torque=.02, work=.05), passed=bool(band_ok)),
        grid_rotation=rotation)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--canonical', action='store_true')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 4))
    ap.add_argument('--resume-from', type=Path,
                    help='reuse series already produced for an identical declared spec; execute() is '
                         'unchanged, so a reused run is the same computation, and rut4_checks.py replays '
                         'a prefix independently to verify that')
    args = ap.parse_args()
    t0 = time.time()
    out_dir = evidence_io.output_dir(args, 'path-memory-rut4')
    series_dir = out_dir/'rut4-series'
    series_dir.mkdir(parents=True, exist_ok=True)
    runs = declared_runs()
    reused = []
    if args.resume_from:
        for s in runs:
            src = args.resume_from/f"{s['name']}.json.gz"
            if not src.is_file():
                continue
            with gzip.open(src, 'rt', encoding='utf-8') as f:
                stored = json.load(f)
            if stored.get('spec') == s and stored.get('status'):
                (series_dir/src.name).write_bytes(src.read_bytes())
                reused.append(s['name'])
    todo = [s for s in runs if s['name'] not in reused]
    code = {f: hashlib.sha256((HERE/f).read_bytes()).hexdigest()
            for f in ('formation.py', 'longrun.py', 'rut4.py', 'rut3.py', 'rut1.py')}
    log = []
    print(f'  reusing {len(reused)} series, running {len(todo)}', flush=True)
    with ProcessPoolExecutor(max_workers=max(1, min(args.workers, len(todo) or 1))) as pool:
        futures = {pool.submit(_worker, s, str(series_dir)): s['name'] for s in todo}
        for fut in as_completed(futures):
            name, path, status, tf, wall = fut.result()
            log.append(dict(name=name, status=status, t_final_periods=tf, wall_seconds=wall))
            print(f'  done {name:<40} {status:<28} t={tf:6.1f} T0  {wall:7.1f}s', flush=True)
    series = {}
    for s in runs:
        with gzip.open(series_dir/f"{s['name']}.json.gz", 'rt', encoding='utf-8') as f:
            series[s['name']] = json.load(f)
    analyses, unc, classes, common, primed, rotation = analyze(series)
    gate = gates(series, analyses, unc, rotation)
    for a in analyses.values():
        a.pop('wall_seconds', None)
    result = dict(
        experiment='RUT-1 stage 4: the instrumented hundred-period campaign',
        protocol='protocol-rut1.md',
        frozen_equation=dict(two_stage='tau_form dE/dt = S - E, dC/dt = E - C/tau_keep', one_stage='tau_form = 0',
                             tau_keep_T0=10., tau_form_T0=3., writing_label=LABEL, jitter=JITTER,
                             horizon_T0=HORIZON),
        code_sha256_at_launch=code, runs_reused=reused,
        execution_note='the campaign ran across two sessions after the first was interrupted with one run '
                       'outstanding; execute() and the physics modules were unchanged between them, only '
                       'rut4.py gained the resume path, and rut4_checks.py replays a 5 T0 prefix of the '
                       'representative run against its archived series to verify that independently',
        runs=runs, analyses=analyses, numerical_uncertainty=unc, classification=classes,
        common_time_comparison=common, primed_challenge=primed, gates=gate,
        thresholds=dict(migration=MIGRATION_MIN, L_drift=L_DRIFT_MIN, dispersion_relative=DISPERSION_REL_MIN,
                        dispersion_absolute=DISPERSION_ABS_MIN, excursion_relative=EXCURSION_REL_MIN,
                        excursion_absolute=EXCURSION_ABS_MIN, resolve_sigma=RESOLVE_SIGMA,
                        resolve_numerical=RESOLVE_NUMERICAL),
        series_sha256={s['name']: hashlib.sha256((series_dir/f"{s['name']}.json.gz").read_bytes()).hexdigest()
                       for s in runs},
        input_sha256={'protocol-rut1.md': hashlib.sha256((HERE/'protocol-rut1.md').read_bytes()).hexdigest()},
        passed=bool(gate['eps_consistency']['passed'] and gate['torque_identity']['passed']
                    and gate['band_attribution']['passed'] and (rotation or {}).get('passed', False)),
        what_this_is_not='one planar ring around a fixed centre, an instantaneous Gaussian kernel, a writing rate '
                         'set from a label, and no field energy or momentum budget. A verdict is about this frozen '
                         'equation under these conditions, not about gravitational memory in general')
    text = json.dumps(result, indent=1, default=float) + '\n'
    (out_dir/'rut4-log.json').write_text(json.dumps(dict(log=log, wall_seconds=round(time.time() - t0, 1)),
                                                    indent=1), encoding='utf-8')
    if args.canonical:
        canon = HERE/'rut4-series'
        canon.mkdir(exist_ok=True)
        for s in runs:
            (canon/f"{s['name']}.json.gz").write_bytes((series_dir/f"{s['name']}.json.gz").read_bytes())
    status = evidence_io.finish(args, 'path-memory-rut4', text, HERE/'rut4-results.json')
    return status if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
