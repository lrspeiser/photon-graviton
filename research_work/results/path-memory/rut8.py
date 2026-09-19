"""RUT-1 stage 8 (protocol-rut8.md): the linear modes of the constructed populations, predicted before they
are measured.

Two phases, run in this order and committed separately:

    --phase predict    parts L and P: the verified calculation and every prediction, archived as
                       rut8-predictions.json. COMMITTED AND PUSHED BEFORE ANYTHING NEW IS SIMULATED.
    --phase measure    part M: the three new populations built, checked, drawn and evolved with stage 7's code,
                       read with the growth-phase rule the protocol fixed beforehand, and compared with the
                       predictions file, whose hash is recorded.

Three statuses are kept separate, as the owner ruled: reproduction and numerical verification set the exit
status; the scientific outcome never does. Every threshold is read from the protocol. rut8_checks.py is the
suite job; this driver is not.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut8_tasks as TASKS  # noqa: E402
from rut8_tasks import CODE, L4_VARIANTS, STAGE7, TH, T0, _find, _safe, np, run_task  # noqa: E402
import rut7 as S7           # noqa: E402  (stage 7's analysis, by import: D_Q, the V6 gate, series i/o)
import warm_modes as WM     # noqa: E402

TASK_CODE = tuple(f for f in CODE if f not in ('rut8.py', 'rut7.py'))
PROTOCOL = 'protocol-rut8.md'
PREDICTIONS = 'rut8-predictions.json'


# ================================================================ part L and part P, from finished tasks
def assemble_spectrum(rows, populations):
    """Per population and m: the located roots, and the L5 checks on every sub-rectangle of both partitions."""
    L = TH['L']
    out, failures = {}, []
    worst = dict(smallest_gap=None, relative_residual=0., polish_movement=0., partition=0., phase_step=0.)
    for name in populations:
        for m in L['m_values']:
            cell = [r for r in rows if r['name'] == name and r['m'] == m]
            parts = {}
            for shifted in (False, True):
                rects = sorted((r for r in cell if r['shifted'] == shifted), key=lambda r: r['rect'][2])
                tag = f"{name}:m{m}:{'shifted' if shifted else 'base'}"
                for r in rects:
                    if not r['counts_agree']:
                        failures.append(f"{tag}: Im [{r['rect'][2]:+.2f}, {r['rect'][3]:+.2f}] located {r['located']} "
                                        f"but the winding number is {r['winding']['count']}")
                    if r['saturated']:
                        failures.append(f'{tag}: moment capacity saturated')
                    if r['winding']['worst_phase_step'] > L['winding_phase_step']*(1 + 1e-12):
                        failures.append(f'{tag}: winding phase step {r["winding"]["worst_phase_step"]:.3g}')
                    worst['phase_step'] = max(worst['phase_step'], r['winding']['worst_phase_step'])
                    for q in r['roots']:
                        worst['smallest_gap'] = r['gap'] if worst['smallest_gap'] is None else min(worst['smallest_gap'], r['gap'])
                        worst['relative_residual'] = max(worst['relative_residual'], q['relative_residual'])
                        worst['polish_movement'] = max(worst['polish_movement'], q['moved'])
                        if q['relative_residual'] >= L['relative_residual'] or q['moved'] >= L['polish_movement']:
                            failures.append(f'{tag}: a root fails its residual or polish bound')
                parts[shifted] = sorted((complex(q['re'], q['im']) for r in rects for q in r['roots']),
                                        key=lambda z: (round(z.imag, 6), z.real))
            pair = S7.TASKS._match(parts[False], parts[True])
            worst['partition'] = max(worst['partition'], pair['worst_relative'])
            if pair['unpaired_first'] or pair['unpaired_second'] or pair['worst_relative'] > L['partition_agreement']:
                failures.append(f'{name}:m{m}: the two partitions return different roots')
            out[f'{name}:m{m}'] = dict(
                roots=[dict(re=z.real, im=z.imag, e_folding_periods=1/(z.real*T0), pattern_speed=-z.imag/m)
                       for z in sorted(parts[False], key=lambda z: -z.real)],
                rectangles=[len([r for r in cell if not r['shifted']]), len([r for r in cell if r['shifted']])],
                partition_agreement=pair['worst_relative'])
    return out, failures, worst


def rectangle_rows(results):
    return [r for strip in results['strip'].values() for r in strip['rectangles']]


def fit_threshold(family):
    """dE at which the m = 2 root reaches the axis: a quadratic through the resolved roots nearest it, with
    the spread under dropping each point in turn."""
    pts = sorted((r['dE'], r['roots'][0]['re']) for r in family if r['roots'] and r['counts_agree']
                 and r['roots'][0]['re'] < TH['P']['threshold_fit_below'])
    if len(pts) < 4:
        return dict(points=pts, note='fewer than four resolved roots below the fitting limit: no fit')
    x, y = np.array(pts).T

    def crossing(xx, yy):
        roots = np.roots(np.polyfit(xx, yy, 2))
        real = [float(z.real) for z in roots if abs(z.imag) < 1e-9 and z.real > xx.max()]
        return min(real) if real else None

    full = crossing(x, y)
    drops = [crossing(np.delete(x, i), np.delete(y, i)) for i in range(len(x))]
    drops = [d for d in drops if d is not None]
    return dict(points=[[float(a), float(b)] for a, b in pts], threshold_dE=full,
                spread_leave_one_out=[min(drops), max(drops)] if drops else None)


def predict(results, populations):
    L, P = TH['L'], TH['P']
    spectrum, failures, worst = assemble_spectrum(rectangle_rows(results), populations)
    control = results['L5_control']['L5_control']
    gates = {}
    l1 = {k: v for k, v in results['L1'].items() if not k.startswith('control')}
    l1c = results['L1']['control:A_cold']
    gates['L1_orbit_library'] = dict(rows=l1, negative_control=dict(
        what='the same comparison with the memory field removed from the potential', row=l1c,
        rejected=bool(not l1c['passed'])), passed=bool(all(r['passed'] for r in l1.values()) and not l1c['passed']))
    l2 = results['L2']
    gates['L2_response_matrix_time_domain'] = dict(
        rows=l2, negative_controls_rejected=bool(all(r['negative_control']['rejected'] for r in l2.values())),
        passed=bool(all(r['passed'] for r in l2.values())))
    l3 = sorted(results['L3'].values(), key=lambda r: -r['radial_width'])
    diffs = [r['growth_relative_difference'] for r in l3]
    ok = all(d is not None for d in diffs)
    falling = ok and all(b < a for a, b in zip(diffs[:-1], diffs[1:]))
    order = (float(np.polyfit(np.log([r['radial_width'] for r in l3]), np.log(diffs), 1)[0]) if ok else None)
    lo, hi = L['ring_limit_order']
    wrong = l3[-1]['growth_relative_difference_wrong_radius']
    gates['L3_ring_limit'] = dict(
        rows=l3, differences=diffs, falling_at_every_step=bool(falling), observed_order=order,
        negative_control=dict(what="the ring evaluated at stage 5's radius, R = 1, instead of the annulus's own",
                              last_difference=wrong, rejected=bool(wrong is not None and wrong > L['ring_limit_last'])),
        passed=bool(falling and diffs[-1] < L['ring_limit_last'] and lo <= order <= hi
                    and wrong is not None and wrong > L['ring_limit_last']))
    l4 = {}
    for name in sorted({k.split(':')[0] for k in results['L4']}):
        base = complex(*results['L4'][f'{name}:base']['root'])
        moved = {}
        for variant in list(L4_VARIANTS) + ['control_l0_only']:
            root = results['L4'][f'{name}:{variant}']['root']
            moved[variant] = None if root is None else float(abs(complex(*root) - base)/abs(base))
        declared = [moved[v] for v in L4_VARIANTS]
        ctl = moved['control_l0_only']
        l4[name] = dict(base_root=[base.real, base.imag], relative_movement=moved,
                        worst=float(max(declared)) if all(d is not None for d in declared) else None,
                        control_rejected=bool(ctl is None or ctl > L['convergence']))
    gates['L4_convergence'] = dict(rows=l4, negative_controls_rejected=bool(all(r['control_rejected'] for r in l4.values())),
                                   passed=bool(all(r['worst'] is not None and r['worst'] < L['convergence']
                                                   and r['control_rejected'] for r in l4.values())))
    gates['L5_completeness'] = dict(failures=failures, worst=worst, sub_rectangles=len(rectangle_rows(results)),
                                    negative_control=control, passed=bool(not failures and control['rejected']))
    family = sorted(results['family'].values(), key=lambda r: r['dE'])
    new = {}
    for name, dE in P['new_populations'].items():
        roots = spectrum[f'{name}:m2']['roots']
        new[name] = dict(dE=dE, predicted_unstable=bool(roots),
                         growth_rate=roots[0]['re'] if roots else None,
                         pattern_speed=roots[0]['pattern_speed'] if roots else None,
                         e_folding_periods=roots[0]['e_folding_periods'] if roots else None,
                         every_m={f'm{m}': spectrum[f'{name}:m{m}']['roots'] for m in L['m_values']})
    return gates, dict(spectrum=spectrum, family_B=[dict(dE=r['dE'], mean_support=r['mean_support'], sigma_r=r['sigma_r'],
                                                         counts_agree=r['counts_agree'], roots=r['roots'])
                                                    for r in family],
                       threshold=fit_threshold(family), new_populations=new)


# ================================================================ part M, from the archived series
def growth_phase(runs, m=2):
    """The declared rule: from the end of the shot-noise imprint to the first time the ensemble's geometric-mean
    amplitude reaches the declared fraction of its peak; one fit per realization."""
    M = TH['M']
    live = [r for r in runs if r['live'] and r['status'] == 'completed']
    t = np.array([row['t_periods'] for row in live[0]['rows']])
    c = np.array([[complex(*row['coefficients']['C'][m]) for row in r['rows']] for r in live])
    c0 = np.array([[abs(complex(*row['coefficients']['C'][0])) for row in r['rows']] for r in live])
    amp = np.abs(c)/c0
    g = np.exp(np.mean(np.log(np.maximum(amp[:, 1:], 1e-300)), axis=0))
    tg = t[1:]
    peak = int(np.argmax(g))
    t_end = float(tg[np.nonzero(g >= M['growth_end_fraction']*g[peak])[0][0]])
    row = dict(realizations=len(live), peak_relative_amplitude=float(g[peak]), peak_at_periods=float(tg[peak]),
               relative_amplitude_at_start=float(g[np.nonzero(tg >= M['growth_start'])[0][0]]))
    if t_end - M['growth_start'] < M['growth_shortest']:
        return dict(row, growth_phase=False)
    sel = (t >= M['growth_start'] - 1e-9) & (t <= t_end + 1e-9)
    rates = [float(np.polyfit(t[sel]*T0, np.log(amp[k, sel]), 1)[0]) for k in range(len(live))]
    speeds = [float(-np.polyfit(t[sel]*T0, np.unwrap(np.angle(c[k, sel])), 1)[0]/m) for k in range(len(live))]
    return dict(row, growth_phase=True, window_periods=[M['growth_start'], t_end],
                growth_rate=S7._ensemble(rates), pattern_speed=S7._ensemble(speeds))


def reading(pred, cells):
    """The declared reading of one new population against its prediction."""
    M = TH['M']
    top = cells[str(max(M['bodies']))]
    if not pred['predicted_unstable']:
        if not any(c['growth_phase'] for c in cells.values()):
            return 'confirmed'
        if top['growth_phase']:
            gr = top['growth_rate']
            if gr['mean'] - M['confirm_sigma']*gr['standard_error'] > M['stable_rate_floor']:
                return 'contradicted'
        return 'unresolved'
    if not top['growth_phase']:
        return 'contradicted'
    gr, ps = top['growth_rate'], top['pattern_speed']
    off = abs(gr['mean'] - pred['growth_rate'])
    rate_ok = off <= M['confirm_sigma']*gr['standard_error'] or off <= M['confirm_relative']*pred['growth_rate']
    poff = abs(ps['mean'] - pred['pattern_speed'])
    speed_ok = poff <= M['confirm_sigma']*ps['standard_error'] or poff <= M['confirm_pattern']
    if all(c['growth_phase'] for c in cells.values()) and rate_ok and speed_ok:
        return 'confirmed'
    if off > M['confirm_sigma']*gr['standard_error'] and off > M['contradict_relative']*pred['growth_rate']:
        return 'contradicted'
    return 'unresolved'


def evaluate_measurement(series, predictions):
    """Everything in part M that is computed from the archived series; the suite job re-derives it."""
    M = TH['M']
    out = {}
    for name, pred in predictions['part_P']['new_populations'].items():
        cells = {str(n): growth_phase(series['X'][(name, n)]) for n in M['bodies']}
        h = M['horizon'][name]
        windows = dict(early=S7.TH['X']['early_window'], late=[h - 5., h], mode=[h - 15., h])
        table = S7.analyze_X({(name, n): series['X'][(name, n)] for n in M['bodies']}, windows)[name]
        out[name] = dict(prediction={k: pred[k] for k in ('dE', 'predicted_unstable', 'growth_rate', 'pattern_speed',
                                                          'e_folding_periods')},
                         measured=cells, reading=reading(pred, cells), horizon_periods=h,
                         D={n: {q: {k: v for k, v in cell['D'][q].items() if k != 'values'} for q in cell['D']}
                            for n, cell in table.items()})
    x_runs = [r for runs in series['X'].values() for r in runs]
    strip = lambda r: {k: v for k, v in r.items() if k not in ('rows', 'seconds')}
    return dict(V6=S7.gate_V6(series['V6']), comparison=out,
                run_index=dict(V6=[strip(r) for r in series['V6']], X=[strip(r) for r in x_runs]),
                all_runs_completed=bool(all(r['status'] == 'completed' for r in x_runs)),
                pairs_identical=bool(all(a['initial_sha256'] == b['initial_sha256'] for a in x_runs for b in x_runs
                                         if (a['name'], a['bodies'], a['realization']) == (b['name'], b['bodies'], b['realization']))))


def load_series(series_dir):
    out = dict(V6=S7._load_gz(series_dir/'V6.json.gz'), X={})
    for path in sorted(series_dir.glob('X_*.json.gz')):
        runs = S7._load_gz(path)
        out['X'][(runs[0]['name'], runs[0]['bodies'])] = runs
    return out


# ================================================================ the campaign
def predict_tasks(new_pops):
    L, P = TH['L'], TH['P']
    pops = dict(TASKS.stage7_populations(), **new_pops)
    tasks = []
    for name, pop in pops.items():
        cost = len(pop['state']['r'])*0 + (3. if pop['state']['dE'] >= .02 else 1.)
        for m in L['m_values']:
            for shifted in (False, True):
                tasks.append((cost*m, 'strip', f"{name}:m{m}:{'shifted' if shifted else 'base'}", (pop, m, shifted)))
        tasks.append((2., 'L1', name, (pop,)))
    tasks.append((2., 'L1', 'control:A_cold', (pops['A_cold'], False)))
    for name, s_pair in L['time_domain_cases'].items():
        tasks.append((60., 'L2', name, (pops[name], s_pair)))
    for cfg in L['ring_limit']:
        tasks.append((10. + cfg[2]/40., 'L3', f'dL{cfg[0]}', (cfg,)))
    for dE in P['family_B_dE']:
        tasks.append((4., 'family', f'{dE:.3f}', (dE,)))
    return pops, sorted(tasks, key=lambda t: -t[0])


def second_wave(results, pops):
    """L4 and L5's control need the sub-rectangle that holds each root, which the first wave finds."""
    tasks = []

    def holding(name):
        rows = [r for r in rectangle_rows(results) if r['name'] == name and r['m'] == 2 and not r['shifted'] and r['roots']]
        best = max(rows, key=lambda r: max(q['re'] for q in r['roots']))
        return best['rect']

    for name in TH['L']['convergence_populations']:
        rect = holding(name)
        for variant in ['base'] + list(L4_VARIANTS) + ['control_l0_only']:
            tasks.append((8. if variant == 'orbit_quadrature' else 2., 'L4', f'{name}:{variant}', (pops[name], variant, rect)))
    tasks.append((2., 'L5_control', 'L5_control', (pops['B_cold'], holding('B_cold'))))
    return sorted(tasks, key=lambda t: -t[0])


def measure_tasks(new_pops):
    M, S = TH['M'], S7.TH['S']
    tasks = []
    for name, pop in new_pops.items():
        for kind in ('V1', 'V2', 'V3', 'V5'):
            tasks.append((5., kind, name, (pop,)))
        for control in (None,) + S7.V6_CONTROLS:
            for k in range(int(S['stationary_realizations'])):
                tasks.append((20., 'V6', f"{name}:{control or 'declared_draw'}:k{k}", (pop, k, control)))
        for n, reps in zip(M['bodies'], M['realizations']):
            for k in range(reps):
                for live in (True, False):
                    tasks.append(((.1 if live else .04)*n*M['horizon'][name]/30., 'X',
                                  f"{name}:N{n}:k{k}:{'live' if live else 'frozen'}", (pop, n, k, live, M['horizon'][name])))
    return sorted(tasks, key=lambda t: -t[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--phase', choices=('predict', 'measure'), required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--resume-from', type=Path)
    args = ap.parse_args()
    t0 = time.time()
    out_dir = args.output_dir
    (out_dir/'tasks').mkdir(parents=True, exist_ok=True)
    code = {f: hashlib.sha256(_find(f).read_bytes()).hexdigest() for f in CODE}
    task_code = {f: code[f] for f in TASK_CODE}
    log = open(out_dir/f'rut8-{args.phase}-progress.log', 'a', encoding='utf-8')

    def say(msg):
        line = f'[{time.time() - t0:7.0f}s] {msg}'
        print(line, flush=True)
        log.write(line + '\n')
        log.flush()

    def path_of(root, kind, key):
        return root/'tasks'/(f'{kind}__{key}'.replace(':', '_').replace('+', 'p').replace('-', 'n') + '.json.gz')

    results, reused, walls, errors = {}, [], {}, []

    def drain(pool, tasks):
        futures = {}
        for _, kind, key, a in tasks:
            src = path_of(args.resume_from, kind, key) if args.resume_from else None
            stored = S7._load_gz(src) if src is not None and src.is_file() else None
            if stored is not None and stored.get('code') == task_code:
                results.setdefault(kind, {})[key] = stored['result']
                S7._dump_gz(path_of(out_dir, kind, key), stored)
                reused.append(f'{kind}:{key}')
            else:
                futures[pool.submit(run_task, kind, key, a)] = (kind, key)
        for n, fut in enumerate(as_completed(futures), 1):
            try:
                kind, key, out, wall = fut.result()
            except Exception as err:
                errors.append(f'{futures[fut]}: {err!r}')
                say(f'{n:4d}/{len(futures)} FAILED {futures[fut]}: {err!r}')
                continue
            results.setdefault(kind, {})[key] = out
            walls[f'{kind}:{key}'] = wall
            S7._dump_gz(path_of(out_dir, kind, key), dict(code=task_code, result=out))
            say(f'{n:4d}/{len(futures)} {kind:<11}{key:<40}{str(out.get("status", out.get("passed", ""))):<11}{wall:8.1f}s')
        if errors:
            raise SystemExit('tasks failed, nothing assembled; finished tasks are cached for --resume-from: ' + '; '.join(errors))

    series_dir = out_dir/'rut8-series'
    series_dir.mkdir(exist_ok=True)
    inputs = {PROTOCOL: hashlib.sha256(_find(PROTOCOL).read_bytes()).hexdigest()}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        say(f'stage 8, phase {args.phase}: {args.workers} workers, output {out_dir}')
        new_names = TH['P']['new_populations']
        if args.phase == 'predict':
            drain(pool, [(1., 'population', name, (name, dE)) for name, dE in new_names.items()])
            new_pops = {name: results['population'][name] for name in new_names}
            pops, tasks = predict_tasks(new_pops)
            say(f'{len(tasks)} tasks in the first wave')
            drain(pool, tasks)
            drain(pool, second_wave(results, pops))
            gates, part_P = predict(results, list(pops))
            S7._dump_gz(series_dir/'new-populations.json.gz', new_pops)
            rows = sorted(rectangle_rows(results), key=lambda r: (r['name'], r['m'], r['shifted'], r['rect'][2]))
            S7._dump_gz(series_dir/'modes.json.gz', rows, 8)
            files = ['new-populations.json.gz', 'modes.json.gz']
            body = dict(part_L=dict(gates=gates, discretization={k: TH['L'][k] for k in
                                    ('n_L', 'n_u', 'n_eta', 'l_max', 'node_spacing', 'node_margin', 'node_clip', 'contour')}),
                        part_P=part_P,
                        populations={n: {k: v for k, v in p.items() if k != 'state'} for n, p in new_pops.items()})
            name_out = PREDICTIONS
        else:
            pred_path = HERE/PREDICTIONS
            predictions = json.loads(pred_path.read_text(encoding='utf-8'))
            inputs[PREDICTIONS] = hashlib.sha256(pred_path.read_bytes()).hexdigest()
            new_pops = S7._load_gz(HERE/'rut8-series'/'new-populations.json.gz')
            drain(pool, measure_tasks(new_pops))
            S7._dump_gz(series_dir/'V6.json.gz', [results['V6'][k] for k in sorted(results['V6'])], S7.SERIES_DIGITS)
            cells = {}
            for k in sorted(results['X']):
                r = results['X'][k]
                cells.setdefault((r['name'], r['bodies']), []).append(r)
            for (name, n), runs in cells.items():
                S7._dump_gz(series_dir/f'X_{name}_N{n}.json.gz', runs, S7.SERIES_DIGITS)
            files = ['V6.json.gz'] + [f'X_{name}_N{n}.json.gz' for name, n in cells]
            derived = evaluate_measurement(load_series(series_dir), predictions)
            gates = dict(V1_the_measure=dict(rows=results['V1'], passed=bool(all(r['passed'] for r in results['V1'].values()))),
                         V2_drawn_moments=dict(rows=results['V2'], passed=bool(all(r['passed'] for r in results['V2'].values()))),
                         V3_support_fits=dict(rows=results['V3'], passed=bool(all(r['fits'] for r in results['V3'].values()))),
                         V5_drawn_source_writes_the_field=dict(rows=results['V5'],
                                                               passed=bool(all(r['passed'] for r in results['V5'].values()))),
                         V6_draw_is_stationary=derived['V6'],
                         V8_ensemble=dict(all_runs_completed=derived['all_runs_completed'],
                                          pairs_identical=derived['pairs_identical'],
                                          passed=bool(derived['all_runs_completed'] and derived['pairs_identical'])))
            body = dict(part_M=dict(gates=gates, comparison=derived['comparison'], run_index=derived['run_index']))
            name_out = 'rut8-results.json'
    failed = sorted(k for k, g in gates.items() if not g['passed'])
    result = dict(experiment='RUT-1 stage 8: the linear modes of the constructed populations, predicted before they are measured',
                  phase=args.phase, protocol=PROTOCOL,
                  statuses=dict(reproduction='decided by rut8_checks.py against this archive',
                                numerical_verification=dict(passed=not failed, failed_gates=failed),
                                scientific_outcome='archived as it fell; never part of the exit status'),
                  deviations=[], thresholds=TH, **body,
                  code_sha256_at_launch=code, tasks_reused=len(reused), input_sha256=inputs,
                  series_sha256={f: hashlib.sha256((series_dir/f).read_bytes()).hexdigest() for f in sorted(files)})
    (out_dir/name_out).write_text(json.dumps(_safe(result), indent=1) + '\n', encoding='utf-8', newline='\n')
    (out_dir/f'rut8-{args.phase}-log.json').write_text(json.dumps(dict(task_wall_seconds=walls,
                                                       wall_seconds=round(time.time() - t0, 1)), indent=1), encoding='utf-8')
    say(f'numerical verification: {"PASSED" if not failed else "FAILED " + ", ".join(failed)}')
    for k, g in sorted(gates.items()):
        say(f"  {'PASS' if g['passed'] else 'FAIL'} {k}")
    return 0 if not failed else 1


if __name__ == '__main__':
    raise SystemExit(main())
