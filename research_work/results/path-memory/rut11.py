"""RUT-1 stage 11 (protocol-rut11.md): the quiet region.

    part A   the analytic rule: it integrates a known integral where stage 9's quadrature cannot, and it
             reproduces stage 9's roots where stage 9's quadrature is valid
    part L   the low box: the search carried from stage 8's floor of 0.006 down to 5e-4, on the populations
             stage 9 found no root in above that floor
    part M   the mode domain: where the response operator is too small to hold a mode at all
    part T   the family's boundary, bracketed by resolution instead of extrapolated
    part D   a controlled disturbance put into each quiet population, and what the live system does with it

Three statuses are kept separate, as the owner ruled: reproduction and numerical verification set the exit
status; the scientific outcome never does. WHICH populations turn out to be quiet is an outcome, so the gates
below ask only whether the search was sound and whether the apparatus would have seen growth had there been
any. Every threshold is read from the protocol. rut11_checks.py is the suite job.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut11_tasks as TASKS    # noqa: E402
from rut11_tasks import CODE, PROTOCOL, REP, T0, TH, _find, np, run_task   # noqa: E402
import rut7 as S7              # noqa: E402  (series i/o, by import)
import rut10 as S10            # noqa: E402  (the mode series of a run, by import)


# ================================================================ reading a disturbed run
def response_fit(run, lo, hi):
    """The disturbed harmonic's amplitude over a declared window, as a fitted rate and as a plain ratio."""
    t, a, ph, _ = S10.mode_series(run, run['m'])
    sel = (t >= lo) & (t <= hi)
    if sel.sum() < 4:
        return None
    g = np.polyfit(t[sel]*T0, np.log(a[sel]), 1)
    p = np.polyfit(t[sel]*T0, ph[sel], 1)
    resid = np.log(a[sel]) - np.polyval(g, t[sel]*T0)
    return dict(window=[float(lo), float(hi)], records=int(sel.sum()), rate=float(g[0]),
                pattern_speed=float(-p[0]/run['m']), rms_residual=float(np.sqrt(np.mean(resid**2))),
                ratio_over_window=float(a[sel][-1]/a[sel][0]), peak_over_initial=float(a.max()/a[0]),
                final_over_initial=float(a[-1]/a[0]), initial_amplitude=float(a[0]))


def _mean_err(values):
    if not values:
        return None, None
    m = float(np.mean(values))
    e = float(np.std(values, ddof=1)/np.sqrt(len(values))) if len(values) > 1 else None
    return m, e


# ================================================================ gates
def gate_A1(row):
    A = TH['A']['toy']
    rows = row['rows']
    worst = max(r['pair_error'] for r in rows)
    control = max(r['gauss_legendre_error'] for r in rows if r['growth_rate'] <= A['control_below'])
    return dict(rows=rows, worst_pair_error=worst, threshold=A['accuracy'],
                negative_control="stage 9's quadrature of the same integral must be seen to miss it by more "
                                 'than %g at the rates at or below %g' % (A['control_min'], A['control_below']),
                control_value=control, control_min=A['control_min'], control_rejected=bool(control > A['control_min']),
                passed=bool(worst < A['accuracy'] and control > A['control_min']))


def gate_A2(rows, old_rule, boxes):
    """Where stage 9's quadrature is valid the two rules must agree, and the rule must not move when refined.
    The control is stage 9's own quadrature used where its poles lie: it must be seen to disagree."""
    A = TH['A']
    out, ok = {}, True
    for key, r in sorted(rows.items()):
        d, f = r.get('difference_from_stage9'), r.get('difference_under_refinement')
        good = d is not None and d < A['stage9_agreement'] and (f is None or f < A['refinement_agreement'])
        ok &= bool(good)
        out[key] = dict(stage9_root=r['stage9_root'], analytic_roots=r['analytic']['roots'],
                        difference_from_stage9=d, difference_under_refinement=f,
                        max_spread=r['analytic']['max_spread'], harmonics=r['harmonics'], passed=bool(good))
    worst, detail = 0., {}
    for key, c in sorted(old_rule.items()):
        mine = boxes.get(c['matches'])
        if mine is None:
            continue
        theirs = max((complex(*s) for s in c['roots']), key=lambda z: z.real, default=None)
        ours = max((complex(*s) for s in mine['roots']), key=lambda z: z.real, default=None)
        diff = None if (theirs is None or ours is None) else float(abs(theirs - ours)/abs(ours))
        broke = (not c['counts_agree']) or (diff is not None and diff > A['control_min'])
        worst = max(worst, diff or 0.)
        detail[key] = dict(stage9_located=c['located'], stage9_winding=c['winding'],
                           counts_agree=c['counts_agree'], stage9_roots=c['roots'], analytic_roots=mine['roots'],
                           relative_difference=diff, breaks=bool(broke))
    control_ok = bool(detail) and all(v['breaks'] for v in detail.values())
    return dict(rows=out, thresholds=dict(stage9_agreement=A['stage9_agreement'],
                                          refinement_agreement=A['refinement_agreement']),
                negative_control="stage 9's quadrature, used in the low box where its poles lie, must be seen "
                                 'to break the count agreement or misplace the root by more than %g' % A['control_min'],
                control_rows=detail, control_worst_difference=worst, control_rejected=control_ok,
                passed=bool(ok and control_ok))


def gate_A3(boxes, ladder):
    """The error the rule declares about itself, on every contour it was used on."""
    rows = list(boxes.values()) + list(ladder.values())
    worst = max((r['max_spread'] for r in rows), default=None)
    where = max(rows, key=lambda r: r['max_spread']) if rows else None
    return dict(contours=len(rows), worst=None if worst is None else float(worst),
                threshold=REP['spread_max'],
                worst_at=None if where is None else dict(name=where.get('name'), m=where.get('m'), rect=where['rect']),
                passed=bool(worst is not None and worst < REP['spread_max']))


def gate_L1(boxes, quiet_and_control):
    """The search is sound: the count agrees with what is located, the declared error is inside its budget,
    and the poles the analytic rule did not remove are the declared margin away."""
    bad = {k: dict(winding=r['winding'], located=r['located'], max_spread=r['max_spread'],
                   pole_margin=r['pole_margin'])
           for k, r in boxes.items()
           if not r['counts_agree'] or r['max_spread'] >= REP['spread_max']
           or r['pole_margin'] < REP['pole_margin_min']}
    covered = {}
    for name in quiet_and_control:
        rows = [r for r in boxes.values() if r['name'] == name]
        covered[name] = dict(rectangles=len(rows), m_values=sorted({r['m'] for r in rows}),
                             partitions=sorted({bool(r['shifted']) for r in rows}),
                             located=int(sum(r['located'] for r in rows)),
                             max_spread=max([r['max_spread'] for r in rows], default=None),
                             min_pole_margin=min([r['pole_margin'] for r in rows], default=None))
    return dict(rectangles=len(boxes), coverage=covered, failures=bad,
                thresholds=dict(spread_max=REP['spread_max'], pole_margin_min=REP['pole_margin_min']),
                note='whether a population holds a root is a reading, not this gate', passed=not bad)


def gate_L2(boxes):
    """The search can see a root that is there: the one the prototype found below stage 8's floor, in both
    partitions. This gate IS the detectability control every quiet reading rests on."""
    L = TH['L']
    rows = [r for r in boxes.values() if r['name'] == L['control_population'] and r['located'] > 0]
    found = {}
    for part in (False, True):
        hits = [complex(*s) for r in rows if r['shifted'] == part for s in r['roots']]
        found[('shifted' if part else 'base')] = None if not hits else max(hits, key=lambda z: z.real)
    a, b = found['base'], found['shifted']
    agree = None if (a is None or b is None) else float(abs(a - b)/abs(a))
    ok = a is not None and b is not None and agree < L['control_agreement'] \
        and REP['gamma_min'] < a.real < REP['gamma_stage8']
    return dict(population=L['control_population'],
                root_base=None if a is None else [a.real, a.imag],
                root_shifted=None if b is None else [b.real, b.imag],
                partitions_agree_to=agree, threshold=L['control_agreement'],
                e_folding_periods=None if a is None else float(1./(a.real*T0)),
                between_the_floors=bool(a is not None and REP['gamma_min'] < a.real < REP['gamma_stage8']),
                passed=bool(ok))


def gate_M1(norms):
    """The grid the operator norm is maximised on is adequate, and the diagnostic discriminates."""
    M = TH['M']
    rows = {(r['name'], r['m']): r for r in norms.values() if not r['refined']}
    fine = {(r['name'], r['m']): r for r in norms.values() if r['refined']}
    rises = {f'{k[0]}:m{k[1]}': float(fine[k]['max_norm'] - rows[k]['max_norm']) for k in fine if k in rows}
    worst = max(rises.values(), default=0.)
    ctrl = rows.get((M['control_population'], M['control_m']))
    cmax = float(ctrl['max_norm']) if ctrl else 0.
    return dict(refinement_rises=rises, worst_rise=worst, threshold=M['refine_rise'],
                by_population={f'{k[0]}:m{k[1]}': v['max_norm'] for k, v in sorted(rows.items())},
                negative_control='%s at m = %d, which holds a root, must be seen to reach 1'
                                 % (M['control_population'], M['control_m']),
                control_value=cmax, control_rejected=bool(cmax >= 1.),
                passed=bool(worst < M['refine_rise'] and cmax >= 1.))


def gate_T1(ladder):
    T = TH['T']
    out, ok = {}, True
    for key, r in sorted(ladder.items()):
        moved = r.get('difference_under_refinement')
        good = bool(r['counts_agree'] and r['max_spread'] < REP['spread_max']
                    and (r['located'] == 0 or (moved is not None and moved < T['refinement_agreement'])))
        ok &= good
        out[key] = dict(dE=r['dE'], located=r['located'], winding=r['winding'], roots=r['roots'],
                        difference_under_refinement=moved, max_spread=r['max_spread'], passed=good)
    return dict(rows=out, threshold=T['refinement_agreement'], passed=bool(ok))


def gate_D1(runs):
    """The disturbance is the declared one, and the sample the bodies were drawn as carries no source of its
    own at the harmonic being disturbed -- so what is watched can only be what was put there."""
    D = TH['D']
    out, ok = {}, True
    for key, r in sorted(runs.items()):
        if r['variant'] == 'ordinary_draw':
            continue
        t, a, ph, C = S10.mode_series(r, r['m'])
        own = float(r['source_harmonics'][r['m'] - 1])
        good = own < D['other_harmonic_max'] and (r['target'] == 0.
                                                  or abs(a[0]/r['target'] - 1.) < D['amplitude_agreement'])
        ok &= bool(good)
        out[key] = dict(initial_amplitude=float(a[0]), target=r['target'], own_source_at_m=own,
                        source_harmonics=r['source_harmonics'][:4], bodies=r['bodies'], fold=r['fold'],
                        passed=bool(good))
    worst = max((float(v['source_harmonics'][v['m'] - 1]) for v in runs.values()
                 if v['variant'] == 'ordinary_draw'), default=0.)
    return dict(rows=out, thresholds=dict(amplitude_agreement=D['amplitude_agreement'],
                                          other_harmonic_max=D['other_harmonic_max']),
                negative_control='an ordinary draw of the same size must be seen to carry a source of its own '
                                 'at that harmonic above %g' % D['control_harmonic_min'],
                control_value=worst, control_rejected=bool(worst > D['control_harmonic_min']),
                passed=bool(ok and worst > D['control_harmonic_min']))


def control_fit(run):
    """The growing control is read over its growth phase, not over the quiet populations' fixed window: it
    starts from a kick, not from its eigenmode, and it reaches the nonlinear ceiling. From the first record
    above the declared factor times the initial amplitude to the last below the ceiling -- stage 8's rule."""
    D = TH['D']
    t, a, ph, _ = S10.mode_series(run, run['m'])
    above = np.nonzero(a > D['control_window_lo_factor']*a[0])[0]
    if not len(above):
        return dict(window=None, note='the amplitude never rose above the lower bound')
    i0 = int(above[0])
    over = np.nonzero(a[i0:] > D['control_ceiling'])[0]
    i1 = i0 + int(over[0]) if len(over) else len(a)
    if t[i1 - 1] - t[i0] < D['control_window_shortest']:
        return dict(window=[float(t[i0]), float(t[i1 - 1])], note='the window is shorter than the declared minimum')
    sel = slice(i0, i1)
    g = np.polyfit(t[sel]*T0, np.log(a[sel]), 1)
    p = np.polyfit(t[sel]*T0, ph[sel], 1)
    resid = np.log(a[sel]) - np.polyval(g, t[sel]*T0)
    return dict(window=[float(t[i0]), float(t[i1 - 1])], records=int(i1 - i0), rate=float(g[0]),
                pattern_speed=float(-p[0]/run['m']), rms_residual=float(np.sqrt(np.mean(resid**2))),
                e_foldings=float(np.log(a[i1 - 1]/a[i0])), initial_amplitude=float(a[0]))


def gate_D2(runs):
    """The apparatus would see growth: the population with a known mode, given the same kind of disturbance,
    grows at the rate stage 9 recalculated. This gate IS the control the quiet readings rest on."""
    D = TH['D']
    rows = {k: r for k, r in sorted(runs.items())
            if r['name'] == D['control_population'] and r['target'] != 0. and r['variant'] == 'base'}
    fits = {k: control_fit(r) for k, r in rows.items()}
    rates = [f['rate'] for f in fits.values() if f and 'rate' in f]
    mean, err = _mean_err(rates)
    pred = D['control_predicted_rate']
    good = bool(rates) and len(rates) == len(rows)         and abs(mean - pred) <= max(D['control_agreement']*pred, 3*(err or 0.))
    return dict(rows=fits, mean_rate=mean, standard_error=err, predicted_rate=pred,
                relative_difference=None if mean is None else float((mean - pred)/pred),
                agreement=D['control_agreement'], realizations=len(rates), runs=len(rows), passed=bool(good))


def gate_D3(runs):
    bad = {k: dict(status=r['status'], reached=r['t_final_periods'], horizon=r['horizon'])
           for k, r in runs.items() if r['status'] != 'completed' or r['t_final_periods'] < r['horizon'] - 1e-9}
    return dict(runs=len(runs), incomplete=bad, passed=not bad)


# ================================================================ the readings
def reading_quiet(boxes, norms, runs):
    D, L, M = TH['D'], TH['L'], TH['M']
    out = {}
    for name in L['quiet'] + [L['control_population']]:
        rows = [r for r in boxes.values() if r['name'] == name]
        found = [s for r in rows for s in r['roots']]
        bounded = sorted({r['m'] for r in norms.values() if r['name'] == name and not r['refined']
                          and r['max_norm'] < M['norm_bound'] and r['m'] >= M['bounded_from']})
        fits = [response_fit(r, *D['window']) for k, r in sorted(runs.items())
                if r['name'] == name and r['target'] != 0. and r['variant'] == 'base' and r['m'] == D['m_primary']]
        fits = [f for f in fits if f]
        rate, err = _mean_err([f['rate'] for f in fits])
        ratio, _ = _mean_err([f['final_over_initial'] for f in fits])
        peak, _ = _mean_err([f['peak_over_initial'] for f in fits])
        out[name] = dict(
            searched_m=sorted({r['m'] for r in rows}), floor=REP['gamma_min'],
            e_folding_bound_periods=float(1./(REP['gamma_min']*T0)),
            rectangles=len(rows), roots_found=found,
            norm_bounded_m=bounded,
            max_norm_where_bounded=max([r['max_norm'] for r in norms.values() if r['name'] == name
                                        and not r['refined'] and r['m'] in bounded], default=None),
            other_windows={f'{lo}-{hi}': _mean_err([response_fit(r, lo, hi)['rate'] for k, r in sorted(runs.items())
                                                    if r['name'] == name and r['target'] != 0.
                                                    and r['variant'] == 'base' and r['m'] == D['m_primary']
                                                    and response_fit(r, lo, hi)])[0]
                           for lo, hi in D['windows_reported']},
            disturbance=dict(realizations=len(fits), m=D['m_primary'], target=D['target'], horizon=D['horizon'],
                             window=D['window'], mean_rate=rate, standard_error=err,
                             mean_final_over_initial=ratio, mean_peak_over_initial=peak))
    return out


def reading_threshold(ladder):
    T = TH['T']
    rows = []
    for r in sorted(ladder.values(), key=lambda d: d['dE']):
        fastest = max(r['roots'], key=lambda s: s[0]) if r['roots'] else None
        rows.append(dict(dE=r['dE'], root=fastest, rate=None if fastest is None else fastest[0],
                         e_folding_periods=r.get('e_folding_periods'), located=r['located']))
    unstable = [r for r in rows if r['rate'] is not None]
    last = max(unstable, key=lambda r: r['dE']) if unstable else None
    quiet = [r for r in rows if r['rate'] is None and (last is None or r['dE'] > last['dE'])]
    first = min(quiet, key=lambda r: r['dE']) if quiet else None
    return dict(ladder=rows, floor=REP['gamma_min'], last_resolved_unstable=last,
                first_with_no_root_above_the_floor=first,
                bracket=None if not (last and first) else [last['dE'], first['dE']],
                stage8_extrapolations=T['stage8_extrapolations'],
                note='a bracket by resolution, not an extrapolation; it says nothing about growth below the '
                     'floor, and a damped mode is invisible to it')


def reading_variants(runs):
    D = TH['D']
    base = [response_fit(r, *D['window']) for k, r in sorted(runs.items())
            if r['name'] == D['variant_population'] and r['variant'] == 'base' and r['realization'] == 0
            and r['target'] == D['target'] and r['m'] == D['m_primary']]
    ref = base[0] if base and base[0] else None
    out = {}
    for key, r in sorted(runs.items()):
        if r['name'] != D['variant_population'] or r['variant'] == 'base' or r['target'] == 0.:
            continue
        f = response_fit(r, *D['window'])
        out[key] = dict(variant=r['variant'], target=r['target'], fit=f,
                        rate_difference=None if not (f and ref) else float(f['rate'] - ref['rate']),
                        ratio_difference=None if not (f and ref) else
                        float(f['final_over_initial'] - ref['final_over_initial']))
    return dict(reference=ref, rows=out)


# ================================================================ the campaign
def disturbance_jobs():
    D = TH['D']
    jobs = []
    for name in D['quiet']:
        for k in range(D['realizations']):
            jobs.append((name, D['m_primary'], D['quarter'], k, D['target'], 'base', D['horizon']))
    for k in range(D['realizations']):        # the control is kicked far smaller: it grows, and must stay
        jobs.append((D['control_population'], D['m_primary'], D['quarter'], k,   # under the nonlinear ceiling
                     D['control_target'], 'base', D['horizon']))
    for m in D['odd_m']:
        for k in range(D['odd_realizations']):
            jobs.append((D['odd_population'], m, D['quarter'], k, D['target'], 'base', D['horizon']))
    for variant in D['variants']:
        jobs.append((D['variant_population'], D['m_primary'], D['quarter'], 0, D['target'], variant, D['horizon']))
    jobs.append((D['variant_population'], D['m_primary'], D['quarter'], 0, D['variant_target_small'], 'base',
                 D['horizon']))
    jobs.append((D['variant_population'], D['m_primary'], D['quarter'], 0, 0., 'base', D['horizon']))
    return jobs


def first_wave():
    # the live runs are submitted first: they are the long pole, and the pool hands work out in this order
    jobs = []
    for name, m, quarter, k, target, variant, horizon in disturbance_jobs():
        jobs.append(('disturb', f'D:{name}:m{m}:q{quarter}:k{k}:t{target:+.0e}:{variant}',
                     (name, m, quarter, k, target, variant, horizon)))
    for name in TH['A']['roots']:
        jobs.append(('A2', f'A2:{name}', (name, 2, TH['A']['rect'])))
    for name in TH['L']['quiet'] + [TH['M']['control_population']]:
        for m in TH['M']['m_values']:
            jobs.append(('norm', f'norm:{name}:m{m}', (name, m, False)))
            if m in TH['M']['refined_m']:
                jobs.append(('norm', f'norm:{name}:m{m}:refined', (name, m, True)))
    for dE in TH['T']['dE_values']:
        jobs.append(('threshold', f'T:{dE:.5f}', (dE, 2)))
    jobs.append(('A1', 'A1', ()))
    return jobs


def box_key(name, m, shifted, i):
    return f"box:{name}:m{m}:{'s' if shifted else 'b'}{i}"


def second_wave():
    jobs = []
    for name in TH['L']['quiet'] + [TH['L']['control_population']]:
        for m in TH['L']['m_values']:
            for shifted in (False, True):
                for i in range(len(TASKS.strip_rectangles(m, shifted))):
                    jobs.append(('box', box_key(name, m, shifted, i), (name, m, i, shifted)))
    for i in TH['L']['control_rectangles']:
        jobs.append(('box_control', f'oldrule:{i}', (TH['L']['control_population'], 2, i, False)))
    return jobs


def _hash_code():
    return {f: hashlib.sha256(_find(f).read_bytes()).hexdigest() for f in CODE}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument('--resume-from')
    ap.add_argument('--only', default='')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cache = Path(args.resume_from) if args.resume_from else out/'cache'
    cache.mkdir(parents=True, exist_ok=True)
    code = _hash_code()
    stamp = hashlib.sha256((json.dumps(code, sort_keys=True)
                            + TASKS.protocol_path().read_text(encoding='utf-8')).encode()).hexdigest()[:16]
    results, t0, reused = {}, time.time(), 0

    errors = []

    def run(jobs, label):
        nonlocal reused
        todo = []
        for kind, key, a in jobs:
            f = cache/f'{stamp}-{kind}-{hashlib.sha256(key.encode()).hexdigest()[:20]}.json.gz'
            if f.is_file():
                results[key] = S7._load_gz(f)
                reused += 1
            else:
                todo.append((kind, key, a, f))
        print(f'{label}: {len(jobs)} tasks, {len(jobs) - len(todo)} reused, {len(todo)} to run', flush=True)
        if not todo:
            return
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_task, kind, key, a): (key, f) for kind, key, a, f in todo}
            done = 0
            for fut in as_completed(futures):
                key, f = futures[fut]
                done += 1
                try:
                    _, _, value, secs = fut.result()
                except Exception as err:                      # a failed task must not be assembled around
                    errors.append(f'{key}: {err!r}')
                    print(f'  [{done}/{len(todo)}] FAILED {key}: {err!r}', flush=True)
                    continue
                results[key] = value
                S7._dump_gz(f, value)
                print(f'  [{done}/{len(todo)}] {key} ({secs:.0f} s, {time.time() - t0:.0f} s elapsed)', flush=True)
        if errors:
            raise SystemExit('tasks failed, nothing assembled; finished tasks are cached for --resume-from: '
                             + '; '.join(errors))

    only = set(args.only.split(',')) if args.only else None
    run([j for j in first_wave() if only is None or j[0] in only], 'first wave')
    run([j for j in second_wave() if only is None or j[0] in only], 'second wave')

    boxes = {k: v for k, v in results.items() if k.startswith('box:')}
    old_rule = {k: v for k, v in results.items() if k.startswith('oldrule:')}
    for k, v in old_rule.items():
        v['matches'] = box_key(v['name'], v['m'], v['shifted'], v['index'])
    norms = {k: v for k, v in results.items() if k.startswith('norm:')}
    ladder = {k: v for k, v in results.items() if k.startswith('T:')}
    runs = {k: v for k, v in results.items() if k.startswith('D:')}
    a2 = {k: v for k, v in results.items() if k.startswith('A2:')}
    quiet_and_control = TH['L']['quiet'] + [TH['L']['control_population']]

    gates = dict(A1_the_rule_on_a_known_integral=gate_A1(results['A1']),
                 A2_the_rule_agrees_with_stage9_where_stage9_is_valid=gate_A2(a2, old_rule, boxes),
                 A3_the_declared_error_on_every_contour=gate_A3(boxes, ladder),
                 L1_the_search_is_sound=gate_L1(boxes, quiet_and_control),
                 L2_the_search_can_see_a_root_below_the_old_floor=gate_L2(boxes),
                 M1_the_operator_norm_grid_is_adequate=gate_M1(norms),
                 T1_every_rung_survives_refinement=gate_T1(ladder),
                 D1_the_disturbance_is_as_declared=gate_D1(runs),
                 D2_the_apparatus_would_see_growth=gate_D2(runs),
                 D3_every_run_completed=gate_D3(runs))
    numerical = all(g['passed'] for g in gates.values())
    archive = dict(experiment='RUT-1 stage 11', protocol=PROTOCOL,
                   protocol_read_from=str(TASKS.protocol_path()),
                   protocol_sha256=hashlib.sha256(TASKS.protocol_path().read_bytes()).hexdigest(),
                   statuses=dict(numerical_verification='passed' if numerical else 'failed',
                                 scientific_outcome='archived as it fell; it never sets the exit status'),
                   thresholds=TH, gates=gates,
                   readings=dict(quiet_region=reading_quiet(boxes, norms, runs),
                                 threshold_bracket=reading_threshold(ladder),
                                 disturbance_variants=reading_variants(runs)),
                   code_sha256_at_launch=code, tasks_reused=reused, seconds=round(time.time() - t0, 1))
    series = out/'rut11-series'
    series.mkdir(exist_ok=True)
    for label, group in (('boxes', {**boxes, **old_rule}), ('norms', norms), ('threshold', ladder),
                         ('disturbance', runs), ('rule', {**a2, 'A1': results['A1']})):
        S7._dump_gz(series/f'{label}.json.gz', {k: group[k] for k in sorted(group)})
    archive['series_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(series.iterdir())}
    (out/'rut11-results.json').write_text(json.dumps(archive, indent=1, sort_keys=True), encoding='utf-8')
    print('numerical verification:', archive['statuses']['numerical_verification'])
    for k, g in gates.items():
        print(f"  {'pass' if g['passed'] else 'FAIL'}  {k}")
    return 0 if numerical else 1


if __name__ == '__main__':
    sys.exit(main())
