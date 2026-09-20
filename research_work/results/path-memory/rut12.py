"""RUT-1 stage 12 (protocol-rut12.md): formation from an empty field, with both reciprocal budgets.

    part A   the angular-momentum ledger this stage derives: it closes, it converges, and L_field is a
             property of the field rather than of the box it is computed in
    part G   formation: an equilibrium of the bare point mass, an empty field, writing switched on, and both
             ledgers carried through the transition
    the reading   what forms, against a prediction fixed before any run and against the verified population

Three statuses are kept separate, as the owner ruled: reproduction and numerical verification set the exit
status; the scientific outcome never does. Every threshold is read from the protocol.
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
import rut12_tasks as TASKS     # noqa: E402
from rut12_tasks import CODE, PROTOCOL, TH, _find, np, run_task   # noqa: E402
import rut7 as S7               # noqa: E402  (series i/o, by import)

DRIVERS = ('rut12.py', 'rut11.py', 'rut10.py', 'rut9.py', 'rut8.py', 'rut7.py')
KEYS = ('mean_L', 'sigma_L', 'mean_rc', 'sigma_rc', 'rms_epicycle')


def _mean_err(values):
    if not values:
        return None, None
    m = float(np.mean(values))
    return m, (float(np.std(values, ddof=1)/np.sqrt(len(values))) if len(values) > 1 else None)


def settled(run):
    """The structure averaged over the records from the declared instant to the horizon."""
    rows = [r for r in run['rows'] if r['t_periods'] >= TH['G']['settled_from']]
    if not rows:
        return None
    out = {k: float(np.mean([r['structure'][k] for r in rows])) for k in run['rows'][0]['structure']}
    out['support'] = float(np.mean([r['support'] for r in rows]))
    out['records'] = len(rows)
    return out


def worst_closure(run):
    return (max(abs(r['energy_balance']) for r in run['rows']),
            max(abs(r['angular_balance']) for r in run['rows']))


# ================================================================ gates
def gate_A1(budgets):
    A = TH['A']
    ref = {r['step']: r for r in budgets.values()
           if r['offset'] == 0. and [r['box'], r['modes']] == list(A['configurations'][0])}
    steps = sorted(ref, reverse=True)
    orders = {}
    for a, b in zip(steps, steps[1:]):
        for kind in ('energy_balance', 'angular_balance'):
            orders[f'{kind}:{a}->{b}'] = float(abs(ref[a][kind])/max(abs(ref[b][kind]), 1e-300))
    fine = ref.get(A['reference_step'])
    # a convergence ratio means nothing once the residual is at round-off, so the order is required only of
    # pairs whose finer member is still above the declared floor; which pairs those are is recorded.
    judged = {k: v for k, v in orders.items()
              if abs(ref[float(k.split('->')[1])][k.split(':')[0]]) > A['order_min_residual']}
    ok = (fine is not None and abs(fine['energy_balance']) < A['energy_closure']
          and abs(fine['angular_balance']) < A['angular_closure']
          and bool(judged)
          and all(A['order_range'][0] <= v <= A['order_range'][1] for v in judged.values()))
    return dict(by_step={str(s): dict(energy=ref[s]['energy_balance'], angular=ref[s]['angular_balance'])
                         for s in steps},
                observed_orders=orders, orders_judged=sorted(judged), thresholds=dict(energy=A['energy_closure'], angular=A['angular_closure'],
                                                        order_range=A['order_range']),
                passed=bool(ok))


def gate_A2(budgets):
    A = TH['A']
    # sorted, not in completion order: an average over tasks must not depend on which worker was quickest
    plain = [r for _, r in sorted(budgets.items()) if r['offset'] == 0. and r['step'] == A['reference_step']]
    vals = {f"box{r['box']}:modes{r['modes']}": r['L_field'] for r in plain}
    ref = float(np.mean([vals[k] for k in sorted(vals)]))
    spread = max(abs(v - ref)/abs(ref) for v in vals.values()) if vals else float('inf')
    worst = max(abs(r['angular_balance']) for r in plain) if plain else float('inf')
    ctrl = [r for _, r in sorted(budgets.items()) if r['offset'] != 0.]
    c_close = max((abs(r['angular_balance']) for r in ctrl), default=0.)
    c_field = max((abs(r['L_field'] - ref)/abs(ref) for r in ctrl), default=0.)
    rejected = c_close > A['control_closure_min'] and c_field > A['control_field_min']
    return dict(L_field=vals, relative_spread=float(spread), worst_angular_closure=float(worst),
                thresholds=dict(invariance=A['field_invariance'], angular=A['angular_closure']),
                negative_control='the position grid displaced by half a box must break the closure and move '
                                 'L_field',
                control_closure=float(c_close), control_field_shift=float(c_field), control_rejected=bool(rejected),
                passed=bool(spread < A['field_invariance'] and worst < A['angular_closure'] and rejected))


def gate_A3(prediction):
    """The spectral field at the declared resolution must be the field the population was built with, or
    neither the mature control nor the comparison below means anything."""
    c = prediction['spectral_field_against_the_population']
    lim = TH['A']['representation_agreement']
    return dict(**c, threshold=lim, passed=bool(c['max_relative'] < lim))


def gate_G1(runs):
    G = TH['G']
    bad = {}
    for k, r in sorted(runs.items()):
        e, a = worst_closure(r)
        if e >= G['energy_closure'] or a >= G['angular_closure'] or not r['rows'] \
                or r['rows'][-1]['t_periods'] < r['horizon'] - 1e-6:
            bad[k] = dict(worst_energy=e, worst_angular=a,
                          reached=r['rows'][-1]['t_periods'] if r['rows'] else None, horizon=r['horizon'])
    worst = {k: dict(zip(('energy', 'angular'), worst_closure(r))) for k, r in sorted(runs.items())}
    return dict(runs=len(runs), worst=worst, failures=bad,
                thresholds=dict(energy=G['energy_closure'], angular=G['angular_closure']), passed=not bad)


def _drift(run):
    a, b = run['rows'][0]['structure'], settled(run)
    return max(abs(b[k] - a[k])/max(abs(a[k]), 1e-300) for k in KEYS)


def gate_G2(runs):
    """With writing frozen the start must hold: whatever the other runs do is caused by the writing."""
    G = TH['G']
    r = next((v for v in runs.values() if v['variant'] == 'no_writing'), None)
    d = None if r is None else _drift(r)
    return dict(present=r is not None, drift=d, threshold=G['no_writing_drift'],
                passed=bool(r is not None and d < G['no_writing_drift']))


def gate_G3(runs):
    G = TH['G']
    r = next((v for v in runs.values() if v['start'] == 'mature'), None)
    d = None if r is None else _drift(r)
    return dict(present=r is not None, drift=d, threshold=G['mature_drift'],
                passed=bool(r is not None and d < G['mature_drift']))


def gate_G4(runs):
    G = TH['G']
    base = next((v for v in runs.values()
                 if v['variant'] == 'base' and v['start'] == G['start'] and v['seed'] == G['seeds'][0]), None)
    out, ok = {}, base is not None
    if base is not None:
        b = settled(base)
        for k, r in sorted(runs.items()):
            if r['variant'] in ('base', 'no_writing') or r['start'] != G['start']:
                continue
            s = settled(r)
            if s is None:
                out[k], ok = dict(moved=None, passed=False), False
                continue
            moved = max(abs(s[q] - b[q])/max(abs(b[q]), 1e-300) for q in KEYS)
            good = moved < G['refinement_agreement'] or r['variant'] == 'long_horizon'
            ok &= good
            out[k] = dict(moved=float(moved), passed=bool(good),
                          note='reported, not bounded' if r['variant'] == 'long_horizon' else '')
    return dict(rows=out, threshold=G['refinement_agreement'], passed=bool(ok))


# ================================================================ the readings
def reading_formation(runs, prediction):
    G = TH['G']
    base = [settled(r) for k, r in sorted(runs.items())
            if r['variant'] == 'base' and r['start'] == G['start']]
    base = [b for b in base if b]
    pred, tgt = prediction['predicted'], prediction['target']
    rows = {}
    for k in KEYS + ('support',):
        vals = [b[k] for b in base]
        m, e = _mean_err(vals)
        p = pred.get(k)
        t = tgt.get(k) if k in tgt else None
        rows[k] = dict(formed=m, standard_error=e, predicted=p, target=t,
                       formed_over_predicted=(None if not (m and p) else float(m/p)),
                       formed_over_target=(None if not (m and t) else float(m/t)))
    return dict(realizations=len(base), settled_from=G['settled_from'], by_quantity=rows,
                note='the prediction was fixed in the protocol before any run; neither it nor the target is a '
                     'gate')


def reading_budget(runs):
    G = TH['G']
    out = {}
    for k, r in sorted(runs.items()):
        if r['variant'] != 'base' or r['start'] != G['start']:
            continue
        last = r['rows'][-1]
        out[k] = dict(dissipated=last['dissipated'], torque_integral=last['torque_integral'],
                      L_field_final=last['L_field'], L_matter_final=last['L_matter'],
                      L_total_initial=last['L_total_0'],
                      angular_through_the_field=float(last['L_field'] + last['torque_integral']),
                      energy_balance=last['energy_balance'], angular_balance=last['angular_balance'])
    return out


# ================================================================ the campaign
def formation_jobs():
    G = TH['G']
    jobs = [(G['population'], G['start'], s, 'base') for s in G['seeds']]
    jobs += [(G['population'], G['start'], G['seeds'][0], v) for v in G['variants']]
    jobs += [(G['population'], 'mature', G['seeds'][0], 'base')]
    jobs += [(G['population'], G['start'], G['seeds'][0], 'no_writing')]
    return jobs


def budget_jobs():
    A = TH['A']
    box, modes = A['configurations'][0]
    jobs = [(s, box, modes, 0., A['orbits']) for s in A['steps']]
    jobs += [(A['reference_step'], b, m, 0., A['orbits']) for b, m in A['configurations'][1:]]
    jobs += [(A['reference_step'], box, modes, A['control_offset'], A['orbits'])]
    return jobs


def first_wave():
    jobs = []
    for name, start, seed, variant in formation_jobs():     # the long runs first
        jobs.append(('formation', f'F:{name}:{start}:s{seed}:{variant}', (name, start, seed, variant)))
    for a in budget_jobs():
        jobs.append(('budget', 'A:h%g:b%g:m%d:o%g' % (a[0], a[1], a[2], a[3]), a))
    jobs.append(('prediction', 'P', (TH['P']['population'], TH['P']['bodies'], TH['P']['seed'])))
    return jobs


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
    code = {f: hashlib.sha256(_find(f).read_bytes()).hexdigest() for f in CODE}
    # the cache is keyed on the code that PRODUCES a task result, not on this driver: correcting an assembly
    # or a gate must not throw away five hours of finished work.
    task_code = {f: h for f, h in code.items() if f not in DRIVERS}
    stamp = hashlib.sha256((json.dumps(task_code, sort_keys=True)
                            + TASKS.protocol_path().read_text(encoding='utf-8')).encode()).hexdigest()[:16]
    results, t0, reused, errors = {}, time.time(), 0, []

    jobs = [j for j in first_wave() if not args.only or j[0] in args.only.split(',')]
    todo = []
    for kind, key, a in jobs:
        f = cache/f'{stamp}-{kind}-{hashlib.sha256(key.encode()).hexdigest()[:20]}.json.gz'
        if f.is_file():
            results[key] = S7._load_gz(f)
            reused += 1
        else:
            todo.append((kind, key, a, f))
    print(f'{len(jobs)} tasks, {reused} reused, {len(todo)} to run', flush=True)
    if todo:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_task, k, key, a): (key, f) for k, key, a, f in todo}
            for n, fut in enumerate(as_completed(futures), 1):
                key, f = futures[fut]
                try:
                    _, _, value, secs = fut.result()
                except Exception as err:
                    errors.append(f'{key}: {err!r}')
                    print(f'  [{n}/{len(todo)}] FAILED {key}: {err!r}', flush=True)
                    continue
                results[key] = value
                S7._dump_gz(f, value)
                print(f'  [{n}/{len(todo)}] {key} ({secs:.0f} s, {time.time() - t0:.0f} s elapsed)', flush=True)
        if errors:
            raise SystemExit('tasks failed, nothing assembled; finished tasks are cached for --resume-from: '
                             + '; '.join(errors))

    runs = {k: v for k, v in results.items() if k.startswith('F:')}
    budgets = {k: v for k, v in results.items() if k.startswith('A:')}
    prediction = results['P']
    gates = dict(A1_both_ledgers_close_and_converge=gate_A1(budgets),
                 A2_the_field_angular_momentum_belongs_to_the_field=gate_A2(budgets),
                 A3_the_spectral_field_is_the_population_s=gate_A3(prediction),
                 G1_every_run_is_sound_throughout=gate_G1(runs),
                 G2_the_start_is_an_equilibrium=gate_G2(runs),
                 G3_the_target_is_a_fixed_point=gate_G3(runs),
                 G4_the_settled_state_survives_refinement=gate_G4(runs))
    numerical = all(g['passed'] for g in gates.values())
    archive = dict(experiment='RUT-1 stage 12', protocol=PROTOCOL,
                   protocol_read_from=str(TASKS.protocol_path()),
                   protocol_sha256=hashlib.sha256(TASKS.protocol_path().read_bytes()).hexdigest(),
                   statuses=dict(numerical_verification='passed' if numerical else 'failed',
                                 scientific_outcome='archived as it fell; it never sets the exit status'),
                   thresholds=TH, gates=gates, prediction=prediction,
                   readings=dict(what_formed=reading_formation(runs, prediction),
                                 budgets_through_the_transition=reading_budget(runs)),
                   code_sha256_at_launch=code, tasks_reused=reused, seconds=round(time.time() - t0, 1))
    series = out/'rut12-series'
    series.mkdir(exist_ok=True)
    for label, group in (('formation', runs), ('budget', budgets), ('prediction', {'P': prediction})):
        S7._dump_gz(series/f'{label}.json.gz', {k: group[k] for k in sorted(group)})
    archive['series_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(series.iterdir())}
    (out/'rut12-results.json').write_text(json.dumps(archive, indent=1, sort_keys=True), encoding='utf-8')
    print('numerical verification:', archive['statuses']['numerical_verification'])
    for k, g in gates.items():
        print(f"  {'pass' if g['passed'] else 'FAIL'}  {k}")
    return 0 if numerical else 1


if __name__ == '__main__':
    sys.exit(main())
