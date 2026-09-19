"""RUT-1 stage 9 (protocol-rut9.md): the corrected response calculation, and what B13's shortfall is not.

    part B   the boundary terms of the sharp truncation, which stage 8 omitted: reproduced without them,
             verified against a resolved taper with them
    part N   the numerical representation: every ring-limit annulus resolved ON ITS OWN before the ring limit
             is read; the same refinements on the populations; the root-free conclusions re-swept
    part F   stage 8's forecasts recalculated with the corrected representation -- BESIDE the originals, which
             stay the forecasts of record
    part S   what the root is sensitive to
    part E   B13's own eigenfunction imposed at its predicted complex rate on bodies of the simulated population

Three statuses are kept separate, as the owner ruled: reproduction and numerical verification set the exit
status; the scientific outcome never does. Every threshold is read from the protocol. Stage 8's archive is not
touched: its L3 and L3b stay failed there whatever this stage finds. rut9_checks.py is the suite job.
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
import rut9_tasks as TASKS     # noqa: E402
from rut9_tasks import CODE, PROTOCOL, REP, TH, _find, _safe, np, run_task   # noqa: E402
import rut7 as S7              # noqa: E402  (series i/o, by import)

TASK_CODE = tuple(f for f in CODE if f not in ('rut9.py', 'rut8.py', 'rut7.py'))
STAGE8 = 'rut8-predictions.json'
KINDS = dict(A='solved for its declared mean support', B="held at family B's coupling")


# ================================================================ what stage 8 archived
def stage8():
    return json.loads(_find(STAGE8).read_text(encoding='utf-8'))


def stage8_root(pred, name):
    """The m = 2 root stage 8 archived for a population or a family member; None where it found none."""
    P = pred['part_P']
    if name in P['new_populations']:
        rows = P['new_populations'][name]['every_m']['m2']
    elif name.startswith('B_dE'):
        member = [r for r in P['family_B'] if abs(r['dE'] - float(name[4:])) < 1e-9][0]
        rows = member['roots']
    else:
        rows = P['spectrum'][f'{name}:m2']['roots'] if isinstance(P['spectrum'][f'{name}:m2'], dict) else P['spectrum'][f'{name}:m2']
    rows = sorted(rows, key=lambda r: -r['re'])
    return [rows[0]['re'], rows[0]['im']] if rows else None


def population_facts(name, pops):
    if name.startswith('B_dE'):
        return 'B', float(name[4:])
    return pops[name]['family'], pops[name]['dE']


# ================================================================ gates, from finished tasks
def _rel(a, b):
    a, b = complex(*a), complex(*b)
    return float(abs(a - b)/abs(b))


def gate_B1(rows):
    B = TH['B']
    out = {}
    for name, r in rows.items():
        ok = (r['orbit_arrays_identical'] and r['nodes_identical'] and r['T_difference'] <= B['reproduction_T']
              and r['root_difference_from_archive'] <= B['reproduction_root'])
        out[name] = dict(r, passed=bool(ok), control_rejected=bool(r['T_difference_with_boundary_rows'] > B['reproduction_control_min']))
    return dict(rows=out, negative_control='the same comparison WITH the boundary rows must be seen to differ from stage 8',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def gate_B2(rows):
    B = TH['B']
    out = {}
    for key, r in rows.items():
        falling = all(a['with_rows'] > b['with_rows'] for a, b in zip(r['rows'][:-1], r['rows'][1:]))
        out[key] = dict(r, falling=bool(falling), passed=bool(falling and r['extrapolated_with_rows'] < B['taper_agreement']),
                        control_rejected=bool(r['extrapolated_without_rows'] > B['taper_control_min']))
    return dict(rows=out, negative_control='the limit of the tapered response against the interior ALONE must miss by the '
                                           'size of the boundary rows',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def convergence_table(variants, key, control_rule):
    """One state's refinement table: every declared variant against the base, the sum, the between-node
    residual, and the control rule put through the same kernel refinements."""
    N = TH['N']
    base = variants[f'{key}:base']
    table = dict(base=base, variants={}, control={})
    ok = base['root'] is not None and base['counts_agree'] and base['located'] == 1
    total = 0.
    for v in REP['variants']:
        row = variants[f'{key}:{v}']
        moved = None if (row['root'] is None or base['root'] is None) else _rel(row['root'], base['root'])
        table['variants'][v] = dict(root=row['root'], relative_movement=moved, kernel_nodes=row['kernel_nodes'],
                                    basis_size=row['basis_size'], counts_agree=row['counts_agree'])
        ok = ok and moved is not None and row['counts_agree'] and moved < N['each_variant']
        total += moved if moved is not None else float('inf')
    table['sum_of_movements'] = total
    table['off_node_residual'] = base.get('off_node_residual')
    ok = ok and total < N['sum_of_variants'] and base.get('off_node_residual', 1.) < N['off_node_residual']
    ctrl = variants[f'{key}:{control_rule}']
    moves = {}
    for v in REP['control_refinement']:
        row = variants[f'{key}:{control_rule}:{v}']
        moves[v] = None if (row['root'] is None or ctrl['root'] is None) else _rel(row['root'], ctrl['root'])
    lost = [v for v, m in moves.items() if m is None]          # a control that loses the root altogether is caught
    worst = max((float('inf') if m is None else m) for m in moves.values())
    table['control'] = dict(rule=control_rule, root=ctrl['root'], kernel_nodes=ctrl['kernel_nodes'], movements=moves,
                            largest_movement=worst, root_lost_under=lost, off_node_residual=ctrl.get('off_node_residual'),
                            against_base=None if ctrl['root'] is None or base['root'] is None else _rel(ctrl['root'], base['root']),
                            rejected=bool(worst >= N['control_min']))
    table['passed'] = bool(ok)
    return table


def gate_N1(variants, annuli):
    N = TH['N']
    tables = {str(k): convergence_table(variants, f'annulus{k}', 'stage8_rule') for k in range(len(N['ring_limit']))}
    must = [str(k) for k in N['control_annuli']]
    return dict(annuli=tables, control_required_on=must,
                negative_control="stage 8's node rule, put through the same kernel refinements, must be caught moving "
                                 'its root by more than the bound on the annuli where it has three or four nodes',
                passed=bool(all(t['passed'] for t in tables.values()) and all(tables[k]['control']['rejected'] for k in must)))


def gate_N2(variants, annuli, n1_passed):
    N = TH['N']
    rows = []
    for k in range(len(N['ring_limit'])):
        a = annuli[f'{k}:base']
        root = variants[f'annulus{k}:base']['root']
        ring = a['ring_root_inertial']
        wrong = a['ring_at_stage5_radius']
        rows.append(dict(index=k, dL=a['dL'], dE=a['dE'], radial_width=a['radial_width'], mean_radius=a['mean_radius'],
                         alpha=a['alpha'], annulus_root=root, ring_root=ring,
                         growth_relative_difference=None if root is None else float(root[0]/ring[0] - 1),
                         frequency_difference=None if root is None else float(root[1] - ring[1]),
                         growth_relative_difference_wrong_radius=None if root is None else float(root[0]/wrong[0] - 1)))
    d = [abs(r['growth_relative_difference']) if r['growth_relative_difference'] is not None else float('inf') for r in rows]
    falling = all(x > y for x, y in zip(d[:-1], d[1:]))
    orders = [float(np.log(d[i]/d[i + 1])/np.log(rows[i]['radial_width']/rows[i + 1]['radial_width'])) for i in range(len(d) - 1)
              if np.isfinite(d[i]) and np.isfinite(d[i + 1])]
    control = abs(rows[-1]['growth_relative_difference_wrong_radius'] or 0.)
    return dict(rows=rows, falling=bool(falling), last=d[-1], observed_orders_diagnostic_only=orders,
                numerical_budget_passed=bool(n1_passed),
                negative_control=dict(what="the ring at R = 1, stage 5's radius, instead of the annulus's own",
                                      last_difference=control, rejected=bool(control > N['ring_control_min'])),
                passed=bool(n1_passed and falling and d[-1] < N['ring_limit_last'] and control > N['ring_control_min']))


def gate_N3(variants):
    N = TH['N']
    tables = {name: convergence_table(variants, name, 'coarse_rule') for name in N['populations']}
    return dict(populations=tables,
                negative_control='a deliberately coarse node rule, put through the same kernel refinement, must be caught '
                                 'on every population',
                passed=bool(all(t['passed'] and t['control']['rejected'] for t in tables.values())))


def gate_N4(strips):
    N = TH['N']
    rows, ok = {}, True
    for name in N['root_free']:
        for config in ['base'] + list(N['root_free_configs']):
            for m in N['m_values']:
                for shifted in ((False, True) if config == 'base' else (False,)):
                    r = strips[f"{name}:m{m}:{config}:{'shifted' if shifted else 'base'}"]
                    agree = all(q['counts_agree'] for q in r['rectangles'])
                    winding = sum(q['winding'] for q in r['rectangles'])
                    located = sum(q['located'] for q in r['rectangles'])
                    rows[f"{name}:m{m}:{config}:{'shifted' if shifted else 'base'}"] = dict(
                        rectangles=len(r['rectangles']), counts_agree=bool(agree), winding=int(winding), located=int(located),
                        kernel_nodes=r['kernel_nodes'], basis_size=r['basis_size'])
                    ok = ok and agree and winding == 0 and located == 0
    name = N['root_free_control']
    found = {}
    for config in ['base'] + list(N['root_free_configs']):
        r = strips[f'{name}:m2:{config}:base']
        roots = [q for rect in r['rectangles'] for q in rect['roots']]
        found[config] = dict(roots=roots, counts_agree=bool(all(q['counts_agree'] for q in r['rectangles'])))
    seen = all(len(f['roots']) == 1 and f['counts_agree'] for f in found.values())
    spread = max((_rel(f['roots'][0], found['base']['roots'][0]) for f in found.values()), default=None) if seen else None
    return dict(rows=rows, negative_control=dict(what=f'the same sweeps on {name}, whose root is the slowest stage 8 resolved, '
                                                      'must find it under every configuration', found=found,
                                                 largest_relative_spread=spread, rejected=bool(seen)),
                passed=bool(ok and seen))


def gate_S1(rows):
    S = TH['S']
    out = {name: dict(first_order_error=r['first_order_error'], control_error=r['control_error'],
                      passed=bool(r['first_order_error'] < S['first_order_agreement']),
                      control_rejected=bool(r['control_error'] > S['control_min'])) for name, r in rows.items()}
    return dict(rows=out, negative_control='the same formula with the memory transfer function held fixed in dT/ds must '
                                           'fail to reproduce the re-solved root',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def gate_E1(reduced):
    E = TH['E']
    gamma = reduced['root'][0]
    base = [r for r in reduced['rows'] if r['variant'] == 'base']
    re = np.array([r['delta_s'][0] for r in base])
    im = np.array([r['delta_s'][1] for r in base])
    ctrl = np.array([r['control_delta_s'][0] for r in base])
    mean, se = float(re.mean()), float(re.std(ddof=1)/np.sqrt(len(re)))
    first = [r for r in reduced['rows'] if r['variant'] == 'refine_base'][0]
    refine = {r['variant']: float(abs(r['delta_s'][0] - first['delta_s'][0])/gamma)
              for r in reduced['rows'] if r['variant'] in E['refinements']}
    return dict(batches=len(base), growth_rate=gamma, implied_growth_shift=mean, standard_error=se,
                implied_growth_shift_relative=mean/gamma, standard_error_relative=se/gamma,
                implied_frequency_shift=float(im.mean()), frequency_standard_error=float(im.std(ddof=1)/np.sqrt(len(im))),
                worst_node_mean=float(np.mean([r['worst_node'] for r in base])),
                regenerated_over_imposed=[float(np.mean([r['regenerated_over_imposed_mean'][0] for r in base])),
                                          float(np.mean([r['regenerated_over_imposed_mean'][1] for r in base]))],
                refinements_relative=refine,
                negative_control=dict(what=f"the prediction scaled by {E['control_scale']}, the uniform change of the response "
                                           "that would account for B13's shortfall",
                                      implied_growth_shift_relative=float(ctrl.mean()/gamma),
                                      rejected=bool(abs(ctrl.mean())/gamma > E['control_min'])),
                passed=bool(abs(mean)/gamma < E['growth_agreement'] and se/gamma < E['standard_error_max']
                            and all(v < E['refinement_agreement'] for v in refine.values())
                            and abs(ctrl.mean())/gamma > E['control_min']))


# ================================================================ readings
def crossings(points):
    """Where a polynomial through the resolved roots reaches zero growth, for each declared degree: an
    EXTRAPOLATION of the resolved branch, never a measured boundary."""
    x, y = np.array(points).T
    out = {}
    for deg in TH['F']['fit_degrees']:
        roots = np.roots(np.polyfit(x, y, deg))
        real = [float(z.real) for z in roots if abs(z.imag) < 1e-9 and z.real > x.max()]
        out[str(deg)] = min(real) if real else None
    return out


def part_F(family, pred, b1):
    F = TH['F']
    rows = []
    for dE in F['family_B_dE']:
        new = family[f'{dE:.3f}']
        old = stage8_root(pred, f'B_dE{dE:.3f}')
        rows.append(dict(dE=dE, stage8=old, recalculated=new['root'], without_boundary_rows=new['root_without_boundary_rows'],
                         difference=None if (old is None or new['root'] is None) else
                         [new['root'][0] - old[0], new['root'][1] - old[1]],
                         growth_relative=None if (old is None or new['root'] is None) else float(new['root'][0]/old[0] - 1),
                         counts_agree=new['counts_agree'], winding=new['winding'], off_node_residual=new.get('off_node_residual')))

    def resolved(which):
        return [[r['dE'], r[which][0]] for r in rows if r[which] is not None and r[which][0] < F['threshold_fit_below']]

    return dict(family_B=rows,
                zero_growth_crossing=dict(note='an extrapolation of the resolved unstable branch in each case; nothing is '
                                               'resolved below the floor, and no boundary has been measured',
                                          stage8_points=resolved('stage8'), stage8=crossings(resolved('stage8')),
                                          recalculated_points=resolved('recalculated'),
                                          recalculated=crossings(resolved('recalculated'))),
                boundary_rows={name: dict(shift=r['shift'], growth_relative=r['growth_relative'],
                                          boundary_over_interior=r['boundary_over_interior']) for name, r in b1.items()})


# ================================================================ the campaign
def first_wave(pops, pred):
    B, N, F, S, E = TH['B'], TH['N'], TH['F'], TH['S'], TH['E']
    tasks = []
    for k in range(len(N['ring_limit'])):
        for refined in (False, True):
            tasks.append((N['ring_limit'][k][2]/20.*(4. if refined else 1.), 'annulus', f"{k}:{'refined' if refined else 'base'}",
                          (k, refined)))
    for name in N['populations']:
        family, dE = population_facts(name, pops)
        tasks.append((40. if family == 'A' else 3., 'refined_population', name, (name, family, dE)))
    for name in B['populations']:
        tasks.append((3., 'B1', name, (pops[name], stage8_root(pred, name))))
    for name in B['taper_populations']:
        for which in B['taper_selections']:
            tasks.append((5., 'B2', f"{name}:{'+'.join(which)}", (pops[name], stage8_root(pred, name), which)))
    for dE in F['family_B_dE']:
        tasks.append((4., 'family', f'{dE:.3f}', (dE,)))
    for name in S['populations']:
        tasks.append((2., 'S', name, (pops[name], stage8_root(pred, name))))
    name = E['population']
    for b in range(int(E['batches'])):
        tasks.append((60., 'E', f'{name}:b{b}:base', (pops[name], stage8_root(pred, name), b, 'base')))
    for variant in ['refine_base'] + E['refinements'] + ['replay']:
        cost = dict(half_step=30., replay=1.).get(variant, 15.)
        tasks.append((cost, 'E', f'{name}:b0:{variant}', (pops[name], stage8_root(pred, name), 0, variant)))
    for name in N['root_free'] + [N['root_free_control']]:
        ms = N['m_values'] if name != N['root_free_control'] else [2]
        for config in ['base'] + list(N['root_free_configs']):
            for m in ms:
                for shifted in ((False, True) if (config == 'base' and name != N['root_free_control']) else (False,)):
                    cost = 12.*m*(2.5 if config == 'quadrature' else 1.)
                    tasks.append((cost, 'strip', f"{name}:m{m}:{config}:{'shifted' if shifted else 'base'}",
                                  (pops[name], m, config, shifted)))
    return sorted(tasks, key=lambda t: -t[0])


def second_wave(results, pops):
    N = TH['N']
    tasks = []
    names = ['base'] + list(REP['variants'])

    def controls(rule):
        return [rule] + [f'{rule}:{v}' for v in REP['control_refinement']]

    for k in range(len(N['ring_limit'])):
        base, fine = results['annulus'][f'{k}:base']['state'], results['annulus'][f'{k}:refined']['state']
        for v in names + controls('stage8_rule'):
            tasks.append((3. if v == 'orbit_quadrature' else 1., 'variant', f'annulus{k}:{v}',
                          (base, f'annulus{k}', v, N['ring_limit_rect'], fine if v == 'population_grid' else None)))
    for name in N['populations']:
        fine = results['refined_population'][name]['state']
        for v in names + controls('coarse_rule'):
            tasks.append((3. if v == 'orbit_quadrature' else 1., 'variant', f'{name}:{v}',
                          (pops[name]['state'], name, v, N['family_rect'], fine if v == 'population_grid' else None)))
    E = TH['E']
    name = E['population']
    batches = [results['E'][k] for k in sorted(results['E'])]
    tasks.append((1., 'E_reduce', name, (pops[name], batches[0]['root'], batches)))
    return sorted(tasks, key=lambda t: -t[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--resume-from', type=Path)
    args = ap.parse_args()
    t0 = time.time()
    out_dir = args.output_dir
    (out_dir/'tasks').mkdir(parents=True, exist_ok=True)
    code = {f: hashlib.sha256(_find(f).read_bytes()).hexdigest() for f in CODE}
    task_code = {f: code[f] for f in TASK_CODE}
    protocol = TASKS.protocol_path()
    task_code['protocol'] = hashlib.sha256(protocol.read_bytes()).hexdigest()   # a toy protocol never shares a cache
    log = open(out_dir/'rut9-progress.log', 'a', encoding='utf-8')

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
            say(f'{n:4d}/{len(futures)} {kind:<19}{key:<44}{wall:8.1f}s')
        if errors:
            raise SystemExit('tasks failed, nothing assembled; finished tasks are cached for --resume-from: ' + '; '.join(errors))

    series_dir = out_dir/'rut9-series'
    series_dir.mkdir(exist_ok=True)
    pred = stage8()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        say(f'stage 9: {args.workers} workers, protocol {protocol}, output {out_dir}')
        pops = TASKS.archived_populations()
        members = sorted({n for n in TH['N']['populations'] + TH['N']['root_free'] + [TH['N']['root_free_control']]
                          if n.startswith('B_dE')})
        drain(pool, [(1., 'family_population', n, (float(n[4:]),)) for n in members])
        pops.update({n: results['family_population'][n] for n in members})
        tasks = first_wave(pops, pred)
        say(f'{len(tasks)} tasks in the first wave')
        drain(pool, tasks)
        tasks = second_wave(results, pops)
        say(f'{len(tasks)} tasks in the second wave')
        drain(pool, tasks)
    V = results['variant']
    n1 = gate_N1(V, results['annulus'])
    gates = dict(B1_reproduces_stage8_without_the_boundary_rows=gate_B1(results['B1']),
                 B2_boundary_rows_are_the_limit_of_a_resolved_taper=gate_B2(results['B2']),
                 N1_every_annulus_resolved_on_its_own=n1,
                 N2_ring_limit=gate_N2(V, results['annulus'], n1['passed']),
                 N3_every_population_root_resolved=gate_N3(V),
                 N4_root_free_conclusions_survive_refinement=gate_N4(results['strip']),
                 S1_sensitivity_formula=gate_S1(results['S']),
                 E1_eigenfunction_in_the_time_domain=gate_E1(results['E_reduce'][TH['E']['population']]))
    S7._dump_gz(series_dir/'strips.json.gz', [results['strip'][k] for k in sorted(results['strip'])])
    S7._dump_gz(series_dir/'variants.json.gz', {k: V[k] for k in sorted(V)})
    S7._dump_gz(series_dir/'time-domain.json.gz', [results['E'][k] for k in sorted(results['E'])])
    files = ['strips.json.gz', 'variants.json.gz', 'time-domain.json.gz']
    strip = lambda d: {k: v for k, v in d.items() if k != 'state'}
    body = dict(gates=gates,
                readings=dict(part_F=part_F(results['family'], pred, results['B1']), part_S=results['S'],
                              part_E=results['E_reduce'][TH['E']['population']]),
                annuli={k: strip(v) for k, v in results['annulus'].items()},
                refined_populations={k: strip(v) for k, v in results['refined_population'].items()},
                family_populations={k: {a: b for a, b in v.items() if a != 'state'} for k, v in results['family_population'].items()})
    failed = sorted(k for k, g in gates.items() if not g['passed'])
    result = dict(experiment='RUT-1 stage 9: the corrected response calculation', protocol=PROTOCOL,
                  statuses=dict(reproduction='decided by rut9_checks.py against this archive',
                                numerical_verification=dict(passed=not failed, failed_gates=failed,
                                                            note="this stage's own status; stage 8's L3 and L3b stay failed in "
                                                                 "stage 8's archive, which this stage does not touch"),
                                scientific_outcome='archived as it fell; never part of the exit status'),
                  deviations=[], thresholds=TH, **body, code_sha256_at_launch=code, tasks_reused=len(reused),
                  input_sha256={PROTOCOL: task_code['protocol'], STAGE8: hashlib.sha256(_find(STAGE8).read_bytes()).hexdigest()},
                  protocol_read_from=str(protocol),
                  series_sha256={f: hashlib.sha256((series_dir/f).read_bytes()).hexdigest() for f in sorted(files)})
    (out_dir/'rut9-results.json').write_text(json.dumps(_safe(result), indent=1) + '\n', encoding='utf-8', newline='\n')
    (out_dir/'rut9-log.json').write_text(json.dumps(dict(task_wall_seconds=walls, wall_seconds=round(time.time() - t0, 1)),
                                                     indent=1), encoding='utf-8')
    say(f'numerical verification: {"PASSED" if not failed else "FAILED " + ", ".join(failed)}')
    for k, g in sorted(gates.items()):
        say(f"  {'PASS' if g['passed'] else 'FAIL'} {k}")
    return 0 if not failed else 1


if __name__ == '__main__':
    raise SystemExit(main())
