"""RUT-1 stage 10 (protocol-rut10.md): the seeded, full-state warm-mode experiment.

    part Q   the quiet sample: it is the population, it carries no m = 1, 2, 3 source, and the live simulator
             keeps it so
    part P   the preparation: the field is seeded as declared, and the two ways of preparing it wrongly are seen
    part G   the growth of the seeded mode: linear in the seed, insensitive to the step and the grid, read over a
             window declared in advance -- and compared with the forecast of record

Three statuses are kept separate, as the owner ruled: reproduction and numerical verification set the exit
status; the scientific outcome never does. Every threshold is read from the protocol. rut10_checks.py is the
suite job.
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
import rut10_tasks as TASKS    # noqa: E402
from rut10_tasks import CODE, PROTOCOL, T0, TH, _find, _safe, np, run_task   # noqa: E402
import rut7 as S7              # noqa: E402  (series i/o, by import)
import rut8 as S8              # noqa: E402  (stage 8's growth-phase rule and its series, by import)
import rut9_tasks as R9        # noqa: E402

TASK_CODE = tuple(f for f in CODE if f not in ('rut10.py', 'rut9.py', 'rut8.py', 'rut7.py'))
STAGE8, STAGE9 = 'rut8-predictions.json', 'rut9-results.json'


# ================================================================ reading a run
def mode_series(run, m=2):
    t = np.array([r['t_periods'] for r in run['rows']])
    C = np.array([[complex(*z) for z in r['coefficients']['C']] for r in run['rows']])
    return t, np.abs(C[:, m])/np.abs(C[:, 0]), np.unwrap(np.angle(C[:, m])), C


def growth_fit(run, lo=None, hi=None):
    """The declared window: from the first record above `lo_factor` times the initial amplitude to the last
    record before the amplitude first exceeds the ceiling. Rate from log amplitude, pattern speed from phase."""
    G = TH['G']
    t, a, ph, _ = mode_series(run)
    lo = G['window_lo_factor']*a[0] if lo is None else lo
    hi = G['window_ceiling'] if hi is None else hi
    above = np.nonzero(a > lo)[0]
    if not len(above):
        return dict(window=None, note='the amplitude never rose above the lower bound')
    i0 = int(above[0])
    over = np.nonzero(a[i0:] > hi)[0]
    i1 = i0 + int(over[0]) if len(over) else len(a)
    if t[i1 - 1] - t[i0] < G['window_shortest']:
        return dict(window=[float(t[i0]), float(t[i1 - 1])], note='the window is shorter than the declared minimum')
    sel = slice(i0, i1)
    g = np.polyfit(t[sel]*T0, np.log(a[sel]), 1)
    w = np.polyfit(t[sel]*T0, ph[sel], 1)
    resid = np.log(a[sel]) - np.polyval(g, t[sel]*T0)
    return dict(window=[float(t[i0]), float(t[i1 - 1])], records=int(i1 - i0), e_foldings=float(np.log(a[i1 - 1]/a[i0])),
                growth_rate=float(g[0]), pattern_speed=float(-w[0]/2.), log_amplitude_rms_residual=float(np.sqrt(np.mean(resid**2))),
                initial_amplitude=float(a[0]))


def _key(name, quarter, k, target, variant):
    return f'{name}:q{quarter}:k{k}:t{target:+.0e}:{variant}'


# ================================================================ gates
def gate_Q1(rows):
    Q = TH['Q']
    out = {}
    for key, r in rows.items():
        quiet, ordinary = r['quiet'], r['ordinary']
        ok = (max(abs(z) for z in quiet['z_scores'].values()) < Q['moment_sigma']
              and max(quiet['source_harmonics'][:3]) < Q['symmetric_harmonic_max'])
        out[key] = dict(z_scores=quiet['z_scores'], source_harmonics=quiet['source_harmonics'][:4], passed=bool(ok),
                        control_m2_source=ordinary['source_harmonics'][1],
                        control_rejected=bool(ordinary['source_harmonics'][1] > Q['control_harmonic_min']))
    return dict(rows=out, negative_control='an ordinary draw of the same size must be seen to carry an m = 2 source',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def gate_Q2(runs):
    Q = TH['Q']
    out = {}
    for name in TH['G']['populations']:
        quiet = runs[_key(name, TH['G']['quarters'][0], 0, 0., 'base')]
        ordinary = runs[_key(name, TH['G']['quarters'][0], 0, 0., 'ordinary_draw')]
        floor = {f'm{m}': float(max(mode_series(quiet, m)[1])) for m in (1, 2, 3)}
        reached = float(max(mode_series(ordinary, 2)[1]))
        out[name] = dict(largest_amplitude_over_the_run=floor, status=quiet['status'],
                         passed=bool(quiet['status'] == 'completed' and max(floor.values()) < Q['live_floor_max']),
                         control_largest_m2=reached, control_rejected=bool(reached > Q['live_control_min']))
    return dict(rows=out, negative_control='an unseeded ORDINARY draw must be seen to grow an m = 2 pattern from its own noise',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def _ratio_at_start(run):
    co = run['rows'][0]['coefficients']
    C2, E2, C0 = complex(*co['C'][2]), complex(*co['E'][2]), complex(*co['C'][0])
    return abs(C2)/abs(C0), E2/C2


def gate_P1(runs, pops):
    P, G = TH['P'], TH['G']
    out = {}
    for name in G['populations']:
        q = G['quarters'][0]
        base = runs[_key(name, q, 0, G['targets'][0], 'base')]
        wrong = runs[_key(name, q, 0, G['targets'][0], 'static_excitation')]
        s = complex(*base['root'])
        expect = s + 1./pops[name]['state']['tau_keep']
        amp, ratio = _ratio_at_start(base)
        _, ratio_wrong = _ratio_at_start(wrong)
        out[name] = dict(target=G['targets'][0], amplitude_at_start=float(amp), excitation_over_field=[ratio.real, ratio.imag],
                         expected=[expect.real, expect.imag],
                         passed=bool(abs(amp/G['targets'][0] - 1) < P['amplitude_relative'] and abs(ratio/expect - 1) < P['excitation_relative']),
                         control_excitation_over_field=[ratio_wrong.real, ratio_wrong.imag],
                         control_rejected=bool(abs(ratio_wrong/expect - 1) > P['excitation_control_min']))
    return dict(rows=out, negative_control='the equilibrium relation delta E = delta C/tau_keep in place of (s + 1/tau_keep) delta C '
                                           'must be seen to violate the mode relation',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def gate_G1(runs):
    """Linear in the seed, read on the SAME records of two runs, so that the wandering of one realization's
    local rate cancels: while the larger run is under the ceiling, the ratio of the two amplitudes keeps its
    initial value, and the two signs keep the same amplitude and phases pi apart. Above the ceiling the ratio
    must be seen to fall: that is the control, and it is why there is a ceiling."""
    G = TH['G']
    out = {}
    for name in G['populations']:
        q = G['quarters'][0]
        t, small, ph_s, _ = mode_series(runs[_key(name, q, 0, G['targets'][0], 'base')])
        _, big, _, _ = mode_series(runs[_key(name, q, 0, G['targets'][1], 'base')])
        _, minus, ph_m, _ = mode_series(runs[_key(name, q, 0, -G['targets'][0], 'base')])
        k = min(len(small), len(big), len(minus))
        ratio = big[:k]/small[:k]/(big[0]/small[0])
        under = np.nonzero(big[:k] > G['window_ceiling'])[0]
        i1 = int(under[0]) if len(under) else k
        over = np.nonzero(big[:k] > G['control_window'][1])[0]
        i2 = int(over[0]) if len(over) else k
        in_small = small[:k] < G['window_ceiling']
        amp = float(np.max(np.abs(ratio[:i1] - 1)))
        sign_amp = float(np.max(np.abs(minus[:k][in_small]/small[:k][in_small] - 1)))
        dphi = (ph_m[:k] - ph_s[:k] + np.pi) % (2*np.pi) - np.pi                 # distance of the phase difference from pi
        sign_phase = float(np.max(np.abs(np.abs(dphi[in_small]) - np.pi)))
        drop = float(1 - np.min(ratio[i1:i2])) if i2 > i1 else None
        out[name] = dict(records_under_the_ceiling=i1, ratio_of_amplitudes_departs_by=amp, two_signs_amplitude=sign_amp,
                         two_signs_phase_from_pi=sign_phase,
                         passed=bool(amp < G['amplitude_agreement'] and sign_amp < G['sign_agreement']
                                     and sign_phase < G['sign_phase_agreement']),
                         control_ratio_falls_by=drop, control_rejected=bool(drop is not None and drop > G['ceiling_control_min']))
    return dict(rows=out, negative_control='the same ratio followed ABOVE the declared ceiling must be seen to fall',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


def gate_G2(fits):
    G = TH['G']
    out = {}
    for name in G['populations']:
        q = G['quarters'][0]
        base = fits[_key(name, q, 0, G['targets'][0], 'base')]

        def change(variant):
            f = fits[_key(name, q, 0, G['targets'][0], variant)]
            return abs(f['growth_rate']/base['growth_rate'] - 1) if 'growth_rate' in f and 'growth_rate' in base else None

        row = {variant: change(variant) for variant in G['refinements']}
        control = change(G['refinement_control'])
        out[name] = dict(relative_change=row, passed=bool(all(v is not None and v < G['refinement_agreement'] for v in row.values())),
                         control=G['refinement_control'], control_relative_change=control,
                         control_rejected=bool(control is None or control > G['refinement_control_min']))
    return dict(rows=out, negative_control='a field grid coarse enough to matter must be seen to change the rate',
                passed=bool(all(r['passed'] and r['control_rejected'] for r in out.values())))


# ================================================================ the declared reading
def cells_of(name):
    """(quarter, realizations) of the seeded base runs of one population."""
    G = TH['G']
    if name in G['populations']:
        return list(zip(G['quarters'], G['realizations']))
    return [(G['reading_only'][name]['quarter'], G['reading_only'][name]['realizations'])]


def all_populations():
    return list(TH['G']['populations']) + list(TH['G']['reading_only'])


def reading(fits, pred8, pred9):
    G, D = TH['G'], TH['D']
    out = {}
    for name in all_populations():
        cells = {}
        for q, reps in cells_of(name):
            rows = [fits[_key(name, q, k, G['targets'][0], 'base')] for k in range(reps)]
            rows = [r for r in rows if 'growth_rate' in r]
            if len(rows) < 2:
                cells[str(4*q)] = dict(realizations=len(rows), note='fewer than two fitted realizations')
                continue
            g = np.array([r['growth_rate'] for r in rows])
            w = np.array([r['pattern_speed'] for r in rows])
            cells[str(4*q)] = dict(realizations=len(rows), growth_rate=dict(mean=float(g.mean()), se=float(g.std(ddof=1)/np.sqrt(len(g))),
                                                                         values=g.tolist()),
                                   pattern_speed=dict(mean=float(w.mean()), se=float(w.std(ddof=1)/np.sqrt(len(w))), values=w.tolist()))
        top = cells[str(4*cells_of(name)[-1][0])]
        p8 = pred8['part_P']['new_populations'][name]
        verdict = 'unresolved'
        if 'growth_rate' in top and top['growth_rate']['se'] > D['power_relative']*p8['growth_rate']:
            verdict = 'unresolved: underpowered'          # a wide error bar is not an agreement
        elif 'growth_rate' in top:
            dg = abs(top['growth_rate']['mean'] - p8['growth_rate'])
            dw = abs(top['pattern_speed']['mean'] - p8['pattern_speed'])
            if (dg <= max(D['confirm_sigma']*top['growth_rate']['se'], D['confirm_relative']*p8['growth_rate'])
                    and dw <= max(D['confirm_sigma']*top['pattern_speed']['se'], D['confirm_pattern_relative']*abs(p8['pattern_speed']))):
                verdict = 'confirmed'
            elif dg > D['confirm_sigma']*top['growth_rate']['se'] and dg > D['contradict_relative']*p8['growth_rate']:
                verdict = 'contradicted'
        out[name] = dict(forecast_of_record=dict(growth_rate=p8['growth_rate'], pattern_speed=p8['pattern_speed']),
                         recalculated=pred9.get(name), measured=cells, reading=verdict)
    return out


# ================================================================ two diagnostics of stage 8's unseeded reading
def diagnostics(runs, background=None):
    """Readings, with no verdict attached. (1) Stage 8's own growth-phase rule applied to the seeded quiet runs:
    is the rule biased when there is no shot noise? (2) Stage 8's archived unseeded runs read again with the
    m = 2 power of B24 -- the root-free population at the same body count, whose m = 2 pattern is shot noise and
    nothing else -- subtracted as a noise floor. B24's noise is a PROXY for the other populations' and the
    reading is rough; it is declared so that it is made once, as written."""
    G = TH['G']
    out = dict(stage8_rule_on_seeded_runs={}, stage8_unseeded_runs_with_a_noise_floor={}, unseeded_quiet_population={})
    for name in G['populations']:                 # does the population heat when m = 2 cannot grow? a reading, no gate
        rows = runs[_key(name, G['quarters'][0], 0, 0., 'base')]['rows']
        t = np.array([r['t_periods'] for r in rows])
        early, late = t <= G['quiet_windows'][0], t >= t[-1] - G['quiet_windows'][1]
        out['unseeded_quiet_population'][name] = {
            q: dict(early=float(np.mean([r[q] for r, e in zip(rows, early) if e])),
                    late=float(np.mean([r[q] for r, l in zip(rows, late) if l])))
            for q in ('radial_spread', 'residual_radial_dispersion', 'residual_azimuthal_dispersion', 'L_mean', 'support_on_bodies')}
    for name in G['populations']:
        q = G['quarters'][0]
        cell = [dict(runs[_key(name, q, k, G['targets'][0], 'base')], live=True) for k in range(G['realizations'][0])]
        if cell[0]['rows'][-1]['t_periods'] <= S8.TH['M']['growth_start']:
            out['stage8_rule_on_seeded_runs'][name] = dict(note="the run ends before stage 8's rule starts")
            continue
        out['stage8_rule_on_seeded_runs'][name] = S8.growth_phase(cell)
    # does the reading depend on where the window starts? the same fit from later and later in the run
    out['window_start'] = {}
    for name in all_populations():
        for q, reps in cells_of(name):
            for factor in G['window_lo_factors_reported']:
                rows = []
                for k in range(reps):
                    run = runs[_key(name, q, k, G['targets'][0], 'base')]
                    f = growth_fit(run, lo=factor*mode_series(run)[1][0])
                    if 'growth_rate' in f:
                        rows.append(f)
                if len(rows) >= 2:
                    g, w = np.array([r['growth_rate'] for r in rows]), np.array([r['pattern_speed'] for r in rows])
                    out['window_start'][f'{name}:N{4*q}:from_{factor:g}_times_the_seed'] = dict(
                        realizations=len(rows), growth_rate=[float(g.mean()), float(g.std(ddof=1)/np.sqrt(len(g)))],
                        pattern_speed=[float(w.mean()), float(w.std(ddof=1)/np.sqrt(len(w)))],
                        window_starts_periods=[r['window'][0] for r in rows])
    # the live axisymmetric field on the ring against its primed value, and what the predicted root's sensitivity
    # to a uniformly scaled C0 would make of it: a rough reading of how much of a run's offset is its background
    out['live_background'] = {}
    for key, run in runs.items():
        if run['target'] == 0. or run['variant'] != 'base' or 'growth_rate' not in growth_fit(run):
            continue
        f = growth_fit(run)
        t_, _, _, C = mode_series(run)
        sel = (t_ >= f['window'][0]) & (t_ <= f['window'][1])
        drift = float(np.mean(np.abs(C[sel, 0]))/abs(C[0, 0]) - 1.)
        row = dict(mean_C0_over_primed_minus_one=drift)
        sens = (background or {}).get(run['name'])
        if sens is not None:
            d_re, d_im = sens['d_root_per_unit_relative_change_of_C0']
            row.update(implied_growth_rate_shift=d_re*drift, implied_pattern_speed_shift=-.5*d_im*drift)
        out['live_background'][key] = row
    series = S8.load_series(_find(STAGE8).parent/'rut8-series')
    M8 = S8.TH['M']

    def power(name, n):
        live = [r for r in series['X'][(name, n)] if r['live']]
        t = np.array([row['t_periods'] for row in live[0]['rows']])
        a = np.array([[abs(complex(*row['coefficients']['C'][2]))/abs(complex(*row['coefficients']['C'][0]))
                       for row in r['rows']] for r in live])
        return t, (a*a).mean(0)

    for name in all_populations():
        for n in M8['bodies']:
            declared = S8.growth_phase(series['X'][(name, n)])
            if not declared.get('growth_phase'):
                continue
            t, p = power(name, n)
            tn, pn = power(G['noise_proxy'], n)
            k = min(len(t), len(tn))
            mode = p[:k] - pn[:k]
            w0, w1 = declared['window_periods']
            sel = (t[:k] >= w0 - 1e-9) & (t[:k] <= w1 + 1e-9) & (mode > 0)
            row = dict(window_periods=[w0, w1], records=int(sel.sum()), declared_rate=declared['growth_rate'])
            if sel.sum() >= 4:
                row['rate_of_total_power'] = float(.5*np.polyfit(t[:k][sel]*T0, np.log(p[:k][sel]), 1)[0])
                row['rate_with_noise_floor_subtracted'] = float(.5*np.polyfit(t[:k][sel]*T0, np.log(mode[sel]), 1)[0])
                row['noise_over_total_at_window_start'] = float(pn[:k][sel][0]/p[:k][sel][0])
            out['stage8_unseeded_runs_with_a_noise_floor'][f'{name}:N{n}'] = row
    return out


# ================================================================ the campaign
def declared_runs():
    G = TH['G']
    runs = []
    for name in G['populations']:
        h = G['horizon'][name]
        q0 = G['quarters'][0]
        for q, reps in zip(G['quarters'], G['realizations']):
            for k in range(reps):
                runs.append((name, q, k, G['targets'][0], 'base', h))
        runs += [(name, q0, 0, -G['targets'][0], 'base', h), (name, q0, 0, G['targets'][1], 'base', h)]
        runs += [(name, q0, 0, G['targets'][0], v, h)
                 for v in list(G['refinements']) + [G['refinement_control'], 'unprepared_bodies', 'static_excitation']]
        runs += [(name, q0, 0, 0., 'base', h), (name, q0, 0, 0., 'ordinary_draw', G['unseeded_ordinary_horizon'])]
    for name, spec in G['reading_only'].items():
        runs += [(name, spec['quarter'], k, G['targets'][0], 'base', spec['horizon']) for k in range(spec['realizations'])]
    return runs


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
    task_code['protocol'] = hashlib.sha256(protocol.read_bytes()).hexdigest()
    task_code['stage9_protocol'] = hashlib.sha256(R9.protocol_path().read_bytes()).hexdigest()
    log = open(out_dir/'rut10-progress.log', 'a', encoding='utf-8')

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
            say(f'{n:4d}/{len(futures)} {kind:<12}{key:<52}{str(out.get("status", "")):<11}{wall:8.1f}s')
        if errors:
            raise SystemExit('tasks failed, nothing assembled; finished tasks are cached for --resume-from: ' + '; '.join(errors))

    series_dir = out_dir/'rut10-series'
    series_dir.mkdir(exist_ok=True)
    pops = R9.archived_populations()
    G = TH['G']
    tasks = []
    for name, q, k, target, variant, h in declared_runs():
        cost = q*h*(2.6 if variant == 'fine_grid' else 2. if variant == 'half_step' else .3 if variant == 'coarse_grid' else 1.)
        tasks.append((cost, 'run', _key(name, q, k, target, variant), (pops[name], q, k, target, variant, h)))
    for name in all_populations():
        for q, reps in cells_of(name):
            tasks += [(1., 'Q', f'{name}:q{q}:k{k}', (pops[name], q, k)) for k in range(reps)]
    for name in G['populations']:
        tasks.append((3., 'background_sensitivity', name, (pops[name],)))
        tasks += [(1., 'field_state', f'{name}:t{t:+.0e}', (pops[name], t))
                  for t in (G['targets'][0], -G['targets'][0], G['targets'][1])]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        say(f'stage 10: {args.workers} workers, protocol {protocol}, output {out_dir}; {len(tasks)} tasks')
        drain(pool, sorted(tasks, key=lambda t: -t[0]))
    runs = results['run']
    fits = {key: growth_fit(r) for key, r in runs.items() if r['target'] != 0.}
    gates = dict(Q1_the_quiet_sample=gate_Q1(results['Q']), Q2_the_symmetry_survives_the_live_run=gate_Q2(runs),
                 P1_the_seeded_field=gate_P1(runs, pops), G1_linear_in_the_seed=gate_G1(runs),
                 G2_step_and_grid=gate_G2(fits),
                 G3_every_run_completed=dict(statuses={k: r['status'] for k, r in runs.items()},
                                             passed=bool(all(r['status'] == 'completed' for r in runs.values()))))
    pred8 = json.loads(_find(STAGE8).read_text(encoding='utf-8'))
    stage9 = json.loads(_find(STAGE9).read_text(encoding='utf-8'))
    pred9 = {}
    for name in all_populations():
        dE = pops[name]['dE']
        row = [r for r in stage9['readings']['part_F']['family_B'] if abs(r['dE'] - dE) < 1e-9][0]
        pred9[name] = dict(growth_rate=row['recalculated'][0], pattern_speed=-row['recalculated'][1]/2.)
    # series: every run's records, its complete initial bodies, and the complete initial fields
    files = []
    for key in sorted(runs):
        f = 'run_' + key.replace(':', '_').replace('+', 'p').replace('-', 'n') + '.json.gz'
        S7._dump_gz(series_dir/f, runs[key])
        files.append(f)
    fields = results['field_state']
    TASKS.dump_arrays(series_dir/'initial-fields.f8.gz', {f'{k}|{part}': np.array(fields[k][part]) for k in fields for part in ('C', 'E', 'axis')})
    S7._dump_gz(series_dir/'initial-fields.json.gz', {k: {a: b for a, b in fields[k].items() if a not in ('C', 'E', 'axis')}
                                                       for k in sorted(fields)})
    S7._dump_gz(series_dir/'quiet-samples.json.gz', {k: results['Q'][k] for k in sorted(results['Q'])})
    files += ['initial-fields.f8.gz', 'initial-fields.json.gz', 'quiet-samples.json.gz']
    failed = sorted(k for k, g in gates.items() if not g['passed'])
    result = dict(experiment='RUT-1 stage 10: the seeded, full-state warm-mode experiment', protocol=PROTOCOL,
                  statuses=dict(reproduction='decided by rut10_checks.py against this archive',
                                numerical_verification=dict(passed=not failed, failed_gates=failed),
                                scientific_outcome='archived as it fell; never part of the exit status'),
                  deviations=[], thresholds=TH, gates=gates, fits=fits, reading=reading(fits, pred8, pred9),
                  diagnostics=diagnostics(runs, results['background_sensitivity']),
                  background_sensitivity=results['background_sensitivity'],
                  run_index={k: dict(status=r['status'], bodies=r['bodies'], seconds=r['seconds'], steps=r['steps'],
                                     sample_sha256=r['sample_sha256'], initial_sha256=r['initial_sha256'],
                                     field_sha256=r['field_sha256']) for k, r in runs.items()},
                  code_sha256_at_launch=code, tasks_reused=len(reused),
                  input_sha256={PROTOCOL: task_code['protocol'], 'protocol-rut9.md': task_code['stage9_protocol'],
                                STAGE8: hashlib.sha256(_find(STAGE8).read_bytes()).hexdigest(),
                                STAGE9: hashlib.sha256(_find(STAGE9).read_bytes()).hexdigest()},
                  protocol_read_from=str(protocol),
                  series_sha256={f: hashlib.sha256((series_dir/f).read_bytes()).hexdigest() for f in sorted(files)})
    (out_dir/'rut10-results.json').write_text(json.dumps(_safe(result), indent=1) + '\n', encoding='utf-8', newline='\n')
    (out_dir/'rut10-log.json').write_text(json.dumps(dict(task_wall_seconds=walls, wall_seconds=round(time.time() - t0, 1)),
                                                      indent=1), encoding='utf-8')
    say(f'numerical verification: {"PASSED" if not failed else "FAILED " + ", ".join(failed)}')
    for k, g in sorted(gates.items()):
        say(f"  {'PASS' if g['passed'] else 'FAIL'} {k}")
    for name, r in result['reading'].items():
        say(f"  reading {name}: {r['reading']}")
    return 0 if not failed else 1


if __name__ == '__main__':
    raise SystemExit(main())
