"""RUT-1 stage 7 (protocol-rut7.md): the cold-ring spectrum certified down to the axis, the population
defined and drawn faithfully, and the same bodies evolved under a live and a frozen field.

This is the campaign driver. It takes of the order of an hour on many cores and is NOT a suite job;
rut7_checks.py is, and verifies in minutes that the committed archive is what the committed code produces.

Three statuses are kept separate, as the owner ruled on 72aef44:

* reproduction             -- does a run reproduce its frozen archive? Sets the exit status.
* numerical verification   -- do the identities, samplers, integrators and solvers pass their correctness
                              checks, and does every check reject its negative control? Sets the exit status.
* scientific outcome       -- what the live-against-frozen experiment shows. Archived however it falls, and
                              never touches the exit status.

Every threshold is read from the machine-readable block of the protocol: this file holds no tolerance of
its own. Anything done differently from the protocol is written into `deviations`, named as a departure.

The sampler of record at every source-sampling call site is the owner's corrected node sampler
(research_work/annulus_sampling, annulus-phase-space-v2). The historical equilibrium.py and rut6.py are
imported where a negative control needs the stage 6 behaviour and are never edited.
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
import rut7_tasks as TASKS  # noqa: E402  (sets one thread per worker before numpy is imported)
from rut7_tasks import (CODE, FM, GRID, OWNER_PACKAGE, PO, POPULATIONS, R_CASES, REFINEMENTS, RUNGS, T0, TH,  # noqa: E402
                        V6_CONTROLS, W, CubicSpline, _annulus, _find, _safe, evidence_io, np, run_task)

TASK_CODE = tuple(f for f in CODE if f != 'rut7.py')       # what a cached task result depends on


def gate_E3_rows(spectra):
    """Completeness, evaluated on every sub-rectangle of every rung and m, in both partitions."""
    E = TH['E']
    failures = []
    worst = dict(smallest_gap=None, residual=0., polish_movement=0., partition=0., phase_step=0.,
                 winding_remainder=0.)
    for row in spectra:
        tag = f"{row['model']}:m{row['m']}"
        for part in ('base', 'shifted'):
            s = row[part]
            if not s['counts_agree']:
                failures.append(f'{tag}:{part}: Beyn count differs from the winding number')
            if s['any_saturated']:
                failures.append(f'{tag}:{part}: moment capacity saturated')
            gap = s['smallest_gap']
            if gap is not None and math.isfinite(gap):
                worst['smallest_gap'] = gap if worst['smallest_gap'] is None else min(worst['smallest_gap'], gap)
                if gap < E['singular_gap']:
                    failures.append(f'{tag}:{part}: singular gap {gap:.3g}')
            if s['worst_residual'] >= E['polished_residual']:
                failures.append(f'{tag}:{part}: residual {s["worst_residual"]:.3g}')
            if s['worst_polish_movement'] >= E['polish_movement']:
                failures.append(f'{tag}:{part}: polish moved a root {s["worst_polish_movement"]:.3g}')
            if s['worst_winding_phase_step'] > E['winding_phase_step']*(1 + 1e-12):
                failures.append(f'{tag}:{part}: winding phase step {s["worst_winding_phase_step"]:.3g}')
            worst['residual'] = max(worst['residual'], s['worst_residual'])
            worst['polish_movement'] = max(worst['polish_movement'], s['worst_polish_movement'])
            worst['phase_step'] = max(worst['phase_step'], s['worst_winding_phase_step'])
            worst['winding_remainder'] = max(worst['winding_remainder'], s['worst_winding_remainder'])
        ag = row['partition_agreement']
        worst['partition'] = max(worst['partition'], ag['worst_relative'])
        if ag['unpaired_first'] or ag['unpaired_second'] or ag['worst_relative'] > E['partition_agreement']:
            failures.append(f'{tag}: the two partitions return different root sets')
    return failures, worst


def gate_V4_rows(pops, refined):
    """How far each refinement moves each population: integrated quantities, and the two profiles."""
    S = TH['S']
    rows, worst_i, worst_p = {}, 0., 0.
    for name, pop in pops.items():
        base = pop['description']
        a = _annulus(pop)
        calibrated = pop['alpha_source'].startswith('solved')
        for variant in REFINEMENTS:
            ref = refined[(name, variant)]
            d = ref['description']
            moved = dict(mean_support=abs(d['support']['mean']/base['support']['mean'] - 1),
                         mean_radius=abs(d['mean_radius']/base['mean_radius'] - 1),
                         radial_width=abs(d['radial_width']/base['radial_width'] - 1),
                         peak_field=abs(d['peak_field']/base['peak_field'] - 1))
            if calibrated:
                moved['alpha'] = abs(d['alpha']/base['alpha'] - 1)
            r = np.asarray(ref['r'])
            inside = (a.r >= r[0]) & (a.r <= r[-1])
            sig = CubicSpline(r, ref['sigma'])(a.r[inside])
            fld = CubicSpline(r, ref['C'])(a.r[inside])
            prof = dict(density=float(np.max(np.abs(sig - a.sigma[inside]))/a.sigma.max()),
                        field=float(np.max(np.abs(fld - a.C[inside]))/a.C.max()))
            rows[f'{name}:{variant}'] = dict(integrated=moved, profiles=prof)
            worst_i = max(worst_i, *moved.values())
            worst_p = max(worst_p, *prof.values())
    return rows, float(worst_i), float(worst_p), bool(worst_i < S['convergence_integrated']
                                                      and worst_p < S['convergence_profile'])


def measure_R(series, window):
    """Growth rate and rotating-frame frequency of the recorded m = 2 coefficient: straight-line fits to
    log|c| and to the unwrapped phase, the declared estimator for a mode that arrives alone."""
    t = np.array(series['t_periods'])*T0
    c = np.array([complex(a, b) for a, b in series['c2']])
    sel = (t >= window[0]*T0) & (t <= window[1]*T0) & (np.abs(c) > 0)
    if sel.sum() < 16:
        return None
    return dict(growth_rate=float(np.polyfit(t[sel], np.log(np.abs(c[sel])), 1)[0]),
                frequency_rotating_frame=float(np.polyfit(t[sel], np.unwrap(np.angle(c[sel])), 1)[0]),
                amplitude_start=float(np.abs(c[sel][0])), amplitude_end=float(np.abs(c[sel][-1])),
                samples=int(sel.sum()))


def gate_R1(series, window=None):
    Rt = TH['R']
    window = Rt['fit_window'] if window is None else window
    rows = {}
    for name, s in series.items():
        m = measure_R(s, window) if s['status'] == FM.COMPLETED else None
        row = dict(status=s['status'], measured=m)
        if m:
            row['growth_relative'] = abs(m['growth_rate']/Rt['predicted_growth'] - 1)
            row['frequency_relative'] = abs(m['frequency_rotating_frame'] - Rt['predicted_frequency'])/Rt['predicted_frequency']
            row['within'] = bool(row['growth_relative'] < Rt['growth_relative']
                                 and row['frequency_relative'] < Rt['frequency_relative'])
        else:
            row['within'] = False
        rows[name] = row
    declared = [k for k in rows if not k.startswith('control')]
    control = rows['control_no_primed_field']
    return dict(rows=rows, predicted=dict(growth_rate=Rt['predicted_growth'], frequency=Rt['predicted_frequency']),
                fit_window_periods=list(window), negative_control_rejected=bool(not control['within']),
                passed=bool(all(rows[k]['within'] for k in declared) and not control['within']))


Q_SCALED_BY_OWN_START = ('mean_radius', 'radial_spread', 'residual_radial_dispersion',
                         'residual_azimuthal_dispersion', 'vr_rms', 'L_mean', 'L_spread', 'support_on_bodies',
                         'support_on_ring')


Q_STREAMING = ('radial_m1_streaming', 'radial_m2_streaming')


def _window_mean(rows, key, window):
    vals = [r[key] for r in rows if window[0] - 1e-9 <= r['t_periods'] <= window[1] + 1e-9]
    return float(np.mean(vals))


def pair_differences(live, frozen, early, late):
    """D_Q for one pair: the difference of the two runs' changes, over a declared scale that is a property
    of the initial population and never another run's drift."""
    first = live['rows'][0]
    sigma0, spread0 = first['residual_radial_dispersion'], first['radial_spread']
    out = {}
    for key in Q_SCALED_BY_OWN_START + Q_STREAMING + ('radial_breathing',):
        if key in Q_SCALED_BY_OWN_START:
            scale = first[key]
        elif key in Q_STREAMING:
            scale = sigma0
        else:
            scale = sigma0/spread0
        change = lambda run: _window_mean(run['rows'], key, late) - _window_mean(run['rows'], key, early)
        dl, df = change(live)/scale, change(frozen)/scale
        out[key] = dict(D=dl - df, live_change=dl, frozen_change=df)
    return out


def mode_fit(live, m, window, field='C'):
    """Growth rate of |c_m| and its pattern speed over the window, from the saved complex coefficients."""
    t = np.array([r['t_periods'] for r in live['rows']])*T0
    c = np.array([complex(*r['coefficients'][field][m]) for r in live['rows']])
    c0 = np.array([complex(*r['coefficients'][field][0]) for r in live['rows']])
    sel = (t >= window[0]*T0 - 1e-9) & (t <= window[1]*T0 + 1e-9) & (np.abs(c) > 0)
    if sel.sum() < 8:
        return None
    rate = float(np.polyfit(t[sel], np.log(np.abs(c[sel])), 1)[0])
    omega = float(np.polyfit(t[sel], np.unwrap(np.angle(c[sel])), 1)[0])
    return dict(growth_rate=rate, pattern_speed=-omega/m,
                relative_amplitude_start=float(np.abs(c[sel][0])/np.abs(c0[sel][0])),
                relative_amplitude_end=float(np.abs(c[sel][-1])/np.abs(c0[sel][-1])))


def _ensemble(values):
    v = np.asarray(values, float)
    return dict(mean=float(v.mean()), standard_error=float(v.std(ddof=1)/math.sqrt(len(v))) if len(v) > 1 else None,
                realizations=int(len(v)), values=[float(x) for x in v])


def _weighted_line(x, y, sy):
    """Slope of y against x with its uncertainty, weights 1/sy^2."""
    x, y, w = np.asarray(x, float), np.asarray(y, float), 1./np.asarray(sy, float)**2
    xm = np.sum(w*x)/np.sum(w)
    sxx = np.sum(w*(x - xm)**2)
    slope = float(np.sum(w*(x - xm)*y)/sxx)
    return slope, float(1./math.sqrt(sxx))


def analyze_X(series, windows=None):
    """Everything the protocol declares for part X, from the archived series alone."""
    X = TH['X']
    windows = windows or dict(early=X['early_window'], late=X['late_window'], mode=X['mode_fit_window'])
    table = {}
    for (name, n), runs in series.items():
        pairs = {}
        for run in runs:
            pairs.setdefault(run['realization'], {})['live' if run['live'] else 'frozen'] = run
        done = [p for p in pairs.values() if len(p) == 2 and all(r['status'] == 'completed' for r in p.values())]
        cell = dict(pairs_declared=len(pairs), pairs_completed=len(done))
        if done:
            diffs = [pair_differences(p['live'], p['frozen'], windows['early'], windows['late']) for p in done]
            cell['D'] = {k: dict(_ensemble([d[k]['D'] for d in diffs]),
                                 live_change_mean=float(np.mean([d[k]['live_change'] for d in diffs])),
                                 frozen_change_mean=float(np.mean([d[k]['frozen_change'] for d in diffs])))
                         for k in diffs[0]}
            cell['modes'] = {}
            for m in range(1, 5):
                fits = [mode_fit(p['live'], m, windows['mode']) for p in done]
                fits = [f for f in fits if f]
                if not fits:
                    continue
                cell['modes'][f'm{m}'] = dict(
                    growth_rate=_ensemble([f['growth_rate'] for f in fits]),
                    pattern_speed=_ensemble([f['pattern_speed'] for f in fits]),
                    relative_amplitude_start=_ensemble([f['relative_amplitude_start'] for f in fits]),
                    relative_amplitude_end=_ensemble([f['relative_amplitude_end'] for f in fits]),
                    excitation_growth_rate=_ensemble([g['growth_rate'] for g in
                                                      (mode_fit(p['live'], m, windows['mode'], 'E') for p in done) if g]))
        table.setdefault(name, {})[str(n)] = cell
    return table


def readings_X(table, populations, bodies=None):
    """The declared readings, applied as declared. `unresolved` is a reading, not a failure."""
    X = TH['X']
    counts = [str(n) for n in (X['bodies'] if bodies is None else bodies)]
    top = counts[-1]
    out = {}
    for name, cells in table.items():
        row = dict(collective_mode={}, quiet={}, discreteness=None)
        for m in range(1, 5):
            key = f'm{m}'
            g = [cells.get(n, {}).get('modes', {}).get(key, {}).get('growth_rate') for n in counts]
            if any(x is None or x['standard_error'] is None for x in g):
                row['collective_mode'][key] = dict(reading='unresolved', why='a body count is missing')
                continue
            above = [x['mean'] - X['mode_sigma']*x['standard_error'] > X['mode_rate_floor'] for x in g]
            agree = all(abs(g[i]['mean'] - g[j]['mean'])
                        <= X['mode_sigma']*math.hypot(g[i]['standard_error'], g[j]['standard_error'])
                        for i in range(len(g)) for j in range(i + 1, len(g)))
            if all(above) and agree:
                reading = 'present'
            elif g[-1]['mean'] + X['mode_sigma']*g[-1]['standard_error'] < X['mode_rate_floor']:
                reading = 'absent'
            else:
                reading = 'unresolved'
            amp = [cells[n]['modes'][key]['relative_amplitude_start']['mean'] for n in counts]
            row['collective_mode'][key] = dict(
                reading=reading, rates={n: [x['mean'], x['standard_error']] for n, x in zip(counts, g)},
                above_floor_at_each_count=above, rates_agree=bool(agree),
                amplitude_exponent_in_N=float(np.polyfit(np.log([float(n) for n in counts]), np.log(amp), 1)[0]),
                pattern_speed_at_largest_count=cells[top]['modes'][key]['pattern_speed']['mean'])
        D = cells.get(top, {}).get('D')
        for q, bound in X['quiet_bounds'].items():
            if not D or D[q]['standard_error'] is None:
                row['quiet'][q] = dict(reading='unresolved', why='the largest body count did not complete')
                continue
            mean, se = abs(D[q]['mean']), D[q]['standard_error']
            reading = ('quiet' if mean + X['quiet_sigma']*se <= bound
                       else 'exceeds the bound' if mean - X['quiet_sigma']*se > bound else 'unresolved')
            row['quiet'][q] = dict(reading=reading, D=D[q]['mean'], standard_error=se, bound=bound)
        q = 'residual_radial_dispersion'
        pts = [cells.get(n, {}).get('D', {}).get(q) for n in counts]
        if all(p and p['standard_error'] for p in pts):
            means, ses = [p['mean'] for p in pts], [p['standard_error'] for p in pts]
            resolved = all(abs(mu) > X['quiet_sigma']*se for mu, se in zip(means, ses))
            same_sign = len({mu > 0 for mu in means}) == 1
            if resolved and same_sign:
                p, sp = _weighted_line(np.log([float(n) for n in counts]), np.log(np.abs(means)),
                                       [se/abs(mu) for mu, se in zip(means, ses)])
                ref = X['discreteness_exponent']
                reading = ('discreteness-dominated' if p + 2*sp < ref else 'collective' if p - 2*sp > ref
                           else 'unresolved')
                row['discreteness'] = dict(reading=reading, exponent=p, exponent_uncertainty=sp, D=means,
                                           standard_errors=ses)
            else:
                row['discreteness'] = dict(reading='unresolved', D=means, standard_errors=ses,
                                           why='the change is not resolved from zero with one sign at every '
                                               'body count, so a power law cannot be fitted to it')
        else:
            row['discreteness'] = dict(reading='unresolved', why='a body count is missing')
        out[name] = row
    family_B = sorted((p['dE'], p['name']) for p in populations.values() if p['family'] == 'B')
    trend = {}
    if all(top in table.get(name, {}) and table[name][top].get('D') for _, name in family_B):
        dEs = [d for d, _ in family_B]
        rate = [table[name][top]['modes'].get('m2', {}).get('growth_rate') for _, name in family_B]
        if all(r and r['standard_error'] for r in rate):
            s, ss = _weighted_line(dEs, [r['mean'] for r in rate], [r['standard_error'] for r in rate])
            trend['m2_growth_rate'] = dict(slope_per_unit_dE=s, uncertainty=ss,
                                           values={name: [r['mean'], r['standard_error']]
                                                   for (_, name), r in zip(family_B, rate)})
        for q in X['quiet_bounds']:
            pts = [table[name][top]['D'][q] for _, name in family_B]
            if all(p['standard_error'] for p in pts):
                s, ss = _weighted_line(dEs, [p['mean'] for p in pts], [p['standard_error'] for p in pts])
                trend[q] = dict(slope_per_unit_dE=s, uncertainty=ss,
                                values={name: [p['mean'], p['standard_error']]
                                        for (_, name), p in zip(family_B, pts)})
    return out, dict(family='B', body_count=int(top), trends=trend,
                     note='family B only: the mass, the coupling, the footprint and the response times are '
                          'fixed, and dE is the only thing that changes')


# ================================================================ gates assembled from task results
def gate_V6(runs):
    """Stationarity of the draw, per realization as declared and on the ensemble as amended; and the three
    broken draws, each of which must be caught."""
    S = TH['S']
    by = {}
    for r in runs:
        by.setdefault((r['name'], r['control']), []).append(r)
    rows, ok, rejected = {}, True, True
    for (name, control), group in sorted(by.items(), key=lambda kv: (kv[0][0], kv[0][1] or '')):
        group = sorted(group, key=lambda r: r['realization'])
        done = [r for r in group if r['status'] == 'completed' and r['statistics']]
        z = {k: [r['statistics']['z_initial_against_time_average'][k] for r in done]
             for k in ('mean_radius', 'radial_spread', 'vr_rms', 'vt_mean')} if done else {}
        ensemble = {k: float(np.mean(v)*math.sqrt(len(v))) for k, v in z.items()}
        each = [bool(r['statistics']['stationary']) for r in done]
        together = bool(done) and all(abs(v) <= S['stationary_sigma'] for v in ensemble.values())
        stationary = bool(len(done) == len(group) and all(each) and together)
        rows[f"{name}:{control or 'declared_draw'}"] = dict(
            realizations=len(group), completed=len(done), each_realization_stationary=each,
            ensemble_z=ensemble, ensemble_stationary=together,
            worst_single_z=float(max([abs(x) for v in z.values() for x in v] or [0.])),
            worst_fluctuation=float(max([r['statistics']['mean_radius_fluctuation_in_standard_errors']
                                         for r in done] or [0.])),
            stationary=stationary)
        if control is None:
            ok = ok and stationary
        else:
            rejected = rejected and not stationary
    return dict(rows=rows, negative_controls_rejected=bool(rejected), passed=bool(ok and rejected))


def gate_V8(v6_runs, r_series, x_runs):
    """Every declared run completes with its status recorded; each pair starts from identical bodies.
    A negative control may end early -- a kicked population that leaves the box has been caught -- so its
    status is recorded and is not required to be `completed`."""
    bad = [f"V6:{r['name']}:{r['realization']}:{r['status']}" for r in v6_runs
           if r['control'] is None and r['status'] != 'completed']
    bad += [f"R:{k}:{s['status']}" for k, s in r_series.items()
            if not k.startswith('control') and s['status'] != FM.COMPLETED]
    controls = {f"V6:{r['name']}:{r['control']}:k{r['realization']}": r['status'] for r in v6_runs if r['control']}
    controls.update({f'R:{k}': s['status'] for k, s in r_series.items() if k.startswith('control')})
    bad += [f"X:{r['name']}:N{r['bodies']}:k{r['realization']}:{'live' if r['live'] else 'frozen'}:{r['status']}"
            for r in x_runs if r['status'] != 'completed']
    pairs = {}
    for r in x_runs:
        pairs.setdefault((r['name'], r['bodies'], r['realization']), {})['live' if r['live'] else 'frozen'] = r
    unpaired = [list(k) for k, p in pairs.items() if len(p) != 2]
    differ = [list(k) for k, p in pairs.items() if len(p) == 2
              and p['live']['initial_sha256'] != p['frozen']['initial_sha256']]
    return dict(runs=len(v6_runs) + len(r_series) + len(x_runs), not_completed=bad,
                negative_control_statuses=controls, pairs=len(pairs),
                unpaired=unpaired, pairs_with_different_bodies=differ,
                passed=bool(not bad and not unpaired and not differ))


def declared_tasks(pops, smoke):
    """(priority, kind, key, args): the whole campaign, heaviest first."""
    S, X, E = TH['S'], TH['X'], TH['E']
    specs = {p['name']: p for p in POPULATIONS}
    alpha_B = pops['B_mid']['description']['alpha']
    tasks = []
    m_max = 2 if smoke else E['m_max']
    for model in RUNGS:
        for m in range(m_max + 1):
            tasks.append((3 if model == 'two_stage' else 1, 'spectrum', f'{model}:m{m}', (model, m)))
    for kind in ('E1', 'E2', 'E3_control', 'E4'):
        tasks.append((4, kind, kind, ()))
    for name, pop in pops.items():
        alpha = None if specs[name]['alpha_from'] is None else alpha_B
        for kind in ('V1', 'V2', 'V3', 'V5'):
            tasks.append((5, kind, name, (pop,)))
        for variant in REFINEMENTS:
            tasks.append((6, 'V4', f'{name}:{variant}', (specs[name], variant, alpha)))
        reps = 1 if smoke else int(S['stationary_realizations'])
        for control in (None,) + V6_CONTROLS:
            for k in range(reps):
                args = (pop, k, control) + ((8., 1024) if smoke else ())
                tasks.append((20, 'V6', f"{name}:{control or 'declared_draw'}:k{k}", args))
        bodies = (16, 32, 64) if smoke else X['bodies']
        counts = (2, 2, 2) if smoke else X['realizations']
        for n, reps in zip(bodies, counts):
            for k in range(reps):
                for live in (True, False):
                    args = (pop, n, k, live) + ((4.,) if smoke else ())
                    tasks.append(((.1 if live else .04)*n, 'X', f"{name}:N{n}:k{k}:{'live' if live else 'frozen'}", args))
    for dE in (.010, .030):
        for n_r in (GRID['n_r'], REFINEMENTS['radial_grid']['n_r']):
            tasks.append((6, 'V4_control', f'dE{dE}:n{n_r}', (dE, n_r)))
    tasks.append((8, 'V7', 'V7', (pops['B_mid'],)))
    for case in R_CASES:
        tasks.append((8, 'R', case['name'], (case,) + ((6.,) if smoke else ())))
    return sorted(tasks, key=lambda t: -t[0])


def _dump_gz(path, obj, digits=None):
    # sorted keys: tasks finish in a different order every time, and a committed file's hash must not
    # depend on which worker was quickest
    text = json.dumps(_safe(obj, digits), separators=(',', ':'), sort_keys=True)
    path.write_bytes(gzip.compress(text.encode('utf-8'), mtime=0))


def _load_gz(path):
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        return json.load(f)


SERIES_DIGITS = 8


def write_series(series_dir, results, pops):
    """The committed series: what every gate and every reading is recomputed from."""
    series_dir.mkdir(parents=True, exist_ok=True)
    _dump_gz(series_dir/'populations.json.gz', pops)
    _dump_gz(series_dir/'spectrum.json.gz', {k: v['detail'] for k, v in results['spectrum'].items()}, SERIES_DIGITS)
    _dump_gz(series_dir/'V6.json.gz', [results['V6'][k] for k in sorted(results['V6'])], SERIES_DIGITS)
    _dump_gz(series_dir/'R.json.gz', results['R'], SERIES_DIGITS)
    cells = {}
    for k in sorted(results['X']):
        r = results['X'][k]
        cells.setdefault((r['name'], r['bodies']), []).append(r)
    for (name, n), runs in cells.items():
        _dump_gz(series_dir/f'X_{name}_N{n}.json.gz', runs, SERIES_DIGITS)
    return sorted(p.name for p in series_dir.iterdir())


def load_series(series_dir):
    out = dict(populations=_load_gz(series_dir/'populations.json.gz'), V6=_load_gz(series_dir/'V6.json.gz'),
               R=_load_gz(series_dir/'R.json.gz'), X={})
    for path in sorted(series_dir.glob('X_*.json.gz')):
        runs = _load_gz(path)
        out['X'][(runs[0]['name'], runs[0]['bodies'])] = runs
    return out


# EXPLORATORY, and labelled so wherever it appears. The declared mode reading fits |C_m| over 15-30 T0. The
# cold annuli grow and SATURATE before 15 T0, so that window sees only the decay and reads "absent" for a
# mode that doubled the population's dispersion. These three numbers were chosen AFTER the series were
# seen, they are not thresholds of any gate, and nothing here is a declared reading.
EXPLORATORY = dict(start_periods=6., end_at_fraction_of_peak=.5, shortest_window_periods=2.,
                   history_every_periods=1.)


def exploratory_growth(series):
    """The growth phase the declared window missed: from the end of the shot-noise imprint to the time the
    ensemble's geometric-mean amplitude first reaches half its peak, fitted per realization."""
    X = EXPLORATORY
    out = {}
    for (name, n), runs in series.items():
        live = [r for r in runs if r['live'] and r['status'] == 'completed']
        if not live:
            continue
        t = np.array([row['t_periods'] for row in live[0]['rows']])
        cell = dict(modes={})
        for m in range(1, 5):
            c = np.array([[complex(*row['coefficients']['C'][m]) for row in r['rows']] for r in live])
            c0 = np.array([[abs(complex(*row['coefficients']['C'][0])) for row in r['rows']] for r in live])
            amp = np.abs(c)/c0
            g = np.exp(np.mean(np.log(np.maximum(amp[:, 1:], 1e-300)), axis=0))      # t = 0 is a smooth field
            tg = t[1:]
            peak = int(np.argmax(g))
            t_end = float(tg[np.nonzero(g >= X['end_at_fraction_of_peak']*g[peak])[0][0]])
            row = dict(peak_relative_amplitude=float(g[peak]), peak_at_periods=float(tg[peak]),
                       relative_amplitude_at_start=float(g[np.nonzero(tg >= X['start_periods'])[0][0]]))
            if t_end - X['start_periods'] < X['shortest_window_periods']:
                row.update(growth_phase=False, why='the amplitude reaches half its peak during the imprint '
                                                   'and does not grow afterwards')
            else:
                sel = (t >= X['start_periods'] - 1e-9) & (t <= t_end + 1e-9)
                rates = [float(np.polyfit(t[sel]*T0, np.log(amp[k, sel]), 1)[0]) for k in range(len(live))]
                speeds = [float(-np.polyfit(t[sel]*T0, np.unwrap(np.angle(c[k, sel])), 1)[0]/m)
                          for k in range(len(live))]
                row.update(growth_phase=True, window_periods=[X['start_periods'], t_end],
                           growth_rate=_ensemble(rates), pattern_speed=_ensemble(speeds))
            cell['modes'][f'm{m}'] = row
        step = int(round(X['history_every_periods']/TH['X']['record_every']))
        mean_of = lambda key: [float(np.mean([r['rows'][i][key] for r in live])) for i in range(0, len(t), step)]
        c2 = np.array([[abs(complex(*row['coefficients']['C'][2]))/abs(complex(*row['coefficients']['C'][0]))
                        for row in r['rows']] for r in live])
        cell['history'] = dict(t_periods=[float(x) for x in t[::step]], radial_spread=mean_of('radial_spread'),
                               residual_radial_dispersion=mean_of('residual_radial_dispersion'),
                               radial_m2_streaming=mean_of('radial_m2_streaming'),
                               C2_over_C0_geometric_mean=[0.] + [float(x) for x in
                                                                 np.exp(np.mean(np.log(c2[:, step::step]), axis=0))])
        out.setdefault(name, {})[str(n)] = cell
    return dict(label='EXPLORATORY: chosen after the series were seen; not a declared reading and not a gate',
                rule=X, cells=out)


def evaluate_series(series, windows=None, r_window=None):
    """Everything that is computed from the archived series: two gates, and the whole of part X. The suite
    job calls this on the committed files and must get the archive back."""
    x_runs = [r for runs in series['X'].values() for r in runs]
    table = analyze_X(series['X'], windows)
    readings, temperature = readings_X(table, series['populations'], sorted({n for _, n in series['X']}))
    strip = lambda r: {k: v for k, v in r.items() if k not in ('rows', 'seconds')}
    return dict(V6=gate_V6(series['V6']), R1=gate_R1(series['R'], r_window),
                V8=gate_V8(series['V6'], series['R'], x_runs),
                run_index=dict(V6=[strip(r) for r in series['V6']], X=[strip(r) for r in x_runs]),
                part_X=dict(differences_and_modes=table, readings=readings, temperature=temperature,
                            exploratory_growth_phase=exploratory_growth(series['X'])))


AMENDMENT = 'protocol-rut7-amendment-1.md'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--canonical', action='store_true')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 4))
    ap.add_argument('--resume-from', type=Path, help='reuse finished tasks of an interrupted campaign that '
                    'ran the same code; rut7_checks.py replays prefixes independently')
    ap.add_argument('--trust-cache', action='store_true', help='reuse finished tasks even though the task code '
                    'has changed since they ran; recorded in the archive, and only for a change that cannot '
                    'touch them')
    ap.add_argument('--smoke', action='store_true', help='a few minutes at toy sizes, to exercise every path; '
                    'never archived')
    args = ap.parse_args()
    if args.smoke and args.canonical:
        raise SystemExit('a smoke run is never archived')
    t0 = time.time()
    out_dir = evidence_io.output_dir(args, 'path-memory-rut7')
    args.output_dir = out_dir
    task_dir = out_dir/'tasks'
    task_dir.mkdir(parents=True, exist_ok=True)
    code = {f: hashlib.sha256(_find(f).read_bytes()).hexdigest() for f in CODE}
    task_code = {f: code[f] for f in TASK_CODE}         # a finished task depends on these, not on rut7.py
    owner_dir = PO.repository_root()/'research_work/annulus_sampling'
    owner = {f: hashlib.sha256((owner_dir/f).read_bytes()).hexdigest() for f in OWNER_PACKAGE}
    log = open(out_dir/'rut7-progress.log', 'a', encoding='utf-8')

    def say(msg):
        line = f'[{time.time() - t0:7.0f}s] {msg}'
        print(line, flush=True)
        log.write(line + '\n')
        log.flush()

    def cached(kind, key):
        if not args.resume_from:
            return None
        src = args.resume_from/'tasks'/(f'{kind}__{key}'.replace(':', '_') + '.json.gz')
        if not src.is_file():
            return None
        stored = _load_gz(src)
        if stored.get('code') != task_code:
            if not args.trust_cache:
                return None
            trusted.append(f'{kind}:{key}')
        return stored['result']

    def keep(kind, key, result):
        _dump_gz(task_dir/(f'{kind}__{key}'.replace(':', '_') + '.json.gz'), dict(code=task_code, result=result))

    results, reused, walls, errors, trusted = {}, [], {}, [], []

    def drain(pool, tasks):
        futures = {}
        for _, kind, key, a in tasks:
            hit = cached(kind, key)
            if hit is not None:
                results.setdefault(kind, {})[key] = hit
                keep(kind, key, hit)
                reused.append(f'{kind}:{key}')
            else:
                futures[pool.submit(run_task, kind, key, a)] = (kind, key)
        for n, fut in enumerate(as_completed(futures), 1):
            try:
                kind, key, out, wall = fut.result()
            except Exception as err:                 # let every other task finish and be cached first
                errors.append(f'{futures[fut]}: {err!r}')
                say(f'{n:4d}/{len(futures)} FAILED {futures[fut]}: {err!r}')
                continue
            results.setdefault(kind, {})[key] = out
            walls[f'{kind}:{key}'] = wall
            keep(kind, key, out)
            status = out.get('status', out.get('passed', ''))
            say(f'{n:4d}/{len(futures)} {kind:<11}{key:<44}{str(status):<12}{wall:8.1f}s')

    specs = {p['name']: p for p in POPULATIONS}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        say(f'stage 7 campaign: {args.workers} workers, smoke={args.smoke}, output {out_dir}')
        drain(pool, [(1, 'population', n, (specs[n],)) for n in specs if specs[n]['alpha_from'] is None])
        alpha_B = results['population']['B_mid']['description']['alpha']
        drain(pool, [(1, 'population', n, (specs[n], alpha_B)) for n in specs if specs[n]['alpha_from']])
        pops = {n: results['population'][n] for n in specs}
        tasks = declared_tasks(pops, args.smoke)
        say(f'{len(tasks)} tasks declared; {len(reused)} reused so far')
        drain(pool, tasks)
    if errors:
        raise SystemExit('tasks failed, nothing assembled; finished tasks are cached for --resume-from: '
                         + '; '.join(errors))

    series_dir = out_dir/'rut7-series'
    files = write_series(series_dir, results, pops)
    series = load_series(series_dir)
    scale = (4./TH['X']['horizon']) if args.smoke else 1.
    windows = None if not args.smoke else {k: [scale*x for x in TH['X'][w]] for k, w in
                                           (('early', 'early_window'), ('late', 'late_window'),
                                            ('mode', 'mode_fit_window'))}
    from_series = evaluate_series(series, windows, [2., 6.] if args.smoke else None)

    spectra = [results['spectrum'][k] for k in sorted(results['spectrum'],
                                                      key=lambda s: (RUNGS.index(s.split(':')[0]), int(s.split(':m')[1])))]
    failures, worst = gate_E3_rows(spectra)
    e3_control = results['E3_control']['E3_control']
    gates = dict(E1_known_spectra=results['E1']['E1'], E2_raw_robustness=results['E2']['E2'],
                 E3_completeness=dict(failures=failures, worst=worst, sub_rectangles=int(sum(
                     s['base']['rectangles'] + s['shifted']['rectangles'] for s in spectra)),
                     negative_control=e3_control, negative_control_rejected=e3_control['rejected'],
                     passed=bool(not failures and e3_control['rejected'])),
                 E4_edge_handling=results['E4']['E4'])
    for g, kind in (('V1_the_measure', 'V1'), ('V2_drawn_moments', 'V2'), ('V5_drawn_source_writes_the_field', 'V5')):
        rows = results[kind]
        gates[g] = dict(rows=rows, negative_controls_rejected=bool(all(r['negative_control']['rejected']
                                                                       for r in rows.values())),
                        passed=bool(all(r['passed'] for r in rows.values())))
    v3 = results['V3']
    warm = [n for n in v3 if specs[n]['dE'] >= .03]
    gates['V3_support_fits'] = dict(rows=v3, negative_control=dict(
        what='the stage 6 five-width support, on the warm annuli', populations=warm,
        rejected=bool(all(not v3[n]['five_width_fits'] for n in warm))),
        passed=bool(all(r['fits'] for r in v3.values()) and all(not v3[n]['five_width_fits'] for n in warm)))
    refined = {tuple(k.split(':')): v for k, v in results['V4'].items()}
    v4_rows, worst_i, worst_p, v4_ok = gate_V4_rows(pops, refined)
    ctl = results['V4_control']
    moved = {f'dE{dE}': abs(ctl[f"dE{dE}:n{REFINEMENTS['radial_grid']['n_r']}"]['alpha']
                            / ctl[f"dE{dE}:n{GRID['n_r']}"]['alpha'] - 1) for dE in (.010, .030)}
    v4_control = dict(what="stage 6's label -- support read at the argmax node, unbracketed secant -- under "
                           'the same radial refinement', alpha_relative_movement=moved,
                      rejected=bool(all(v > TH['S']['convergence_integrated'] for v in moved.values())))
    gates['V4_discretization'] = dict(rows=v4_rows, worst_integrated=worst_i, worst_profile=worst_p,
                                      negative_control=v4_control, passed=bool(v4_ok and v4_control['rejected']))
    gates['V6_draw_is_stationary'] = from_series['V6']
    gates['V7_restart'] = results['V7']['V7']
    gates['V8_ensemble'] = from_series['V8']
    gates['R1_simulator_against_linear_theory'] = from_series['R1']
    failed = sorted(k for k, g in gates.items() if not g['passed'])

    unstable = {}
    for s in spectra:
        rows = [r for r in s['roots'] if r['kind'] == 'unstable']
        if rows:
            best = max(rows, key=lambda r: r['re'])
            unstable.setdefault(s['model'], {})[f"m{s['m']}"] = dict(
                unstable_roots=len(rows), fastest=[best['re'], best['im']],
                e_folding_periods=best['e_folding_periods'],
                pattern_speed_over_omega=best['pattern_speed_over_omega'])
    result = dict(
        experiment='RUT-1 stage 7: the spectrum certified to the axis, the population drawn faithfully, and '
                   'live compared against frozen with a measured uncertainty',
        protocol='protocol-rut7.md', amendment=AMENDMENT, smoke=bool(args.smoke),
        statuses=dict(
            reproduction='decided by rut7_checks.py against this archive, and by this driver when re-run',
            numerical_verification=dict(passed=not failed, failed_gates=failed),
            scientific_outcome='part_X below: archived as it fell, and never part of the exit status'),
        deviations=[],
        thresholds=TH,
        frozen_equation=dict(two_stage='tau_form dE/dt = S - E, dC/dt = E - C/tau_keep', tau_keep_T0=10.,
                             tau_form_T0=3., footprint_w=W),
        part_E=dict(strip=dict(re=TH['E']['strip_re'], im=TH['E']['strip_im'], m_max=TH['E']['m_max']),
                    unstable_summary=unstable,
                    spectra={f"{s['model']}:m{s['m']}": {k: v for k, v in s.items() if k != 'detail'}
                             for s in spectra}),
        part_S=dict(populations={n: {k: v for k, v in p.items() if k != 'state'} for n, p in pops.items()},
                    grid=GRID, refinements=REFINEMENTS, sampler_of_record='research_work/annulus_sampling '
                    '(annulus-phase-space-v2)'),
        gates=gates, run_index=from_series['run_index'], part_X=from_series['part_X'],
        code_sha256_at_launch=code, owner_package_sha256=owner, tasks_reused=sorted(reused),
        tasks_reused_across_a_code_change=sorted(trusted),
        series_sha256={f: hashlib.sha256((series_dir/f).read_bytes()).hexdigest() for f in files},
        input_sha256={f: hashlib.sha256(_find(f).read_bytes()).hexdigest() for f in ('protocol-rut7.md', AMENDMENT)},
        what_this_is_not='thirty periods, one planar annulus about a fixed centre, an instantaneous Gaussian '
                         'kernel, a declared mass and a labelled writing rate. A population quiet for thirty '
                         'periods may hold a mode with a growth time of hundreds; only the linear mode '
                         'calculation around these populations can exclude one, and it is not run here.')
    text = json.dumps(_safe(result), indent=1) + '\n'
    (out_dir/'rut7-log.json').write_text(json.dumps(dict(task_wall_seconds=walls,
                                                         wall_seconds=round(time.time() - t0, 1)), indent=1),
                                         encoding='utf-8')
    say(f'numerical verification: {"PASSED" if not failed else "FAILED " + ", ".join(failed)}')
    for k, g in sorted(gates.items()):
        say(f"  {'PASS' if g['passed'] else 'FAIL'} {k}")
    canonical = HERE/'rut7-results.json'
    if args.smoke:
        (out_dir/'rut7-results.json').write_text(text, encoding='utf-8', newline='\n')
        return 0 if not failed else 1
    if args.canonical:
        canon = HERE/'rut7-series'
        canon.mkdir(exist_ok=True)
        for f in files:
            (canon/f).write_bytes((series_dir/f).read_bytes())
        if not canonical.is_file():
            canonical.write_text(text, encoding='utf-8', newline='\n')
    if not canonical.is_file():
        (out_dir/'rut7-results.json').write_text(text, encoding='utf-8', newline='\n')
        say('no archive yet: results written to the output directory only')
        return 0 if not failed else 1
    status = evidence_io.finish(args, 'path-memory-rut7', text, canonical)
    return status if not failed else 1


if __name__ == '__main__':
    raise SystemExit(main())
