"""Suite job for RUT-1 stage 8 (protocol-rut8.md). rut8.py is the campaign driver, in two phases; this job
checks in a few minutes that its two archives are what the committed code produces, and reports the three
statuses separately.

REPRODUCTION
- the committed code, the protocol and every series file hash to the values recorded at launch, in both
  archives; the measurement archive records the hash of the predictions file it was compared with, and that
  file is the committed one; the thresholds in both archives are the ones the protocol declares today;
- recomputed from scratch and compared with the predictions archive: gate L1 on every population and its
  negative control, L5's negative control, the sub-rectangle holding the m = 2 root of `A_cold`, of `B_cold`
  and of each new population predicted unstable, one member of the family scan, and the first of the ring-limit
  annuli with its exact-ring target;
- recomputed from the archive: the threshold fit, and BOTH ring-limit gates -- L3 exactly as declared, which
  failed, and L3b as amended, its replacement -- re-derived from the archived rows;
- a deterministic 2 T0 prefix of one live and frozen pair of a new population is replayed against the series;
- every part M reading is re-derived from the committed series and the committed predictions.

NUMERICAL VERIFICATION
- every gate in both archives passed and rejected its negative control, and so did every gate recomputed here.

SCIENTIFIC OUTCOME
- the predictions, the measurements and the declared readings: reported, never part of the exit status.
"""
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut8 as R8  # noqa: E402
from rut8_tasks import TH, _find, _safe, run_task, stage7_populations, thresholds  # noqa: E402
import rut7 as S7  # noqa: E402
from rut7_tasks import evidence_io  # noqa: E402

PREFIX = 2.


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _holding(modes_rows, name, m=2):
    rows = [r for r in modes_rows if r['name'] == name and r['m'] == m and not r['shifted'] and r['roots']]
    return max(rows, key=lambda r: max(q['re'] for q in r['roots'])) if rows else None


def main():
    t0 = time.time()
    pred = json.loads((HERE/'rut8-predictions.json').read_text(encoding='utf-8'))
    meas = json.loads((HERE/'rut8-results.json').read_text(encoding='utf-8'))
    series_dir = HERE/'rut8-series'
    out = {}

    mism = []
    for label, archive in (('predictions', pred), ('results', meas)):
        mism += [f'{label}:{f}' for f, h in archive['code_sha256_at_launch'].items() if _sha(_find(f)) != h]
        mism += [f'{label}:{f}' for f, h in archive['series_sha256'].items() if _sha(series_dir/f) != h]
        mism += [f'{label}:{f}' for f, h in archive['input_sha256'].items() if _sha(_find(f)) != h]
    out['code_inputs_and_series_are_what_ran'] = dict(mismatched=mism, passed=not mism)
    out['thresholds_are_the_protocols'] = dict(passed=bool(thresholds() == pred['thresholds'] == meas['thresholds'] == TH))

    modes_rows = S7._load_gz(series_dir/'modes.json.gz')
    new_pops = S7._load_gz(series_dir/'new-populations.json.gz')
    pops = dict(stage7_populations(), **new_pops)
    unstable_new = [n for n, p in pred['part_P']['new_populations'].items() if p['predicted_unstable']]
    jobs = [('L1', n, (pops[n],)) for n in pops] + [('L1', 'control:A_cold', (pops['A_cold'], False))]
    rect_jobs = {}
    for name in ['A_cold', 'B_cold'] + unstable_new:
        row = _holding(modes_rows, name)
        rect_jobs[name] = row
        jobs.append(('rectangle', f'{name}:m2', (pops[name], 2, row['rect'], False)))
    jobs.append(('L5_control', 'L5_control', (pops['B_cold'], rect_jobs['B_cold']['rect'])))
    family_dE = TH['P']['family_B_dE'][len(TH['P']['family_B_dE'])//2]
    jobs.append(('family', f'{family_dE:.3f}', (family_dE,)))
    jobs.append(('L3', 'first', (TH['L']['ring_limit'][0],)))
    jobs.append(('population', unstable_new[0], (unstable_new[0], TH['P']['new_populations'][unstable_new[0]])))
    n_small = TH['M']['bodies'][0]
    for live in (True, False):
        jobs.append(('X', f"prefix:{'live' if live else 'frozen'}", (new_pops[unstable_new[0]], n_small, 0, live, PREFIX)))
    fresh = {}
    with ProcessPoolExecutor(max_workers=6) as pool:
        for kind, key, result, _ in pool.map(run_task, *zip(*jobs)):
            fresh[(kind, key)] = result

    L = pred['part_L']['gates']
    checks = {}
    drop = lambda d: {k: v for k, v in d.items() if k != 'seconds'}
    for n in pops:
        checks[f'L1_{n}'] = evidence_io.compare(fresh[('L1', n)], L['L1_orbit_library']['rows'][n], rtol=1e-6, atol=1e-10)
    checks['L1_control'] = evidence_io.compare(fresh[('L1', 'control:A_cold')], L['L1_orbit_library']['negative_control']['row'],
                                               rtol=1e-6, atol=1e-10)
    for name, row in rect_jobs.items():
        keep = lambda r: dict(roots=[dict(re=q['re'], im=q['im']) for q in r['roots']], located=r['located'],
                              winding=r['winding']['count'], counts_agree=r['counts_agree'])
        checks[f'root_{name}'] = evidence_io.compare(keep(fresh[('rectangle', f'{name}:m2')]), keep(row), rtol=1e-7, atol=1e-9)
    checks['L5_control'] = evidence_io.compare(fresh[('L5_control', 'L5_control')], L['L5_completeness']['negative_control'])
    stored_family = next(r for r in pred['part_P']['family_B'] if abs(r['dE'] - family_dE) < 1e-12)
    got_family = fresh[('family', f'{family_dE:.3f}')]
    checks['family_member'] = evidence_io.compare(
        dict(roots=[dict(re=q['re'], im=q['im']) for q in got_family['roots']], sigma_r=got_family['sigma_r'],
             mean_support=got_family['mean_support']),
        dict(roots=[dict(re=q['re'], im=q['im']) for q in stored_family['roots']], sigma_r=stored_family['sigma_r'],
             mean_support=stored_family['mean_support']), rtol=1e-7, atol=1e-9)
    first = sorted(L['L3_ring_limit']['rows'], key=lambda r: -r['radial_width'])[0]
    checks['ring_limit_first'] = evidence_io.compare(drop(fresh[('L3', 'first')]), drop(first), rtol=1e-6, atol=1e-9)
    name0 = unstable_new[0]
    got_pop = fresh[('population', name0)]
    checks['new_population'] = (evidence_io.compare({k: v for k, v in got_pop.items() if k != 'state'}, pred['populations'][name0])
                                + evidence_io.compare(got_pop['state'], new_pops[name0]['state']))
    checks['threshold_fit'] = evidence_io.compare(_safe(R8.fit_threshold(pred['part_P']['family_B'])), pred['part_P']['threshold'])
    # both ring-limit gates re-derived from the archived rows: L3 exactly as declared, L3b as amended
    import numpy as np
    for gate, last, lo, hi in (('L3_ring_limit', TH['L']['ring_limit_last'], *TH['L']['ring_limit_order']),
                               ('L3b_ring_limit_amended', R8.TH_A['L3b']['last'], R8.TH_A['L3b']['order_min'], float('inf'))):
        rows = L[gate]['rows']
        d = [r['growth_relative_difference'] for r in rows]
        order = float(np.polyfit(np.log([r['radial_width'] for r in rows]), np.log(d), 1)[0])
        passed = bool(all(b < a for a, b in zip(d[:-1], d[1:])) and d[-1] < last and lo <= order <= hi
                      and L[gate]['negative_control']['rejected'])
        checks[f'{gate}_rederived'] = evidence_io.compare(dict(order=order, passed=passed),
                                                         dict(order=L[gate]['observed_order'], passed=L[gate]['passed']))

    series = R8.load_series(series_dir)
    cell = series['X'][(name0, n_small)]
    for live in (True, False):
        tag = 'live' if live else 'frozen'
        stored = next(r for r in cell if r['realization'] == 0 and r['live'] == live)
        got = fresh[('X', f'prefix:{tag}')]
        diffs = []
        for i, (a, b) in enumerate(zip(_safe(got['rows'], S7.SERIES_DIGITS), stored['rows'])):
            diffs += evidence_io.compare(a, b, f'/rows[{i}]', rtol=1e-6, atol=1e-12)
        if got['initial_sha256'] != stored['initial_sha256']:
            diffs.append('the replay drew different bodies')
        checks[f'prefix_X_{tag}'] = diffs
    again = _safe(R8.evaluate_measurement(series, pred))
    M = meas['part_M']
    checks['V6_from_series'] = evidence_io.compare(again['V6'], M['gates']['V6_draw_is_stationary'])
    checks['comparison_from_series'] = evidence_io.compare(again['comparison'], M['comparison'])
    checks['run_index'] = evidence_io.compare(again['run_index'], M['run_index'])
    out['recomputed'] = {k: dict(differences=v[:5], count=len(v), passed=not v) for k, v in checks.items()}
    out['predictions_file_is_the_one_measured_against'] = dict(
        passed=bool(meas['input_sha256'].get('rut8-predictions.json') == _sha(HERE/'rut8-predictions.json')))

    reproduction = bool(all(v['passed'] for k, v in out.items() if k != 'recomputed')
                        and all(c['passed'] for c in out['recomputed'].values()))
    recomputed = {f'L1_{n}': bool(fresh[('L1', n)]['passed']) for n in pops}
    recomputed['L1_control_rejected'] = bool(not fresh[('L1', 'control:A_cold')]['passed'])
    recomputed['L5_control_rejected'] = bool(fresh[('L5_control', 'L5_control')]['rejected'])
    recomputed['V6'] = bool(again['V6']['passed'])
    numerical = bool(pred['statuses']['numerical_verification']['passed'] and meas['statuses']['numerical_verification']['passed']
                     and all(recomputed.values()))
    out['statuses'] = dict(
        reproduction=dict(passed=reproduction),
        numerical_verification=dict(passed=numerical,
                                    archive_failed_gates=pred['statuses']['numerical_verification']['failed_gates']
                                    + meas['statuses']['numerical_verification']['failed_gates'],
                                    failed_as_declared=pred['statuses']['numerical_verification']['failed_as_declared'],
                                    superseded_by_amendment_1=pred['statuses']['numerical_verification']['superseded_by_amendment_1'],
                                    recomputed_here={k: v for k, v in recomputed.items() if not v} or 'all passed'),
        scientific_outcome=dict(affects_exit_status=False,
                                threshold_dE=pred['part_P']['threshold'].get('threshold_dE'),
                                readings={n: c['reading'] for n, c in M['comparison'].items()},
                                predicted_growth_rate={n: c['prediction']['growth_rate'] for n, c in M['comparison'].items()},
                                measured_growth_rate_largest_count={
                                    n: (c['measured'][str(max(TH['M']['bodies']))].get('growth_rate') or {}).get('mean')
                                    for n, c in M['comparison'].items()}))
    out['deviations'] = pred['deviations'] + meas['deviations']
    out['passed'] = bool(reproduction and numerical)
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
