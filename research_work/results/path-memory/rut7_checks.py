"""Suite job for RUT-1 stage 7 (protocol-rut7.md, amendment 1). rut7.py is the campaign driver and takes of
the order of an hour on many cores; this job checks, in a couple of minutes, that its archive is what the
committed code produces, and reports the three statuses separately:

REPRODUCTION
- the committed code, the owner's sampler package, the protocol and its amendment hash to the values
  recorded when the campaign was launched, and every committed series file hashes to its recorded value;
- the thresholds in the archive are the ones the protocol declares today;
- the cheap gates are recomputed from scratch and must reproduce the archive: E1, E2, E4 and the E3 control,
  the two-stage spectrum at m = 1 and 2, all five populations, V1 and V3 on all five, V2 and one V4
  refinement on one population, one V5 population, V7, and the declared part R run;
- deterministic 2 T0 prefixes of one V6 run and of one live and frozen pair are replayed and must
  reproduce the archived series;
- everything that is computed from the series -- gates V6, V8 and R1, and the whole of part X with its
  declared readings -- is recomputed from the committed files and must reproduce the archive.

NUMERICAL VERIFICATION
- every gate in the archive passed and rejected its negative control, and so did every gate recomputed here.

SCIENTIFIC OUTCOME
- reported, and never part of the exit status.

The exit status is non-zero if reproduction or numerical verification fails.
"""
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut7 as R7  # noqa: E402
from rut7_tasks import PO, POPULATIONS, R_CASES, TH, _find, _safe, evidence_io, run_task, thresholds  # noqa: E402

PREFIX = 2.                 # periods replayed
SERIES_RTOL = 1e-6          # the archived series are rounded to eight significant figures
REPLAY_POPULATION = 'A_cold'


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _rows_match(fresh, stored, skip=('seconds',)):
    """A replayed prefix against the archived rows, to the archive's rounding."""
    diffs = []
    for i, (a, b) in enumerate(zip(fresh, stored)):
        a = {k: v for k, v in a.items() if k not in skip}
        b = {k: v for k, v in b.items() if k not in skip}
        diffs += evidence_io.compare(a, b, f'/rows[{i}]', rtol=SERIES_RTOL, atol=1e-12)
    if len(fresh) > len(stored):
        diffs.append(f'replay has {len(fresh)} rows, archive {len(stored)}')
    return diffs


def _light(d):
    """Drop timings, which are not results."""
    if isinstance(d, dict):
        return {k: _light(v) for k, v in d.items() if k not in ('seconds', 'wall_seconds')}
    if isinstance(d, list):
        return [_light(v) for v in d]
    return d


def main():
    t0 = time.time()
    archive = json.loads((HERE/'rut7-results.json').read_text(encoding='utf-8'))
    series_dir = HERE/'rut7-series'
    out = {}

    # ------------------------------------------------------------ reproduction: provenance
    code = {f: _sha(_find(f)) for f in archive['code_sha256_at_launch']}
    owner_dir = PO.repository_root()/'research_work/annulus_sampling'
    owner = {f: _sha(owner_dir/f) for f in archive['owner_package_sha256']}
    inputs = {f: _sha(_find(f)) for f in archive['input_sha256']}
    mism = ([f for f in code if code[f] != archive['code_sha256_at_launch'][f]]
            + [f'annulus_sampling/{f}' for f in owner if owner[f] != archive['owner_package_sha256'][f]]
            + [f for f in inputs if inputs[f] != archive['input_sha256'][f]])
    out['code_and_inputs_are_what_ran'] = dict(mismatched=mism, passed=not mism)
    bad = [f for f, h in archive['series_sha256'].items() if _sha(series_dir/f) != h]
    out['series_hashes'] = dict(files=len(archive['series_sha256']), mismatched=bad, passed=not bad)
    out['thresholds_are_the_protocols'] = dict(passed=bool(thresholds() == archive['thresholds'] == TH))

    # ------------------------------------------------------------ reproduction: recomputed from scratch
    series = R7.load_series(series_dir)
    pops = series['populations']
    specs = {p['name']: p for p in POPULATIONS}
    alpha_B = pops['B_mid']['description']['alpha']
    pop = pops[REPLAY_POPULATION]
    base_case = next(c for c in R_CASES if c['name'] == 'w5_h0.01')
    jobs = [('E1', 'E1', ()), ('E2', 'E2', ()), ('E3_control', 'E3_control', ()), ('E4', 'E4', ()),
            ('spectrum', 'two_stage:m1', ('two_stage', 1)), ('spectrum', 'two_stage:m2', ('two_stage', 2)),
            ('V2', REPLAY_POPULATION, (pop,)), ('V5', REPLAY_POPULATION, (pop,)),
            ('V4', f'{REPLAY_POPULATION}:radial_grid', (specs[REPLAY_POPULATION], 'radial_grid', None)),
            ('V7', 'V7', (pops['B_mid'],)), ('R', base_case['name'], (base_case,)),
            ('V6', 'prefix', (pop, 0, None, PREFIX)),
            ('X', 'prefix:live', (pop, TH['X']['bodies'][0], 0, True, PREFIX)),
            ('X', 'prefix:frozen', (pop, TH['X']['bodies'][0], 0, False, PREFIX))]
    jobs += [('population', n, (specs[n],) if specs[n]['alpha_from'] is None else (specs[n], alpha_B))
             for n in specs]
    jobs += [(kind, n, (pops[n],)) for n in specs for kind in ('V1', 'V3')]
    fresh = {}
    with ProcessPoolExecutor(max_workers=6) as pool:
        for kind, key, result, _ in pool.map(run_task, *zip(*jobs)):
            fresh[(kind, key)] = result

    gates = archive['gates']
    checks = {}
    for kind, name in (('E1', 'E1_known_spectra'), ('E2', 'E2_raw_robustness'), ('E4', 'E4_edge_handling'),
                       ('V7', 'V7_restart')):
        checks[name] = evidence_io.compare(_light(fresh[(kind, kind)]), _light(gates[name]))
    checks['E3_control'] = evidence_io.compare(fresh[('E3_control', 'E3_control')],
                                               gates['E3_completeness']['negative_control'])
    for m in (1, 2):
        got = {k: v for k, v in fresh[('spectrum', f'two_stage:m{m}')].items() if k != 'detail'}
        checks[f'spectrum_two_stage_m{m}'] = evidence_io.compare(got, archive['part_E']['spectra'][f'two_stage:m{m}'],
                                                                 rtol=1e-7, atol=1e-9)
    for n in specs:
        got = fresh[('population', n)]
        checks[f'population_{n}'] = (
            evidence_io.compare({k: v for k, v in got.items() if k != 'state'}, archive['part_S']['populations'][n])
            + evidence_io.compare(got['state'], pops[n]['state']))
        checks[f'V1_{n}'] = evidence_io.compare(fresh[('V1', n)], gates['V1_the_measure']['rows'][n])
        checks[f'V3_{n}'] = evidence_io.compare(fresh[('V3', n)], gates['V3_support_fits']['rows'][n])
    checks['V2'] = evidence_io.compare(fresh[('V2', REPLAY_POPULATION)],
                                       gates['V2_drawn_moments']['rows'][REPLAY_POPULATION])
    checks['V5'] = evidence_io.compare(fresh[('V5', REPLAY_POPULATION)],
                                       gates['V5_drawn_source_writes_the_field']['rows'][REPLAY_POPULATION])
    ref = fresh[('V4', f'{REPLAY_POPULATION}:radial_grid')]['description']
    base = pop['description']
    moved = dict(alpha=abs(ref['alpha']/base['alpha'] - 1), mean_radius=abs(ref['mean_radius']/base['mean_radius'] - 1),
                 radial_width=abs(ref['radial_width']/base['radial_width'] - 1),
                 peak_field=abs(ref['peak_field']/base['peak_field'] - 1),
                 mean_support=abs(ref['support']['mean']/base['support']['mean'] - 1))
    stored = gates['V4_discretization']['rows'][f'{REPLAY_POPULATION}:radial_grid']['integrated']
    checks['V4_radial_grid'] = evidence_io.compare(moved, stored, rtol=1e-6, atol=1e-10)
    fresh_R = R7.measure_R(_safe(fresh[('R', base_case['name'])], R7.SERIES_DIGITS), TH['R']['fit_window'])
    checks['R_base_case'] = evidence_io.compare(
        fresh_R, gates['R1_simulator_against_linear_theory']['rows'][base_case['name']]['measured'], rtol=1e-6)

    # ------------------------------------------------------------ reproduction: deterministic prefixes
    v6 = next(r for r in series['V6'] if r['name'] == REPLAY_POPULATION and r['control'] is None
              and r['realization'] == 0)
    checks['prefix_V6'] = _rows_match(_safe(fresh[('V6', 'prefix')]['rows'], R7.SERIES_DIGITS), v6['rows'])
    cell = series['X'][(REPLAY_POPULATION, TH['X']['bodies'][0])]
    for live in (True, False):
        tag = 'live' if live else 'frozen'
        stored_run = next(r for r in cell if r['realization'] == 0 and r['live'] == live)
        got = fresh[('X', f'prefix:{tag}')]
        checks[f'prefix_X_{tag}'] = (_rows_match(_safe(got['rows'], R7.SERIES_DIGITS), stored_run['rows'])
                                     + ([] if got['initial_sha256'] == stored_run['initial_sha256']
                                        else ['the replay drew different bodies']))

    # ------------------------------------------------------------ reproduction: everything from the series
    again = _safe(R7.evaluate_series(series))
    checks['gates_from_series'] = (evidence_io.compare(again['V6'], gates['V6_draw_is_stationary'])
                                   + evidence_io.compare(again['V8'], gates['V8_ensemble'])
                                   + evidence_io.compare(again['R1'], gates['R1_simulator_against_linear_theory']))
    checks['run_index'] = evidence_io.compare(again['run_index'], archive['run_index'])
    checks['part_X'] = evidence_io.compare(again['part_X'], archive['part_X'])
    out['recomputed'] = {k: dict(differences=v[:5], count=len(v), passed=not v) for k, v in checks.items()}

    # ------------------------------------------------------------ the three statuses
    reproduction = bool(out['code_and_inputs_are_what_ran']['passed'] and out['series_hashes']['passed']
                        and out['thresholds_are_the_protocols']['passed']
                        and all(c['passed'] for c in out['recomputed'].values()))
    recomputed_gates = {k: bool(fresh[(k, k)]['passed']) for k in ('E1', 'E2', 'E4', 'V7')}
    recomputed_gates['E3_control'] = bool(fresh[('E3_control', 'E3_control')]['rejected'])
    recomputed_gates.update({f'V1_{n}': bool(fresh[('V1', n)]['passed']) for n in specs})
    recomputed_gates.update({f'V3_{n}': bool(fresh[('V3', n)]['fits']) for n in specs})
    recomputed_gates.update(V2=bool(fresh[('V2', REPLAY_POPULATION)]['passed']),
                            V5=bool(fresh[('V5', REPLAY_POPULATION)]['passed']),
                            V6=bool(again['V6']['passed']), V8=bool(again['V8']['passed']),
                            R1=bool(again['R1']['passed']))
    archived = archive['statuses']['numerical_verification']
    numerical = bool(archived['passed'] and all(recomputed_gates.values()))
    readings = archive['part_X']['readings']
    out['statuses'] = dict(
        reproduction=dict(passed=reproduction),
        numerical_verification=dict(passed=numerical, archive_failed_gates=archived['failed_gates'],
                                    recomputed_here={k: v for k, v in recomputed_gates.items() if not v} or 'all passed',
                                    archive_gates={k: bool(g['passed']) for k, g in gates.items()}),
        scientific_outcome=dict(
            affects_exit_status=False,
            collective_m2_declared_reading={n: r['collective_mode'].get('m2', {}).get('reading')
                                            for n, r in readings.items()},
            caution='the declared window, 15-30 T0, lies AFTER the cold annuli have grown and saturated, so '
                    '"absent" there describes the decay, not the mode. The growth phase below is exploratory: '
                    'its rule was chosen after the series were seen.',
            exploratory_m2_growth_rate={
                n: {count: (c['modes']['m2']['growth_rate']['mean'] if c['modes']['m2']['growth_phase'] else None)
                    for count, c in cells.items()}
                for n, cells in archive['part_X']['exploratory_growth_phase']['cells'].items()},
            quiet={n: {q: v['reading'] for q, v in r['quiet'].items()} for n, r in readings.items()},
            discreteness={n: (r['discreteness'] or {}).get('reading') for n, r in readings.items()}))
    out['deviations'] = archive['deviations']
    out['passed'] = bool(reproduction and numerical)
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
