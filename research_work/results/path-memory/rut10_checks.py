"""Suite job for RUT-1 stage 10 (protocol-rut10.md). rut10.py is the campaign driver; this job checks in a few
minutes that its archive is what the committed code produces, and reports the three statuses separately.

REPRODUCTION
- the committed code, both protocols, stage 8's predictions, stage 9's archive and every series file hash to the
  values recorded at launch; the archive was made from the committed protocol; the thresholds are today's;
- THE EXACT PREPARATION IS REPLAYED for one seeded run: the quiet sample is drawn again, moved again by the
  predicted mode's field, and the field is seeded again; the bodies and the complete field must hash to the
  archived values, and the bodies must equal the archived initial state; a declared prefix of that run and of
  the unseeded quiet run is then evolved and compared with the series, record for record;
- the saved complete initial fields hash to the values every run recorded;
- recomputed from the series: every fit, every gate and the declared reading.

NUMERICAL VERIFICATION
- every gate in the archive passed and rejected its negative control, and so did every gate recomputed here.

SCIENTIFIC OUTCOME
- the seeded growth rates and pattern speeds against the forecast of record: reported, never part of the exit
  status.
"""
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut10 as R10  # noqa: E402
from rut10_tasks import PROTOCOL, TH, _find, _safe, _sha, load_arrays, np, run_task, thresholds  # noqa: E402
import rut7 as S7  # noqa: E402
import rut9_tasks as R9  # noqa: E402
from rut7_tasks import evidence_io  # noqa: E402


def _file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load_runs(series_dir, index):
    runs = {}
    for key in index:
        f = 'run_' + key.replace(':', '_').replace('+', 'p').replace('-', 'n') + '.json.gz'
        runs[key] = S7._load_gz(series_dir/f)
    return runs


def main():
    t0 = time.time()
    arc = json.loads((HERE/'rut10-results.json').read_text(encoding='utf-8'))
    series_dir = HERE/'rut10-series'
    out = {}
    mism = [f for f, h in arc['code_sha256_at_launch'].items() if _file_sha(_find(f)) != h]
    mism += [f for f, h in arc['series_sha256'].items() if _file_sha(series_dir/f) != h]
    mism += [f for f, h in arc['input_sha256'].items() if _file_sha(_find(f)) != h]
    out['code_inputs_and_series_are_what_ran'] = dict(mismatched=mism, passed=not mism)
    out['made_from_the_committed_protocol'] = dict(read_from=Path(arc['protocol_read_from']).name,
                                                   passed=bool(Path(arc['protocol_read_from']).name == PROTOCOL))
    out['thresholds_are_the_protocols'] = dict(passed=bool(thresholds() == arc['thresholds'] == TH))

    G = TH['G']
    pops = R9.archived_populations()
    runs = _load_runs(series_dir, arc['run_index'])
    name, q = G['populations'][0], G['quarters'][0]
    seeded_key = R10._key(name, q, 0, G['targets'][0], 'base')
    quiet_key = R10._key(name, q, 0, 0., 'base')
    prefix = TH['C']['replay_prefix']
    jobs = [('run', seeded_key, (pops[name], q, 0, G['targets'][0], 'base', prefix)),
            ('run', quiet_key, (pops[name], q, 0, 0., 'base', prefix)),
            ('Q', f'{name}:q{q}:k0', (pops[name], q, 0)),
            ('field_state', f"{name}:t{G['targets'][0]:+.0e}", (pops[name], G['targets'][0]))]
    fresh = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        for kind, key, result, _ in pool.map(run_task, *zip(*jobs)):
            fresh[(kind, key)] = result

    checks = {}
    for key in (seeded_key, quiet_key):
        got, stored = fresh[('run', key)], runs[key]
        diffs = [f'{k} differs' for k in ('sample_sha256', 'initial_sha256', 'field_sha256') if got[k] != stored[k]]
        diffs += evidence_io.compare(got['initial_positions'], stored['initial_positions'], '/initial_positions', rtol=0, atol=0)
        diffs += evidence_io.compare(got['initial_velocities'], stored['initial_velocities'], '/initial_velocities', rtol=0, atol=0)
        for i, (a, b) in enumerate(zip(got['rows'], stored['rows'])):
            diffs += evidence_io.compare(a, b, f'/rows[{i}]', rtol=1e-6, atol=1e-13)
        checks[f'replay_{key}'] = diffs
    fields = S7._load_gz(series_dir/'initial-fields.json.gz')
    whole = load_arrays(series_dir/'initial-fields.f8.gz')
    fkey = f"{name}:t{G['targets'][0]:+.0e}"
    got = fresh[('field_state', fkey)]
    checks['complete_field_state'] = ([] if all(np.array_equal(np.array(got[part]), whole[f'{fkey}|{part}']) for part in ('C', 'E', 'axis'))
                                      else ['the regenerated field differs from the saved one'])
    wrong_hash = [k for k, f in fields.items() if _sha(whole[f'{k}|C'], whole[f'{k}|E']) != f['field_sha256']]
    used = {r['field_sha256'] for k, r in runs.items() if r['target'] != 0. and r['variant'] in ('base', 'half_step', 'unprepared_bodies')
            and r['quarter'] == q}
    saved = {f['field_sha256'] for f in fields.values()}
    checks['saved_fields_are_the_fields_used'] = wrong_hash + sorted(used - saved)
    quiet_samples = S7._load_gz(series_dir/'quiet-samples.json.gz')
    checks['quiet_sample'] = evidence_io.compare(fresh[('Q', f'{name}:q{q}:k0')], quiet_samples[f'{name}:q{q}:k0'], rtol=1e-9, atol=1e-15)

    fits = {key: R10.growth_fit(r) for key, r in runs.items() if r['target'] != 0.}
    again = dict(Q1_the_quiet_sample=R10.gate_Q1(quiet_samples), Q2_the_symmetry_survives_the_live_run=R10.gate_Q2(runs),
                 P1_the_seeded_field=R10.gate_P1(runs, pops), G1_linear_in_the_seed=R10.gate_G1(runs),
                 G2_step_and_grid=R10.gate_G2(fits),
                 G3_every_run_completed=dict(statuses={k: r['status'] for k, r in runs.items()},
                                             passed=bool(all(r['status'] == 'completed' for r in runs.values()))))
    checks['fits'] = evidence_io.compare(_safe(fits), arc['fits'], rtol=1e-9, atol=1e-13)
    for k, g in again.items():
        checks[f'gate_{k}'] = evidence_io.compare(_safe(g), arc['gates'][k], rtol=1e-9, atol=1e-13)
    pred8 = json.loads(_find(R10.STAGE8).read_text(encoding='utf-8'))
    checks['diagnostics'] = evidence_io.compare(_safe(R10.diagnostics(runs, arc['background_sensitivity'])), arc['diagnostics'], rtol=1e-9, atol=1e-13)
    checks['reading'] = evidence_io.compare(_safe(R10.reading(fits, pred8, {n: r['recalculated'] for n, r in arc['reading'].items()})),
                                            arc['reading'], rtol=1e-9, atol=1e-13)
    out['recomputed'] = {k: dict(differences=v[:5], count=len(v), passed=not v) for k, v in checks.items()}

    reproduction = bool(all(v['passed'] for k, v in out.items() if k != 'recomputed')
                        and all(c['passed'] for c in out['recomputed'].values()))
    numerical = bool(arc['statuses']['numerical_verification']['passed'] and all(g['passed'] for g in again.values()))
    out['statuses'] = dict(reproduction=dict(passed=reproduction),
                           numerical_verification=dict(passed=numerical,
                                                       archive_failed_gates=arc['statuses']['numerical_verification']['failed_gates'],
                                                       rederived_here={k: bool(g['passed']) for k, g in again.items()}),
                           scientific_outcome=dict(affects_exit_status=False,
                                                   reading={n: r['reading'] for n, r in arc['reading'].items()},
                                                   measured={n: r['measured'] for n, r in arc['reading'].items()},
                                                   forecast_of_record={n: r['forecast_of_record'] for n, r in arc['reading'].items()}))
    out['deviations'] = arc['deviations']
    out['passed'] = bool(reproduction and numerical)
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
