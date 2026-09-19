"""Suite job for RUT-1 stage 9 (protocol-rut9.md). rut9.py is the campaign driver; this job checks in a few
minutes that its archive is what the committed code produces, and reports the three statuses separately.

REPRODUCTION
- the committed code, the protocol, stage 8's predictions file and every series file hash to the values
  recorded at launch; the archive was made from the committed protocol and not from a toy one; the thresholds
  in the archive are the ones the protocol declares today;
- recomputed from scratch and compared with the archive: gate B1 on one population; the widest ring-limit
  annulus, built, with its exact-ring target, and its root under the base representation and under stage 8's
  node rule; one member of the family scan; the sensitivity of B13's root; and the small declared replay batch
  of the time-domain check, body for body;
- recomputed from the archive: every gate, re-derived from the archived task results and series, and every
  reading of part F.

NUMERICAL VERIFICATION
- every gate in the archive passed and rejected its negative control, and so did every gate recomputed here.

SCIENTIFIC OUTCOME
- what the corrected calculation says about stage 8's forecasts and about B13's shortfall: reported, never part
  of the exit status.

Stage 8's own job, rut8_checks.py, is not touched by any of this: its L3 and L3b stay failed and it stays red.
"""
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rut9 as R9  # noqa: E402
from rut9_tasks import PROTOCOL, TH, _find, _safe, archived_populations, run_task, thresholds  # noqa: E402
import rut7 as S7  # noqa: E402
from rut7_tasks import evidence_io  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    t0 = time.time()
    arc = json.loads((HERE/'rut9-results.json').read_text(encoding='utf-8'))
    series_dir = HERE/'rut9-series'
    out = {}
    mism = [f for f, h in arc['code_sha256_at_launch'].items() if _sha(_find(f)) != h]
    mism += [f for f, h in arc['series_sha256'].items() if _sha(series_dir/f) != h]
    mism += [f for f, h in arc['input_sha256'].items() if _sha(_find(f)) != h]
    out['code_inputs_and_series_are_what_ran'] = dict(mismatched=mism, passed=not mism)
    out['made_from_the_committed_protocol'] = dict(
        read_from=Path(arc['protocol_read_from']).name,
        passed=bool(Path(arc['protocol_read_from']).name == PROTOCOL and arc['input_sha256'][PROTOCOL] == _sha(HERE/PROTOCOL)))
    out['thresholds_are_the_protocols'] = dict(passed=bool(thresholds() == arc['thresholds'] == TH))

    pred = R9.stage8()
    pops = archived_populations()
    V = S7._load_gz(series_dir/'variants.json.gz')
    strips = {f"{r['name']}:m{r['m']}:{r['config']}:{'shifted' if r['shifted'] else 'base'}": r
              for r in S7._load_gz(series_dir/'strips.json.gz')}
    batches = S7._load_gz(series_dir/'time-domain.json.gz')
    name_B1 = TH['B']['populations'][0]
    name_E = TH['E']['population']
    resolved = [r['dE'] for r in arc['readings']['part_F']['family_B'] if r['recalculated'] is not None]
    family_dE = max(resolved) if resolved else TH['F']['family_B_dE'][0]      # the slowest root resolved
    jobs = [('B1', name_B1, (pops[name_B1], R9.stage8_root(pred, name_B1))),
            ('annulus', '0:base', (0, False)),
            ('family', f'{family_dE:.3f}', (family_dE,)),
            ('S', name_E, (pops[name_E], R9.stage8_root(pred, name_E))),
            ('E', f'{name_E}:b0:replay', (pops[name_E], R9.stage8_root(pred, name_E), 0, 'replay'))]
    fresh = {}
    with ProcessPoolExecutor(max_workers=5) as pool:
        for kind, key, result, _ in pool.map(run_task, *zip(*jobs)):
            fresh[(kind, key)] = result
    annulus = fresh[('annulus', '0:base')]
    jobs2 = [('variant', 'annulus0:base', (annulus['state'], 'annulus0', 'base', TH['N']['ring_limit_rect'], None)),
             ('variant', 'annulus0:stage8_rule', (annulus['state'], 'annulus0', 'stage8_rule', TH['N']['ring_limit_rect'], None))]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for kind, key, result, _ in pool.map(run_task, *zip(*jobs2)):
            fresh[(kind, key)] = result

    drop = lambda d: {k: v for k, v in d.items() if k not in ('seconds', 'state')}
    G = arc['gates']
    checks = {}
    checks['B1'] = evidence_io.compare(fresh[('B1', name_B1)],
                                       {k: v for k, v in G['B1_reproduces_stage8_without_the_boundary_rows']['rows'][name_B1].items()
                                        if k not in ('passed', 'control_rejected')}, rtol=1e-6, atol=1e-12)
    checks['annulus'] = evidence_io.compare(drop(annulus), drop(arc['annuli']['0:base']), rtol=1e-7, atol=1e-10)
    for key in ('annulus0:base', 'annulus0:stage8_rule'):
        checks[f'variant_{key}'] = evidence_io.compare(drop(fresh[('variant', key)]), drop(V[key]), rtol=1e-6, atol=1e-9)
    stored = [r for r in arc['readings']['part_F']['family_B'] if abs(r['dE'] - family_dE) < 1e-12][0]
    got = fresh[('family', f'{family_dE:.3f}')]
    checks['family_member'] = evidence_io.compare(dict(recalculated=got['root'], without=got['root_without_boundary_rows']),
                                                  dict(recalculated=stored['recalculated'], without=stored['without_boundary_rows']),
                                                  rtol=1e-7, atol=1e-9)
    checks['sensitivity'] = evidence_io.compare(fresh[('S', name_E)], arc['readings']['part_S'][name_E], rtol=1e-5, atol=1e-10)
    replay = [b for b in batches if b['variant'] == 'replay'][0]
    checks['time_domain_replay'] = evidence_io.compare(drop(fresh[('E', f'{name_E}:b0:replay')]), drop(replay), rtol=1e-7, atol=1e-10)

    # every gate re-derived from the archived task results
    n1 = R9.gate_N1(V, arc['annuli'])
    again = dict(B1_reproduces_stage8_without_the_boundary_rows=R9.gate_B1(
                     {k: {a: b for a, b in v.items() if a not in ('passed', 'control_rejected')}
                      for k, v in G['B1_reproduces_stage8_without_the_boundary_rows']['rows'].items()}),
                 B2_boundary_rows_are_the_limit_of_a_resolved_taper=R9.gate_B2(
                     {k: {a: b for a, b in v.items() if a not in ('passed', 'control_rejected', 'falling')}
                      for k, v in G['B2_boundary_rows_are_the_limit_of_a_resolved_taper']['rows'].items()}),
                 N1_every_annulus_resolved_on_its_own=n1,
                 N2_ring_limit=R9.gate_N2(V, arc['annuli'], n1['passed']),
                 N3_every_population_root_resolved=R9.gate_N3(V),
                 N4_root_free_conclusions_survive_refinement=R9.gate_N4(strips),
                 S1_sensitivity_formula=R9.gate_S1(arc['readings']['part_S']),
                 E1_eigenfunction_in_the_time_domain=R9.gate_E1(arc['readings']['part_E']))
    for k, g in again.items():
        checks[f'gate_{k}'] = evidence_io.compare(_safe(g), G[k], rtol=1e-9, atol=1e-12)
    out['recomputed'] = {k: dict(differences=v[:5], count=len(v), passed=not v) for k, v in checks.items()}

    reproduction = bool(all(v['passed'] for k, v in out.items() if k != 'recomputed')
                        and all(c['passed'] for c in out['recomputed'].values()))
    numerical = bool(arc['statuses']['numerical_verification']['passed'] and all(g['passed'] for g in again.values()))
    F, E = arc['readings']['part_F'], G['E1_eigenfunction_in_the_time_domain']
    out['statuses'] = dict(
        reproduction=dict(passed=reproduction),
        numerical_verification=dict(passed=numerical, archive_failed_gates=arc['statuses']['numerical_verification']['failed_gates'],
                                    rederived_here={k: bool(g['passed']) for k, g in again.items()}),
        scientific_outcome=dict(affects_exit_status=False,
                                boundary_rows_growth_relative={n: r['growth_relative'] for n, r in F['boundary_rows'].items()},
                                ring_limit_differences=[r['growth_relative_difference'] for r in G['N2_ring_limit']['rows']],
                                zero_growth_crossing_extrapolated=F['zero_growth_crossing']['recalculated'],
                                implied_growth_shift_relative=E['implied_growth_shift_relative'],
                                implied_growth_shift_standard_error=E['standard_error_relative']))
    out['deviations'] = arc['deviations']
    out['passed'] = bool(reproduction and numerical)
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
