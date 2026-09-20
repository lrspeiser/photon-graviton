"""RUT-1 stage 11: the suite job. It re-derives everything the archive claims from the archive and the series,
and it recomputes from scratch the pieces that can be recomputed in seconds — the analytic rule's closed forms
against a reference integral, the disturbance's exact preparation, and every gate, reading and fit.

Reproduction and numerical verification set the exit status; the scientific outcome never does.
"""
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np                                  # noqa: E402
import rut11 as S11                                 # noqa: E402
import rut11_tasks as TASKS                         # noqa: E402
from rut11_tasks import REP, TH, _find, _sha        # noqa: E402
import resonant_response as RR                      # noqa: E402
import rut7 as S7                                   # noqa: E402

FAILURES = []


def check(name, ok, detail=''):
    print(('  ok   ' if ok else '  FAIL ') + name + (f'  {detail}' if detail else ''))
    if not ok:
        FAILURES.append(name)
    return ok


def main():
    archive = json.loads(_find('rut11-results.json').read_text(encoding='utf-8'))
    series = _find('rut11-results.json').parent/'rut11-series'
    print('RUT-1 stage 11: the quiet region')

    # ---------------------------------------------------------------- the archive is of this protocol and code
    proto = _find(TASKS.PROTOCOL)
    check('the archive is of the committed protocol',
          archive['protocol_sha256'] == hashlib.sha256(proto.read_bytes()).hexdigest())
    check('the archive was not made from an override protocol',
          Path(archive['protocol_read_from']).name == TASKS.PROTOCOL,
          archive['protocol_read_from'])
    check('the thresholds in the archive are the thresholds in the protocol',
          archive['thresholds'] == TASKS.thresholds())
    check('every series file has the hash the archive records',
          all(hashlib.sha256((series/n).read_bytes()).hexdigest() == h
              for n, h in archive['series_sha256'].items()),
          '%d files' % len(archive['series_sha256']))

    groups = {n.split('.')[0]: S7._load_gz(series/n) for n in archive['series_sha256']}
    boxes = {k: v for k, v in groups['boxes'].items() if k.startswith('box:')}
    old_rule = {k: v for k, v in groups['boxes'].items() if k.startswith('oldrule:')}
    for v in old_rule.values():
        v['matches'] = S11.box_key(v['name'], v['m'], v['shifted'], v['index'])
    norms, ladder, runs = groups['norms'], groups['threshold'], groups['disturbance']
    a2 = {k: v for k, v in groups['rule'].items() if k.startswith('A2:')}

    # ---------------------------------------------------------------- the campaign is the declared one
    want1 = {key for _, key, _ in S11.first_wave()}
    want2 = {key for _, key, _ in S11.second_wave()}
    have = set(boxes) | set(old_rule) | set(norms) | set(ladder) | set(runs) | set(groups['rule'])
    check('every declared task is in the series and no other', want1 | want2 == have,
          'declared %d, present %d' % (len(want1 | want2), len(have)))

    # ---------------------------------------------------------------- the rule, recomputed here
    A = TH['A']['toy']
    worst_pair, worst_gauss = 0., 0.
    for row in groups['rule']['A1']['rows']:
        z = complex(A['omega'], row['growth_rate'])
        ref = RR.reference_integral(z, A['a'], A['b'], A['c'], A['dL'], A['dE'], A['reach'])
        check('the reference integral at growth rate %g is the archived one' % row['growth_rate'],
              abs(complex(*row['reference']) - ref) <= 1e-12*abs(ref))
        vals = {k: RR.reference_by_rule(z, A['a'], A['b'], A['c'], A['dL'], A['dE'], A['reach'],
                                        A['n_L'], A['n_u'], k) for k in REP['refines']}
        rich = vals[REP['refines'][1]] + (vals[REP['refines'][1]] - vals[REP['refines'][0]])/3.
        worst_pair = max(worst_pair, abs(rich - ref)/abs(ref))
        if row['growth_rate'] <= A['control_below']:
            g = RR.reference_by_rule(z, A['a'], A['b'], A['c'], A['dL'], A['dE'], A['reach'], A['n_L'], A['n_u'],
                                     1, gauss=True)
            worst_gauss = max(worst_gauss, abs(g - ref)/abs(ref))
    check('recomputed here, the rule holds its declared accuracy', worst_pair < A['accuracy'], '%.2e' % worst_pair)
    check("recomputed here, stage 9's quadrature fails the same integral", worst_gauss > A['control_min'],
          '%.2e' % worst_gauss)

    # the closed forms where they do not cancel, and the two branches where they meet. The identities
    # themselves lose digits for small x -- x - log(1+x) is x^2/2 computed as a difference of two numbers of
    # size x -- so they are only checked where the right-hand side is well conditioned.
    big = np.array([.3, 1., 7., 30.])*np.exp(.7j)
    b0, b1, b2 = RR.logs(big)
    check('L0 satisfies its identity where it does not cancel', np.allclose(b0*big, np.log1p(big), rtol=1e-14, atol=0))
    check('L1 satisfies its identity where it does not cancel',
          np.allclose(b1*big*big, big - np.log1p(big), rtol=1e-13, atol=0))
    check('L2 satisfies its identity where it does not cancel',
          np.allclose(b2*big**3, np.log1p(big) - big + .5*big*big, rtol=1e-12, atol=0))
    eps = 1e-9
    lo = RR.logs(np.array([RR._SMALL*(1 - eps)])*np.exp(.7j))            # the series branch
    hi = RR.logs(np.array([RR._SMALL*(1 + eps)])*np.exp(.7j))            # the closed-form branch
    worst = max(float(abs(a - b)/abs(b)) for a, b in zip(lo, hi))
    # the seam is set by the closed form, not the series: at |x| = 0.01 the numerator of L2 is x^3/3 computed
    # as a difference of terms of size x, which leaves about 10^-11 of it. The rule's own target is 10^-8.
    check('the series and the closed form agree where the code switches between them', worst < 1e-9,
          '%.2e' % worst)
    # and the rule integrates a constant exactly: sum of the weights = area / (nu - z) when nu is constant
    xa, xb = np.linspace(0., 1., 17), np.linspace(0., 2., 9)
    flat = RR.Triangles(xa, xb, np.full((17, 9), .7))
    z = .3 + .002j
    check('the rule integrates a cell of constant frequency exactly',
          abs(complex(flat.weights(z).sum()) - 2./(.7 - z)) < 1e-14)

    # ---------------------------------------------------------------- one disturbance, prepared again
    D = TH['D']
    key = next(k for k in sorted(runs) if runs[k]['variant'] == 'base' and runs[k]['target'] != 0.)
    r = runs[key]
    run, info, (x2, v2) = TASKS.build_disturbance(r['name'], r['m'], r['quarter'], r['realization'], r['target'],
                                                  r['variant'])
    check('the disturbed run %s draws the same bodies' % key, info['initial_sha256'] == r['initial_sha256'])
    check('the disturbed run %s builds the same field' % key, info['field_sha256'] == r['field_sha256'])
    check('the disturbed run %s has the same amplitude' % key, abs(info['eps']/r['eps'] - 1.) < 1e-12)
    pre = run.advance(TH['C']['replay_prefix']*TASKS.T0).result()
    for a, b in zip(pre['rows'], r['rows'][:len(pre['rows'])]):
        if abs(a['t_periods'] - b['t_periods']) > 1e-12 or a['coefficients']['C'] != b['coefficients']['C']:
            check('the replayed prefix of %s matches the series record for record' % key, False,
                  't = %.3f' % a['t_periods'])
            break
    else:
        check('the replayed prefix of %s matches the series record for record' % key, True,
              '%d records' % len(pre['rows']))

    # ---------------------------------------------------------------- every gate and reading, recomputed
    quiet_and_control = TH['L']['quiet'] + [TH['L']['control_population']]
    gates = dict(A1_the_rule_on_a_known_integral=S11.gate_A1(groups['rule']['A1']),
                 A2_the_rule_agrees_with_stage9_where_stage9_is_valid=S11.gate_A2(a2, old_rule, boxes),
                 A3_the_declared_error_on_every_contour=S11.gate_A3(boxes, ladder),
                 L1_the_search_is_sound=S11.gate_L1(boxes, quiet_and_control),
                 L2_the_search_can_see_a_root_below_the_old_floor=S11.gate_L2(boxes),
                 M1_the_operator_norm_grid_is_adequate=S11.gate_M1(norms),
                 T1_every_rung_survives_refinement=S11.gate_T1(ladder),
                 D1_the_disturbance_is_as_declared=S11.gate_D1(runs),
                 D2_the_apparatus_would_see_growth=S11.gate_D2(runs),
                 D3_every_run_completed=S11.gate_D3(runs))
    for name, g in sorted(gates.items()):
        check('gate %s is as archived' % name, g == archive['gates'][name])
        check('gate %s passed' % name, g['passed'])
    readings = dict(quiet_region=S11.reading_quiet(boxes, norms, runs),
                    threshold_bracket=S11.reading_threshold(ladder),
                    disturbance_variants=S11.reading_variants(runs))
    check('every reading is as archived', readings == archive['readings'])
    check('the archived status is the one the gates give',
          archive['statuses']['numerical_verification'] == ('passed' if all(g['passed'] for g in gates.values())
                                                            else 'failed'))

    # ---------------------------------------------------------------- what the readings may not say
    text = json.dumps(archive['readings']).lower()
    check('no reading calls a population stable', re.search(r'(?<!un)stable', text) is None)

    print('\n%d checks failed' % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
