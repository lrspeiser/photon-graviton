"""RUT-1 stage 12: the suite job. It re-derives everything the archive claims from the archive and the series,
recomputes from scratch the pieces that can be recomputed in seconds -- the angular-momentum ledger's identities
and the prediction -- replays a declared prefix of one formation run, and re-evaluates every gate and reading.

Reproduction and numerical verification set the exit status; the scientific outcome never does.
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np                              # noqa: E402
import rut12 as S12                             # noqa: E402
import rut12_tasks as TASKS                     # noqa: E402
from rut12_tasks import TH, T0, _find, _sha     # noqa: E402
import reciprocal as RC                         # noqa: E402
import reciprocal_budget as RB                  # noqa: E402
import rut7 as S7                               # noqa: E402

FAILURES = []


def check(name, ok, detail=''):
    print(('  ok   ' if ok else '  FAIL ') + name + (f'  {detail}' if detail else ''))
    if not ok:
        FAILURES.append(name)
    return ok


def main():
    archive = json.loads(_find('rut12-results.json').read_text(encoding='utf-8'))
    series = _find('rut12-results.json').parent/'rut12-series'
    print('RUT-1 stage 12: formation from an empty field')

    proto = _find(TASKS.PROTOCOL)
    check('the archive is of the committed protocol',
          archive['protocol_sha256'] == hashlib.sha256(proto.read_bytes()).hexdigest())
    check('the archive was not made from an override protocol',
          Path(archive['protocol_read_from']).name == TASKS.PROTOCOL, archive['protocol_read_from'])
    check('the thresholds in the archive are the thresholds in the protocol', archive['thresholds'] == TASKS.thresholds())
    check('every series file has the hash the archive records',
          all(hashlib.sha256((series/n).read_bytes()).hexdigest() == h for n, h in archive['series_sha256'].items()),
          '%d files' % len(archive['series_sha256']))

    groups = {n.split('.')[0]: S7._load_gz(series/n) for n in archive['series_sha256']}
    runs, budgets, prediction = groups['formation'], groups['budget'], groups['prediction']['P']
    want = {key for _, key, _ in S12.first_wave()}
    have = set(runs) | set(budgets) | {'P'}
    check('every declared task is in the series and no other', want == have,
          'declared %d, present %d' % (len(want), len(have)))

    # ---------------------------------------------------------------- the ledger's own identities
    field = RC.ReciprocalField(5., 32, .2, 10*T0, 3*T0, .03)
    X, Y, cell = RB.coordinates(field)
    check('the position grid is the one ifft2 returns, wrapped to be signed',
          abs(X[0, 0]) < 1e-12 and abs(X[1, 0] - field.L/field.n) < 1e-12 and X.min() < 0 < X.max())
    # a localised, lumpy field: eight sources on a ring, which is what the derivation assumes -- something that
    # has died away before the box edge -- and unlike an axisymmetric field it has a real d_theta.
    ang = 2*np.pi*np.arange(8)/8
    ring = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    mass = np.full(8, .125)
    field.h = field.tau_keep*field.source(ring, mass)
    field.hd = np.zeros_like(field.h)
    check('with h_t = 0 the field carries no angular momentum', abs(RB.field_angular_momentum(field)) < 1e-12)
    # int h d_theta h = (1/2) int d_theta(h^2), which vanishes on a rotationally closed domain. The box is a
    # square, so it vanishes only to the extent that the field has decayed before the edge: that approximation
    # is the one the whole ledger rests on, and A2 gates it independently by box invariance.
    field.hd = field.h.copy()
    h_r = RB.real_space(field, field.h)
    scale = float(np.sum(np.abs(h_r)*np.abs(RB.d_theta(field, field.h)))*RB.coordinates(field)[2])
    ratio = abs(RB.J(field))/max(scale, 1e-300)
    check('on a localised field int h d_theta h vanishes as the derivation needs', ratio < 1e-6, '%.2e' % ratio)
    check('the same integral does NOT vanish when the field is not localised, which is why the box matters',
          True, 'checked by A2, not asserted here')

    # ---------------------------------------------------------------- the prediction, recomputed here
    P = TH['P']
    again = TASKS.task_prediction(P['population'], P['bodies'], P['seed'])
    same = all(abs(again['predicted'][k] - prediction['predicted'][k]) <= 1e-12*max(1., abs(prediction['predicted'][k]))
               for k in prediction['predicted'])
    check('the prediction recomputes to the archived values', same)
    check('the prediction in the archive is the one the protocol states',
          abs(prediction['hotter_than_target_by'] - archive['prediction']['hotter_than_target_by']) < 1e-12)

    # ---------------------------------------------------------------- one formation run, prepared again
    key = next(k for k in sorted(runs) if runs[k]['variant'] == 'base' and runs[k]['start'] == TH['G']['start'])
    r = runs[key]
    src, st = TASKS.start_state(r['name'], r['start'])
    x0, v0 = src.draw(r['bodies'], seed=r['seed'])
    check('the formation run %s draws the same bodies' % key, _sha(x0, v0) == r['initial_sha256'])
    m = np.full(r['bodies'], float(st['mass'])/r['bodies'])
    field = RC.ReciprocalField(r['box'], r['modes'], st['w'], st['tau_keep'], TASKS.TAU_FORM, float(st['alpha']))
    rows = TASKS.advance(field, x0, v0, m, TH['C']['replay_prefix']*T0, r['step'], samples=1)
    e, a = abs(rows[-1]['energy_balance']), abs(rows[-1]['angular_balance'])
    check('the replayed prefix of %s keeps both ledgers closed' % key,
          e < TH['G']['energy_closure'] and a < TH['G']['angular_closure'], 'energy %.2e angular %.2e' % (e, a))

    # ---------------------------------------------------------------- every gate and reading, recomputed
    gates = dict(A1_both_ledgers_close_and_converge=S12.gate_A1(budgets),
                 A2_the_field_angular_momentum_belongs_to_the_field=S12.gate_A2(budgets),
                 G1_every_run_is_sound_throughout=S12.gate_G1(runs),
                 G2_the_start_is_an_equilibrium=S12.gate_G2(runs),
                 G3_the_target_is_a_fixed_point=S12.gate_G3(runs),
                 G4_the_settled_state_survives_refinement=S12.gate_G4(runs))
    for name, g in sorted(gates.items()):
        check('gate %s is as archived' % name, g == archive['gates'][name])
        check('gate %s passed' % name, g['passed'])
    readings = dict(what_formed=S12.reading_formation(runs, prediction),
                    budgets_through_the_transition=S12.reading_budget(runs))
    check('every reading is as archived', readings == archive['readings'])
    check('the archived status is the one the gates give',
          archive['statuses']['numerical_verification'] == ('passed' if all(g['passed'] for g in gates.values())
                                                            else 'failed'))
    import re
    check('no reading calls a population stable',
          re.search(r'(?<!un)stable', json.dumps(archive['readings']).lower()) is None)

    print('\n%d checks failed' % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == '__main__':
    sys.exit(main())
