"""Round 24: the suite's full-tier scoreboards for the Casimir-EFT candidates side by side with the adopted law's
baseline (regression/baseline.json), by group, with the key numbers each group is graded on.

    python code/eft_suite_table_v24.py --runs run-eft-v24/suite --output run-eft-v24/suite_table_v24.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
KEYS = ['galaxies.sparc_rms_all', 'galaxies.sparc_bulges', 'clusters.xcop_rms', 'clusters.xcop_radial_trend',
        'clusters.u_preferred', 'lensing.kids_red', 'lensing.kids_gap', 'lensing.mistele_etg', 'lensing.slacs_light_equals_matter',
        'mw.v_sun', 'precision.planets', 'precision.cassini_q2', 'bullet.outer_stars', 'bullet.kappa_sub', 'bullet.gas_main',
        'bullet.peak_main', 'bullet.peak_sub', 'bullet.m250_main', 'bullet.m250_sub', 'stack.beta_gas', 'macs0025.peak_nw',
        'macs0025.sigma', 'abell520.m710', 'abell520.sigma_peaks', 'elgordo.m500', 'elgordo.sigma_nw']
ORDER = ['eft_r12', 'eft_r12_refit', 'eft_nohold_r12', 'eft_memory_r12', 'eft_excess_local_r12', 'eft_excess_r12',
         'eft_scalar_local_r12', 'eft_scalar_r12']


def tally(checks, groups=None):
    t = dict(pass_=0, close=0, fail=0)
    for c in checks:
        if groups and c['group'] not in groups:
            continue
        if c['status'] in ('pass', 'close', 'fail'):
            t['pass_' if c['status'] == 'pass' else c['status']] += 1
    return f"{t['pass_']} / {t['close']} / {t['fail']}"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--runs', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    base = json.loads((HERE.parent / 'regression/baseline.json').read_text())
    runs = {'round12 (the law)': base}
    for name in ORDER:
        f = args.runs / f'{name}-full' / 'results.json'
        if f.exists():
            runs[name] = json.loads(f.read_text())
    groups = ['machinery', 'galaxies', 'clusters', 'lensing', 'milky_way', 'dwarfs', 'precision', 'collisions']
    out = dict(runs={}, keys=KEYS)
    for name, r in runs.items():
        ch = r['checks']; by = {c['id']: c for c in ch}
        law = r['law']
        out['runs'][name] = dict(
            constants=dict(a_SI=law['a_SI'], u_kms=law['u_kms'], g_d_SI=law['g_d_SI']),
            hot_geometry=law.get('hot_geometry', 'two_way'), collision_field=law.get('collision_field', 'law'),
            total=tally(ch), by_group={g: tally(ch, {g}) for g in groups},
            regressed=sorted(c['id'] for c in ch if c.get('move') == 'regressed'),
            improved=sorted(c['id'] for c in ch if c.get('move') == 'improved'),
            key_values={k: dict(value=by[k]['value'], status=by[k]['status']) for k in KEYS if k in by})
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    names = list(runs)
    print('| group | ' + ' | '.join(names) + ' |')
    print('|---|' + '---|' * len(names))
    print('| **all 77 graded checks** | ' + ' | '.join(out['runs'][n]['total'] for n in names) + ' |')
    for g in groups:
        print(f'| {g} | ' + ' | '.join(out['runs'][n]['by_group'][g] for n in names) + ' |')
    print()
    for k in KEYS:
        vals = []
        for n in names:
            kv = out['runs'][n]['key_values'].get(k)
            vals.append(f"{kv['value']:.3g} ({kv['status'][0]})" if kv else '-')
        print(f'| {k} | ' + ' | '.join(vals) + ' |')


if __name__ == '__main__':
    main()
