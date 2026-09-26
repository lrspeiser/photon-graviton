"""Round 25: the suite's full-tier scoreboards for the trapped-companion candidates (round25-trapped-suite.md) side by
side with the adopted law's baseline (regression/baseline.json): tallies by group, every change of grade, the key
numbers, and what the trapping did (escaping shares and zones, from each run's trapping_stats.json).

    python code/trapped_suite_table_v25.py --runs run-trapped-v25/suite --output run-trapped-v25/suite_table_v25.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ORDER = ['trapped_own_heat_r12', 'trapped_cluster_heat_r12']
GROUPS = ['machinery', 'galaxies', 'clusters', 'lensing', 'milky_way', 'dwarfs', 'precision', 'collisions']
KEYS = ['law.u_kms', 'galaxies.sparc_rms_all', 'galaxies.sparc_median', 'galaxies.sparc_bulges', 'clusters.xcop_rms',
        'clusters.xcop_radial_trend', 'clusters.u_preferred', 'clusters.xcop_crossval', 'lensing.kids_all', 'lensing.kids_red',
        'lensing.kids_blue', 'lensing.kids_disc', 'lensing.kids_bulge', 'lensing.kids_gama', 'lensing.kids_gap',
        'lensing.kids_gap_sersic', 'lensing.mistele_ltg', 'lensing.mistele_etg', 'lensing.slacs_light_equals_matter',
        'mw.v_sun', 'mw.rc_outer', 'mw.vertical_pull', 'mw.mass_20', 'mw.mass_50', 'mw.mass_100', 'mw.mass_200', 'mw.v_escape',
        'mw.inner_share', 'dwarfs.chi2', 'precision.double_pulsar', 'bullet.outer_stars', 'bullet.kappa_main', 'bullet.kappa_sub',
        'bullet.gas_main', 'bullet.gas_sub', 'bullet.peak_main', 'bullet.peak_sub', 'bullet.m250_main', 'bullet.m250_sub',
        'stack.beta_gas', 'macs0025.m300_se', 'macs0025.m300_nw', 'macs0025.peak_se', 'macs0025.peak_nw', 'macs0025.sigma',
        'abell520.m710', 'abell520.sigma_peaks', 'elgordo.m500', 'elgordo.m1000', 'elgordo.sigma_nw', 'elgordo.sigma_se']


def tally(checks, groups=None):
    t = dict(p=0, c=0, f=0)
    for c in checks:
        if groups and c['group'] not in groups:
            continue
        if c['status'] in ('pass', 'close', 'fail'):
            t[c['status'][0]] += 1
    return f"{t['p']} / {t['c']} / {t['f']}"


def trapping_summary(stats):
    """The escaping shares and zones the run's last calls found."""
    if not stats:
        return {}
    out = {}
    if 'xcop' in stats:
        sh = [v['escaping_share'] for v in stats['xcop'].values()]
        out['xcop_escaping_share'] = dict(mean=float(np.mean(sh)), min=float(np.min(sh)), max=float(np.max(sh)),
                                          u_last=float(list(stats['xcop'].values())[-1]['u']))
    if 'lenses' in stats:
        out['lenses'] = {k: dict(escaping_share=v['escaping_share'], k=v['k'], size_kpc=v['size_kpc'], profile=v['profile'],
                                 zone_outer_kpc=v['zone_outer_kpc']) for k, v in stats['lenses'].items()}
    if 'sparc' in stats:
        held = [v['bulge_share_held'] for v in stats['sparc'].values()]
        out['sparc'] = dict(galaxies_with_a_zone=len(held), bulge_share_held_median=float(np.median(held)),
                            radii_in_zone_median=float(np.median([v['radii_in_zone'] for v in stats['sparc'].values()])),
                            band_median=float(np.median([v['band'] for v in stats['sparc'].values()])))
    if 'milky_way' in stats:
        out['milky_way'] = stats['milky_way']
    if 'maps' in stats:
        out['maps'] = dict(n=len(stats['maps']), escaping_share=[round(m['escaping_share'], 4) for m in stats['maps']],
                           iterations=[m['iterations'] for m in stats['maps']])
    for k in ('slacs_last', 'jeans_last'):
        if k in stats:
            out[k] = stats[k]
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--runs', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    base = json.loads((HERE.parent / 'regression/baseline.json').read_text())
    runs = {'round12 (the law)': (base, None)}
    for name in ORDER:
        d = args.runs / f'{name}-full'
        if (d / 'results.json').exists():
            st = d / 'trapping_stats.json'
            runs[name] = (json.loads((d / 'results.json').read_text()), json.loads(st.read_text()) if st.exists() else None)
    bid = {c['id']: c for c in base['checks']}
    out = dict(registration='round25-trapped-suite.md', runs={}, keys=KEYS)
    for name, (r, st) in runs.items():
        ch = r['checks']; by = {c['id']: c for c in ch}
        law = r['law']
        moved = [dict(id=c['id'], title=c['title'], value=c['value'], status=c['status'], baseline_value=bid[c['id']]['value'],
                      baseline_status=bid[c['id']]['status'])
                 for c in ch if c['id'] in bid and c['status'] != bid[c['id']]['status']]
        out['runs'][name] = dict(
            constants=dict(a_SI=law['a_SI'], u_kms=law['u_kms'], g_d_SI=law['g_d_SI'], refit_log=law.get('refit_log')),
            companion_trapping=law.get('companion_trapping', 'none'), commit=r.get('commit'), seconds=r.get('seconds'),
            total=tally(ch), by_group={g: tally(ch, {g}) for g in GROUPS}, grade_changes=moved,
            key_values={k: dict(value=by[k]['value'], status=by[k]['status'], unit=by[k].get('unit', '')) for k in KEYS if k in by},
            trapping=trapping_summary(st))
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    names = list(runs)
    print('| group | ' + ' | '.join(names) + ' |')
    print('|---|' + '---|' * len(names))
    print('| **all 77 graded checks** | ' + ' | '.join(out['runs'][n]['total'] for n in names) + ' |')
    for g in GROUPS:
        print(f'| {g} | ' + ' | '.join(out['runs'][n]['by_group'][g] for n in names) + ' |')
    print()
    for k in KEYS:
        vals = []
        for n in names:
            kv = out['runs'][n]['key_values'].get(k)
            vals.append(f"{kv['value']:.4g} ({kv['status'][0]})" if kv and kv['value'] is not None else '-')
        print(f'| {k} | ' + ' | '.join(vals) + ' |')
    for n in names[1:]:
        print(f'\n{n}: grade changes against the baseline')
        for m in out['runs'][n]['grade_changes']:
            print(f"  {m['id']}: {m['baseline_value']:.4g} ({m['baseline_status']}) -> {m['value']:.4g} ({m['status']})")
        print(json.dumps(out['runs'][n]['trapping'], indent=1, default=float)[:3000])


if __name__ == '__main__':
    main()
