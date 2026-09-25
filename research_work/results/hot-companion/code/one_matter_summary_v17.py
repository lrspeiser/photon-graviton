"""Summaries of code/one_matter_v17.py's runs.

Per case (structure, one-way or not, motion, heat weight, collision rate, source size): the mean over arrangements of
the sources' output (also relative to the same arrangement at rest), the wave reaching the receivers, the pulls
(net, quiet channel, from the sources' wave alone, radiators), keeping step, the sources' rhythm spread and the energy
books. For runs with receivers on several shells, per shell: the receivers' net pull F, their keeping step, the
sources' wave intensity I at them, and the review's ratio
    R(sigma, r) = (F_extra(sigma, r)/F_extra(0, r)) * sqrt(I(0, r)/I(sigma, r))
with F_extra the pull on the receivers' quiet channel from the sources' wave (the part the heat should raise), and the
same with the net pull.

    python code/one_matter_summary_v17.py run-one-matter-v17/parity.json [more.json ...]
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

KEYS = ('src_out', 'pull', 'pull_s', 'pull_s_src', 'pull_B', 'lead_src', 'E_src', 'fed_src', 'rhythm_spread_sources',
        'rhythm_sources', 'rhythm_receivers', 'amp_s', 'w_s', 'src_fam', 'src_chamber')


def medium(c):
    # round 18: runs with the absorbing stream are grouped by its kappa; otherwise by the one-way leak (None: two-way)
    return ('absorb', c['absorb']) if c.get('absorb') is not None else c.get('one_way')


def key(c):
    return (c['structure'], medium(c), c['tag'], round(c.get('k', 0.0), 3), c.get('nu', 0.0), c['Ns'])


def summarize(path):
    d = json.loads(Path(path).read_text())
    groups = {}
    for r in d['runs']:
        groups.setdefault(key(r['cfg']), []).append(r)
    cold = {(k[0], k[1], k[5], r['cfg']['seed']): r for k, rs in groups.items() if k[2] == 'cold' for r in rs}
    order = sorted(groups, key=lambda k: (k[0], str(k[1]), k[5], {'cold': 0, 'free': 1, 'collisional': 2, 'stopped': 3}[k[2]], k[3], k[4]))
    rows = []
    for kk in order:
        rs = groups[kk]
        ab = isinstance(kk[1], tuple)
        row = dict(structure=kk[0], one_way=None if ab else kk[1], tag=kk[2], k=kk[3], nu=kk[4], Ns=kk[5], arrangements=len(rs))
        if ab:
            row['absorb'] = kk[1][1]
            bk = [r.get('absorb_booking') for r in rs if r.get('absorb_booking')]
            if bk:
                row['absorbed_fraction'] = float(np.mean([b['absorbed_fraction'] for b in bk]))
                row['passivity_min'] = float(min(b['passivity_min'] for b in bk))
        for name in KEYS:
            v = np.array([r.get(name, np.nan) for r in rs], float)
            row[name] = float(np.nanmean(v)); row[name + '_sd'] = float(np.nanstd(v, ddof=1)) if len(v) > 1 else None
        base = [cold.get((kk[0], kk[1], kk[5], r['cfg']['seed'])) for r in rs]
        if all(b is not None for b in base):
            row['output_rel_rest'] = float(np.mean([r['src_out'] / b['src_out'] for r, b in zip(rs, base)]))
            row['wave_rel_rest'] = float(np.mean([r['E_src'] / b['E_src'] for r, b in zip(rs, base)]))
            row['pull_from_sources_rel_rest'] = float(np.mean([r['pull_s_src'] / b['pull_s_src'] for r, b in zip(rs, base)]))
        row['energy_residual_max'] = float(max(abs(r['energy']['residual']) for r in rs))
        if 'lead_src_each' in rs[0] and 'receiver_distance' in rs[0]:
            dist = np.round(np.array(rs[0]['receiver_distance']), 3)
            shells = sorted(set(dist.tolist()))
            per = []
            for sh in shells:
                m = dist == sh
                F = np.mean([np.array(r['pull_each'])[m].mean() for r in rs])
                Fq = np.mean([np.array(r['pull_s_src_each'])[m].mean() for r in rs])
                L = np.mean([np.array(r['lead_src_each'])[m].mean() for r in rs])
                I = np.mean([np.array(r['intensity_src_each'])[m].mean() for r in rs])
                ent = dict(r=sh, pull=float(F), pull_quiet_from_sources=float(Fq), keeping_step=float(L), intensity=float(I))
                if all(b is not None and 'lead_src_each' in b for b in base):
                    F0q = np.mean([np.array(b['pull_s_src_each'])[m].mean() for b in base])
                    F0 = np.mean([np.array(b['pull_each'])[m].mean() for b in base])
                    I0 = np.mean([np.array(b['intensity_src_each'])[m].mean() for b in base])
                    ent['R_quiet'] = float((Fq / F0q) * np.sqrt(I0 / I)) if F0q != 0 else None
                    ent['R_net'] = float((F / F0) * np.sqrt(I0 / I)) if F0 != 0 else None
                    ent['intensity_rel_rest'] = float(I / I0)
                per.append(ent)
            row['shells'] = per
        rows.append(row)
    return rows


def main():
    for p in sys.argv[1:]:
        rows = summarize(p)
        print(p)
        for r in rows:
            med = f"absorb {r['absorb']:g}" if r.get('absorb') is not None else ('one-way' if r['one_way'] is not None else 'two-way')
            head = (f"{r['structure']:>36} {med:>10} Ns {r['Ns']:3d} {r['tag']:>11} "
                    f"k {r['k']:5.1f} nu {r['nu']:4.0f} | out {r['src_out']:.3e} (x{r.get('output_rel_rest', 1):.2f}) wave x{r.get('wave_rel_rest', 1):.2f} | "
                    f"pull {r['pull']:+.2e} quiet {r['pull_s']:+.2e} from src {r['pull_s_src']:+.2e} (x{r.get('pull_from_sources_rel_rest', 1):.2f}) "
                    f"rad {r['pull_B']:+.2e} | step {r['lead_src']:+.2f} | spread {r['rhythm_spread_sources']:.1e} | E res {r['energy_residual_max']:.0e}")
            print(head)
            for sh in r.get('shells', []):
                print(f"{'':>60} r {sh['r']:5.1f}: pull {sh['pull']:+.2e}, quiet from sources {sh['pull_quiet_from_sources']:+.2e}, "
                      f"step {sh['keeping_step']:+.2f}, I x{sh.get('intensity_rel_rest', 1):.2f}, R quiet {sh.get('R_quiet') or 0:+.2f}, R net {sh.get('R_net') or 0:+.2f}")
        out = Path(p).with_name(Path(p).stem + '_summary.json')
        out.write_text(json.dumps(rows, indent=1) + '\n')
        print('wrote', out)


if __name__ == '__main__':
    main()
