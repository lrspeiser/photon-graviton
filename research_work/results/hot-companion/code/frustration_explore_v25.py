"""Round 25, after the registered test: what does organise the law's misses? Exploratory (not registered): read the
test points of code/frustration_v25.py and ask
  1. how one-sided the misses are, passes included (sign test over the graded checks);
  2. whether the misses follow the Newtonian pull where they are made (low pull -> larger miss?);
  3. for the dwarfs, whether they follow how strongly the Galaxy's pull outweighs their own (the proposal's "a larger
     outside current overwhelms a smaller one"), which the frustration measure does not see (one dominant stream
     flows cleanly, so its chi is near zero);
  4. the size of the miss where nothing is frustrated (chi < 0.01).
Any lead found here needs its own registered test on new data.

    python code/frustration_explore_v25.py run-frustration-v25/frustration_v25.json
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu, binomtest


def main():
    f = Path(sys.argv[1])
    d = json.loads(f.read_text())
    P, ck = d['points'], d['checks']['rows']
    out = {}
    # 1. one-sidedness: the law below (res > 0) or above (res < 0) each graded measurement
    below = [c['id'] for c in ck if c['res'] > 0]; above = [c['id'] for c in ck if c['res'] < 0]
    bp = [c for c in ck if c['status'] == 'pass']
    out['sign'] = dict(n=len(ck), law_below=len(below), law_above=len(above), above_ids=above,
                       p_two_sided=float(binomtest(len(below), len(ck)).pvalue),
                       passes=dict(n=len(bp), law_below=sum(c['res'] > 0 for c in bp)),
                       not_passing=dict(n=len(ck) - len(bp), law_below=sum(c['res'] > 0 for c in ck if c['status'] != 'pass')),
                       median_res_pass=float(np.median([c['res'] for c in bp])),
                       median_res_not=float(np.median([c['res'] for c in ck if c['status'] != 'pass'])))
    # 2. the miss against the Newtonian pull, check by check (mean over a check's points)
    g = {}
    for p in P:
        if np.isfinite(p['gN_SI']):
            g.setdefault(p['id'], []).append(p['gN_SI'])
    rows = [(np.log10(np.mean(g[c['id']])), c['res'], c['status'], c['family']) for c in ck if c['id'] in g]
    lg, rs = np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
    rho = spearmanr(lg, rs)
    out['vs_gN_checks'] = dict(n=len(rows), spearman=float(rho.statistic), p=float(rho.pvalue),
                               table=[dict(log10_gN=float(a), res=float(b), status=s, family=fm) for a, b, s, fm in rows])
    fam = {}
    for p in P:
        fam.setdefault(p['family'], []).append(p)
    out['vs_gN_families'] = {}
    for k, rr in fam.items():
        x = np.log10([p['gN_SI'] for p in rr]); y = np.array([p['res'] for p in rr])
        if np.isfinite(x).all() and np.ptp(x) > 0 and len(rr) > 5:
            r = spearmanr(x, y)
            out['vs_gN_families'][k] = dict(n=len(rr), spearman=float(r.statistic), p=float(r.pvalue))
    # 3. the dwarfs: the Galaxy's pull over the dwarf's own at the half-light radius
    dw = [p for p in P if p['family'] == 'dwarfs']
    eta = np.array([p['ext_over_int_at_rh'] for p in dw]); res = np.array([p['res'] for p in dw])
    npass = np.array([p['status'] != 'pass' for p in dw])
    r = spearmanr(eta, res)
    u = mannwhitneyu(eta[npass], eta[~npass], alternative='greater')
    out['dwarfs_outside_pull'] = dict(spearman=float(r.statistic), p=float(r.pvalue), mannwhitney_p=float(u.pvalue),
                                      table=sorted([dict(name=p['id'].split('.')[1], ext_over_int=p['ext_over_int_at_rh'], res=p['res'],
                                                         status=p['status'], chi=p['chi']) for p in dw], key=lambda t: t['ext_over_int']))
    # 4. the miss where nothing is frustrated
    out['unfrustrated'] = {k: dict(n=int(np.sum([p['chi'] < 0.01 for p in rr])),
                                   median_res=float(np.median([p['res'] for p in rr if p['chi'] < 0.01])))
                           for k, rr in fam.items() if any(p['chi'] < 0.01 for p in rr)}
    (f.parent / 'explore_v25.json').write_text(json.dumps(out, indent=1) + '\n')
    s = out['sign']
    print(f"graded checks with a pull-like measurement: {s['n']}; law below {s['law_below']}, above {s['law_above']} "
          f"(p = {s['p_two_sided']:.2g}); passes {s['passes']['law_below']}/{s['passes']['n']} below; "
          f"not passing {s['not_passing']['law_below']}/{s['not_passing']['n']} below")
    print('  law above:', ', '.join(s['above_ids']))
    print(f"miss vs Newtonian pull (checks): rho = {out['vs_gN_checks']['spearman']:+.2f} (p = {out['vs_gN_checks']['p']:.2g})")
    for k, v in out['vs_gN_families'].items():
        print(f"  {k:10s} n = {v['n']:5d}  rho = {v['spearman']:+.3f} (p = {v['p']:.2g})")
    dd = out['dwarfs_outside_pull']
    print(f"dwarfs: miss vs outside/own pull rho = {dd['spearman']:+.2f} (p = {dd['p']:.2g}); non-passes have the larger ratio: p = {dd['mannwhitney_p']:.3g}")
    for t in dd['table']:
        print(f"  {t['name']:12s} outside/own {t['ext_over_int']:7.2f}  miss {t['res']:+.2f}  {t['status']:5s}  chi {t['chi']:.3f}")
    for k, v in out['unfrustrated'].items():
        print(f"unfrustrated (chi < 0.01) {k:10s} n = {v['n']:5d}  median miss {v['median_res']:+.3f}  (x{np.exp(v['median_res']):.2f})")


if __name__ == '__main__':
    main()
