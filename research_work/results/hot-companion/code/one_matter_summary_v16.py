"""Summaries of code/one_matter_v16.py's runs: per case, the mean over arrangements (and the spread), with the sources'
output and the receivers' pull also given relative to the same arrangement at rest.

    python code/one_matter_summary_v16.py run-one-matter-v16/one_matter_v16.json [more.json ...]
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from one_matter_v16 import k_nominal                 # noqa: E402

KEYS = ('src_out', 'src_bright', 'pull', 'pull_s', 'pull_s_src', 'pull_B', 'lead', 'lead_src', 'fed_src', 'E_src', 'amp_s',
        'w_s', 'sync', 'rhythm_sources', 'rhythm_receivers', 'rhythm_spread_sources', 'rhythm_spread_receivers')


def summarize(path, baseline=None):
    d = json.loads(Path(path).read_text())
    runs = d['runs']
    groups = {}
    for r in runs:
        c = r['cfg']
        groups.setdefault((c['tag'], round(c.get('q', 0.0), 4), c.get('nu', 0.0), c.get('Delta0', 0.0)), []).append(r)
    cold = {r['cfg']['seed']: r for r in groups.get(('cold', 0.0, 0.0, 0.0), [])}
    if not cold and baseline is not None:                  # e.g. the detuned set, measured against the main run's rest
        cold = {r['cfg']['seed']: r for r in json.loads(Path(baseline).read_text())['runs'] if r['cfg']['tag'] == 'cold'}
    rows = []
    order = sorted(groups, key=lambda k: ({'cold': 0, 'cold_detuned': 1, 'free': 2, 'collisional': 3}[k[0]], k[1], k[2], k[3]))
    for key in order:
        rs = groups[key]
        row = dict(tag=key[0], q=key[1], nu=key[2], Delta0=key[3], k=k_nominal(rs[0]['cfg']), arrangements=len(rs))
        for kk in KEYS:
            v = np.array([r.get(kk, np.nan) for r in rs], float)
            row[kk] = float(np.nanmean(v)); row[kk + '_sd'] = float(np.nanstd(v, ddof=1)) if len(v) > 1 else None
        rel_out = [r['src_out'] / cold[r['cfg']['seed']]['src_out'] for r in rs if r['cfg']['seed'] in cold]
        rel_pull = [r['pull_s_src'] / cold[r['cfg']['seed']]['pull_s_src'] for r in rs if r['cfg']['seed'] in cold]
        row['output_rel_rest'] = float(np.mean(rel_out)) if rel_out else None
        row['E_src_rel_rest'] = float(np.mean([r['E_src'] / cold[r['cfg']['seed']]['E_src'] for r in rs if r['cfg']['seed'] in cold])) if cold else None
        row['quiet_pull_from_sources_rel_rest'] = float(np.mean(rel_pull)) if rel_pull else None
        row['sqrt_output_rel_rest'] = float(np.sqrt(np.mean(rel_out))) if rel_out else None
        row['energy_residual_max'] = float(max(abs(r['energy']['residual']) for r in rs))
        row['momentum_balance_max'] = float(max(m['balance'] for r in rs for m in r['momentum']))
        rows.append(row)
    return rows


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--baseline=')]
    base = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--baseline=')), None)
    for p in args:
        rows = summarize(p, base)
        print(p)
        print(f"{'case':>12} {'q':>5} {'nu':>5} {'k':>5} | {'output':>9} {'×rest':>6} {'√':>5} | {'net pull':>10} {'quiet':>10} "
              f"{'from src':>10} {'radiators':>10} | {'lead':>5} {'l_src':>5} | {'fed_src':>9} | {'src rhythm spread':>17} | {'E res':>7}")
        for r in rows:
            print(f"{r['tag']:>12} {r['q']:5.2f} {r['nu']:5.0f} {r['k']:5.2f} | {r['src_out']:9.3e} {r['output_rel_rest'] or 1:6.2f} "
                  f"{r['sqrt_output_rel_rest'] or 1:5.2f} | {r['pull']:+10.3e} {r['pull_s']:+10.3e} {r['pull_s_src']:+10.3e} {r['pull_B']:+10.3e} | "
                  f"{r['lead']:5.2f} {r['lead_src']:5.2f} | {r['fed_src']:9.3e} | {r['rhythm_spread_sources']:17.2e} | {r['energy_residual_max']:7.1e}")
        out = Path(p).with_name(Path(p).stem + '_summary.json')
        out.write_text(json.dumps(rows, indent=1) + '\n')
        print('wrote', out)


if __name__ == '__main__':
    main()
