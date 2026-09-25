"""Round 18: why a warm model source glows about x4.75 of its cold output at heat weight k = 8, not the law's 1 + k = 9.

k is defined for an isolated piece in linear response (code/rhythm_protect_v17.py: the mixing q is chosen so that one
piece's released glow is k times its cold glow, at fixed quiet amplitude). In the full model the quiet amplitude is not
fixed: each piece is an inverted, self-sustained oscillator whose steady amplitude falls as its total loss rises, and
in a cloud the pieces also share one wave. This separates the two: the glow of one isolated piece (Ns = 1, one far
receiver) and of 48 pieces in a dense ball (radius 3, as in rounds 16-17) and in a dilute one (radius 30, spacing many
wavelengths), at rest and at k = 2, 8 and 16, for round 16's matter (f = 0.2) and for the f = 0.05 matter.
"""
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from one_matter_v17 import case, simulate, shells       # noqa: E402

LOWF = dict(piece=dict(G=40.0, W0=16.0, g_par=0.01), dt=0.01)


def configs():
    runs = []
    for stn, extra in (('single', {}), ('single, odd, f = 0.05', LOWF)):
        for label, Ns, Rb, rec in (('isolated piece', 1, 0.0, shells((30.0,), per=1)),
                                   ('dense ball (radius 3)', 48, 3.0, shells((12.0,), per=2)),
                                   ('dilute ball (radius 30)', 48, 30.0, shells((60.0,), per=2))):
            for k in (0.0, 2.0, 8.0, 16.0):
                kw = dict(Ns=Ns, Rb=Rb, receivers=rec, T=4000.0, burn=2000.0, momentum=False, seed=1, **extra)
                c = case(stn, 'cold', **kw) if k == 0 else case(stn, 'free', k=k, **kw)
                c['label'] = label
                runs.append(c)
    return runs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=3)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic(); runs = configs(); res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            c = r['cfg']
            res.append(dict(structure=c['structure'], label=c['label'], Ns=c['Ns'], Rb=c['Rb'], k=c['k'], q=c['q'],
                            src_out=r['src_out'], amp_s=r.get('amp_s'), w_s=r.get('w_s'),
                            energy_residual=r['energy']['residual']))
            print(f"[{time.monotonic() - t0:5.0f} s] {c['structure']:>24} {c['label']:>24} k {c['k']:5.1f}: output {r['src_out']:.4e}", flush=True)
    table = []
    for stn in sorted({r['structure'] for r in res}):
        for label in ('isolated piece', 'dense ball (radius 3)', 'dilute ball (radius 30)'):
            rows = {r['k']: r for r in res if r['structure'] == stn and r['label'] == label}
            if 0.0 not in rows:
                continue
            base = rows[0.0]['src_out']
            table.append(dict(structure=stn, label=label, glow_rel_rest={f'{k:g}': rows[k]['src_out'] / base for k in sorted(rows)},
                              law_1_plus_k={f'{k:g}': 1 + k for k in sorted(rows)},
                              quiet_amplitude_rel_rest={f'{k:g}': (rows[k]['amp_s'] / rows[0.0]['amp_s']) if rows[0.0].get('amp_s') else None for k in sorted(rows)}))
    for t in table:
        print(t['structure'], '|', t['label'], '| glow ÷ rest', {k: round(v, 2) for k, v in t['glow_rel_rest'].items()},
              '| quiet amplitude ÷ rest', {k: (round(v, 3) if v else None) for k, v in t['quiet_amplitude_rel_rest'].items()})
    args.output.write_text(json.dumps(dict(experiment='round 18: the heat factor, isolated piece against dense and dilute balls',
                                           runs=res, table=table, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
