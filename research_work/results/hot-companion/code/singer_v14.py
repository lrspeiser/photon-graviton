"""Round 14, step 2: where "feed quiet waves, absorb loud ones" (round 13's rule) comes from.

Two physical versions of the emitter, dropped into the joint experiment of round 13 (code/coherent_force_v13.py):

A. A singer that also absorbs (c_abs): a locked self-sustained oscillation a quarter cycle ahead of its wave, plus a
   passive part i c0 E driven by the wave (energy taken in proportion to |E|^2). A lone one in a wave of strength E is
   pulled by (k/2)(|E| - c0 |E|^2): fed quiet waves, absorbed loud ones. With gamma_abs, the passive part is resonant.
B. A singer with a power balance (supply, loss): a fixed supply P_s (the law's l per kilogram) and a loss c_L |E|^2 to
   the wave around it; it must feed the wave the difference, so its timing is set like a synchronous machine's load
   angle: |E| sin(lead) = P_s - c_L |E|^2. A generator (ahead) in quiet waves, a motor (behind) in loud ones.

Parts (each writes its runs to the output JSON):
  a_rest      A at rest: three densities, four absorption strengths
  a_motion    A in the densest cloud (c0 = 1), moving freely, colliding (nu = 20 Gamma max(1, q)) or rotating
  a_stop      A: moving freely (q = 1/2) until t = 600, then at rest (does the cloud keep what motion gave it?)
  a_resonant  A with resonant passive parts (gamma_abs = Gamma): at rest and moving, up to q = 8
  b_rest      B at rest: three densities, four (P_s, E_s = sqrt(P_s/c_L)) pairs
  b_motion    B in the two dense clouds (P_s = 0.1, E_s = 0.1), the three motions, q = 1/32 ... 1
  b_stop      B: moving freely (q = 1/4) until t = 600, then at rest

    python code/singer_v14.py --output run-coherent-force-v13/singer_v14.json [--processes 4] [--parts a_rest,...]
"""
from __future__ import annotations
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_force_v13 import K0, simulate   # noqa: E402

BASE = dict(N=100, Gamma=0.1, T=900.0, burn=350.0, rp=6.0, P=32, w=0.01, delta=-np.pi / 2)
QS = (1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0)
MOTIONS = lambda q: (('free', 0.0), ('collisional', 20 * BASE['Gamma'] * max(1.0, q)), ('ordered', 0.0))


def configs(parts):
    runs = []
    def add(part, **kw):
        runs.append(dict(BASE, part=part, **kw))
    if 'a_rest' in parts:
        for Rb in (1.0, 0.5, 0.35):
            for c0 in (0.25, 0.5, 1.0, 2.0):
                for sd in (1, 2):
                    add('a_rest', Rb=Rb, c_abs=c0, kind='rest', sigma=0.0, nu=0.0, seed=sd, ratio=0.0)
    if 'a_motion' in parts:
        for sd in (1, 2, 3):
            add('a_motion', Rb=0.35, c_abs=1.0, kind='rest', sigma=0.0, nu=0.0, seed=sd, ratio=0.0)
            for q in QS:
                for kind, nu in MOTIONS(q):
                    add('a_motion', Rb=0.35, c_abs=1.0, kind=kind, sigma=q * BASE['Gamma'] / K0, nu=nu, seed=sd, ratio=q)
    if 'a_stop' in parts:
        for sd in (1, 2, 3, 4):
            add('a_stop', Rb=0.35, c_abs=1.0, kind='free', sigma=0.5 * BASE['Gamma'] / K0, nu=0.0, seed=sd, ratio=0.5,
                stop_at=600.0, T=4000.0, burn=3000.0)
    if 'a_resonant' in parts:
        for sd in (1, 2):
            add('a_resonant', Rb=0.35, c_abs=1.0, gamma_abs=0.1, kind='rest', sigma=0.0, nu=0.0, seed=sd, ratio=0.0, burn=450.0)
            for q in (1 / 4, 1.0, 2.0, 4.0, 8.0):
                for kind, nu in MOTIONS(q):
                    add('a_resonant', Rb=0.35, c_abs=1.0, gamma_abs=0.1, kind=kind, sigma=q * BASE['Gamma'] / K0, nu=nu, seed=sd,
                        ratio=q, burn=450.0)
    if 'b_rest' in parts:
        for Rb in (1.0, 0.5, 0.35):
            for Ps, Es in ((0.05, 0.1), (0.1, 0.1), (0.05, 0.2), (0.2, 0.3)):
                for sd in (1, 2):
                    add('b_rest', Rb=Rb, supply=Ps, loss=Ps / Es ** 2, E_switch=Es, kind='rest', sigma=0.0, nu=0.0, seed=sd, ratio=0.0)
    if 'b_motion' in parts:
        for Rb in (0.35, 0.5):
            for sd in (1, 2, 3, 4):
                add('b_motion', Rb=Rb, supply=0.1, loss=10.0, E_switch=0.1, kind='rest', sigma=0.0, nu=0.0, seed=sd, ratio=0.0)
                for q in QS:
                    for kind, nu in MOTIONS(q):
                        add('b_motion', Rb=Rb, supply=0.1, loss=10.0, E_switch=0.1, kind=kind, sigma=q * BASE['Gamma'] / K0, nu=nu,
                            seed=sd, ratio=q)
    if 'b_stop' in parts:
        for sd in (1, 2, 3, 4):
            add('b_stop', Rb=0.35, supply=0.1, loss=10.0, E_switch=0.1, kind='free', sigma=0.25 * BASE['Gamma'] / K0, nu=0.0, seed=sd,
                ratio=0.25, stop_at=600.0, T=4000.0, burn=3000.0)
    return runs


def cost(c):
    rates = [0.1, 0.05 / c['Gamma']] + ([0.2 / c['nu']] if c['nu'] > 0 else []) + ([0.1 / (K0 * c['sigma'])] if c['sigma'] > 0 else [])
    return c['T'] / min(rates) * (2.0 if c.get('c_abs') else 1.0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    ap.add_argument('--parts', default='a_rest,a_motion,a_stop,a_resonant,b_rest,b_motion,b_stop')
    args = ap.parse_args(); t0 = time.monotonic()
    runs = sorted(configs(set(args.parts.split(','))), key=lambda c: -cost(c))
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] {c['part']:10s} Rb {c['Rb']:.2f} {c['kind']:11s} q {c['ratio']:.4f} seed {c['seed']}:"
                  f" radiated/indep {r['radiated_over_independent']:.4f}, pull/indep {r['pull'] / r['independent_pull']:+.4f}, lock {r['locking']:.3f}", flush=True)
    for r in res:
        r['cfg']['delta_over_pi'] = r['cfg']['delta'] / np.pi
    args.output.write_text(json.dumps(dict(experiment='round 14, step 2: the singer that absorbs, and the singer with a power balance',
                                           runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
