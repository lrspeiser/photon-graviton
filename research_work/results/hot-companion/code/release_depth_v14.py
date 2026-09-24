"""Round 14, step 3: how the released power grows with the speed spread, and when test bodies lose step.

The singer with a power balance (code/singer_v14.py, part B: supply P_s = 0.1, loss c_L = 10, switch at |E| = 0.1),
in three clouds of N = 100 of increasing density (Rb = 0.5, 0.35, 0.25 wavelengths: 191, 557 and 1,528 per cubic
wavelength), moving freely or colliding (nu = 20 Gamma max(1, q)), q = k0 sigma / Gamma = 1/64 ... 2, three seeds.
For each: the power released by the motion (radiated - radiated at rest), its growth per doubling of sigma, the wave at
the test bodies, their feeding (how well they keep in step: 1 = perfectly) and the pull.

    python code/release_depth_v14.py --output run-coherent-force-v13/release_depth_v14.json [--processes 4]
"""
from __future__ import annotations
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_force_v13 import K0, simulate   # noqa: E402

BASE = dict(N=100, Gamma=0.1, T=900.0, burn=350.0, rp=6.0, P=32, w=0.01, delta=-np.pi / 2, supply=0.1, loss=10.0)
QS = (1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0, 2.0)


def configs():
    runs = []
    for Rb in (0.5, 0.35, 0.25):
        for sd in (1, 2, 3):
            runs.append(dict(BASE, Rb=Rb, kind='rest', sigma=0.0, nu=0.0, seed=sd, ratio=0.0))
            for q in QS:
                for kind, nu in (('free', 0.0), ('collisional', 20 * BASE['Gamma'] * max(1.0, q))):
                    runs.append(dict(BASE, Rb=Rb, kind=kind, sigma=q * BASE['Gamma'] / K0, nu=nu, seed=sd, ratio=q))
    return runs


def cost(c):
    rates = [0.1, 0.05 / c['Gamma']] + ([0.2 / c['nu']] if c['nu'] > 0 else []) + ([0.1 / (K0 * c['sigma'])] if c['sigma'] > 0 else [])
    return c['T'] / min(rates)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    args = ap.parse_args(); t0 = time.monotonic()
    runs = sorted(configs(), key=lambda c: -cost(c)); res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] Rb {c['Rb']:.2f} {c['kind']:11s} q {c['ratio']:.4f} seed {c['seed']}: radiated/indep "
                  f"{r['radiated_over_independent']:.4f}, pull/indep {r['pull'] / r['independent_pull']:+.4f}, test feeding {r['feeding_test_bodies']:+.3f}", flush=True)
    for r in res:
        r['cfg']['delta_over_pi'] = r['cfg']['delta'] / np.pi
    args.output.write_text(json.dumps(dict(experiment='round 14, step 3: the release against the speed spread, in three densities (power-balance singers)',
                                           runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
