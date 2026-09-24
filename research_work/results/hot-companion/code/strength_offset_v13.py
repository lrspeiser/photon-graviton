"""Round 13, part 2: the joint experiment with an offset that depends on the strength of the local wave.

Part 1 (code/coherent_force_v13.py, code/joint_checks_v13.py) showed that with one fixed offset, a locked test body
is pulled exactly when the cloud is brighter than independent emitters, so random motion can only weaken the pull.
A single rule escapes this if an emitter's offset differs between the strong waves inside a body and the weak ones
arriving from far away. The simplest such rule, for sources and test bodies alike:
    delta(|E|) = -pi/2 + pi x^4 / (1 + x^4),   x = |E| / E_s
a quarter cycle ahead (feeding) in weak waves, a quarter cycle behind (absorbing) in strong ones.

1. rest: the cloud at rest for several E_s and densities (dark cloud, pulled test bodies?).
2. motion: the dark dense clouds (N = 100 in a ball of radius 0.5 or 0.35 wavelengths) with E_s = 0.1, at rest and
   moving freely, colliding (each velocity redrawn 20 times per locking time or per Doppler time, whichever is
   shorter: nu = 20 Gamma max(1, q)) or rotating, at q = k0 sigma / Gamma = 1/32 ... 1, four seeds.

    python code/strength_offset_v13.py --output run-coherent-force-v13/strength_offset_v13.json [--processes 4]
"""
from __future__ import annotations
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_force_v13 import K0, simulate   # noqa: E402

BASE = dict(N=100, Gamma=0.1, T=900.0, burn=350.0, rp=6.0, P=32, w=0.01, delta=-np.pi / 2, delta_strong=np.pi / 2, n_s=4)
QS = (1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0)


def configs():
    runs = []
    for Rb in (1.0, 0.5, 0.35):
        for Es in (0.03, 0.1, 0.3, 1.0, 3.0):
            for sd in (1, 2):
                runs.append(dict(BASE, part='rest', Rb=Rb, E_s=Es, seed=sd, kind='rest', sigma=0.0, nu=0.0, ratio=0.0))
    for Rb in (0.5, 0.35):
        for sd in (1, 2, 3, 4):
            runs.append(dict(BASE, part='motion', Rb=Rb, E_s=0.1, seed=sd, kind='rest', sigma=0.0, nu=0.0, ratio=0.0))
            for q in QS:
                sig = q * BASE['Gamma'] / K0
                for kind, nu in (('free', 0.0), ('collisional', 20 * BASE['Gamma'] * max(1.0, q)), ('ordered', 0.0)):
                    runs.append(dict(BASE, part='motion', Rb=Rb, E_s=0.1, seed=sd, kind=kind, sigma=sig, nu=nu, ratio=q))
    return runs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    args = ap.parse_args(); t0 = time.monotonic()
    runs = configs(); res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] {c['part']:6s} Rb {c['Rb']:.2f} E_s {c['E_s']:.2f} {c['kind']:11s} q {c['ratio']:.4f} seed {c['seed']}:"
                  f" radiated/indep {r['radiated_over_independent']:.4f}, pull/indep {r['pull'] / r['independent_pull']:+.4f},"
                  f" feeding src {r['feeding_sources']:+.2f} test {r['feeding_test_bodies']:+.2f}", flush=True)
    for r in res:
        r['cfg']['delta_over_pi'] = r['cfg']['delta'] / np.pi
    args.output.write_text(json.dumps(dict(experiment='round 13, part 2: an offset that depends on the local wave strength',
                                           runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
