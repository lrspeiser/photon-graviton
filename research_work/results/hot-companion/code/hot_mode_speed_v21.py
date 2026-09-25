"""Round 21, step 1: does a faster-spreading hot mode close the Bullet Cluster's gap?

A proposal supplied to the project (25 September 2026) keeps the law's u = 169.4 km/s for the heat weight and for the
cold companion's memory, but lets the heat generated in a collision spread at its own group speed v_h. It reports that
v_h of about 600 km/s puts both halves of the Bullet Cluster inside their measured lensing masses.

Round 16's crossing-heat model (code/crossing_heat_v16.py) is where that speed enters: the heat a star picked up a time
t ago now sits at a distance d from it, and the model sets d = u t. Here d = v_h t instead, in two versions:

  'spread'   only the timing changes (the proposal as stated);
  'energy'   the timing changes AND the energy is booked: a glow emitted at the same power but spreading v_h/u times
             faster is v_h/u times thinner (energy density = power / (4 pi d^2 speed)), so its heat weight carries a
             factor u / v_h. This is what the law's own energy balance (a = 2 l / u) implies for a mode that travels
             at v_h, unless the collision emits v_h/u times more power.
  'power'    the energy-booked version with the collision's power into the fast glow multiplied by m (--power m):
             weight x m u / v_h. 'spread' is the case m = v_h / u.

Everything else is round 16's model: the round-12 law in the static distances, the straight-pass history (receding at
3,900 km/s since pericentre, approaching at 3,000, impact parameter 150 kpc), the shares p of the other system's
companion, the settled heat of both systems unchanged.

    python code/hot_mode_speed_v21.py --grid coarse --output run-hot-mode-v21/hot_mode_speed_coarse.json
    python code/hot_mode_speed_v21.py --grid fine --speeds 600 --output run-hot-mode-v21/hot_mode_speed_fine.json
    python code/hot_mode_speed_v21.py --grid coarse --versions power --speeds 600,800 --power 2.5,3 \
        --output run-hot-mode-v21/hot_mode_speed_power.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import crossing_heat_v16 as CH                    # noqa: E402
import collisions_v10 as C10                      # noqa: E402

KPC_PER_KMS_GYR = 1.0227121650537077
SETTINGS = dict(v_h=None, energy=False, power=1.0)
_ORIGINAL = CH.make_kernel


def make_kernel_vh(u, I_other, I_own_at_shell, D_now, v_out, v_in, reach, kmax_d, b=0.0):
    """round 16's kernel with the heat's spreading speed v_h (and, in the 'energy' version, the factor u / v_h)."""
    vh = SETTINGS['v_h'] or u
    r_o, I_o = I_other
    scale = (SETTINGS['power'] * u / vh) if SETTINGS['energy'] else 1.0
    dmax = vh * 13.0 * KPC_PER_KMS_GYR

    def k_of_d(d):
        t = d / vh * CH.MYR_PER_KPC_PER_KMS
        D, v, cos, inside = CH.history(t, D_now, v_out, v_in, reach, b)
        Io = np.exp(np.interp(np.log(np.maximum(D, r_o[0])), np.log(r_o), np.log(I_o)))
        p = np.where(inside, Io / (Io + I_own_at_shell), 0.0)
        k = scale * p * np.maximum(v ** 2 - 2 * u * v * cos, 0.0) / u ** 2
        return np.where(d <= dmax, k, 0.0)
    return k_of_d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--grid', choices=['coarse', 'fine'], default='coarse')
    ap.add_argument('--speeds', type=str, default='169.4,300,450,600,800,1200')
    ap.add_argument('--versions', type=str, default='spread,energy')
    ap.add_argument('--power', type=str, default='1', help="for version 'power': comma-separated multipliers m")
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    n, dx = (128, 22.5) if args.grid == 'coarse' else (192, 15.0)
    CH.make_kernel = make_kernel_vh
    cross = dict(v_out=3900.0, v_in=3000.0, b=150.0)
    out = dict(experiment='round 21: the Bullet Cluster with a separate spreading speed for the collision heat', law=law['name'],
               grid=dict(n=n, dx_kpc=dx), crossing=cross, runs=[])
    SETTINGS['v_h'] = None
    base = CH.run_bullet_cross(law, dict(f), cross=None, n=n, dx=dx)
    out['no_crossing_heat'] = {k: base[k] for k in ('m250_sub', 'm250_main', 'peak_main', 'peak_sub', 'gas_main', 'gas_sub', 'kappa_sub', 'kappa_main', 'targets')}
    print(f"[{time.monotonic() - t0:5.0f} s] no crossing heat: sub {base['m250_sub']:.3e}, main {base['m250_main']:.3e}", flush=True)
    combos = []
    for version in args.versions.split(','):
        for vh in [float(s) for s in args.speeds.split(',')]:
            for m in ([float(x) for x in args.power.split(',')] if version == 'power' else [1.0]):
                combos.append((version, vh, m))
    for version, vh, m in combos:
        SETTINGS['v_h'] = vh; SETTINGS['energy'] = (version in ('energy', 'power')); SETTINGS['power'] = m
        r = CH.run_bullet_cross(law, dict(f), cross=cross, n=n, dx=dx)
        row = dict(version=version, v_h=vh, power=m, m250_sub=r['m250_sub'], m250_main=r['m250_main'], peak_main=r['peak_main'],
                   peak_sub=r['peak_sub'], gas_main=r['gas_main'], gas_sub=r['gas_sub'], kappa_sub=r['kappa_sub'],
                   kappa_main=r['kappa_main'], S_cross_share=r['report'].get('S_cross_over_S_mid_plane'), targets=r['targets'])
        out['runs'].append(row)
        print(f"[{time.monotonic() - t0:5.0f} s] {version:6s} v_h {vh:6.1f}" + (f" power x{m:g}" if version == 'power' else '') + f": sub {r['m250_sub']:.3e} (target {r['targets']['m250_sub'][0]:.2e}-"
              f"{r['targets']['m250_sub'][1]:.2e}), main {r['m250_main']:.3e} (target {r['targets']['m250_main'][0]:.2e}-{r['targets']['m250_main'][1]:.2e}); "
              f"peaks {r['peak_main']:.0f}, {r['peak_sub']:.0f} kpc; gas {r['gas_main']:.3f}, {r['gas_sub']:.3f}; crossing share {row['S_cross_share']:.2f}", flush=True)
        args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    CH.make_kernel = _ORIGINAL
    SETTINGS.update(v_h=None, energy=False, power=1.0)
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
