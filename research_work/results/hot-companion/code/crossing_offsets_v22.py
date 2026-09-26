"""Round 22, step 3: which way the heat of crossing moves the Bullet's lensing peaks, as the collision stack measures.

Harvey et al. (2015) stacked 72 collisions and measured how far the lensing moves toward the gas as the gas falls
behind its galaxies: beta = -0.04 +- 0.07 (the suite's test, regression/t_collisions.py; round 4 found beta = 0.03
without the heat of crossing). The suite's stack family holds the galaxies' separation fixed and only moves the gas, so
a heat of crossing with one history adds the same lensing to every member and cannot change its beta. What can change
it is heat that moves toward the gas as time goes on. That is measured here directly on the Bullet: the lensing peaks'
signed offsets from their galaxies along the collision axis (positive toward each system's own gas), with and without
the heat of crossing, for round 16's rule, the fast glow of round 21 (as proposed and with the energy booked), and the
frame versions of step 1 (code/crossing_frame_v22.py). Since the gas lag and the heat's drift both grow with the time
since the crossing, (offset with heat - offset without) / (galaxy-gas separation) estimates the change in beta.

    python code/crossing_offsets_v22.py --output run-crossing-frame-v22/crossing_offsets_v22.json
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
import crossing_frame_v22 as CF                   # noqa: E402
import hot_mode_speed_v21 as HM                   # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--cases', type=str, default='none,u,spread600,energy600')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    import bullet_static_v11 as BS
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    n, dx = 128, 22.5
    cross = dict(v_out=3900.0, v_in=3000.0, b=150.0)
    orig_make, orig_map = CH.make_kernel, CH.kappa_map_cross
    captured = {}

    def capture(*a, **k):
        res = orig_map(*a, **k)
        captured['pos'] = a[3]; captured['x'], captured['y'], captured['Se'] = res[0], res[1], res[2]
        return res
    out = dict(experiment='round 22: signed lensing-peak offsets toward the gas, Bullet, with and without the heat of crossing',
               law=law['name'], grid=dict(n=n, dx_kpc=dx), crossing=cross, runs=[])
    CH.kappa_map_cross = capture
    try:
        for case in args.cases.split(','):
            if case == 'none':
                CH.make_kernel = orig_make; c = None
            elif case == 'u':
                CH.make_kernel = orig_make; c = cross
            else:
                CH.make_kernel = HM.make_kernel_vh; c = cross
                HM.SETTINGS['v_h'] = float(case.replace('spread', '').replace('energy', ''))
                HM.SETTINGS['energy'] = case.startswith('energy'); HM.SETTINGS['power'] = 1.0
            r = CH.run_bullet_cross(law, dict(f), cross=c, n=n, dx=dx)
            kap = captured['Se'] / BS.static_sigma_crit()
            off = CF.signed_offsets(captured['x'], captured['y'], kap, captured['pos'], f['size'])
            pos = captured['pos']
            sep = {w: float(np.linalg.norm(np.asarray(pos[f'{w}_plasma']) - np.asarray(pos[f'{w}_bcg']))) for w in ('main', 'sub')}
            row = dict(case=case, m250_sub=r['m250_sub'], m250_main=r['m250_main'], peak_toward_gas_main=off['main'],
                       peak_toward_gas_sub=off['sub'], galaxy_gas_separation_kpc=sep, gas_main=r['gas_main'], gas_sub=r['gas_sub'])
            out['runs'].append(row)
            print(f"[{time.monotonic() - t0:5.0f} s] {case:10s}: sub {r['m250_sub']:.3e}, main {r['m250_main']:.3e}; peaks toward gas "
                  f"{off['main']:+.1f} (main, gas {sep['main']:.0f} kpc away), {off['sub']:+.1f} (sub, {sep['sub']:.0f})", flush=True)
            args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    finally:
        CH.make_kernel, CH.kappa_map_cross = orig_make, orig_map
        HM.SETTINGS.update(v_h=None, energy=False, power=1.0)
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
