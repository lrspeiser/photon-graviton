"""Round 23: the heat of crossing driven by the interference (beat) of the two systems' companions.

Round 16's heat of crossing (code/crossing_heat_v16.py) weighs the other system's flow by its share of the intensity a
star sits in, p = I_other / (I_other + I_own): the two companions are added as if they did not interfere. If they are
waves that do interfere, a star moving through both sees the combined intensity rise and fall at the rate it crosses
the other flow's crests (round 22, section 33.3: hundreds of times the companion's own rate). The depth of that
flicker, as a fraction of the mean intensity, is
    m = 2 sqrt(I_other I_own) / (I_other + I_own),
which is 1 when the two are equally strong and 2 sqrt(I_other / I_own) when the other is faint: interference is linear
in the fainter wave's amplitude, so a wave with 1% of the intensity makes the total flicker by 20% (the principle of
heterodyne detection). Here the heat of crossing is driven by that flicker, k = m (v^2 - 2 u v cos theta)/u^2 in place
of p (...), everything else round 16's: the straight-pass history, the other system's reach, the heat spreading from
the stars at u (the companion's memory) or at a faster v_h with its energy booked (round 21). No constant is refitted.

    python code/beat_heat_v23.py --cases share,beat --output run-beat-heat-v23/beat_bullet.json
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

KPC_PER_KMS_GYR = 1.0227121650537077
SETTINGS = dict(weight='share', v_h=None, energy=True, power=1.0)
_ORIG_MAKE, _ORIG_MAP = CH.make_kernel, CH.kappa_map_cross


def weight(Io, I_own):
    if SETTINGS['weight'] == 'beat':
        return 2 * np.sqrt(Io * I_own) / (Io + I_own)
    return Io / (Io + I_own)


def make_kernel_beat(u, I_other, I_own_at_shell, D_now, v_out, v_in, reach, kmax_d, b=0.0):
    """round 16's k(d) with the share p replaced by the beat depth m (SETTINGS['weight'] = 'beat'), spreading at v_h
    (energy booked: weight x u / v_h) when SETTINGS['v_h'] is set."""
    vh = SETTINGS['v_h'] or u
    r_o, I_o = I_other
    scale = (SETTINGS['power'] * u / vh) if SETTINGS['energy'] else 1.0
    dmax = vh * 13.0 * KPC_PER_KMS_GYR

    def k_of_d(d):
        t = d / vh * CH.MYR_PER_KPC_PER_KMS
        D, v, cos, inside = CH.history(t, D_now, v_out, v_in, reach, b)
        Io = np.exp(np.interp(np.log(np.maximum(D, r_o[0])), np.log(r_o), np.log(I_o)))
        w = np.where(inside, weight(Io, I_own_at_shell), 0.0)
        k = scale * w * np.maximum(v ** 2 - 2 * u * v * cos, 0.0) / u ** 2
        return np.where(d <= dmax, k, 0.0)
    return k_of_d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--cases', type=str, default='share,beat,beat:300,beat:600',
                    help="'share' (round 16), 'beat', 'beat:<v_h>' (spreading at v_h with the energy booked), or 'beat:<v_h>:<power>' (the power into the glow multiplied)")
    ap.add_argument('--grid', choices=['coarse', 'fine'], default='coarse')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    import bullet_static_v11 as BS
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    n, dx = (128, 22.5) if args.grid == 'coarse' else (192, 15.0)
    cross = dict(v_out=3900.0, v_in=3000.0, b=150.0)
    captured = {}

    def capture(*a, **k):
        res = _ORIG_MAP(*a, **k)
        captured['pos'] = a[3]; captured['x'], captured['y'], captured['Se'] = res[0], res[1], res[2]
        return res
    out = dict(experiment='round 23: the heat of crossing driven by the beat of the two companions (Bullet)', law=law['name'],
               grid=dict(n=n, dx_kpc=dx), crossing=cross, runs=[])
    CH.make_kernel, CH.kappa_map_cross = make_kernel_beat, capture
    try:
        for case in args.cases.split(','):
            parts = case.split(':')
            name, vh, pw = parts[0], (parts[1] if len(parts) > 1 else ''), (parts[2] if len(parts) > 2 else '1')
            SETTINGS.update(weight=name, v_h=float(vh) if vh else None, energy=True, power=float(pw))
            r = CH.run_bullet_cross(law, dict(f), cross=cross, n=n, dx=dx)
            kap = captured['Se'] / BS.static_sigma_crit()
            off = CF.signed_offsets(captured['x'], captured['y'], kap, captured['pos'], f['size'])
            row = dict(case=case, weight=name, v_h=SETTINGS['v_h'] or law['u_kms'], power=SETTINGS['power'], m250_sub=r['m250_sub'], m250_main=r['m250_main'],
                       peak_main=r['peak_main'], peak_sub=r['peak_sub'], peak_toward_gas_main=off['main'], peak_toward_gas_sub=off['sub'],
                       gas_main=r['gas_main'], gas_sub=r['gas_sub'], kappa_main=r['kappa_main'], kappa_sub=r['kappa_sub'],
                       S_cross_share=r['report'].get('S_cross_over_S_mid_plane'), shells=r['crossing'].get('shells'),
                       targets=r['targets'], peak_limits=r['peak_limits'])
            out['runs'].append(row)
            print(f"[{time.monotonic() - t0:5.0f} s] {case:9s}: sub {r['m250_sub']:.3e} (target {r['targets']['m250_sub'][0]:.2e}-"
                  f"{r['targets']['m250_sub'][1]:.2e}), main {r['m250_main']:.3e} (target {r['targets']['m250_main'][0]:.2e}-"
                  f"{r['targets']['m250_main'][1]:.2e}); peaks {r['peak_main']:.0f}, {r['peak_sub']:.0f} kpc (toward gas {off['main']:+.0f}, "
                  f"{off['sub']:+.0f}); gas {r['gas_main']:.3f}, {r['gas_sub']:.3f}; crossing share {row['S_cross_share']:.2f}", flush=True)
            args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    finally:
        CH.make_kernel, CH.kappa_map_cross = _ORIG_MAKE, _ORIG_MAP
        SETTINGS.update(weight='share', v_h=None, energy=True, power=1.0)
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
