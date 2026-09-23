#!/usr/bin/env python3
"""The collision family again, now with the companion's memory (round 4).

    python collisions_v4.py --output-dir ../run-collisions-v4

Same 20 pieces of substructure as collisions_v3.py (the Bullet as a template; the gas lagging
its galaxies by 40-300 kpc; subcluster mass x1 and x3). New in round 4:
* the old companion, emitted before the collision, rides with the galaxies (the ghost of each
  settled pre-collision cluster);
* a fresh sphere of radius u t is rebuilt around the gas where it is now, where t is the time
  since the gas was stopped, t = lag / (gas-galaxy separation speed, 1,240 km/s as in the Bullet:
  190 kpc in 150 Myr).
Because u = 197 km/s is far below separation speeds of order 1,000 km/s, the fresh sphere never
catches up with the lag: the lensing should stay with the galaxies at every stage.

Measurement. Peaks are found as in round 3 (30 kpc smoothing) and then refined below the pixel
with a quadratic fit. A lag-0 run (gas left on its galaxies) is the reference: the peak of a
summed map is pulled toward the other cluster whether or not the gas has moved, and Harvey et
al. fit each clump separately, which removes that. The number to compare with their
beta = -0.04 +- 0.07 is how far the lensing moves when the gas moves: beta_gas = ds / d(lag).
"""
from __future__ import annotations
import argparse, copy, json, time
from pathlib import Path
import numpy as np
import bullet_v3 as B
import bullet_v4 as V
import collisions_v3 as C3

SEP_SPEED = 1240.0      # km/s
LAGS = (0, 40, 80, 150, 209, 300)


def refined_peaks(x, y, kap, pos, smooth_kpc=30.0):
    """B.peaks, with each peak moved below the pixel by a quadratic fit to its 3x3 neighbourhood."""
    from scipy.ndimage import gaussian_filter, maximum_filter
    dx = x[1] - x[0]
    ks = gaussian_filter(kap, smooth_kpc / dx)
    mx = (ks == maximum_filter(ks, size=7)) & (ks > 0.25 * ks.max())
    out = []
    for i, j in zip(*np.nonzero(mx)):
        if not (0 < i < len(x) - 1 and 0 < j < len(y) - 1): continue
        w = ks[i - 1:i + 2, j - 1:j + 2]
        gx, gy = (w[2, 1] - w[0, 1]) / 2, (w[1, 2] - w[1, 0]) / 2
        hxx, hyy = w[2, 1] - 2 * w[1, 1] + w[0, 1], w[1, 2] - 2 * w[1, 1] + w[1, 0]
        hxy = (w[2, 2] - w[2, 0] - w[0, 2] + w[0, 0]) / 4
        try: d = -np.linalg.solve([[hxx, hxy], [hxy, hyy]], [gx, gy])
        except np.linalg.LinAlgError: d = np.zeros(2)
        d = np.clip(d, -1, 1)
        pt = np.array([x[i] + d[0] * dx, y[j] + d[1] * dx]); pix = np.array([x[i], y[j]])
        out.append(dict(x=float(pt[0]), y=float(pt[1]), kappa=float(ks[i, j]),
                        **{f'dist_{k}': float(np.linalg.norm(pt - v)) for k, v in pos.items()},
                        **{f'pixdist_{k}': float(np.linalg.norm(pix - v)) for k, v in pos.items()}))
    return sorted(out, key=lambda d: -d['kappa'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=160); ap.add_argument('--dx', type=float, default=18.0)
    ap.add_argument('--ml', type=float, default=1.5)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    u = consts['u_kms']
    base_pos = B.positions(); comps0, _ = B.fit_components(base_pos)
    comps = copy.deepcopy(comps0)
    for k in ('st_main', 'st_sub'): comps[k]['M'] *= args.ml / 2.0
    rows = []
    axis = (base_pos['sub_bcg'] - base_pos['main_bcg']) / np.linalg.norm(base_pos['sub_bcg'] - base_pos['main_bcg'])
    for lag, ms in [(lg, m) for lg in LAGS for m in (1.0, 3.0)]:
        c, pos = C3.scenario(comps, base_pos, lag, lag, sub_mass_scale=ms)
        gg, sats, pre = V.pre_collision_models(c, ratio=8.0 / ms)
        rm, sm, _ = V.own_sigma_multi([c['gas_main']], [c['st_main']], consts)
        rs, ss, _ = V.own_sigma_multi([dict(c['gas_sub']), gg['gas_atm_ghost']], [c['st_sub'], sats], consts)
        gst = [(dict(c['st_main']), (rm, sm)), (dict(c['st_sub']), (rs, ss))] + ([(sats, (rs, ss))] if sats['M'] > 0 else [])
        t_myr = lag / SEP_SPEED * 977.8                     # kpc/(km/s) -> Myr
        fresh = u * t_myr / 977.8
        x, y, Se, Sb = V.kappa_map_v4(c, gg, gst, pos, consts, n=args.n, dx=args.dx, fresh_kpc=fresh)
        kap = Se / B.sigma_crit(); pk = refined_peaks(x, y, kap, pos, smooth_kpc=30.0)
        row = dict(lag_kpc=lag, sub_mass_scale=ms, time_since_stop_myr=t_myr, fresh_sphere_kpc=fresh)
        for side, st, gs, sgn in (('main', 'main_bcg', 'main_plasma', 1.0), ('sub', 'sub_bcg', 'sub_plasma', -1.0)):
            near = [p for p in pk if p[f'dist_{st}'] < lag + 60]
            if near:
                p = min(near, key=lambda d: d[f'dist_{st}'])
                rel = np.array([p['x'], p['y']]) - pos[st]
                row[f'along_kpc_{side}'] = float(sgn * rel @ axis)            # + toward its own gas (and the other cluster)
                row[f'across_kpc_{side}'] = float(abs(rel[0] * axis[1] - rel[1] * axis[0]))
                row[f'offset_kpc_{side}'] = p[f'dist_{st}']; row[f'pixel_offset_kpc_{side}'] = p[f'pixdist_{st}']
                row[f'beta_raw_{side}'] = row[f'along_kpc_{side}'] / lag if lag > 0 else None
            else:
                for k in ('along_kpc', 'across_kpc', 'offset_kpc', 'pixel_offset_kpc', 'beta_raw'): row[f'{k}_{side}'] = None
            row[f'separate_peak_on_{side}_gas'] = bool(lag > 0 and [p for p in pk if p[f'dist_{gs}'] < 40 and p[f'dist_{st}'] > 60])
        rows.append(row)
        print(f"lag {lag:3d} kpc (t = {t_myr:4.0f} Myr, fresh sphere {fresh:4.1f} kpc), sub x{ms:.0f}: peak along the axis toward the gas"
              f" main {row['along_kpc_main']:+.1f}, sub {row['along_kpc_sub']:+.1f} kpc (pixel offsets {row['pixel_offset_kpc_main']:.1f}, {row['pixel_offset_kpc_sub']:.1f});"
              f" separate peak on the gas: {row['separate_peak_on_main_gas']}, {row['separate_peak_on_sub_gas']}", flush=True)
    # how far the lensing moves when the gas moves, measured from the lag-0 reference
    ref = {(r['sub_mass_scale'], sd): r[f'along_kpc_{sd}'] for r in rows if r['lag_kpc'] == 0 for sd in ('main', 'sub')}
    beta_gas, slopes = [], {}
    for r in rows:
        for sd in ('main', 'sub'):
            if r['lag_kpc'] > 0 and r[f'along_kpc_{sd}'] is not None:
                r[f'beta_gas_{sd}'] = (r[f'along_kpc_{sd}'] - ref[(r['sub_mass_scale'], sd)]) / r['lag_kpc']
                beta_gas.append(r[f'beta_gas_{sd}'])
    for ms in (1.0, 3.0):
        for sd in ('main', 'sub'):
            L = np.array([r['lag_kpc'] for r in rows if r['sub_mass_scale'] == ms]); S = np.array([r[f'along_kpc_{sd}'] for r in rows if r['sub_mass_scale'] == ms])
            slopes[f'{sd}_x{ms:.0f}'] = dict(slope=float(np.polyfit(L, S, 1)[0]), offset_at_lag0_kpc=float(ref[(ms, sd)]),
                                           spread_kpc=float(S.max() - S.min()))
    offs = [r[k] for r in rows if r['lag_kpc'] > 0 for k in ('offset_kpc_main', 'offset_kpc_sub') if r[k] is not None]
    pix = [r[k] for r in rows if r['lag_kpc'] > 0 for k in ('pixel_offset_kpc_main', 'pixel_offset_kpc_sub') if r[k] is not None]
    betas = [r[k] for r in rows for k in ('beta_raw_main', 'beta_raw_sub') if r[k] is not None]
    summary = dict(n=len(offs), offset_kpc_median=float(np.median(offs)), offset_kpc_range=[float(min(offs)), float(max(offs))],
                   pixel_offset_kpc_median=float(np.median(pix)), pixel_offset_kpc_range=[float(min(pix)), float(max(pix))],
                   beta_raw_median=float(np.median(betas)),
                   beta_gas_median=float(np.median(beta_gas)), beta_gas_mean=float(np.mean(beta_gas)),
                   beta_gas_range=[float(min(beta_gas)), float(max(beta_gas))], slopes=slopes,
                   gas_peaks=sum(1 for r in rows for k in ('separate_peak_on_main_gas', 'separate_peak_on_sub_gas') if r[k]),
                   observed='lensing within 5.8 +- 8.2 kpc of the galaxies; beta = -0.04 +- 0.07 (Harvey et al. 2015)')
    (out / 'collisions_v4.json').write_text(json.dumps(dict(rows=rows, summary=summary, mass_to_light=args.ml, sep_speed=SEP_SPEED,
                                                             constants=consts, seconds=time.monotonic() - t0), indent=2, default=float) + '\n')
    print(summary)


if __name__ == '__main__':
    main()
