#!/usr/bin/env python3
"""Round 4: galaxy speeds in the Bullet Cluster predicted by our law, as the observers measure them.

    python bullet_speeds_v4.py --output-dir ../run-bullet-speeds-v4

* Main cluster: the line-of-sight speed spread averaged over a projected aperture (Barrena et al.
  2002 measured 1,249 +109/-100 km/s from 71 galaxies; their virial mass uses about 1.5 Mpc), for
  star mass-to-light ratios 1, 1.5 and 2 (Clowe et al.'s masses assume 2).
* Subcluster: its stars' speeds before the collision (they keep them), for pre-collision mass
  ratios 1:6, 1:8 and 1:10. Barrena et al.'s 7 galaxies give 212 +67/-52 km/s; from its X-ray
  temperature and luminosity they infer about 700 km/s and a pre-merger ratio of about 1:6.
Isotropic Jeans equation in the gravity our law makes from each cluster's own visible matter.
"""
from __future__ import annotations
import argparse, copy, json, time
from pathlib import Path
import numpy as np
import bullet_v3 as B
import bullet_v4 as V


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    pos = B.positions(); comps, _ = B.fit_components(pos)
    main_rows, sub_rows = [], []
    for ml in (1.0, 1.5, 2.0):
        cm = copy.deepcopy(comps)
        for k in ('st_main', 'st_sub'): cm[k]['M'] = comps[k]['M'] * ml / 2.0
        r, s, _ = V.own_sigma_multi([cm['gas_main']], [cm['st_main']], consts)
        dms = V.shells_on(r, np.gradient(r), cm['st_main'])
        ap_avg = {}
        for R in (500, 1000, 1500):
            f = np.where(r <= R, 1.0, 1 - np.sqrt(np.clip(1 - (R / r) ** 2, 0, 1)))    # share of each shell inside the cylinder
            ap_avg[f'{R}kpc'] = float(np.sqrt(np.sum(dms * f * s ** 2) / np.sum(dms * f)))
        main_rows.append(dict(mass_to_light=ml, sigma_r_kms={f'{x}kpc': float(np.interp(x, r, s)) for x in (100, 400, 800)},
                              line_of_sight_inside_aperture_kms=ap_avg))
        print(f'main, M/L {ml}: line-of-sight average inside 0.5/1/1.5 Mpc = ' + ' / '.join(f'{v:.0f}' for v in ap_avg.values())
              + ' km/s (measured 1,249 +109/-100)', flush=True)
        if ml == 1.5:
            for ratio in (6.0, 8.0, 10.0):
                gg, sats, pre = V.pre_collision_models(cm, ratio=ratio)
                rs, ss, _ = V.own_sigma_multi([dict(cm['gas_sub']), gg['gas_atm_ghost']], [cm['st_sub'], sats], consts)
                sub_rows.append(dict(ratio=ratio, mass_to_light=ml, sigma_before_kms={f'{x}kpc': float(np.interp(x, rs, ss)) for x in (20, 50, 100, 200)}))
                print(f'subcluster before the collision, 1:{ratio:.0f}: ' + ', '.join(f'{v:.0f}' for v in sub_rows[-1]['sigma_before_kms'].values())
                      + ' km/s at 20/50/100/200 kpc', flush=True)
    (out / 'speeds.json').write_text(json.dumps(dict(main=main_rows, sub_before=sub_rows,
        observed=dict(main='1249 +109/-100 km/s, 71 galaxies (Barrena et al. 2002)', sub='212 +67/-52 km/s, 7 galaxies; ~700 km/s implied by its X-ray temperature and luminosity'),
        note='the 1.5 Mpc aperture is pulled down by the model galaxies being cut off at 1.5 Mpc',
        constants=consts, seconds=time.monotonic() - t0), indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
