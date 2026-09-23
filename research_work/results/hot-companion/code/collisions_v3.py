#!/usr/bin/env python3
"""A family of colliding clusters under round-3 hot-companion gravity, against the stack of
72 collisions of Harvey et al. (2015, Science 347, 1462).

    python collisions_v3.py --output-dir ../run-collisions-v3

Harvey et al. measure, for each piece of substructure, the separation from its galaxies to its
gas (delta_SG) and from its galaxies to its lensing mass along that line (delta_SI). Their
stack gives <delta_SI> = -5.8 +- 8.2 kpc and a fractional offset beta = delta_SI/delta_SG =
-0.04 +- 0.07 (68%), after correcting the lensing position for the pull of the gas (on
average 4.3 +- 1.6 kpc); without that correction their limits loosen only slightly.

We use the Bullet Cluster's fitted gas and stars (bullet_v3.py) as a template and vary the
stage of the collision (how far each gas cloud lags its galaxies) and the mass ratio. Star
speeds come from our law and the visible matter (the self-consistent Jeans profiles), so
nothing is tuned. For each piece of substructure we report beta_raw, the lensing peak's
position between the galaxies (0) and the gas (1); the comparable observed value is beta
before the gas correction, about 0.0 +- 0.07.
"""
from __future__ import annotations
import argparse, json, time, copy
from pathlib import Path
import numpy as np
import bullet_v3 as B


def scenario(comps, base_pos, lag_main, lag_sub, sub_mass_scale=1.0, separation=None):
    """Move each gas cloud to lag its galaxies by the given distance along the collision axis."""
    c = copy.deepcopy(comps)
    mb, sb = base_pos['main_bcg'], base_pos['sub_bcg']
    if separation is not None:
        sb = mb + (sb - mb) / np.linalg.norm(sb - mb) * separation
    axis = (sb - mb) / np.linalg.norm(sb - mb)
    pos = dict(main_bcg=mb, sub_bcg=sb, main_plasma=mb + axis * lag_main, sub_plasma=sb - axis * lag_sub)
    for k in ('gas_sub', 'st_sub'):
        c[k]['M'] *= sub_mass_scale
    return c, pos


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=160); ap.add_argument('--dx', type=float, default=18.0)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    base_pos = B.positions(); comps, _ = B.fit_components(base_pos)
    Mmain_b = comps['gas_main']['M'] + comps['st_main']['M']
    rows = []
    grid = [(lm, ls, ms) for lm, ls in ((40, 40), (80, 80), (150, 150), (209, 193), (300, 300)) for ms in (1.0, 3.0)]
    for lag_main, lag_sub, ms in grid:
        c, pos = scenario(comps, base_pos, lag_main, lag_sub, sub_mass_scale=ms)
        prof_main = B.own_sigma_profile(c['gas_main'], c['st_main'], consts)
        lost = dict(kind='beta', M=max(ms * Mmain_b / 8 - c['gas_sub']['M'] - c['st_sub']['M'], 1e10), scale=150.0, rt=1000.0)
        prof_sub = B.own_sigma_profile(c['gas_sub'], c['st_sub'], consts, extra_gas=lost)
        res = B.kappa_map(c, pos, consts, prof_main, prof_sub, n=args.n, dx=args.dx, direction='mix')
        x, y, Sig_eff = res[0], res[1], res[2]
        kap = Sig_eff / B.sigma_crit()
        pk = B.peaks(x, y, kap, pos, smooth_kpc=30.0)
        row = dict(lag_main_kpc=lag_main, lag_sub_kpc=lag_sub, sub_mass_scale=ms)
        for side, st, gs in (('main', 'main_bcg', 'main_plasma'), ('sub', 'sub_bcg', 'sub_plasma')):
            near = [p for p in pk if p[f'dist_{st}'] < np.linalg.norm(pos[gs] - pos[st]) + 60]   # matched as Harvey et al.: nearest clump
            if near:
                p = min(near, key=lambda d: d[f'dist_{st}'])
                row[f'beta_raw_{side}'] = B.along_axis_fraction(p, pos[st], pos[gs])
                row[f'offset_kpc_{side}'] = float(p[f'dist_{st}'])
            else:
                row[f'beta_raw_{side}'] = None; row[f'offset_kpc_{side}'] = None
            on_gas = [p for p in pk if p[f'dist_{gs}'] < 40]
            row[f'separate_peak_on_{side}_gas'] = bool(on_gas)
        rows.append(row)
        print(f"lag main {lag_main:3d} / sub {lag_sub:3d} kpc, sub mass x{ms:.0f}: beta_raw main {row['beta_raw_main']}, sub {row['beta_raw_sub']}"
              f"; separate gas peaks main {row['separate_peak_on_main_gas']}, sub {row['separate_peak_on_sub_gas']}", flush=True)
    betas = [r[k] for r in rows for k in ('beta_raw_main', 'beta_raw_sub') if r[k] is not None]
    summary = dict(n_substructures=len(betas), beta_raw_median=float(np.median(betas)), beta_raw_mean=float(np.mean(betas)),
                   beta_raw_range=[float(np.min(betas)), float(np.max(betas))],
                   missing_peaks=sum(1 for r in rows for k in ('beta_raw_main', 'beta_raw_sub') if r[k] is None),
                   observed='beta = -0.04 +- 0.07 after the gas correction; about 0.0 +- 0.07 before it (Harvey et al. 2015)')
    (out / 'collisions.json').write_text(json.dumps(dict(rows=rows, summary=summary, constants=consts,
                                                          seconds=time.monotonic() - t0), indent=2, default=float) + '\n')
    print(summary)


if __name__ == '__main__':
    main()
