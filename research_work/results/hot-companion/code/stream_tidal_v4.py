#!/usr/bin/env python3
"""Round 4: two effects examined for the Bullet's subcluster before the memory was found.

    python stream_tidal_v4.py --output-dir ../run-stream-tidal-v4

1. Stream heat. The heat weight is the mean-square speed of free-streaming matter about its
   LOCAL mean motion. Where the two galaxy streams overlap at relative speed V, that adds
   rho_m rho_s V^2 / (rho_m + rho_s) to rho <|v - v_mean|^2> (streaming along one axis counts
   V^2, not 3 V^2). Evaluated on the instantaneous field (no memory), as an upper bound on what
   the collision could have switched on so far.
2. Tidal shaking. As the subcluster crossed the main cluster, the main cluster's tides (under
   our law's pull) kicked its stars. The impulse approximation along a straight path, with a
   150 kpc miss distance, gives the kick at 50-200 kpc from the subcluster's centre. A fast pass
   gives every star at one place the same kick: ordered motion that becomes heat only after
   phase mixing (about an orbital time).
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
import bullet_v3 as B
import law as L

G = L.G


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    pos = B.positions(); comps, _ = B.fit_components(pos)
    Mmain_b = comps['gas_main']['M'] + comps['st_main']['M']
    pm = B.own_sigma_profile(comps['gas_main'], comps['st_main'], consts)
    lost = dict(kind='beta', M=Mmain_b / 8 - comps['gas_sub']['M'] - comps['st_sub']['M'], scale=150.0, rt=1000.0)
    ps = B.own_sigma_profile(comps['gas_sub'], comps['st_sub'], consts, extra_gas=lost)

    stream = []
    for V in (None, 2600., 3000., 4700.):
        x, y, Se = B.kappa_map(comps, pos, consts, pm, ps, direction='mix', v_rel=V)[:3]
        kap = Se / B.sigma_crit(); dec = B.clowe_decomposition(x, y, kap, pos)
        pk = B.peaks(x, y, kap, pos)
        sub = min(pk, key=lambda d: d['dist_sub_bcg']); main_ = min(pk, key=lambda d: d['dist_main_bcg'])
        row = dict(v_rel_kms=V, main=dec['main_bcg'], sub=dec['sub_bcg'], main_gas=dec['main_plasma'], sub_gas=dec['sub_plasma'],
                   main_peak_from_bcg_kpc=main_['dist_main_bcg'], nearest_peak_to_sub_bcg_kpc=sub['dist_sub_bcg'],
                   that_peak_from_sub_gas_kpc=sub['dist_sub_plasma'])
        stream.append(row)
        print(f"stream heat V = {V}: main {row['main']:.3f}  sub {row['sub']:.3f}  gas {row['main_gas']:.3f}/{row['sub_gas']:.3f};"
              f" nearest peak to the sub's galaxies {row['nearest_peak_to_sub_bcg_kpc']:.0f} kpc away", flush=True)

    # tidal kicks from the main cluster's pull under our law (spherical, its own gas and stars)
    r = np.geomspace(1.0, 5000.0, 800); dr = np.gradient(r)
    def shells(c):
        prof = B.rho_beta if c['kind'] == 'beta' else B.rho_nfw
        return c['M'] / B.mass3d(prof, (c['scale'],), c['rt']) * prof(r, c['scale'], c['rt']) * 4 * np.pi * r ** 2 * dr
    dmg, dms = shells(comps['gas_main']), shells(comps['st_main'])
    rm, sm = pm
    k = 3 * np.interp(r, rm, sm) ** 2 / u ** 2
    W = L.shell_weights(r, r); S = G * (W @ (k * dms)) / r ** 2; gN = G * np.cumsum(dmg + dms) / r ** 2
    g = L.total(gN, S, a, lam)
    gi = lambda X: np.interp(X, r, g)
    d_now = float(np.linalg.norm(pos['sub_bcg'] - pos['main_bcg']))
    tidal = []
    for V in (2600., 3000., 4700.):
        b = 150.0
        f = lambda s: gi(np.hypot(s, b)) / np.hypot(s, b)            # (dPhi/dR)/R: the squeeze across the path
        Cperp = (quad(f, -6000, 0, limit=400)[0] + quad(f, 0, d_now, limit=400)[0]) / V
        cpar = gi(np.hypot(d_now, b)) * d_now / np.hypot(d_now, b) / V
        kicks = {f'{x}kpc': float(np.sqrt((Cperp ** 2 * 2 / 3 + cpar ** 2 / 3) * x ** 2 / 3)) for x in (50, 100, 200)}
        tidal.append(dict(v_kms=V, miss_kpc=b, squeeze_kms_per_kpc=Cperp, stretch_kms_per_kpc=cpar, kick_kms=kicks))
        print(f'tidal V = {V:.0f} km/s, miss 150 kpc: kicks at 50/100/200 kpc ' + ', '.join(f'{v:.0f}' for v in kicks.values()) + ' km/s', flush=True)

    (out / 'stream_tidal_v4.json').write_text(json.dumps(dict(
        stream_heat=stream, tidal=tidal,
        main_pull_law_code_units={f'{x}kpc': float(gi(x)) for x in (50, 200, 700, 1500)},
        notes='stream heat evaluated on the instantaneous field (no memory); tidal kicks are ordered, heat only after phase mixing',
        constants=consts, seconds=time.monotonic() - t0), indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
