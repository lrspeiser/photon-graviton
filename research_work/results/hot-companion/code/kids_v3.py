#!/usr/bin/env python3
"""KiDS-1000 early- versus late-type lensing under round-3 hot-companion gravity.

    python kids_v3.py --output-dir ../run-kids-v3

Brouwer et al. 2021 (A&A 650, A113, Sect. 5.4): at equal stellar mass, early types show a
higher lensing acceleration than late types; the mean ratio log10(g_obs,E / g_obs,L) is
0.17 dex (split by Sersic index) and 0.27 dex (split by colour), at >= 5.7 sigma. g_bar
counts stars and cold gas only.

Round 3 changes the story from round 2: gas no longer feeds the heat term (it collides), so a
hot gas halo counts only as extra mass. The heat now comes from the early types' own stars,
which move randomly and never collide. Late types' disk stars move in near-circular orbits
with small random speeds.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import law as L

G = L.G
Mstar = 10 ** 10.8
RADII = np.array([30., 100., 300., 1000.])


def offsets(consts, sig_e, sig_l, f_halo, gas_frac_late=0.15):
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    s = np.geomspace(0.5, 3000, 3000)
    rho = (1 + (s / 5.0) ** 2) ** -0.75 * (s <= 300)
    dmh = 4 * np.pi * s ** 2 * rho * np.gradient(s); dmh *= f_halo * Mstar / max(dmh.sum(), 1e-30)
    ke, kl = L.heat_weight(sig_e, u), L.heat_weight(sig_l, u)
    out = {}
    for R in RADII:
        # equal g_bar: stars (+ cold gas for the late type) are what the survey counts
        gbar = G * Mstar * (1 + gas_frac_late) / R ** 2
        gN_e = gbar + G * np.sum(dmh[s < R]) / R ** 2          # early type: halo is extra ordinary mass
        g_e = L.total(gN_e, ke * G * Mstar / R ** 2, a, lam)    # stars far away act as a point for S
        g_l = L.total(gbar, kl * G * Mstar / R ** 2, a, lam)
        out[int(R)] = float(np.log10(g_e / g_l))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True)
    consts = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    grid = []
    print('log10(g_E / g_L) at 30 / 100 / 300 / 1000 kpc;  observed mean 0.17 (Sersic split) and 0.27 (colour split)')
    for sig_e in (150., 175., 200.):
        for f_halo in (0.0, 0.3, 1.0):
            o = offsets(consts, sig_e, 30., f_halo)
            grid.append(dict(sigma_early=sig_e, sigma_late_disk=30., hot_halo_mass_over_Mstar=f_halo, offsets=o))
            print(f'  early-type stars {sig_e:.0f} km/s, hot halo {f_halo:.1f} M*: ' + '  '.join(f'{o[r]:+.3f}' for r in o))
    (out / 'kids.json').write_text(json.dumps(dict(constants=consts, grid=grid,
        observed=dict(sersic_split=0.17, colour_split=0.27, source='Brouwer et al. 2021, Sect. 5.4')), indent=2) + '\n')


if __name__ == '__main__':
    main()
