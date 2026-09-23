#!/usr/bin/env python3
"""KiDS-1000: why do ellipticals bend light more than spirals of the same stellar mass?

    python kids_haloes.py --output-dir ../run-kids

Brouwer et al. 2021 (A&A 650, A113) measure, for isolated lenses of equal stellar mass,
at least a factor 1.5 (~0.2 dex) higher lensing acceleration around early types, at
>= 5.7 sigma. g_bar there counts stars + cold gas only.

Model (both types have the same stars + cold gas, M_gal, so the same g_bar):
  spiral     : cold stars and gas, no hot halo
  elliptical : stars with random speed sigma_*, plus a hot corona of mass f_c M_* at
               temperature T_c (beta-model, beta = 0.5, core 5 kpc, out to 300 kpc)
Offset = log10 g_elliptical - log10 g_spiral at the same radius.
Ours: g = g_N + exp(-g_N/g_d) sqrt(a (g_N + S_hot)), with the corona's heat counted.
MOND (context): the corona counts only as extra mass.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import law as L

G, K = L.G, L.KMS2_PER_KPC
A = 6.5797e-11 / K; GD = 4.281 * A; U = 874.1
A_MOND = 8.423e-11 / K                     # MOND simple, fitted on the same SPARC data
MU, MP_KEV = 0.6, 938272.088
Mstar = 10 ** 10.8; Mgal = 1.15 * Mstar    # cold gas ~15% of stars for both (enters g_bar)
RADII = np.array([30., 100., 300., 1000.])


def corona(fc, s):
    rho = (1 + (s / 5.0) ** 2) ** -0.75 * (s <= 300)
    dm = 4 * np.pi * s ** 2 * rho * np.gradient(s)
    return dm * fc * Mstar / dm.sum()


def offsets(fc, Tc_keV, sigma_star=200.):
    s = np.geomspace(0.5, 3000, 3000)
    dmc = corona(fc, s) if fc > 0 else np.zeros_like(s)
    sig_c = np.sqrt(Tc_keV / (MU * MP_KEV)) * L.C_KMS
    k_c = L.heat_weight(sig_c, U); k_s = L.heat_weight(sigma_star, U)
    out = {}
    for R in RADII:
        gNs = G * Mgal / R ** 2
        gNc = G * np.sum(dmc[s < R]) / R ** 2
        Sc = G * (L.shell_weights([R], s) @ dmc)[0] / R ** 2
        g_sp = L.total(gNs, 0., A, 4.281)
        g_el = L.total(gNs + gNc, k_s * G * Mstar / R ** 2 + k_c * Sc, A, 4.281)
        m_sp = L.mond_simple(gNs, A_MOND); m_el = L.mond_simple(gNs + gNc, A_MOND)
        out[int(R)] = dict(ours=float(np.log10(g_el / g_sp)), mond=float(np.log10(m_el / m_sp)))
    return out, float(sig_c), float(k_c)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True)
    grid = {}
    print('offset (dex) at 30 / 100 / 300 / 1000 kpc;  observed >= 0.2')
    for Tc in (0.3, 0.6, 1.0):
        for fc in (0.0, 0.25, 0.5, 1.0):
            o, sc, kc = offsets(fc, Tc)
            grid[f'T{Tc}_f{fc}'] = dict(T_keV=Tc, f_corona=fc, sigma_corona=sc, k_corona=kc, offsets=o)
            print(f'  T_c={Tc:.1f} keV (sigma {sc:3.0f}, k {kc:.2f})  f_c={fc:4.2f}   ours: ' +
                  ' '.join(f"{o[r]['ours']:+.3f}" for r in o) + '   MOND: ' + ' '.join(f"{o[r]['mond']:+.3f}" for r in o))
    # corona mass needed for 0.2 dex at 100 kpc
    need = {}
    for Tc in (0.3, 0.6, 1.0):
        f_ours = brentq(lambda f: offsets(f, Tc)[0][100]['ours'] - 0.2, 0.0, 20.)
        f_mond = brentq(lambda f: offsets(f, Tc)[0][100]['mond'] - 0.2, 0.0, 20.)
        need[Tc] = dict(ours=f_ours, mond=f_mond)
        print(f'  corona mass needed for 0.2 dex at 100 kpc, T_c = {Tc} keV:  ours {f_ours:.2f} M*   MOND {f_mond:.2f} M*')
    (out / 'kids.json').write_text(json.dumps(dict(grid=grid, corona_mass_needed_for_0p2_dex=need), indent=2) + '\n')


if __name__ == '__main__':
    main()
