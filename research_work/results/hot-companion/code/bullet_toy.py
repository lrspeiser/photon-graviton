#!/usr/bin/env python3
"""Bullet Cluster toy: ballistic companions after a cluster collision.

    python bullet_toy.py --output-dir ../run-bullet

Idea (proposed here): companion energy leaves each piece of matter carrying that matter's
velocity at the moment of emission, then spreads at speed u. In steady motion the
companion simply travels with its source. When the gas of two colliding clusters is
stopped by ram pressure, the companion it emitted BEFORE the collision flies on with the
original velocity -- along with the galaxies, which pass through each other. Only
companion emitted AFTER the collision is centred on the stopped gas, and it has spread
only to u * t_since.

Geometry is a representative toy, NOT a fit to the Bullet Cluster maps:
  150 Myr after core passage (relative speed 4700 km/s): galaxy groups 720 kpc apart,
  each cluster's gas lagging its galaxies. Gas outweighs stars ~7:1, as observed.
Lensing (light feels the same pull, eta = 1): surface density of ordinary matter plus the
effective ('phantom') density of the extra pull, projected along the line of sight.
Compared with the same law without memory (companion follows the matter as it is now).
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import field_equation as FE

G = FE.G; A = FE.A; GD = FE.GD; U = FE.U
MYR_KPC_PER_KMS = 1.0227e-3                 # kpc per Myr per (km/s)
MU, MP_KEV, C = 0.6, 938272.088, 299792.458

COMPONENTS = [
    # name,           kind,  mass,   core,  x_now, x_ghost, pre-T or sigma, post-T or sigma
    ('main gas',      'gas', 1.0e14, 150., -80.,  -200.,  8.0, 12.0),
    ('main galaxies', 'gal', 1.3e13, 150., -200., -200., 1100., 1100.),
    ('sub gas',       'gas', 2.0e13,  70., 320.,   520.,  5.0, 15.0),
    ('sub galaxies',  'gal', 3.0e12, 100., 520.,   520.,  800., 800.),
]
SMOOTH_KPC = 60.                            # typical resolution of a lensing mass map


def sigma_of(kind, v):
    return np.sqrt(v / (MU * MP_KEV)) * C if kind == 'gas' else v


def profile(X, Y, Z, kind, M, core, x0):
    """Gas: isothermal beta-model (beta = 2/3). Galaxies: the same form with a galaxy-count
    core -- galaxy groups in clusters are extended, not point-like."""
    r = np.sqrt((X - x0) ** 2 + Y ** 2 + Z ** 2)
    rho = (1 + (r / core) ** 2) ** -1.0 * (r < 1500)
    return rho * M / (rho.sum())                            # mass per cell


class MaskedConv:
    """Radial-kernel FFT convolutions, full / only inside d<R / only outside d>R."""
    def __init__(self, shape, dx):
        self.shape = shape; self.dx = dx
        ax = [(np.arange(2 * n) - (np.arange(2 * n) >= n) * 2 * n) * dx for n in shape]
        X, Y, Z = np.meshgrid(*ax, indexing='ij')
        self.d = np.sqrt(X ** 2 + Y ** 2 + Z ** 2); self.dvec = (X, Y, Z)
        self.cache = {}

    def kernels(self, R, part):
        key = (R, part)
        if key not in self.cache:
            d = np.maximum(self.d, 0.5 * self.dx)
            m = np.ones_like(d) if part == 'full' else ((d <= R) if part == 'in' else (d > R))
            vec = [np.fft.rfftn((c / d ** 3 * m).astype(np.float32)) for c in self.dvec]
            sca = np.fft.rfftn((1 / d ** 2 * m).astype(np.float32))
            self.cache[key] = (vec, sca)
        return self.cache[key]

    def conv(self, f, K):
        s = [2 * n for n in self.shape]
        pad = np.zeros(s, np.float32); pad[:self.shape[0], :self.shape[1], :self.shape[2]] = f
        F = np.fft.rfftn(pad)
        return np.fft.irfftn(F * K, s=s, axes=(0, 1, 2))[:self.shape[0], :self.shape[1], :self.shape[2]]

    def fields(self, mass, R, part):
        vec, sca = self.kernels(R, part)
        # pull toward the mass: g = -G sum m (x - x')/|x - x'|^3
        return np.array([-G * self.conv(mass, k) for k in vec]), G * self.conv(mass, sca)


def run(t_myr, memory, grid):
    X, Y, Z, dx, mc = grid
    R = U * t_myr * MYR_KPC_PER_KMS
    rho_now = np.zeros(X.shape, np.float32)
    F = np.zeros((3,) + X.shape); S = np.zeros(X.shape)
    for name, kind, M, sc, xnow, xghost, pre, post in COMPONENTS:
        m_now = profile(X, Y, Z, kind, M, sc, xnow); rho_now += m_now
        k_pre = 3 * sigma_of(kind, pre) ** 2 / U ** 2; k_post = 3 * sigma_of(kind, post) ** 2 / U ** 2
        if memory and xghost != xnow:
            m_ghost = profile(X, Y, Z, kind, M, sc, xghost)
            f_out, s_out = mc.fields(m_ghost, R, 'out')      # emitted before: centred on the ghost
            f_in, s_in = mc.fields(m_now, R, 'in')           # emitted after: centred on the gas now
            F += f_out + f_in; S += k_pre * s_out + k_post * s_in
        else:
            f, s = mc.fields(m_now, R, 'full')
            F += f; S += (k_post if memory else k_pre) * s if kind == 'gas' else k_pre * s
    gN, _ = mc.fields(rho_now, R, 'full')
    q = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fm = np.sqrt(np.sum(F ** 2, axis=0)) + 1e-30
    gc = (np.exp(-q / GD) * np.sqrt(A * (Fm + S)) / Fm) * F
    rho_ph = -FE.div(gc, dx) / (4 * np.pi * G) * dx ** 3        # phantom mass per cell
    from scipy.ndimage import gaussian_filter
    Sigma_b = gaussian_filter(rho_now.sum(axis=2), SMOOTH_KPC / dx)
    Sigma_t = Sigma_b + gaussian_filter(rho_ph.sum(axis=2), SMOOTH_KPC / dx)
    j0 = X.shape[1] // 2
    xs = X[:, 0, 0]
    prof_b = Sigma_b[:, j0 - 1:j0 + 1].mean(1); prof_t = Sigma_t[:, j0 - 1:j0 + 1].mean(1)
    from scipy.signal import argrelmax
    def maxima(prof):
        return [float(xs[i]) for i in argrelmax(prof, order=3)[0]]
    ph = gaussian_filter(rho_ph.sum(axis=2), SMOOTH_KPC / dx)[:, j0 - 1:j0 + 1].mean(1)
    def centroid(w, lo, hi):
        sel = (xs > lo) & (xs < hi); return float(np.sum(xs[sel] * w[sel]) / np.sum(w[sel]))
    return dict(t_myr=t_myr, spread_kpc=R, memory=memory,
                lens_maxima=maxima(prof_t), baryon_maxima=maxima(prof_b),
                extra_mass_centre_main=centroid(ph, -500, 60), extra_mass_centre_sub=centroid(ph, 260, 800),
                baryon_centre_main=centroid(prof_b, -500, 60), baryon_centre_sub=centroid(prof_b, 260, 800),
                phantom_fraction=float(rho_ph.sum() / (rho_ph.sum() + rho_now.sum())),
                profile_x=xs.tolist(), profile_total=prof_t.tolist(), profile_baryons=prof_b.tolist())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    dx = 20.; shape = (150, 80, 80)
    xs = (np.arange(shape[0]) - shape[0] / 2 + 0.5) * dx + 160; ys = (np.arange(shape[1]) - shape[1] / 2 + 0.5) * dx
    X, Y, Z = np.meshgrid(xs, ys, ys, indexing='ij')
    mc = MaskedConv(shape, dx); grid = (X, Y, Z, dx, mc)
    rows = [run(150, False, grid)] + [run(t, True, grid) for t in (35, 70, 105, 150, 300, 600, 1200)]
    (out / 'bullet.json').write_text(json.dumps(rows, indent=1) + '\n')
    print('galaxies at x = -200 (main) and +520 (sub); gas at -80 (main) and +320 (sub)   [kpc]')
    print('                                   lensing maxima            ordinary-matter maxima   extra-mass centre main | sub')
    for r in rows:
        lab = f"memory, {r['t_myr']:4.0f} Myr (spread {r['spread_kpc']:4.0f} kpc)" if r['memory'] else 'no memory (follows matter now)'
        print(f"  {lab:36s} {str([round(v) for v in r['lens_maxima']]):24s} {str([round(v) for v in r['baryon_maxima']]):22s}"
              f"  {r['extra_mass_centre_main']:+6.0f} | {r['extra_mass_centre_sub']:+6.0f}")
    print(f'{time.monotonic() - t0:.0f} s')


if __name__ == '__main__':
    main()
