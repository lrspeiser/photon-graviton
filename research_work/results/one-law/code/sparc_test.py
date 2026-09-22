#!/usr/bin/env python3
"""JR-13b: the one law against 3152 raw SPARC rotation measurements plus the clusters.

JR-13 used R10 evaluated at the Einstein radius as the galaxy observable, so it
inherited R10. This removes that entirely: raw measured rotation speeds, raw
catalogue baryonic components, no companion model anywhere in the inputs.

    g_obs = v_obs^2 / r          from the measured rotation curve
    g_N   = (v_gas|v_gas| + Y v_disk^2 + Y_b v_bulge^2) / r     from the catalogue

and the law is

    g = tau(F) * [ g_N/2 + sqrt(g_N^2/4 + g_N a0) ],   tau = 1 + 1/(1+(F0/F)^m)

F = L_3.6 / (4 pi r^2) is the local starlight flux density, which varies WITH RADIUS
inside every galaxy. That makes the crossover testable within the rotation curves
themselves, not just between galaxies and clusters, and 149 galaxies spanning three
decades in surface brightness pin F0 in a way six lenses never could.

    python sparc_test.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, csv, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

import model as M

G = M.J.G
KMS2_PER_KPC_TO_SI = 1e6 / 3.0856775814913673e19
ROOT = Path(__file__).resolve().parents[4]


def law(gN, a0):
    return gN / 2.0 + np.sqrt(gN * gN / 4.0 + gN * a0)


def tau_of_flux(F, F0, m):
    return 1.0 + 1.0 / (1.0 + (F0 / F) ** m)


def load_sparc():
    """Raw measured rotation curves and catalogue baryonic components. No model."""
    gals = M.J.load_sparc()
    out = []
    for d in gals:
        r = d['r']
        gN = (d['gas_v2'] + d['star_v2']) / r          # catalogue baryonic force, km^2/s^2/kpc
        gobs = d['observed'] ** 2 / r
        ok = (gN > 0) & (gobs > 0) & np.isfinite(gN) & np.isfinite(gobs)
        if ok.sum() < 3:
            continue
        out.append(dict(name=d['name'], split=d['split'], r=r[ok], gN=gN[ok], gobs=gobs[ok],
                        v=d['observed'][ok], err=d['error'][ok], L36=d['L36'],
                        flux=d['L36'] / (4 * np.pi * r[ok] ** 2)))
    return out


def load_clusters():
    cc = {x['cluster']: x for x in csv.DictReader(
        (ROOT / 'companion_wave_test/cluster_comparison.csv').open())}
    rows = []
    for p in csv.DictReader((ROOT / 'companion_deposition_fit/cluster_predictions.csv').open()):
        c = cc[p['cluster']]
        M500 = float(c['M500_1e14Msun']) * 1e14
        Mstar = 0.02 * M500
        Mb = (float(c['fgas500']) + 0.02) * M500
        R = float(c['R500_Mpc']) * 1000.0
        Mdyn = float(p['required_extra_Msun']) + Mb
        rows.append(dict(name=p['cluster'], r=R, gN=G * Mb / R ** 2, gobs=G * Mdyn / R ** 2,
                         flux=(Mstar / 3.0) / (4 * np.pi * R ** 2)))
    return rows


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)
    t0 = time.monotonic()

    gals = load_sparc()
    clus = load_clusters()
    gN = np.concatenate([g['gN'] for g in gals])
    gobs = np.concatenate([g['gobs'] for g in gals])
    flux = np.concatenate([g['flux'] for g in gals])
    gid = np.concatenate([np.full(len(g['r']), i) for i, g in enumerate(gals)])
    n_gal_pts = len(gN)
    cN = np.array([c['gN'] for c in clus])
    cobs = np.array([c['gobs'] for c in clus])
    cflux = np.array([c['flux'] for c in clus])
    print(f'{len(gals)} SPARC galaxies, {n_gal_pts} rotation points; {len(clus)} clusters',
          flush=True)
    print(f'starlight flux spans {flux.min():.3g} to {flux.max():.3g} Lsun/kpc^2 in galaxies, '
          f'{cflux.min():.3g} to {cflux.max():.3g} in clusters', flush=True)

    # Per galaxy, weight all points equally so big curves do not dominate.
    w = 1.0 / np.bincount(gid)[gid] / len(gals)
    wc = np.full(len(clus), 1.0 / len(clus))

    def residuals(x):
        a0, F0, m = 10 ** x[0], 10 ** x[1], np.exp(x[2])
        pg = tau_of_flux(flux, F0, m) * law(gN, a0)
        pc = tau_of_flux(cflux, F0, m) * law(cN, a0)
        return np.log(gobs / pg), np.log(cobs / pc)

    def cost(x):
        try:
            rg, rc = residuals(x)
        except FloatingPointError:
            return 1e9
        if not (np.all(np.isfinite(rg)) and np.all(np.isfinite(rc))):
            return 1e9
        return float(np.sum(w * rg ** 2) + np.sum(wc * rc ** 2))

    best = min((minimize(cost, [np.log10(a), np.log10(f), np.log(mm)], method='Nelder-Mead',
                         options=dict(maxiter=30000, maxfev=30000, fatol=1e-13, xatol=1e-10))
                for a in (4e3, 1.5e4, 4e4) for f in (1e5, 1e6, 1e7) for mm in (0.3, 1.0)),
               key=lambda r: r.fun)
    a0, F0, m = 10 ** best.x[0], 10 ** best.x[1], float(np.exp(best.x[2]))
    rg, rc = residuals(best.x)
    tg = tau_of_flux(flux, F0, m)
    tc = tau_of_flux(cflux, F0, m)

    # Control: the same law with tau frozen at 1 everywhere (no stress state at all).
    def cost_notau(x):
        a0_ = 10 ** x[0]
        r1 = np.log(gobs / law(gN, a0_))
        r2 = np.log(cobs / law(cN, a0_))
        return float(np.sum(w * r1 ** 2) + np.sum(wc * r2 ** 2))
    b2 = min((minimize(cost_notau, [np.log10(a)], method='Nelder-Mead',
                       options=dict(maxiter=8000, fatol=1e-13, xatol=1e-10))
              for a in (4e3, 1.5e4, 4e4)), key=lambda r: r.fun)
    a0_nt = 10 ** b2.x[0]
    rg_nt = np.log(gobs / law(gN, a0_nt))
    rc_nt = np.log(cobs / law(cN, a0_nt))

    per_gal = []
    for i, g in enumerate(gals):
        sel = gid == i
        vpred = np.sqrt(np.maximum(g['r'] * tau_of_flux(g['flux'], F0, m) * law(g['gN'], a0), 0))
        per_gal.append(dict(name=g['name'], split=g['split'], n=int(sel.sum()),
                            log_scatter=float(np.std(rg[sel], ddof=1)) if sel.sum() > 1 else 0.0,
                            mean_log_ratio=float(np.mean(rg[sel])),
                            v_RMSE_kms=float(np.sqrt(np.mean((vpred - g['v']) ** 2))),
                            v_fractional_RMS=float(np.sqrt(np.mean(((vpred - g['v']) / g['v']) ** 2))),
                            tau_min=float(tau_of_flux(g['flux'], F0, m).min()),
                            tau_max=float(tau_of_flux(g['flux'], F0, m).max())))
    rmse = np.array([p['v_RMSE_kms'] for p in per_gal])

    payload = dict(
        experiment='JR-13b',
        scope='The one law against raw SPARC rotation measurements and the X-COP clusters. '
              'No companion model, no R10, appears anywhere in the inputs.',
        law='g = tau(F) * [ g_N/2 + sqrt(g_N^2/4 + g_N a0) ], tau = 1 + 1/(1+(F0/F)^m)',
        constants=dict(a0_kms2_per_kpc=float(a0), a0_SI=float(a0 * KMS2_PER_KPC_TO_SI),
                       F0_Lsun_per_kpc2=float(F0), sharpness=m, n_free_numbers=3),
        data=dict(n_galaxies=len(gals), n_rotation_points=int(n_gal_pts), n_clusters=len(clus),
                  galaxy_flux_range=[float(flux.min()), float(flux.max())],
                  cluster_flux_range=[float(cflux.min()), float(cflux.max())]),
        fit=dict(galaxy_log_scatter=float(np.std(rg, ddof=1)),
                 cluster_log_scatter=float(np.std(rc, ddof=1)),
                 galaxy_mean_log_ratio=float(np.mean(rg)),
                 cluster_mean_log_ratio=float(np.mean(rc)),
                 mean_galaxy_v_RMSE_kms=float(rmse.mean()),
                 median_galaxy_v_RMSE_kms=float(np.median(rmse)),
                 galaxies_within_20pct=int(sum(p['v_fractional_RMS'] < 0.20 for p in per_gal)),
                 galaxies_within_10pct=int(sum(p['v_fractional_RMS'] < 0.10 for p in per_gal)),
                 tau_range_galaxies=[float(tg.min()), float(tg.max())],
                 tau_range_clusters=[float(tc.min()), float(tc.max())]),
        control_tau_frozen_at_one=dict(
            a0_SI=float(a0_nt * KMS2_PER_KPC_TO_SI),
            galaxy_log_scatter=float(np.std(rg_nt, ddof=1)),
            cluster_log_scatter=float(np.std(rc_nt, ddof=1)),
            cluster_mean_log_ratio=float(np.mean(rc_nt)),
            galaxy_mean_log_ratio=float(np.mean(rg_nt)),
            note='One constant, no stress state. The offset between the two classes is what '
                 'tau has to explain.'),
        per_galaxy=per_gal,
        clusters=[dict(name=c['name'], tau=float(t), ratio=float(np.exp(r)))
                  for c, t, r in zip(clus, tc, rc)],
        seconds=time.monotonic() - t0,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'sparc-one-law.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    f = payload['fit']
    print(f"\na0 = {a0 * KMS2_PER_KPC_TO_SI:.3e} m/s^2   F0 = {F0:.3g} Lsun/kpc^2   m = {m:.3f}")
    print(f"  galaxies: log scatter {f['galaxy_log_scatter']:.1%}, mean v RMSE "
          f"{f['mean_galaxy_v_RMSE_kms']:.2f} km/s, {f['galaxies_within_20pct']}/{len(gals)} "
          f"within 20%, {f['galaxies_within_10pct']} within 10%")
    print(f"  clusters: log scatter {f['cluster_log_scatter']:.1%}, mean log ratio "
          f"{f['cluster_mean_log_ratio']:+.3f}")
    print(f"  tau spans {f['tau_range_galaxies'][0]:.3f}-{f['tau_range_galaxies'][1]:.3f} in "
          f"galaxies, {f['tau_range_clusters'][0]:.3f}-{f['tau_range_clusters'][1]:.3f} in clusters")
    c = payload['control_tau_frozen_at_one']
    print(f"\ncontrol, tau frozen at 1 (one constant, no stress state):")
    print(f"  a0 = {c['a0_SI']:.3e}   galaxy scatter {c['galaxy_log_scatter']:.1%}, "
          f"cluster scatter {c['cluster_log_scatter']:.1%}")
    print(f"  offset between classes: galaxies {c['galaxy_mean_log_ratio']:+.3f}, "
          f"clusters {c['cluster_mean_log_ratio']:+.3f} "
          f"-> factor {np.exp(c['galaxy_mean_log_ratio']-c['cluster_mean_log_ratio']):.2f}")
    print('WROTE', out / 'sparc-one-law.json')


if __name__ == '__main__':
    main()
