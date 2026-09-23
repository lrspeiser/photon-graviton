#!/usr/bin/env python3
"""Round 6: a star count of the Bullet Cluster from the Legacy Survey (DR10).

    python bullet_light_v6.py --output-dir ../run-bullet-light-v6

Round 5 left two predictions about where the Bullet Cluster's stars are:
1. **The main cluster's outskirts** hold more stars than a profile fitted to two central 100 kpc
   apertures: the measured galaxy speed needs an extra 6.2 x 10^12 Msun of stars spread to 3 Mpc
   (NFW scale 800 kpc).
2. **The subcluster** held about a third of the main cluster's matter before the collision, so
   about 7 x 10^12 Msun of its stars should still travel with it. Barrena et al.'s velocities show
   no such galaxies (bullet_members_v6.py), so any such stars would have to be spread out or
   diffuse.

Data: the DR10 Tractor catalogue (DECam g, r, i, z; photometric redshifts of Zhou et al.) in a
1.6 x 0.9 degree box around the Bullet, galaxies with z < 22.5, from the NOIRLab Astro Data Lab.
The query is in QUERY below, and the result is stored in ../data/ls_dr10_bullet.csv.gz.

Method.
* Cluster galaxies: z-band magnitude below 21.5 (absolute about -19.4) and photometric redshift
  within 0.04 of 0.296. The field level of the same selection is measured 4-6 Mpc out and
  subtracted.
* Light: z-band flux (rest frame about I). Stellar mass at M/L_I = 2, as Clowe et al. assumed.
  The comparison uses shapes (light beyond R relative to light inside 250 kpc), so the M/L,
  the magnitude limit and the photo-z completeness cancel.
* Main cluster: the eastern half, away from the subcluster (the western half is also reported).
* Subcluster: the excess light within 250 and 500 kpc of its brightest galaxy, after removing
  the field and the main cluster's own light at the same distance from the main cluster (taken
  from the eastern side).
"""
from __future__ import annotations
import argparse, csv, gzip, json, time
from pathlib import Path
import numpy as np
import bullet_v3 as B

QUERY = ("SELECT t.ra, t.dec, t.type, t.flux_g, t.flux_r, t.flux_i, t.flux_z, t.flux_w1, t.mw_transmission_g, t.mw_transmission_r, "
         "t.mw_transmission_i, t.mw_transmission_z, t.maskbits, t.shape_r, p.z_phot_mean, p.z_phot_std FROM ls_dr10.tractor AS t "
         "LEFT JOIN ls_dr10.photo_z AS p ON t.ls_id = p.ls_id WHERE t.ra BETWEEN 103.85 AND 105.45 AND t.dec BETWEEN -56.40 AND -55.50 "
         "AND t.type <> 'PSF' AND t.flux_z > 1.0")
Z_CL, DZ, ZLIM = 0.296, 0.04, 21.5
DM = 40.91                      # distance modulus at z = 0.296 (H0 = 70, Om = 0.3)
MSUN_Z = 4.50                   # the Sun's absolute AB magnitude in z
ML = 2.0                        # M/L_I, as Clowe et al.
ROUND5 = dict(main_inner=dict(kind='nfw', M=9.022e12, scale=466.5, rt=1500.0),
              main_outer=dict(kind='nfw', M=6.23e12, scale=800.0, rt=3000.0),
              sub_core=dict(kind='nfw', M=0.986e12, scale=12.1, rt=600.0))
LOST = {'1:3': 7.1e12, '1:8': 2.05e12}


def proj_mass(c, R):
    """Projected mass inside radius R (unit normalisation times M)."""
    prof = B.rho_nfw if c['kind'] == 'nfw' else B.rho_beta
    r = np.geomspace(0.5, c['rt'], 20000); dr = np.gradient(r)
    dm = c['M'] / B.mass3d(prof, (c['scale'],), c['rt']) * prof(r, c['scale'], c['rt']) * 4 * np.pi * r ** 2 * dr
    R = np.atleast_1d(R)[:, None]
    f = np.where(r[None, :] <= R, 1.0, 1 - np.sqrt(np.clip(1 - (R / r[None, :]) ** 2, 0, 1)))
    return (f * dm[None, :]).sum(1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    here = Path(__file__).resolve().parent
    rows = list(csv.DictReader(gzip.open(here.parent / 'data/ls_dr10_bullet.csv.gz', 'rt')))
    f = lambda k: np.array([float(r[k]) if r[k] not in ('', 'None') else np.nan for r in rows])
    ra, dec, zp = f('ra'), f('dec'), f('z_phot_mean')
    zmag = 22.5 - 2.5 * np.log10(f('flux_z') / f('mw_transmission_z'))
    mask = f('maskbits').astype(int)
    starmask = (mask & ((1 << 1) | (1 << 11))) != 0          # bright and medium star halos
    pos = B.positions()
    ref = B.radec_to_kpc(*B.POS_RAW['main_bcg'])
    x = -((ra * 3600 * np.cos(np.radians(B.DEC0))) - ref[0]) * B.SCALE      # west
    y = (-(np.abs(dec) * 3600) - ref[1]) * B.SCALE                          # north
    Rm = np.hypot(x - pos['main_bcg'][0], y - pos['main_bcg'][1]); Rs = np.hypot(x - pos['sub_bcg'][0], y - pos['sub_bcg'][1])
    Lz = 10 ** (-0.4 * (zmag - DM - MSUN_Z))                                # Lsun (z band, no k-correction)

    results = {}
    for dz in (0.04, 0.06):
        sel = (zmag < ZLIM) & (np.abs(zp - Z_CL) < dz) & ~starmask
        bg_ring = sel & (Rm > 4000) & (Rm < 6000)
        bg_density = Lz[bg_ring].sum() / (np.pi * (6000 ** 2 - 4000 ** 2))  # Lsun per kpc^2
        edges = np.array([0, 100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000])
        east, west = x < pos['main_bcg'][0], x >= pos['main_bcg'][0]
        prof = {}
        for side, half in (('east', east), ('west', west)):
            cum = []
            for Rr in edges[1:]:
                area = np.pi * Rr ** 2 / 2
                cum.append(2 * (Lz[sel & half & (Rm < Rr)].sum() - bg_density * area))    # full-circle equivalent
            prof[side] = np.array(cum)
        # shapes: light inside R relative to light inside 250 kpc (full circle there: no subcluster within 400 kpc)
        full250 = Lz[sel & (Rm < 250)].sum() - bg_density * np.pi * 250 ** 2
        shape_e = prof['east'] / full250
        mA = proj_mass(ROUND5['main_inner'], edges[1:])
        mL = mA + proj_mass(ROUND5['main_outer'], edges[1:])
        shapeA, shapeL = mA / mA[1], mL / mL[1]
        # subcluster: light within 250 / 500 kpc of its BCG, minus field and the main cluster's light at the same R_main (east side)
        sub = {}
        for Rap in (250.0, 500.0):
            inside = sel & (Rs < Rap)
            # main cluster's surface brightness at each R_main, from the east side, averaged over the aperture
            ann = np.arange(0, 3100, 100.0)
            sb_e = []
            for lo in ann[:-1]:
                ring = sel & east & (Rm >= lo) & (Rm < lo + 100)
                sb_e.append(Lz[ring].sum() / (np.pi * ((lo + 100) ** 2 - lo ** 2) / 2) - bg_density)
            sb_e = np.array(sb_e)
            g = np.random.default_rng(1).uniform(-Rap, Rap, (200000, 2)); g = g[np.hypot(*g.T) < Rap]
            rm_g = np.hypot(g[:, 0] + pos['sub_bcg'][0] - pos['main_bcg'][0], g[:, 1] + pos['sub_bcg'][1] - pos['main_bcg'][1])
            main_there = np.interp(rm_g, ann[:-1] + 50, sb_e).mean() * np.pi * Rap ** 2
            excess = Lz[inside].sum() - bg_density * np.pi * Rap ** 2 - main_there
            # model: the core, and the lost galaxies if they were still where they started (NFW scale 150 kpc)
            core_m = proj_mass(ROUND5['sub_core'], [Rap])[0]
            lost = {k: proj_mass(dict(kind='nfw', M=M, scale=150.0, rt=1000.0), [Rap])[0] for k, M in LOST.items()}
            sub[f'{int(Rap)}kpc'] = dict(excess_Lsun=float(excess), excess_Msun=float(ML * excess), main_light_removed_Lsun=float(main_there),
                                         model_core_Msun=float(core_m), model_lost_Msun={k: float(v) for k, v in lost.items()})
        results[f'dz_{dz}'] = dict(n_selected=int(sel.sum()), bg_Lsun_per_kpc2=float(bg_density),
                                   radii_kpc=edges[1:].tolist(), light_inside_east_Lsun=prof['east'].tolist(), light_inside_west_Lsun=prof['west'].tolist(),
                                   stars_inside_east_Msun=(ML * prof['east']).tolist(), shape_east=shape_e.tolist(),
                                   model_shape_round4=shapeA.tolist(), model_shape_round5=shapeL.tolist(),
                                   model_Msun_round4=mA.tolist(), model_Msun_round5=mL.tolist(), subcluster=sub)
        print(f"photo-z window +-{dz}: {sel.sum()} galaxies; field level {bg_density:.3e} Lsun/kpc^2", flush=True)
        print('  R (kpc)          ' + ' '.join(f'{r:7.0f}' for r in edges[1:]), flush=True)
        print('  stars inside R, east half x2 (1e12 Msun) ' + ' '.join(f'{v / 1e12:7.2f}' for v in ML * prof['east']), flush=True)
        print('  shape (east)     ' + ' '.join(f'{v:7.2f}' for v in shape_e), flush=True)
        print('  round-4 model    ' + ' '.join(f'{v:7.2f}' for v in shapeA), flush=True)
        print('  round-5 model    ' + ' '.join(f'{v:7.2f}' for v in shapeL), flush=True)
        for k, v in sub.items():
            print(f"  subcluster within {k}: excess stars {v['excess_Msun'] / 1e12:.2f} x 1e12 Msun (after removing {ML * v['main_light_removed_Lsun'] / 1e12:.2f} of main-cluster light);"
                  f" model core {v['model_core_Msun'] / 1e12:.2f}, lost galaxies if still in place: 1:3 {v['model_lost_Msun']['1:3'] / 1e12:.2f}, 1:8 {v['model_lost_Msun']['1:8'] / 1e12:.2f}", flush=True)

    (out / 'light_v6.json').write_text(json.dumps(dict(
        experiment='Bullet Cluster star count from Legacy Survey DR10 (round 6)', query=QUERY,
        selection=dict(z_cluster=Z_CL, zmag_limit=ZLIM, distance_modulus=DM, msun_z=MSUN_Z, mass_to_light=ML),
        results=results, seconds=time.monotonic() - t0), indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
