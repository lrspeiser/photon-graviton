#!/usr/bin/env python3
"""Round 8: three more colliding clusters under our law -- MACS J0025.4-1222, Abell 520, El Gordo.

    python collisions_v8.py --output-dir ../run-collisions-v8

The machinery is the Bullet's (bullet_v4.py, rounds 4-5), with nothing refitted:
* today's gas and stars set the ordinary pull and the release factor;
* the companion emitted before the collision rides with the galaxies: each pre-collision
  subcluster is a settled cluster whose gas was centred on its galaxies (its 'ghost');
* inside a fresh sphere of radius u t around the gas where it is now (t = time since the gas was
  stopped), the companion has been rebuilt around today's gas;
* the stars' random speeds come from the Jeans equation in the pull our law makes from each
  pre-collision subcluster's own gas and stars, and set the heat weight k = 3 sigma^2 / u^2.
Lensing is the projected 'as if' density -div h / 4 pi G; masses in apertures are compared
directly (they do not depend on the source redshift). Distances and masses follow the papers'
convention (flat LCDM, H0 = 70, Om = 0.3), as for the Bullet.

Inputs (all published; see data/collisions_v8.json, written by this script, for every number and
its source):
* MACS J0025.4-1222 (z = 0.586; Bradac et al. 2008): two near-equal subclusters 540 kpc apart that
  collided in the plane of the sky a few 10^8 yr ago at about 2,000 km/s; one gas cloud between them
  (King core 400 kpc; 3.6e13 Msun inside 500 kpc); stars 2.7 and 1.9 x 10^12 inside 300 kpc of the
  brightest galaxies; galaxy speed spread 835 km/s. Lensing: 2.5 and 2.6 x 10^14 inside 300 kpc,
  peaks with the galaxies (> 4 sigma from the gas).
* Abell 520 (z = 0.201; Mahdavi et al. 2007; Jee et al. 2012, 2014; Clowe et al. 2012; Wang et al.
  2016; Girardi et al. 2008): a 'train wreck' along NE-SW, about 0.5-1 Gyr after the main passage;
  galaxies in clumps P1, P2 (NE), P4 (SW), P5 (E), P6; the gas stuck in the middle (P3), where two
  teams disagree on the lensing mass: 3.4-4.0 (Jee) or 2.3 (Clowe) x 10^13 inside 150 kpc.
  Gas and light per clump from Clowe et al. 2012; galaxy speed spread 1,066 km/s.
* El Gordo (z = 0.870; Menanteau et al. 2012; Jee et al. 2014; Ng et al. 2015; Kim et al. 2021):
  NW and SE subclusters 700 kpc apart, 0.46 (outgoing) or 0.91 Gyr (returning) after pericentre
  at 2,400 km/s; gas 2.2e14 inside R500 (14.5 keV), with the SE cool core 8 arcsec beyond the SE
  lensing peak; stars 7.5 (NW) and 5.6 (SE) x 10^12; galaxy speeds 1,290 and 1,089 km/s; lensing
  M200c 9.9 and 6.5 x 10^14 (Kim et al. 2021, NFW fits consistent with their model-free masses).
"""
from __future__ import annotations
import argparse, copy, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import bullet_v3 as B
import bullet_v4 as V
import bullet_main_v5 as BM
import collisions_v4 as C4

G = B.G
HERE = Path(__file__).resolve().parent
U0 = 197.41001228101223


# ----------------------------------------------------------------------------- geometry
def kpc_per_arcsec(z, h=0.7, om=0.3):
    c = 299792.458
    chi = quad(lambda x: 1 / np.sqrt(om * (1 + x) ** 3 + 1 - om), 0, z)[0] * c / (100 * h)     # Mpc
    return chi / (1 + z) * 1e3 / 206264.80624709636


def sky(ra, dec):
    h, m, s = [float(v) for v in ra.split(':')]
    sgn = -1.0 if dec.strip().startswith('-') else 1.0
    d, dm, ds = [abs(float(v)) for v in dec.split(':')]
    return (h + m / 60 + s / 3600) * 15.0, sgn * (d + dm / 60 + ds / 3600)


def to_kpc(coords, ref, scale):
    """x to the west, y to the north (kpc), as in bullet_v3."""
    r0, d0 = sky(*ref)
    out = {}
    for k, (ra, dec) in coords.items():
        r, d = sky(ra, dec)
        out[k] = np.array([-(r - r0) * np.cos(np.radians(d0)) * 3600 * scale, (d - d0) * 3600 * scale])
    return out


def beta(M, rc, rt, centre):
    return dict(kind='beta', M=float(M), scale=float(rc), rt=float(rt), centre=centre)


def nfw(M, rs, rt, centre):
    return dict(kind='nfw', M=float(M), scale=float(rs), rt=float(rt), centre=centre)


def projected(c, R):
    return BM.projected_mass(c, R)


def sphere(c, R):
    return float(BM.cum_mass(c, np.array([R]))[0])


def with_projected(c, R, target):
    c = dict(c); c['M'] = 1.0; c['M'] = target / projected(c, R); return c


def column(c, centre_xy, at_xy, R=150.0, n=64):
    """Projected mass of a component inside a circle of radius R centred at at_xy."""
    d = float(np.linalg.norm(np.asarray(at_xy) - np.asarray(centre_xy)))
    prof = B.rho_beta if c['kind'] == 'beta' else B.rho_nfw
    unit = B.aperture(prof, (c['scale'],), c['rt'], d, R)
    return c['M'] / B.mass3d(prof, (c['scale'],), c['rt']) * unit


# ----------------------------------------------------------------------------- the clusters
def macs0025():
    s = 6.61                                           # kpc/arcsec (Bradac et al. 2008, Sect. 1)
    am = 60 * s                                        # kpc per arcmin
    pos = dict(gas_peak=np.array([0.0, 0.0]),
               se_gal=np.array([-0.79, -0.50]) * am, nw_gal=np.array([0.39, 0.18]) * am,
               se_lens=np.array([-0.70, -0.43]) * am, nw_lens=np.array([0.46, 0.20]) * am)
    # gas: one beta model about the X-ray peak, ending at 1.5 Mpc (the paper's 1.1-2 Mpc), with the core
    # that gives both published gas masses: 3.6e13 inside a 500 kpc sphere and 5.5e13 projected inside
    # 500 kpc (the quoted King core, 400 kpc, gives 7.5e13 projected)
    def icm_for(rc):
        c = beta(1.0, rc, 1500.0, 'gas_peak'); c['M'] = 3.6e13 / sphere(c, 500.0); return c
    rc = brentq(lambda rc: projected(icm_for(rc), 500.0) - 5.5e13, 20.0, 400.0)
    icm = icm_for(rc)
    # stars: NFW-shaped galaxy populations normalised to the published masses inside 300 kpc
    st_se = with_projected(nfw(1.0, 250.0, 1500.0, 'se_gal'), 300.0, 2.7e12)
    st_nw = with_projected(nfw(1.0, 250.0, 1500.0, 'nw_gal'), 300.0, 1.9e12)
    gas_total = icm['M']
    subs = [dict(name='SE', galaxies='se_gal', stars=['st_se'], gas_share=2.5 / 5.1, rc_pre=150.0, rt_pre=1500.0),
            dict(name='NW', galaxies='nw_gal', stars=['st_nw'], gas_share=2.6 / 5.1, rc_pre=150.0, rt_pre=1500.0)]
    return dict(name='MACS J0025.4-1222', z=0.586, kpc_per_arcsec=s, pos=pos,
                current=dict(gas_icm=icm, st_se=st_se, st_nw=st_nw), gas_total=gas_total, subs=subs,
                t_gyr=(0.26, 0.15, 0.40), centre=(-80.0, -60.0), dx=15.0,
                checks=dict(gas_projected_500=[5.5e13, 0.6e13], stars_projected_500=[5.0e12, 1.0e12]),
                sources=dict(geometry='Bradac et al. 2008 (arXiv:0806.2320) Table 2; 6.61 kpc/arcsec',
                             gas='Sect. 5-6: King core ~60 arcsec (400 kpc); 3.6 +- 0.4e13 in a 500 kpc sphere; 5.5 +- 0.6e13 projected',
                             stars='Table 3: 0.027 +- 0.008 and 0.019 +- 0.006 x 1e14 inside 300 kpc (M/L_K 0.74)',
                             lensing='Table 3: 2.5 +1.0/-1.7 and 2.6 +0.5/-1.4 x 1e14 inside 300 kpc of BCG1 and BCG3; 6.2 +1.2/-4.0 inside 500 kpc of the gas peak',
                             speeds='Sect. 3: 835 +58/-59 km/s (108 galaxies); merger ~2000 km/s in the sky plane, closest approach a few 1e8 yr ago'))


def abell520(ml_z=2.0):
    z = 0.201; s = kpc_per_arcsec(z)
    # Clowe et al. 2012, Table 1 (r < 150 kpc): positions, z-band light (1e11 Lsun), gas column (1e13)
    clowe = dict(P1=('04:54:19.60', '+02:57:49.09', 2.43, 0.25), P2=('04:54:14.84', '+02:57:06.25', 4.16, 0.40),
                 P3=('04:54:11.25', '+02:55:37.28', 1.38, 0.69), P4=('04:54:04.57', '+02:53:58.60', 3.11, 0.50),
                 P5=('04:54:17.11', '+02:55:30.09', 2.66, 0.44), P6=('04:54:09.61', '+02:53:55.90', 1.15, 0.65))
    pos = to_kpc({k: v[:2] for k, v in clowe.items()}, clowe['P3'][:2], s)
    pos.update(to_kpc({'P3_jee': ('04:54:11.07', '+02:55:35.3'), 'P3p_jee': ('04:54:07.51', '+02:54:41.3')}, clowe['P3'][:2], s))
    # gas: one beta model about P3 fitted to the six gas columns
    from scipy.optimize import least_squares
    cols = np.array([clowe[k][3] for k in clowe]) * 1e13
    def res(q):
        c = beta(np.exp(q[0]), np.exp(q[1]), 1500.0, 'P3')
        return np.log([column(c, pos['P3'], pos[k]) for k in clowe]) - np.log(cols)
    q = least_squares(res, [np.log(1e14), np.log(250.0)], bounds=([np.log(1e12), np.log(30.0)], [np.log(1e16), np.log(1500.0)])).x
    icm = beta(np.exp(q[0]), np.exp(q[1]), 1500.0, 'P3')
    # stars: each clump's light (z band) at M/L_z = 2 inside 150 kpc, NFW-shaped (scale 100 kpc, cut at 1 Mpc)
    current = dict(gas_icm=icm)
    for k, v in clowe.items():
        current[f'st_{k.lower()}'] = with_projected(nfw(1.0, 100.0, 1000.0, k), 150.0, ml_z * v[2] * 1e11)
    Ltot = sum(v[2] for v in clowe.values())
    subs = [dict(name=k, galaxies=k, stars=[f'st_{k.lower()}'], gas_share=v[2] / Ltot, rc_pre=120.0, rt_pre=1500.0)
            for k, v in clowe.items()]
    return dict(name='Abell 520', z=z, kpc_per_arcsec=s, pos=pos, current=current, gas_total=icm['M'], subs=subs,
                t_gyr=(0.5, 0.3, 1.0), centre=(0.0, 0.0), dx=15.0, mass_to_light_z=ml_z,
                gas_fit=dict(M=icm['M'], rc=icm['scale'], columns_model=[column(icm, pos['P3'], pos[k]) for k in clowe],
                             columns_obs=cols.tolist()),
                sources=dict(geometry='Clowe et al. 2012 (arXiv:1209.2143) Table 1; Jee et al. 2014 (arXiv:1401.3356) Table 1',
                             gas='Clowe et al. 2012 Table 1 gas columns inside 150 kpc (fitted by one beta model about P3)',
                             stars='Clowe et al. 2012 Table 1 z-band light inside 150 kpc, M/L_z = 2 (as M/L_I = 2 for the Bullet)',
                             lensing='Jee et al. 2014 Table 1 and Clowe et al. 2012 Table 1 (masses inside 150 kpc)',
                             timing='Markevitch et al. 2005 (shock 2,300 km/s); Mahdavi et al. 2007 (about 1 Gyr)',
                             speeds='Girardi et al. 2008: 1066 +67/-61 km/s (167 members)'))


def el_gordo():
    z = 0.870; s = kpc_per_arcsec(z)
    coords = {'nw_gal': ('01:02:51.23', '-49:15:02.56'), 'se_gal': ('01:02:56.95', '-49:16:21.86'),
              'com': ('01:02:53.49', '-49:15:33.96')}
    pos = to_kpc(coords, coords['com'], s)
    axis = (pos['se_gal'] - pos['nw_gal']) / np.linalg.norm(pos['se_gal'] - pos['nw_gal'])
    pos['cool_core'] = pos['se_gal'] + axis * 8.0 * s              # 8 arcsec beyond the SE lensing peak
    pos['xray_c'] = 0.5 * (pos['cool_core'] + pos['com'])          # the bulk of the X-ray gas: SE half, wake toward NW
    icm = beta(1.0, 250.0, 2500.0, 'xray_c'); icm['M'] = 2.1e14 / sphere(icm, 1300.0)
    core = beta(1.0e13, 30.0, 300.0, 'cool_core')
    st_nw = nfw(7.5e12, 350.0, 2100.0, 'nw_gal'); st_se = nfw(5.6e12, 350.0, 2100.0, 'se_gal')
    gas_total = icm['M'] + core['M']
    subs = [dict(name='NW', galaxies='nw_gal', stars=['st_nw'], gas_share=9.9 / 16.4, rc_pre=200.0, rt_pre=2500.0),
            dict(name='SE', galaxies='se_gal', stars=['st_se'], gas_share=6.5 / 16.4, rc_pre=100.0, rt_pre=2500.0)]
    return dict(name='El Gordo (ACT-CL J0102-4915)', z=z, kpc_per_arcsec=s, pos=pos,
                current=dict(gas_icm=icm, gas_core=core, st_nw=st_nw, st_se=st_se), gas_total=gas_total, subs=subs,
                t_gyr=(0.46, 0.91), centre=tuple(pos['com']), dx=25.0,
                sources=dict(geometry='Kim et al. 2021 (arXiv:2106.00031) Table 2 centroids; cool core 8 arcsec beyond the SE peak (Jee et al. 2014; Ng et al. 2015)',
                             gas='Menanteau et al. 2012 (arXiv:1109.0953): Mgas = 2.2 +- 0.1e14 (R500 about 1.3 Mpc), T = 14.5 keV',
                             stars='Menanteau et al. 2012 Sect. 3.2: 7.5 +- 1.4 (NW) and 5.6 +- 1.3 (SE) x 1e12 inside r200',
                             lensing='Kim et al. 2021 Table 2: M200c 9.9 +2.1/-2.2 (c 2.54) and 6.5 +1.9/-1.4 (c 3.20) x 1e14; total 21.3 +2.5/-2.3',
                             timing='Ng et al. 2015 (arXiv:1412.1826): 0.46 (outgoing) or 0.91 Gyr (returning) after pericentre at 2,400 km/s',
                             speeds='Menanteau et al. 2012: 1290 +- 134 (NW), 1089 +- 200 (SE), 1321 +- 106 km/s (all)'))


# ----------------------------------------------------------------------------- the law on a cluster
def nfw_projected(M200, c, z, R, h=0.7, om=0.3):
    """Projected mass inside R (kpc) of an NFW halo with M200c, concentration c at redshift z."""
    H = 100 * h * np.sqrt(om * (1 + z) ** 3 + 1 - om)                  # km/s/Mpc
    rho_c = 3 * (H / 1e3) ** 2 / (8 * np.pi * G)                       # Msun/kpc^3 (H in km/s/kpc)
    r200 = (3 * M200 / (4 * np.pi * 200 * rho_c)) ** (1 / 3); rs = r200 / c
    m = lambda x: np.log(1 + x) - x / (1 + x)
    rho_s = M200 / (4 * np.pi * rs ** 3 * m(c))
    x = R / rs
    g = np.log(x / 2) + (np.arccosh(1 / x) / np.sqrt(1 - x ** 2) if x < 1 else np.arccos(1 / x) / np.sqrt(x ** 2 - 1))
    return float(4 * np.pi * rho_s * rs ** 3 * g), float(r200)


def solve(spec, law, t_gyr, n=192, dx=None, heat=True, memory=True):
    dx = dx or spec['dx']
    pos, cur = spec['pos'], spec['current']
    ghost_gas, ghost_stars, speeds = {}, [], {}
    r = np.geomspace(1.0, 3000.0, 500); dr = np.gradient(r)
    for sb in spec['subs']:
        g = beta(spec['gas_total'] * sb['gas_share'], sb['rc_pre'], sb['rt_pre'], sb['galaxies'])
        ghost_gas[f"gas_{sb['name'].lower()}_pre"] = g
        stars = [cur[k] for k in sb['stars']]
        rr, ss, gg = V.own_sigma_multi([g], stars, law)
        for st in stars:
            ghost_stars.append((dict(st), (rr, ss)))
        dms = sum(V.shells_on(r, dr, st) for st in stars)
        speeds[sb['name']] = {f'{int(Ra)}kpc': BM.sigma_los_aperture(r, dms, np.interp(r, rr, ss), 0.0, Ra) for Ra in (500.0, 1000.0, 1500.0)}
        speeds[sb['name']]['sigma_r_100kpc'] = float(np.interp(100.0, rr, ss))
    fresh = law['u_kms'] * t_gyr * 1.0227121650537077
    x, y, Se, Sb = V.kappa_map_v4(cur, ghost_gas, ghost_stars, pos, law, n=n, dx=dx, centre=spec['centre'],
                                  fresh_kpc=fresh, heat=heat, memory=memory)
    return dict(x=x, y=y, Se=Se, Sb=Sb, fresh_kpc=fresh, speeds=speeds, dx=dx)


def aperture(sol, xy, R):
    X, Y = np.meshgrid(sol['x'], sol['y'], indexing='ij')
    m = np.hypot(X - xy[0], Y - xy[1]) <= R
    return float(sol['Se'][m].sum() * sol['dx'] ** 2), float(sol['Sb'][m].sum() * sol['dx'] ** 2)


def peaks(sol, pos, smooth=30.0):
    return C4.refined_peaks(sol['x'], sol['y'], sol['Se'] / sol['Se'].max(), pos, smooth_kpc=smooth)


def measure(name, spec, sol):
    pos = spec['pos']; out = dict(fresh_kpc=sol['fresh_kpc'], speeds=sol['speeds'])
    pk = peaks(sol, pos)
    if name == 'macs0025':
        for w in ('se', 'nw'):
            out[f'M300_{w}'], out[f'M300_{w}_baryons'] = aperture(sol, pos[f'{w}_gal'], 300.0)
            p = min(pk, key=lambda d: d[f'dist_{w}_gal'])
            out[f'peak_{w}'] = dict(to_galaxies=p[f'dist_{w}_gal'], to_gas=p['dist_gas_peak'],
                                    galaxies_to_gas=float(np.linalg.norm(pos[f'{w}_gal'] - pos['gas_peak'])),
                                    to_published_lens_peak=p[f'dist_{w}_lens'])
        out['M500_gas_peak'], out['M500_gas_peak_baryons'] = aperture(sol, pos['gas_peak'], 500.0)
        out['sigma_los_1Mpc'] = float(np.mean([sol['speeds'][k]['1000kpc'] for k in sol['speeds']]))
    elif name == 'abell520':
        for k in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P3_jee', 'P3p_jee'):
            out[f'M150_{k}'], out[f'M150_{k}_baryons'] = aperture(sol, pos[k], 150.0)
        out['M710_P3'], _ = aperture(sol, pos['P3'], 710.0)
        near3 = min(pk, key=lambda d: d['dist_P3'])
        out['peak_near_P3'] = dict(dist=near3['dist_P3'], height_rel=near3['kappa'])
        out['peaks'] = [dict(x=p['x'], y=p['y'], h=p['kappa'], nearest=min(('P1', 'P2', 'P3', 'P4', 'P5', 'P6'), key=lambda k: p[f'dist_{k}'])) for p in pk[:8]]
        L = {k: spec['current'][f'st_{k.lower()}']['M'] for k in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6')}
        out['sigma_los_1Mpc'] = float(sum(sol['speeds'][k]['1000kpc'] * L[k] for k in L) / sum(L.values()))
    else:
        for w in ('nw', 'se'):
            for R in (500.0, 1000.0):
                out[f'M{int(R)}_{w}'], _ = aperture(sol, pos[f'{w}_gal'], R)
        for R in (500.0, 1000.0, 1500.0):
            out[f'M{int(R)}_com'], out[f'M{int(R)}_com_baryons'] = aperture(sol, pos['com'], R)
        p = min(pk, key=lambda d: d['dist_se_gal'])
        out['peak_se'] = dict(to_galaxies=p['dist_se_gal'], to_cool_core=p['dist_cool_core'],
                              galaxies_to_cool_core=float(np.linalg.norm(pos['se_gal'] - pos['cool_core'])))
        pn = min(pk, key=lambda d: d['dist_nw_gal']); out['peak_nw_to_galaxies'] = pn['dist_nw_gal']
    return out


def observed_el_gordo(spec):
    """Kim et al. 2021's two NFW haloes, projected about each subcluster and about the centre of mass."""
    z = spec['z']; pos = spec['pos']
    halos = dict(nw=(9.9e14, 2.54, 2.15e14), se=(6.5e14, 3.20, 1.65e14))
    out = {}
    for R in (500.0, 1000.0, 1500.0):
        tot = 0.0
        for w, (M, c, e) in halos.items():
            # projected mass of each halo inside a circle about the centre of mass (numerical, the halo is off-centre)
            d = float(np.linalg.norm(pos[f'{w}_gal'] - pos['com']))
            Rg = np.linspace(0.5, R + d + 1, 3000)
            Mcyl = np.array([nfw_projected(M, c, z, rr)[0] for rr in Rg])
            Sig = np.gradient(Mcyl, Rg) / (2 * np.pi * Rg)
            with np.errstate(invalid='ignore', divide='ignore'):
                cosang = (Rg ** 2 + d ** 2 - R ** 2) / (2 * Rg * d)
            ang = np.where(cosang <= -1, np.pi, np.where(cosang >= 1, 0.0, np.arccos(np.clip(cosang, -1, 1))))
            tot += float(np.trapezoid(2 * ang * Rg * Sig, Rg))
        out[f'M{int(R)}_com'] = tot
    for w, (M, c, e) in halos.items():
        out[f'M500_{w}'] = nfw_projected(M, c, z, 500.0)[0]
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=192)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    law = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    specs = dict(macs0025=macs0025(), abell520=abell520(), el_gordo=el_gordo())
    res = {}
    for name, spec in specs.items():
        res[name] = dict(runs=[])
        for t in spec['t_gyr']:
            sol = solve(spec, law, t, n=args.n)
            m = measure(name, spec, sol); m['t_gyr'] = t
            res[name]['runs'].append(m)
            np.savez_compressed(out / f"{name}_t{t:.2f}.npz", x=sol['x'], y=sol['y'], Se=sol['Se'].astype(np.float32), Sb=sol['Sb'].astype(np.float32))
            print(f"{spec['name']}, t = {t} Gyr (fresh sphere {sol['fresh_kpc']:.0f} kpc): " +
                  json.dumps({k: (round(v / 1e13, 2) if isinstance(v, float) and v > 1e9 else v) for k, v in m.items() if k not in ('speeds', 'peaks')}, default=float), flush=True)
            print('   galaxy speeds (line of sight, 1 Mpc): ' + ', '.join(f"{k} {v['1000kpc']:.0f}" for k, v in sol['speeds'].items()), flush=True)
        if name == 'el_gordo':
            res[name]['observed_nfw'] = observed_el_gordo(spec)
            print('   El Gordo, Kim et al. 2021 NFW haloes projected: ' + json.dumps({k: round(v / 1e14, 2) for k, v in res[name]['observed_nfw'].items()}), flush=True)
    data = {k: dict(name=s['name'], z=s['z'], kpc_per_arcsec=s['kpc_per_arcsec'], positions_kpc={p: v.tolist() for p, v in s['pos'].items()},
                    components=s['current'], subclusters=s['subs'], gas_total=s['gas_total'], t_gyr=s['t_gyr'], sources=s['sources'],
                    **({'gas_fit': s['gas_fit']} if 'gas_fit' in s else {})) for k, s in specs.items()}
    (out / 'collisions_v8.json').write_text(json.dumps(dict(experiment='Three more colliding clusters under our law (round 8)', constants=law,
                                                            inputs=data, results=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'{time.monotonic() - t0:.0f} s')


if __name__ == '__main__':
    main()
