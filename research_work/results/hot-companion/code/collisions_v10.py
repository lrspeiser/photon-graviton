"""Round 10: MACS J0025.4-1222, Abell 520 and El Gordo in the project's own (static) distance law.

Round 8 (collisions_v8.py) turned angles into kiloparsecs, and light into masses, with flat LCDM
distances (H0 = 70, Omega_m = 0.3), as the published measurements it compared with had done. The
project's rules exclude an expanding universe, and its own distance law is static and Euclidean:
    1 + z = exp(alpha D),  alpha = 2.488993e-4 per Mpc,   D_A = D,  D_L = (1 + z) D.
Every published quantity was measured at fixed angles and fluxes, so it converts with a factor
(static / LCDM, at the same angle):
    sizes             D_A
    stellar masses    D_L^2                      (light times a mass-to-light ratio)
    X-ray gas masses  D_L D_A^1.5                (emission measure of a resolved profile)
    lensing masses    D_l D_s / D_ls             (critical density times area; source redshift z_s)
    galaxy speeds     unchanged                  (redshift differences)
The model clusters of round 8 are rescaled accordingly: every length (positions, core radii, cut
radii, the grid) by the size factor, every gas and star mass by its factor. The physics of the
law (u, the fresh sphere u t, the release, the reach) is not rescaled. Apertures that were quoted
in kpc are applied at the same angle (radius times the size factor), and the measured masses
are multiplied by the lensing factor.

    python code/collisions_v10.py --output-dir run-collisions-v10
"""
from __future__ import annotations
import argparse, copy, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
import collisions_v8 as V8
import bullet_v4 as V
import bullet_main_v5 as BM

HERE = Path(__file__).resolve().parent
C_KMS = 299792.458
ALPHA = 2.488993286382367e-4          # per Mpc: the project's static redshift law (research_work/results/joint-light-forward)


PLANCK15 = (67.74, 0.3089)             # Kim et al. 2021 (El Gordo lensing): "Planck 2016", D_A = 1636 Mpc


def lcdm(z, H0=70.0, Om=0.3):
    chi = quad(lambda x: 1 / np.sqrt(Om * (1 + x) ** 3 + 1 - Om), 0, z)[0] * C_KMS / H0
    return dict(DA=chi / (1 + z), DL=chi * (1 + z), chi=chi)


# Which static geometry turns angles into sizes (round 11). Both share the brightness law D_L = (1 + z) D
# (energy loss and arrival-rate stretching), with D = ln(1 + z) / ALPHA the path length; they differ in D_A:
#   'fixed'  (the fixed-material transport branch, adopted since round 10): D_A = D;
#   'metric' (the material-coasting geometry of five_candidate_tests/prior/conformal_action_derivation_derivation.md,
#             eq. 13, where rulers co-scale): D_A = D / (1 + z), so that D_L = (1 + z)^2 D_A.
VARIANT = 'fixed'


def static(z):
    D = np.log1p(z) / ALPHA
    DA = D if VARIANT == 'fixed' else D / (1 + z)
    return dict(DA=DA, DL=(1 + z) * D, chi=D)


def static_Dls(zl, zs):
    Sl, Ss = static(zl), static(zs)
    return (Ss['DA'] - Sl['DA']) if VARIANT == 'fixed' else (Ss['chi'] - Sl['chi']) / (1 + zs)


def factors(z, zs, ref=(70.0, 0.3)):
    """static / LCDM(ref) at the same angle and flux; zs is the effective source redshift in that LCDM."""
    L, S = lcdm(z, *ref), static(z)
    Ls, Ss = lcdm(zs, *ref), static(zs)
    lens_L = L['DA'] * Ls['DA'] / ((Ls['chi'] - L['chi']) / (1 + zs))      # flat LCDM: D_ls = (chi_s - chi_l)/(1 + z_s)
    lens_S = S['DA'] * Ss['DA'] / static_Dls(z, zs)                         # Euclidean (or the metric variant)
    return dict(size=S['DA'] / L['DA'], stars=(S['DL'] / L['DL']) ** 2,
                gas=(S['DL'] / L['DL']) * (S['DA'] / L['DA']) ** 1.5, lens=lens_S / lens_L,
                kpc_per_arcsec_lcdm=L['DA'] * 1e3 / 206264.806, kpc_per_arcsec_static=S['DA'] * 1e3 / 206264.806, z=z, z_source=zs,
                reference_cosmology=dict(H0=ref[0], Om=ref[1]))


def z_for_beta(zl, beta, ref=(70.0, 0.3)):
    """The source redshift whose D_ls/D_s equals beta in LCDM(ref) (to express a published <beta> as z_eff)."""
    from scipy.optimize import brentq
    cl = lcdm(zl, *ref)['chi']
    return brentq(lambda zs: (lcdm(zs, *ref)['chi'] - cl) / lcdm(zs, *ref)['chi'] - beta, zl + 1e-3, 20.0)


# How each measurement converts (the papers' own conventions; literature check of round 10):
#  MACS J0025: Bradac et al. 2008, flat 0.3/0.7/70; strong + weak lensing, weak-lensing sources at an assumed z = 1.4
#              (one spectroscopic arc system at z = 2.38); stars = F814W light x M/L_K 0.74 (no age fitting).
#  Abell 520:  Jee et al. 2014 and Clowe et al. 2012 clumps, <D_ls/D_s> = 0.73 (z_eff about 0.85, UDF photo-z);
#              Mahdavi et al. 2007, 710 kpc about P3, <D_ls/D_s> = 0.59 (CFHTLS photo-z); light x M/L 2 (ours).
#  El Gordo:   Kim et al. 2021, Planck 2015 cosmology, <D_ls/D_s> = 0.254 (z_eff 1.31); aperture masses about the
#              centre of mass read from their Fig. 11 (5.8, 14.9, 20.1 x 1e14 inside 0.5, 1, 1.5 Mpc, +-9-14%);
#              stars from SED fits with ages up to 7.0 Gyr (Menanteau et al. 2012, 0.27/0.73/70).
CONVENTIONS = dict(
    macs0025=dict(ref=(70.0, 0.3), zs=1.4),
    abell520=dict(ref=(70.0, 0.3), zs=None, beta_clumps=0.73, beta_710=0.59),
    el_gordo=dict(ref=PLANCK15, zs=1.31, aperture_Mpc=(0.5, 1.0, 1.5), aperture_M=(5.8e14, 14.9e14, 20.1e14), aperture_err=0.12))


# The stars' IMF basis (round 11). The law's u is calibrated on X-COP, whose stellar masses assume a Chabrier IMF
# (Ghizzardi et al. 2021), so every star mass the law is compared with should be on that basis. MACS J0025's
# masses use M/L_K = 0.74 after Drory et al. (2004), whose masses assume a Salpeter IMF (the IMF is not stated by
# Bradac et al. 2008; inferred in the round-11 literature audit): x 10^-0.25 to Chabrier, the convention of the
# SLACS code (Auger et al. 2009). El Gordo (Menanteau et al. 2012) is Chabrier; Abell 520 and the Bullet use
# M/L = 2 after Clowe et al. 2006, who cite Kauffmann et al. 2003 (Kroupa, close to Chabrier).
CHABRIER_BASIS = dict(macs0025=10 ** -0.25, abell520=1.0, el_gordo=1.0)


def chabrier_basis(name):
    return CHABRIER_BASIS.get(name, 1.0)


def lens_factor(name, which='main'):
    z = dict(macs0025=0.586, abell520=0.201, el_gordo=0.870)[name]
    c = CONVENTIONS[name]
    if name == 'abell520':
        zs = z_for_beta(z, c['beta_clumps'] if which == 'main' else c['beta_710'], c['ref'])
    else:
        zs = c['zs']
    return factors(z, zs, c['ref'])


def rescale(spec, f, star_extra=1.0):
    """Lengths times f['size']; gas masses times f['gas']; stars times f['stars'] * star_extra."""
    s = copy.deepcopy(spec); fs = f['size']
    s['pos'] = {k: np.asarray(v) * fs for k, v in spec['pos'].items()}
    for k, c in s['current'].items():
        c['scale'] *= fs; c['rt'] *= fs
        c['M'] *= f['gas'] if k.startswith('gas') else f['stars'] * star_extra
    s['gas_total'] = spec['gas_total'] * f['gas']
    for sb in s['subs']:
        sb['rc_pre'] *= fs; sb['rt_pre'] *= fs
    s['centre'] = tuple(np.asarray(spec['centre']) * fs)
    s['dx'] = spec['dx'] * fs
    s['kpc_per_arcsec'] = f['kpc_per_arcsec_static']
    return s


def solve(spec, law, t_gyr, fs, n=192):
    """collisions_v8.solve with the speed apertures at the same angles (500, 1000, 1500 kpc in LCDM)."""
    dx = spec['dx']; pos, cur = spec['pos'], spec['current']
    ghost_gas, ghost_stars, speeds = {}, [], {}
    r = np.geomspace(1.0, 3000.0 * fs, 500); dr = np.gradient(r)
    for sb in spec['subs']:
        g = V8.beta(spec['gas_total'] * sb['gas_share'], sb['rc_pre'], sb['rt_pre'], sb['galaxies'])
        ghost_gas[f"gas_{sb['name'].lower()}_pre"] = g
        stars = [cur[k] for k in sb['stars']]
        rr, ss, gg = V.own_sigma_multi([g], stars, law, rmax=3000.0 * max(fs, 1.0))
        for st in stars:
            ghost_stars.append((dict(st), (rr, ss)))
        dms = sum(V.shells_on(r, dr, st) for st in stars)
        speeds[sb['name']] = {f'{int(Ra)}kpc': BM.sigma_los_aperture(r, dms, np.interp(r, rr, ss), 0.0, Ra * fs) for Ra in (500.0, 1000.0, 1500.0)}
    fresh = law['u_kms'] * t_gyr * 1.0227121650537077
    x, y, Se, Sb = V.kappa_map_v4(cur, ghost_gas, ghost_stars, pos, law, n=n, dx=dx, centre=spec['centre'],
                                  fresh_kpc=fresh, heat=True, memory=True)
    return dict(x=x, y=y, Se=Se, Sb=Sb, fresh_kpc=fresh, speeds=speeds, dx=dx)


def measure(name, spec, sol, fs):
    pos = spec['pos']; out = dict(fresh_kpc=sol['fresh_kpc'], speeds=sol['speeds'])
    ap = lambda xy, R: V8.aperture(sol, xy, R * fs)
    pk = V8.peaks(sol, pos)
    if name == 'macs0025':
        for w in ('se', 'nw'):
            out[f'M300_{w}'], out[f'M300_{w}_baryons'] = ap(pos[f'{w}_gal'], 300.0)
            p = min(pk, key=lambda d: d[f'dist_{w}_gal'])
            out[f'peak_{w}'] = dict(to_galaxies=p[f'dist_{w}_gal'], to_gas=p['dist_gas_peak'],
                                    galaxies_to_gas=float(np.linalg.norm(pos[f'{w}_gal'] - pos['gas_peak'])))
        out['M500_gas_peak'], _ = ap(pos['gas_peak'], 500.0)
        out['sigma_los_1p5Mpc'] = float(np.mean([sol['speeds'][k]['1500kpc'] for k in sol['speeds']]))
    elif name == 'abell520':
        for k in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
            out[f'M150_{k}'], out[f'M150_{k}_baryons'] = ap(pos[k], 150.0)
        out['M710_P3'], _ = ap(pos['P3'], 710.0)
        out['sigma_los_500kpc'] = {k: sol['speeds'][k]['500kpc'] for k in ('P1', 'P2', 'P4', 'P5')}
    else:
        for R in (500.0, 1000.0, 1500.0):
            out[f'M{int(R)}_com'], out[f'M{int(R)}_com_baryons'] = ap(pos['com'], R)
        # apertures at Kim et al.'s Planck-2015 radii, in the static law
        fk = lens_factor('el_gordo')['size']
        for R in CONVENTIONS['el_gordo']['aperture_Mpc']:
            out[f'Mkim{int(R * 1000)}_com'], _ = V8.aperture(sol, pos['com'], R * 1000 * fk)
        p = min(pk, key=lambda d: d['dist_se_gal'])
        out['peak_se'] = dict(to_galaxies=p['dist_se_gal'], to_cool_core=p['dist_cool_core'],
                              galaxies_to_cool_core=float(np.linalg.norm(pos['se_gal'] - pos['cool_core'])))
    return out


# published measurements (LCDM conventions of the papers), as in round 8
OBS = dict(
    macs0025=dict(M300_se=(2.5e14, 1.0e14, 1.7e14), M300_nw=(2.6e14, 0.5e14, 1.4e14), M500_gas_peak=(6.2e14, 1.2e14, 4.0e14),
                  sigma_los_1p5Mpc=(835.0, 59.0, 59.0)),
    abell520=dict(M150_P1=((2.10e13, 0.43e13), (2.81e13, 0.67e13)), M150_P2=((4.05e13, 0.28e13), (4.16e13, 0.67e13)),
                  M150_P3=((3.35e13, 0.34e13), (2.84e13, 0.64e13)), M150_P4=((4.23e13, 0.28e13), (5.59e13, 0.68e13)),
                  M150_P5=((2.93e13, 0.39e13), (3.17e13, 0.66e13)), M150_P6=(None, (3.68e13, 0.68e13)),
                  M710_P3=(5.0e14, 0.55e14)),
    el_gordo=dict(sigma_NW=(1290.0, 134.0), sigma_SE=(1089.0, 200.0)))


def compare(name, m, f, obs_gordo=None):
    """Model against the measurements converted to the static law (lensing masses times each measurement's
    lensing factor, from the paper's own cosmology and source redshifts)."""
    fl = lens_factor(name)['lens']; rows = []
    if name == 'macs0025':
        for k in ('M300_se', 'M300_nw', 'M500_gas_peak'):
            v, up, lo = OBS[name][k]; v, up, lo = v * fl, up * fl, lo * fl
            z = (m[k] - v) / (up if m[k] > v else lo)
            rows.append(dict(check=k, model=m[k], measured=v, z=z))
        v, e, _ = OBS[name]['sigma_los_1p5Mpc']; rows.append(dict(check='sigma_los_1p5Mpc', model=m['sigma_los_1p5Mpc'], measured=v, z=(m['sigma_los_1p5Mpc'] - v) / e))
    elif name == 'abell520':
        for k in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
            jee, clowe = OBS[name][f'M150_{k}']
            vals = [x for x in (jee, clowe) if x is not None]
            lo = min(x[0] for x in vals) * fl; hi = max(x[0] for x in vals) * fl; err = max(x[1] for x in vals) * fl
            mv = m[f'M150_{k}']
            z = 0.0 if lo <= mv <= hi else ((mv - hi) / err if mv > hi else (mv - lo) / err)
            rows.append(dict(check=f'M150_{k}', model=mv, measured_range=[lo, hi], z=z))
        f7 = lens_factor(name, '710')['lens']
        v, e = OBS[name]['M710_P3']; rows.append(dict(check='M710_P3', model=m['M710_P3'], measured=v * f7, z=(m['M710_P3'] - v * f7) / (e * f7)))
    else:
        c = CONVENTIONS['el_gordo']
        for R, M in zip(c['aperture_Mpc'], c['aperture_M']):
            k = f'Mkim{int(R * 1000)}_com'; v = M * fl
            rows.append(dict(check=f'{k} (aperture)', model=m[k], measured=v, ratio=m[k] / v, z=(m[k] - v) / (c['aperture_err'] * v)))
        f70 = factors(0.870, c['zs'], (70.0, 0.3))['lens']        # the NFW sums were evaluated in 0.3/0.7/70
        for R in (500, 1000, 1500):
            v = obs_gordo[f'M{R}_com'] * f70
            rows.append(dict(check=f'M{R}_com (NFW sum)', model=m[f'M{R}_com'], measured=v, ratio=m[f'M{R}_com'] / v, z=(m[f'M{R}_com'] - v) / (0.15 * v)))
        for w in ('NW', 'SE'):
            v, e = OBS[name][f'sigma_{w}']; mv = m['speeds'][w]['1000kpc']
            rows.append(dict(check=f'sigma_{w}', model=mv, measured=v, z=(mv - v) / e))
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=192)
    args = ap.parse_args(); out = args.output_dir
    out.mkdir(parents=True, exist_ok=True); t0 = time.monotonic()
    law = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    base = dict(macs0025=V8.macs0025(), abell520=V8.abell520(), el_gordo=V8.el_gordo())
    t_main = dict(macs0025=0.5, abell520=0.3, el_gordo=0.46)
    res = dict(experiment='round 10: the three collisions in the static distance law', distance_law=dict(alpha_per_Mpc=ALPHA),
               law_constants=law, results={})
    for name, spec in base.items():
        f = factors(spec['z'], 1.0)                        # inputs: positions from coordinates, masses as used in round 8 (0.3/0.7/70)
        f.update(lens_conventions={w: lens_factor(name, w) for w in (('main', '710') if name == 'abell520' else ('main',))})
        obs_gordo = V8.observed_el_gordo(spec) if name == 'el_gordo' else None
        entry = dict(factors=f, runs={})
        variants = [('published stars', 1.0)] + ([('stars x 2', 2.0)] if name == 'el_gordo' else [])
        for label, sx in variants:
            st = rescale(spec, f, star_extra=sx)
            sol = solve(st, law, t_main[name], f['size'], n=args.n)
            m = measure(name, st, sol, f['size'])
            rows = compare(name, m, f, obs_gordo)
            entry['runs'][label] = dict(measures={k: v for k, v in m.items() if k != 'speeds'}, speeds=m['speeds'], comparison=rows)
            print(f"{spec['name']} (static; size x{f['size']:.2f}, stars x{f['stars']:.2f}, gas x{f['gas']:.2f}, lensing x{lens_factor(name)['lens']:.2f}), {label}:", flush=True)
            for r in rows:
                meas = r.get('measured', r.get('measured_range'))
                ms = f"{meas:.3g}" if isinstance(meas, float) else f"{meas[0]:.3g}-{meas[1]:.3g}"
                print(f"   {r['check']:18s} model {r['model']:.3g}  measured {ms}  z {r['z']:+.2f}", flush=True)
        res['results'][name] = entry
    res['seconds'] = time.monotonic() - t0
    (out / 'collisions_v10.json').write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {out / 'collisions_v10.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
