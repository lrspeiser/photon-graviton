#!/usr/bin/env python3
"""Round 5: the Bullet main cluster's galaxy speeds, with its outskirts taken from independent data.

    python bullet_main_v5.py --output-dir ../run-bullet-main-v5

Round 4 found our law gives the main cluster's galaxies about 1,000 km/s (line of sight,
averaged over the region Barrena et al. 2002 surveyed), against 1,249 +109/-100 km/s from 71
galaxies. That main-cluster model was built only from Clowe et al.'s two 100 kpc apertures:
its stars were an NFW profile fitted there and cut off at 1.5 Mpc.

What the other data say (checked below):
* Gas. The fitted beta model agrees with the ROSAT beta model of the whole system (Ota &
  Mitsuda 2004, as used by Paraficz et al. 2016) at 1-2 Mpc, and with the gas share of the
  lensing mass inside 250 kpc (Paraficz et al.: 9 +- 3% of 2.5e14). Gas is not the problem.
* Stars inside 340 kpc. Barrena et al.'s R-band light, 1.0e12 Lsun in 5.4 arcmin^2 (0.37 Mpc^2)
  excluding the subcluster, is 2-3e12 Msun of stars for M/L_R 2-3. The model has 2.1e12
  (M/L_I 1.5) to 2.7e12 (M/L_I 2). Consistent.
* Stars beyond 340 kpc: no counts exist. In the 7 X-COP clusters with measured stellar profiles
  (Ghizzardi et al. 2021) -- the clusters that fixed u = 197 km/s -- the cumulative star/gas
  ratio is 0.074 (median) at 0.5 R500 and 0.062 at R500 (range 0.035-0.073). The round-4
  model has 0.045 and 0.035 at M/L 2, 0.034 and 0.026 at M/L 1.5: about half the usual stars
  in the outskirts, an artefact of extrapolating two central apertures.

The test keeps the central stars (Clowe et al.) and adds an extended population of galaxies in
the outskirts (NFW shape, scale 800 kpc, cut at 3 Mpc, typical of galaxy counts in clusters):
1. at the X-COP median star/gas ratio (a template, nothing fitted);
2. at the amount the measured galaxy speed requires under our law (one number inferred from
   one measurement), then checked against two independent constraints: Barrena et al.'s
   light inside 340 kpc, and the X-COP range of star/gas ratios;
and recomputes the lensing with the companion's memory (bullet_v4). Also tested: the 1.5 Mpc
cut-off on its own, the outer scale (500 / 1200 kpc), central M/L 1.0-2.0, and orbits that are
partly radial (constant anisotropy 0.3), common in cluster outskirts.

Lensing is compared two ways:
* Clowe et al. 2006's mean kappa in 100 kpc apertures. Their paper states that the method
  "systematically underestimate[s] kappa in the cores of massive clusters ... our measurements
  of kappa in the peaks of the components are only lower bounds". So 0.36 and 0.20 are floors.
* The projected mass inside 250 kpc of each brightest galaxy from strong (+ weak) lensing:
  Bradac et al. 2006, main 2.8 +- 0.2 and sub 2.3 +- 0.2 x 1e14 Msun; Paraficz et al. 2016
  (strong lensing, 14 multiple-image systems), main 2.5 +- 0.1 and sub 2.0 +- 0.2.
The model box has to hold the line of sight: --dx 25 gives a 4.8 Mpc box (default 2.9 Mpc).
The subcluster's pre-collision size is also varied (1:8, 1:4, 1:3, 1:2 of the main's baryons).
"""
from __future__ import annotations
import argparse, copy, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import bullet_v3 as B
import bullet_v4 as V
import collisions_v4 as C4
import law as L
import run as R

G = B.G
R500 = 1500.0                  # kpc; X-ray and SZ estimates for the main cluster are 1.4-1.6 Mpc
BARRENA_RAP = np.sqrt(1.8e6 / np.pi)   # members lie in about 1.8 Mpc^2: equivalent radius 757 kpc
OBS_SIGMA = (1249.0, 109.0, 100.0)
OBS_M250 = dict(bradac2006=dict(main=(2.8e14, 0.2e14), sub=(2.3e14, 0.2e14)),
                paraficz2016=dict(main=(2.5e14, 0.1e14), sub=(2.0e14, 0.2e14)))


def mass_within(x, y, kap, centre, R, scrit):
    X, Y = np.meshgrid(x, y, indexing='ij'); dx = x[1] - x[0]
    m = np.hypot(X - centre[0], Y - centre[1]) <= R
    return float(scrit * kap[m].sum() * dx * dx)


def xcop_ratio():
    """Cumulative M*(<r) / Mgas(<r) of the 7 X-COP clusters with measured stellar profiles."""
    D = json.loads((R.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    XS = np.geomspace(0.05, 1.5, 30); fr = []
    for c in D.values():
        if 'stellar_mass' in c:
            R5 = c['header']['R500_kpc']
            fr.append(R.logi(XS * R5, c['stellar_mass']['radius_kpc'], c['stellar_mass']['Mstar']) /
                      R.logi(XS, c['gas_mass']['RADIUS'], c['gas_mass']['MGAS']))
    fr = np.array(fr)
    return XS, np.median(fr, 0), fr.min(0), fr.max(0)


def cum_mass(c, r):
    prof = B.rho_beta if c['kind'] == 'beta' else B.rho_nfw
    rr = np.geomspace(1e-2, max(c['rt'], r.max()), 6000)
    m = np.concatenate([[0], np.cumsum(np.diff(rr) * 0.5 * ((4 * np.pi * rr ** 2 * prof(rr, c['scale'], c['rt']))[1:] +
                                                           (4 * np.pi * rr ** 2 * prof(rr, c['scale'], c['rt']))[:-1]))])
    return c['M'] / B.mass3d(prof, (c['scale'],), c['rt']) * np.interp(r, rr, m)


def projected_mass(c, Rap):
    r = np.geomspace(0.5, c['rt'], 20000); dr = np.gradient(r)
    dm = V.shells_on(r, dr, c)
    f = np.where(r <= Rap, 1.0, 1 - np.sqrt(np.clip(1 - (Rap / r) ** 2, 0, 1)))
    return float(np.sum(dm * f))


def outer_component(st_inner, gas, target_x, target_ratio, r500, rt=3000.0):
    """Extended stellar component (NFW shape, cut at rt) that brings the cumulative star/gas
    ratio up to the target between 0.3 and 1 R500."""
    xs = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0]); rr = xs * r500
    want = np.interp(np.log(xs), np.log(target_x), target_ratio) * cum_mass(gas, rr)
    have = cum_mass(st_inner, rr)
    def resid(q):
        c = dict(kind='nfw', M=float(np.exp(q[0])), scale=float(np.exp(q[1])), rt=rt, centre=st_inner['centre'])
        return np.log((have + cum_mass(c, rr)) / want)
    q = least_squares(resid, np.log([5e12, 800.0]), bounds=(np.log([1e10, 100.0]), np.log([1e15, 1e4]))).x
    return dict(kind='nfw', M=float(np.exp(q[0])), scale=float(np.exp(q[1])), rt=rt, centre=st_inner['centre'])


def fixed_outer(st_inner, M, rs=800.0, rt=3000.0):
    return dict(kind='nfw', M=float(M), scale=float(rs), rt=rt, centre=st_inner['centre'])


def jeans_aniso(r, dr, rho, g, beta):
    """rho sigma_r^2 = r^(-2 beta) * integral_r^inf r'^(2 beta) rho g dr'  (constant anisotropy)."""
    w = r ** (2 * beta)
    tail = np.cumsum((w * rho * g * dr)[::-1])[::-1]
    return np.where(rho > 1e-30, tail / np.maximum(w * rho, 1e-300), 0.0)


def own_sigma(gas_list, star_list, consts, beta=0.0, n=600, rmax=4000.0, iters=600):
    """Stars' radial dispersion in the gravity our law makes from the cluster's own gas and
    stars, iterated to a fixed point. With anisotropy the heat weight is <v^2>/u^2 =
    (3 - 2 beta) sigma_r^2 / u^2."""
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    r = np.geomspace(1.0, rmax, n); dr = np.gradient(r)
    dmg = sum(V.shells_on(r, dr, c) for c in gas_list); dms = sum(V.shells_on(r, dr, c) for c in star_list)
    gN = G * np.cumsum(dmg + dms) / r ** 2
    W = L.shell_weights(r, r); rho_s = dms / (4 * np.pi * r ** 2 * dr)
    sig2 = jeans_aniso(r, dr, rho_s, gN + np.sqrt(a * gN), beta)
    for _ in range(iters):
        S = G * (W @ (L.k_from_sig2(sig2, u, pref=3 - 2 * beta) * dms)) / r ** 2
        g = gN + np.exp(-gN / (lam * a)) * np.sqrt(a * (gN + S))
        new = jeans_aniso(r, dr, rho_s, g, beta)
        if np.max(np.abs(new - sig2) * dms) < 1e-8 * np.max(new * dms): sig2 = new; break
        sig2 = 0.5 * sig2 + 0.5 * new
    return r, np.sqrt(np.maximum(sig2, 0.0)), g, dms


def sigma_los_aperture(r, dms, sig_r, beta, Rap):
    """Line-of-sight dispersion averaged over a projected circle of radius Rap. A shell of radius
    r contributes the caps inside the circle; on them <sin^2 theta> sets the anisotropy term."""
    x = np.clip(Rap / r, 0.0, 1.0)
    c0 = np.where(r <= Rap, 0.0, np.sqrt(1 - x ** 2))
    f = 1 - c0
    sin2 = 1 - (1 + c0 + c0 ** 2) / 3
    w = dms * f
    return float(np.sqrt(np.sum(w * sig_r ** 2 * (1 - beta * sin2)) / np.sum(w)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=192); ap.add_argument('--dx', type=float, default=15.0)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    pos = B.positions(); comps0, _ = B.fit_components(pos); scrit = B.sigma_crit()
    XS, med, lo, hi = xcop_ratio()
    gas = dict(comps0['gas_main'], centre='main_bcg')      # before the collision, centred on its galaxies

    # ---- data checks on the round-4 model
    r_chk = np.array([250.0, 340.0, 0.5 * R500, R500])
    checks = dict(
        gas_projected_250kpc=projected_mass(comps0['gas_main'], 250.0),
        gas_paraficz_250kpc=[0.09 * 2.5e14, 0.03 * 2.5e14],
        gas_3d_1_1p5_2Mpc=[float(cum_mass(gas, np.array([x]))[0]) for x in (1000.0, 1500.0, 2000.0)],
        stars_projected_340kpc_ml2=projected_mass(comps0['st_main'], 340.0),
        stars_barrena_340kpc=[2.0e12, 3.0e12],
        xcop_ratio_at=dict(x=[0.17, 0.23, 0.5, 1.0], median=[float(np.interp(np.log(x), np.log(XS), med)) for x in (0.17, 0.23, 0.5, 1.0)],
                           low=[float(np.interp(np.log(x), np.log(XS), lo)) for x in (0.17, 0.23, 0.5, 1.0)],
                           high=[float(np.interp(np.log(x), np.log(XS), hi)) for x in (0.17, 0.23, 0.5, 1.0)]),
        round4_ratio_ml2_at=[float(cum_mass(comps0['st_main'], np.array([x]))[0] / cum_mass(gas, np.array([x]))[0]) for x in r_chk])
    print('checks:', json.dumps(checks, default=float), flush=True)

    rows = []

    def build(ml, outer, rs_out, rt_inner):
        cm = copy.deepcopy(comps0)
        for k in ('st_main', 'st_sub'): cm[k]['M'] = comps0[k]['M'] * ml / 2.0
        inner = dict(cm['st_main'])
        if rt_inner != inner['rt']:           # same density law, cut off further out
            inner['M'] *= B.mass3d(B.rho_nfw, (inner['scale'],), rt_inner) / B.mass3d(B.rho_nfw, (inner['scale'],), inner['rt'])
            inner['rt'] = rt_inner
        stars = [inner]
        if outer == 'median':
            stars.append(outer_component(inner, gas, XS, med, R500))
        elif isinstance(outer, (int, float)) and outer > 0:
            stars.append(fixed_outer(inner, outer, rs_out))
        return cm, stars

    def speeds(stars, beta):
        r, sr, g, dms = own_sigma([gas], stars, consts, beta=beta)
        return r, sr, dms, {f'{int(Ra)}kpc': sigma_los_aperture(r, dms, sr, beta, Ra) for Ra in (500.0, BARRENA_RAP, 1000.0)}

    def solve_outer(ml, target, rs_out, beta):
        """Outer stellar mass for which our law gives the measured line-of-sight speed."""
        key = f'{int(BARRENA_RAP)}kpc'
        f = lambda lm: speeds(build(ml, 10 ** lm, rs_out, 1500.0)[1], beta)[3][key] - target
        lo_, hi_ = 11.0, 13.8
        if f(lo_) > 0: return 0.0
        for _ in range(40):
            mid = 0.5 * (lo_ + hi_)
            lo_, hi_ = (mid, hi_) if f(mid) < 0 else (lo_, mid)
        return 10 ** (0.5 * (lo_ + hi_))

    def variant(label, ml=1.5, outer='none', rs_out=800.0, beta=0.0, rt_inner=1500.0, lens=False, sub_ratio=8.0):
        cm, stars = build(ml, outer, rs_out, rt_inner)
        r, sr, dms, los = speeds(stars, beta)
        rr = np.array([0.3, 0.5, 1.0]) * R500
        ratio = [float(sum(cum_mass(st, np.array([x]))[0] for st in stars) / cum_mass(gas, np.array([x]))[0]) for x in rr]
        row = dict(case=label, mass_to_light_centre=ml, outskirts=outer, outer_scale_kpc=rs_out, anisotropy=beta, rt_inner_kpc=rt_inner,
                   outer_component=(stars[1] if len(stars) > 1 else None),
                   stars_projected_340kpc=float(sum(projected_mass(st, 340.0) for st in stars)),
                   stars_3d_r500=float(sum(cum_mass(st, np.array([R500]))[0] for st in stars)),
                   star_gas_ratio_03_05_1_R500=ratio,
                   sigma_r_kms={f'{x}kpc': float(np.interp(x, r, sr)) for x in (100, 400, 800, 1500)},
                   sigma_los_kms=los)
        bar = los[f'{int(BARRENA_RAP)}kpc']
        row['gap_sigma'] = float((OBS_SIGMA[0] - bar) / (OBS_SIGMA[2] if bar < OBS_SIGMA[0] else OBS_SIGMA[1]))
        msg = (f"{label}: stars in 340 kpc {row['stars_projected_340kpc']:.2e} (Barrena 2-3e12); star/gas at 0.3/0.5/1 R500 " +
               '/'.join(f'{v:.3f}' for v in ratio) + f"; line-of-sight speed in Barrena's region {bar:.0f} km/s"
               f" (0.5/1 Mpc: {los['500kpc']:.0f}/{los['1000kpc']:.0f}); gap {row['gap_sigma']:+.1f} sigma")
        if lens:
            ghost_gas, sub_sats, pre = V.pre_collision_models(cm, sub_ratio, 0.07, 150.0)
            rs_s, ss_s, _ = V.own_sigma_multi([dict(cm['gas_sub']), ghost_gas['gas_atm_ghost']], [cm['st_sub'], sub_sats], consts)
            s_eff = sr * np.sqrt((3 - 2 * beta) / 3)      # heat weight (3 - 2 beta) sigma_r^2 / u^2
            cur = dict(cm); cur['st_main'] = stars[0]
            ghost_stars = [(dict(stars[0]), (r, s_eff))]
            if len(stars) > 1:
                cur['st_main_outer'] = stars[1]; ghost_stars.append((dict(stars[1]), (r, s_eff)))
            ghost_stars += [(dict(cm['st_sub']), (rs_s, ss_s))] + ([(sub_sats, (rs_s, ss_s))] if sub_sats['M'] > 0 else [])
            x, y, Se, Sb = V.kappa_map_v4(cur, ghost_gas, ghost_stars, pos, consts, n=args.n, dx=args.dx, fresh_kpc=30.0)
            kap = Se / scrit
            dec = B.clowe_decomposition(x, y, kap, pos)
            pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0)
            mp = min(pk, key=lambda d: d['dist_main_bcg']); sp = min(pk, key=lambda d: d['dist_sub_bcg'])
            m250 = {w: mass_within(x, y, kap, pos[f'{w}_bcg'], 250.0, scrit) for w in ('main', 'sub')}
            row['lensing'] = dict(main=dec['main_bcg'], sub=dec['sub_bcg'], main_gas=dec['main_plasma'], sub_gas=dec['sub_plasma'],
                                  main_peak_kpc=mp['dist_main_bcg'], sub_peak_kpc=sp['dist_sub_bcg'], kappa_min=float(kap.min()),
                                  mass_250kpc=m250, sub_ratio=sub_ratio, sub_sigma_before_kms={f'{x0}kpc': float(np.interp(x0, rs_s, ss_s)) for x0 in (50, 100, 200)})
            np.savez_compressed(out / f"kappa_{len(rows) + 1}.npz", x=x, y=y, kappa=kap.astype(np.float32))
            msg += (f"\n    lensing (sub 1:{sub_ratio:.0f}): kappa main {dec['main_bcg']:.2f} (floor 0.36)  sub {dec['sub_bcg']:.2f} (floor 0.20)  gas "
                    f"{dec['main_plasma']:.2f} / {dec['sub_plasma']:.2f} (0.05 / 0.02); peaks {mp['dist_main_bcg']:.0f} / {sp['dist_sub_bcg']:.0f} kpc from the galaxies;"
                    f" mass inside 250 kpc: main {m250['main'] / 1e14:.2f} (2.5-2.8), sub {m250['sub'] / 1e14:.2f} (2.0-2.3) x 1e14;"
                    f" sub stars before {np.interp(100, rs_s, ss_s):.0f} km/s at 100 kpc")
        rows.append(row); print(msg, flush=True)
        return row

    variant('A. round 4: central stars only, cut off at 1.5 Mpc (M/L 1.5)', lens=True)
    variant('B. the same stars, not cut off (to 3 Mpc)', rt_inner=3000.0)
    variant('C. outskirts at the X-COP median star/gas ratio', outer='median')
    solved = {}
    for tgt, tag in ((OBS_SIGMA[0], 'measured'), (OBS_SIGMA[0] - OBS_SIGMA[2], 'measured - 1 sigma'), (OBS_SIGMA[0] + OBS_SIGMA[1], 'measured + 1 sigma')):
        solved[tag] = solve_outer(1.5, tgt, 800.0, 0.0)
        print(f'outer stars needed for {tgt:.0f} km/s: {solved[tag]:.3e} Msun', flush=True)
    variant('D. outskirts inferred from the galaxy speed (1,249 km/s)', outer=solved['measured'], lens=True)
    variant('E. the same, speed 1 sigma lower (1,149 km/s)', outer=solved['measured - 1 sigma'])
    variant('F. the same, speed 1 sigma higher (1,358 km/s)', outer=solved['measured + 1 sigma'])
    variant('G. D with orbits partly radial (anisotropy 0.3)', outer=solved['measured'], beta=0.3)
    for rs_o in (500.0, 1200.0):
        m_rs = solve_outer(1.5, OBS_SIGMA[0], rs_o, 0.0)
        variant(f'H. outer scale {rs_o:.0f} kpc, inferred from the speed', outer=m_rs, rs_out=rs_o)
    m_ml1 = solve_outer(1.0, OBS_SIGMA[0], 800.0, 0.0)
    variant('I. central stars at M/L 1.0, outskirts inferred from the speed', ml=1.0, outer=m_ml1, lens=True)
    m_ml2 = solve_outer(2.0, OBS_SIGMA[0], 800.0, 0.0)
    variant('J. central stars at M/L 2.0 (published), outskirts inferred from the speed', ml=2.0, outer=m_ml2)
    for rat in (4.0, 3.0, 2.0):
        variant(f'K. D with the subcluster 1:{rat:.0f} of the main before the collision', outer=solved['measured'], lens=True, sub_ratio=rat)
    # Clowe et al.'s kappas are floors, so the published star masses (M/L 2) need no lowering
    for rat in (8.0, 4.0, 3.0):
        variant(f'L. published star masses (M/L 2), outskirts from the speed, subcluster 1:{rat:.0f}', ml=2.0, outer=m_ml2, lens=True, sub_ratio=rat)

    (out / 'bullet_main_v5.json').write_text(json.dumps(dict(
        experiment='Bullet main cluster: galaxy speeds with outskirts from independent data (round 5)',
        observed=dict(sigma_los='1249 +109/-100 km/s, 71 main-cluster galaxies in about 1.8 Mpc^2 (Barrena et al. 2002)',
                      kappa='Clowe et al. 2006: 0.36 +- 0.06, 0.20 +- 0.05 (stated as lower bounds), gas residuals 0.05 +- 0.06, 0.02 +- 0.06',
                      mass_250kpc=OBS_M250),
        r500_kpc=R500, barrena_aperture_kpc=BARRENA_RAP, checks=checks, outer_mass_needed=solved, runs=rows, constants=consts,
        grid=dict(n=args.n, dx_kpc=args.dx), seconds=time.monotonic() - t0), indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
