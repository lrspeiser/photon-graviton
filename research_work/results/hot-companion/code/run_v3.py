#!/usr/bin/env python3
"""Round 3: hot-companion gravity with the collision rule, against the same real data.

    python run_v3.py --output-dir ../run-v3

The one change from round 1 (run.py): only matter that moves freely between collisions
feeds the heat term. Stars and galaxies do; gas and plasma do not, because their particles
change direction long before they travel a companion wavelength (Dicke narrowing; see
dicke_toy.py). Gas still counts in full as ordinary mass (it is in g_N) and adds coherently,
exactly like cold matter.

    g = g_N + exp(-|g_N|/g_d) sqrt(a (|g_N| + S_free)),   S_free = G sum k m_free / d^2,
    k = 3 sigma^2 / u^2 for free-streaming matter, 0 for colliding matter.

Galaxies: 149 SPARC rotation curves (the hot part is the bulge, as in round 1).
Clusters: 12 X-COP clusters; the hot part is now the stars (measured profiles for seven,
Ghizzardi et al. 2021; the median star/gas ratio for the other five). The stars' random
speed comes from the Jeans equation in the measured X-ray acceleration; a self-consistent
variant uses the law's own acceleration instead, and a third simply sets it to the gas's.
Three universal constants again; nothing is fitted per object.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run as R          # round-1 loaders, galaxy model, lenses, benchmarks
import law as L
import formula_guard as FG

G, K = L.G, L.KMS2_PER_KPC
XR = R.XR


def loglog(x, xp, yp):
    """log-log interpolation with power-law extrapolation at both ends."""
    lx, lxp, lyp = np.log(x), np.log(xp), np.log(yp)
    y = np.interp(lx, lxp, lyp)
    lo, hi = lx < lxp[0], lx > lxp[-1]
    if lo.any():
        s = (lyp[3] - lyp[0]) / (lxp[3] - lxp[0]); y[lo] = lyp[0] + s * (lx[lo] - lxp[0])
    if hi.any():
        s = (lyp[-1] - lyp[-4]) / (lxp[-1] - lxp[-4]); y[hi] = lyp[-1] + s * (lx[hi] - lxp[-1])
    return np.exp(y)


def prepare_clusters():
    """Round-1 cluster inputs plus the stars' Jeans dispersion from the X-ray acceleration."""
    import json as _j
    raw = _j.loads((R.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    cls = R.load_xcop()
    for c in cls:
        s = c['s']; ds = np.gradient(s)
        hm = raw[c['name']]['hydro_mass']
        g_obs = G * loglog(s, np.array(hm['RADIUS']), np.array(hm['M_FORW'])) / s ** 2
        c['g_obs_grid'] = g_obs
        c['sig2_star_hse'] = jeans_sigma2(s, ds, c['dms'], g_obs)
        # coarse grid for the self-consistent variant
        blk = 10; n = len(s) // blk * blk
        c['sc'] = s[:n].reshape(-1, blk).mean(1)
        c['dms_c'] = c['dms'][:n].reshape(-1, blk).sum(1)
        c['dmg_c'] = c['dmg'][:n].reshape(-1, blk).sum(1)
        c['Wc'] = L.shell_weights(c['sc'], c['sc'])
        c['Wk_c'] = L.shell_weights(c['Rk'], c['sc'])
        c['Mb_c'] = np.cumsum(c['dms_c'] + c['dmg_c'])
        c['g_obs_c'] = G * loglog(c['sc'], np.array(hm['RADIUS']), np.array(hm['M_FORW'])) / c['sc'] ** 2
    return cls


def jeans_sigma2(s, ds, dm, g):
    """Isotropic Jeans equation: sigma^2(r) = (1/rho) * integral_r^inf rho g ds'."""
    rho = dm / (4 * np.pi * s ** 2 * ds)
    tail = np.cumsum((rho * g * ds)[::-1])[::-1]
    with np.errstate(divide='ignore', invalid='ignore'):
        sig2 = np.where(rho > 1e-30, tail / rho, 0.0)
    return sig2


def cluster_M3(c, a, u, lam, sigma='hse', gas_weight=0.0, stars_weight=1.0):
    """Predicted enclosed mass at the six radii under the collision rule."""
    gN = G * c['Mb'] / c['Rk'] ** 2
    if sigma == 'sc':                       # self-consistent dispersion lives on the coarse grid
        S = G * (c['Wk_c'] @ (stars_weight * L.heat_weight(np.sqrt(c['sig2_star_sc_c']), u) * c['dms_c'])) / c['Rk'] ** 2
        return L.total(gN, S, a, lam) * c['Rk'] ** 2 / G
    sig2_star = c['sig2_star_hse'] if sigma == 'hse' else c['sig2']
    src = stars_weight * L.heat_weight(np.sqrt(sig2_star), u) * c['dms']
    if gas_weight:
        src = src + gas_weight * L.heat_weight(np.sqrt(c['sig2']), u) * c['dmg']
    S = G * (c['W'] @ src) / c['Rk'] ** 2
    if L.HOT_GEOMETRY == 'stream':          # round 19: the cold glow heard through the absorbing stream too
        heard = G * (L.stream_cold_weights(c['Rk'], c['s']) @ (c['dmg'] + c['dms'])) / c['Rk'] ** 2
        return L.total_heard(gN, heard, S, a, lam) * c['Rk'] ** 2 / G
    return L.total(gN, S, a, lam) * c['Rk'] ** 2 / G


def self_consistent_sigma(c, a, u, lam, iters=200, tol=1e-6):
    """Stars' dispersion from the Jeans equation in the law's OWN acceleration (baryons
    only), iterated to a fixed point on the coarse grid.  Started from the X-ray estimate."""
    sc = c['sc']; ds_c = np.gradient(sc)
    sig2_c = jeans_sigma2(sc, ds_c, c['dms_c'], c['g_obs_c'])
    gN = G * c['Mb_c'] / sc ** 2
    for i in range(iters):
        k = L.heat_weight(np.sqrt(sig2_c), u)
        S = G * (c['Wc'] @ (k * c['dms_c'])) / sc ** 2
        g = L.total(gN, S, a, lam)
        new = jeans_sigma2(sc, ds_c, c['dms_c'], g)
        tail_old, tail_new = sig2_c * c['dms_c'], new * c['dms_c']
        if np.max(np.abs(tail_new - tail_old)) <= tol * np.max(np.abs(tail_new)):
            sig2_c = new; break
        sig2_c = 0.5 * sig2_c + 0.5 * new
    c['sig2_star_sc_c'] = sig2_c; c['sc_iterations'] = i + 1
    c['g_sc_over_obs_c'] = g / c['g_obs_c']
    return sig2_c


def resid(cls, fn):
    return np.array([np.log(c['Mh'] / fn(c)) for c in cls])


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


def fit_u3(cls, a, lam, **kw):
    r = minimize_scalar(lambda lu: np.mean(resid(cls, lambda c: cluster_M3(c, a, 10 ** lu, lam, **kw)) ** 2),
                        bounds=(1.3, 4.0), method='bounded')
    return 10 ** r.x


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    gals = R.load_sparc(); cls = prepare_clusters()
    v1 = json.loads((HERE.parent / 'run-v1/results.json').read_text())
    print(f'{len(gals)} galaxies; {len(cls)} clusters ({sum(c["own_stars"] for c in cls)} with measured stellar profiles)', flush=True)

    # --- stars' dispersion, for the record (Jeans in the X-ray acceleration)
    disp = []
    for c in cls:
        rr = np.array([30., 100., 300.]) ; row = dict(name=c['name'])
        for r0 in rr:
            wdw = (c['s'] > r0 / 1.2) & (c['s'] < r0 * 1.2)
            row[f'sigma_star_{int(r0)}kpc'] = float(np.sqrt(np.sum((c['sig2_star_hse'] * c['dms'])[wdw]) / max(np.sum(c['dms'][wdw]), 1e-30)))
            row[f'sigma_gas_{int(r0)}kpc'] = float(np.sqrt(np.interp(r0, c['s'], c['sig2'])))
        disp.append(row)

    # --- constants: (a, lambda) on galaxies, u on clusters, alternated
    a, lam = v1['constants']['a_code'], v1['constants']['lam']
    u = fit_u3(cls, a, lam)
    for _ in range(4):
        a, lam = R.fit_a_lam(gals, u, (a, lam))
        u = fit_u3(cls, a, lam)
    print(f'constants: a = {a * K:.4e} m/s^2; lambda = {lam:.3f} (g_d = {lam * a * K:.3e}); u = {u:.1f} km/s', flush=True)

    # --- galaxies
    ours = lambda g: R.galaxy_g(g, a, u, lam)
    a_mond = v1['mond_context_a_SI'] / K
    mond = lambda g: L.mond_simple(g['gN'], a_mond)
    splits = (None, 'train', 'validation', 'test')
    gal = {'ours (round 3, 3 constants)': {str(s): R.sparc_score(gals, ours, s)[0] for s in splits},
           'ours (round 1, for reference)': {str(k): v for k, v in v1['galaxies']['hot companion (ours, 3 universal constants)'].items()},
           'MOND simple (context)': {str(s): R.sparc_score(gals, mond, s)[0] for s in splits},
           'Newton (context)': {str(s): R.sparc_score(gals, lambda g: g['gN'], s)[0] for s in splits},
           'NFW dark matter (context, 2 per galaxy)': {'None': v1['galaxies']['NFW dark matter (context, 2 parameters per galaxy)']['null']
                                                       if 'null' in v1['galaxies']['NFW dark matter (context, 2 parameters per galaxy)'] else
                                                       list(v1['galaxies']['NFW dark matter (context, 2 parameters per galaxy)'].values())[0]}}
    # bulge-dominated galaxies: does a stronger heat term hurt them?
    bf = np.array([np.max(g['vb2'] / np.maximum(g['vbar2'], 1e-9)) for g in gals])
    bulgy = [g for g, f in zip(gals, bf) if f > 0.5]
    gal['bulge-dominated subset'] = dict(n=len(bulgy), ours=R.sparc_score(bulgy, ours)[0], mond=R.sparc_score(bulgy, mond)[0],
                                         ours_round1=R.sparc_score(bulgy, lambda g: R.galaxy_g(g, v1['constants']['a_code'], v1['constants']['u_kms'], v1['constants']['lam']))[0])

    # --- clusters
    res = {}
    R3 = resid(cls, lambda c: cluster_M3(c, a, u, lam))
    res['ours (round 3: stars carry the heat)'] = R3
    R1 = resid(cls, lambda c: R.cluster_M(c, v1['constants']['a_code'], v1['constants']['u_kms'], v1['constants']['lam']))
    res['ours (round 1: gas and stars carry the heat)'] = R1
    res['ours without any heat'] = resid(cls, lambda c: R.cluster_M(c, a, None, lam))
    res['MOND simple'] = resid(cls, lambda c: L.mond_simple(G * c['Mb'] / c['Rk'] ** 2, a_mond) * c['Rk'] ** 2 / G)
    res['Newton'] = resid(cls, lambda c: c['Mb'])
    res['NFW dark matter (release fit, 2 per cluster)'] = resid(cls, lambda c: c['Mnfw'])
    # variants of the stars' dispersion
    u_gas = fit_u3(cls, a, lam, sigma='gas')
    res[f'variant: stars at the gas dispersion (u refit {u_gas:.0f})'] = resid(cls, lambda c: cluster_M3(c, a, u_gas, lam, sigma='gas'))
    u_sc = u
    for _ in range(12):                      # damped fixed point: dispersion depends on u
        for c in cls: self_consistent_sigma(c, a, u_sc, lam)
        u_new = fit_u3(cls, a, lam, sigma='sc')
        if abs(u_new / u_sc - 1) < 2e-3: u_sc = u_new; break
        u_sc = np.sqrt(u_sc * u_new)
    for c in cls: self_consistent_sigma(c, a, u_sc, lam)
    res[f'variant: self-consistent Jeans, law only (u refit {u_sc:.0f})'] = resid(cls, lambda c: cluster_M3(c, a, u_sc, lam, sigma='sc'))
    own = [c for c in cls if c['own_stars']]
    u_own = fit_u3(own, a, lam)
    res_own = resid(own, lambda c: cluster_M3(c, a, u_own, lam))
    # what do the data say about a residual gas term? fit a gas weight alongside u
    best = None
    for fg in (0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
        uu = fit_u3(cls, a, lam, gas_weight=fg)
        rr = rms(resid(cls, lambda c: cluster_M3(c, a, uu, lam, gas_weight=fg)))
        best = (best or []) + [dict(gas_weight=fg, u=float(uu), rms=rr)]
    clus = {k: dict(mean_by_radius=v.mean(0).tolist(), rms=rms(v)) for k, v in res.items()}
    clus['own stellar profiles only (7), u refit'] = dict(u=float(u_own), rms=rms(res_own), mean_by_radius=res_own.mean(0).tolist())
    clus['gas weight scan'] = best
    clus['per_cluster_round3'] = [dict(name=c['name'], own_stars=c['own_stars'], ln_MHSE_over_Mpred=r.tolist()) for c, r in zip(cls, R3)]
    # cross-validation of u (924 half splits)
    cv_u, cv_rms = [], []
    for tr in itertools.combinations(range(len(cls)), len(cls) // 2):
        te = [i for i in range(len(cls)) if i not in tr]
        uu = fit_u3([cls[i] for i in tr], a, lam); cv_u.append(uu)
        cv_rms.append(rms(resid([cls[i] for i in te], lambda c: cluster_M3(c, a, uu, lam))))
    clus['cross_validation'] = dict(n_splits=len(cv_u), u_median=float(np.median(cv_u)),
                                    u_5_95=[float(np.percentile(cv_u, 5)), float(np.percentile(cv_u, 95))],
                                    heldout_rms_median=float(np.median(cv_rms)),
                                    heldout_rms_5_95=[float(np.percentile(cv_rms, 5)), float(np.percentile(cv_rms, 95))])

    # --- lenses (stars are collisionless: same rule as round 1, new constants)
    print('lenses...', flush=True)
    lenses = R.lens_tests(a, u, lam)

    # --- KiDS: early types' stars are hot and collisionless; late types' disks are cold
    kids = {f'sigma {s} km/s, all baryons in hot stars': float(0.5 * np.log10(1 + L.heat_weight(s, u))) for s in (120, 150, 200, 250)}

    # --- the guard
    guards = [FG.guard('hot companion, round 3', point_mass=L.point_mass(a, u, lam),
                       sparc=[(g['gN'], ours(g)) for g in gals], a0_si=a * K)]

    # --- Solar System and energy budget
    GM_SUN, AU = 1.32712440018e20, 1.495978707e11
    solar = {pl: float(np.exp(-(GM_SUN / (d * AU) ** 2) / (lam * a * K)) * np.sqrt(a * K * GM_SUN / (d * AU) ** 2))
             for pl, d in (('Mercury', 0.387), ('Earth', 1.0), ('Saturn', 9.58), ('Neptune', 30.1))}
    ell = a * K * u * 1e3 / 2                       # W/kg, a = 2 l / u
    c_si, msun, yr = 2.998e8, 1.98847e30, 3.15576e7
    energy = dict(emission_W_per_kg=ell,
                  sun_power_W=ell * msun, sun_power_over_Lsun=ell * msun / 3.828e26,
                  sun_extra_mass_loss_per_year=ell / c_si ** 2 * yr,
                  note='the Sun is collisional plasma: no heat term, so only the base emission counts')

    payload = dict(experiment='hot-companion round 3: only free-streaming matter feeds the heat term',
                   constants=dict(a_SI=a * K, a_code=a, lam=lam, g_d_SI=lam * a * K, u_kms=u),
                   stars_dispersion_kms=disp, galaxies=gal, clusters=clus, lenses=lenses,
                   kids_offset_dex=kids, guard=guards, solar_system_anomalous_pull_SI=solar, energy=energy,
                   seconds=time.monotonic() - t0,
                   code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                (HERE / 'run_v3.py', HERE / 'run.py', HERE / 'law.py')})
    (out / 'results.json').write_text(json.dumps(payload, indent=2, default=float) + '\n')

    print('\nstars vs gas dispersion (km/s) at 30 / 100 / 300 kpc')
    for d in disp:
        print(f"  {d['name']:8s} stars {d['sigma_star_30kpc']:5.0f} {d['sigma_star_100kpc']:5.0f} {d['sigma_star_300kpc']:5.0f}"
              f"   gas {d['sigma_gas_30kpc']:5.0f} {d['sigma_gas_100kpc']:5.0f} {d['sigma_gas_300kpc']:5.0f}")
    print('\nGALAXIES mean velocity error (km/s): all / train / validation / test')
    for k, v in gal.items():
        print(f'  {k:42s} {v}')
    print('\nCLUSTERS ln(M_HSE/M_pred) at', XR, 'R500')
    for k, v in clus.items():
        if isinstance(v, dict) and 'mean_by_radius' in v:
            print(f"  {k:62s} {np.round(v['mean_by_radius'], 2)}  rms {v['rms']:.3f}")
    print('  gas weight scan:', [(b['gas_weight'], round(b['u']), round(b['rms'], 3)) for b in clus['gas weight scan']])
    cv = clus['cross_validation']
    print(f"  cross-validated u {cv['u_median']:.0f} ({cv['u_5_95'][0]:.0f}-{cv['u_5_95'][1]:.0f}); held-out rms {cv['heldout_rms_median']:.3f} ({cv['heldout_rms_5_95'][0]:.3f}-{cv['heldout_rms_5_95'][1]:.3f})")
    print('\nLENSES stellar-mass offset needed (dex): cold | hot;  in errors: cold | hot')
    for l in lenses:
        print(f"  {l['name']}  {l['dm_needed_cold_companion']:+.3f} | {l['dm_needed_hot_companion']:+.3f}   "
              f"{l['dm_needed_sigma_cold_companion']:+.2f} | {l['dm_needed_sigma_hot_companion']:+.2f}  (sigma {l['sigma_kms']:.0f}, k {l['heat_weight']:.3f})")
    print('\nKiDS offset (dex):', {k: round(v, 3) for k, v in kids.items()}, ' observed >= 0.2')
    print('guard:', guards[0]['verdict'])
    print('Solar System pull (m/s^2):', {k: f'{v:.1e}' for k, v in solar.items()})
    print('energy:', {k: (f'{v:.3e}' if isinstance(v, float) else v) for k, v in energy.items()})
    print(f'{time.monotonic() - t0:.0f} s')


if __name__ == '__main__':
    main()
