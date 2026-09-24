#!/usr/bin/env python3
"""Hot-companion gravity against real data, with the standing formula guard.

    python run.py --output-dir ../run-v1

Galaxies: 149 SPARC rotation curves (3,150+ points).  Clusters: 12 X-COP hydrostatic
mass profiles with measured gas density and temperature.  Lenses: 6 SLACS/KCWI strong
lenses with measured Einstein radii and resolved kinematics.  Plus the KiDS-1000
early/late-type lensing offset (Brouwer et al. 2021) as an out-of-sample check.

Three constants: a and lambda (SPARC) and u (X-COP). Nothing is fitted per object.
"""
from __future__ import annotations
import argparse, hashlib, io, json, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar, minimize, least_squares, brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / 'research_work/tools'))
sys.path.insert(0, str(ROOT / 'research_work/results/cross-prediction-response/code'))
import law as L
import formula_guard as FG

G, K = L.G, L.KMS2_PER_KPC
KPC_CM = 3.0856775814913673e21; MP = 1.67262192e-24; MSUN = 1.98847e33
MU, MUE, MP_KEV = 0.6, 1.14, 938272.088
XR = np.array([0.1, 0.2, 0.35, 0.5, 0.7, 1.0])          # cluster radii in R500


# ----------------------------------------------------------------------------- data
def load_sparc(alpha=None, h0_catalogue=73.0):
    """The 149 SPARC galaxies. With alpha (per Mpc) the Hubble-flow galaxies (f_D = 1 in Lelli et al. 2016,
    D = cz/H0 with H0 = 73) take the static law's distance D = ln(1 + z)/alpha instead, D x 73/(alpha c):
    radii scale with D, the mass models' speeds with sqrt(D) (round 12). Other distances are kept."""
    splits = json.loads((ROOT / 'companion_wave_test/data/sparc_frozen.json').read_text())['split']
    cat = {}
    for line in (ROOT / 'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
        f = line.split()
        if len(f) == 19:
            try: cat[f[0]] = dict(L=float(f[7]) * 1e9, MHI=float(f[13]) * 1e9, Rd=float(f[11]), f_D=int(f[4]))
            except ValueError: pass
    gals = []
    with zipfile.ZipFile(ROOT / 'temporal_candidate_audit/data/Rotmod_LTG.zip') as z:
        for split, names in splits.items():
            for n in names:
                s = h0_catalogue / (alpha * 299792.458) if (alpha is not None and cat[n]['f_D'] == 1) else 1.0
                a = np.atleast_2d(np.loadtxt(io.BytesIO(z.read(n + '_rotmod.dat')))); a = a[a[:, 0] > 0]
                r, v, e, vg, vd, vb = a[:, 0] * s, a[:, 1], a[:, 2], a[:, 3] * np.sqrt(s), a[:, 4] * np.sqrt(s), a[:, 5] * np.sqrt(s)
                vb2 = 0.7 * vb ** 2
                gN = (vg * np.abs(vg) + 0.5 * vd ** 2 + vb2) / r
                # bulge (the only hot component in disk galaxies): spherical mass profile
                sfine = np.geomspace(r.min() / 20, r.max(), 400)
                Mb = np.maximum.accumulate(np.interp(sfine, r, vb2 * r / G, left=0.))
                dmb = np.diff(np.r_[0., Mb])
                sigb = 0.65 * np.sqrt(vb2.max()) if vb2.max() > 0 else 0.
                ok = gN > 0
                Msys = (0.5 * cat[n]['L'] + 1.33 * cat[n]['MHI']) * s ** 2
                gals.append(dict(name=n, split=split, r=r[ok], v=v[ok], err=e[ok], gN=gN[ok],
                                 sfine=sfine, dmb=dmb, sigb=sigb, Msys=Msys, Rd=cat[n]['Rd'] * s,
                                 vb2=vb2[ok], vbar2=(vg * np.abs(vg) + 0.5 * vd ** 2 + vb2)[ok], f_D=cat[n]['f_D']))
    assert len(gals) == 149
    return gals


def logi(x, xp, yp):
    return np.exp(np.interp(np.log(x), np.log(xp), np.log(yp)))


def load_xcop(rmax=3.0):
    D = json.loads((ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    XS = np.geomspace(.05, 1.5, 30); fr = []
    for c in D.values():
        if 'stellar_mass' in c:
            R5 = c['header']['R500_kpc']
            fr.append(logi(XS * R5, c['stellar_mass']['radius_kpc'], c['stellar_mass']['Mstar']) /
                      logi(XS, c['gas_mass']['RADIUS'], c['gas_mass']['MGAS']))
    FRAC = np.median(fr, axis=0)
    out = []
    for n in sorted(D):
        c = D[n]; R5 = c['header']['R500_kpc']; T5 = c['temperature']['T500_keV']
        de = c['density']; rm = np.sqrt(np.array(de['r_in_kpc']) * np.array(de['r_out_kpc']))
        ne = np.array(de['ne_cm3']); ok = ne > 0; rm, ne = rm[ok], ne[ok]
        s = np.geomspace(rm[0], rmax * R5, 3000); ds = np.gradient(s)
        slope = np.polyfit(np.log(rm[-6:]), np.log(ne[-6:]), 1)[0]
        nes = np.where(s <= rm[-1], np.exp(np.interp(np.log(s), np.log(rm), np.log(ne))), ne[-1] * (s / rm[-1]) ** slope)
        dmg = 4 * np.pi * s ** 2 * MUE * MP * nes * KPC_CM ** 3 / MSUN * ds
        Mg = np.cumsum(dmg)
        own_stars = 'stellar_mass' in c
        Ms = logi(s, c['stellar_mass']['radius_kpc'], c['stellar_mass']['Mstar']) if own_stars else \
            Mg * np.exp(np.interp(np.log(s / R5), np.log(XS), np.log(FRAC)))
        dms = np.maximum(np.gradient(Ms, s) * ds, 0)
        tx, ts = c['temperature']['xray'], c['temperature']['xsz']
        tr = np.r_[tx['r_over_R500'], [x for x in ts['r_over_R500'] if x > max(tx['r_over_R500'])]] * R5
        tt = np.r_[tx['T_over_T500'], [t for x, t in zip(ts['r_over_R500'], ts['T_over_T500'])
                                       if x > max(tx['r_over_R500'])]] * T5
        o = np.argsort(tr); T = np.interp(np.log(s), np.log(tr[o]), tt[o])
        sig2 = T / (MU * MP_KEV) * L.C_KMS ** 2
        hm = c['hydro_mass']; Rk = XR * R5
        W = L.shell_weights(Rk, s)
        out.append(dict(name=n, R5=R5, T5=T5, s=s, dmg=dmg, dms=dms, sig2=sig2, W=W, Rk=Rk,
                        Mb=np.array([np.sum((dmg + dms)[s < R]) for R in Rk]),
                        Mh=logi(Rk, hm['RADIUS'], hm['M_FORW']),
                        eMh=np.interp(Rk, hm['RADIUS'], hm['EM_FORW']),
                        Mnfw=logi(Rk, hm['RADIUS'], hm['M_NFW']), own_stars=own_stars))
    return out


# ----------------------------------------------------------------------------- law on data
def galaxy_g(g, a, u, lam=np.inf):
    gN = g['gN']
    if u is None or g['sigb'] == 0:
        S = 0.
    else:
        S = L.scalar_sum(g['r'], g['sfine'], g['dmb'], L.heat_weight(g['sigb'], u))
    return L.total(gN, S, a, lam)


def sparc_score(gals, fn, split=None):
    rm, ll = [], []
    for g in gals:
        if split and g['split'] != split: continue
        p = fn(g); rm.append(np.sqrt(np.mean((np.sqrt(g['r'] * p) - g['v']) ** 2)))
        ll.append(np.mean(np.log(g['v'] ** 2 / g['r'] / p) ** 2))
    return float(np.mean(rm)), float(np.mean(ll))


def fit_a(gals, make_fn):
    r = minimize_scalar(lambda la: sparc_score(gals, make_fn(10 ** la))[1], bounds=(0.5, 5.5), method='bounded')
    return 10 ** r.x


def cluster_M(cl, a, u, lam=np.inf, stars_hot=True):
    gN = G * cl['Mb'] / cl['Rk'] ** 2
    if u is None:
        S = 0.
    else:
        k = L.heat_weight(np.sqrt(cl['sig2']), u)
        S = G * (cl['W'] @ (k * (cl['dmg'] + (cl['dms'] if stars_hot else 0)))) / cl['Rk'] ** 2
    return L.total(gN, S, a, lam) * cl['Rk'] ** 2 / G


def cluster_resid(cls, fn):
    return np.array([np.log(c['Mh'] / fn(c)) for c in cls])


def fit_u(cls, a, lam=np.inf):
    r = minimize_scalar(lambda lu: np.mean(cluster_resid(cls, lambda c: cluster_M(c, a, 10 ** lu, lam)) ** 2),
                        bounds=(2, 4.5), method='bounded')
    return 10 ** r.x


def fit_a_lam(gals, u, start):
    cost = lambda p: sparc_score(gals, lambda g: galaxy_g(g, 10 ** p[0], u, 10 ** p[1]))[1]
    best = min((minimize(cost, [np.log10(start[0]), l0], method='Nelder-Mead',
                         options=dict(xatol=1e-5, fatol=1e-9, maxiter=4000)) for l0 in (0.3, 0.7)),
               key=lambda r: r.fun)
    return 10 ** best.x[0], 10 ** best.x[1]


# ----------------------------------------------------------------------------- lenses
def lens_tests(a, u, lam=np.inf):
    import model as M
    out = []
    for name in M.LENSES:
        lens = M.make_lens(name)
        r = lens.r; dm_star = lens.Mstar * np.diff(np.r_[0., lens.frac])
        sigma = float(np.mean(lens.y))                       # stars' own dispersion
        Wr = L.shell_weights(r, r)
        S_unit = G * (Wr @ dm_star) / r ** 2                 # scalar sum of the stars
        gN_unit = G * lens.Mstar * lens.frac / r ** 2

        def g_on_r(dm, uu):
            f = 10 ** dm; k = 0. if uu is None else L.heat_weight(sigma, uu)
            return L.total(f * gN_unit, f * k * S_unit, a, lam)

        def theta(dm, uu):
            gr = g_on_r(dm, uu)
            gi = lambda x: np.exp(np.interp(np.log(x), np.log(r), np.log(gr)))
            def f(b):
                rr = b / np.cos(lens.lens_t)
                return lens.ratio * 4 / L.C_KMS ** 2 * np.dot(lens.lens_w, gi(rr) * rr) - b / lens.Dl
            grid = lens.Re * np.geomspace(1e-4, 1e3, 300); vals = np.array([f(b) for b in grid])
            ii = np.flatnonzero((vals[:-1] > 0) & (vals[1:] <= 0))
            if not len(ii): return 0.
            i = ii[-1]; return float(brentq(f, grid[i], grid[i + 1], xtol=1e-10) / lens.Dl * 206264.80624709636)

        def newton_theta(dm):
            gr = 10 ** dm * gN_unit
            gi = lambda x: np.exp(np.interp(np.log(x), np.log(r), np.log(gr)))
            def f(b):
                rr = b / np.cos(lens.lens_t)
                return lens.ratio * 4 / L.C_KMS ** 2 * np.dot(lens.lens_w, gi(rr) * rr) - b / lens.Dl
            grid = lens.Re * np.geomspace(1e-4, 1e3, 300); vals = np.array([f(b) for b in grid])
            ii = np.flatnonzero((vals[:-1] > 0) & (vals[1:] <= 0))
            return float(brentq(f, grid[ii[-1]], grid[ii[-1] + 1]) / lens.Dl * 206264.80624709636) if len(ii) else 0.

        def vrms_chi2(dm, uu, beta=0.):
            v = lens.moments(g_on_r(dm, uu), {'beta': beta}, False)
            e = v - lens.y
            from scipy.linalg import solve_triangular
            w = solve_triangular(lens.chol, e, lower=True)
            return float(w @ w), v

        def dm_needed(uu):
            f = lambda dm: theta(dm, uu) - lens.theta
            try: return float(brentq(f, -1.5, 1.5, xtol=1e-8))
            except ValueError: return float('nan')

        row = dict(name=name, sigma_kms=sigma, theta_obs=float(lens.theta),
                   mass_log_error=float(lens.mass_log_error),
                   heat_weight=float(L.heat_weight(sigma, u)))
        row['theta_newton'] = newton_theta(0.)
        for tag, uu in (('cold_companion', None), ('hot_companion', u)):
            row['theta_' + tag] = theta(0., uu)
            dm = dm_needed(uu); row['dm_needed_' + tag] = dm
            row['dm_needed_sigma_' + tag] = dm / row['mass_log_error'] if np.isfinite(dm) else float('nan')
            row['dm_needed_sigma_salpeter_' + tag] = (dm - 0.25) / row['mass_log_error'] if np.isfinite(dm) else float('nan')
            c2, v = vrms_chi2(0., uu); row['vrms_chi2_' + tag] = c2; row['vrms_' + tag] = v.tolist()
        row['vrms_obs'] = lens.y.tolist()
        out.append(row)
    return out


# ----------------------------------------------------------------------------- benchmarks
def nfw_fit_sparc(gals):
    rm = []
    for g in gals:
        def res(p):
            v2 = g['vbar2'] + L.nfw_velocity2(g['r'], p[0], p[1])
            return (np.sqrt(np.maximum(v2, 1e-6)) - g['v']) / g['err']
        best = min((least_squares(res, [v0, c0], bounds=([5, 1], [600, 60])) for v0 in (60, 150, 300) for c0 in (5, 15)),
                   key=lambda s: s.cost)
        v2 = g['vbar2'] + L.nfw_velocity2(g['r'], *best.x)
        rm.append(np.sqrt(np.mean((np.sqrt(np.maximum(v2, 0)) - g['v']) ** 2)))
    return float(np.mean(rm))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    gals = load_sparc(); cls = load_xcop()
    print(f'{len(gals)} galaxies, {sum(len(g["r"]) for g in gals)} points; {len(cls)} clusters', flush=True)

    # --- fit the three constants, alternating (a, lambda on galaxies; u on clusters)
    a = fit_a(gals, lambda aa: (lambda g: galaxy_g(g, aa, None)))
    u = fit_u(cls, a); lam = 3.0
    for _ in range(4):
        a, lam = fit_a_lam(gals, u, (a, lam))
        u = fit_u(cls, a, lam)
    print(f'constants: a = {a * K:.4e} m/s^2 ; lambda = {lam:.3f} (release scale {lam * a * K:.3e} m/s^2) ; u = {u:.1f} km/s', flush=True)

    # --- galaxies
    ours = lambda g: galaxy_g(g, a, u, lam)
    a_mond = fit_a(gals, lambda aa: (lambda g: L.mond_simple(g['gN'], aa)))
    mond = lambda g: L.mond_simple(g['gN'], a_mond)
    newton = lambda g: g['gN']
    gal_scores = {
        'hot companion (ours, 3 universal constants)': {s: sparc_score(gals, ours, s)[0] for s in (None, 'train', 'validation', 'test')},
        'MOND simple (context, 1 constant)': {s: sparc_score(gals, mond, s)[0] for s in (None, 'train', 'validation', 'test')},
        'Newton, visible matter only (context)': {s: sparc_score(gals, newton, s)[0] for s in (None, 'train', 'validation', 'test')},
    }
    print('fitting an NFW dark-matter halo to every galaxy (298 parameters)...', flush=True)
    gal_scores['NFW dark matter (context, 2 parameters per galaxy)'] = {None: nfw_fit_sparc(gals)}

    # --- clusters
    R_ours = cluster_resid(cls, lambda c: cluster_M(c, a, u, lam))
    R_cold = cluster_resid(cls, lambda c: cluster_M(c, a, None, lam))
    R_mond = cluster_resid(cls, lambda c: L.mond_simple(G * c['Mb'] / c['Rk'] ** 2, a_mond) * c['Rk'] ** 2 / G)
    R_newt = cluster_resid(cls, lambda c: c['Mb'])
    R_nfw = cluster_resid(cls, lambda c: c['Mnfw'])
    bias = 0.08 + 0.10 * XR                                  # context only, not used in any fit
    clus = dict(radii_R500=XR.tolist(),
                ours=dict(mean_by_radius=R_ours.mean(0).tolist(), rms=float(np.sqrt(np.mean(R_ours ** 2))),
                          mean_by_radius_with_context_bias=(R_ours - np.log(1 - bias)).mean(0).tolist()),
                ours_without_heat=dict(mean_by_radius=R_cold.mean(0).tolist(), rms=float(np.sqrt(np.mean(R_cold ** 2)))),
                mond_simple=dict(mean_by_radius=R_mond.mean(0).tolist(), rms=float(np.sqrt(np.mean(R_mond ** 2)))),
                newton=dict(mean_by_radius=R_newt.mean(0).tolist(), rms=float(np.sqrt(np.mean(R_newt ** 2)))),
                nfw_release_fit=dict(rms=float(np.sqrt(np.mean(R_nfw ** 2))), note='the X-COP release NFW fit, 2 parameters per cluster'),
                per_cluster=[dict(name=c['name'], T500_keV=c['T5'], ln_MHSE_over_Mpred=r.tolist(), own_stellar_profile=c['own_stars'])
                             for c, r in zip(cls, R_ours)])

    # --- robustness: MOND allowed its own cluster constant; cross-validation of u
    import itertools
    rm = minimize_scalar(lambda la: np.mean(cluster_resid(cls, lambda c: L.mond_simple(G * c['Mb'] / c['Rk'] ** 2, 10 ** la) * c['Rk'] ** 2 / G) ** 2),
                         bounds=(1, 6), method='bounded')
    a_mc = 10 ** rm.x
    R_mc = cluster_resid(cls, lambda c: L.mond_simple(G * c['Mb'] / c['Rk'] ** 2, a_mc) * c['Rk'] ** 2 / G)
    cv_u, cv_rms = [], []
    for tr in itertools.combinations(range(len(cls)), len(cls) // 2):
        te = [i for i in range(len(cls)) if i not in tr]
        uu = fit_u([cls[i] for i in tr], a, lam); cv_u.append(uu)
        Rt = cluster_resid([cls[i] for i in te], lambda c: cluster_M(c, a, uu, lam))
        cv_rms.append(float(np.sqrt(np.mean(Rt ** 2))))
    clus['mond_simple_refit_on_clusters'] = dict(a0_SI=a_mc * K, times_galaxy_value=a_mc / a_mond,
                                                 mean_by_radius=R_mc.mean(0).tolist(), rms=float(np.sqrt(np.mean(R_mc ** 2))))
    clus['cross_validation'] = dict(n_splits=len(cv_u), u_median=float(np.median(cv_u)),
                                    u_5_95=[float(np.percentile(cv_u, 5)), float(np.percentile(cv_u, 95))],
                                    heldout_rms_median=float(np.median(cv_rms)),
                                    heldout_rms_5_95=[float(np.percentile(cv_rms, 5)), float(np.percentile(cv_rms, 95))])

    # --- lenses
    print('lenses...', flush=True)
    lenses = lens_tests(a, u, lam)

    # --- KiDS-1000 early/late offset (Brouwer et al. 2021: >= 0.2 dex, early types higher)
    kids = {f'sigma {s} km/s': float(0.5 * np.log10(1 + L.heat_weight(s, u))) for s in (150, 200, 250, 300)}

    # --- the guard
    print('guard...', flush=True)
    guards = []
    guards.append(FG.guard('MOND simple (control: must be flagged)',
                           point_mass=lambda r, M, s: L.mond_simple(G * M / r ** 2, a_mond),
                           sparc=[(g['gN'], mond(g)) for g in gals], a0_si=a_mond * K))
    guards.append(FG.guard('hot companion (ours)', point_mass=L.point_mass(a, u, lam),
                           sparc=[(g['gN'], ours(g)) for g in gals], a0_si=a * K))
    # the two failed non-MOND candidates from this session, for the record
    pool = lambda g: g['gN'] * (1 + g['r'] / np.sqrt(G * g['Msys'] / a))
    guards.append(FG.guard('companion pool (failed on SPARC)',
                           point_mass=lambda r, M, s: (G * M / r ** 2) * (1 + r / np.sqrt(G * M / a)),
                           sparc=[(g['gN'], pool(g)) for g in gals], a0_si=a * K))
    for gr in guards: print('  ', gr['name'], '->', gr['verdict'], flush=True)

    # --- Solar System: anomalous pull at each planet from the released companion
    GM_SUN, AU = 1.32712440018e20, 1.495978707e11
    solar = {}
    for pl, d in (('Mercury', 0.387), ('Earth', 1.0), ('Jupiter', 5.20), ('Saturn', 9.58), ('Neptune', 30.1)):
        gS = GM_SUN / (d * AU) ** 2; aS = a * K
        solar[pl] = float(np.exp(-gS / (lam * aS)) * np.sqrt(aS * gS))
    payload = dict(experiment='hot-companion v1', constants=dict(a_SI=a * K, a_code=a, lam=lam,
                   release_scale_SI=lam * a * K, u_kms=u,
                   note='a and lambda fitted on 149 SPARC galaxies; u fitted on 12 X-COP clusters; alternated to convergence'),
                   solar_system_anomalous_pull_SI=solar,
                   mond_context_a_SI=a_mond * K, galaxies=gal_scores, clusters=clus, lenses=lenses,
                   kids_early_late_predicted_offset_dex=kids,
                   kids_observed='>= 0.2 dex (factor >= 1.5), early types higher, >= 5.7 sigma (Brouwer et al. 2021, A&A 650, A113)',
                   guard=guards, seconds=time.monotonic() - t0,
                   code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                (HERE / 'run.py', HERE / 'law.py', ROOT / 'research_work/tools/formula_guard.py')})
    (out / 'results.json').write_text(json.dumps(payload, indent=2, default=float) + '\n')

    print('\nGALAXIES (mean velocity error, km/s; all / train / validation / test)')
    for k, v in gal_scores.items():
        print(f'  {k:52s} ' + '  '.join(f'{x:6.2f}' for x in v.values()))
    print('\nCLUSTERS ln(M_HSE / M_pred) mean at', XR, 'R500')
    for k in ('ours', 'ours_without_heat', 'mond_simple', 'mond_simple_refit_on_clusters', 'newton'):
        print(f"  {k:30s} {np.round(clus[k]['mean_by_radius'], 2)}  rms {clus[k]['rms']:.3f}")
    cv = clus['cross_validation']
    print(f"  cross-validated u: {cv['u_median']:.0f} km/s ({cv['u_5_95'][0]:.0f}-{cv['u_5_95'][1]:.0f}); held-out rms {cv['heldout_rms_median']:.3f} ({cv['heldout_rms_5_95'][0]:.3f}-{cv['heldout_rms_5_95'][1]:.3f}) over {cv['n_splits']} splits")
    print(f"  MOND refit on clusters needs a0 x{clus['mond_simple_refit_on_clusters']['times_galaxy_value']:.1f} its galaxy value")
    print(f"  ours + context bias  {np.round(clus['ours']['mean_by_radius_with_context_bias'], 2)}")
    print(f"  NFW (release fit, 24 params) rms {clus['nfw_release_fit']['rms']:.3f}")
    print('\nLENSES  theta_E obs | Newton | cold | hot (arcsec);  stellar-mass offset needed (in published errors) cold | hot;  V_rms chi2 cold | hot')
    for l in lenses:
        print(f"  {l['name']}  {l['theta_obs']:.3f} | {l['theta_newton']:.3f} | {l['theta_cold_companion']:.3f} | {l['theta_hot_companion']:.3f}"
              f"   {l['dm_needed_sigma_cold_companion']:+.2f} | {l['dm_needed_sigma_hot_companion']:+.2f}"
              f"  [Salpeter {l['dm_needed_sigma_salpeter_cold_companion']:+.2f} | {l['dm_needed_sigma_salpeter_hot_companion']:+.2f}]"
              f"   {l['vrms_chi2_cold_companion']:.1f} | {l['vrms_chi2_hot_companion']:.1f}   (sigma {l['sigma_kms']:.0f}, k {l['heat_weight']:.3f})")
    print('\nSolar System anomalous pull (m/s^2):', {k: f'{v:.1e}' for k, v in solar.items()})
    print('\nKiDS early/late predicted offset (dex):', {k: round(v, 3) for k, v in kids.items()}, ' observed >= 0.2')
    print(f'\n{time.monotonic() - t0:.0f} s')


if __name__ == '__main__':
    main()
