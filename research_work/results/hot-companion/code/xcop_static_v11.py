"""Round 11: the 12 X-COP clusters in the project's static distance law, and the companion speed u.

The X-COP release gives every profile in physical units for a flat LCDM with H0 = 70, Omega_m = 0.3.
At each cluster's redshift, at fixed angle and flux (collisions_v10.factors):
    radii and hydrostatic masses (M ~ kT r)       x  D_A(static) / D_A(LCDM)
    gas masses (from the X-ray emission measure)   x  (D_L ratio) (D_A ratio)^1.5
    star masses (from light)                       x  (D_L ratio)^2
    temperatures, and the stars' Jeans speeds       unchanged (sigma^2 ~ g r ~ M / r)
The law's a and g_d are held at their SPARC values (SPARC's distances are nearby: our law's scale
corresponds to a Hubble-constant-like 74.6, against SPARC's 73 for its Hubble-flow galaxies, a 2% change),
and u is refitted.

    python code/xcop_static_v11.py --output-dir run-xcop-static-v11
"""
from __future__ import annotations
import argparse, copy, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import run as RUN                                 # noqa: E402
import run_v3 as R3                               # noqa: E402
import collisions_v10 as C10                      # noqa: E402

G = L.G


# ---------------------------------------------------------------------------- the stars, deprojected
# The release's cumulative stellar mass profiles (Ghizzardi et al. 2021, arXiv:2007.01084, Sect. 4.1) are masses
# inside a PROJECTED radius (a cylinder). Their Table 3 gives the mass inside the R500 sphere as 0.75 x the
# projected value, from a gNFW galaxy distribution (c = 0.72, alpha = 1.64; van der Burg et al. 2015). Rounds 1-10
# used the projected profiles as if they were spherical, which overstates the stars by 1.3 at R500 and more inside.
# Here each profile is fitted with a BCG (Hernquist, a = 15 kpc) plus satellites (that gNFW, R200 = R500 / 0.65,
# truncated along the line of sight at R_t, chosen so that the satellites' sphere/cylinder ratio at R500 is the
# paper's 0.75), and the spherical profile of the fit replaces the projected one.
R500_OVER_R200 = 0.65


def _shapes(R500, r):
    R200 = R500 / R500_OVER_R200
    rs = R200 / 0.72
    x = r / rs
    sat = x ** -1.64 * (1 + x) ** (1.64 - 3)
    a = 15.0
    bcg = a / (2 * np.pi * r * (r + a) ** 3)
    return sat, bcg, R200


def _proj_and_3d(rho, r, dr, R, Rt):
    """cumulative projected mass inside radii R and spherical mass inside radii R, for rho on the grid r (cut at Rt)."""
    dm = 4 * np.pi * r ** 2 * rho * dr * (r <= Rt)
    m3 = np.array([dm[r <= Ri].sum() for Ri in R])
    w = np.where(r[None, :] <= R[:, None], 1.0, 1 - np.sqrt(np.clip(1 - (R[:, None] / r[None, :]) ** 2, 0, 1)))
    return (w * dm[None, :]).sum(1), m3


def deproject_profile(radius_kpc, Mstar_proj, R500):
    """Fit BCG + satellites to a projected cumulative stellar mass profile; return a function M3D(<r) and the fit."""
    from scipy.optimize import brentq, least_squares
    r = np.geomspace(0.05, 6 * R500 / R500_OVER_R200, 6000); dr = np.gradient(r)
    sat, bcg, R200 = _shapes(R500, r)
    # line-of-sight cut-off: the satellites' sphere/cylinder ratio at R500 equals the paper's 0.75
    def ratio(Rt):
        p, m3 = _proj_and_3d(sat, r, dr, np.array([R500]), Rt)
        return m3[0] / p[0] - 0.75
    Rt = brentq(ratio, 1.02 * R500, 6 * R200) if ratio(6 * R200) < 0 else 6 * R200
    R = np.asarray(radius_kpc, float); M = np.asarray(Mstar_proj, float)
    ok = (R >= 10.0) & (M > 0)
    ps, m3s = _proj_and_3d(sat, r, dr, R[ok], Rt); pb, m3b = _proj_and_3d(bcg, r, dr, R[ok], Rt)
    fit = least_squares(lambda q: np.log(np.exp(q[0]) * pb + np.exp(q[1]) * ps) - np.log(M[ok]),
                        np.log([0.3 * M[ok].max() / max(pb[-1], 1e-30), 0.7 * M[ok].max() / max(ps[-1], 1e-30)]))
    Mb, Ms = np.exp(fit.x)
    cum_s = np.cumsum(4 * np.pi * r ** 2 * sat * dr * (r <= Rt)); cum_b = np.cumsum(4 * np.pi * r ** 2 * bcg * dr * (r <= Rt))
    M3 = lambda rr: Mb * np.interp(rr, r, cum_b) + Ms * np.interp(rr, r, cum_s)
    P5 = float(np.exp(np.interp(np.log(R500), np.log(R[R > 0]), np.log(np.maximum(M[R > 0], 1e-30)))))
    return M3, dict(R_t_over_R200=float(Rt / R200), M_bcg=float(Mb * cum_b[-1]), M_sat=float(Ms * cum_s[-1]),
                    rms_log_fit=float(np.sqrt(np.mean(fit.fun ** 2))), M3D_R500=float(M3(R500)), M2D_R500=P5,
                    sphere_over_cylinder_R500=float(M3(R500) / P5))


def deproject_xcop(cls, raw):
    """Replace each cluster's projected stellar profile by its spherical deprojection (the 7 clusters with
    profiles) or by the median spherical star/gas ratio of those 7 (the other 5), as run.load_xcop does."""
    import copy as _copy
    R = RUN                          # (a bare 'import run' can pick up another project's run.py once the SLACS code has loaded it)
    XSg = np.geomspace(.05, 1.5, 30)
    out, fits, fr = [], {}, []
    M3s = {}
    for c in cls:
        rc = raw[c['name']]
        if 'stellar_mass' in rc:
            M3, info = deproject_profile(rc['stellar_mass']['radius_kpc'], rc['stellar_mass']['Mstar'], rc['header']['R500_kpc'])
            M3s[c['name']] = M3; fits[c['name']] = info
            R5 = rc['header']['R500_kpc']
            Mg_at = R.logi(XSg, rc['gas_mass']['RADIUS'], rc['gas_mass']['MGAS'])
            fr.append(M3(XSg * R5) / Mg_at)
    FRAC = np.median(fr, axis=0)
    for c in cls:
        d = _copy.copy(c); s_ = c['s']; ds = np.gradient(s_)
        if c['name'] in M3s:
            Ms = M3s[c['name']](s_)
        else:
            Mg = np.cumsum(c['dmg'])
            Ms = Mg * np.exp(np.interp(np.log(s_ / c['R5']), np.log(XSg), np.log(FRAC)))
        d['dms'] = np.maximum(np.gradient(Ms, s_) * ds, 0)
        d['Mb'] = np.array([np.sum((d['dmg'] + d['dms'])[s_ < R_]) for R_ in c['Rk']])
        d['sig2_star_hse'] = R3.jeans_sigma2(s_, ds, d['dms'], c['g_obs_grid'])
        d['stars_deprojected'] = True
        out.append(d)
    return out, dict(fits=fits, median_star_over_gas_at=dict(zip([float(x) for x in XSg[[9, 18, 23, 29]]], [float(x) for x in FRAC[[9, 18, 23, 29]]])))


def to_static(c, z):
    f = C10.factors(z, 1.0)                       # the source redshift only matters for lensing, not used here
    fs, fg, fst = f['size'], f['gas'], f['stars']
    d = copy.copy(c)
    d['s'] = c['s'] * fs; d['Rk'] = c['Rk'] * fs; d['R5'] = c['R5'] * fs
    d['dmg'] = c['dmg'] * fg; d['dms'] = c['dms'] * fst
    d['Mb'] = np.array([np.sum((d['dmg'] + d['dms'])[d['s'] < R]) for R in d['Rk']])
    d['Mh'] = c['Mh'] * fs; d['eMh'] = c['eMh'] * fs
    d['W'] = L.shell_weights(d['Rk'], d['s'])     # depends on ratios only; recomputed for clarity
    # the stars' Jeans speeds from the observed (hydrostatic) pull: g ~ M/r^2 scales as 1/fs, r as fs
    d['g_obs_grid'] = c['g_obs_grid'] / fs
    d['sig2_star_hse'] = R3.jeans_sigma2(d['s'], np.gradient(d['s']), d['dms'], d['g_obs_grid'])
    d['factors'] = dict(size=fs, gas=fg, stars=fst, z=z)
    return d


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law('round9'); a, lam = law['a_code'], law['lam']
    raw = json.loads((RUN.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    ctx = C.Context(tier='quick', verbose=False)
    cls = ctx.xcop()
    dep, dep_info = deproject_xcop(cls, raw)
    zof = lambda c: raw[c['name']]['header']['z']
    samples = {'standard distances, stars as used in rounds 1-10 (projected)': cls,
               'static distances, stars as used in rounds 1-10 (projected)': [to_static(c, zof(c)) for c in cls],
               'standard distances, stars deprojected': dep,
               'static distances, stars deprojected': [to_static(c, zof(c)) for c in dep]}
    res = dict(experiment='round 11: X-COP in the static distance law, with the stars deprojected', law=dict(a_SI=law['a_SI'], lam=lam, u_kms=law['u_kms']),
               deprojection=dep_info, clusters=[dict(name=c['name'], **to_static(c, zof(c))['factors']) for c in cls])
    rows = {}
    for tag, sample in samples.items():
        u = R3.fit_u3(sample, a, lam)
        r = R3.resid(sample, lambda c: R3.cluster_M3(c, a, u, lam))
        rows[tag] = dict(u_best=float(u), rms=float(R3.rms(r)), mean_by_radius=[float(x) for x in r.mean(0)])
        print(f"{tag:62s} best u {u:6.1f} km/s, rms {R3.rms(r):.3f}; by radius " + ' '.join(f'{x:+.2f}' for x in r.mean(0)), flush=True)
    res['fits'] = rows
    st_sph = samples['static distances, stars deprojected']
    ratio = np.array([[c['Mb'][-1]] for c in st_sph])
    gals = RUN.load_sparc()
    for key, sample in (('joint_refit_static_projected_stars', samples['static distances, stars as used in rounds 1-10 (projected)']),
                        ('joint_refit_static', st_sph)):
        aa, ll, uu = a, lam, R3.fit_u3(sample, a, lam)
        log = []
        for it in range(3):
            aa, ll = RUN.fit_a_lam(gals, uu, (aa, ll))
            uu = R3.fit_u3(sample, aa, ll)
            sc = RUN.sparc_score(gals, lambda g: RUN.galaxy_g(g, aa, uu, ll))[0]
            r = R3.resid(sample, lambda c: R3.cluster_M3(c, aa, uu, ll))
            log.append(dict(iteration=it + 1, a_SI=float(aa * L.KMS2_PER_KPC), lam=float(ll), g_d_SI=float(aa * ll * L.KMS2_PER_KPC),
                            u_kms=float(uu), sparc_rms_kms=float(sc), xcop_rms=float(R3.rms(r))))
        res[key] = log
        print(f"{key}: a {aa * L.KMS2_PER_KPC:.4e} m/s2, g_d {aa * ll * L.KMS2_PER_KPC:.4e}, u {uu:.1f} km/s; SPARC {sc:.2f} km/s, X-COP rms {R3.rms(r):.3f}", flush=True)
    f = np.array([[d['factors']['size'], d['factors']['gas'], d['factors']['stars']] for d in st_sph])
    res['factor_ranges'] = dict(size=[float(f[:, 0].min()), float(f[:, 0].max())], gas=[float(f[:, 1].min()), float(f[:, 1].max())],
                                stars=[float(f[:, 2].min()), float(f[:, 2].max())])
    star_gas = [float(np.sum(c['dms'][c['s'] < c['R5']]) / np.sum(c['dmg'][c['s'] < c['R5']])) for c in dep]
    res['star_over_gas_R500'] = dict(deprojected=star_gas, projected=[float(np.sum(c['dms'][c['s'] < c['R5']]) / np.sum(c['dmg'][c['s'] < c['R5']])) for c in cls])
    print('star/gas at R500, deprojected: %.3f-%.3f (median %.3f); as used before: %.3f-%.3f' % (min(star_gas), max(star_gas), np.median(star_gas),
          min(res['star_over_gas_R500']['projected']), max(res['star_over_gas_R500']['projected'])))
    # star/gas inside spheres of 0.02-1 R500 (the seven clusters with optical data), before and after the deprojection
    withdata = [c['name'] for c in dep if raw[c['name']].get('stellar_mass')]
    res['star_over_gas_by_radius'] = {}
    for tag, sample in (('projected', cls), ('deprojected', dep)):
        rows_r = {}
        for x in (0.02, 0.05, 0.1, 1.0):
            v = [float(np.sum(c['dms'][c['s'] < x * c['R5']]) / np.sum(c['dmg'][c['s'] < x * c['R5']])) for c in sample if c['name'] in withdata]
            rows_r[f'{x:g} R500'] = [min(v), float(np.median(v)), max(v)]
        res['star_over_gas_by_radius'][tag] = rows_r
        print(f'star/gas ({tag}, {len(withdata)} clusters with optical data), min/median/max: ' +
              '; '.join(f"{k}: {v[0]:.3f}/{v[1]:.3f}/{v[2]:.3f}" for k, v in rows_r.items()), flush=True)
    # the adopted sample's u from each half of the clusters (924 splits), a and g_d from the joint refit
    j = res['joint_refit_static'][-1]; aj, lj = j['a_SI'] / L.KMS2_PER_KPC, j['lam']
    import itertools
    cv_u = [R3.fit_u3([st_sph[i] for i in tr], aj, lj) for tr in itertools.combinations(range(len(st_sph)), len(st_sph) // 2)]
    res['cross_validation_u'] = dict(n_splits=len(cv_u), u_p5=float(np.percentile(cv_u, 5)), u_median=float(np.median(cv_u)),
                                     u_p95=float(np.percentile(cv_u, 95)))
    print('u from half the clusters (%d splits): 5-95%% %.1f-%.1f km/s, median %.1f' % (len(cv_u), res['cross_validation_u']['u_p5'],
          res['cross_validation_u']['u_p95'], res['cross_validation_u']['u_median']), flush=True)
    res['seconds'] = time.monotonic() - t0
    (out / 'xcop_static_v11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'xcop_static_v11.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
