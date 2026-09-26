"""Round 25, DRT-1 against the data. Two of DRT-1's results (drt1-transport/README.md, main at ab1c027) turned into
numbers from the suite's own models:

1. Where the energy completion fails. DRT-1 section 7: the simplest energy completion of the repaired field equation
   has no lower bound where the companion's gradient |grad phi| exceeds a/2. In round systems |grad phi|^2 =
   a (|g_N| + S) (round 24), so the condition is |g_N| + S > a/4 = 1.57e-11 m/s^2. Share of each family's test points
   (run-frustration-v25/frustration_v25.json, the same points as the registered test) above it, and of the Bullet
   Cluster's aperture walls.
2. How much of its own companion cluster gas could keep. DRT-1 section 4: scattering does not pile an outside bath
   onto gas, but gas that scatters its own cold emission with optical depth tau keeps (1/2 + tau/12)/(1/2) of it, an
   excess eps = tau/6 of its transparent brightness. In the law's terms that adds eps G int rho_gas / d^2 to S (energy
   that does not flow). Bounded by the 12 X-COP clusters (u fixed and refitted) and by the Bullet Cluster's lensing
   (the leftover lensing on the gas, the peaks' distance from the galaxies, the masses), with the gas's own brightness
   kept around today's gas (the most dangerous case: fully built up and stopped with the gas).

    python code/drt1_bearing_v25.py --output run-frustration-v25/drt1_bearing_v25.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(RESULTS / 'regression'))
T0 = time.monotonic()
EPS = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0)                      # gas brightness kept, in units of its transparent brightness


def log(msg):
    print(f'[{time.monotonic() - T0:5.0f} s] {msg}', flush=True)


# ------------------------------------------------------------------------------------------------ 1. the threshold
def threshold_coverage(law, ctx):
    import law as L
    import mw_dwarfs_v7 as D
    import t_dwarfs as TD
    import t_milky_way as TMW
    import milky_way_v7 as MW
    import mw_model as M
    import kids_heat_v12 as KH
    import kids_static_v11 as KS
    import collisions_v10 as C10
    K = L.KMS2_PER_KPC
    thr = law['a_SI'] / 4                                    # |g_N| + S above this: |grad phi| > a/2
    P = json.loads((RESULTS / 'run-frustration-v25/frustration_v25.json').read_text())['points']
    out = {}
    tot = lambda gs: np.asarray(gs, float)
    sp = [p for p in P if p['family'] == 'sparc']
    v = tot([p['gN_SI'] * (1 + p['S_over_gN']) for p in sp])
    out['sparc'] = dict(n=len(sp), share_above=float(np.mean(v > thr)))
    xc = [p for p in P if p['family'] == 'xcop']
    v = tot([p['gN_SI'] * (1 + p['S_over_gN']) for p in xc])
    out['xcop'] = dict(n=len(xc), share_above=float(np.mean(v > thr)), min_over_threshold=float(v.min() / thr))
    # KiDS: point lenses, |g_N| + S = g_bar (1 + 0.85 k) (lensing_census_v7.pull_profile), reliable bins
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    tabs = KS.tables(); k = KH.heat(law['u_kms'])
    gb = tabs['all']['gbar'] * f['stars'] / f['size'] ** 2
    rel = tabs['all']['gbar'] >= 1e-13
    kk = dict(red=k['red'], blue=k['blue'], disc=k['disc'], bulge=k['bulge'])
    out['kids'] = {s: dict(n=int(rel.sum()), share_above=float(np.mean(gb[rel] * (1 + 0.85 * kv) > thr))) for s, kv in kk.items()}
    # the Milky Way: the rotation-curve radii, in the plane
    comps, grid, F = TMW.setup(ctx)
    ev = MW.Evaluator(comps, grid, F, law)
    parts = MW.models(comps)['M17']
    ev.run(parts, 'ours', reach=law['reach_kpc'])
    gR = sum(ff * F[kx][0] for kx, ff in parts.items()); gz = sum(ff * F[kx][1] for kx, ff in parts.items())
    hot = [MW.Scaled(comps[kx], ff) for kx, ff in parts.items() if kx in ('bulge', 'halo')]
    S, _, _ = M.heat_fields(hot, grid, law['u_kms'])
    from scipy.interpolate import RegularGridInterpolator as RGI
    fI = RGI((grid.R, grid.z), (np.hypot(gR, gz) + S) * K, bounds_error=False, fill_value=None)
    Rs = np.array([p['r_kpc'] for p in P if p['family'] == 'mw' and p['id'].startswith('mw.rc')])
    v = fI(np.c_[Rs, 0 * Rs])
    Rgrid = np.geomspace(1, 300, 400); vg = fI(np.c_[Rgrid, 0 * Rgrid])
    out['mw'] = dict(n=int(Rs.size), share_above=float(np.mean(v > thr)),
                     crossing_radius_kpc=float(Rgrid[np.argmax(vg < thr)]) if (vg < thr).any() else None)
    # dwarfs: at the half-light radius, averaged over directions
    data = json.loads((RESULTS / 'data/mw_dwarfs.json').read_text())
    mw = TD.galaxy_profile(law, ctx)
    mu, w = np.polynomial.legendre.leggauss(256)
    rows = []
    for d in data['dwarfs']:
        env = D.galaxy_env(d['D_gc_kpc'], mw, law)
        _, gi = D.plummer(2.0 * d['L_V'], d['r_h_pc'] / 1000.0)
        g_i = gi(d['r_h_pc'] / 1000.0); k_d = 3 * d['sigma_obs'] ** 2 / law['u_kms'] ** 2
        tot_ = np.hypot(env['gN'] * mu - g_i, env['gN'] * np.sqrt(1 - mu ** 2)) + env['S'] + k_d * g_i
        rows.append(dict(name=d['name'], over_threshold=float(np.sum(w * tot_) / 2 * K / thr)))
    out['dwarfs'] = dict(n=len(rows), share_above=float(np.mean([r['over_threshold'] > 1 for r in rows])), rows=rows)
    out['threshold_m_s2'] = thr
    return out


# ------------------------------------------------------------------------------------------------ 2a. X-COP
def xcop_gas(law, ctx):
    import law as L
    import run_v3 as R3
    import t_clusters as TC
    from scipy.optimize import minimize_scalar
    G = L.G
    a, lam, u0 = law['a_code'], law['lam'], law['u_kms']
    cls = TC.static_clusters(ctx)
    W = {c['name']: L.shell_weights(c['Rk'], c['s']) for c in cls}

    def M(c, u, eps):
        gN = G * c['Mb'] / c['Rk'] ** 2
        S = G * (W[c['name']] @ (L.heat_weight(np.sqrt(c['sig2_star_hse']), u) * c['dms'] + eps * c['dmg'])) / c['Rk'] ** 2
        return L.total(gN, S, a, lam) * c['Rk'] ** 2 / G

    rows = []
    for eps in EPS:
        res = np.array([np.log(c['Mh'] / M(c, u0, eps)) for c in cls])
        fit = minimize_scalar(lambda lu: np.mean(np.array([np.log(c['Mh'] / M(c, 10 ** lu, eps)) for c in cls]) ** 2),
                              bounds=(1.3, 4.0), method='bounded')
        uf = 10 ** fit.x
        rf = np.array([np.log(c['Mh'] / M(c, uf, eps)) for c in cls])
        rows.append(dict(eps=eps, tau=6 * eps, rms_u_fixed=float(np.sqrt(np.mean(res ** 2))), worst_radius_u_fixed=float(np.max(np.abs(res.mean(0)))),
                         mean_u_fixed=float(res.mean()), u_refit=float(uf), rms_refit=float(np.sqrt(np.mean(rf ** 2))),
                         worst_radius_refit=float(np.max(np.abs(rf.mean(0))))))
        log(f"X-COP eps {eps:g}: rms {rows[-1]['rms_u_fixed']:.3f} (u fixed), refit u {uf:.1f} -> rms {rows[-1]['rms_refit']:.3f}, "
            f"worst radius {rows[-1]['worst_radius_refit']:.3f}")
    return rows


# ------------------------------------------------------------------------------------------------ 2b. the Bullet Cluster
def _map(args, kw, eps, keep=None):
    """bullet_v4.kappa_map_v4 with eps G int rho_gas(now) / d^2 added to S (the gas's own kept brightness)."""
    import bullet_v3 as B
    import bullet_v4 as V
    import law as L
    current, ghost_gas, ghost_stars, pos, consts = args
    n, dx, centre = kw.get('n', 192), kw.get('dx', 15.0), kw.get('centre', (360., 50.))
    fresh_kpc, heat, memory = kw.get('fresh_kpc', 30.0), kw.get('heat', True), kw.get('memory', True)
    G = L.G
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
    rho_now = rho_gas_now + B.build_density(current, pos, x, y, z, 'st')
    rho_ghost_gas = B.build_density(ghost_gas, pos, x, y, z, 'gas') if memory else rho_gas_now
    rho_ghost_st = np.zeros_like(rho_now); krho = np.zeros_like(rho_now)
    for comp, prof in ghost_stars:
        rho_c = B.build_density({'c': comp}, pos, x, y, z, 'c')
        cx, cy = pos[comp['centre']]
        r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        rho_ghost_st += rho_c
        if heat:
            krho += L.k_from_sig2(np.interp(r3, prof[0], prof[1]) ** 2, u).astype(np.float32) * rho_c
        del r3, rho_c
    inv_r = B.Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    gN = -np.array(np.gradient(-G * inv_r(rho_now * dV), dx))
    Fv = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx))
    del inv_r
    if memory and fresh_kpc > 0:
        Kin = V.vec_conv(n, dx, rmax=fresh_kpc)
        Fv = Fv + G * dV * (V.apply_vec(Kin, rho_gas_now, n) - V.apply_vec(Kin, rho_ghost_gas, n))
        del Kin
    inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    S_st = G * inv_r2(krho * dV)
    S_gas = G * inv_r2(rho_gas_now * dV)
    del inv_r2
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fmag = np.sqrt(np.sum(Fv ** 2, axis=0)) + 1e-30
    hmag = np.sqrt(np.sum(g_hot ** 2, axis=0))
    out = []
    for e in (eps if isinstance(eps, (list, tuple)) else [eps]):
        S = S_st + e * S_gas
        extra = np.exp(-mag / gd) * np.sqrt(a * (Fmag + S))
        h = gN + extra * (Fv + g_hot) / (Fmag + hmag + 1e-30)
        divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
        Sig_eff = (-divh / (4 * np.pi * G)).sum(axis=2) * dx
        out.append(Sig_eff)
        del extra, h, divh
    if keep is not None:
        keep.update(x=x, y=y, z=z, dx=dx, level=(Fmag + S_st), S_gas_over_st=None)
    Sig_b = rho_now.sum(axis=2) * dx
    return x, y, out, Sig_b


def bullet_gas(law):
    import bullet_v4 as V
    import bullet_v3 as B
    import bullet_main_v5 as BM
    import bullet_static_v11 as BS
    import collisions_v4 as C4
    import collisions_v10 as C10
    import law as L
    f = C10.factors(0.296, 1.0, (70.0, 0.3)); fs, fl = f['size'], f['lens']
    got = {}
    orig = V.kappa_map_v4

    def capture(*a, **k):
        got['args'], got['kw'] = a, dict(k)
        return orig(*a, **k)
    V.kappa_map_v4 = capture
    try:
        r = BS.run_bullet(law, dict(f), n=192, dx=15.0)
    finally:
        V.kappa_map_v4 = orig
    scrit, pos = r['sigma_crit'], got['args'][3]
    keep = {}
    x, y, maps, _ = _map(got['args'], got['kw'], list(EPS), keep=keep)
    sep = {w: float(np.linalg.norm(pos[f'{w}_plasma'] - pos[f'{w}_bcg'])) for w in ('main', 'sub')}
    rows = []
    for eps, Se in zip(EPS, maps):
        kap = Se / scrit
        dec = B.clowe_decomposition(x, y, kap, pos, rmax=1200.0 * fs)
        pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0 * fs)
        mp = min(pk, key=lambda d: d['dist_main_bcg']); spk = min(pk, key=lambda d: d['dist_sub_bcg'])
        m250 = {w: BM.mass_within(x, y, kap, pos[f'{w}_bcg'], 250.0 * fs, scrit) for w in ('main', 'sub')}
        row = dict(eps=eps, tau=6 * eps, kappa_main=dec['main_bcg'], kappa_sub=dec['sub_bcg'], gas_main=dec['main_plasma'],
                   gas_sub=dec['sub_plasma'], peak_main_kpc=mp['dist_main_bcg'], peak_sub_kpc=spk['dist_sub_bcg'],
                   peak_limit_main=0.25 * sep['main'], peak_limit_sub=0.25 * sep['sub'], m250_main=m250['main'], m250_sub=m250['sub'])
        z_gm, z_gs = (row['gas_main'] - 0.05) / 0.06, (row['gas_sub'] - 0.02) / 0.06
        row['gas_checks_pass'] = bool(abs(z_gm) <= 2 and abs(z_gs) <= 2)
        row['peaks_pass'] = bool(row['peak_main_kpc'] <= row['peak_limit_main'] and row['peak_sub_kpc'] <= row['peak_limit_sub'])
        row['m250_main_in_range'] = bool(2.5e14 * fl - 2 * 0.15e14 * fl <= m250['main'] <= 2.8e14 * fl + 2 * 0.15e14 * fl)
        rows.append(row)
        log(f"Bullet eps {eps:g}: gas lensing {row['gas_main']:.3f}/{row['gas_sub']:.3f}, peaks {row['peak_main_kpc']:.0f}/{row['peak_sub_kpc']:.0f} kpc "
            f"(limits {row['peak_limit_main']:.0f}/{row['peak_limit_sub']:.0f}), m250 {m250['main']:.3g}/{m250['sub']:.3g}")
    # the threshold on the Bullet's aperture walls (250 kpc about each galaxy centre)
    thr = law['a_SI'] / 4 / L.KMS2_PER_KPC
    from scipy.ndimage import map_coordinates
    cover = {}
    for w in ('main', 'sub'):
        ph = np.linspace(0, 2 * np.pi, 180, endpoint=False)
        cx, cy = pos[f'{w}_bcg']
        X = (cx + 250.0 * fs * np.cos(ph))[:, None] + 0 * keep['z'][None, :]
        Y = (cy + 250.0 * fs * np.sin(ph))[:, None] + 0 * keep['z'][None, :]
        Z = 0 * X + keep['z'][None, :]
        idx = np.array([(X.ravel() - keep['x'][0]) / keep['dx'], (Y.ravel() - keep['y'][0]) / keep['dx'], (Z.ravel() - keep['z'][0]) / keep['dx']])
        lv = map_coordinates(keep['level'], idx, order=1, mode='nearest').reshape(X.shape)
        inner = np.abs(Z) <= 250.0 * fs
        cover[w] = float(np.mean(lv[inner] > thr))
    return rows, cover


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    log('threshold coverage')
    thr = threshold_coverage(law, ctx)
    log('X-COP with the gas keeping its own brightness')
    xc = xcop_gas(law, ctx)
    log('the Bullet Cluster with the gas keeping its own brightness')
    bl, cover = bullet_gas(law)
    thr['bullet_aperture_walls'] = cover
    out = dict(source='code/drt1_bearing_v25.py', drt1='drt1-transport/README.md (main at ab1c027)', threshold=thr, xcop=xc, bullet=bl)
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(thr, indent=1, default=float)[:3000])
    log(f'wrote {args.output}')


if __name__ == '__main__':
    main()
