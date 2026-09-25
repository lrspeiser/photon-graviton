"""Round 19, step 3: the Milky Way, Cassini and wide binaries from one Galaxy fitted to independent constraints.

The review: the Milky Way test used McMillan's (2017) matter, which was itself fitted together with a dark halo; the
law's 209 km/s at 20 kpc against 229-234 measured (20-25% too little pull) may be the matter model, not the law. So
here the Galaxy's ordinary matter is varied within independent (non-rotation-curve) constraints, and fitted to the
radial and vertical pull together, with the law's constants held at their adopted (SPARC + X-COP) values:

  rotation curve  Feng, Huang, Zhang & Liu 2026 (MNRAS 546, stag011; arXiv:2512.21780): 903 classical Cepheids,
                  Gaia DR3, 6.6-17.6 kpc, R0 = 8.275 kpc, V_phi,sun = 250.2 km/s; statistical errors (1-2 km/s) plus
                  a 3 km/s floor for the features no axisymmetric model has (their dip and bump) and a common
                  scale (1 +- 0.02: the Sun's own speed); the other Gaia curves (Eilers 2019, Zhou 2023, Ou 2024) as
                  alternates, each in its own convention
  vertical pull   K_z(1.1 kpc)/2 pi G = 70 +- 5 Msun/pc^2 at the Sun (Bland-Hawthorn & Gerhard 2016's consensus of
                  Kuijken & Gilmore 1991, 71 +- 6; Holmberg & Flynn 2004, 74 +- 6; Bovy & Rix 2013, 67-68)
  priors          stars at the Sun 33.4 +- 3 Msun/pc^2 (McKee, Parravano & Hollenbach 2015); thick disk's share of
                  them 0.12 +- 0.04 (Bland-Hawthorn & Gerhard 2016); gas 13.7 +- 1.6 (McKee et al.); bulge
                  (1.2 +- 0.35) x 10^10 Msun (Licquia & Newman 2015: 0.91; Bland-Hawthorn & Gerhard: 1.4-1.7);
                  thin disk scale length 2.6 +- 0.5 kpc, thick 2.0 +- 0.6 kpc (Bland-Hawthorn & Gerhard; wide
                  for the thick disk, whose published lengths run 1.8-4.9 kpc)

Shapes (scale lengths) on a grid, normalisations fitted continuously. The same fit with Newton's pull alone (the
same matter) for comparison. From every fit, the Galaxy's pull and heat at the Sun, which set Cassini's Q2 and the
wide-binary boost for each release length L (regression/t_precision.sun_in_galaxy).

    python code/mw_joint_v19.py --output run-mw-joint-v19/mw_joint_v19.json [--processes 4]
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, itertools, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import mw_model as M                     # noqa: E402
import milky_way_v7 as MW                # noqa: E402

PC2, G, K_SI = M.PC2, M.G, M.K_SI
FENG = dict(ref='Feng, Huang, Zhang & Liu 2026, MNRAS 546, stag011 (arXiv:2512.21780), Table 1', R0=8.275, V_phi_sun=250.2,
            R=[6.58, 7.49, 8.48, 9.50, 10.52, 11.42, 12.50, 13.50, 14.46, 15.21, 15.88, 17.58],
            v=[242.92, 239.31, 236.54, 232.36, 229.97, 233.57, 233.38, 237.05, 237.11, 237.00, 231.03, 224.03],
            err=[1.25, 1.33, 0.86, 1.32, 1.13, 1.39, 1.30, 1.53, 1.31, 2.22, 1.95, 1.91])
R0_OF = dict(Feng2026=8.275, Eilers2019=8.122, Zhou2023=8.122, Ou2024=8.178)
FLOOR = 3.0                              # km/s: non-axisymmetric features
SCALE_PRIOR = 0.02                       # the curve's common scale (the Sun's speed, +-~5 km/s)
PRIORS = dict(Sigma_star=(33.4, 3.0), thick_share=(0.12, 0.04), Sigma_gas=(13.7, 1.6), M_bulge=(1.2e10, 0.35e10),
              L_thin=(2.6, 0.5), L_thick=(2.0, 0.6), Kz11=(70.0, 5.0))
L_THIN = (2.1, 2.6, 3.1)
L_THICK = (2.0, 3.0)
Z_THIN, Z_THICK = 0.300, 0.900

_STATE = {}


def setup():
    if 'grid' in _STATE:
        return _STATE
    import common as C
    ctx = C.Context(tier='quick', verbose=False)
    comps = MW.components(); grid = M.Grid()
    cache = ctx.cache
    F = {}
    def field(key, comp):
        p = cache / f'mw_joint_v19_{key}.npz'
        if p.exists():
            d = np.load(p); return (d['R'], d['z'])
        gR, gz = M._one_field(comp, grid, 60.0, 10, 40.0, 0.005)
        np.savez(p, R=gR, z=gz); return (gR, gz)
    for k in ('bulge', 'HI', 'H2', 'halo'):
        F[k] = field(k, comps[k])
    disks = {}
    for L in L_THIN:
        c = M.Disk('thin', 1.0 * PC2 * np.exp(8.2 / L), L, Z_THIN)            # 1 Msun/pc^2 at 8.2 kpc
        disks[('thin', L)] = c; F[('thin', L)] = field(f'thin_L{L}', c)
    for L in L_THICK:
        c = M.Disk('thick', 1.0 * PC2 * np.exp(8.2 / L), L, Z_THICK)
        disks[('thick', L)] = c; F[('thick', L)] = field(f'thick_L{L}', c)
    # the heat of the bulge (at its McMillan mass) and of the stellar halo, once: S and g_hot are sums of k dm, and the
    # bulge's k grows with its mass (the SPARC rule sigma = 0.65 max v_bulge), so its part scales as the factor squared
    from law_config import load_law
    u = load_law('round12')['u_kms']
    vb = np.sqrt(np.maximum(-F['bulge'][0][:, grid.nzh] * grid.R, 0))
    comps['bulge'].sigma = 0.65 * vb[grid.R < 10].max()
    Hb = M.heat_fields([MW.Scaled(comps['bulge'], 1.0)], grid, u)
    Hh = M.heat_fields([MW.Scaled(comps['halo'], 1.0)], grid, u)
    # ring sums need target points on cell edges (a target inside a cell meets a near-singular ring)
    Rt = np.unique([grid.Re[np.argmin(abs(grid.Re - x))] for x in np.arange(5.0, 26.01, 0.25)])
    _STATE.update(comps=comps, grid=grid, F=F, disks=disks, Hb=Hb, Hh=Hh, Rt=Rt, sigma_bulge_1=float(comps['bulge'].sigma))
    return _STATE


def evaluate(p, shapes, law, laws='ours', R_eval=None, R0=8.275):
    """Circular speed at R_eval, K_z(1.1)/2piG at R0, and the Sun's pull and heat, for parameters p = (Sigma_star,
    thick_share, gas factor, bulge factor) and shapes (L_thin, L_thick)."""
    st = setup(); comps, grid, F, disks = st['comps'], st['grid'], st['F'], st['disks']
    Sst, share, fg, fb = p
    Lt, Lk = shapes
    thin, thick = disks[('thin', Lt)], disks[('thick', Lk)]
    # 1 Msun/pc^2 at 8.2 kpc in the library; the priors are at R0: rescale by the exponential
    s_thin = (1 - share) * Sst * np.exp((R0 - 8.2) / Lt); s_thick = share * Sst * np.exp((R0 - 8.2) / Lk)
    parts = {('thin', Lt): s_thin, ('thick', Lk): s_thick, 'HI': fg, 'H2': fg, 'bulge': fb, 'halo': 1.0}
    cmap = {('thin', Lt): thin, ('thick', Lk): thick, 'HI': comps['HI'], 'H2': comps['H2'], 'bulge': comps['bulge'], 'halo': comps['halo']}
    gR = sum(f * F[k][0] for k, f in parts.items()); gz = sum(f * F[k][1] for k, f in parts.items())
    rho_b = sum(f * cmap[k].rho(grid.RR, grid.ZZ) for k, f in parts.items())
    jz = grid.nzh
    consts = dict(a_code=law['a_code'], lam=law['lam'], u_kms=law['u_kms'])
    if laws == 'ours':
        Hb, Hh = st['Hb'], st['Hh']
        S, hR, hz = (fb ** 2 * Hb[i] + Hh[i] for i in range(3))
        eR, ez = M.law_extra(gR, gz, S, hR, hz, consts, law='ours', reach=law['reach_kpc'], grid=grid)
        rho = rho_b - M.divergence(grid, eR, ez) / (4 * np.pi * G)
    else:
        S = np.zeros_like(gR); hR = hz = np.zeros_like(gR)
        rho = rho_b
    Rt = st['Rt']; Re = np.asarray(R_eval, float)
    lo, hi = np.searchsorted(Rt, Re.min()) - 2, np.searchsorted(Rt, Re.max()) + 2
    Ru = Rt[max(lo, 0):hi]
    g1, _ = M.ring_sum(Ru, 0 * Ru, grid, rho)
    v = np.interp(Re, Ru, np.sqrt(np.maximum(-g1 * Ru, 0)))
    R0e = grid.Re[np.argmin(abs(grid.Re - R0))]; z11 = grid.ze[np.argmin(abs(grid.ze - 1.1))]
    _, gzz = M.ring_sum(np.array([R0e]), np.array([z11]), grid, rho)
    Kz = float(-gzz[0] / (2 * np.pi * G) / PC2)
    iR = np.argmin(abs(grid.R - R0))
    sun = dict(g_N_SI=float(np.hypot(gR[iR, jz], gz[iR, jz]) * K_SI), S_SI=float(S[iR, jz] * K_SI),
               g_hot_SI=float(np.hypot(hR[iR, jz], hz[iR, jz]) * K_SI))
    M_b = fb * comps['bulge'].mass(); Sgas = fg * (comps['HI'].Sigma(R0) + comps['H2'].Sigma(R0)) / PC2
    masses = dict(thin=float(s_thin * thin.mass()), thick=float(s_thick * thick.mass()), bulge=float(M_b),
                  gas=float(fg * (comps['HI'].mass() + comps['H2'].mass())))
    return v, Kz, sun, dict(M_bulge=float(M_b), Sigma_gas=float(Sgas), masses=masses)


def chi2(p, shapes, law, laws, data, return_parts=False):
    Sst, share, fg, fb = p
    if Sst <= 0 or not (0 <= share <= 0.9) or fg <= 0 or fb <= 0:
        return 1e9
    R, vd, e = np.array(data['R']), np.array(data['v']), np.array(data['err'])
    v, Kz, sun, extra = evaluate(p, shapes, law, laws, R, data.get('R0', 8.275))
    et2 = e ** 2 + FLOOR ** 2
    # the common scale s of the data (the Sun's speed), marginalised in closed form: data -> s data
    A = np.sum(vd * vd / et2) + 1 / SCALE_PRIOR ** 2; B = np.sum(v * vd / et2) + 1 / SCALE_PRIOR ** 2
    s = B / A
    c_rc = float(np.sum((v - s * vd) ** 2 / et2) + ((s - 1) / SCALE_PRIOR) ** 2)
    pri = {'Sigma_star': Sst, 'thick_share': share, 'Sigma_gas': extra['Sigma_gas'], 'M_bulge': extra['M_bulge'],
           'L_thin': shapes[0], 'L_thick': shapes[1], 'Kz11': Kz}
    c_pri = {k: float(((pri[k] - PRIORS[k][0]) / PRIORS[k][1]) ** 2) for k in PRIORS}
    tot = c_rc + sum(c_pri.values())
    if return_parts:
        return dict(total=tot, rotation_curve=c_rc, scale=float(s), priors=c_pri, values=pri, v=v.tolist(), R=R.tolist(),
                    residual=(v - s * vd).tolist(), Kz11=Kz, sun=sun, extra=extra, n_points=len(R))
    return tot


def fit_one(args):
    shapes, laws, dname, data = args
    from law_config import load_law
    law = load_law('round12')
    x0 = np.array([33.4, 0.12, 13.7 / 12.23, 1.2e10 / 8.878e9])
    f = lambda x: chi2(x, shapes, law, laws, data)
    r = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-3, fatol=1e-3, maxiter=400, initial_simplex=np.array(
        [x0, x0 * [1.15, 1, 1, 1], x0 * [1, 1.6, 1, 1], x0 * [1, 1, 1.2, 1], x0 * [1, 1, 1, 1.3]])))
    parts = chi2(r.x, shapes, law, laws, data, return_parts=True)
    return dict(shapes=dict(L_thin=shapes[0], L_thick=shapes[1]), law=laws, data=dname, params=dict(Sigma_star=float(r.x[0]),
                thick_share=float(r.x[1]), gas_factor=float(r.x[2]), bulge_factor=float(r.x[3])), nfev=int(r.nfev), **parts)


def curves():
    rc = json.loads((HERE.parent / 'data/mw_rotation_curves.json').read_text())
    out = {'Feng2026': dict(R=FENG['R'], v=FENG['v'], err=FENG['err'], R0=FENG['R0'])}
    for k in ('Eilers2019', 'Zhou2023', 'Ou2024'):
        d = rc[k]; R = np.array(d['R']); m = (R >= 5.0) & (R <= 25.0)
        out[k] = dict(R=R[m].tolist(), v=np.array(d['v'])[m].tolist(), err=np.array(d['err'])[m].tolist(), R0=R0_OF[k])
    return out


def cassini_binaries(sun, law, Ls_pc=(0.0, 0.05, 0.10, 0.15, 0.19, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0)):
    from t_precision import sun_in_galaxy, MSUN
    PC_AU = 206264.806
    rows = []
    for L in Ls_pc:
        Q2, boost = sun_in_galaxy(MSUN, sun['g_N_SI'], sun['S_SI'], sun['g_hot_SI'], L * PC_AU, law['a_SI'], law['g_d_SI'])
        rows.append(dict(L_pc=L, Q2=Q2, z_2026=(Q2 - 1.6e-27) / 1.8e-27, boost_7000=boost[7000], boost_20000=boost[20000]))
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    setup()                                                           # builds the field library once (cached)
    data = curves()
    jobs = [(sh, laws, 'Feng2026', data['Feng2026']) for sh in itertools.product(L_THIN, L_THICK) for laws in ('ours', 'newton')]
    jobs += [(sh, 'ours', k, data[k]) for sh in itertools.product(L_THIN, L_THICK) for k in ('Eilers2019', 'Zhou2023', 'Ou2024')]
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(fit_one, jobs):
            res.append(r)
            print(f"[{time.monotonic() - t0:5.0f} s] {r['data']:10s} {r['law']:6s} L_thin {r['shapes']['L_thin']} L_thick {r['shapes']['L_thick']}: chi2 {r['total']:7.2f} "
                  f"(curve {r['rotation_curve']:6.2f} on {r['n_points']}, scale {r['scale']:.3f}; K_z {r['Kz11']:.1f}) Sigma* {r['params']['Sigma_star']:.1f}, "
                  f"thick {r['params']['thick_share']:.2f}, gas x{r['params']['gas_factor']:.2f}, bulge x{r['params']['bulge_factor']:.2f}; "
                  f"priors " + ', '.join(f"{k} {v:.1f}" for k, v in r['priors'].items()), flush=True)
            args.output.write_text(json.dumps(dict(experiment='round 19: the Milky Way fitted to independent constraints', fits=res,
                                                   seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    from law_config import load_law
    law = load_law('round12')
    best = {}
    for laws in ('ours', 'newton'):
        rs = [r for r in res if r['data'] == 'Feng2026' and r['law'] == laws]
        b = min(rs, key=lambda r: r['total']); best[laws] = b
    # Cassini and wide binaries from the best fit and from every shape's fit (the spread they give)
    cas = dict(best=cassini_binaries(best['ours']['sun'], law),
               spread=[dict(shapes=r['shapes'], total=r['total'], sun=r['sun'], rows=cassini_binaries(r['sun'], law))
                       for r in res if r['data'] == 'Feng2026' and r['law'] == 'ours'])
    out = dict(experiment='round 19: the Milky Way fitted to independent constraints, with Cassini and wide binaries', data=FENG, priors=PRIORS,
               floor_kms=FLOOR, scale_prior=SCALE_PRIOR, fits=res, best=best, cassini_binaries=cas, seconds=time.monotonic() - t0)
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    for laws, b in best.items():
        print(f"best {laws}: chi2 {b['total']:.2f} (curve {b['rotation_curve']:.2f} on 12), shapes {b['shapes']}, params {b['params']}, K_z {b['Kz11']:.1f}, sun {b['sun']}")
    for row in cas['best']:
        print(f"  L = {row['L_pc']:.2f} pc: Q2 = {row['Q2']:.2e} ({row['z_2026']:+.2f} sigma), binaries +{100 * (row['boost_7000'] - 1):.1f}% at 7000 AU, +{100 * (row['boost_20000'] - 1):.1f}% at 20000 AU")


if __name__ == '__main__':
    main()
