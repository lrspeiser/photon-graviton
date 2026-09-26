"""Round 25: flow frustration. Do the adopted law's misses grow with the share of the companion's energy that is not
flowing? The test registered in round25-flow-frustration.md, run as registered.

The proposal's measure, with the law's own fields at each test point:

    chi_law = 1 - |g_N + g_hot| / (|g_N| + S)

(|g_N| + S: the companion's total intensity in the law; |g_N + g_hot|: the part that flows as one net stream; on the
collision grids the flowing companion is the memory field F, so chi_law = 1 - |F + g_hot| / (|F| + S)). Residual:
ln(observed / law) of the pull or of a quantity proportional to it. Secondary: chi_all = 1 - |g_N + g_hot| / (S_cold + S),
the ordered matter's energy also counted as brightness (S_cold = G int rho / d^2; spherical stand-ins, exact on the
collision grids), and the collisions' galaxy speeds with chi of today's collision field.

    python code/frustration_v25.py --output-dir run-frustration-v25            # everything (collisions: ~15 min)
    python code/frustration_v25.py --output-dir run-frustration-v25 --skip-collisions
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

BASE = {c['id']: c for c in json.loads((RESULTS / 'regression/baseline.json').read_text())['checks']}
T0 = time.monotonic()


def log(msg):
    print(f'[{time.monotonic() - T0:5.0f} s] {msg}', flush=True)


def status(cid):
    return BASE[cid]['status'] if cid in BASE else 'info'


def row(family, cid, res, err, chi, chi_all, gN, **kw):
    """One test point. gN in m/s^2 (the Newtonian pull there, for the control)."""
    return dict(family=family, id=cid, status=kw.pop('status', status(cid)), res=float(res), err=float(err), chi=float(chi),
                chi_all=float(chi_all) if chi_all is not None else None, gN_SI=float(gN), **kw)


def _shells_mid(s, M):
    """Shells between the nodes of a mass profile M(s): midpoints (geometric) and masses (never on a node)."""
    return np.sqrt(s[1:] * s[:-1]), np.diff(M)


# ------------------------------------------------------------------------------------------------ SPARC
def sparc(law, ctx):
    import run as R
    import law as L
    a, lam, u, G, K = law['a_code'], law['lam'], law['u_kms'], L.G, L.KMS2_PER_KPC
    out = []
    for g in ctx.sparc():
        r, gN = g['r'], g['gN']
        gl = R.galaxy_g(g, a, u, lam)
        s, dmb = g['sfine'], g['dmb']
        if g['sigb'] > 0:                                   # the law's hot matter: the bulge (spherical shells)
            k = L.heat_weight(g['sigb'], u)
            S = L.scalar_sum(r, s, dmb, k)
            ghot = G * k * np.array([dmb[s < x].sum() for x in r]) / r ** 2
        else:
            S = ghot = np.zeros_like(r)
        chi = 1 - (gN + ghot) / (gN + S)
        # secondary: the ordered matter's brightness, from spherical shells of the rotation curve's baryonic mass
        Mr = np.maximum(g['vbar2'], 0.0) * r / G
        sg = np.geomspace(r.min() / 20, r.max(), 401)
        Ms = np.maximum.accumulate(np.where(sg < r[0], Mr[0] * (sg / r[0]) ** 2, np.interp(sg, r, Mr)))
        sm, dm = _shells_mid(sg, Ms)
        Scold = G * (L.shell_weights(r, sm) @ dm) / r ** 2
        chi_all = 1 - (gN + ghot) / (Scold + S)
        res = np.log(g['v'] ** 2 / r / gl)
        err = 2 * g['err'] / g['v']
        bulgy = bool(np.max(g['vb2'] / np.maximum(g['vbar2'], 1e-9)) > 0.5)
        for i in range(r.size):
            out.append(row('sparc', 'galaxies', res[i], err[i], chi[i], chi_all[i], gN[i] * K, status='pass',
                           galaxy=g['name'], r_kpc=float(r[i]), bulge=bool(g['sigb'] > 0), bulge_dominated=bulgy,
                           S_over_gN=float(S[i] / gN[i]), ghot_over_gN=float(ghot[i] / gN[i])))
    return out


# ------------------------------------------------------------------------------------------------ X-COP
def xcop(law, ctx):
    import law as L
    import run_v3 as R3
    import t_clusters as TC
    a, lam, u, G, K = law['a_code'], law['lam'], law['u_kms'], L.G, L.KMS2_PER_KPC
    out = []
    for c in TC.static_clusters(ctx):
        Rk, s = c['Rk'], c['s']
        W = L.shell_weights(Rk, s)
        gN = G * c['Mb'] / Rk ** 2
        src = L.heat_weight(np.sqrt(c['sig2_star_hse']), u) * c['dms']
        S = G * (W @ src) / Rk ** 2
        ghot = G * np.array([src[s < x].sum() for x in Rk]) / Rk ** 2
        Scold = G * (W @ (c['dmg'] + c['dms'])) / Rk ** 2
        M = R3.cluster_M3(c, a, u, lam)
        res, err = np.log(c['Mh'] / M), c['eMh'] / c['Mh']
        for i in range(Rk.size):
            out.append(row('xcop', 'clusters', res[i], err[i], 1 - (gN[i] + ghot[i]) / (gN[i] + S[i]),
                           1 - (gN[i] + ghot[i]) / (Scold[i] + S[i]), gN[i] * K, status='pass', cluster=c['name'],
                           r_kpc=float(Rk[i]), x_R500=float(Rk[i] / c['R5']), S_over_gN=float(S[i] / gN[i]),
                           ghot_over_gN=float(ghot[i] / gN[i]), w_max=float(W[i].max())))
    return out


# ------------------------------------------------------------------------------------------------ KiDS and Mistele
def lensing(law):
    """The suite's point lenses: S = k g_N and g_hot = k g_N along g_N at every radius, so chi = 0 exactly (and chi_all:
    a point's brightness equals its pull). The residuals per bin, from the suite's own comparison."""
    import kids_static_v11 as KS
    import kids_heat_v12 as KH
    import lensing_census_v7 as LC
    import collisions_v10 as C10
    import common as C
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
    tabs = KS.tables(); k = KH.heat(law['u_kms'])
    fg, fl, fm = f['stars'] / f['size'] ** 2, f['lens'] / f['size'] ** 2, f['stars']
    gb = tabs['all']['gbar'] * fg
    Mtyp = 10 ** 10.6 * fm
    g = {s: LC.gconv_at(gb, Mtyp, 'ours', consts, k=k[s], reach=law['reach_kpc']) for s in KH.PAIRS}
    fr = tabs['red']['w'] / (tabs['red']['w'] + tabs['blue']['w'])
    mix = fr * g['red'] + (1 - fr) * g['blue']
    rel = tabs['all']['gbar'] >= 1e-13
    out = []
    for s, pred in (('all', mix), ('blue', g['blue']), ('red', g['red']), ('disc', g['disc']), ('bulge', g['bulge']), ('gama', mix)):
        obs, err = tabs[s]['gobs'] * fl, tabs[s]['err'] * fl
        for i in np.where(rel)[0]:
            out.append(row('kids', f'lensing.kids_{s}', np.log(obs[i] / pred[i]), err[i] / obs[i], 0.0, 0.0, gb[i], sample=s))
    # Mistele et al. 2024: circular speeds 50-300 kpc (speed squared is proportional to the pull)
    m = KH.mistele(consts, law['u_kms'], law['reach_kpc'], f)
    for s in ('LTG', 'ETG'):
        for q in m[s]['observed_over_predicted']:
            out.append(row('mistele', f'lensing.mistele_{s.lower()}', 2 * np.log(q), np.nan, 0.0, 0.0, np.nan, sample=s))
    return out


# ------------------------------------------------------------------------------------------------ Milky Way
def milky_way(law, ctx):
    import milky_way_v7 as MW
    import mw_model as M
    import law as L
    import t_milky_way as TMW
    import t_dwarfs as TD
    from scipy.interpolate import RegularGridInterpolator as RGI
    G, K = L.G, L.KMS2_PER_KPC
    comps, grid, F = TMW.setup(ctx)
    ev = MW.Evaluator(comps, grid, F, law)
    parts = MW.models(comps)['M17']
    res_law = ev.run(parts, 'ours', reach=law['reach_kpc'])          # also sets the bulge's heat (the SPARC rule)
    Rt, vt = np.array(res_law['R']), np.array(res_law['v'])
    gR = sum(f * F[k][0] for k, f in parts.items()); gz = sum(f * F[k][1] for k, f in parts.items())
    hot = [MW.Scaled(comps[k], f) for k, f in parts.items() if k in ('bulge', 'halo')]
    S, hR, hz = M.heat_fields(hot, grid, law['u_kms'])
    eR, ez = M.law_extra(gR, gz, S, hR, hz, law, law='ours', reach=law['reach_kpc'], grid=grid)
    gmag = np.hypot(gR, gz)
    chi_f = 1 - np.hypot(gR + hR, gz + hz) / (gmag + S + 1e-300)
    I = lambda A: RGI((grid.R, grid.z), A, bounds_error=False, fill_value=None)
    fchi, fg = I(chi_f), I(gmag)
    # secondary: spherical stand-ins for the whole Galaxy (the dwarfs' profile): cold and hot brightness, pulls
    prof = TD.galaxy_profile(law, ctx)
    s, dm, dmk = prof['s'], prof['dm'], prof['dmk']
    sm, dms = _shells_mid(s, np.cumsum(dm)); _, dmks = _shells_mid(s, np.cumsum(dmk))

    def chi_all_sph(r):
        r = np.atleast_1d(r)
        W = L.shell_weights(r, sm)
        gNs = G * np.array([dms[sm < x].sum() for x in r]) / r ** 2
        ghs = G * np.array([dmks[sm < x].sum() for x in r]) / r ** 2
        return 1 - (gNs + ghs) / (G * (W @ dms) / r ** 2 + G * (W @ dmks) / r ** 2)

    out = []
    plane = lambda R, z=0.0: np.c_[np.atleast_1d(R), np.full(np.size(R), z)]
    # rotation curves, 5-27 kpc
    rc = json.loads((RESULTS / 'data/mw_rotation_curves.json').read_text())
    for name, d in rc.items():
        R, v, e = np.array(d['R']), np.array(d['v']), np.array(d['err'])
        m = (R >= 5.0) & (R <= 27.5)
        vl = np.interp(R[m], Rt, vt); et = np.hypot(e[m], d['sys_frac'] * v[m])
        ch, gg, ca = fchi(plane(R[m])), fg(plane(R[m])), chi_all_sph(R[m])
        for i in range(m.sum()):
            out.append(row('mw', 'mw.rc_outer' if R[m][i] >= 15 else 'mw.rc_rms_eilers', 2 * np.log(v[m][i] / vl[i]), 2 * et[i] / v[m][i],
                           ch[i], ca[i], gg[i] * K, analysis=name, r_kpc=float(R[m][i])))
    # the other graded radii
    v0 = float(np.interp(MW.R0, Rt, vt))
    out.append(row('mw', 'mw.v_sun', 2 * np.log(231.5 / v0), 2 * 7.0 / 231.5, fchi(plane(MW.R0))[0], chi_all_sph(MW.R0)[0],
                   fg(plane(MW.R0))[0] * K, r_kpc=MW.R0))
    share = BASE['mw.inner_share']['value']
    Rin = np.linspace(2.0, 3.0, 11)
    out.append(row('mw', 'mw.inner_share', 2 * np.log(share / 0.88), 2 * 0.07 / 0.88, float(np.mean(fchi(plane(Rin)))),
                   float(np.mean(chi_all_sph(Rin))), float(np.mean(fg(plane(Rin)))) * K, r_kpc=2.5))
    out.append(row('mw', 'mw.vertical_pull', np.log(69.8 / BASE['mw.vertical_pull']['value']), 3.3 / 69.8,
                   fchi(plane(MW.R0, 1.1))[0], chi_all_sph(np.hypot(MW.R0, 1.1))[0], fg(plane(MW.R0, 1.1))[0] * K, r_kpc=MW.R0, z_kpc=1.1))
    # enclosed masses: chi on the sphere, weighted by the companion's inward flux through it
    fe = (I(eR), I(ez))
    mu, w = np.polynomial.legendre.leggauss(400)
    for rr, obs, err in ((20.0, 1.91e11, 0.18e11), (50.0, 4.5e11, 0.4e11), (100.0, 6.9e11, 0.7e11), (200.0, 1.10e12, 0.25e12)):
        st = np.sqrt(1 - mu ** 2); pts = np.c_[rr * st, rr * mu]
        inward = np.maximum(-(fe[0](pts) * st + fe[1](pts) * mu), 0.0) * w
        cid = f'mw.mass_{int(rr)}'
        out.append(row('mw', cid, np.log(obs / BASE[cid]['value']), err / obs, float(np.sum(inward * fchi(pts)) / inward.sum()),
                       chi_all_sph(rr)[0], float(np.sum(inward * fg(pts)) / inward.sum()) * K, r_kpc=rr))
    return out


# ------------------------------------------------------------------------------------------------ dwarfs
def dwarfs(law, ctx):
    import mw_dwarfs_v7 as D
    import law as L
    import t_dwarfs as TD
    G, K, u = L.G, L.KMS2_PER_KPC, law['u_kms']
    data = json.loads((RESULTS / 'data/mw_dwarfs.json').read_text())
    mw = TD.galaxy_profile(law, ctx)
    s = mw['s']; sm, dmm = _shells_mid(s, np.cumsum(mw['dm']))
    hold = law['external_hold']
    mu, w = np.polynomial.legendre.leggauss(256)
    out = []
    for d in data['dwarfs']:
        cid = 'dwarfs.' + d['name'].lower().replace(' ', '_')
        env = D.galaxy_env(d['D_gc_kpc'], mw, law)
        ge, Se = hold * env['gN'], hold * env['S']
        Mdw, b = 2.0 * d['L_V'], d['r_h_pc'] / 1000.0
        rho, gi = D.plummer(Mdw, b)
        k_d = 3 * d['sigma_obs'] ** 2 / u ** 2
        r = np.geomspace(0.02 * b, 30 * b, 600)
        g_i = gi(r)[:, None]
        gN_par, gN_perp = ge * mu[None, :] - g_i, ge * np.sqrt(1 - mu ** 2)[None, :] + 0 * g_i
        gh_par, gh_perp = Se * mu[None, :] - k_d * g_i, Se * np.sqrt(1 - mu ** 2)[None, :] + 0 * g_i
        gmag = np.hypot(gN_par, gN_perp)
        flow = np.hypot(gN_par + gh_par, gN_perp + gh_perp)
        chi = 1 - flow / (gmag + Se + k_d * g_i)
        # secondary: the Galaxy's cold brightness at the dwarf (spherical shells, both ways) and the dwarf's own
        S_mw_cold = G * float(L.shell_weights(np.array([d['D_gc_kpc']]), sm)[0] @ dmm) / d['D_gc_kpc'] ** 2
        sd = np.geomspace(1e-3 * b, 300 * b, 3001)
        Md = Mdw * sd ** 3 / (sd ** 2 + b ** 2) ** 1.5
        sdm, dmd = _shells_mid(sd, Md)
        S_dw_cold = G * (L.shell_weights(r, sdm) @ dmd) / r ** 2
        chi_all = 1 - flow / (hold * S_mw_cold + S_dw_cold[:, None] + Se + k_d * g_i)
        wr = rho(r) * r ** 3                                  # Plummer mass per ln r
        av = lambda X: float(np.sum(wr[:, None] * X * w[None, :] / 2) / np.sum(wr))
        sig = BASE[cid]['value']
        i_b = np.argmin(np.abs(r - b))
        out.append(row('dwarfs', cid, 2 * np.log(d['sigma_obs'] / sig), 2 * d['sigma_err'] / d['sigma_obs'], av(chi), av(chi_all),
                       float(np.hypot(ge, gi(b))) * K, D_kpc=d['D_gc_kpc'], ext_over_int_at_rh=float(ge / gi(b)),
                       chi_at_rh=float(np.sum(chi[i_b] * w) / 2), sigma_law=sig, sigma_obs=d['sigma_obs']))
    return out


# ------------------------------------------------------------------------------------------------ collisions
FIELDS = []          # one entry per collision map computed while the capture is installed


def _map_with_fields(current, ghost_gas, ghost_stars, pos, consts, n=192, dx=15.0, centre=(360., 50.),
                     fresh_kpc=30.0, heat=True, memory=True):
    """bullet_v4.kappa_map_v4, line for line, that also keeps its fields (and the ordered matter's brightness)."""
    import bullet_v3 as B
    import bullet_v4 as V
    import law as L
    G = L.G
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_now = B.build_density(current, pos, x, y, z, 'gas') + B.build_density(current, pos, x, y, z, 'st')
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
    rho_st_now = B.build_density(current, pos, x, y, z, 'st')
    rho_ghost_gas = B.build_density(ghost_gas, pos, x, y, z, 'gas') if memory else rho_gas_now
    rho_ghost_st = np.zeros_like(rho_now); krho = np.zeros_like(rho_now)
    for comp, prof in ghost_stars:
        cname = {'c': comp}
        rho_c = B.build_density(cname, pos, x, y, z, 'c')
        cx, cy = pos[comp['centre']]
        r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        rho_ghost_st += rho_c
        if heat:
            krho += L.k_from_sig2(np.interp(r3, prof[0], prof[1]) ** 2, u).astype(np.float32) * rho_c
        del r3, rho_c
    inv_r = B.Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    gN = -np.array(np.gradient(-G * inv_r(rho_now * dV), dx))
    F = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx)) if heat else 0.0
    del inv_r
    if memory and fresh_kpc > 0:
        Kin = V.vec_conv(n, dx, rmax=fresh_kpc)
        F = F + G * dV * (V.apply_vec(Kin, rho_gas_now, n) - V.apply_vec(Kin, rho_ghost_gas, n))
        del Kin
    o2 = B.cube_average(lambda r: 1 / r ** 2) / dx ** 2
    inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, o2)
    S = G * inv_r2(krho * dV) if heat else 0.0
    Scold = G * inv_r2((rho_ghost_gas + rho_ghost_st) * dV)             # secondary: the ordered matter's brightness
    del inv_r2
    if memory and fresh_kpc > 0:
        fr2 = B.Conv(n, dx, lambda r: np.where(r < fresh_kpc, 1.0 / r ** 2, 0.0), o2)
        Scold = Scold + G * dV * (fr2(rho_gas_now) - fr2(rho_ghost_gas))
        del fr2
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fmag = np.sqrt(np.sum(F ** 2, axis=0)) + 1e-30
    hmag = np.sqrt(np.sum(g_hot ** 2, axis=0)) if heat else 0.0
    extra = np.exp(-mag / gd) * np.sqrt(a * (Fmag + S))
    h = gN + extra * (F + g_hot) / (Fmag + hmag + 1e-30)
    divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
    Sig_eff = (-divh / (4 * np.pi * G)).sum(axis=2) * dx
    Sig_b = rho_now.sum(axis=2) * dx
    flow = np.sqrt(np.sum((F + g_hot) ** 2, axis=0))
    FIELDS.append(dict(x=x, y=y, z=z, dx=dx, pos={k: np.asarray(v, float) for k, v in pos.items()},
                       chi=(1 - flow / (Fmag + S)).astype(np.float32), chi_all=(1 - flow / (Scold + S)).astype(np.float32),
                       ex=(extra * (F + g_hot) / (Fmag + hmag + 1e-30)).astype(np.float32), gN=mag.astype(np.float32),
                       h=h.astype(np.float32), rho_st=rho_st_now, G=G))
    return x, y, Sig_eff, Sig_b


def _interp(fd, A, pts):
    from scipy.ndimage import map_coordinates
    x0, y0, z0, dx = fd['x'][0], fd['y'][0], fd['z'][0], fd['dx']
    idx = np.array([(pts[:, 0] - x0) / dx, (pts[:, 1] - y0) / dx, (pts[:, 2] - z0) / dx])
    return map_coordinates(A, idx, order=1, mode='nearest')


def wall(fd, centre, R, nphi=240):
    """chi on the wall of the cylinder of radius R about `centre` (all z on the grid), weighted by the companion's
    inward flux through it; also the lensing mass the wall's flux gives (a cross-check of the aperture mass)."""
    ph = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    X = centre[0] + R * np.cos(ph)[:, None] + 0 * fd['z'][None, :]
    Y = centre[1] + R * np.sin(ph)[:, None] + 0 * fd['z'][None, :]
    Z = 0 * X + fd['z'][None, :]
    pts = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    nx, ny = np.repeat(np.cos(ph), fd['z'].size), np.repeat(np.sin(ph), fd['z'].size)
    ein = -(_interp(fd, fd['ex'][0], pts) * nx + _interp(fd, fd['ex'][1], pts) * ny)
    hin = -(_interp(fd, fd['h'][0], pts) * nx + _interp(fd, fd['h'][1], pts) * ny)
    wgt = np.maximum(ein, 0.0)
    dA = R * (2 * np.pi / nphi) * fd['dx']
    return dict(chi=float(np.sum(wgt * _interp(fd, fd['chi'], pts)) / wgt.sum()),
                chi_all=float(np.sum(wgt * _interp(fd, fd['chi_all'], pts)) / wgt.sum()),
                gN_SI=float(np.sum(wgt * _interp(fd, fd['gN'], pts)) / wgt.sum()),
                M_wall=float(np.sum(hin) * dA / (4 * np.pi * fd['G'])), companion_share=float(np.sum(ein) / np.sum(hin)))


def ball(fd, centre, R):
    """chi of today's collision field in the sphere of radius R about `centre` (in the sky plane, z = 0), weighted by
    today's stars (secondary, for the galaxy-speed checks)."""
    X, Y, Z = np.meshgrid(fd['x'], fd['y'], fd['z'], indexing='ij')
    m = (X - centre[0]) ** 2 + (Y - centre[1]) ** 2 + Z ** 2 <= R ** 2
    w = fd['rho_st'][m]
    return dict(chi=float(np.sum(w * fd['chi'][m]) / w.sum()), chi_all=float(np.sum(w * fd['chi_all'][m]) / w.sum()))


def settled(sub_gas, stars, law, rmax, Rap):
    """chi in the settled pre-collision cluster that predicts the galaxies' speed spread (bullet_v4.own_sigma_multi),
    weighted by the stars inside the aperture (bullet_main_v5.sigma_los_aperture's weights)."""
    import bullet_v4 as V
    import law as L
    G, u = L.G, law['u_kms']
    r, sig, g = V.own_sigma_multi(sub_gas, stars, law, rmax=rmax)
    dr = np.gradient(r)
    dmg = sum(V.shells_on(r, dr, c) for c in sub_gas); dms = sum(V.shells_on(r, dr, c) for c in stars)
    gN = G * np.cumsum(dmg + dms) / r ** 2
    W = L.shell_weights(r, r)                               # as the prediction computes it
    src = L.k_from_sig2(sig ** 2, u) * dms
    S = G * (W @ src) / r ** 2
    ghot = G * np.cumsum(src) / r ** 2
    Scold = G * (W @ (dmg + dms)) / r ** 2
    chi, chi_all = 1 - (gN + ghot) / (gN + S), 1 - (gN + ghot) / (Scold + S)
    x = np.clip(Rap / r, 0.0, 1.0); c0 = np.where(r <= Rap, 0.0, np.sqrt(1 - x ** 2))
    w = dms * (1 - c0)
    return dict(chi=float(np.sum(w * chi) / w.sum()), chi_all=float(np.sum(w * chi_all) / w.sum()),
                gN_SI=float(np.sum(w * gN) / w.sum()) * L.KMS2_PER_KPC)


def collisions(law, ctx):
    import bullet_v4 as V
    import bullet_static_v11 as BS
    import collisions_v10 as C10
    import collisions_v8 as V8
    import t_new_collisions as NC
    out = []
    orig = V.kappa_map_v4
    V.kappa_map_v4 = _map_with_fields
    try:
        # ---- the Bullet Cluster (t_collisions.bullet)
        log('Bullet')
        f = C10.factors(0.296, 1.0, (70.0, 0.3)); fs, fl = f['size'], f['lens']
        FIELDS.clear()
        r = BS.run_bullet(law, dict(f), n=192, dx=15.0)
        fd = FIELDS[-1]
        for w, lo, hi in (('main', 2.5e14, 2.8e14), ('sub', 2.0e14, 2.3e14)):
            cid = f'bullet.m250_{w}'
            wl = wall(fd, fd['pos'][f'{w}_bcg'], 250.0 * fs)
            mid = 0.5 * (lo + hi) * fl
            out.append(row('collisions', cid, np.log(mid / r['checks'][f'm250_{w}']['value']), 0.15e14 * fl / mid if w == 'main' else 0.2e14 * fl / mid,
                           wl['chi'], wl['chi_all'], wl['gN_SI'] * 3.2408e-14, system='Bullet', aperture_kpc=250.0 * fs,
                           M_law=r['checks'][f'm250_{w}']['value'], M_wall=wl['M_wall'], companion_share=wl['companion_share']))
        # ---- MACS J0025, Abell 520, El Gordo (t_new_collisions)
        for name, t in (('macs0025', 0.3), ('abell520', 0.3), ('el_gordo', 0.46)):
            log(name)
            spec = dict(macs0025=V8.macs0025, abell520=V8.abell520, el_gordo=V8.el_gordo)[name]()
            ff = C10.factors(spec['z'], 1.0); fs = ff['size']
            st = C10.rescale(spec, ff, star_extra=C10.chabrier_basis(name))
            FIELDS.clear()
            sol = C10.solve(st, law, t, fs)
            m = C10.measure(name, st, sol, fs)
            fd = FIELDS[-1]
            pos = fd['pos']
            sets = {}
            for sb in st['subs']:                              # the settled clusters behind the galaxy speeds
                g = V8.beta(st['gas_total'] * sb['gas_share'], sb['rc_pre'], sb['rt_pre'], sb['galaxies'])
                stars = [st['current'][k] for k in sb['stars']]
                sets[sb['name']] = (g, stars, sb['galaxies'])
            if name == 'macs0025':
                L_ = C10.lens_factor('macs0025')['lens']
                for w, (o, lo, hi) in (('se', (2.5e14, 1.7e14, 1.0e14)), ('nw', (2.6e14, 1.4e14, 0.5e14))):
                    wl = wall(fd, pos[f'{w}_gal'], 300.0 * fs)
                    mv = m[f'M300_{w}']
                    e = (hi if mv > o * L_ else lo) * L_ / (o * L_)
                    out.append(row('collisions', f'macs0025.m300_{w}', np.log(o * L_ / mv), e, wl['chi'], wl['chi_all'], wl['gN_SI'] * 3.2408e-14,
                                   system='MACS J0025', aperture_kpc=300.0 * fs, M_law=mv, M_wall=wl['M_wall'], companion_share=wl['companion_share']))
                st_ = [settled([sets[k][0]], sets[k][1], law, 3000.0 * max(fs, 1.0), 1500.0 * fs) for k in sets]
                bl = [ball(fd, pos[sets[k][2]], 1500.0 * fs) for k in sets]
                out.append(row('collisions', 'macs0025.sigma', 2 * np.log(835.0 / m['sigma_los_1p5Mpc']), 2 * 59.0 / 835.0,
                               np.mean([q['chi'] for q in st_]), np.mean([q['chi_all'] for q in st_]), np.mean([q['gN_SI'] for q in st_]),
                               system='MACS J0025', sigma_law=m['sigma_los_1p5Mpc'], chi_today=np.mean([q['chi'] for q in bl]),
                               chi_all_today=np.mean([q['chi_all'] for q in bl])))
            elif name == 'abell520':
                L_ = C10.lens_factor('abell520')['lens']; L7 = C10.lens_factor('abell520', '710')['lens']
                jee = dict(P1=(2.10, 0.43), P2=(4.05, 0.28), P3=(3.35, 0.34), P4=(4.23, 0.28), P5=(2.93, 0.39))
                clowe = dict(P1=(2.81, 0.67), P2=(4.16, 0.67), P3=(2.84, 0.64), P4=(5.59, 0.68), P5=(3.17, 0.66), P6=(3.68, 0.68))
                for k in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
                    vals = [dd[k] for dd in (jee, clowe) if k in dd]
                    lo, hi = min(v[0] for v in vals) * 1e13 * L_, max(v[0] for v in vals) * 1e13 * L_
                    err = max(v[1] for v in vals) * 1e13 * L_
                    wl = wall(fd, pos[k], 150.0 * fs)
                    mv = m[f'M150_{k}']
                    out.append(row('collisions', f'abell520.m150_{k.lower()}', np.log(0.5 * (lo + hi) / mv), err / (0.5 * (lo + hi)), wl['chi'],
                                   wl['chi_all'], wl['gN_SI'] * 3.2408e-14, system='Abell 520', aperture_kpc=150.0 * fs, M_law=mv,
                                   M_wall=wl['M_wall'], companion_share=wl['companion_share']))
                wl = wall(fd, pos['P3'], 710.0 * fs)
                out.append(row('collisions', 'abell520.m710', np.log(5.0e14 * L7 / m['M710_P3']), 0.55 / 5.0, wl['chi'], wl['chi_all'],
                               wl['gN_SI'] * 3.2408e-14, system='Abell 520', aperture_kpc=710.0 * fs, M_law=m['M710_P3'], M_wall=wl['M_wall'],
                               companion_share=wl['companion_share']))
                girardi = dict(P1=(811.0, 71.0, 278.0), P2=(749.0, 88.0, 186.0), P4=(579.0, 151.0, 523.0), P5=(668.0, 187.0, 570.0))
                for k, (o, lo, hi) in girardi.items():
                    v = m['speeds'][k]['500kpc']
                    q = settled([sets[k][0]], sets[k][1], law, 3000.0 * max(fs, 1.0), 500.0 * fs)
                    bl = ball(fd, pos[sets[k][2]], 500.0 * fs)
                    out.append(row('collisions', 'abell520.sigma_peaks', 2 * np.log(o / v), 2 * (hi if v > o else lo) / o, q['chi'], q['chi_all'],
                                   q['gN_SI'], system='Abell 520', peak=k, sigma_law=v, sigma_obs=o, chi_today=bl['chi'], chi_all_today=bl['chi_all']))
            else:
                conv = C10.CONVENTIONS['el_gordo']; L_ = C10.lens_factor('el_gordo')['lens']; fk = C10.lens_factor('el_gordo')['size']
                for R, M in zip(conv['aperture_Mpc'][:2], conv['aperture_M'][:2]):
                    kk = f'Mkim{int(R * 1000)}_com'
                    wl = wall(fd, pos['com'], R * 1000 * fk)
                    out.append(row('collisions', f'elgordo.m{int(R * 1000)}', np.log(M * L_ / m[kk]), conv['aperture_err'], wl['chi'], wl['chi_all'],
                                   wl['gN_SI'] * 3.2408e-14, system='El Gordo', aperture_kpc=R * 1000 * fk, M_law=m[kk], M_wall=wl['M_wall'],
                                   companion_share=wl['companion_share']))
                for w, (o, e) in (('NW', (1290.0, 134.0)), ('SE', (1089.0, 200.0))):
                    v = m['speeds'][w]['1000kpc']
                    q = settled([sets[w][0]], sets[w][1], law, 3000.0 * max(fs, 1.0), 1000.0 * fs)
                    bl = ball(fd, pos[sets[w][2]], 1000.0 * fs)
                    out.append(row('collisions', f'elgordo.sigma_{w.lower()}', 2 * np.log(o / v), 2 * e / o, q['chi'], q['chi_all'], q['gN_SI'],
                                   system='El Gordo', sigma_law=v, chi_today=bl['chi'], chi_all_today=bl['chi_all']))
    finally:
        V.kappa_map_v4 = orig
        FIELDS.clear()
    return out


# ------------------------------------------------------------------------------------------------ statistics
def _ranks(x):
    from scipy.stats import rankdata
    return rankdata(x)


def spearman(x, y):
    from scipy.stats import spearmanr
    if np.ptp(x) == 0 or np.ptp(y) == 0 or len(x) < 4:
        return None, None
    r = spearmanr(x, y)
    return float(r.statistic), float(r.pvalue)


def partial_spearman(x, y, controls):
    """Rank correlation of x and y after removing (linearly, in ranks) the controls; p from Student's t (n - 2 - k)."""
    from scipy.stats import t as T
    if np.ptp(x) == 0 or len(x) < 6:
        return None, None
    Z = np.c_[np.ones(len(x)), np.column_stack([_ranks(c) for c in controls])]
    rx, ry = _ranks(x), _ranks(y)
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    r = float(np.corrcoef(ex, ey)[0, 1])
    dof = len(x) - 2 - len(controls)
    tt = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * T.sf(abs(tt), dof))


def stratified(fams, key='chi', n_perm=20000, seed=25, partial=False):
    """Sum over families of the within-family rank correlation of residual and chi (or its partial, at fixed ln g_N);
    p: the share of permutations (chi shuffled within each family) with a sum at least as large."""
    rng = np.random.default_rng(seed)
    prep = []
    for f, rows in fams.items():
        x = np.array([r[key] for r in rows]); y = np.array([r['res'] for r in rows])
        if np.ptp(x) == 0 or len(rows) < 6:
            continue
        rx, ry = _ranks(x), _ranks(y)
        if partial:
            ctl = [np.log(np.array([r['gN_SI'] for r in rows]))]
            if f == 'xcop':
                ctl.append(np.array([r['x_R500'] for r in rows]))
            Z = np.c_[np.ones(len(x)), np.column_stack([_ranks(c) for c in ctl])]
            rx = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ry = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        rx = (rx - rx.mean()) / (rx.std() + 1e-300); ry = (ry - ry.mean()) / (ry.std() + 1e-300)
        prep.append((f, rx, ry))
    obs = sum(float(np.mean(rx * ry)) for _, rx, ry in prep)
    perm = np.zeros(n_perm)
    for _, rx, ry in prep:
        n = len(rx)
        idx = np.argsort(rng.random((n_perm, n)), axis=1)
        perm += (rx[idx] * ry[None, :]).mean(1)
    return dict(families=[f for f, _, _ in prep], sum_of_correlations=obs, p_one_sided=float((np.sum(perm >= obs) + 1) / (n_perm + 1)))


CHECK_LEVEL = ('lensing.kids_', 'lensing.mistele_', 'mw.v_sun', 'mw.mass_', 'mw.inner_share', 'mw.vertical_pull', 'mw.rc_outer',
               'dwarfs.', 'bullet.m250_', 'macs0025.', 'abell520.', 'elgordo.')


def check_level(rows):
    """One entry per graded check whose measurement is proportional to the pull: the mean residual and mean chi of its
    points (a check with several points, e.g. KiDS bins or the outer rotation curve, is one entry). Abell 520's four
    galaxy speeds form one graded check (close) and enter as one entry."""
    groups = {}
    for r in rows:
        if r['status'] not in ('pass', 'close', 'fail') or not r['id'].startswith(CHECK_LEVEL):
            continue
        groups.setdefault(r['id'], []).append(r)
    out = []
    for cid, rr in groups.items():
        res = np.array([r['res'] for r in rr])
        out.append(dict(id=cid, status=rr[0]['status'], family=rr[0]['family'], n=len(rr), res=float(np.mean(res)),
                        chi=float(np.mean([r['chi'] for r in rr])), chi_all=float(np.mean([r['chi_all'] for r in rr]))))
    # SPARC and X-COP (all pass): their offset checks, with the chi of the points they grade
    return out


def analyse(rows):
    from scipy.stats import mannwhitneyu
    fams = {}
    for r in rows:
        fams.setdefault(r['family'], []).append(r)
    per = {}
    for f, rr in fams.items():
        x = np.array([r['chi'] for r in rr]); xa = np.array([r['chi_all'] for r in rr], float); y = np.array([r['res'] for r in rr])
        lg = np.log(np.array([r['gN_SI'] for r in rr], float))
        ctl = [lg] + ([np.array([r['x_R500'] for r in rr])] if f == 'xcop' else [])
        ok = np.isfinite(lg).all()
        d = dict(n=len(rr), chi_min=float(x.min()), chi_median=float(np.median(x)), chi_max=float(x.max()),
                 chi_all_median=float(np.nanmedian(xa)), res_median=float(np.median(y)),
                 spearman=spearman(x, y), partial=partial_spearman(x, y, ctl) if ok else (None, None),
                 spearman_all=spearman(xa, y) if np.isfinite(xa).all() else (None, None),
                 partial_all=partial_spearman(xa, y, ctl) if ok and np.isfinite(xa).all() else (None, None),
                 spearman_lngN=spearman(lg, y) if ok else (None, None))
        if f == 'sparc':
            b = np.array([r['bulge'] for r in rr])
            d['bulge_points'] = int(b.sum())
            d['spearman_bulge_galaxies'] = spearman(x[b], y[b])
            d['partial_bulge_galaxies'] = partial_spearman(x[b], y[b], [lg[b]])
            d['res_median_chi0'] = float(np.median(y[x == 0])); d['res_median_chi_pos'] = float(np.median(y[x > 0]))
        per[f] = d
    fam_ok = {f: rr for f, rr in fams.items() if np.ptp([r['chi'] for r in rr]) >= 0.1}
    pooled = dict(chi=stratified(fams), chi_partial=stratified({f: rr for f, rr in fams.items() if np.isfinite([r['gN_SI'] for r in rr]).all()}, partial=True),
                  chi_all=stratified(fams, key='chi_all'), families_spanning_0p1=sorted(fam_ok))
    ck = check_level(rows)
    pas = np.array([c['chi'] for c in ck if c['status'] == 'pass']); non = np.array([c['chi'] for c in ck if c['status'] != 'pass'])
    pas_a = np.array([c['chi_all'] for c in ck if c['status'] == 'pass']); non_a = np.array([c['chi_all'] for c in ck if c['status'] != 'pass'])
    mw = mannwhitneyu(non, pas, alternative='greater'); mwa = mannwhitneyu(non_a, pas_a, alternative='greater')
    checks = dict(n_pass=int(pas.size), n_not=int(non.size), chi_median_pass=float(np.median(pas)), chi_median_not=float(np.median(non)),
                  mannwhitney_p=float(mw.pvalue), chi_all_median_pass=float(np.median(pas_a)), chi_all_median_not=float(np.median(non_a)),
                  mannwhitney_p_all=float(mwa.pvalue), rows=ck)
    a_ok = pooled['chi']['sum_of_correlations'] > 0 and pooled['chi']['p_one_sided'] < 0.01
    b_ok = all((per[f]['spearman'][0] or 0) > 0 for f in fam_ok) and pooled['chi_partial']['sum_of_correlations'] > 0 \
        and pooled['chi_partial']['p_one_sided'] < 0.05
    c_ok = checks['mannwhitney_p'] < 0.05
    verdict = dict(a_pooled=bool(a_ok), b_every_family_and_control=bool(b_ok), c_nonpasses_higher=bool(c_ok),
                   supported=bool(a_ok and b_ok and c_ok))
    return dict(per_family=per, pooled=pooled, checks=checks, verdict=verdict)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--skip-collisions', action='store_true')
    ap.add_argument('--reuse', action='store_true', help='reuse the collision rows of an earlier run in the same folder')
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    rows = []
    log('SPARC'); rows += sparc(law, ctx)
    log('X-COP'); rows += xcop(law, ctx)
    log('KiDS and Mistele'); rows += lensing(law)
    log('Milky Way'); rows += milky_way(law, ctx)
    log('dwarfs'); rows += dwarfs(law, ctx)
    prev = out / 'collision_rows.json'
    if args.reuse and prev.exists():
        rows += json.loads(prev.read_text())
    elif not args.skip_collisions:
        cr = collisions(law, ctx)
        prev.write_text(json.dumps(cr, indent=1) + '\n')
        rows += cr
    log('statistics')
    res = analyse(rows)
    payload = dict(source='code/frustration_v25.py', registration='round25-flow-frustration.md', law=law['name'],
                   constants=dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=law['u_kms']), **res, points=rows)
    (out / 'frustration_v25.json').write_text(json.dumps(payload, indent=1, default=float) + '\n')
    v = res['verdict']
    print(json.dumps(dict(verdict=v, pooled=res['pooled'], checks={k: res['checks'][k] for k in res['checks'] if k != 'rows'}), indent=1))
    for f, d in res['per_family'].items():
        print(f"{f:10s} n={d['n']:5d} chi {d['chi_min']:.3f}/{d['chi_median']:.3f}/{d['chi_max']:.3f}  res med {d['res_median']:+.3f}  "
              f"rho {d['spearman']}  partial {d['partial']}  | chi_all med {d['chi_all_median']:.3f} rho {d['spearman_all']} partial {d['partial_all']}")
    log(f"wrote {out / 'frustration_v25.json'}")


if __name__ == '__main__':
    main()
