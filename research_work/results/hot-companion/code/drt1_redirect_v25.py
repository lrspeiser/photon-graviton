"""Round 25, DRT-1's redirection against the data: how much room do the clusters leave for it?

DRT-1 (drt1-transport/README.md, main at ab1c027) gives a conservative way for opposing companion streams to turn one
another, (+x) + (-x) <-> (+y) + (-y), and shows that scattering keeps energy a system makes itself (a slab of optical
depth tau holds (1/2 + tau/12) of its own production against 1/2 when transparent) without piling an outside bath onto
gas. In the law's terms (DRT-1 section 1: S_hot = (4 pi G u / ell) U_hot, g_hot = -(4 pi G / ell) F_hot):
  * redirection alone keeps U_hot, so it leaves the pull's size sqrt(a (|g_N| + S)) unchanged and only moves its
    direction toward the ordered flow (the factor (g_N + g_hot)/(|g_N| + |g_hot|) tends to g_N/|g_N|);
  * redirection that slows the escape of a system's own hot glow raises U_hot in the steady state: S -> S (1 + tau/6),
    and it acts where streams oppose, so take tau = tau0 chi, chi = 1 - |g_N + g_hot| / (|g_N| + S) (round 25).
Two exploratory bounds (not adopted, nothing refitted in the suite):
  1. X-COP: fit tau0 together with u on the 12 clusters' hydrostatic masses; does holding back the hot glow where streams
     oppose help or hurt?
  2. What refitting u for that does where streams do not oppose: galaxy lensing (KiDS, its red/blue gap) and the
     bulge-dominated SPARC galaxies.
  3. The Bullet Cluster: the most redirection alone can do (the hot glow's direction fully merged into the ordered
     flow), and frustration-scaled retention with tau0 = 3, 10, 30: the smaller half's lensing mass (measured
     2.47-2.85 x 10^14 in the static law), the main cluster's, the leftover lensing on the gas and the peaks.

    python code/drt1_redirect_v25.py --output run-frustration-v25/drt1_redirect_v25.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, pickle, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(RESULTS / 'regression'))
T0 = time.monotonic()
TAU0 = (0.0, 3.0, 10.0, 30.0)


def log(msg):
    print(f'[{time.monotonic() - T0:5.0f} s] {msg}', flush=True)


def xcop(law, ctx):
    import law as L
    import t_clusters as TC
    from scipy.optimize import minimize
    G = L.G
    a, lam = law['a_code'], law['lam']
    cls = TC.static_clusters(ctx)
    W = {c['name']: L.shell_weights(c['Rk'], c['s']) for c in cls}

    def M(c, u, tau0):
        Rk, s = c['Rk'], c['s']
        gN = G * c['Mb'] / Rk ** 2
        src = L.heat_weight(np.sqrt(c['sig2_star_hse']), u) * c['dms']
        S = G * (W[c['name']] @ src) / Rk ** 2
        gh = G * np.array([src[s < x].sum() for x in Rk]) / Rk ** 2
        chi = 1 - (gN + gh) / (gN + S)
        return L.total(gN, S * (1 + tau0 * chi / 6), a, lam) * Rk ** 2 / G

    res = lambda u, t: np.array([np.log(c['Mh'] / M(c, u, t)) for c in cls])
    rows = []
    for t in TAU0:
        from scipy.optimize import minimize_scalar
        f = minimize_scalar(lambda lu: np.mean(res(10 ** lu, t) ** 2), bounds=(1.3, 4.0), method='bounded')
        r = res(10 ** f.x, t)
        rows.append(dict(tau0=t, u_refit=float(10 ** f.x), rms=float(np.sqrt(np.mean(r ** 2))), worst_radius=float(np.max(np.abs(r.mean(0)))),
                         mean_by_radius=r.mean(0).tolist()))
        log(f"X-COP tau0 {t:g}: u {10 ** f.x:.1f}, rms {rows[-1]['rms']:.3f}, worst radius {rows[-1]['worst_radius']:.3f}")
    best = minimize(lambda p: np.mean(res(10 ** p[0], max(p[1], -6.0)) ** 2), x0=[np.log10(law['u_kms']), 1.0], method='Nelder-Mead',
                    options=dict(xatol=1e-4, fatol=1e-7, maxiter=400))
    r = res(10 ** best.x[0], max(best.x[1], -6.0))
    return dict(rows=rows, joint_fit=dict(u=float(10 ** best.x[0]), tau0=float(best.x[1]), rms=float(np.sqrt(np.mean(r ** 2))),
                                          note='tau0 below zero would mean the clusters want less held-back glow where streams oppose'))


def side_effects(law, ctx, xrows):
    """What refitting u to suit the clusters does where the streams do not oppose (chi = 0 at the lensing radii of the
    suite's point lenses): the KiDS samples and the red/blue gap, and the bulge-dominated SPARC galaxies (their bulge
    heat held back by tau0 chi as well)."""
    import kids_heat_v12 as KH
    import kids_static_v11 as KS
    import collisions_v10 as C10
    import run as R
    import law as L
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
    a, lam = law['a_code'], law['lam']
    gals = ctx.sparc()
    bulgy = [g for g in gals if np.max(g['vb2'] / np.maximum(g['vbar2'], 1e-9)) > 0.5]

    def gal(g, u, tau0):
        if g['sigb'] == 0:
            return L.total(g['gN'], 0.0, a, lam)
        k = L.heat_weight(g['sigb'], u)
        S = L.scalar_sum(g['r'], g['sfine'], g['dmb'], k)
        gh = L.G * k * np.array([g['dmb'][g['sfine'] < x].sum() for x in g['r']]) / g['r'] ** 2
        chi = 1 - (g['gN'] + gh) / (g['gN'] + S)
        return L.total(g['gN'], S * (1 + tau0 * chi / 6), a, lam)

    rows = []
    for r in xrows:
        u, t = r['u_refit'], r['tau0']
        k = KH.kids(consts, u, law['reach_kpc'], f)
        rb = R.sparc_score(bulgy, lambda g: gal(g, u, t))[0]
        rows.append(dict(tau0=t, u=u, kids_all=k['all'], kids_red=k['red'], kids_blue=k['blue'], kids_gap=k['gap_model'],
                         kids_gap_observed=k['gap_observed'], sparc_bulges_rms_kms=rb))
        log(f"u {u:.1f} (tau0 {t:g}): KiDS all/red/blue {k['all']:+.3f}/{k['red']:+.3f}/{k['blue']:+.3f} dex, gap {k['gap_model']:.3f} "
            f"(observed {k['gap_observed']:.3f} +- 0.04); SPARC bulges {rb:.2f} km/s (MOND 30.35)")
    return rows


def bullet_args(law):
    """The Bullet Cluster's map inputs as the suite builds them (bullet_static_v11.run_bullet), cached."""
    import bullet_v4 as V
    import bullet_static_v11 as BS
    import collisions_v10 as C10
    cache = RESULTS / 'regression/cache/bullet_map_args_r12.pkl'
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
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
    out = dict(args=got['args'], kw=got['kw'], scrit=r['sigma_crit'], f=f, m250=dict(main=r['checks']['m250_main']['value'], sub=r['checks']['m250_sub']['value']))
    cache.write_bytes(pickle.dumps(out))
    return out


def bullet(law):
    import bullet_v3 as B
    import bullet_v4 as V
    import bullet_main_v5 as BM
    import collisions_v4 as C4
    import law as L
    ba = bullet_args(law)
    current, ghost_gas, ghost_stars, pos, consts = ba['args']
    kw = ba['kw']; fs, fl = ba['f']['size'], ba['f']['lens']; scrit = ba['scrit']
    n, dx, centre = kw.get('n', 192), kw.get('dx', 15.0), kw.get('centre', (360., 50.))
    fresh_kpc = kw.get('fresh_kpc', 30.0)
    G = L.G
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
    rho_now = rho_gas_now + B.build_density(current, pos, x, y, z, 'st')
    rho_ghost_gas = B.build_density(ghost_gas, pos, x, y, z, 'gas')
    rho_ghost_st = np.zeros_like(rho_now); krho = np.zeros_like(rho_now)
    for comp, prof in ghost_stars:
        rho_c = B.build_density({'c': comp}, pos, x, y, z, 'c')
        cx, cy = pos[comp['centre']]
        r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        rho_ghost_st += rho_c
        krho += L.k_from_sig2(np.interp(r3, prof[0], prof[1]) ** 2, u).astype(np.float32) * rho_c
        del r3, rho_c
    inv_r = B.Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    gN = -np.array(np.gradient(-G * inv_r(rho_now * dV), dx))
    Fv = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx))
    del inv_r
    Kin = V.vec_conv(n, dx, rmax=fresh_kpc)
    Fv = Fv + G * dV * (V.apply_vec(Kin, rho_gas_now, n) - V.apply_vec(Kin, rho_ghost_gas, n))
    del Kin
    inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    S = G * inv_r2(krho * dV)
    del inv_r2
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fmag = np.sqrt(np.sum(Fv ** 2, axis=0)) + 1e-30
    hmag = np.sqrt(np.sum(g_hot ** 2, axis=0))
    chi = 1 - np.sqrt(np.sum((Fv + g_hot) ** 2, axis=0)) / (Fmag + S)
    sep = {w: float(np.linalg.norm(pos[f'{w}_plasma'] - pos[f'{w}_bcg'])) for w in ('main', 'sub')}

    def measure(label, h):
        divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
        kap = (-divh / (4 * np.pi * G)).sum(axis=2) * dx / scrit
        dec = B.clowe_decomposition(x, y, kap, pos, rmax=1200.0 * fs)
        pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0 * fs)
        mp = min(pk, key=lambda d: d['dist_main_bcg']); spk = min(pk, key=lambda d: d['dist_sub_bcg'])
        m = {w: BM.mass_within(x, y, kap, pos[f'{w}_bcg'], 250.0 * fs, scrit) for w in ('main', 'sub')}
        row = dict(case=label, m250_main=m['main'], m250_sub=m['sub'], kappa_main=dec['main_bcg'], kappa_sub=dec['sub_bcg'],
                   gas_main=dec['main_plasma'], gas_sub=dec['sub_plasma'], peak_main_kpc=mp['dist_main_bcg'], peak_sub_kpc=spk['dist_sub_bcg'],
                   peak_limits_kpc=[0.25 * sep['main'], 0.25 * sep['sub']],
                   sub_target=[2.0e14 * fl, 2.3e14 * fl], main_target=[2.5e14 * fl, 2.8e14 * fl])
        log(f"Bullet {label}: m250 main {m['main']:.3g} sub {m['sub']:.3g} (targets {2.5e14 * fl:.3g}-{2.8e14 * fl:.3g}, {2.0e14 * fl:.3g}-{2.3e14 * fl:.3g}); "
            f"gas {dec['main_plasma']:.3f}/{dec['sub_plasma']:.3f}; peaks {mp['dist_main_bcg']:.0f}/{spk['dist_sub_bcg']:.0f} kpc")
        return row

    rows = []
    fr = np.exp(-mag / gd)
    for t in TAU0:
        St = S * (1 + t * chi / 6)
        h = gN + fr * np.sqrt(a * (Fmag + St)) * (Fv + g_hot) / (Fmag + hmag + 1e-30)
        rows.append(measure(f'held back where streams oppose, tau0 = {t:g}', h)); del h
    # redirection alone: the hot glow's energy kept, its direction merged into the ordered flow
    h = gN + fr * np.sqrt(a * (Fmag + S)) * Fv / Fmag
    rows.append(measure('redirection alone (hot direction merged into the ordered flow)', h)); del h
    return rows, dict(chi_median_within_300kpc=float(np.median(chi[(np.hypot(x[:, None, None] - pos['sub_bcg'][0], y[None, :, None] - pos['sub_bcg'][1]) + 0 * z[None, None, :] < 300 * fs)])))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    log('X-COP: hot glow held back where streams oppose')
    xc = xcop(law, ctx)
    log('side effects of refitting u: galaxy lensing and bulges')
    se = side_effects(law, ctx, xc['rows'])
    log('the Bullet Cluster')
    bl, extra = bullet(law)
    out = dict(source='code/drt1_redirect_v25.py', drt1='drt1-transport/README.md (main at ab1c027)', xcop=xc, side_effects=se,
               bullet=bl, bullet_extra=extra)
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    log(f'wrote {args.output}')


if __name__ == '__main__':
    main()
