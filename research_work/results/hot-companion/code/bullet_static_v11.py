"""Round 11: the Bullet Cluster (the suite's round-5 case) in the project's static distance law.

The round-5 model starts from Clowe et al. 2006's eight aperture masses of gas and stars (inside 100 kpc
of the two galaxy concentrations and the two gas peaks), their kpc-per-arcsec scale, and Barrena et al.
2002's galaxy speed spread. At z = 0.296, with the sources at z = 1 as in the round-5 model, at fixed
angle and flux (collisions_v10.factors):
    every length (positions, apertures, truncation and core radii, the speed aperture, the grid)  x size
    gas masses                                                                                    x gas
    star masses, including the outer stars' target from the Legacy Survey star count              x stars
    measured lensing masses inside a fixed angle                                                   x lens
    the convergence kappa (an observable)                                                          unchanged
The fresh companion sphere (u t) scales with the time since the crossing, which is the separation over
the (spectroscopic) shock speed: x size. Galaxy speeds are unchanged. The law's constants are those
refitted in the static law (run-xcop-static-v11); the round-9 constants are shown for comparison.

    python code/bullet_static_v11.py --output-dir run-bullet-static-v11
"""
from __future__ import annotations
import argparse, copy, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import bullet_v3 as B                             # noqa: E402
import bullet_v4 as V                             # noqa: E402
import bullet_main_v5 as BM                       # noqa: E402
import collisions_v4 as C4                        # noqa: E402
import collisions_v10 as C10                      # noqa: E402

U0 = 197.41001228101223


def static_sigma_crit(zl=0.296, zs=1.0):
    S, Ss = C10.static(zl), C10.static(zs)
    c = 299792.458
    return (c ** 2) / (4 * np.pi * B.G) * Ss['DA'] / (S['DA'] * C10.static_Dls(zl, zs)) / 1e3


def fit_components(pos, fs, fg, fst):
    """bullet_v3.fit_components with every length times fs and the eight aperture masses converted."""
    keys = ['main_bcg', 'main_plasma', 'sub_bcg', 'sub_plasma']
    dist = lambda a, b: float(np.linalg.norm(pos[a] - pos[b]))
    RT = dict(gas_main=2000. * fs, gas_sub=400. * fs, st_main=1500. * fs, st_sub=600. * fs)
    rap = 100.0 * fs

    def model_gas(q):
        Mm, rcm, Ms, rcs = np.exp(q)
        um = Mm / B.mass3d(B.rho_beta, (rcm,), RT['gas_main']); us = Ms / B.mass3d(B.rho_beta, (rcs,), RT['gas_sub'])
        return {k: um * B.aperture(B.rho_beta, (rcm,), RT['gas_main'], dist(k, 'main_plasma'), rap) +
                   us * B.aperture(B.rho_beta, (rcs,), RT['gas_sub'], dist(k, 'sub_plasma'), rap) for k in keys}

    def model_st(q):
        Mm, rsm, Ms, rss = np.exp(q)
        um = Mm / B.mass3d(B.rho_nfw, (rsm,), RT['st_main']); us = Ms / B.mass3d(B.rho_nfw, (rss,), RT['st_sub'])
        return {k: um * B.aperture(B.rho_nfw, (rsm,), RT['st_main'], dist(k, 'main_bcg'), rap) +
                   us * B.aperture(B.rho_nfw, (rss,), RT['st_sub'], dist(k, 'sub_bcg'), rap) for k in keys}

    rg = least_squares(lambda q: [np.log(model_gas(q)[k] / (B.OBS['gas'][k] * 1e12 * fg)) for k in keys],
                       np.log([1.5e14 * fg, 250. * fs, 2e13 * fg, 60. * fs]),
                       bounds=(np.log([1e12, 5 * fs, 1e11, 30 * fs]), np.log([1e16, 2000 * fs, 1e15, 400 * fs])))
    rs = least_squares(lambda q: [np.log(model_st(q)[k] / (B.OBS['stars'][k] * 1e12 * fst)) for k in keys],
                       np.log([3e12 * fst, 300. * fs, 1e12 * fst, 60. * fs]),
                       bounds=(np.log([1e10, 5 * fs, 1e10, 5 * fs]), np.log([1e15, 3000 * fs, 1e15, 3000 * fs])))
    return dict(gas_main=dict(kind='beta', M=float(np.exp(rg.x[0])), scale=float(np.exp(rg.x[1])), rt=RT['gas_main'], centre='main_plasma'),
                gas_sub=dict(kind='beta', M=float(np.exp(rg.x[2])), scale=float(np.exp(rg.x[3])), rt=RT['gas_sub'], centre='sub_plasma'),
                st_main=dict(kind='nfw', M=float(np.exp(rs.x[0])), scale=float(np.exp(rs.x[1])), rt=RT['st_main'], centre='main_bcg'),
                st_sub=dict(kind='nfw', M=float(np.exp(rs.x[2])), scale=float(np.exp(rs.x[3])), rt=RT['st_sub'], centre='sub_bcg'))


def pre_collision_models(comps, fs, ratio=8.0, star_frac=0.07, atm_core=150.0):
    """bullet_v4.pre_collision_models with its lengths times fs."""
    Mmain_b = comps['gas_main']['M'] + comps['st_main']['M']
    Msub_b = Mmain_b / ratio
    gas_pre = Msub_b / (1 + star_frac); stars_pre = Msub_b - gas_pre
    atm = gas_pre - comps['gas_sub']['M']
    sats = max(stars_pre - comps['st_sub']['M'], 0.0)
    main_gas = dict(comps['gas_main']); main_gas['M'] = comps['gas_main']['M'] - atm; main_gas['centre'] = 'main_bcg'
    ghost_gas = dict(gas_main_ghost=main_gas, gas_bullet_ghost=dict(comps['gas_sub'], centre='sub_bcg'),
                     gas_atm_ghost=dict(kind='beta', M=atm, scale=atm_core * fs, rt=1000.0 * fs, centre='sub_bcg'))
    sub_sats = dict(kind='nfw', M=sats, scale=atm_core * fs, rt=1000.0 * fs, centre='sub_bcg')
    return ghost_gas, sub_sats


def run_bullet(law, f, n=192, dx=15.0):
    fs, fg, fst, fl = f['size'], f['gas'], f['stars'], f['lens']
    pos = {k: v * fs for k, v in B.positions().items()}
    comps0 = fit_components(pos, fs, fg, fst)
    scrit = static_sigma_crit() if fs != 1.0 else B.sigma_crit()
    gas = dict(comps0['gas_main'], centre='main_bcg')
    cm = copy.deepcopy(comps0); inner = dict(cm['st_main'])
    target = BM.OBS_SIGMA[0]; rap_speed = BM.BARRENA_RAP * fs

    def los(M_out):
        stars = [inner] + ([BM.fixed_outer(inner, M_out, 800.0 * fs, 3000.0 * fs)] if M_out > 0 else [])
        r, sr, g, dms = BM.own_sigma([gas], stars, law, rmax=4000.0 * max(fs, 1.0))
        return BM.sigma_los_aperture(r, dms, sr, 0.0, rap_speed), (r, sr, stars)

    lo, hi = 11.0, 13.8
    if los(10 ** lo)[0] > target:
        M_out = 0.0
    else:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if los(10 ** mid)[0] < target else (lo, mid)
        M_out = 10 ** (0.5 * (lo + hi))
    _, (r, sr, stars) = los(M_out)
    ghost_gas, sub_sats = pre_collision_models(cm, fs)
    rs_s, ss_s, _ = V.own_sigma_multi([dict(cm['gas_sub']), ghost_gas['gas_atm_ghost']], [cm['st_sub'], sub_sats], law, rmax=3000.0 * max(fs, 1.0))
    cur = dict(cm); cur['st_main'] = stars[0]
    ghost_stars = [(dict(stars[0]), (r, sr))]
    if len(stars) > 1:
        cur['st_main_outer'] = stars[1]; ghost_stars.append((dict(stars[1]), (r, sr)))
    ghost_stars += [(dict(cm['st_sub']), (rs_s, ss_s))] + ([(sub_sats, (rs_s, ss_s))] if sub_sats['M'] > 0 else [])
    fresh = 30.0 * law['u_kms'] / U0 * fs
    x, y, Se, Sb = V.kappa_map_v4(cur, ghost_gas, ghost_stars, pos, law, n=n, dx=dx * fs, fresh_kpc=fresh, centre=tuple(np.array((360., 50.)) * fs))
    kap = Se / scrit
    dec = B.clowe_decomposition(x, y, kap, pos, rmax=1200.0 * fs)
    pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0 * fs)
    mp = min(pk, key=lambda d: d['dist_main_bcg']); sp = min(pk, key=lambda d: d['dist_sub_bcg'])
    m250 = {w: BM.mass_within(x, y, kap, pos[f'{w}_bcg'], 250.0 * fs, scrit) for w in ('main', 'sub')}
    sep = {w: float(np.linalg.norm(pos[f'{w}_plasma'] - pos[f'{w}_bcg'])) for w in ('main', 'sub')}
    zc = lambda v, o, e: (v - o) / e
    rng = lambda v, lo_, hi_, e: 0.0 if lo_ <= v <= hi_ else ((v - lo_) / e if v < lo_ else (v - hi_) / e)
    checks = dict(
        outer_stars=dict(value=M_out, target=f'{3.9 * fst:.2f}-{6.7 * fst:.2f} x 10^12', z=rng(M_out, 3.9e12 * fst, 6.7e12 * fst, 0.7e12 * fst)),
        kappa_main=dict(value=dec['main_bcg'], target='>= 0.36 +- 0.06', z=min(0.0, zc(dec['main_bcg'], 0.36, 0.06))),
        kappa_sub=dict(value=dec['sub_bcg'], target='>= 0.20 +- 0.05', z=min(0.0, zc(dec['sub_bcg'], 0.20, 0.05))),
        gas_main=dict(value=dec['main_plasma'], target='0.05 +- 0.06', z=zc(dec['main_plasma'], 0.05, 0.06)),
        gas_sub=dict(value=dec['sub_plasma'], target='0.02 +- 0.06', z=zc(dec['sub_plasma'], 0.02, 0.06)),
        peak_main=dict(value=mp['dist_main_bcg'], target=f"<= {0.25 * sep['main']:.0f} kpc", ok=bool(mp['dist_main_bcg'] <= 0.25 * sep['main'])),
        peak_sub=dict(value=sp['dist_sub_bcg'], target=f"<= {0.25 * sep['sub']:.0f} kpc", ok=bool(sp['dist_sub_bcg'] <= 0.25 * sep['sub'])),
        m250_main=dict(value=m250['main'], target=f'{2.5 * fl:.2f}-{2.8 * fl:.2f} x 10^14', z=rng(m250['main'], 2.5e14 * fl, 2.8e14 * fl, 0.15e14 * fl)),
        m250_sub=dict(value=m250['sub'], target=f'{2.0 * fl:.2f}-{2.3 * fl:.2f} x 10^14', z=rng(m250['sub'], 2.0e14 * fl, 2.3e14 * fl, 0.2e14 * fl)))
    return dict(checks=checks, masses=dict(gas_main=comps0['gas_main']['M'], gas_sub=comps0['gas_sub']['M'],
                                           st_main=comps0['st_main']['M'], st_sub=comps0['st_sub']['M'], outer_stars=M_out),
                sigma_crit=scrit, factors=dict(size=fs, gas=fg, stars=fst, lens=fl))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law, with_constants
    law9 = load_law('round9')
    st = json.loads((HERE.parent / 'run-xcop-static-v11/xcop_static_v11.json').read_text())['joint_refit_static'][-1]
    law11 = with_constants(law9, a_code=st['a_SI'] / L.KMS2_PER_KPC, lam=st['lam'], u_kms=st['u_kms'])
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    one = dict(size=1.0, gas=1.0, stars=1.0, lens=1.0)
    res = dict(experiment='round 11: the Bullet Cluster in the static distance law', factors=dict(size=f['size'], gas=f['gas'], stars=f['stars'], lens=f['lens']), runs={})
    for tag, lw, ff in (('standard distances, round-9 constants', law9, one), ('static distances, round-9 constants', law9, f),
                        ('static distances, static-refit constants', law11, f)):
        r = run_bullet(lw, dict(ff))
        res['runs'][tag] = r
        print(tag + ':', flush=True)
        for k, v in r['checks'].items():
            zz = v.get('z'); flag = ('pass' if abs(zz) <= 2 else 'close' if abs(zz) <= 3 else 'fail') if zz is not None else ('pass' if v['ok'] else 'fail')
            val = f"{v['value']:.3g}"
            print(f"   {k:12s} {val:>10s}  target {v['target']:24s} {('z %+.2f' % zz) if zz is not None else '':10s} {flag}", flush=True)
    res['seconds'] = time.monotonic() - t0
    (out / 'bullet_static_v11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'bullet_static_v11.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
