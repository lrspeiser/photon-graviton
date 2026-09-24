"""Round 16, step 2b: the crossing heat. If the quiet store is opened by matter's velocity relative to the companion
flowing through it (code/stream_store_v16.py), then matter moving through ANOTHER system's flow is heated by that
relative motion. For a piece that sees its own system's flow (share 1 - p) and another system's (share p, the other
flow moving outward from the other system's centre at u, and carrying that system's velocity):
    k = 3 sigma^2/u^2  +  p (v^2 - 2 u v cos theta)/u^2,
v the piece's speed relative to the other system and theta the angle between its motion and the other flow
(receding: cos = +1; approaching: -1). Gas collides, so it takes no crossing heat (the law's rule).

The Bullet Cluster (round 5's case in the static distances, code/bullet_static_v11.py) is rerun with this term, and
nothing else changed. Each system's stars carry the heat they picked up; the heat emitted a time t ago now sits at a
distance u t from them (the companion's memory), so it enters through a distance-dependent kernel k(d):
    S_cross(x) = G sum_shells int rho_shell(x') k_shell(|x - x'|) / |x - x'|^2 dV'   (and the matching vector pull).
The shares p are the other system's companion intensity (|g_N| + S, from each system's own spherical model before
the collision) over the sum, at each shell of stars; the history is a straight pass: receding now at v_out since
pericentre (t_p = D_now / v_out), approaching before it at v_in, from where the other system's flow ends (its reach,
u x 13 Gyr).

    python code/crossing_heat_v16.py --output run-crossing-heat-v16/crossing_heat_v16.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, copy, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import bullet_v3 as B                             # noqa: E402
import bullet_v4 as V                             # noqa: E402
import bullet_main_v5 as BM                       # noqa: E402
import bullet_static_v11 as BS                    # noqa: E402
import collisions_v4 as C4                        # noqa: E402
import collisions_v10 as C10                      # noqa: E402

G = B.G
MYR_PER_KPC_PER_KMS = 977.8                       # 1 kpc at 1 km/s takes 977.8 Myr


def intensity_profile(gas_list, star_list, law, rmax=6000.0, n=600):
    """A system's own companion intensity I(r) = |g_N| + S (the law's cold plus hot parts), spherical, with its stars'
    dispersion from the Jeans equation in its own gravity (bullet_v4.own_sigma_multi)."""
    r, sig, g = V.own_sigma_multi(gas_list, star_list, law, rmax=rmax, n=n)
    dr = np.gradient(r)
    dmg = sum(V.shells_on(r, dr, c) for c in gas_list); dms = sum(V.shells_on(r, dr, c) for c in star_list)
    gN = G * np.cumsum(dmg + dms) / r ** 2
    S = G * (L.shell_weights(r, r) @ (L.k_from_sig2(sig ** 2, law['u_kms']) * dms)) / r ** 2
    return r, gN + S


def history(t_back, D_now, v_out, v_in, reach, b=0.0):
    """Separation D (kpc), speed v (km/s) and cos theta (+1 receding, -1 approaching) t_back Myr ago; b is the impact
    parameter (the closest approach), so D = sqrt(path^2 + b^2)."""
    vo, vi = v_out / MYR_PER_KPC_PER_KMS, v_in / MYR_PER_KPC_PER_KMS          # kpc per Myr
    tp = np.sqrt(max(D_now ** 2 - b ** 2, 0.0)) / vo
    path = np.where(t_back <= tp, vo * (tp - t_back), vi * (t_back - tp))
    D = np.sqrt(path ** 2 + b ** 2)
    v = np.where(t_back <= tp, v_out, v_in)
    cos = np.where(t_back <= tp, 1.0, -1.0)
    return D, v, cos, (D <= reach)


def make_kernel(u, I_other, I_own_at_shell, D_now, v_out, v_in, reach, kmax_d, b=0.0):
    """k(d) for a shell of stars whose own-system intensity is I_own_at_shell: the heat emitted a time d/u ago."""
    r_o, I_o = I_other
    def k_of_d(d):
        t = d / u * MYR_PER_KPC_PER_KMS
        D, v, cos, inside = history(t, D_now, v_out, v_in, reach, b)
        Io = np.exp(np.interp(np.log(np.maximum(D, r_o[0])), np.log(r_o), np.log(I_o)))
        p = np.where(inside, Io / (Io + I_own_at_shell), 0.0)
        k = p * np.maximum(v ** 2 - 2 * u * v * cos, 0.0) / u ** 2
        return np.where(d <= kmax_d, k, 0.0)
    return k_of_d


def vec_conv_k(n, dx, kfun):
    """Three convolutions giving G-less sum m k(d) (x' - x)/|x' - x|^3 (bullet_v4.vec_conv with a distance factor)."""
    m = 2 * n
    kk = ((np.arange(m) - (np.arange(m) >= n) * m) * dx).astype(np.float32)
    X = kk[:, None, None]; Y = kk[None, :, None]; Z = kk[None, None, :]
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2, dtype=np.float32)
    kr = np.where(r > 0, kfun(np.maximum(r, 1e-6).astype(np.float64)), 0.0).astype(np.float32)
    out = []
    for comp in (X, Y, Z):
        ker = np.where(r > 0, -comp * kr / np.maximum(r, 1e-6) ** 3, 0.0).astype(np.float32)
        out.append(np.fft.rfftn(ker)); del ker
    del r, kr
    return out


def kappa_map_cross(current, ghost_gas, ghost_stars, pos, consts, crossing=(), n=192, dx=15.0, centre=(360., 50.),
                    fresh_kpc=30.0, heat=True, memory=True, report=None):
    """bullet_v4.kappa_map_v4 plus the crossing heat: `crossing` is a list of (component, shell (r0, r1) about its
    centre, k_of_d). With crossing = () it is kappa_map_v4 line for line."""
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_now = B.build_density(current, pos, x, y, z, 'gas') + B.build_density(current, pos, x, y, z, 'st')
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
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
    F = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx)) if heat else 0.0
    del inv_r
    if memory and fresh_kpc > 0:
        Kin = V.vec_conv(n, dx, rmax=fresh_kpc)
        F = F + G * dV * (V.apply_vec(Kin, rho_gas_now, n) - V.apply_vec(Kin, rho_ghost_gas, n))
        del Kin
    inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    S = G * inv_r2(krho * dV) if heat else 0.0
    del inv_r2
    S_cross_total = 0.0
    for comps_, (r0, r1), kfun in crossing:
        rho_c = np.zeros_like(rho_now)
        for comp in comps_:
            rc = B.build_density({'c': comp}, pos, x, y, z, 'c')
            cx, cy = pos[comp['centre']]
            r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
            rho_c += (rc * ((r3 >= r0) & (r3 < r1))).astype(np.float32); del r3, rc
        conv = B.Conv(n, dx, lambda r: kfun(r) / r ** 2, float(kfun(np.array([0.5 * dx]))[0]) * B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
        Sx = G * conv(rho_c * dV); del conv
        Ks = vec_conv_k(n, dx, kfun)
        gx = G * dV * V.apply_vec(Ks, rho_c, n); del Ks
        S = S + Sx; g_hot = g_hot + gx
        S_cross_total = S_cross_total + Sx
        del rho_c
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fmag = np.sqrt(np.sum(F ** 2, axis=0)) + 1e-30
    hmag = np.sqrt(np.sum(g_hot ** 2, axis=0)) if heat else 0.0
    extra = np.exp(-mag / gd) * np.sqrt(a * (Fmag + S))
    h = gN + extra * (F + g_hot) / (Fmag + hmag + 1e-30)
    divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
    Sig_eff = (-divh / (4 * np.pi * G)).sum(axis=2) * dx
    Sig_b = rho_now.sum(axis=2) * dx
    if report is not None and crossing:
        k0 = n // 2
        report['S_cross_over_S_mid_plane'] = float(np.sum(S_cross_total[:, :, k0]) / max(np.sum(S[:, :, k0]), 1e-30))
    return x, y, Sig_eff, Sig_b


def run_bullet_cross(law, f, cross=None, n=192, dx=15.0):
    """bullet_static_v11.run_bullet with the crossing heat (cross = dict(v_out, v_in, shells)); cross=None reproduces it."""
    fs, fg, fst, fl = f['size'], f['gas'], f['stars'], f['lens']
    pos = {k: v * fs for k, v in B.positions().items()}
    comps0 = BS.fit_components(pos, fs, fg, fst)
    scrit = BS.static_sigma_crit() if fs != 1.0 else B.sigma_crit()
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
    ghost_gas, sub_sats = BS.pre_collision_models(cm, fs)
    sub_gas_list = [dict(cm['gas_sub']), ghost_gas['gas_atm_ghost']]
    sub_star_list = [cm['st_sub']] + ([sub_sats] if sub_sats['M'] > 0 else [])
    rs_s, ss_s, _ = V.own_sigma_multi(sub_gas_list, sub_star_list, law, rmax=3000.0 * max(fs, 1.0))
    cur = dict(cm); cur['st_main'] = stars[0]
    ghost_stars = [(dict(stars[0]), (r, sr))]
    if len(stars) > 1:
        cur['st_main_outer'] = stars[1]; ghost_stars.append((dict(stars[1]), (r, sr)))
    ghost_stars += [(dict(cm['st_sub']), (rs_s, ss_s))] + ([(sub_sats, (rs_s, ss_s))] if sub_sats['M'] > 0 else [])
    fresh = 30.0 * law['u_kms'] / BS.U0 * fs
    crossing = []; info = {}
    if cross:
        u = law['u_kms']; reach = u * 1.0227 * 13.0 * 1e3 / 1e3 * 1e3 / 1e3     # u x 13 Gyr in kpc (1 km/s x 1 Gyr = 1.0227 kpc)
        reach = u * 1.0227121650537077 * 13.0
        D_now = float(np.linalg.norm(pos['sub_bcg'] - pos['main_bcg']))
        main_gas_pre = dict(ghost_gas['gas_main_ghost'])
        I_main = intensity_profile([main_gas_pre], [dict(s) for s in stars], law)
        I_sub = intensity_profile(sub_gas_list, sub_star_list, law)
        bb = cross.get('b', 0.0) * fs
        info = dict(D_now_kpc=D_now, reach_kpc=reach, impact_parameter_kpc=bb,
                    t_since_pericentre_Myr=float(np.sqrt(max(D_now ** 2 - bb ** 2, 0.0)) / (cross['v_out'] / MYR_PER_KPC_PER_KMS)), shells=[])
        edges = [e * fs for e in cross.get('shells', (0.0, 50.0, 150.0, 400.0, 1e5))]
        for who, own_I, other_I, comp_list in (('sub', I_sub, I_main, [c for c in sub_star_list]),
                                               ('main', I_main, I_sub, [dict(s) for s in stars])):
            for r0, r1 in zip(edges[:-1], edges[1:]):
                rm = np.sqrt(r0 * r1) if r0 > 0 else r1 / 2
                I_own = float(np.exp(np.interp(np.log(rm), np.log(own_I[0]), np.log(own_I[1]))))
                kfun = make_kernel(u, other_I, I_own, D_now, cross['v_out'], cross['v_in'], reach, kmax_d=u * 13.0 * 1.0227121650537077,
                                   b=cross.get('b', 0.0) * fs)
                dd = np.array([10.0, 30.0, 100.0, 300.0]) * fs
                info['shells'].append(dict(system=who, stars_M=float(sum(c['M'] for c in comp_list)), shell_kpc=(r0, r1),
                                           k_at_d={f'{d:.0f}': float(kfun(np.array([d]))[0]) for d in dd}))
                crossing.append((comp_list, (r0, r1), kfun))
    rep = {}
    x, y, Se, Sb = kappa_map_cross(cur, ghost_gas, ghost_stars, pos, law, crossing=crossing, n=n, dx=dx * fs, fresh_kpc=fresh,
                                   centre=tuple(np.array((360., 50.)) * fs), report=rep)
    kap = Se / scrit
    dec = B.clowe_decomposition(x, y, kap, pos, rmax=1200.0 * fs)
    pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0 * fs)
    mp = min(pk, key=lambda d: d['dist_main_bcg']); sp = min(pk, key=lambda d: d['dist_sub_bcg'])
    m250 = {w: BM.mass_within(x, y, kap, pos[f'{w}_bcg'], 250.0 * fs, scrit) for w in ('main', 'sub')}
    sep = {w: float(np.linalg.norm(pos[f'{w}_plasma'] - pos[f'{w}_bcg'])) for w in ('main', 'sub')}
    return dict(outer_stars=M_out, kappa_main=dec['main_bcg'], kappa_sub=dec['sub_bcg'], gas_main=dec['main_plasma'],
                gas_sub=dec['sub_plasma'], peak_main=mp['dist_main_bcg'], peak_sub=sp['dist_sub_bcg'],
                peak_limits=dict(main=0.25 * sep['main'], sub=0.25 * sep['sub']),
                m250_main=m250['main'], m250_sub=m250['sub'], targets=dict(m250_main=(2.5e14 * fl, 2.8e14 * fl), m250_sub=(2.0e14 * fl, 2.3e14 * fl)),
                crossing=info, report=rep)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--quick', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    out = dict(experiment='round 16: the crossing heat in the Bullet Cluster', law=law['name'], factors=f, runs={})
    cases = [('the law (no crossing heat)', None),
             ('crossing heat: receding 3900 km/s since pericentre, approaching at 3000, impact parameter 150 kpc', dict(v_out=3900.0, v_in=3000.0, b=150.0)),
             ('crossing heat: slower (2700, 2000), b 150', dict(v_out=2700.0, v_in=2000.0, b=150.0)),
             ('crossing heat: faster (4700, 4000), b 150', dict(v_out=4700.0, v_in=4000.0, b=150.0)),
             ('crossing heat: 3900, 3000, head-on (b 0)', dict(v_out=3900.0, v_in=3000.0, b=0.0)),
             ('crossing heat: 3900, 3000, b 300', dict(v_out=3900.0, v_in=3000.0, b=300.0))]
    if args.quick: cases = cases[:2]
    for tag, cross in cases:
        r = run_bullet_cross(law, dict(f), cross=cross, n=192 if not args.quick else 128, dx=15.0 if not args.quick else 22.5)
        out['runs'][tag] = r
        print(f"[{time.monotonic() - t0:5.0f} s] {tag}: M250 sub {r['m250_sub']:.3e} (target {r['targets']['m250_sub'][0]:.2e}-{r['targets']['m250_sub'][1]:.2e}), "
              f"main {r['m250_main']:.3e} (target {r['targets']['m250_main'][0]:.2e}-{r['targets']['m250_main'][1]:.2e}); kappa sub {r['kappa_sub']:.3f}, main {r['kappa_main']:.3f}; "
              f"gas residuals {r['gas_main']:.3f}, {r['gas_sub']:.3f}; peaks {r['peak_main']:.0f}, {r['peak_sub']:.0f} kpc; report {r['report']}", flush=True)
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
