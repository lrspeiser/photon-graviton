"""The owner's self-energy repair (26 September 2026), checked and followed one step further with the suite's own models.

The repair: add V_self(T) = sqrt(a) T^(3/2) / (6 pi G) to the repaired field energy
    V0 = g^3 / (12 pi G a) - [(T - |d|) g + d . grad phi] / (4 pi G),
so that V0 + V_self = (g - h)^2 (g + 2h) / (12 pi G a) + (|d| g - d . grad phi) / (4 pi G) >= 0, h = sqrt(a T); T = S_hot,
d its directional part, g = |grad phi|. With DRT-1's normalisation (T = (4 pi G u / ell) U, d = (4 pi G / ell) F,
ell = a u / 2) the packets move with epsilon = dH/df, and for a packet of free energy eps0 and velocity v:
    epsilon = eps0 [1 + 2 (h - g) / a + (2 / (a u)) (g d_hat - grad phi) . v],
so in the aligned (round) case epsilon = n eps0 with n = 1 + 2 (h - g) / a. The total energy is bounded, but n < 0
wherever g - h > a/2, i.e. where the ordinary-source term Q exceeds a/4 + sqrt(a T). In a static field a packet keeps
its energy, so one born where n < 0 can never reach a region where n > 0 (there its energy would have to be positive):
the hot glow born there stays there, and, since adding packets there lowers the energy, the region fills until n = 0,
at T_eq = (Q - a/4)^2 / a. This script:
  1. checks the bound and the per-packet energy numerically;
  2. maps n over the suite's test points and finds each system's trapping radius;
  3. consequence A (trapping): galaxy lensing with the lenses' hot glow kept inside them; the clusters with their
     galaxies' glow kept where each galaxy's own pull exceeds a/4 + sqrt(a T);
  4. consequence B (fill-up): the condensate T_eq added inside the trapping zone, on SPARC and the Milky Way.
Exploratory; nothing adopted.

    python code/self_energy_check_v25.py --output run-self-energy-v25/self_energy_check_v25.json
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


def log(msg):
    print(f'[{time.monotonic() - T0:5.0f} s] {msg}', flush=True)


def n_factor(Q, T, a):
    """Per-packet energy factor, aligned case: 1 - 2 (sqrt(Q + T) - sqrt(T)) / sqrt(a)."""
    Q, T = np.asarray(Q, float), np.maximum(np.asarray(T, float), 0.0)
    return 1 - 2 * (np.sqrt(Q + T) - np.sqrt(T)) / np.sqrt(a)


def T_eq(Q, a):
    """The fill-up level: n = 0 at T = (Q - a/4)^2 / a where Q > a/4."""
    return np.maximum(np.asarray(Q, float) - a / 4, 0.0) ** 2 / a


# ------------------------------------------------------------------------------------------------ 1. checks
def algebra(rng):
    a, fg = 1.0, 1.0                                     # normalised: a = 1, 4 pi G = 1
    err = 0.0; low = np.inf
    for _ in range(100000):
        T = rng.exponential(3.0); dm = T * rng.uniform(0, 1)
        dh = rng.normal(size=3); dh /= np.linalg.norm(dh); d = dm * dh
        gv = rng.normal(size=3) * rng.exponential(2.0); g = np.linalg.norm(gv); h = np.sqrt(a * T)
        V = g ** 3 / (3 * fg * a) - ((T - dm) * g + d @ gv) / fg + 2 * np.sqrt(a) * T ** 1.5 / (3 * fg)
        F = (g - h) ** 2 * (g + 2 * h) / (3 * fg * a) + (dm * g - d @ gv) / fg
        err = max(err, abs(V - F) / max(1.0, abs(F))); low = min(low, V)
    rows = []
    for _ in range(200):                                 # per-packet energy: finite difference of H against the formula
        U = rng.exponential(1.0); nh = rng.normal(size=3); nh /= np.linalg.norm(nh)
        Fl = U * rng.uniform(0, 0.9) * nh; gv = rng.exponential(2.0) * rng.normal(size=3)
        v = rng.normal(size=3); v /= np.linalg.norm(v)

        def H(U_, F_):
            T_, d_ = 2 * U_, 2 * F_; g_ = np.linalg.norm(gv); dm_ = np.linalg.norm(d_)
            return U_ + g_ ** 3 / 3 - ((T_ - dm_) * g_ + d_ @ gv) + 2 * T_ ** 1.5 / 3
        e = 1e-7
        fd = (H(U + e, Fl + e * v) - H(U, Fl)) / e
        T, d = 2 * U, 2 * Fl; g = np.linalg.norm(gv); h = np.sqrt(T)
        pred = 1 + 2 * (h - g) + 2 * (g * d / np.linalg.norm(d) - gv) @ v
        rows.append(abs(fd - pred) / max(1.0, abs(pred)))
    return dict(factorisation_max_rel_error=float(err), smallest_completed_energy=float(low), packet_energy_max_rel_error=float(max(rows)),
                previous_minimum=[float(y - 4 / 3 * y ** 1.5) for y in (1, 10, 100)], repaired_minimum=[1.0, 10.0, 100.0])


# ------------------------------------------------------------------------------------------------ 2. n over the data
def data_map(law, ctx):
    import law as L
    import mw_dwarfs_v7 as D
    import t_dwarfs as TD
    import t_milky_way as TMW
    import milky_way_v7 as MW
    import mw_model as M
    import kids_heat_v12 as KH
    import collisions_v10 as C10
    import kids_static_v11 as KS
    K, aS, G = L.KMS2_PER_KPC, law['a_SI'], L.G
    P = json.loads((RESULTS / 'run-frustration-v25/frustration_v25.json').read_text())['points']
    out = {}
    for fam in ('sparc', 'xcop'):
        pts = [p for p in P if p['family'] == fam]
        Q = np.array([p['gN_SI'] for p in pts]); T = Q * np.array([p['S_over_gN'] for p in pts])
        n = n_factor(Q, T, aS)
        out[fam] = dict(n=len(pts), share_negative=float(np.mean(n < 0)), n_min=float(n.min()), n_median=float(np.median(n)))
    # SPARC bulges: the share of each bulge's mass born where n < 0 (its own glow counted in T)
    fr = []
    for g in ctx.sparc():
        if g['sigb'] <= 0:
            continue
        s, dmb = g['sfine'], g['dmb']
        Qs = np.interp(s, g['r'], g['gN']) * K
        k = L.heat_weight(g['sigb'], law['u_kms'])
        sm = np.r_[s[0] / 2, np.sqrt(s[1:] * s[:-1])]           # shells between the nodes (never on a node)
        Ts = L.scalar_sum(s, sm, dmb, k) * K
        fr.append(float(np.sum(dmb[n_factor(Qs, Ts, aS) < 0]) / max(np.sum(dmb), 1e-30)))
    out['sparc_bulge_mass_born_where_n_negative'] = dict(galaxies=len(fr), median=float(np.median(fr)), min=float(np.min(fr)))
    # the Milky Way, in the plane
    comps, grid, F = TMW.setup(ctx)
    ev = MW.Evaluator(comps, grid, F, law)
    parts = MW.models(comps)['M17']
    ev.run(parts, 'ours', reach=law['reach_kpc'])
    gR = sum(ff * F[kx][0] for kx, ff in parts.items()); gz = sum(ff * F[kx][1] for kx, ff in parts.items())
    hot = [MW.Scaled(comps[kx], ff) for kx, ff in parts.items() if kx in ('bulge', 'halo')]
    S, _, _ = M.heat_fields(hot, grid, law['u_kms'])
    jz = grid.nzh
    Rp = grid.R; Qp = np.hypot(gR[:, jz], gz[:, jz]) * K; Tp = S[:, jz] * K
    npl = n_factor(Qp, Tp, aS)
    Rz = float(Rp[np.where(npl < 0)[0].max()]) if (npl < 0).any() else 0.0
    out['milky_way'] = dict(zone_radius_kpc=Rz, n_at_sun=float(np.interp(8.2, Rp, npl)), n_at_3kpc=float(np.interp(3.0, Rp, npl)))
    # KiDS lenses: point lenses with S = 0.85 k g_N; the radius inside which n < 0, and n at the graded bins
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    kk = KH.heat(law['u_kms'])
    Mt = 10 ** 10.6 * f['stars']
    tabs = KS.tables(); gb = tabs['all']['gbar'] * f['stars'] / f['size'] ** 2; rel = tabs['all']['gbar'] >= 1e-13
    lens = {}
    for s_, k in kk.items():
        c = np.sqrt(1 + 0.85 * k) - np.sqrt(0.85 * k)
        gz0 = aS / (4 * c * c)                                  # g_N where n = 0
        r0 = np.sqrt(G * Mt * K / gz0)                          # kpc
        lens[s_] = dict(k=float(k), zone_radius_kpc=float(r0), n_at_bins_min=float(n_factor(gb[rel], 0.85 * k * gb[rel], aS).min()))
    out['kids_lenses'] = lens
    # dwarfs at the half-light radius
    data = json.loads((RESULTS / 'data/mw_dwarfs.json').read_text())
    mw = TD.galaxy_profile(law, ctx)
    mu, w = np.polynomial.legendre.leggauss(128)
    dn = []
    for d in data['dwarfs']:
        env = D.galaxy_env(d['D_gc_kpc'], mw, law)
        _, gi = D.plummer(2.0 * d['L_V'], d['r_h_pc'] / 1000.0)
        g_i = gi(d['r_h_pc'] / 1000.0); k_d = 3 * d['sigma_obs'] ** 2 / law['u_kms'] ** 2
        Q = np.hypot(env['gN'] * mu - g_i, env['gN'] * np.sqrt(1 - mu ** 2))
        dn.append(float(np.sum(w * n_factor(Q * K, (env['S'] + k_d * g_i) * K, aS)) / 2))
    out['dwarfs'] = dict(n_min=float(min(dn)), n_max=float(max(dn)))
    return out


# ------------------------------------------------------------------------------------------------ 3. consequence A: trapping
def kids_trapped(law):
    """The lenses' hot glow kept inside them: no hot brightness at the graded radii (k = 0)."""
    import kids_heat_v12 as KH
    import kids_static_v11 as KS
    import collisions_v10 as C10
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
    base = KH.kids(consts, law['u_kms'], law['reach_kpc'], f)
    trap = KH.kids(consts, law['u_kms'], law['reach_kpc'], f, k={s: 0.0 for s in KH.PAIRS})
    keys = ('all', 'red', 'blue', 'disc', 'bulge', 'gama', 'gap_model', 'gap_sersic_model')
    return dict(law={k: base[k] for k in keys}, glow_kept_inside={k: trap[k] for k in keys}, gap_observed=base['gap_observed'],
                gap_sersic_observed=base['gap_sersic_observed'])


def kids_extended(law):
    """The same, for extended lenses: stars outside the trapping zone send their glow out. Early types (red, bulge):
    Hernquist, R_e = 3.5 kpc; late types (blue, disc): exponential disks, R_d = 3.0 kpc (sizes typical of 10^10.6 Msun
    galaxies, times the static size factor). The zone (n = 0 for a point of the lens's mass, with the escaping share of
    its heat) and the escaping share are solved together; the lensing then uses k x (escaping share)."""
    import kids_heat_v12 as KH
    import kids_static_v11 as KS
    import collisions_v10 as C10
    import law as L
    G, K, aS = L.G, L.KMS2_PER_KPC, law['a_SI']
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    consts = dict(a_SI=aS, g_d_SI=law['g_d_SI'])
    kk = KH.heat(law['u_kms'])
    Mt = 10 ** 10.6 * f['stars']
    early = lambda r, Re=3.5 * f['size']: (r / (r + Re / 1.815)) ** 2
    late = lambda r, Rd=3.0 * f['size']: 1 - (1 + r / Rd) * np.exp(-r / Rd)
    prof = dict(red=early, bulge=early, blue=late, disc=late)
    keff, info = {}, {}
    for s_, k in kk.items():
        fe = 1.0
        for _ in range(100):
            ke = 0.85 * k * fe
            c = np.sqrt(1 + ke) - np.sqrt(ke)
            r0 = np.sqrt(G * Mt * K / (aS / (4 * c * c)))
            fe_new = 1 - prof[s_](r0)
            if abs(fe_new - fe) < 1e-10:
                break
            fe = 0.5 * fe + 0.5 * fe_new
        keff[s_] = k * fe; info[s_] = dict(k=float(k), escaping_share=float(fe), zone_radius_kpc=float(r0))
    r = KH.kids(consts, law['u_kms'], law['reach_kpc'], f, k=keff)
    keys = ('all', 'red', 'blue', 'disc', 'bulge', 'gama', 'gap_model', 'gap_sersic_model')
    return dict(lenses=info, result={k_: r[k_] for k_ in keys})


def xcop_trapped(law, ctx):
    """Cluster galaxies keep their own glow where their own pull exceeds a/4 + sqrt(a T): a Hernquist galaxy of mass Mg
    and scale r_H (R_e / 1.815) traps the share (r_t / (r_t + r_H))^2 of its stars, r_t = sqrt(G Mg / Qc) - r_H. The
    clusters' brightness then scales by the escaping share, solved together with T (which sets Qc)."""
    import law as L
    import run_v3 as R3
    import t_clusters as TC
    from scipy.optimize import minimize_scalar
    G, K = L.G, L.KMS2_PER_KPC
    a, lam = law['a_code'], law['lam']
    cls = TC.static_clusters(ctx)
    W = {c['name']: L.shell_weights(c['Rk'], c['s']) for c in cls}
    rows = []
    for Mg, Re in ((3e10, 3.0), (1e11, 5.0)):
        rH = Re / 1.815

        def esc(Tc):
            Qc = a / 4 + np.sqrt(a * np.maximum(Tc, 0.0))          # code units
            rt = np.maximum(np.sqrt(G * Mg / Qc) - rH, 0.0)
            return 1 - (rt / (rt + rH)) ** 2

        def M(c, u):
            gN = G * c['Mb'] / c['Rk'] ** 2
            S0 = G * (W[c['name']] @ (L.heat_weight(np.sqrt(c['sig2_star_hse']), u) * c['dms'])) / c['Rk'] ** 2
            S = S0.copy()
            for _ in range(60):                                    # T sets the trapping; the escaping glow sets T
                S = 0.5 * S + 0.5 * S0 * esc(S)
            return L.total(gN, S, a, lam) * c['Rk'] ** 2 / G, float(np.mean(S / S0))

        res = lambda u: np.array([np.log(c['Mh'] / M(c, u)[0]) for c in cls])
        r0 = res(law['u_kms'])
        fit = minimize_scalar(lambda lu: np.mean(res(10 ** lu) ** 2), bounds=(1.0, 4.0), method='bounded')
        rf = res(10 ** fit.x)
        share = float(np.mean([M(c, law['u_kms'])[1] for c in cls]))
        rows.append(dict(galaxy_mass=Mg, R_e_kpc=Re, glow_escaping_share=share, rms_u_fixed=float(np.sqrt(np.mean(r0 ** 2))),
                         mean_u_fixed=float(r0.mean()), u_refit=float(10 ** fit.x), rms_refit=float(np.sqrt(np.mean(rf ** 2))),
                         worst_radius_refit=float(np.max(np.abs(rf.mean(0))))))
        log(f"X-COP, galaxies {Mg:.0e} Msun: glow escaping {share:.2f}; rms {rows[-1]['rms_u_fixed']:.3f} (u fixed), "
            f"refit u {10 ** fit.x:.1f} -> rms {rows[-1]['rms_refit']:.3f}")
    return rows


# ------------------------------------------------------------------------------------------------ 4. consequence B: fill-up
def sparc_fill(law, ctx):
    import run as R
    import law as L
    from scipy.optimize import minimize, minimize_scalar
    gals = ctx.sparc()
    a0, lam0, u = law['a_code'], law['lam'], law['u_kms']
    gd = lam0 * a0

    def gfun(mode, a, lam):
        def f(g):
            gN = g['gN']
            S = 0.0
            if g['sigb'] > 0:
                S = L.scalar_sum(g['r'], g['sfine'], g['dmb'], L.heat_weight(g['sigb'], u))
            if mode == 'law':
                return L.total(gN, S, a, lam)
            Tc = T_eq(gN, a)                                       # the fill-up level (n = 0) inside the zone
            if mode == 'condensate':                               # added; the hot glow otherwise as in the law
                return L.total(gN, np.maximum(S, Tc), a, lam)
            inside = gN > a / 4 + np.sqrt(a * np.maximum(S, 0.0))  # 'condensate+trapping': no escaped bulge glow outside
            return L.total(gN, np.where(inside, np.maximum(S, Tc), 0.0), a, lam)
        return f
    rows = []
    for mode in ('law', 'condensate', 'condensate+trapping'):
        rms, msq = R.sparc_score(gals, gfun(mode, a0, lam0))
        ra = minimize_scalar(lambda la: R.sparc_score(gals, gfun(mode, 10 ** la, gd / 10 ** la))[1],
                             bounds=(np.log10(a0) - 0.5, np.log10(a0) + 0.5), method='bounded', options=dict(xatol=1e-5))
        a1 = 10 ** ra.x
        rms1, msq1 = R.sparc_score(gals, gfun(mode, a1, gd / a1))
        rb = minimize(lambda p: R.sparc_score(gals, gfun(mode, 10 ** p[0], 10 ** p[1]))[1], [np.log10(a0), np.log10(lam0)],
                      method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-7, maxiter=300))
        a2, l2 = 10 ** rb.x[0], 10 ** rb.x[1]
        rms2, msq2 = R.sparc_score(gals, gfun(mode, a2, l2))
        K = L.KMS2_PER_KPC
        med = lambda a_, l_: float(np.median(np.concatenate([np.log10(g['v'] ** 2 / g['r'] / gfun(mode, a_, l_)(g)) for g in gals])))
        rows.append(dict(mode=mode, rms_fixed=rms, msq_fixed=msq, median_dex_fixed=med(a0, lam0),
                         a_refit_SI=a1 * K, rms_a_refit=rms1, msq_a_refit=msq1,
                         a_both_SI=a2 * K, gd_both_SI=a2 * l2 * K, rms_both=rms2, msq_both=msq2, median_dex_both=med(a2, l2)))
        log(f"SPARC {mode}: fixed rms {rms:.2f} msq {msq:.4f}; a refit {a1 * K:.3e} rms {rms1:.2f}; a, g_d refit "
            f"{a2 * K:.3e}, {a2 * l2 * K:.3e}: rms {rms2:.2f} msq {msq2:.4f}")
    return rows


def milky_way_fill(law, ctx, sparc_rows):
    """The Milky Way (McMillan 2017 matter; the local census for the pull above the Sun, as the suite grades it) with the
    law's constants and with the constants SPARC prefers once the fill-up is in."""
    import milky_way_v7 as MW
    import mw_model as M
    import t_milky_way as TMW
    from law_config import with_constants
    import law as L
    comps, grid, F = TMW.setup(ctx)
    rc = json.loads((RESULTS / 'data/mw_rotation_curves.json').read_text())
    mods = MW.models(comps)
    orig = M.law_extra
    refit = next(r for r in sparc_rows if r['mode'] == 'condensate')
    a2 = refit['a_both_SI'] / L.KMS2_PER_KPC; lam2 = refit['gd_both_SI'] / refit['a_both_SI']
    cases = [('law', law, False), ('condensate, law constants', law, True),
             ('condensate, SPARC-refitted a and g_d', with_constants(law, a_code=a2, lam=lam2, u_kms=law['u_kms']), True)]
    out = {}
    for label, cst, fill in cases:
        def patched(gR, gz, S, hR, hz, consts, law='ours', a0_SI=1.2e-10, reach=None, grid=None, _fill=fill):
            if _fill and law == 'ours':
                S = np.maximum(S, T_eq(np.sqrt(gR ** 2 + gz ** 2), consts['a_code']))
            return orig(gR, gz, S, hR, hz, consts, law=law, a0_SI=a0_SI, reach=reach, grid=grid)
        M.law_extra = patched
        try:
            ev = MW.Evaluator(comps, grid, F, cst)
            r = ev.run(mods['M17'], 'ours', reach=cst['reach_kpc'])
            rl = ev.run(mods['local census'], 'ours', reach=cst['reach_kpc'])
        finally:
            M.law_extra = orig
        cmp = MW.compare_rc(r, rc)
        out[label] = dict(v_sun=float(np.interp(MW.R0, r['R'], r['v'])), rc_outer_mean_offset=float(np.mean([cmp[k]['outer_15_27'] for k in cmp])),
                          rc_inner_5_10=float(np.mean([cmp[k]['inner_5_10'] for k in cmp if cmp[k]['inner_5_10'] is not None])),
                          M20=r['M_eff']['20.0'], M50=r['M_eff']['50.0'], M100=r['M_eff']['100.0'], vertical_1p1_local_census=rl['Sigma_z']['1.1'])
        o = out[label]
        log(f"Milky Way, {label}: v_sun {o['v_sun']:.1f} (229-234 +- 7), outer offset {o['rc_outer_mean_offset']:+.1f} km/s, "
            f"M20 {o['M20']:.3g} (1.91e11), M50 {o['M50']:.3g} (4.5e11), vertical {o['vertical_1p1_local_census']:.1f} (69.8 +- 3.3)")
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    rng = np.random.default_rng(26)
    log('algebra'); alg = algebra(rng)
    log('n over the data'); dm = data_map(law, ctx)
    log('A: KiDS with the lenses\' glow kept inside'); kt = kids_trapped(law)
    log('A: KiDS with extended lenses'); ke = kids_extended(law)
    from law_config import with_constants
    ke['with_u_113'] = kids_extended(with_constants(law, u_kms=113.3))   # the speed X-COP prefers once its galaxies trap their glow
    log('A: X-COP with the galaxies\' glow kept inside them'); xt = xcop_trapped(law, ctx)
    log('B: SPARC with the fill-up'); sf = sparc_fill(law, ctx)
    log('B: the Milky Way with the fill-up'); mf = milky_way_fill(law, ctx, sf)
    out = dict(source='code/self_energy_check_v25.py', proposal='the owner\'s self-energy repair, 26 September 2026 (message; code not in the repository)',
               threshold_note='n < 0 where Q > a/4 + sqrt(a T); fill-up T_eq = (Q - a/4)^2 / a', algebra=alg, n_map=dm,
               A_kids=kt, A_kids_extended=ke, A_xcop=xt, B_sparc=sf, B_milky_way=mf)
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(dict(algebra=alg, n_map=dm, A_kids=kt), indent=1, default=float)[:4000])
    log(f'wrote {args.output}')


if __name__ == '__main__':
    main()
