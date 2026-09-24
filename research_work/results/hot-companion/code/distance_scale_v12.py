"""Round 12: the distance law's scale alpha (1 + z = e^(alpha D)), fitted jointly.

alpha sets every distance the project infers from a redshift. Four data sets respond to it here:
* SPARC: the 'Hubble-flow' galaxies (f_D = 1 in Lelli et al. 2016, D = cz/H0 with H0 = 73) have
  D = ln(1 + z)/alpha in our law, so D x 73/(alpha c). Radii scale with D, the mass-model speeds with
  sqrt(D), and the measured pull v^2/r as 1/D. TRGB, Cepheid, Ursa Major and supernova distances do not
  depend on alpha.
* X-COP (z = 0.05-0.09): sizes, gas and stars converted to the static law at each alpha (xcop_static_v11).
* KiDS and Mistele (lenses at z = 0.25): the conversion of kids_heat_v12, with the lenses' heat measured.
* the local distance ladder, which fixes alpha c = H0 directly: SH0ES Cepheids 73.0 +- 1.0 (Riess et al.
  2022), CCHP TRGB 69.8 +- 1.7 (Freedman 2021), JAGB 67.8 +- 2.7 (Lee et al. 2024); the project's own
  164-group calibration gives the adopted alpha (74.6).
At each alpha the constants are refitted as in round 11 (a and g_d on SPARC, u on X-COP, alternated three
times), then the lensing tests are graded. The adopted geometry (D_A = D) is kept.

    python code/distance_scale_v12.py --output-dir run-distance-scale-v12
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import run as RUN                                 # noqa: E402
import run_v3 as R3                               # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import xcop_static_v11 as XS                      # noqa: E402
import kids_static_v11 as KS                      # noqa: E402
import kids_heat_v12 as KH                        # noqa: E402

C_KMS = 299792.458
LADDER = dict(SH0ES_cepheids=(73.04, 1.04), CCHP_TRGB=(69.8, 1.7), JAGB=(67.8, 2.7))


def load_sparc_scaled(alpha):
    """SPARC with the Hubble-flow galaxies' distances in the static law at this alpha (run.load_sparc)."""
    return RUN.load_sparc(alpha)


def xcop_static(scale):
    import common as C
    raw = json.loads((RUN.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    dep, _ = XS.deproject_xcop(C.Context(tier='quick', verbose=False).xcop(), raw)
    return [XS.to_static(c, raw[c['name']]['header']['z']) for c in dep]


def point(scale, law0):
    alpha0 = C10.ALPHA
    C10.ALPHA = alpha0 * scale
    try:
        gals = load_sparc_scaled(C10.ALPHA)
        cls = xcop_static(scale)
        a, lam, u = law0['a_code'], law0['lam'], R3.fit_u3(cls, law0['a_code'], law0['lam'])
        for _ in range(3):
            a, lam = RUN.fit_a_lam(gals, u, (a, lam))
            u = R3.fit_u3(cls, a, lam)
        sparc = RUN.sparc_score(gals, lambda g: RUN.galaxy_g(g, a, u, lam))[0]
        xr = float(R3.rms(R3.resid(cls, lambda c: R3.cluster_M3(c, a, u, lam))))
        consts = dict(a_SI=a * L.KMS2_PER_KPC, g_d_SI=a * lam * L.KMS2_PER_KPC)
        reach = u * 1.0227121650537077 * 13.0
        f = C10.factors(0.25, 0.75, KS.WMAP9)
        k = KH.kids(consts, u, reach, f); mi = KH.mistele(consts, u, reach, f)
    finally:
        C10.ALPHA = alpha0
    return dict(alpha_scale=scale, H0_like=scale * alpha0 * C_KMS, a_SI=consts['a_SI'], g_d_SI=consts['g_d_SI'], u_kms=float(u),
                sparc_rms_kms=float(sparc), xcop_rms=xr, kids={s: k[s] for s in ('all', 'blue', 'red', 'disc', 'bulge', 'gama', 'gap_model', 'gap_sersic_model')},
                kids_heat=k['k'], mistele=dict(LTG=mi['LTG']['rms_z'], ETG=mi['ETG']['rms_z']),
                n_hubble_flow=int(sum(g['f_D'] == 1 for g in gals)))


def sparc_scale_bootstrap(law0, scales, B=400, seed=20260924):
    """What alpha do the galaxies alone prefer? At each alpha refit a and g_d on SPARC (u from X-COP, as above),
    keep each galaxy's mean squared log residual (the statistic the constants are fitted with) and its rms in
    km/s, then resample the 149 galaxies and take the best alpha by each."""
    alpha0 = C10.ALPHA
    per_gal, per_kms, fits = [], [], []
    for sc in scales:
        C10.ALPHA = alpha0 * sc
        try:
            gals = load_sparc_scaled(C10.ALPHA); cls = xcop_static(sc)
            a, lam, u = law0['a_code'], law0['lam'], R3.fit_u3(cls, law0['a_code'], law0['lam'])
            for _ in range(3):
                a, lam = RUN.fit_a_lam(gals, u, (a, lam))
                u = R3.fit_u3(cls, a, lam)
            per_gal.append([float(np.mean(np.log(g['v'] ** 2 / g['r'] / RUN.galaxy_g(g, a, u, lam)) ** 2)) for g in gals])
            per_kms.append([float(np.sqrt(np.mean((np.sqrt(g['r'] * RUN.galaxy_g(g, a, u, lam)) - g['v']) ** 2))) for g in gals])
            hf = np.array([g['f_D'] == 1 for g in gals])
            fits.append(dict(alpha_scale=sc, a_SI=float(a * L.KMS2_PER_KPC), g_d_SI=float(a * lam * L.KMS2_PER_KPC), u_kms=float(u),
                             sparc_msq=float(np.mean(per_gal[-1])), sparc_rms_kms=float(np.mean(per_kms[-1])),
                             msq_hubble_flow=float(np.mean(np.array(per_gal[-1])[hf])), msq_other=float(np.mean(np.array(per_gal[-1])[~hf]))))
            print(f"  SPARC at alpha x{sc:.2f}: msq {fits[-1]['sparc_msq']:.5f}, {fits[-1]['sparc_rms_kms']:.3f} km/s, u {u:.1f}", flush=True)
        finally:
            C10.ALPHA = alpha0
    out = dict(fits=fits)
    for tag, M in (('msq', np.array(per_gal)), ('kms', np.array(per_kms))):
        rng = np.random.default_rng(seed); n = M.shape[1]
        best = [scales[int(np.argmin(M[:, rng.integers(0, n, n)].mean(1)))] for _ in range(B)]
        lo, med, hi = (float(x) for x in np.percentile(best, [16, 50, 84]))
        out[tag] = dict(best_scale=float(scales[int(np.argmin(M.mean(1)))]), bootstrap_16_50_84=[lo, med, hi],
                        H0_like_16_50_84=[x * alpha0 * C_KMS for x in (lo, med, hi)],
                        at_grid_edge=bool(np.argmin(M.mean(1)) in (0, len(scales) - 1)))
    return out


def joint_scale(out):
    """Combine the scales the supernovae (sn_scale_v12.py, all redshifts, Cepheid calibration) and SPARC's
    Hubble-flow galaxies (the fit statistic's bootstrap) prefer, by inverse variance."""
    sn = next(w for w in json.loads((out / 'sn_scale_v12.json').read_text())['windows'] if w['window'].startswith('all'))
    sp = json.loads((out / 'distance_scale_v12.json').read_text())['sparc_preferred']['msq']
    s1, e1 = sn['best_scale'], (sn['scale_1sigma'][1] - sn['scale_1sigma'][0]) / 2
    s2, e2 = sp['best_scale'], (sp['bootstrap_16_50_84'][2] - sp['bootstrap_16_50_84'][0]) / 2
    w1, w2 = e1 ** -2, e2 ** -2
    return dict(supernovae=[s1, e1], sparc=[s2, e2], combined=[(w1 * s1 + w2 * s2) / (w1 + w2), (w1 + w2) ** -0.5])


def adopt(scale, law0, out, iterations=8):
    """The constants at one scale, a and g_d on SPARC (in the static law) and u on X-COP, alternated to convergence."""
    alpha0 = C10.ALPHA
    C10.ALPHA = alpha0 * scale
    try:
        gals = load_sparc_scaled(C10.ALPHA); cls = xcop_static(scale)
        a, lam = law0['a_code'], law0['lam']; u = R3.fit_u3(cls, a, lam)
        log = []
        for it in range(iterations):
            a, lam = RUN.fit_a_lam(gals, u, (a, lam)); u_new = R3.fit_u3(cls, a, lam)
            rms_kms, msq = RUN.sparc_score(gals, lambda g: RUN.galaxy_g(g, a, u_new, lam))
            log.append(dict(iteration=it + 1, a_SI=float(a * L.KMS2_PER_KPC), lam=float(lam), g_d_SI=float(a * lam * L.KMS2_PER_KPC), u_kms=float(u_new),
                            sparc_rms_kms=float(rms_kms), sparc_msq=float(msq), xcop_rms=float(R3.rms(R3.resid(cls, lambda c: R3.cluster_M3(c, a, u_new, lam))))))
            print(f"  iteration {it + 1}: a {log[-1]['a_SI']:.5e}, g_d {log[-1]['g_d_SI']:.5e}, u {u_new:.3f}; SPARC {rms_kms:.3f} km/s, X-COP {log[-1]['xcop_rms']:.4f}", flush=True)
            done = abs(u_new - u) < 1e-3 and it > 0 and abs(log[-1]['a_SI'] / log[-2]['a_SI'] - 1) < 1e-5
            u = u_new
            if done: break
    finally:
        C10.ALPHA = alpha0
    res = dict(experiment="round 12: the adopted distance scale and the constants refitted there", alpha_scale=scale,
               alpha_per_Mpc=alpha0 * scale, H0_like=alpha0 * scale * C_KMS, basis=joint_scale(out), joint_refit=log,
               note='SPARC\'s Hubble-flow galaxies (f_D = 1) at D = ln(1 + z)/alpha; X-COP deprojected and converted at this alpha')
    (out / 'adopted_v12.json').write_text(json.dumps(res, indent=1) + '\n')
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--adopt', type=float, default=None, help='only refit the constants at this scale (x alpha0) and write adopted_v12.json')
    ap.add_argument('--scales', default='0.80,0.85,0.90,0.95,1.00,1.05')
    ap.add_argument('--fine', default='0.60:1.05:0.01', help='the grid for the SPARC-alone preference, lo:hi:step')
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    law0 = load_law('round11')
    if args.adopt is not None:
        b = joint_scale(out)
        print(f"supernovae x{b['supernovae'][0]:.4f} +- {b['supernovae'][1]:.4f}; SPARC x{b['sparc'][0]:.3f} +- {b['sparc'][1]:.3f};"
              f" combined x{b['combined'][0]:.4f} +- {b['combined'][1]:.4f}; adopting x{args.adopt}", flush=True)
        adopt(args.adopt, law0, out)
        return
    rows = []
    for s in [float(x) for x in args.scales.split(',')]:
        r = point(s, law0); rows.append(r)
        kk = r['kids']
        print(f"alpha x{s:.2f} (H0-like {r['H0_like']:.1f}): a {r['a_SI']:.3e}, g_d {r['g_d_SI']:.3e}, u {r['u_kms']:.1f}; SPARC {r['sparc_rms_kms']:.2f},"
              f" X-COP {r['xcop_rms']:.3f}; KiDS " + ', '.join(f"{x} {kk[x]:+.3f}" for x in ('all', 'blue', 'red', 'disc', 'bulge', 'gama')) +
              f"; gaps {kk['gap_model']:.3f} / {kk['gap_sersic_model']:.3f}; Mistele {r['mistele']['LTG']:.2f} / {r['mistele']['ETG']:.2f}", flush=True)
    # where KiDS centres: zero of the 'all' median offset, by linear interpolation
    xs = np.array([r['alpha_scale'] for r in rows]); ys = np.array([r['kids']['all'] for r in rows])
    order = np.argsort(xs); xs, ys = xs[order], ys[order]
    cross = [float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i])) for i in range(len(xs) - 1) if ys[i] * ys[i + 1] <= 0]
    lo, hi, step = (float(x) for x in args.fine.split(':'))
    fine = [round(x, 3) for x in np.arange(lo, hi + step / 2, step)]
    boot = sparc_scale_bootstrap(law0, fine)
    for tag, what in (('msq', 'mean squared log residual (the fit statistic)'), ('kms', 'rms in km/s')):
        bb = boot[tag]; b = bb['H0_like_16_50_84']
        print(f"SPARC alone, by {what}: best alpha x{bb['best_scale']:.2f}{' (grid edge)' if bb['at_grid_edge'] else ''};"
              f" bootstrap 16-50-84% x{bb['bootstrap_16_50_84']}, H0-like {b[0]:.1f}-{b[1]:.1f}-{b[2]:.1f}", flush=True)
    res = dict(experiment="round 12: the distance law's scale, with SPARC's Hubble-flow distances, X-COP and KiDS moved together",
               alpha0_per_Mpc=C10.ALPHA, H0_like_adopted=C10.ALPHA * C_KMS, ladder=LADDER, rows=rows,
               kids_all_zero_at_scale=cross, sparc_preferred=boot, seconds=time.monotonic() - t0)
    (out / 'distance_scale_v12.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"KiDS (all) centres at alpha x {cross}; wrote {out / 'distance_scale_v12.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
