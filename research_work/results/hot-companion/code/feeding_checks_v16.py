"""Round 16, step 1b: is the galaxies' preference for the feeding feedback real? Splits, bootstraps, the
surface-brightness trend, and the clusters.

Uses galaxy_g_feed from code/feeding_feedback_v16.py (beta = v_phase / u; beta = 0 is the law).
  1. SPARC: (a, g_d, beta) fitted on the training galaxies only (the round-3 split); the typical miss on the
     validation and test galaxies, against the adopted law.
  2. The outer-point residual against W/(aM) (how much of a galaxy's mass sits where the companion pulls), and
     against surface density M_sys / R_d^2: slopes with 400 bootstrap resamples of the galaxies, at beta = 0 and
     at the fitted beta.
  3. X-COP (static, deprojected; the suite's sample): the same feeding bookkeeping on each cluster's fine grid
     (gas and stars), a and g_d from SPARC at each beta, u refitted; the clusters' typical miss.

    python code/feeding_checks_v16.py --output run-feeding-feedback-v16/feeding_checks_v16.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                    # noqa: E402
import run as R                                    # noqa: E402
import common as C                                 # noqa: E402
from law_config import load_law                    # noqa: E402
from feeding_feedback_v16 import galaxy_g_feed     # noqa: E402

G = L.G


def score(gals, a, gd, u, beta, split=None):
    return R.sparc_score(gals, lambda g: galaxy_g_feed(g, a, u, gd / a, beta), split)


def fit_train(gals, u, a0, gd0, beta0, free_beta=True):
    tr = [g for g in gals if g['split'] == 'train']
    def stat(p):
        a, gd = np.exp(p[:2]); beta = p[2] if free_beta else beta0
        if beta < 0 or beta > 3: return 1e9
        return score(tr, a, gd, u, beta)[1]
    x0 = np.r_[np.log([a0, gd0]), [beta0]] if free_beta else np.log([a0, gd0])
    res = minimize(stat, x0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-10, maxiter=800))
    a, gd = np.exp(res.x[:2]); beta = float(res.x[2]) if free_beta else beta0
    return a, gd, beta


def outer_residuals(gals, a, gd, u, beta):
    w, res, sd = [], [], []
    for g in gals:
        gp, gc, ww = galaxy_g_feed(g, a, u, gd / a, beta, detail=True)
        w.append(ww[-1]); res.append(np.log10(g['v'][-1] ** 2 / g['r'][-1] / gp[-1]))
        sd.append(np.log10(g['Msys'] / (2 * np.pi * g['Rd'] ** 2)))
    return np.array(w), np.array(res), np.array(sd)


def slope_boot(x, y, n=400, seed=5):
    rng = np.random.default_rng(seed)
    A = np.c_[np.ones_like(x), x]
    b0 = np.linalg.lstsq(A, y, rcond=None)[0][1]
    bs = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        bs.append(np.linalg.lstsq(A[i], y[i], rcond=None)[0][1])
    return float(b0), float(np.std(bs))


def cluster_pred_feed(c, a, u, lam, beta, iters=80):
    """Enclosed mass at the six radii with the feeding bookkeeping on the cluster's own fine grid."""
    s = c['s']; dm = c['dmg'] + c['dms']; M = np.cumsum(dm)
    gN = G * M / s ** 2
    if 'W_ss' not in c:
        c['W_ss'] = L.shell_weights(s, s)
    src = L.heat_weight(np.sqrt(c['sig2_star_hse']), u) * c['dms']
    S = G * (c['W_ss'] @ src) / s ** 2
    f = np.exp(-gN / (lam * a))
    gc = f * np.sqrt(a * (gN + S))
    if beta > 0:
        for _ in range(iters):
            W = np.cumsum(np.r_[gc[0], 0.5 * (gc[1:] + gc[:-1])] * dm)
            new = f * np.sqrt(a * (gN * (1 + 2 * beta * W / (a * M)) + S))
            if np.max(np.abs(new - gc) / np.maximum(gc, 1e-30)) < 1e-9:
                gc = new; break
            gc = new
    g = gN + gc
    return np.interp(c['Rk'], s, g) * c['Rk'] ** 2 / G


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    gals = ctx.sparc(); u = law['u_kms']; a0, gd0 = law['a_code'], law['lam'] * law['a_code']
    out = dict(experiment='round 16: robustness of the feeding feedback', u_kms=u)
    # 1. train-only fits
    fits = {}
    for tag, free, b in (('law (beta = 0), refitted on the training galaxies', False, 0.0), ('beta free, fitted on the training galaxies', True, 0.2)):
        a, gd, beta = fit_train(gals, u, a0, gd0, b, free_beta=free)
        row = dict(a_SI=a * C.K_SI, g_d_SI=gd * C.K_SI, beta=beta)
        for split in (None, 'train', 'validation', 'test'):
            rms, msq = score(gals, a, gd, u, beta, split)
            row[f'rms_{split or "all"}'] = rms; row[f'msq_{split or "all"}'] = msq
        fits[tag] = row
        print(tag, {k: (round(v, 4) if isinstance(v, float) and v < 100 else f'{v:.4g}') for k, v in row.items()}, flush=True)
    out['train_fits'] = fits
    # 2. residual trends
    beta_fit = fits['beta free, fitted on the training galaxies']['beta']
    trends = {}
    for tag, beta in (('beta = 0', 0.0), (f'beta = {beta_fit:.2f} (fitted)', beta_fit), ('beta = 0.5', 0.5), ('beta = 1', 1.0)):
        # constants refitted on all galaxies at this beta
        def stat(p):
            a, gd = np.exp(p); return score(gals, a, gd, u, beta)[1]
        rr = minimize(stat, np.log([a0, gd0]), method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-10, maxiter=400))
        a, gd = np.exp(rr.x)
        w, res, sd = outer_residuals(gals, a, gd, u, beta)
        sw, ew = slope_boot(w, res); ss, es = slope_boot(sd, res)
        trends[tag] = dict(beta=beta, a_SI=a * C.K_SI, g_d_SI=gd * C.K_SI, rms_all=score(gals, a, gd, u, beta)[0],
                           slope_on_W=sw, slope_on_W_err=ew, slope_on_surface_density=ss, slope_on_surface_density_err=es,
                           median_outer_residual=float(np.median(res)))
        print(tag, {k: round(v, 4) if isinstance(v, float) and abs(v) < 100 else v for k, v in trends[tag].items()}, flush=True)
    out['trends'] = trends
    # 3. clusters
    import t_clusters as TC, run_v3 as R3
    cls = TC.static_clusters(ctx)
    crow = []
    for tag, v in trends.items():
        a, lam, beta = v['a_SI'] / C.K_SI, v['g_d_SI'] / v['a_SI'], v['beta']
        fn = lambda uu: np.mean(R3.resid(cls, lambda c: cluster_pred_feed(c, a, uu, lam, beta)) ** 2)
        r = minimize_scalar(lambda lu: fn(10 ** lu), bounds=(1.3, 4.0), method='bounded')
        ub = 10 ** r.x
        res_u = R3.resid(cls, lambda c: cluster_pred_feed(c, a, ub, lam, beta))
        res_law_u = R3.resid(cls, lambda c: cluster_pred_feed(c, a, u, lam, beta))
        crow.append(dict(tag=tag, beta=beta, u_refit=float(ub), rms_u_refit=float(R3.rms(res_u)), rms_u_adopted=float(R3.rms(res_law_u)),
                         worst_radius_mean=float(np.max(np.abs(res_u.mean(0))))))
        print('clusters', crow[-1], flush=True)
    out['clusters'] = crow
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
