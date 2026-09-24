"""Round 16, step 1: the power a pulled body feeds into the companion, and what the galaxies allow.

A body pulled by the companion is pulled by the recoil of the power it feeds the passing wave (round 13, §23.3):
F = P / v_phase, since a wave carries momentum k/omega = 1/v_phase per unit energy. So every kilogram pulled with
the companion's extra acceleration g_c feeds P/m = g_c v_phase into the companion (round 10, §20.1, condition 3).
That power joins the outgoing stream: the stream's energy grows by l per kilogram (every kilogram's own feed,
a = 2 l / u) plus g_c v_phase per kilogram (the feeding). With beta = v_phase / u, the ordered companion's
intensity at radius r becomes (spherical bookkeeping, as for the law's |g_N| = G M(r)/r^2)
    I(r) = |g_N(r)| (1 + 2 beta W(r) / (a M(r))),      W(r) = integral_0^r g_c dM,
and the extra pull g_c = exp(-|g_N|/g_d) sqrt(a (I + S)) is solved self-consistently from the centre out. beta = 0
is the law. Round 10 estimated that beta = 1 would spoil the tight v^4 = G M a and asked for beta <~ 0.05; here the
149 SPARC galaxies measure it: for each beta, a and g_d are refitted (the round-3 statistic, the mean squared log
residual) with u held, and the typical speed miss and the residual trend with how much of a galaxy's mass sits in
the companion's regime are reported.

    python code/feeding_feedback_v16.py --output run-feeding-feedback-v16/feeding_feedback_v16.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                    # noqa: E402
import run as R                                    # noqa: E402
import common as C                                 # noqa: E402
from law_config import load_law                    # noqa: E402

G = L.G


def galaxy_g_feed(g, a, u, lam, beta, iters=60, detail=False):
    gN = g['gN']; r = g['r']
    S = 0.0 if g['sigb'] == 0 else L.scalar_sum(r, g['sfine'], g['dmb'], L.heat_weight(g['sigb'], u))
    f = np.exp(-np.abs(gN) / (lam * a))
    M = np.maximum.accumulate(np.maximum(g['vbar2'] * r / G, 0.0))     # mass inside r, spherical bookkeeping
    dM = np.diff(np.r_[0.0, M])
    gc = f * np.sqrt(a * (np.abs(gN) + S))
    if beta > 0:
        for _ in range(iters):
            gc_mid = np.r_[gc[0], 0.5 * (gc[1:] + gc[:-1])]           # the shell between r_{i-1} and r_i
            W = np.cumsum(gc_mid * dM)
            I = np.abs(gN) * (1 + 2 * beta * W / (a * np.maximum(M, 1e-30)))
            new = f * np.sqrt(a * (I + S))
            if np.max(np.abs(new - gc) / np.maximum(gc, 1e-30)) < 1e-10:
                gc = new; break
            gc = new
    if detail:
        W = np.cumsum(np.r_[gc[0], 0.5 * (gc[1:] + gc[:-1])] * dM)
        return gN + gc, gc, W / (a * np.maximum(M, 1e-30))
    return gN + gc


def fit(gals, u, beta, a0, lam0):
    gd0 = lam0 * a0
    def stat(p):
        a, gd = np.exp(p)
        return R.sparc_score(gals, lambda g: galaxy_g_feed(g, a, u, gd / a, beta))[1]
    res = minimize(stat, np.log([a0, gd0]), method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-9, maxiter=400))
    a, gd = np.exp(res.x)
    rms, msq = R.sparc_score(gals, lambda g: galaxy_g_feed(g, a, u, gd / a, beta))
    return dict(beta=beta, a_SI=float(a * C.K_SI), g_d_SI=float(gd * C.K_SI), rms_kms=rms, msq=msq, a_code=float(a), lam=float(gd / a))


def trend(gals, fitres, u):
    """Outermost-point residual log10(v_obs^2 / v_pred^2) against the share of the galaxy's mass whose companion
    feeds (W/(aM) at the last point, which is what the feeding adds): the slope of the residual on it."""
    a, lam, beta = fitres['a_code'], fitres['lam'], fitres['beta']
    x, y = [], []
    for g in gals:
        gp, gc, w = galaxy_g_feed(g, a, u, lam, max(beta, 0.0), detail=True)
        x.append(w[-1]); y.append(np.log10(g['v'][-1] ** 2 / g['r'][-1] / gp[-1]))
    x, y = np.array(x), np.array(y)
    A = np.c_[np.ones_like(x), x]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return dict(slope=float(coef[1]), intercept=float(coef[0]), w_median=float(np.median(x)), w_p10=float(np.percentile(x, 10)),
                w_p90=float(np.percentile(x, 90)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    gals = ctx.sparc(); u = law['u_kms']
    base = R.sparc_score(gals, lambda g: R.galaxy_g(g, law['a_code'], u, law['lam']))
    same = R.sparc_score(gals, lambda g: galaxy_g_feed(g, law['a_code'], u, law['lam'], 0.0))
    print(f'adopted law: rms {base[0]:.4f} km/s, statistic {base[1]:.6f}; this code at beta = 0: {same[0]:.4f}, {same[1]:.6f}', flush=True)
    rows = []
    a0, lam0 = law['a_code'], law['lam']
    for beta in (0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
        fr = fit(gals, u, beta, a0, lam0)
        fr['trend'] = trend(gals, fr, u)
        held = R.sparc_score(gals, lambda g: galaxy_g_feed(g, law['a_code'], u, law['lam'], beta))
        fr['rms_constants_held'] = held[0]; fr['msq_constants_held'] = held[1]
        rows.append(fr); a0, lam0 = fr['a_code'], fr['lam']
        print(f"beta {beta:5.2f}: a {fr['a_SI']:.3e}, g_d {fr['g_d_SI']:.3e}; rms {fr['rms_kms']:.3f} km/s (statistic {fr['msq']:.5f}); "
              f"constants held: {held[0]:.3f}; outer residual slope on W/(aM) {fr['trend']['slope']:+.3f} (W/aM median {fr['trend']['w_median']:.2f})", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 16: the feeding feedback on the 149 SPARC galaxies', u_kms=u,
                                           adopted=dict(rms=base[0], msq=base[1]), rows=rows, seconds=time.monotonic() - t0), indent=1) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
