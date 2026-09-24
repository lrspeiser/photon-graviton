"""Round 14: is the law's release factor the same thing as "absorb loud waves"?

The law holds the companion back where the Newtonian pull is strong: extra = R(|g_N|) sqrt(a (|g_N| + S)) with
R = exp(-|g_N|/g_d). A body that both sings (a self-sustained oscillation locked a quarter cycle ahead: it feeds a
passing wave in proportion to its amplitude |E|) and absorbs (a passive part: it takes energy in proportion to |E|^2)
is pulled by (k/2)(a0 |E| - c0 |E|^2) = (k/2) a0 |E| (1 - |E|/E_s): the same pull, held back by 1 - |E|/E_s. The
ordered companion's intensity is |E|^2 ~ |g_N| (round 11: J = (l/4 pi G)(-g_N)), so the singer-absorber gives
R = 1 - sqrt(|g_N|/g_s) (and a push beyond g_s). Here the 149 SPARC galaxies judge the shapes, each with its two
constants (a and the scale) refitted, u and the heat held at the round-12 law:
    exp          R = exp(-g/g_d)                        (the law)
    absorber     R = 1 - sqrt(g/g_s), pushing beyond g_s (the lone singer-absorber, as derived)
    absorber+    R = max(0, 1 - sqrt(g/g_s))            (the same, stopped at zero)
    saturating   R = 1 - x / (1 + x^2),  x = sqrt(g/g_s) (an absorber that saturates in loud waves: no push; the
                                                          pull comes back in very strong waves)
    balance      R = clip(A (1/x - x), 0, 1),  x = sqrt(g/g_s) (the singer with a power balance, round 14: its pull is
                                                          the supply minus the loss, (k/2)(P_s - c_L |E|^2), and the
                                                          full round-10 pull (k/2)|E| in quiet waves; A = sqrt(P_s c_L),
                                                          a third constant, stopped at zero)
    python code/release_shape_v14.py --output run-coherent-force-v13/release_shape_v14.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import run as R                                   # noqa: E402
from law_config import load_law, K_SI             # noqa: E402

SHAPES = {
    'exp': lambda g, s: np.exp(-g / s),
    'absorber': lambda g, s: 1.0 - np.sqrt(g / s),
    'absorber+': lambda g, s: np.maximum(0.0, 1.0 - np.sqrt(g / s)),
    'saturating': lambda g, s: 1.0 - np.sqrt(g / s) / (1.0 + g / s),
    'balance': lambda g, s, A=1.0: np.clip(A * (np.sqrt(s / np.maximum(g, 1e-300)) - np.sqrt(g / s)), 0.0, 1.0),
}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    law = load_law('round12')
    u = law['u_kms']; a0 = law['a_code']; gd0 = law['lam'] * a0
    gals = R.load_sparc(alpha=law['alpha_per_Mpc'])
    pre = []
    for g in gals:
        S = 0.0 if g['sigb'] == 0 else L.scalar_sum(g['r'], g['sfine'], g['dmb'], L.heat_weight(g['sigb'], u))
        pre.append((np.abs(g['gN']), g['gN'], S))

    def score(shape, a, s, A=None):
        rm, ll = [], []
        for g, (ag, gN, S) in zip(gals, pre):
            Rf = SHAPES[shape](ag, s) if A is None else SHAPES[shape](ag, s, A)
            p = gN + Rf * np.sqrt(a * (ag + S))
            p = np.maximum(p, 1e-30)
            rm.append(np.sqrt(np.mean((np.sqrt(g['r'] * p) - g['v']) ** 2)))
            ll.append(np.mean(np.log(g['v'] ** 2 / g['r'] / p) ** 2))
        return float(np.mean(rm)), float(np.mean(ll))

    rows = []
    for shape in SHAPES:
        three = shape == 'balance'
        starts = [(np.log(a0), np.log(gd0 * f)) + ((np.log(A),) if three else ()) for f in (0.3, 1.0, 3.0, 10.0, 30.0)
                  for A in ((0.3, 1.0, 3.0) if three else (1.0,))]
        best = None
        for x0 in starts:
            fn = (lambda x: score(shape, np.exp(x[0]), np.exp(x[1]), np.exp(x[2]))[1]) if three else (lambda x: score(shape, np.exp(x[0]), np.exp(x[1]))[1])
            r = minimize(fn, x0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-9, maxiter=3000))
            if best is None or r.fun < best.fun: best = r
        a, s = np.exp(best.x[:2]); A = float(np.exp(best.x[2])) if three else None
        rm, ll = score(shape, a, s, A)
        # where the hold halves the pull, as a Newtonian acceleration
        gg = np.geomspace(1e-3, 1e3, 20001) * s
        Rg = SHAPES[shape](gg, s) if A is None else SHAPES[shape](gg, s, A)
        half = float(gg[np.argmin(np.abs(Rg - 0.5))] * K_SI)
        rows.append(dict(shape=shape, a_SI=float(a * K_SI), scale_SI=float(s * K_SI), A=A, typical_miss_kms=rm, fit_statistic=ll, half_hold_at_SI=half))
        print(f"{shape:10s}: a {a * K_SI:.4e} m/s^2, scale {s * K_SI:.4e} m/s^2" + (f", A {A:.3f}" if A else '') + f" (halves the pull at g_N = {half:.3e}); "
              f"typical miss {rm:.3f} km/s, fit statistic {ll:.5f}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 14: the release factor against the singer-absorber hold, on SPARC (round-12 law, u and heat held)',
                                           rows=rows, law_u_kms=u), indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
