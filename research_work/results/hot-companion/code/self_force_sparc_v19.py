"""Round 19: how large a constant self-pull can the galaxies carry?

code/full_wave_v19.py finds that in a stream that absorbs counter-moving waves every emitter is pulled against the
stream (toward the source of the stream, a galaxy's centre) by a fraction of its emitted power over the wave speed:
3-22% of P/c for absorption strengths of 0.5-20 per wavelength. Per kilogram that is a constant acceleration
eta a (a = 2 l/u; with the wave speed c <= u/2 and a few per cent of P/c, eta ~ 0.03-0.1), present where the companion
is released. This asks what SPARC allows:
    g = g_N + exp(-|g_N|/g_d) [ sqrt(a (|g_N| + S)) + eta a ]
with a refitted at each eta (g_d and u held), on all 149 galaxies and on the held-out split.

    python code/self_force_sparc_v19.py --output run-full-wave-v19/self_force_sparc_v19.json
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); t0 = time.monotonic()
    import common as C, law as L, run as R
    from law_config import load_law
    law = load_law('round12'); ctx = C.Context(tier='quick', verbose=False); C.apply_distances(law, ctx)
    gals = ctx.sparc(); a0, lam, u = law['a_code'], law['lam'], law['u_kms']; gd = lam * a0
    def g_of(g, a, eta):
        S = 0.0 if g['sigb'] == 0 else L.scalar_sum(g['r'], g['sfine'], g['dmb'], L.heat_weight(g['sigb'], u))
        f = np.exp(-np.abs(g['gN']) / gd)
        return g['gN'] + f * (np.sqrt(a * (np.abs(g['gN']) + S)) + eta * a)
    rows = []
    for eta in (0.0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2):
        r = minimize_scalar(lambda la: R.sparc_score(gals, lambda g: g_of(g, 10 ** la, eta))[1],
                            bounds=(np.log10(a0) - 0.3, np.log10(a0) + 0.3), method='bounded', options=dict(xatol=1e-5))
        a = 10 ** r.x
        rm, st = R.sparc_score(gals, lambda g: g_of(g, a, eta))
        out = {s: R.sparc_score(gals, lambda g: g_of(g, a, eta), split=s) for s in ('train', 'validation', 'test')}
        # the outer points: log(v_obs^2 / v_pred^2) at each galaxy's last point, against radius
        last = np.array([np.log(g['v'][-1] ** 2 / g['r'][-1] / g_of(g, a, eta)[-1]) for g in gals])
        rows.append(dict(eta=eta, a_SI=a * C.K_SI, rms_kms=rm, statistic=st, splits={k: dict(rms_kms=v[0], statistic=v[1]) for k, v in out.items()},
                         outer_mean=float(last.mean()), outer_se=float(last.std(ddof=1) / np.sqrt(len(last)))))
        print(f"eta {eta:4.2f}: a = {a * C.K_SI:.3e}, SPARC {rm:.2f} km/s (statistic {st:.4f}); held-out validation {out['validation'][0]:.2f}, "
              f"test {out['test'][0]:.2f}; last points mean log(v_obs^2/v_pred^2) {rows[-1]['outer_mean']:+.3f} +- {rows[-1]['outer_se']:.3f}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 19: a constant self-pull on SPARC', law=law['name'], rows=rows,
                                           seconds=time.monotonic() - t0), indent=1) + '\n')


if __name__ == '__main__':
    main()
