"""Round 20: the clusters with the stars' random speeds taken from the law's own gravity (no X-ray input to the heat).

The suite's cluster test (regression/t_clusters.py) computes the stars' velocity dispersion from the Jeans equation in
the MEASURED hydrostatic acceleration, which is the quantity the law is asked to predict (README §10.2, audit §1.6).
Here the stars' dispersion comes from the Jeans equation in the law's OWN acceleration (visible matter only, with the
heat term the stars themselves supply), iterated to a fixed point, exactly as round 3's self_consistent_sigma did, but
on today's inputs: the round-12 law, the static distance law, the deprojected stars. Nothing else changes.

Three questions:
  1. how much the typical miss grows at the adopted constants (u = 169.4 km/s);
  2. which u the self-consistent clusters prefer, and how good they are there;
  3. where the self-consistent speeds differ from the X-ray ones (radius by radius), i.e. where the law's own gravity
     falls short of the measured one.

    python code/xcop_selfconsistent_v20.py --output run-xcop-selfconsistent-v20/xcop_selfconsistent_v20.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                     # noqa: E402
import run_v3 as R3                 # noqa: E402

G = L.G
BLK = 10                            # coarse grid: blocks of 10 fine shells (300 points), as in round 3


def coarse(c):
    """The coarse-grid fields round 3's self-consistent solver needs, rebuilt from the static, deprojected sample."""
    s = c['s']; n = len(s) // BLK * BLK
    c['sc'] = s[:n].reshape(-1, BLK).mean(1)
    c['dms_c'] = c['dms'][:n].reshape(-1, BLK).sum(1)
    c['dmg_c'] = c['dmg'][:n].reshape(-1, BLK).sum(1)
    c['Wc'] = L.shell_weights(c['sc'], c['sc'])
    c['Wk_c'] = L.shell_weights(c['Rk'], c['sc'])
    c['Mb_c'] = np.cumsum(c['dms_c'] + c['dmg_c'])
    c['g_obs_c'] = np.interp(c['sc'], s, c['g_obs_grid'])
    return c


def solve(cls, a, lam, u0, fit_u=True):
    """Damped outer loop over u (the dispersions depend on u), inner fixed point per cluster (round 3's solver)."""
    u = u0
    for _ in range(20):
        for c in cls:
            R3.self_consistent_sigma(c, a, u, lam)
        if not fit_u:
            break
        u_new = R3.fit_u3(cls, a, lam, sigma='sc')
        if abs(u_new / u - 1) < 1e-3:
            u = u_new
            break
        u = np.sqrt(u * u_new)
    for c in cls:
        R3.self_consistent_sigma(c, a, u, lam)
    return u


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    import t_clusters as TC
    law = load_law('round12')
    ctx = C.Context(tier='quick', verbose=False)
    C.apply_distances(law, ctx)
    cls = [coarse(dict(c)) for c in TC.static_clusters(ctx)]
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    out = dict(experiment='round 20: X-COP with the stars\' dispersions from the law\'s own gravity (static distances, deprojected stars)',
               law=dict(name=law['name'], a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=u))

    def score(tag, fn, uu):
        r = R3.resid(cls, fn)
        row = dict(u_kms=float(uu), rms=R3.rms(r), typical_miss_percent=float(100 * (np.exp(R3.rms(r)) - 1)),
                   mean_by_radius=[float(x) for x in r.mean(0)], worst_radius=float(np.max(np.abs(r.mean(0)))),
                   per_cluster_rms={c['name']: float(np.sqrt(np.mean(ri ** 2))) for c, ri in zip(cls, r)})
        out[tag] = row
        print(f"{tag:58s} u {uu:6.1f}  rms {row['rms']:.3f} ({row['typical_miss_percent']:.0f}%)  by radius "
              + ' '.join(f'{x:+.2f}' for x in row['mean_by_radius']), flush=True)
        return row

    score('X-ray speeds, adopted u', lambda c: R3.cluster_M3(c, a, u, lam), u)
    u_x = R3.fit_u3(cls, a, lam)
    score('X-ray speeds, u refitted', lambda c: R3.cluster_M3(c, a, u_x, lam), u_x)
    solve(cls, a, lam, u, fit_u=False)
    score('law\'s own speeds, adopted u', lambda c: R3.cluster_M3(c, a, u, lam, sigma='sc'), u)
    ratio_adopted = {c['name']: c['g_sc_over_obs_c'].copy() for c in cls}
    u_sc = solve(cls, a, lam, u, fit_u=True)
    score('law\'s own speeds, u refitted', lambda c: R3.cluster_M3(c, a, u_sc, lam, sigma='sc'), u_sc)

    # where the law's own gravity differs from the measured one, and the stars' speeds with it (at the refitted u)
    xs = np.array([0.1, 0.2, 0.3, 0.5, 0.7, 1.0])
    g_ratio, sig_ratio = [], []
    for c in cls:
        x = c['sc'] / c['R5']
        g_ratio.append(np.interp(xs, x, c['g_sc_over_obs_c']))
        sig_hse = R3.jeans_sigma2(c['sc'], np.gradient(c['sc']), c['dms_c'], c['g_obs_c'])
        sig_ratio.append(np.sqrt(np.interp(xs, x, c['sig2_star_sc_c']) / np.maximum(np.interp(xs, x, sig_hse), 1e-30)))
    g_ratio, sig_ratio = np.array(g_ratio), np.array(sig_ratio)
    out['radii_R500'] = xs.tolist()
    out['law_gravity_over_measured'] = dict(median=np.median(g_ratio, 0).tolist(), p16=np.percentile(g_ratio, 16, 0).tolist(),
                                            p84=np.percentile(g_ratio, 84, 0).tolist())
    out['star_speed_law_over_xray'] = dict(median=np.median(sig_ratio, 0).tolist(), p16=np.percentile(sig_ratio, 16, 0).tolist(),
                                           p84=np.percentile(sig_ratio, 84, 0).tolist())
    out['iterations'] = {c['name']: int(c['sc_iterations']) for c in cls}
    print('law gravity / measured (median) at ' + ', '.join(f'{x:g}' for x in xs) + ' R500: '
          + ' '.join(f'{v:.2f}' for v in out['law_gravity_over_measured']['median']))
    print('stars\' speed, law / X-ray (median): ' + ' '.join(f'{v:.2f}' for v in out['star_speed_law_over_xray']['median']))
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
