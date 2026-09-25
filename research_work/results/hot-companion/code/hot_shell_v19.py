"""Round 19, step 2b: the hot-shell benchmark on the clusters and galaxies the law was fitted on.

The law's heat term S counts the hot matter's extra glow at a receiver from every shell, inside and outside it
(code/law.py: S = G sum k dm <1/d^2>_shell). A medium whose wave only travels outward (round 17's rule, round 18's
absorbing stream) does not let a receiver hear hot matter farther out than itself. For round sources this is the only
difference between the two (README 29.3: a stream that absorbs inward waves hears exactly half of every inner shell,
in the scalar sum and in the net flux alike, which only rescales a). This script asks what the X-COP clusters and the
SPARC galaxies say, first with every constant held at the adopted (round-12) values, then with u refitted on X-COP:

  two_way          the law (every shell)
  one_way          inner shells only, each with the law's weight <1/d^2> (artanh(x)/x)
  one_way_vector   inner shells only, each as if at the centre (the net outward flux; Gauss)

    python code/hot_shell_v19.py --output run-hot-shell-v19/hot_shell_v19.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))

GEOMETRIES = ('two_way', 'one_way', 'one_way_vector')


def outer_share(cls, law):
    """The fraction of the law's two-way S that comes from shells outside the receiver, at each of the six radii."""
    import law as L
    rows = []
    for c in cls:
        x = c['s'][None, :] / c['Rk'][:, None]
        w_in = np.where(x < 1, L.w_inside(x), 0.0); w_out = np.where(x >= 1, L.w_outside(x), 0.0)
        src = L.heat_weight(np.sqrt(c['sig2_star_hse']), law['u_kms']) * c['dms']
        Si, So = w_in @ src, w_out @ src
        rows.append(So / (Si + So))
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    import common as C, law as L, run_v3 as R3, run as R, t_clusters as TC
    from law_config import load_law
    law = load_law('round12')
    ctx = C.Context(tier='quick', verbose=False)
    C.apply_distances(law, ctx)
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    out = dict(experiment='round 19: the hot-shell benchmark (how the hot matter is heard) on X-COP and SPARC',
               law=dict(name=law['name'], a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=u), geometries={})
    L.HOT_GEOMETRY = 'two_way'
    cls0 = TC.static_clusters(ctx)
    share = outer_share(cls0, law)
    out['outer_share_of_S'] = dict(radii_over_R500=[0.1, 0.2, 0.3, 0.5, 0.7, 1.0] if len(cls0[0]['Rk']) == 6 else None,
                                   per_cluster={c['name']: s.tolist() for c, s in zip(cls0, share)},
                                   median=np.median(share, 0).tolist())
    print('share of the law\'s S from shells outside the receiver (median over clusters, at the six radii):',
          ' '.join(f'{v:.2f}' for v in np.median(share, 0)), flush=True)
    gals = ctx.sparc()
    for geo in GEOMETRIES:
        L.HOT_GEOMETRY = geo
        ctx.shared.pop('xcop_static', None)
        cls = TC.static_clusters(ctx)                     # rebuilt: their shell weights follow the geometry
        res = R3.resid(cls, lambda c: R3.cluster_M3(c, a, u, lam))
        u_fit = R3.fit_u3(cls, a, lam)
        res_fit = R3.resid(cls, lambda c: R3.cluster_M3(c, a, u_fit, lam))
        gscore = R.sparc_score(gals, lambda g: R.galaxy_g(g, a, u, lam))
        gscore_fit = R.sparc_score(gals, lambda g: R.galaxy_g(g, a, u_fit, lam))
        row = dict(xcop_rms=R3.rms(res), xcop_mean_by_radius=res.mean(0).tolist(), xcop_per_cluster={c['name']: r.tolist() for c, r in zip(cls, res)},
                   u_refit_kms=float(u_fit), xcop_rms_refit=R3.rms(res_fit), xcop_mean_by_radius_refit=res_fit.mean(0).tolist(),
                   sparc_rms_kms=gscore[0], sparc_statistic=gscore[1], sparc_rms_kms_at_refit_u=gscore_fit[0], sparc_statistic_at_refit_u=gscore_fit[1])
        out['geometries'][geo] = row
        print(f"[{time.monotonic() - t0:5.0f} s] {geo:15s} X-COP rms {row['xcop_rms']:.3f} (worst radius {np.max(np.abs(res.mean(0))):.3f}); "
              f"u refit {u_fit:.1f} km/s -> rms {row['xcop_rms_refit']:.3f}; SPARC {gscore[0]:.2f} km/s (stat {gscore[1]:.4f}), "
              f"at the refitted u {gscore_fit[0]:.2f} km/s", flush=True)
    L.HOT_GEOMETRY = 'two_way'
    print('the absorbing stream at finite strength (straight rays; hot glow scalar, cold net flux):', flush=True)
    out['absorbing_stream_scan'] = stream_scan()
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('wrote', args.output)



# ---------------------------------------------------------------------------------------------------------------
# The absorbing stream at finite strength: how far inward must the hot glow travel for the clusters to be happy?
def stream_weights(R, s, kappa):
    """The absorbing stream's scalar and net-flux shell weights (code/law.py: stream_weights)."""
    import law as L
    return L.stream_weights(R, s, kappa)


def stream_scan(kappas_per_Mpc=(0.0, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)):
    """X-COP with the stream's hearing (hot glow: scalar; cold: net flux; both with the absorption), a held at its
    SPARC value (galaxies are far smaller than the absorption lengths scanned, so they hear two-way), u refitted."""
    import common as C, law as L, run_v3 as R3, t_clusters as TC
    from law_config import load_law
    from scipy.optimize import minimize_scalar
    law = load_law('round12'); ctx = C.Context(tier='quick', verbose=False); C.apply_distances(law, ctx)
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    L.HOT_GEOMETRY = 'two_way'; cls = TC.static_clusters(ctx)
    G = L.G
    rows = []
    for kap in kappas_per_Mpc:
        for c in cls:
            ws, wv = stream_weights(c['Rk'], c['s'], kap / 1000.0)
            c['_ws'], c['_wv'] = ws, wv
        def pred(c, uu):
            gN = G * c['Mb'] / c['Rk'] ** 2
            cold = G * (c['_wv'] @ (c['dmg'] + c['dms'])) / c['Rk'] ** 2
            S = G * (c['_ws'] @ (L.heat_weight(np.sqrt(c['sig2_star_hse']), uu) * c['dms'])) / c['Rk'] ** 2
            f = L.released(gN, a, lam)
            return (gN + f * np.sqrt(a * np.maximum(cold + S, 0.0))) * c['Rk'] ** 2 / G
        res = lambda uu: np.array([np.log(c['Mh'] / pred(c, uu)) for c in cls])
        r0 = res(u)
        opt = minimize_scalar(lambda lu: np.mean(res(10 ** lu) ** 2), bounds=(1.3, 4.0), method='bounded')
        uf = 10 ** opt.x; rf = res(uf)
        cold_frac = float(np.median([(c['_wv'] @ (c['dmg'] + c['dms'])) / (c['Mb']) for c in cls]))
        rows.append(dict(kappa_per_Mpc=kap, absorption_length_kpc=(1000.0 / kap if kap else None),
                         rms_fixed=R3.rms(r0), mean_by_radius_fixed=r0.mean(0).tolist(),
                         u_refit_kms=float(uf), rms_refit=R3.rms(rf), mean_by_radius_refit=rf.mean(0).tolist(),
                         cold_heard_fraction_median=cold_frac))
        print(f"kappa {kap:6.1f}/Mpc (absorption length {'inf' if not kap else f'{1000 / kap:6.0f} kpc'}): rms {R3.rms(r0):.3f} at u = {u:.0f}; "
              f"u refit {uf:6.1f} -> rms {R3.rms(rf):.3f}; mean by radius " + ' '.join(f'{x:+.2f}' for x in rf.mean(0)) +
              f"; cold heard {cold_frac:.2f}", flush=True)
    return rows


if __name__ == '__main__':
    main()
