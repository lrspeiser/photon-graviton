"""Clusters: the 12 X-COP hydrostatic mass profiles (the data u was fitted on).

Since round 11 the profiles are converted from the release's flat-LCDM units (H0 = 70, Omega_m = 0.3) to
the project's static distance law at each cluster's redshift (code/xcop_static_v11.py): radii and
hydrostatic masses x D_A ratio, gas x D_L D_A^1.5 ratio, stars x D_L^2 ratio. The release's stellar
profiles are masses inside a projected radius (Ghizzardi et al. 2021, Sect. 4.1); since round 11 they are
deprojected first (xcop_static_v11.deproject_xcop). The 'round11' constants were refitted this way; the
round-3/9 constants were fitted in the release's units with the projected stars."""
from __future__ import annotations
import itertools
import numpy as np
from checks import make, at_most

GROUP = 'clusters'


def static_clusters(ctx):
    def build():
        import json
        import xcop_static_v11 as XS
        R = XS.RUN                        # the hot-companion run.py (a bare 'import run' can resolve to another project's)
        raw = json.loads((R.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
        dep, _ = XS.deproject_xcop(ctx.xcop(), raw)            # round 11: the release's stellar profiles are projected
        return [XS.to_static(c, raw[c['name']]['header']['z']) for c in dep]
    return ctx.get('xcop_static', build)


def run(law, ctx):
    import run_v3 as R3
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    cls = static_clusters(ctx)
    ctx.log('X-COP masses')
    res = R3.resid(cls, lambda c: R3.cluster_M3(c, a, u, lam))
    rms = R3.rms(res)
    out = [make('clusters.xcop_rms', GROUP, 'X-COP: typical mass miss, rms ln(M_hydrostatic/M_predicted) (static distances)', rms,
                crit=at_most(rms, 0.30, 0.45), target='<= 0.30 (hydrostatic masses carry 10-20% bias and error); MOND 1.06, NFW fit 0.10',
                refs='Ettori et al. 2019; Ghizzardi et al. 2021 (X-COP)')]
    mb = res.mean(0)
    worst = float(np.max(np.abs(mb)))
    out.append(make('clusters.xcop_radial_trend', GROUP, 'X-COP: worst mean miss at any of the six radii (0.1-1 R500)', worst,
                    crit=at_most(worst, 0.25, 0.35), target='<= 0.25 (the usual hydrostatic bias at R500)',
                    detail=dict(mean_by_radius=mb.tolist())))
    u_best = R3.fit_u3(cls, a, lam)
    out.append(make('clusters.u_preferred', GROUP, 'the companion speed X-COP prefers for these a and g_d', u_best,
                    unit='km/s', target=f"law uses {u:.1f}", crit=('info', abs(np.log(u_best / u)))))
    if ctx.full:
        ctx.log('X-COP cross-validation (924 half splits)')
        cv = []
        for tr in itertools.combinations(range(len(cls)), len(cls) // 2):
            te = [i for i in range(len(cls)) if i not in tr]
            uu = R3.fit_u3([cls[i] for i in tr], a, lam)
            cv.append(R3.rms(R3.resid([cls[i] for i in te], lambda c: R3.cluster_M3(c, a, uu, lam))))
        m = float(np.median(cv))
        out.append(make('clusters.xcop_crossval', GROUP, 'X-COP: held-out typical miss, u fitted on the other half', m,
                        crit=at_most(m, 0.30, 0.45), target='<= 0.30', detail=dict(p5=float(np.percentile(cv, 5)), p95=float(np.percentile(cv, 95)))))
    return out
