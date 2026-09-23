"""Shared plumbing for the suite: paths, the run context, cached data and the constant refits."""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
CODE = RESULTS / 'code'
ROOT = RESULTS.parents[2]
for p in (CODE, ROOT / 'research_work/tools', ROOT / 'research_work/results/cross-prediction-response/code'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

K_SI = 1e6 / 3.0856775814913673e19


class Context:
    """What every test module receives: the tier, a cache folder, a logger and shared data."""
    def __init__(self, tier='quick', cache_dir=None, verbose=True):
        self.tier = tier
        self.full = tier == 'full'
        self.cache = Path(cache_dir) if cache_dir else HERE / 'cache'
        self.cache.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose
        self.shared = {}
        self.t0 = time.monotonic()

    def log(self, msg):
        if self.verbose:
            print(f'  [{time.monotonic() - self.t0:6.0f} s] {msg}', flush=True)

    def get(self, key, fn):
        if key not in self.shared:
            self.shared[key] = fn()
        return self.shared[key]

    def sparc(self):
        import run as R
        return self.get('sparc', R.load_sparc)

    def xcop(self):
        import run_v3 as R3
        return self.get('xcop', R3.prepare_clusters)


def refit_constants(law, ctx):
    """Refit the constants named in law['refit'] on their home data, with the amendments in place:
    'a' on the 149 SPARC galaxies with g_d held fixed (the fit statistic of round 3, the mean squared
    log residual), 'u' on the 12 X-COP clusters. Alternated three times when both are named."""
    from scipy.optimize import minimize_scalar
    import run as R, run_v3 as R3
    from law_config import with_constants
    want = set(law['refit'])
    if not want:
        return law
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    gd = lam * a
    log = []
    for _ in range(3 if want >= {'a', 'u'} else 1):
        if 'a' in want:
            gals = ctx.sparc()
            r = minimize_scalar(lambda la: R.sparc_score(gals, lambda g: R.galaxy_g(g, 10 ** la, u, gd / 10 ** la))[1],
                                bounds=(np.log10(a) - 0.3, np.log10(a) + 0.3), method='bounded', options=dict(xatol=1e-5))
            a = 10 ** r.x; lam = gd / a
            log.append(dict(step='a on SPARC (g_d held)', a_SI=a * K_SI, lam=lam))
        if 'u' in want:
            u = R3.fit_u3(ctx.xcop(), a, lam)
            log.append(dict(step='u on X-COP', u_kms=u))
    out = with_constants(law, a_code=a, lam=lam, u_kms=u)
    out['refit_log'] = log
    return out
