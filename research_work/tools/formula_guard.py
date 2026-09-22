"""Formula guard: check a candidate gravity law against MOND, Newton and dark matter.

Standing rule (RULES.md, 22 September 2026): no candidate is reported as ours until it
has been checked against all three. This module does the checking; it never fits the
candidate. Every published MOND interpolation function below is labelled with its source.

A candidate supplies
    point_mass(r_kpc, M_sun, sigma_kms) -> g   (km^2/s^2/kpc), an isolated point mass;
    predictions on SPARC: arrays gN (catalogue Newtonian) and g (candidate), per galaxy;
    metadata: number of per-object free parameters for an extra component, and whether
    its extra gravity is an independent substance with its own conserved amount.

Usage:  from formula_guard import guard, MOND_CATALOGUE
"""
from __future__ import annotations
import numpy as np

G = 4.30091727003628e-6                     # kpc (km/s)^2 / Msun
KMS2_PER_KPC = 1e6 / 3.0856775814913673e19  # (km/s)^2/kpc -> m/s^2

# Published MOND interpolation functions, nu(y) with y = g_N / a0.
MOND_CATALOGUE = {
    'simple (Famaey & Binney 2005)': lambda y: 0.5 + np.sqrt(0.25 + 1 / y),
    'standard (Milgrom 1983)': lambda y: np.sqrt(0.5 + np.sqrt(0.25 + 1 / y ** 2)),
    'Bekenstein toy / FB eq. 5 (g_N + sqrt(a0 g_N))': lambda y: 1 + 1 / np.sqrt(y),
    'RAR exponential (McGaugh, Lelli & Schombert 2016)': lambda y: 1 / (1 - np.exp(-np.sqrt(y))),
    'delta-family, delta=2 (Famaey & McGaugh 2012)': lambda y: np.sqrt(1 + 1 / y),
    'n-family, n=3 (Famaey & McGaugh 2012)': lambda y: (0.5 + np.sqrt(0.25 + y ** -3.)) ** (1 / 3),
}


def point_mass_test(point_mass, a0_si=1.2e-10, sigma_kms=0.0):
    """Does the law, for isolated point masses, collapse to one function of g_N/a0?"""
    a = a0_si / KMS2_PER_KPC
    y = np.geomspace(1e-4, 1e3, 120)
    curves = []
    for M in np.geomspace(1e7, 1e15, 17):
        r = np.sqrt(G * M / (y * a))
        g = point_mass(r, M, sigma_kms)
        curves.append(np.log10(g / (y * a)))              # log nu(y) for this mass
    curves = np.array(curves)
    spread = float(np.max(np.ptp(curves, axis=0)))
    out = dict(sigma_kms=sigma_kms, max_spread_dex=spread, collapses=spread < 1e-6)
    if out['collapses']:
        nu = curves[0]
        best = []
        for name, f in MOND_CATALOGUE.items():
            # allow the candidate its own acceleration constant: scan a scale factor
            scales = np.geomspace(0.05, 20, 4001)
            dev = [np.max(np.abs(nu - np.log10(f(y * s)))) for s in scales]
            i = int(np.argmin(dev))
            best.append((dev[i], name, scales[i]))
        best.sort()
        out['closest_published'] = [dict(function=n, max_dev_dex=float(d), a0_scale=float(s))
                                    for d, n, s in best[:3]]
        out['identical_to'] = best[0][1] if best[0][0] < 1e-6 else None
    return out


def local_function_test(gN, g, k=12):
    """How much of log g cannot be written as a smooth function of log g_N alone?

    For every point, fit a straight line to log g against log g_N among its k nearest
    neighbours in log g_N (excluding itself) and measure the point's departure from it.
    A law that is any smooth function of g_N alone scores ~1e-4 dex; a law that also
    depends on something else (radius, size, temperature) scores more.
    """
    x = np.log10(gN); z = np.log10(g)
    o = np.argsort(x); xs, zs = x[o], z[o]; n = len(xs)
    resid = np.empty(n)
    for i in range(n):
        lo = max(0, min(i - k // 2, n - k - 1)); idx = np.r_[lo:i, i + 1:lo + k + 1]
        A = np.vstack([xs[idx], np.ones(len(idx))]).T
        c = np.linalg.lstsq(A, zs[idx], rcond=None)[0]
        resid[i] = zs[i] - (c[0] * xs[i] + c[1])
    return dict(rms_beyond_local_gN_dex=float(np.sqrt(np.mean(resid ** 2))),
                max_abs_dex=float(np.max(np.abs(resid))), residual=resid)


def temperature_test(point_mass, a0_si=1.2e-10, sigmas=(0., 100., 300., 1000.)):
    """At the SAME Newtonian acceleration, does the pull depend on the matter's random speed?"""
    a = a0_si / KMS2_PER_KPC
    y = np.geomspace(1e-3, 1e2, 60); M = 1e12; r = np.sqrt(G * M / (y * a))
    logs = np.array([np.log10(point_mass(r, M, s)) for s in sigmas])
    return dict(sigmas_kms=list(sigmas), max_spread_dex=float(np.max(np.ptp(logs, axis=0))),
                depends_on_temperature=bool(np.max(np.ptp(logs, axis=0)) > 1e-6))


def newton_test(per_galaxy):
    """Is the boost g/g_N constant within every galaxy (Newton with rescaled mass)?"""
    spreads = [float(np.std(np.log10(g / gN))) for gN, g in per_galaxy]
    return dict(median_within_galaxy_boost_spread_dex=float(np.median(spreads)),
                reduces_to_newton=bool(np.max(spreads) < 1e-3))


def guard(name, point_mass=None, sparc=None, per_object_params=0,
          independent_substance=False, a0_si=1.2e-10):
    """Run every check and return a report with a plain-language verdict."""
    rep = dict(name=name)
    if point_mass is not None:
        rep['cold_point_mass'] = point_mass_test(point_mass, a0_si, 0.0)
        rep['temperature'] = temperature_test(point_mass, a0_si)
    if sparc is not None:
        gN = np.concatenate([p[0] for p in sparc]); g = np.concatenate([p[1] for p in sparc])
        lf = local_function_test(gN, g)
        rep['sparc_local_function'] = {k: v for k, v in lf.items() if k != 'residual'}
        rep['newton'] = newton_test(sparc)
    rep['dark_matter'] = dict(per_object_free_parameters=per_object_params,
                              independent_substance=independent_substance,
                              is_dark_matter=bool(per_object_params > 0 or independent_substance))
    # verdict
    v = []
    cp = rep.get('cold_point_mass')
    tp = rep.get('temperature')
    lf = rep.get('sparc_local_function')
    mond_everywhere = (cp is None or cp['collapses']) and (tp is None or not tp['depends_on_temperature']) and \
        (lf is None or lf['rms_beyond_local_gN_dex'] < 0.005)
    if mond_everywhere:
        v.append('MOND-CLASS: a function of the local Newtonian acceleration everywhere tested')
    else:
        if cp is not None and cp['collapses']:
            ident = cp.get('identical_to')
            v.append('cold isolated mass: reduces to a MOND form' +
                     (f' (identical to {ident})' if ident else
                      f" (closest published: {cp['closest_published'][0]['function']},"
                      f" {cp['closest_published'][0]['max_dev_dex']:.3f} dex)"))
        if tp is not None and tp['depends_on_temperature']:
            v.append(f"same Newtonian pull, different random speeds 0-1000 km/s: predictions differ by "
                     f"{tp['max_spread_dex']:.3f} dex -- NOT a function of g_N alone")
        if lf is not None:
            v.append(f"SPARC: {lf['rms_beyond_local_gN_dex']:.4f} dex beyond local g_N")
    if rep.get('newton', {}).get('reduces_to_newton'):
        v.append('NEWTON-CLASS')
    if rep['dark_matter']['is_dark_matter']:
        v.append('DARK-MATTER-CLASS')
    rep['verdict'] = '; '.join(v)
    return rep
