"""Round 11: why the ordered companion adds up the way Newton's pulls do.

A. Energy bookkeeping. In a steady state the companion's energy flux J obeys div J = l rho: every watt
   emitted flows out (nothing is lost on the way). If, at each point, the companion is one stream
   moving at speed u (the condition for an undiluted pull, from the local form of round 10), its
   energy density is |J|/u. If the flow has no whirlpools (curl J = 0), J is unique:
       J = (l / 4 pi G) (-g_N),
   so the stream's intensity is |g_N| and the pull is proportional to sqrt|g_N|: round 3's rule for
   ordered matter. The alternative "one stream carrying everything" (round 10's rule 2, intensity S_N)
   is not a consistent flow: through a closed surface it would carry out more power than the matter
   inside emits. The ratio is computed here for spheres (exact, by symmetry) and for the Milky Way
   model (a real disk galaxy, by integrating over spheres).
   Independent waves (rule 1) do conserve energy (their net flux is J), but their pull is diluted.
   So round 3 is the one rule that is both lossless and undiluted.

B. How precisely do the data demand it? A three-parameter family that contains all three rules:
       I = w |g_N| + (1 - w) S_N + h S + (1 - h) |g_hot|        (intensity)
       D = |g_N + g_hot| / (S_N + S)                             (directedness of free streaming)
       extra = release x sqrt(a I) x D^gamma,  along the net flow
   Round 3 is (w, h, gamma) = (1, 1, 0); independent waves (0, 1, 1); one stream of everything
   (0, 1, 0); "the heat cancels too, like Newton's arrows" (1, 0, 0). Scans through round 3 along
   each parameter, with a and g_d refitted on SPARC and u on X-COP at every point; the preferred
   value's spread from bootstrap resampling of galaxies and of clusters. X-COP is used as the round-11
   suite uses it: in the project's static distances, with the stellar profiles deprojected
   (xcop_static_v11); the scans start from the adopted (round-11) constants.

    python code/companion_flow_v11.py --output-dir run-companion-flow-v11
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import run as RUN                                 # noqa: E402

G = L.G


# --------------------------------------------------------------------------- A. energy bookkeeping
def sphere_ratio(profile, x):
    """S_N / |g_N| at radius x (units of the scale) for a spherical density profile rho(r): the power a
    'one stream of everything' would carry out through the sphere of radius x, over the power emitted inside."""
    r = np.geomspace(1e-4, 200.0, 20000); dr = np.gradient(r)
    dm = 4 * np.pi * r ** 2 * profile(r) * dr
    W = L.shell_weights(np.atleast_1d(x), r)
    SN = (W @ dm) / x ** 2
    gN = dm[r < x].sum() / x ** 2
    return float(SN[0] / gN)


def milky_way_power_ratio(radii=(4.0, 8.2, 12.0, 20.0, 30.0), n_theta=24):
    """For the Milky Way model (McMillan 2017 visible matter): power carried out through spheres by a
    stream of intensity S_N along the net pull, over the power emitted inside (= the flux of g_N)."""
    import common as C
    import t_milky_way as TMW
    import milky_way_v7 as MW
    import first_principles_v10 as FP
    from scipy.interpolate import RegularGridInterpolator
    ctx = C.Context(tier='quick', verbose=False)
    comps, grid, F = TMW.setup(ctx)
    fac = MW.models(comps)['M17']
    rho = lambda R, z: sum(f * comps[n].rho(R, z) for n, f in fac.items())
    gR = sum(f * F[n][0] for n, f in fac.items()); gz = sum(f * F[n][1] for n, f in fac.items())
    iR = RegularGridInterpolator((grid.R, grid.z), gR); iz = RegularGridInterpolator((grid.R, grid.z), gz)
    rows = []
    mu = np.cos(np.linspace(0, np.pi / 2, n_theta + 1)[:-1] + np.pi / (4 * n_theta))   # upper hemisphere, midpoints in angle
    th = np.arccos(mu); w = np.sin(th) * (np.pi / 2 / n_theta)
    for r0 in radii:
        R = r0 * np.sin(th); z = r0 * np.cos(th)
        gr_ = iR(np.c_[R, z]); gz_ = iz(np.c_[R, z])
        g = np.hypot(gr_, gz_)
        gn_out = -(gr_ * np.sin(th) + gz_ * np.cos(th))                 # inward component (positive = inward)
        SN = G * np.array([FP.mw_scalar_sum(rho, Ri, zi) for Ri, zi in zip(R, z)])
        emitted = np.sum(gn_out * w)                                     # flux of g_N through the sphere / (2 pi r^2), = 4 pi G M(<r) / (2 pi r^2) / 2
        carried = np.sum(SN * gn_out / np.maximum(g, 1e-30) * w)         # a stream of intensity S_N along the net pull
        g0 = float(abs(iR([[r0, 0.0]])[0]))
        rows.append(dict(r_kpc=r0, carried_over_emitted=float(carried / emitted),
                         SN_over_gN_midplane=float(G * FP.mw_scalar_sum(rho, r0, 0.0) / g0)))
    return rows


# --------------------------------------------------------------------------- B. the family, on SPARC and X-COP
def galaxy_family(g, a, u, lam, w, h, gamma):
    k = L.heat_weight(g['sigb'], u) if g['sigb'] > 0 else 0.0
    S = L.scalar_sum(g['r'], g['sfine'], g['dmb'], k * np.ones_like(g['dmb'])) if g['sigb'] > 0 else 0.0
    ghot = k * g['vb2'] / g['r']
    I = w * g['gN'] + (1 - w) * g['SN'] + h * S + (1 - h) * ghot
    D = np.clip((g['gN'] + ghot) / (g['SN'] + S), 0, 1)
    return g['gN'] + L.released(g['gN'], a, lam) * np.sqrt(a * I) * D ** gamma


def cluster_parts(c, u):
    k = L.heat_weight(np.sqrt(c['sig2_star_hse']), u)
    gN = G * c['Mb'] / c['Rk'] ** 2
    S = G * (c['W'] @ (k * c['dms'])) / c['Rk'] ** 2
    SN = G * (c['W'] @ (c['dms'] + c['dmg'])) / c['Rk'] ** 2
    ghot = G * np.array([np.sum((k * c['dms'])[c['s'] < R]) for R in c['Rk']]) / c['Rk'] ** 2
    return gN, S, SN, ghot


def cluster_family(c, a, u, lam, w, h, gamma):
    gN, S, SN, ghot = cluster_parts(c, u)
    I = w * gN + (1 - w) * SN + h * S + (1 - h) * ghot
    D = np.clip((gN + ghot) / (SN + S), 0, 1)
    return (gN + L.released(gN, a, lam) * np.sqrt(a * I) * D ** gamma) * c['Rk'] ** 2 / G


def fit_point(gals, cls, R3, law, w, h, gamma, a0, lam0):
    u = law['u_kms']
    cost = lambda p: RUN.sparc_score(gals, lambda g: galaxy_family(g, 10 ** p[0], u, 10 ** p[1], w, h, gamma))[1]
    best = min((minimize(cost, [np.log10(a0), l0], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-8))
                for l0 in (np.log10(lam0), 0.3, 1.2)), key=lambda r: r.fun)
    a, lam = 10 ** best.x[0], 10 ** best.x[1]
    rms = RUN.sparc_score(gals, lambda g: galaxy_family(g, a, u, lam, w, h, gamma))[0]
    r = minimize_scalar(lambda lu: np.mean(R3.resid(cls, lambda c: cluster_family(c, a, 10 ** lu, lam, w, h, gamma)) ** 2),
                        bounds=(0.5, 3.5), method='bounded')
    uu = 10 ** r.x
    res = R3.resid(cls, lambda c: cluster_family(c, a, uu, lam, w, h, gamma))
    return dict(w=w, h=h, gamma=gamma, a_SI=float(a * L.KMS2_PER_KPC), g_d_SI=float(a * lam * L.KMS2_PER_KPC), lam=float(lam),
                sparc_rms_kms=float(rms), sparc_msq=float(best.fun), u_kms=float(uu), xcop_rms=float(np.sqrt(np.mean(res ** 2))),
                per_galaxy_msq=[float(np.mean(np.log(g['v'] ** 2 / g['r'] / galaxy_family(g, a, u, lam, w, h, gamma)) ** 2)) for g in gals],
                per_cluster_resid=res.tolist())


def bootstrap_best(values, per_item, rng, B=400):
    """per_item[i][j] = the statistic of item j at scan value i (refitted constants held); resample items."""
    M = np.array(per_item)                      # [n_values, n_items]
    n = M.shape[1]; best = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        best.append(values[int(np.argmin(M[:, idx].mean(1)))])
    return [float(x) for x in np.percentile(best, [16, 50, 84])]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic(); rng = np.random.default_rng(20260926)
    from law_config import load_law
    import common as C
    import run_v3 as R3
    import combination_rules_v10 as CR
    import xcop_static_v11 as XS
    law = load_law('round11')
    res = dict(experiment='round 11: why the ordered companion adds up the way Newton\'s pulls do')

    # A
    prof = dict(uniform=lambda r: (r < 1.0).astype(float), plummer=lambda r: (1 + r ** 2) ** -2.5,
                hernquist=lambda r: 1 / (np.maximum(r, 1e-6) * (1 + r) ** 3))
    A_sph = {name: [dict(x=x, carried_over_emitted=sphere_ratio(f, x)) for x in (0.25, 0.5, 1.0, 2.0, 5.0)] for name, f in prof.items()}
    res['A_spheres'] = A_sph
    print('A. a stream of intensity S_N along the net pull: power carried out through a sphere / power emitted inside it')
    for name, rows in A_sph.items():
        print(f"   {name:9s} " + ', '.join(f"r = {r['x']:g}: {r['carried_over_emitted']:.2f}" for r in rows))
    A_mw = milky_way_power_ratio()
    res['A_milky_way'] = A_mw
    print('   Milky Way model: ' + ', '.join(f"{r['r_kpc']:g} kpc: {r['carried_over_emitted']:.2f}" for r in A_mw), flush=True)

    # B
    ctx = C.Context(tier='quick', verbose=False)
    raw = json.loads((RUN.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    dep, _ = XS.deproject_xcop(ctx.xcop(), raw)
    cls = [XS.to_static(c, raw[c['name']]['header']['z']) for c in dep]
    gals = RUN.load_sparc()
    info = CR.sparc_plain_totals(gals, 0.2, 1.0)
    res['sparc_plain_totals'] = info
    a0, lam0 = law['a_code'], law['lam']
    scans = dict(w=[0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0],
                 gamma=[0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0],
                 h=[0.0, 0.25, 0.5, 0.75, 0.9, 1.0])
    base = dict(w=1.0, h=1.0, gamma=0.0)
    res['B_scans'] = {}
    for par, vals in scans.items():
        rows = []
        for v in vals:
            p = dict(base); p[par] = v
            rows.append(fit_point(gals, cls, R3, law, p['w'], p['h'], p['gamma'], a0, lam0))
            r = rows[-1]
            print(f"B. {par} = {v:5.2f}: SPARC {r['sparc_rms_kms']:.2f} km/s (a {r['a_SI']:.3g}, g_d {r['g_d_SI']:.3g});"
                  f" X-COP rms {r['xcop_rms']:.3f} (u {r['u_kms']:.0f} km/s)", flush=True)
        boot_g = bootstrap_best(vals, [r['per_galaxy_msq'] for r in rows], rng)
        boot_c = bootstrap_best(vals, [np.mean(np.array(r['per_cluster_resid']) ** 2, axis=1) for r in rows], rng)
        for r in rows:
            r.pop('per_galaxy_msq'); r['per_cluster_resid'] = [float(np.sqrt(np.mean(np.square(x)))) for x in r['per_cluster_resid']]
        res['B_scans'][par] = dict(values=vals, rows=rows, sparc_preferred_16_50_84=boot_g, xcop_preferred_16_50_84=boot_c)
        print(f"   {par}: SPARC prefers {boot_g[1]:g} (16-84%: {boot_g[0]:g}-{boot_g[2]:g}); X-COP prefers {boot_c[1]:g} ({boot_c[0]:g}-{boot_c[2]:g})", flush=True)
    res['seconds'] = time.monotonic() - t0
    (out / 'companion_flow_v11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'companion_flow_v11.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
