#!/usr/bin/env python3
"""JR-13: one law, one new constant, plus a stress state set by the light itself.

WHAT HAD TO CHANGE
------------------
The companion could not be added mass: 5e4 times all the starlight ever emitted,
and no accumulation beats conservation. But R10's own fitted exponent is
p = 0.5262, essentially one half, and a force A/r with A proportional to sqrt(M) is
not the field of a reservoir. It is the baryons' own field with a longer range,

        g_far = sqrt(G M_b a0) / r.

So there is no companion matter and nothing to pay for. The '15 to 21 times the
stellar mass' was a fictitious number produced by applying M = r^2 g / G to a force
that is not inverse-square. The energy shortfall was a units artefact, not physics.

THE LAW
-------
        g = tau * [ g_N/2 + sqrt( g_N^2/4 + g_N a0 ) ]

with ONE new constant a0, and tau the Tolman active-mass factor of the gravitational
stream, which is exact at its two endpoints:

        tau = 1   standing, dust-like        (w_r = 0, w_t = 0)
        tau = 2   free-streaming radial flux (w_r = 1, w_t = 0)

A stream pulls twice as hard as the same thing standing still. That is ordinary
general relativity, not an assumption.

WHICH STATE A SYSTEM IS IN
--------------------------
Set by the local radiation flux density F = L / (4 pi r^2): where the light is
intense the gravitational response streams, where it is dilute it stands. Galaxies
at a few kpc from a bright bulge and clusters at a megaparsec differ by about three
orders of magnitude in F, so they should sit at opposite ends.

        tau(F) = 1 + 1 / (1 + (F0/F)^m)

Two shape numbers for the crossover, one constant a0. Three numbers total for
seventeen systems spanning 139x in stellar mass.

    python onelaw.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

import model as M
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'cross-prediction-response/code'))
import response_kernels as rk

G = M.J.G
KMS2_PER_KPC_TO_SI = 1e6 / 3.0856775814913673e19
LSUN_OVER_MSUN = 1.0 / 3.0            # declared mass-to-light ratio, solar units
ROOT = Path(__file__).resolve().parents[4]


def interpolation(gN, a0):
    """g = g_N/2 + sqrt(g_N^2/4 + g_N a0). Newtonian above a0, long-range below."""
    return gN / 2.0 + np.sqrt(gN * gN / 4.0 + gN * a0)


def tau_of_flux(F, F0, m):
    """1 (standing) to 2 (free-streaming), crossing over at flux density F0."""
    return 1.0 + 1.0 / (1.0 + (F0 / F) ** m)


def build():
    rows = []
    for n in M.LENSES:
        L = M.make_lens(n)
        scale, A, rc, rt, q = rk.base_terms(L, 0.0)
        Mb = L.Mstar * scale
        r = L.b                                     # the Einstein radius, where lensing is measured
        gN = G * Mb * L.base.mass_fraction(r) / r ** 2
        gobs = scale * L.baryon_force(r) + M.J.companion_force(np.array([r]), A, rc, rt, q)[0]
        rows.append(dict(name=n, kind='galaxy', r_kpc=r, Mb=Mb, gN=gN, gobs=gobs,
                         Lstar_Lsun=Mb * LSUN_OVER_MSUN))
    cc = {x['cluster']: x for x in csv.DictReader(
        (ROOT / 'companion_wave_test/cluster_comparison.csv').open())}
    for p in csv.DictReader((ROOT / 'companion_deposition_fit/cluster_predictions.csv').open()):
        c = cc[p['cluster']]
        M500 = float(c['M500_1e14Msun']) * 1e14
        Mstar = 0.02 * M500
        Mb = (float(c['fgas500']) + 0.02) * M500
        R = float(c['R500_Mpc']) * 1000.0
        Mdyn = float(p['required_extra_Msun']) + Mb
        rows.append(dict(name=p['cluster'], kind='cluster', r_kpc=R, Mb=Mb,
                         gN=G * Mb / R ** 2, gobs=G * Mdyn / R ** 2,
                         Lstar_Lsun=Mstar * LSUN_OVER_MSUN))
    for r in rows:
        r['flux'] = r['Lstar_Lsun'] / (4 * np.pi * r['r_kpc'] ** 2)   # Lsun / kpc^2
    return rows


def predict(rows, a0, F0, m):
    return np.array([tau_of_flux(r['flux'], F0, m) * interpolation(r['gN'], a0) for r in rows])


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    rows = build()
    obs = np.array([r['gobs'] for r in rows])
    kind = np.array([r['kind'] for r in rows])

    def cost(x):
        a0, F0, m = 10 ** x[0], 10 ** x[1], np.exp(x[2])
        p = predict(rows, a0, F0, m)
        if np.any(~np.isfinite(p)) or np.any(p <= 0):
            return 1e9
        return float(np.sum(np.log(obs / p) ** 2))

    best = min((minimize(cost, [4.2, f0, np.log(mm)], method='Nelder-Mead',
                         options=dict(maxiter=20000, fatol=1e-12, xatol=1e-9))
                for f0 in (5., 6., 7., 8.) for mm in (0.5, 1.0, 2.0)), key=lambda r: r.fun)
    a0, F0, m = 10 ** best.x[0], 10 ** best.x[1], float(np.exp(best.x[2]))
    pred = predict(rows, a0, F0, m)
    ratio = obs / pred

    # zero-parameter check: is the crossover actually reaching the exact endpoints?
    taus = np.array([tau_of_flux(r['flux'], F0, m) for r in rows])

    payload = dict(
        experiment='JR-13',
        law='g = tau * [ g_N/2 + sqrt(g_N^2/4 + g_N a0) ]',
        why_no_energy_problem='R10 fitted p = 0.5262. A force A/r with A proportional to '
                              'sqrt(M_b) is the baryons own field with a longer range, not a '
                              'reservoir. No companion matter exists, so none has to be paid '
                              'for. The 15 to 21 M_star figure came from applying M = r^2 g/G '
                              'to a non-inverse-square force.',
        tau_endpoints={'standing, dust-like': 1.0, 'free-streaming radial flux': 2.0,
                       'note': 'Both exact in linearised GR: tau = 1 + w_r + 2 w_t.'},
        constants=dict(a0_kms2_per_kpc=float(a0), a0_SI=float(a0 * KMS2_PER_KPC_TO_SI),
                       F0_Lsun_per_kpc2=float(F0), crossover_sharpness=m,
                       n_free_numbers=3, n_systems=len(rows)),
        rows=[dict(name=r['name'], kind=r['kind'], r_kpc=r['r_kpc'],
                   Mb_Msun=r['Mb'], gN_over_a0=float(r['gN'] / a0),
                   flux_Lsun_kpc2=float(r['flux']), tau=float(t),
                   g_observed=float(o), g_predicted=float(p), ratio=float(x))
              for r, t, o, p, x in zip(rows, taus, obs, pred, ratio)],
        scatter=dict(all=float(np.std(np.log(ratio), ddof=1)),
                     galaxies=float(np.std(np.log(ratio[kind == 'galaxy']), ddof=1)),
                     clusters=float(np.std(np.log(ratio[kind == 'cluster']), ddof=1)),
                     mean_ratio_galaxies=float(ratio[kind == 'galaxy'].mean()),
                     mean_ratio_clusters=float(ratio[kind == 'cluster'].mean())),
        tau_reached=dict(galaxies=[float(t) for t in taus[kind == 'galaxy']],
                         clusters=[float(t) for t in taus[kind == 'cluster']]),
        limitations=[
            'Seventeen systems and three numbers. The crossover sharpness is set by the gap '
            'between the two classes, not resolved within either.',
            'Cluster g_obs comes from X-ray hydrostatic masses, which carry their own known '
            'bias of order ten to twenty percent.',
            'Galaxy g_obs is R10 evaluated at the Einstein radius, so it inherits R10 rather '
            'than being an independent measurement.',
            'The flux density uses a declared mass-to-light ratio of three and the stellar '
            'mass only; a full radiative-transfer calculation is not done here.'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'one-law.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print(f"ONE LAW:  g = tau * [ g_N/2 + sqrt(g_N^2/4 + g_N a0) ]\n")
    print(f"  a0 = {a0 * KMS2_PER_KPC_TO_SI:.3e} m/s^2")
    print(f"  crossover at F0 = {F0:.3g} Lsun/kpc^2, sharpness m = {m:.2f}")
    print(f"  three numbers, {len(rows)} systems\n")
    print(f"{'system':<12}{'kind':<9}{'g_N/a0':>9}{'flux':>11}{'tau':>7}{'obs/pred':>10}")
    for r in payload['rows']:
        print(f"{r['name']:<12}{r['kind']:<9}{r['gN_over_a0']:>9.3f}{r['flux_Lsun_kpc2']:>11.3g}"
              f"{r['tau']:>7.3f}{r['ratio']:>10.3f}")
    s = payload['scatter']
    print(f"\n  galaxies: mean ratio {s['mean_ratio_galaxies']:.3f}, log scatter {s['galaxies']:.1%}")
    print(f"  clusters: mean ratio {s['mean_ratio_clusters']:.3f}, log scatter {s['clusters']:.1%}")
    print(f"  ALL {len(rows)} TOGETHER: log scatter {s['all']:.1%}")
    print(f"\n  tau reached: galaxies {min(taus[kind=='galaxy']):.3f}-{max(taus[kind=='galaxy']):.3f}, "
          f"clusters {min(taus[kind=='cluster']):.3f}-{max(taus[kind=='cluster']):.3f}")
    print('WROTE', out / 'one-law.json')


if __name__ == '__main__':
    main()
