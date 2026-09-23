"""Round 9: what the adopted release length allows, and what the Galaxy's hold does to the dwarfs.

Three scans, with the round-3 constants (run-v3) and the suite's own test code:
1. Cassini's Q2 and the wide-binary boosts (3,000 / 7,000 / 20,000 AU) for release lengths
   L = 0 ... 2 x 10^6 AU and holds c = 1, 0.3, 0.1, 0 (t_precision.sun_in_galaxy).
2. The ten Milky Way dwarfs' predicted speed spreads for c = 1 ... 0 (t_dwarfs).
3. Which part of the Galaxy's influence slows the dwarfs: its pull g_N, its heat S, or both.
4. How much of each dwarf's own companion the release length still holds back at its half-light
   radius (conservative: every direction counted as a scalar flux, Plummer light profile).

The hold c scales the outside galaxy's pull and heat as a separate system (a dwarf, the Sun, a
binary) feels them in its own law; c = 1 is the law as it stands. The release length L makes the
companion leaving an emitter build up as R(r) = 1 - exp(-r/L).

    python code/release_hold_v9.py --output-dir run-release-hold-v9
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'regression'))
import common as C                      # noqa: E402  (sets the paths to the model code)
from law_config import load_law         # noqa: E402
import t_precision as P                 # noqa: E402
import t_dwarfs as TD                   # noqa: E402
import mw_dwarfs_v7 as D                # noqa: E402

LENGTHS_AU = (0, 30000, 50000, 100000, 200000, 500000, 1000000, 2000000)
HOLDS_SUN = (1.0, 0.3, 0.1, 0.0)
HOLDS_DWARFS = (1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.0)
AU_PER_PC = 206264.8


def unreleased_share(b_pc, L_pc):
    """1 - R_eff at r = b for a Plummer sphere of scale b: int rho e^{-d/L}/d^2 dV / int rho/d^2 dV,
    integrated in spherical shells of radius d about the point."""
    d = np.geomspace(1e-6 * b_pc, 200 * b_pc, 4000)
    mu, w = np.polynomial.legendre.leggauss(200)
    r2 = b_pc ** 2 + d[:, None] ** 2 + 2 * b_pc * d[:, None] * mu[None, :]
    shell = 2 * np.pi * ((1 + r2 / b_pc ** 2) ** -2.5 * w[None, :]).sum(1)
    return float(np.trapezoid(shell * np.exp(-d / L_pc), d) / np.trapezoid(shell, d))


def zscore(v, obs, err):
    return (v - obs) / err


def grade(z):
    return 'pass' if abs(z) <= 2 else 'close' if abs(z) <= 3 else 'fail'


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    law = load_law('round3')
    a, gd = law['a_SI'], law['g_d_SI']
    ge, gobs = P.galactic_pull(law)

    # 1. Cassini and wide binaries
    sun = []
    print(f"Galaxy's Newtonian pull at the Sun {ge:.3e} m/s^2 (observed {gobs:.3e})")
    print(f"{'L (AU)':>9} {'L (pc)':>7} {'hold':>5} {'Q2 (1e-27)':>11} {'z':>6} {'b3k':>6} {'b7k':>6} {'b20k':>6}")
    for c in HOLDS_SUN:
        for L in LENGTHS_AU:
            Q2, b = P.sun_in_galaxy(P.MSUN, c * ge, 0.0, 0.0, L, a, gd)
            z = zscore(Q2, 3e-27, 3e-27)
            row = dict(L_au=L, L_pc=L / AU_PER_PC, hold=c, Q2=Q2, Q2_z=z, Q2_grade=grade(z),
                       boost_3000=b[3000], boost_7000=b[7000], boost_20000=b[20000],
                       binaries_grade='pass' if 1.0 - 0.2 <= b[20000] <= 1.5 + 0.2 else 'close' if b[20000] <= 1.5 + 0.3 else 'fail')
            sun.append(row)
            print(f"{L:9.0f} {L / AU_PER_PC:7.2f} {c:5.2f} {Q2 / 1e-27:11.2f} {z:6.2f} {b[3000]:6.3f} {b[7000]:6.3f} {b[20000]:6.3f}", flush=True)

    # 2 and 3. The dwarfs
    ctx = C.Context(tier='quick', verbose=False)
    mw = TD.galaxy_profile(law, ctx)
    data = json.loads((C.RESULTS / 'data/mw_dwarfs.json').read_text())
    variants = dict(full=(1.0, 1.0), pull_only=(1.0, 0.0), heat_only=(0.0, 1.0), neither=(0.0, 0.0))
    dwarfs, chi_hold, pass_hold, chi_var = [], {c: 0.0 for c in HOLDS_DWARFS}, {c: 0 for c in HOLDS_DWARFS}, {k: 0.0 for k in variants}
    for d in data['dwarfs']:
        env = D.galaxy_env(d['D_gc_kpc'], mw, law)
        b = d['r_h_pc'] / 1000.0
        g_int = D.G * 2.0 * d['L_V'] * b / (2 * b * b) ** 1.5
        sig = lambda cg, cs: D.sigma_los(2.0 * d['L_V'], b, dict(gN=cg * env['gN'], S=cs * env['S']), 'ours', law,
                                         3 * d['sigma_obs'] ** 2 / law['u_kms'] ** 2)[0]
        row = dict(name=d['name'], sigma_obs=d['sigma_obs'], sigma_err=d['sigma_err'],
                   galaxy_pull_over_own=env['gN'] / g_int, galaxy_heat_over_own=env['S'] / g_int, by_hold={}, by_part={},
                   unreleased_at_r_h={str(L): unreleased_share(d['r_h_pc'], L) for L in (0.15, 0.5, 1.0, 10.0)})
        for c in HOLDS_DWARFS:
            s = sig(c, c); z = zscore(s, d['sigma_obs'], d['sigma_err'])
            row['by_hold'][str(c)] = dict(sigma=s, z=z, grade=grade(z))
            chi_hold[c] += z * z; pass_hold[c] += grade(z) == 'pass'
        for k, (cg, cs) in variants.items():
            s = sig(cg, cs); z = zscore(s, d['sigma_obs'], d['sigma_err'])
            row['by_part'][k] = dict(sigma=s, z=z, grade=grade(z))
            chi_var[k] += z * z
        dwarfs.append(row)
    print('\ndwarf'.ljust(13) + 'obs'.rjust(11) + ''.join(f'{c:>12}' for c in HOLDS_DWARFS))
    for r in dwarfs:
        print(r['name'].ljust(12) + f"{r['sigma_obs']:6.2f}+-{r['sigma_err']:4.2f}"
              + ''.join(f"{r['by_hold'][str(c)]['sigma']:7.2f} {r['by_hold'][str(c)]['grade'][:4]:4}" for c in HOLDS_DWARFS))
    print('chi2'.ljust(23) + ''.join(f'{chi_hold[c]:12.1f}' for c in HOLDS_DWARFS))
    print('passes'.ljust(23) + ''.join(f'{pass_hold[c]:12d}' for c in HOLDS_DWARFS))
    print('\nthe Galaxy in each dwarf\'s law, chi-squared: ' + ', '.join(f'{k} {v:.1f}' for k, v in chi_var.items()))
    print('\nshare of each dwarf\'s companion still held at its half-light radius, L = 0.15 / 0.5 / 1 / 10 pc:')
    for r in dwarfs:
        print(f"  {r['name']:12s} " + '  '.join(f'{v * 100:5.2f}%' for v in r['unreleased_at_r_h'].values()))

    res = dict(experiment='round 9: release length and the Galaxy\'s hold',
               constants=dict(a_SI=a, g_d_SI=gd, u_kms=law['u_kms']), galaxy_pull_at_sun_SI=ge, observed_pull_at_sun_SI=gobs,
               cassini=dict(measured='(3 +- 3) x 10^-27 s^-2 (Hees et al. 2014)', rows=sun),
               wide_binaries=dict(range='1.0 (Banik et al. 2024) to 1.5 (Chae 2023-24) at 20,000 AU'),
               dwarfs=dict(model='Plummer spheres, M/L_V = 2, the Galaxy from the McMillan 2017 visible components (as in the suite)',
                           rows=dwarfs, chi2_by_hold={str(c): chi_hold[c] for c in HOLDS_DWARFS},
                           passes_by_hold={str(c): pass_hold[c] for c in HOLDS_DWARFS}, chi2_by_part=chi_var,
                           mond_same_stars_chi2=119),
               seconds=time.monotonic() - t0)
    (out / 'release_hold_v9.json').write_text(json.dumps(res, indent=1))
    print(f"\nwrote {out / 'release_hold_v9.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
