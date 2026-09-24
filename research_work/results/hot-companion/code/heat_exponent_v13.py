"""Round 13: what exponent the data give the heat weight, k = 3 (sigma / u)^p (the law has p = 2).

In the heat-dominated regime the extra pull goes as sqrt(k), so p = 2 means "doubling the velocity dispersion doubles
the extra pull". Two data sets at very different speeds test p:
  * X-COP's 12 clusters (galaxies at sigma ~ 500-1200 km/s): for each p, u is refitted with a and g_d held at the
    round-12 values, in the round-12 (static, deprojected) cluster sample.
  * KiDS's early/late lensing gap (lenses' stars at sigma ~ 100-190 km/s, SDSS-measured, round 12): with u from the
    clusters at the same p, the heat weights of red and blue lenses follow, and so does the gap (measured 0.153).
    <f sigma_e^p> is approximated as <f> (<f sigma_e^2> / <f>)^(p/2) within each 0.1-dex mass bin.
Clusters alone fix only a combination of u and p; the lenses, 7 times slower, pin p.

    python code/heat_exponent_v13.py --output run-coherent-force-v13/heat_exponent_v13.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import common                                     # noqa: E402,F401
import law as L                                   # noqa: E402
import run_v3 as R3                               # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import kids_static_v11 as KS                      # noqa: E402
import kids_heat_v12 as KH                        # noqa: E402
import distance_scale_v12 as DS                   # noqa: E402
from law_config import load_law                   # noqa: E402

HEAT = json.loads((HERE.parent / 'data/lens_heat_sdss_v12.json').read_text())


def cluster_M(c, a, u, lam, p):
    gN = R3.G * c['Mb'] / c['Rk'] ** 2
    k = 3.0 * (np.sqrt(c['sig2_star_hse']) / u) ** p
    S = R3.G * (c['W'] @ (k * c['dms'])) / c['Rk'] ** 2
    return L.total(gN, S, a, lam) * c['Rk'] ** 2 / R3.G


def fit_u(cls, a, lam, p):
    f = lambda lu: np.mean(R3.resid(cls, lambda c: cluster_M(c, a, 10 ** lu, lam, p)) ** 2)
    r = minimize_scalar(f, bounds=(0.5, 4.0), method='bounded', options=dict(xatol=1e-6))
    u = 10 ** r.x
    return u, float(np.sqrt(f(r.x)))


def k_lens(sample, u, p, lo=10.3, hi=10.9):
    rows = [r for r in HEAT['samples'][sample] if lo - 1e-9 <= r['logM_lo'] < hi - 1e-9]
    ks = []
    for r in rows:
        f = 1 - r['mean_1_minus_f']
        sig = np.sqrt(r['mean_f_sigma2'] / max(f, 1e-9))
        ks.append(3 * (f * (sig / u) ** p + r['mean_1_minus_f'] * (HEAT['sigma_disk_kms'] / u) ** p))
    return float(np.mean(ks))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    law = load_law('round12')
    C10.ALPHA = law['alpha_per_Mpc']
    cls = DS.xcop_static(0.95)
    a, lam = law['a_code'], law['lam']
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    rows = []
    for p in (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0):
        u, rms = fit_u(cls, a, lam, p)
        k = {s: k_lens(s, u, p) for s in ('red', 'blue', 'bulge', 'disc')}
        kk = KH.kids(consts, u, u * 1.0227121650537077 * 13.0, f, k=k)
        sig_mean = float(np.mean([np.sqrt(np.average(c['sig2_star_hse'], weights=c['dms'][:len(c['sig2_star_hse'])] if len(c['dms']) == len(c['sig2_star_hse']) else None)) for c in cls]))
        rows.append(dict(p=p, u_kms=u, xcop_rms=rms, k_red=k['red'], k_blue=k['blue'], gap_colour=kk['gap_model'], gap_sersic=kk['gap_sersic_model'],
                         kids_all=kk['all'], kids_red=kk['red'], kids_blue=kk['blue'], typical_cluster_sigma=sig_mean))
        print(f"p {p:4.2f}: u {u:7.1f} km/s, X-COP rms {rms:.4f}; k red {k['red']:.2f}, blue {k['blue']:.2f}; KiDS gap {kk['gap_model']:.3f} (colour, obs 0.153),"
              f" {kk['gap_sersic_model']:.3f} (Sersic, obs 0.154); all {kk['all']:+.3f}", flush=True)
    obs, err = 0.153, 0.04
    ps = np.array([r['p'] for r in rows]); gaps = np.array([r['gap_colour'] for r in rows]); rms = np.array([r['xcop_rms'] for r in rows])
    chi2 = ((gaps - obs) / err) ** 2
    res = dict(experiment='round 13: the exponent of the heat weight, from clusters and lenses together', rows=rows,
               kids_gap_observed=[obs, err], chi2_gap=chi2.tolist(),
               p_range_gap_within_1sigma=[float(ps[chi2 <= 1].min()) if (chi2 <= 1).any() else None, float(ps[chi2 <= 1].max()) if (chi2 <= 1).any() else None])
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    print('p within 1 sigma of the KiDS gap:', res['p_range_gap_within_1sigma'], '; wrote', args.output)


if __name__ == '__main__':
    main()
