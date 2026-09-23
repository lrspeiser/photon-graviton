#!/usr/bin/env python3
"""Round 7: does a later switch-off suit the SPARC galaxies as well as the Milky Way?

    python transition_check_v7.py --output-dir ../run-transition-check-v7

Round 7 found our law 8% low at the Sun (211 against 229-234 km/s) with McMillan's (2017) visible
matter, and 0.03 dex low at the median SPARC point of the same pull. The Milky Way curve rises to
218-221 km/s if the release scale g_d is 1.5-2 times larger (milky_way_v7.py). This script asks
whether the 149 SPARC galaxies allow that: the typical miss (km/s), the fit statistic used in
round 3 (mean squared log residual) and the median residual by pull, for g_d x 1, 1.5, 2 and 3,
with a held fixed and with a refitted. Nothing is adopted here; it maps what would have to change.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import run as R
import law as L


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    c = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    a, lam, u = c['a_code'], c['lam'], c['u_kms']
    gals = R.load_sparc(); K = L.KMS2_PER_KPC

    def stats(a_, lam_):
        rows = [np.c_[g['gN'] * K, (g['v'] ** 2 / g['r']) / R.galaxy_g(g, a_, u, lam_)] for g in gals]
        A = np.vstack(rows); lr = np.log10(A[:, 1]); lg = np.log10(A[:, 0])
        rms_kms, msq = R.sparc_score(gals, lambda g: R.galaxy_g(g, a_, u, lam_))
        bins = {f'{lo}..{hi}': float(np.median(lr[(lg >= lo) & (lg < hi)])) for lo, hi in ((-12, -11), (-11, -10.25), (-10.25, -9.5), (-9.5, -8.75))}
        return dict(rms_kms=float(rms_kms), mean_sq_log=float(msq), median=float(np.median(lr)), median_by_pull=bins)

    res = []
    for gs in (1.0, 1.5, 2.0, 3.0):
        fixed = stats(a, lam * gs)
        gd = lam * gs * a
        r = minimize_scalar(lambda la: R.sparc_score(gals, lambda g: R.galaxy_g(g, 10 ** la, u, gd / 10 ** la))[1], bounds=(3.0, 3.6), method='bounded')
        a2 = 10 ** r.x; refit = stats(a2, gd / a2)
        res.append(dict(g_d_scale=gs, g_d_SI=c['g_d_SI'] * gs, a_fixed=fixed, a_refit=dict(a_over_a0=a2 / a, **refit)))
        print(f"g_d x{gs}: a fixed: {fixed['rms_kms']:.2f} km/s, median {fixed['median']:+.3f}, by pull {', '.join(f'{k} {v:+.3f}' for k, v in fixed['median_by_pull'].items())}; "
              f"a refit x{a2 / a:.3f}: {refit['rms_kms']:.2f} km/s, median {refit['median']:+.3f}", flush=True)
    (out / 'transition_check_v7.json').write_text(json.dumps(dict(
        experiment='SPARC typical miss and median residuals for a later switch-off (g_d x 1-3), round 7', constants=c, results=res,
        seconds=time.monotonic() - t0), indent=1) + '\n')


if __name__ == '__main__':
    main()
