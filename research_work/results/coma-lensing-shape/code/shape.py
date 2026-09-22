#!/usr/bin/env python3
"""JR-14: the law makes a sharp shape prediction for cluster weak lensing. Test it.

In the low-acceleration branch the force from a compact baryonic mass is

    g = sqrt(G M_b a0) / r

so the effective lensing mass grows LINEARLY with radius, the effective density goes
as 1/r^2, and the projected excess surface density that weak lensing measures is

    Delta Sigma (R)  =  (1/pi - 1/4) * w * sqrt(M_b a0 / G) / R     proportional to 1/R

with w = (1 + 1/tau)/2 the lensing weight of the extra term. The amplitude carries
tau and the unknown critical surface density together, so it cannot separate them.
The SHAPE can be tested on its own, and 1/R is a strong, parameter-free prediction.

Compared against an NFW shape, which is what dark matter predicts and what the
source paper fitted. NFW falls much faster than 1/R in the outskirts, so the two
separate exactly where these measurements are.

Data: six tangential-shear points reconstructed from Kubo et al. 2007 figure 2,
archived in this repository. Their limits are recorded there and repeated below.

    python shape.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar, minimize

ROOT = Path(__file__).resolve().parents[4]
KUBO = ROOT / 'research_work/results/cluster-observation-readiness/kubo-figure-data.json'


def nfw_delta_sigma(R, rs):
    """NFW excess surface density shape, Wright and Brainerd, up to amplitude."""
    x = np.asarray(R, float) / rs
    out = np.empty_like(x)
    lo, hi, mid = x < 1 - 1e-8, x > 1 + 1e-8, np.abs(x - 1) <= 1e-8
    if np.any(lo):
        a = x[lo]
        c = np.arccosh(1 / a)
        out[lo] = (8 * c / (a ** 2 * np.sqrt(1 - a ** 2)) + 4 / a ** 2 * np.log(a / 2)
                   - 2 / (a ** 2 - 1) + 4 * c / ((a ** 2 - 1) * np.sqrt(1 - a ** 2)))
    if np.any(hi):
        a = x[hi]
        c = np.arccos(1 / a)
        out[hi] = (8 * c / (a ** 2 * np.sqrt(a ** 2 - 1)) + 4 / a ** 2 * np.log(a / 2)
                   - 2 / (a ** 2 - 1) + 4 * c / (a ** 2 - 1) ** 1.5)
    if np.any(mid):
        out[mid] = 10.0 / 3.0 + 4 * np.log(0.5)
    return out


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    d = json.loads(KUBO.read_text())
    rows = d['rows']
    R = np.array([r['published_radius_h_inverse_Mpc'] for r in rows])
    g = np.array([r['shear_t'] for r in rows])
    s = np.array([r['plotted_sigma_t'] for r in rows])
    gx = np.array([r['shear_cross'] for r in rows])
    sx = np.array([r['plotted_sigma_cross'] for r in rows])

    def chi2_scaled(model):
        """Amplitude is a free scale; it carries tau and Sigma_crit together."""
        a = float(np.sum(model * g / s ** 2) / np.sum(model ** 2 / s ** 2))
        return float(np.sum(((a * model - g) / s) ** 2)), a

    one_over_R = 1.0 / R
    chi_law, amp_law = chi2_scaled(one_over_R)

    best_nfw = None
    for rs0 in (0.3, 0.7, 1.5, 3.0, 6.0):
        r = minimize_scalar(lambda lrs: chi2_scaled(nfw_delta_sigma(R, np.exp(lrs)))[0],
                            bracket=(np.log(rs0 * 0.5), np.log(rs0 * 2.0)))
        if best_nfw is None or r.fun < best_nfw.fun:
            best_nfw = r
    rs = float(np.exp(best_nfw.x))
    chi_nfw, amp_nfw = chi2_scaled(nfw_delta_sigma(R, rs))

    # free power law, to see what slope the data actually prefers
    def chi_pow(p):
        return chi2_scaled(R ** (-p))[0]
    pr = minimize(lambda x: chi_pow(x[0]), [1.0], method='Nelder-Mead',
                  options=dict(xatol=1e-6, fatol=1e-9))
    slope = float(pr.x[0])
    chi_pow_best, _ = chi2_scaled(R ** (-slope))

    chi_null = float(np.sum((g / s) ** 2))

    payload = dict(
        experiment='JR-14',
        scope='Shape-only test of the low-acceleration branch against reconstructed Coma '
              'tangential shear. The amplitude carries tau and Sigma_crit together and is '
              'marginalised, so this tests the radial shape and nothing else.',
        prediction='Delta Sigma proportional to 1/R, from g = sqrt(G M_b a0)/r with no '
                   'free shape parameter.',
        data=dict(source=d['source'], n_points=len(R),
                  radius_hinv_Mpc=[float(x) for x in R],
                  shear=[float(x) for x in g], error=[float(x) for x in s],
                  archived_limits=d.get('restrictions', [])),
        cross_component_control=dict(
            chi2=float(np.sum((gx / sx) ** 2)), paper_value=d['checks']['paper_null_chi2_cross'],
            note='The cross component should be consistent with zero if the signal is real '
                 'lensing rather than a systematic.'),
        results=dict(
            one_over_R=dict(chi2=chi_law, dof=len(R) - 1, shape_parameters=0,
                            reduced_chi2=chi_law / (len(R) - 1)),
            nfw=dict(chi2=chi_nfw, dof=len(R) - 2, shape_parameters=1,
                     scale_radius_hinv_Mpc=rs, reduced_chi2=chi_nfw / (len(R) - 2)),
            free_power_law=dict(chi2=chi_pow_best, fitted_slope=slope,
                                dof=len(R) - 2, shape_parameters=1),
            no_signal=dict(chi2=chi_null, dof=len(R),
                           paper_value=d['checks']['paper_null_chi2_t'])),
        reading='The free power law tells us the slope the data prefers on its own. If that '
                'sits near 1, the prediction is confirmed as a shape; it does not confirm tau, '
                'which needs the critical surface density.',
        limitations=[
            'Figure reconstruction, not author machine-readable measurements.',
            'No bin edges, source weights, full covariance or source geometry, so the errors '
            'are the plotted bars and correlations between bins are not modelled.',
            'The amplitude is marginalised, so nothing here measures tau. Separating tau needs '
            'the critical surface density, which is not recoverable from the figure.',
            'Coma baryonic mass is not used; only the shape is tested.'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'coma-shape.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    r_ = payload['results']
    print('Coma tangential shear, six reconstructed points\n')
    print(f"  no signal at all        chi2 {r_['no_signal']['chi2']:7.2f}  "
          f"({len(R)} dof)   [paper quotes {r_['no_signal']['paper_value']}]")
    print(f"  1/R, our prediction     chi2 {r_['one_over_R']['chi2']:7.2f}  "
          f"({r_['one_over_R']['dof']} dof)   reduced {r_['one_over_R']['reduced_chi2']:.2f}   "
          f"ZERO shape parameters")
    print(f"  NFW (dark matter)       chi2 {r_['nfw']['chi2']:7.2f}  "
          f"({r_['nfw']['dof']} dof)   reduced {r_['nfw']['reduced_chi2']:.2f}   "
          f"1 shape parameter, r_s = {rs:.2f}")
    print(f"  free power law R^-p     chi2 {r_['free_power_law']['chi2']:7.2f}   "
          f"fitted p = {slope:.3f}")
    print(f"\n  the data on its own prefers p = {slope:.3f}; the law predicts p = 1 with nothing fitted")
    print('WROTE', out / 'coma-shape.json')


if __name__ == '__main__':
    main()
