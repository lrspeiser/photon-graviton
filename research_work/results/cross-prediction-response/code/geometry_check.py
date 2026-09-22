#!/usr/bin/env python3
"""JR-9 C: is the residual under-bending new physics, or the distance convention?

The cross-prediction angle error turned out to separate the six lenses almost
perfectly by foreground redshift. Before any of it is attributed to companion
structure, this script asks the ordinary question: what single change in the
lensing distance ratio D_ls/D_s would each lens need, and is that pattern what
a standard FLRW geometry would have given instead of the project's frozen
static-Euclidean convention?

JR-1 and everything downstream use "static Euclidean distances inferred using a
frozen exponential redshift rule". That is a declared model-dependent choice, not
a measurement. Substituting a different convention is not allowed as a silent
improvement, so this script only *reports* the comparison. It does not re-fit
anything and it does not change the project's primary distance convention.

    python geometry_check.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

import model as M
import response_kernels as rk

CKMS = 299792.458
OBS = {x['Name']: x for x in json.loads(
    (M.ROOT / 'research_work/results/lensing-data-readiness/'
     'lens-observations-and-image-models.json').read_text())}
CROSS = M.ROOT / 'research_work/results/cross-prediction-response/run-v1/cross-prediction-results.json'


def flrw_comoving(z, H0, Om, Ol):
    """Comoving distance in Mpc for flat or curved FLRW, by direct quadrature."""
    E = lambda zz: np.sqrt(Om * (1 + zz) ** 3 + (1 - Om - Ol) * (1 + zz) ** 2 + Ol)
    return CKMS / H0 * quad(lambda zz: 1 / E(zz), 0, z, epsabs=1e-10, epsrel=1e-10)[0]


def flrw_ratio(zl, zs, H0=70.0, Om=0.3, Ol=0.7):
    """D_ls / D_s for a flat FLRW cosmology (comoving distances subtract)."""
    dl, ds = flrw_comoving(zl, H0, Om, Ol), flrw_comoving(zs, H0, Om, Ol)
    return (ds - dl) / ds


def required_ratio_factor(lens, dm):
    """Multiplier on D_ls/D_s that makes the predicted angle equal the catalog angle."""
    original = lens.ratio
    try:
        def f(factor):
            lens.ratio = original * factor
            return M.einstein_angle(lens, dm) - lens.theta
        lo, hi = 0.2, 8.0
        if f(lo) * f(hi) > 0:
            return None
        return float(brentq(f, lo, hi, xtol=1e-12, rtol=1e-13))
    finally:
        lens.ratio = original


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    cross = json.loads(CROSS.read_text())
    kin_dm = {r['name']: r['K_to_L']['dm_MAP_dex'] for r in cross['companion_R10']}

    rows = []
    for name in M.LENSES:
        lens = M.make_lens(name)
        o = OBS[name]
        dm = kin_dm[name]
        factor = required_ratio_factor(lens, dm)
        project_ratio = float(lens.ratio)
        std = flrw_ratio(o['zFG'], o['zBG'])
        rows.append(dict(
            name=name, zFG=o['zFG'], zBG=o['zBG'],
            dm_from_kinematics_dex=dm,
            angle_observed_arcsec=float(lens.theta),
            angle_predicted_arcsec=M.einstein_angle(lens, dm),
            angle_fractional_error=float(M.einstein_angle(lens, dm) / lens.theta - 1),
            project_Dls_over_Ds=project_ratio,
            flrw_Dls_over_Ds_Om0p3_h0p70=std,
            flrw_over_project=float(std / project_ratio),
            required_ratio_factor=factor,
            required_fractional_change=None if factor is None else float(factor - 1)))

    need = np.array([r['required_fractional_change'] for r in rows], float)
    z = np.array([r['zFG'] for r in rows], float)
    flrw = np.array([r['flrw_over_project'] for r in rows], float) - 1
    # Descriptive linear trend of the required correction with lens redshift.
    Adesign = np.vstack([np.ones_like(z), z]).T
    coef, *_ = np.linalg.lstsq(Adesign, need, rcond=None)
    resid = need - Adesign @ coef
    # What a single constant rescaling would leave behind, for comparison.
    flat_resid = need - need.mean()

    payload = dict(
        experiment='JR-9C', protocol='../PROTOCOL.md',
        scope='Diagnostic only. The project distance convention is NOT changed here; '
              'this reports how much of the residual deflection gap is degenerate with it.',
        distance_convention_in_use='Static Euclidean, frozen exponential redshift rule '
                                   '(research_work/results/lensing-data-readiness/conditional-geometry.json)',
        comparison_cosmology='Flat FLRW, Om=0.3, OL=0.7, H0=70 km/s/Mpc, comoving-difference D_ls',
        rows=rows,
        required_change_vs_z_linear_fit=dict(intercept=float(coef[0]), slope_per_unit_z=float(coef[1]),
                                             residual_rms=float(np.sqrt(np.mean(resid ** 2))),
                                             raw_rms=float(np.sqrt(np.mean(need ** 2))),
                                             rms_after_constant_rescaling=float(np.sqrt(np.mean(flat_resid ** 2)))),
        correlation_required_change_with_zFG=float(np.corrcoef(need, z)[0, 1]),
        correlation_required_change_with_flrw_offset=float(np.corrcoef(need, flrw)[0, 1]),
        mean_flrw_over_project=float(np.mean(flrw + 1)),
        interpretation=[
            'A redshift-dependent distance error and a redshift-independent companion '
            'correction are not distinguishable by these six systems alone.',
            'Reporting both is mandatory; choosing whichever is smaller is not.'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'geometry-degeneracy.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print('lens          zFG     err%   need_Dls/Ds   project  FLRW(0.3,70)  FLRW/project')
    for r in rows:
        print(f"{r['name']:<12}{r['zFG']:>6.4f}{r['angle_fractional_error']*100:>8.2f}"
              f"{r['required_fractional_change']*100:>12.2f}%{r['project_Dls_over_Ds']:>10.4f}"
              f"{r['flrw_Dls_over_Ds_Om0p3_h0p70']:>12.4f}{r['flrw_over_project']:>13.4f}")
    f = payload['required_change_vs_z_linear_fit']
    print(f"\nrequired change vs zFG: {f['intercept']:+.4f} {f['slope_per_unit_z']:+.4f}*z, "
          f"residual RMS {f['residual_rms']*100:.2f}% (raw {f['raw_rms']*100:.2f}%, "
          f"constant rescale {f['rms_after_constant_rescaling']*100:.2f}%)")
    print(f"corr(required, zFG) = {payload['correlation_required_change_with_zFG']:+.3f}")
    print(f"corr(required, FLRW offset) = {payload['correlation_required_change_with_flrw_offset']:+.3f}")
    print('WROTE', out / 'geometry-degeneracy.json')


if __name__ == '__main__':
    main()
