"""Post-hoc diagnostic for CR-2, requested in the owner's review; not part of the declared protocol.

The CR-2 report attributed the condensate's Milky Way shortfall to its edge near 75 kpc. The model adds
G M(<r)/r to the baryonic v^2 at the 38 Eilers radii (5–25 kpc), so material beyond 25 kpc cannot change it.
This script compares, inside the measured range, three things:
- the extra enclosed mass each fitted model supplies;
- the extra mass the observed speeds require, M_req(<r) = r (v_obs^2 - v_b^2)/G, with the archive's baseline-I
  baryon speeds;
- the signed velocity residuals.
The models are the CR-2 condensate at its lens-fixed R_TF and the shared-scale NFW at its lens-fixed r_s, each
with only its mass fitted to the Milky Way.

    python milky-way-interior.py

Writes milky-way-interior.json next to this script.
"""
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cr2  # noqa: E402
import tf  # noqa: E402
import inputs as I  # noqa: E402  (capture-to-orbit, on the path via cr2)


def main():
    res = json.loads((HERE/'cr2-results.json').read_text(encoding='utf-8'))
    uni = json.loads((HERE/'universal-halo.json').read_text(encoding='utf-8'))
    base = next(r for r in I.milky_way_runs() if r['baryons'] == 'I' and abs(r['rd'] - 2.6) < 1e-9 and abs(r['lf'] - 1) < 1e-9)
    R, y, vb = base['R'], base['y'], base['vb']
    Mb = I.milky_way_receivers('I')['mass']
    M_req = R*(y**2 - vb**2)/cr2.G
    out = dict(scope='Post-hoc diagnostic of the CR-2 Milky Way comparison, requested in review.', R_kpc=R.tolist(),
               required_extra_mass_Msun=M_req.tolist(), geometries={})
    bands = [(5, 10), (10, 15), (15, 20), (20, 25.5)]
    for which in cr2.GEOMETRIES:
        mw = res['geometries'][which]['milky_way']
        fam = tf.Family(I.GRID, Mb/Mb[-1], tf.k_of(mw['R_tf_kpc']))
        s = fam.solve(mw['lam'])
        M_c = Mb[-1]*s['m'](fam.k*R)/mw['lam']
        nb = uni['geometries'][which]['best']
        rs, A = nb['rs_kpc'], nb['milky_way']['halo_normalization_Msun']
        x = R/rs
        M_n = A*(np.log1p(x) - x/(1 + x))
        rows = {}
        for label, M in (('baryons_only', np.zeros_like(R)), ('condensate', M_c), ('shared_nfw', M_n)):
            v = np.sqrt(vb**2 + cr2.G*M/R)
            res_v = v - y
            sel = (R >= 8) & (R <= 20)
            slope = float(np.polyfit(np.log(R[sel]), np.log(np.maximum(M[sel], 1e-30)), 1)[0]) if M.any() else None
            rows[label] = dict(
                rmse_kms=float(np.sqrt(np.mean(res_v**2))),
                mean_signed_residual_by_band={f'{a:g}-{b:g} kpc': float(np.mean(res_v[(R >= a) & (R < b)])) for a, b in bands},
                extra_mass_over_required={f'{r0:g} kpc': float(np.interp(r0, R, M)/np.interp(r0, R, M_req)) for r0 in (6, 10, 15, 20, 24)},
                enclosed_mass_log_slope_8_20kpc=slope,
                signed_residual_kms=res_v.tolist())
        sel = (R >= 8) & (R <= 20)
        rows['required_log_slope_8_20kpc'] = float(np.polyfit(np.log(R[sel]), np.log(M_req[sel]), 1)[0])
        rows['condensate_edge_kpc'] = s['x_edge']/fam.k
        rows['condensate_mass_inside_25kpc_over_total'] = float(Mb[-1]*s['m'](np.array([fam.k*25.]))[0]/mw['lam']/mw['condensate_mass_Msun'])
        out['geometries'][which] = rows
        print(which, json.dumps({k: (v if not isinstance(v, dict) else {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != 'signed_residual_kms'}) for k, v in rows.items()}, indent=1, default=float))
    (HERE/'milky-way-interior.json').write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8', newline='\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
