"""PF-1 T4: conditional exposed-data brightness test (protocol.md).

    python brightness.py [--output-dir DIR] [--canonical]

PF-1 predicts D_L = (1+z) ln(1+z)/alpha, the flux law already scored by
shared_interaction_test/brightness.py and brightness-distance-consistency/rate.py; the archived
frame convention is kept (heliocentric z for the energy/time factor, zHD for the path). The
comparator is flat FLRW with Omega_m = 0.3 and free H0. Both models get the same rows, the same
full covariance and the same freedom: one scale parameter, with the absolute magnitude taken from
the Cepheid hosts. The scale enters each prediction as a redshift-independent magnitude shift, so
the best fit equals the free-offset (shape-only) fit, and the calibrators only turn the offset into
a value of alpha or H0. These data are exposed; this is not a blind test.

Regenerates brightness-results.json into a fresh directory and compares it with the archived
copy; --canonical overwrites the archive.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
DATA = ROOT/'shared_interaction_test/data'
C_KMS = 299792.458
ALPHA0 = 0.0002488993286382367
BINS = [(.1, .3), (.3, .6), (.6, 1.), (1., 3.)]


def main():
    args = evidence_io.parse(__doc__)
    source, covsource = DATA/'Pantheon+SH0ES.dat', DATA/'Pantheon+SH0ES_STAT+SYS.cov'
    df = pd.read_csv(source, sep=r'\s+')
    raw = np.loadtxt(covsource)
    N = int(raw[0])
    C = raw[1:].reshape(N, N)
    C = (C + C.T)/2
    cal = np.flatnonzero(df.IS_CALIBRATOR.to_numpy() == 1)
    hi = np.flatnonzero((df.IS_CALIBRATOR.to_numpy() == 0) & (df.zHD.to_numpy() >= .1))
    Cc, Ch, Hc = C[np.ix_(cal, cal)], C[np.ix_(hi, hi)], C[np.ix_(hi, cal)]
    w = cho_solve(cho_factor(Cc), np.ones(len(cal)))
    w /= w.sum()
    mu_cal = df.CEPH_DIST.to_numpy()[cal]
    Mbase = float(w@(df.m_b_corr.to_numpy()[cal] - mu_cal))
    Dbar = float(w@10**((mu_cal - 25)/5))                  # GLS-weighted Cepheid distance, Mpc
    var = float(w@Cc@w)
    cross = Hc@w
    V = Ch + var - cross[:, None] - cross[None, :]           # calibrator-mean uncertainty propagated
    vf = cho_factor(V)
    ones = np.ones(len(hi))
    vi1 = cho_solve(vf, ones)
    z, zh, obs = df.zHD.to_numpy()[hi], df.zHEL.to_numpy()[hi], df.m_b_corr.to_numpy()[hi]

    shapes = dict(
        PF1=(1 + zh)*np.log1p(z),                            # alpha * D_L, Mpc per (1/Mpc)
        FLRW=(1 + zh)*C_KMS*np.array([quad(lambda s: 1/np.sqrt(.3*(1 + s)**3 + .7), 0, zz, epsabs=1e-11)[0]
                                      for zz in z]))        # H0 * D_L
    # Magnitude offset implied by each scale value, calibrators included:
    # PF-1 hosts are dimmed by their own propagation, M = Mbase - 5 alpha Dbar/ln10 (as in rate.py);
    # FLRW takes the Cepheid moduli as luminosity-distance moduli.
    offset_of = dict(PF1=lambda a: Mbase - 5*a*Dbar/np.log(10) - 5*np.log10(a),
                     FLRW=lambda H: Mbase - 5*np.log10(H))
    bounds = dict(PF1=(1/C_KMS, 1000/C_KMS), FLRW=(1., 1000.))

    models = {}
    for name, shape in shapes.items():
        r0 = obs - (25 + 5*np.log10(shape))
        off = float(vi1@r0/(ones@vi1))
        sig = float((ones@vi1)**-.5)
        r = r0 - off
        chi2 = float(r@cho_solve(vf, r))
        f = offset_of[name]
        inv = lambda target: brentq(lambda s: f(s) - target, *bounds[name], xtol=1e-16)
        scale = inv(off)
        interval = sorted([inv(off - sig), inv(off + sig)])
        bins, B = [], []
        for lo, up in BINS:
            ix = np.flatnonzero((z >= lo) & (z < up))
            sub = V[np.ix_(ix, ix)]
            bw = cho_solve(cho_factor(sub), np.ones(len(ix)))
            bw /= bw.sum()
            full = np.zeros(len(hi))
            full[ix] = bw
            B.append(full)
            bins.append(dict(z_min=lo, z_max=up, light_curves=len(ix), mean_residual_mag=float(bw@r[ix]),
                             sigma_mag=float(np.sqrt(bw@sub@bw))))
        B = np.array(B)
        bc = B@V@B.T
        means = B@r
        contrasts = [dict(bin_index=i, minus_first_bin_mag=float(means[i] - means[0]),
                          sigma_mag=float(np.sqrt(bc[i, i] + bc[0, 0] - 2*bc[i, 0]))) for i in range(1, len(BINS))]
        scale_key = 'c_alpha_km_s_Mpc' if name == 'PF1' else 'H0_km_s_Mpc'
        conv = C_KMS if name == 'PF1' else 1.
        models[name] = {'chi_squared': chi2, 'rows': len(hi), 'fitted_parameters': 1,
                        'best_offset_mag': off, 'offset_sigma_mag': sig,
                        scale_key: scale*conv, scale_key + '_offset_sigma_interval': [v*conv for v in interval],
                        'redshift_bins': bins, 'contrasts_to_first_bin': contrasts}
    result = dict(
        test='PF-1 T4, conditional exposed-data brightness test',
        exposed_data=True, blind=False,
        pf1_luminosity_distance='D_L = (1+zHEL) ln(1+zHD)/alpha', comparator='flat FLRW, Omega_m = 0.3 fixed, H0 free',
        calibrators=dict(light_curves=len(cal), Mbase=Mbase, Mbase_sigma=var**.5, weighted_distance_Mpc=Dbar),
        models=models, delta_chi_squared_PF1_minus_FLRW=models['PF1']['chi_squared'] - models['FLRW']['chi_squared'],
        archived_alpha_c_km_s_Mpc=C_KMS*ALPHA0,
        caveats=['Standardized SALT2 magnitudes and the released covariance inherit Pantheon+ reduction assumptions, '
                 'including a fiducial cosmology in some corrections.',
                 'Cepheid distances are stipulated geometric distances with their own ladder assumptions.',
                 'The same flux law was scored before (shared_interaction_test/brightness.py, rate.py); this run '
                 'repeats it under PF-1 with a matched FLRW comparator and is not new evidence.'],
        input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, covsource, HERE/'protocol.md')})
    text = json.dumps(result, indent=1) + '\n'
    print(text)
    return evidence_io.finish(args, 'propagation-field-brightness', text, HERE/'brightness-results.json')


if __name__ == '__main__':
    raise SystemExit(main())
