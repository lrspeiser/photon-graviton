"""Round 12: the scale alpha of the static distance law (1 + z = e^(alpha D)), read from supernovae.

Pantheon+SH0ES standardized magnitudes with their full STAT+SYS covariance (the release the project already
uses in research_work/results/brightness-distance-consistency), calibrated on the 77 Cepheid-host rows exactly
as there:
* at a calibrator the geometric modulus mu (CEPH_DIST) holds, and the static law adds its two factors of
  S = e^(alpha D) (energy loss and arrival-rate stretching): M = m - mu - 5 log10(S);
* in the Hubble flow D = ln(1 + z_HD)/alpha and m = M + 25 + 5 log10(D (1 + z_hel)), D in Mpc.
The generalized least squares score keeps the calibrators' covariance with the far rows. Only alpha is fitted,
in windows of redshift: the SH0ES Hubble-flow range (0.0233-0.15), the lenses' range (0.1-0.3), both, and
all. A window's best alpha that drifts with depth measures the law's shape against the supernovae, not its scale.
The optional bounded beam-area term of bounded-area.py (flux x 1/(1 + eta f), f = 1 - 1/(1 + z)) is fitted
jointly with alpha as a comparison.

    python code/sn_scale_v12.py --output run-distance-scale-v12/sn_scale_v12.json
"""
from __future__ import annotations
import argparse, gzip, hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DATA = ROOT / 'shared_interaction_test/data'
BD = ROOT / 'research_work/results/brightness-distance-consistency'
ALPHA0 = 2.488993286382367e-4          # per Mpc, the adopted scale (collisions_v10.ALPHA)
C_KMS = 299792.458
WINDOWS = {'SH0ES Hubble flow, 0.0233-0.15': (0.0233, 0.15), 'the lenses, 0.1-0.3': (0.1, 0.3),
           '0.0233-0.3': (0.0233, 0.3), 'all, 0.0233-2.3': (0.0233, 2.3)}


def load():
    rec = json.loads((BD / 'opacity-results.json').read_text())['hashes']
    src = DATA / 'Pantheon+SH0ES.dat'
    assert hashlib.sha256(src.read_bytes()).hexdigest() == rec['Pantheon+SH0ES.dat']
    raw = gzip.decompress((DATA / 'Pantheon+SH0ES_STAT+SYS.cov.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest() == rec['Pantheon+SH0ES_STAT+SYS.cov']
    lines = src.read_text().split('\n'); head = lines[0].split()
    rows = [l.split() for l in lines[1:] if l.strip()]
    col = lambda c: np.array([float(r[head.index(c)]) for r in rows])
    df = {c: col(c) for c in ('zHD', 'zHEL', 'm_b_corr', 'CEPH_DIST', 'IS_CALIBRATOR')}
    v = np.array(raw.split(), dtype=float); n = int(v[0]); cov = v[1:].reshape(n, n)
    assert len(rows) == n
    return df, (cov + cov.T) / 2


class Fit:
    def __init__(self, df, cov, lo, hi):
        self.cal = np.flatnonzero(df['IS_CALIBRATOR'] == 1)
        z = df['zHD']
        self.ev = np.flatnonzero((df['IS_CALIBRATOR'] == 0) & (z >= lo) & (z < hi))
        Cc = cov[np.ix_(self.cal, self.cal)]
        w = cho_solve(cho_factor(Cc), np.ones(len(self.cal))); self.w = w / w.sum()
        cross = cov[np.ix_(self.ev, self.cal)] @ self.w
        V = cov[np.ix_(self.ev, self.ev)] + float(self.w @ Cc @ self.w) - cross[:, None] - cross[None, :]
        self.fac = cho_factor(V)
        self.mu_cal = df['CEPH_DIST'][self.cal]; self.m_cal = df['m_b_corr'][self.cal]
        self.D_cal = 10 ** ((self.mu_cal - 25) / 5)
        self.z = z[self.ev]; self.zh = df['zHEL'][self.ev]; self.m = df['m_b_corr'][self.ev]

    def resid(self, alpha, eta=0.0):
        area = lambda f: 1 + eta * f
        M = float(self.w @ (self.m_cal - self.mu_cal - 5 * alpha * self.D_cal / np.log(10) - 2.5 * np.log10(area(-np.expm1(-alpha * self.D_cal)))))
        D = np.log1p(self.z) / alpha
        return self.m - (M + 25 + 5 * np.log10(D * (1 + self.zh)) + 2.5 * np.log10(area(1 - 1 / (1 + self.zh)))), M

    def chi2(self, alpha, eta=0.0):
        r = self.resid(alpha, eta)[0]
        return float(r @ cho_solve(self.fac, r))


def fit_scale(F):
    """Best alpha / alpha0 and its Delta chi^2 = 1 range."""
    g = lambda s: F.chi2(ALPHA0 * s)
    best = minimize_scalar(g, bounds=(0.6, 1.4), method='bounded', options=dict(xatol=1e-7)).x
    c0 = g(best)
    lo = minimize_scalar(lambda s: (g(s) - c0 - 1) ** 2, bounds=(best - 0.2, best), method='bounded', options=dict(xatol=1e-8)).x
    hi = minimize_scalar(lambda s: (g(s) - c0 - 1) ** 2, bounds=(best, best + 0.2), method='bounded', options=dict(xatol=1e-8)).x
    return float(best), float(lo), float(hi), c0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    df, cov = load()
    res = dict(experiment='round 12: the static law\'s scale alpha from Pantheon+ supernovae, Cepheid-calibrated',
               alpha0_per_Mpc=ALPHA0, H0_like_adopted=ALPHA0 * C_KMS, windows=[])
    for name, (lo, hi) in WINDOWS.items():
        F = Fit(df, cov, lo, hi)
        s, sl, sh, c0 = fit_scale(F)
        r, M = F.resid(ALPHA0 * s)
        c_adopted = F.chi2(ALPHA0)
        row = dict(window=name, z_range=[lo, hi], n=len(F.ev), best_scale=s, scale_1sigma=[sl, sh], H0_like=s * ALPHA0 * C_KMS,
                   H0_like_1sigma=[sl * ALPHA0 * C_KMS, sh * ALPHA0 * C_KMS], chi2_best=c0, chi2_at_adopted=c_adopted, M=M)
        # the bounded beam-area term, fitted jointly (one more parameter)
        if hi > 0.2:
            opt = minimize(lambda p: F.chi2(ALPHA0 * p[0], p[1]) if p[1] >= 0 else 1e9, [s, 0.5], method='Nelder-Mead',
                           options=dict(xatol=1e-7, fatol=1e-7, maxiter=4000))
            row['with_area_term'] = dict(scale=float(opt.x[0]), eta=float(opt.x[1]), chi2=float(opt.fun), H0_like=float(opt.x[0] * ALPHA0 * C_KMS))
        res['windows'].append(row)
        extra = f"; with the area term: x{row['with_area_term']['scale']:.4f}, eta {row['with_area_term']['eta']:.2f}, chi2 {row['with_area_term']['chi2']:.1f}" if 'with_area_term' in row else ''
        print(f"{name:32s} n {len(F.ev):4d}: alpha x{s:.4f} (x{sl:.4f}-x{sh:.4f}), H0-like {s * ALPHA0 * C_KMS:.2f}; chi2 {c0:.1f} (adopted {c_adopted:.1f}){extra}", flush=True)
    # the calibration's zero point: other ladders move M, hence alpha, by 10^(dM/5)
    res['ladder_note'] = ('TRGB (Freedman 2021, 69.8) and JAGB (Lee et al. 2024, 67.8) calibrations sit 0.097 and 0.160 mag fainter '
                          'than SH0ES (73.04) in M; in our law they lower alpha by the same factors, x0.956 and x0.928.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
