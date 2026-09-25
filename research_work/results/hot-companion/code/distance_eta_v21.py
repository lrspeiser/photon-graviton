"""Round 21, step 3a: the supernovae and the proposed distance law's extra factor.

A proposal supplied to the project (25 September 2026) rereads round 12's bounded beam-area term, flux x 1/(1 + eta b)
with b = z/(1 + z), as a correction to the path distance itself,
    D*(z) = (ln(1 + z) / alpha) sqrt(1 + eta z/(1 + z)),   D_A = D*/(1 + z),   D_L = (1 + z) D*,
which gives the same supernova magnitudes (D_L enters as 5 log10 D_L) while keeping distance duality,
D_L = (1 + z)^2 D_A. It notes that b = z/(1 + z) = 1 - e^(-y), y = ln(1 + z), is the solution of db/dy = 1 - b (a state
relaxing at the redshift's own rate), and that eta = 1/2 may be allowed.

Here, with round 12's Pantheon+ reduction (code/sn_scale_v12.py: the Cepheid-calibrated standardized magnitudes, the
full STAT+SYS covariance, generalized least squares):
  1. the profile chi^2(eta), alpha refitted at each eta, for all 1,365 supernovae at z > 0.023 (and by depth);
  2. eta = 1/2 against the best eta;
  3. for reference only, flat LCDM (Omega_m and the scale fitted) in the same reduction; that model is excluded as an
     answer by the project's rules and is shown only to calibrate what a chi^2 difference means;
  4. the proposal's operation as written, a_out(w) = sqrt(s) a_in(s w), conserves each pulse's energy (and multiplies
     the photon number by s = 1 + z), so the flux is dimmed by the stretching alone: D_L = (1 + z)^(1/2) D*. The
     distances above assume D_L = (1 + z) D*, which needs a_out(w) = a_in(s w) (photons conserved, each redshifted).
     Both are fitted, and the power p of (1 + z) is also left free.

    python code/distance_eta_v21.py --output run-distance-eta-v21/distance_eta_v21.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import cho_solve
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sn_scale_v12 as SN                # noqa: E402

C_KMS = SN.C_KMS


def profile_eta(F, etas):
    rows = []
    s_prev = 0.98
    for eta in etas:
        r = minimize_scalar(lambda s: F.chi2(SN.ALPHA0 * s, eta), bounds=(0.8, 1.2), method='bounded', options=dict(xatol=1e-8))
        rows.append(dict(eta=float(eta), scale=float(r.x), H0_like=float(r.x * SN.ALPHA0 * C_KMS), chi2=float(r.fun)))
    return rows


def chi2_power(F, alpha, eta, p):
    """D_L = (1 + z)^p D*. p = 1: the photons are conserved and each arrives with 1/(1 + z) of its energy (what the
    distances need); p = 1/2: the pulse's energy is conserved, as in the proposal's unitary operation
    a_out(w) = sqrt(s) a_in(s w), so only the stretching (1 + z) dims the flux (the photon number grows by 1 + z)."""
    area = lambda f: 1 + eta * f
    M = float(F.w @ (F.m_cal - F.mu_cal - 5 * p * alpha * F.D_cal / np.log(10) - 2.5 * np.log10(area(-np.expm1(-alpha * F.D_cal)))))
    D = np.log1p(F.z) / alpha
    r = F.m - (M + 25 + 5 * np.log10(D * (1 + F.zh) ** p) + 2.5 * np.log10(area(1 - 1 / (1 + F.zh))))
    return float(r @ cho_solve(F.fac, r))


def lcdm_chi2(F, H0, Om):
    """flat LCDM in the same reduction: calibrators take the geometric modulus as it is; far rows D_L = (1 + z_hel) chi(z_HD)."""
    M = float(F.w @ (F.m_cal - F.mu_cal))
    zg = np.linspace(0, max(F.z.max(), 2.5), 6000)
    E = np.sqrt(Om * (1 + zg) ** 3 + 1 - Om)
    chi = C_KMS / H0 * cumulative_trapezoid(1 / E, zg, initial=0.0)
    D = np.interp(F.z, zg, chi) * (1 + F.zh)
    r = F.m - (M + 25 + 5 * np.log10(D))
    return float(r @ cho_solve(F.fac, r)), M


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    df, cov = SN.load()
    out = dict(experiment='round 21: the supernovae and the path-distance factor sqrt(1 + eta z/(1 + z))', windows={})
    for name, (lo, hi) in (('all, 0.0233-2.3', (0.0233, 2.3)), ('0.0233-0.15', (0.0233, 0.15)), ('0.15-2.3', (0.15, 2.3))):
        F = SN.Fit(df, cov, lo, hi)
        etas = np.round(np.arange(0.0, 1.201, 0.01), 3)
        prof = profile_eta(F, etas)
        c = np.array([p['chi2'] for p in prof]); i = int(np.argmin(c))
        # refine the minimum and the Delta chi^2 = 1 interval
        best = minimize(lambda q: F.chi2(SN.ALPHA0 * q[0], q[1]), [prof[i]['scale'], prof[i]['eta']], method='Nelder-Mead',
                        options=dict(xatol=1e-8, fatol=1e-8, maxiter=6000))
        c_best = float(best.fun); eta_best = float(best.x[1])
        within = [p['eta'] for p in prof if p['chi2'] - c_best <= 1.0]
        at_half = [p for p in prof if abs(p['eta'] - 0.5) < 1e-9][0]
        at_zero = prof[0]
        row = dict(n=len(F.ev), eta_best=eta_best, scale_best=float(best.x[0]), H0_like_best=float(best.x[0] * SN.ALPHA0 * C_KMS),
                   chi2_best=c_best, eta_1sigma=[min(within), max(within)] if within else None,
                   eta_half=dict(chi2=at_half['chi2'], delta_chi2=at_half['chi2'] - c_best, H0_like=at_half['H0_like']),
                   eta_zero=dict(chi2=at_zero['chi2'], delta_chi2=at_zero['chi2'] - c_best, H0_like=at_zero['H0_like']),
                   profile=prof)
        # the reference model, same reduction
        opt = minimize(lambda q: lcdm_chi2(F, q[0], q[1])[0] if 0.01 < q[1] < 1.0 else 1e12, [70.0, 0.3], method='Nelder-Mead',
                       options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
        row['reference_flat_lcdm'] = dict(H0=float(opt.x[0]), Om=float(opt.x[1]), chi2=float(opt.fun))
        # the operation's normalization: energy-conserving (p = 1/2) against photon-conserving (p = 1), and p free
        ops = {}
        assert abs(chi2_power(F, SN.ALPHA0 * best.x[0], eta_best, 1.0) - c_best) < 1e-6
        for label, p_fixed in (('pulse energy conserved, p = 1/2', 0.5), ('p free', None)):
            fun = (lambda q: chi2_power(F, SN.ALPHA0 * q[0], q[1], p_fixed if p_fixed is not None else q[2]) if q[1] > -0.99 else 1e12)
            start = [best.x[0], eta_best] + ([] if p_fixed is not None else [1.0])
            o = minimize(fun, start, method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-8, maxiter=20000))
            if p_fixed is not None:     # a second start, far from the first
                o2 = minimize(fun, [best.x[0], 3.0], method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-8, maxiter=20000))
                o = o2 if o2.fun < o.fun else o
            ops[label] = dict(scale=float(o.x[0]), eta=float(o.x[1]), p=float(p_fixed if p_fixed is not None else o.x[2]),
                              H0_like=float(o.x[0] * SN.ALPHA0 * C_KMS), chi2=float(o.fun), delta_chi2=float(o.fun - c_best))
        row['operation_normalization'] = ops
        out['windows'][name] = row
        print(f"{name:18s} n {len(F.ev):4d}: best eta {eta_best:.3f} (1 sigma {row['eta_1sigma']}), chi2 {c_best:.1f}, H0-like {row['H0_like_best']:.2f}; "
              f"eta 0.5: Delta chi2 {row['eta_half']['delta_chi2']:.2f}; eta 0: Delta chi2 {row['eta_zero']['delta_chi2']:.1f}; "
              f"reference flat LCDM chi2 {opt.fun:.1f} (Om {opt.x[1]:.3f}, H0 {opt.x[0]:.2f})", flush=True)
        for label, o in ops.items():
            print(f"    {label:32s}: p {o['p']:.3f}, eta {o['eta']:.3f}, H0-like {o['H0_like']:.2f}, chi2 {o['chi2']:.1f} (Delta {o['delta_chi2']:+.1f})", flush=True)
    # the relaxing state: b(y) solving db/dy = 1 - b from b(0) = 0 is 1 - e^-y = z/(1 + z) exactly
    y = np.linspace(0, np.log(3.3), 200); b = np.zeros_like(y)
    for j in range(1, len(y)):
        h = y[j] - y[j - 1]; bj = b[j - 1]
        k1 = 1 - bj; k2 = 1 - (bj + h * k1 / 2); k3 = 1 - (bj + h * k2 / 2); k4 = 1 - (bj + h * k3)
        b[j] = bj + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    z = np.expm1(y)
    out['relaxing_state_check'] = dict(max_abs_difference_from_z_over_1_plus_z=float(np.max(np.abs(b - z / (1 + z)))))
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    print('relaxing state: max |b - z/(1+z)| =', out['relaxing_state_check']['max_abs_difference_from_z_over_1_plus_z'])
    print('wrote', args.output)


if __name__ == '__main__':
    main()
