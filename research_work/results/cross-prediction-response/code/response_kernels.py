#!/usr/bin/env python3
"""JR-9 D: where would extra gravitating structure have to live, and what does it cost?

For each lens we add a thin spherical shell of mass m at 3D radius r_s on top of
the frozen R10 response, and ask two questions with numbers:

  1. How much shell mass at r_s closes the Einstein-angle gap exactly?
  2. What chi-square does that cost against the measured V_rms bins, after the
     ordinary nuisance pair is allowed to re-adjust in its own defence?

The deflection integral and the Jeans pressure integral are both linear in the
force, so the shell response is the linear-response kernel of each observable.
Any spherically symmetric redistribution is a superposition of these shells, so
the minimum cost over r_s is a lower bound on the kinematic price of repairing
that lens with *any* spherical rearrangement of the companion.

    python response_kernels.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.linalg import cho_factor, cho_solve

import model as M

G = M.J.G
SHELL_RADII_PER_DECADE = 24
SHELL_RANGE_IN_RE = (0.02, 400.0)


def shell_force(r, mass, r_s):
    """Newtonian force from a thin shell: zero inside, G m / r^2 outside."""
    r = np.asarray(r, float)
    return np.where(r >= r_s, G * mass / np.maximum(r, 1e-30) ** 2, 0.0)


def base_terms(lens, dm):
    """R10 baryonic + companion force pieces for this lens at one stellar offset."""
    p = M.params(dm, 0.)
    local = lens.local_parameters(p, M.SPEC)
    q = local['q_sph']
    A, rc, rt = M.J.source_parameters(lens.Mstar, 0., lens.Re, local, M.SPEC['saturation'], 1.)
    scale = 10 ** local['logu']
    return scale, float(A), float(rc), float(rt), float(q)


def force_on(lens, r, dm):
    scale, A, rc, rt, q = base_terms(lens, dm)
    return scale * lens.baryon_force(r) + M.J.companion_force(np.asarray(r, float), A, rc, rt, q)


def vrms_with_shell(lens, dm, beta, mass, r_s):
    """Seeing-convolved annular V_rms including a shell perturbation."""
    scale, A, rc, rt, q = base_terms(lens, dm)
    g = scale * lens.gb + M.J.companion_force(lens.r, A, rc, rt, q) + shell_force(lens.r, mass, r_s)
    p = M.params(dm, beta)
    local = lens.local_parameters(p, M.SPEC)
    try:
        return lens.moments(g, local, M.SPEC['radial_beta'])
    except ValueError:
        return None


def angle_with_shell(lens, dm, mass, r_s):
    """Einstein angle including a shell perturbation, by the same root find as JR-1."""
    def f(b):
        r = b / np.cos(lens.lens_t)
        g = force_on(lens, r, dm) + shell_force(r, mass, r_s)
        return lens.ratio * 4 / M.J.CLIGHT ** 2 * np.dot(lens.lens_w, g * r) - b / lens.Dl
    grid = lens.Re * np.geomspace(1e-6, 1e4, 200)
    vals = np.array([f(b) for b in grid])
    ii = np.flatnonzero((vals[:-1] > 0) & (vals[1:] <= 0))
    if not len(ii):
        return None
    i = ii[-1]
    return float(brentq(f, grid[i], grid[i + 1], xtol=1e-10) / lens.Dl * M.J.ARCSEC)


def kinematic_map(lens, chol):
    """Best ordinary nuisance pair from kinematics alone (no shell). Reference point."""
    def cost(x):
        dm, beta = float(x[0]), float(np.clip(x[1], *M.BETA_BOUNDS))
        v = M.vrms(lens, dm, beta)
        if v is None:
            return 1e9
        d = v - lens.y
        return float(d @ cho_solve(chol, d)) + (dm / lens.mass_log_error) ** 2 \
            + (beta / M.BETA_PRIOR_SIGMA) ** 2
    best = min((minimize(cost, [d0, b0], method='Nelder-Mead',
                         options=dict(xatol=1e-6, fatol=1e-8, maxiter=800))
                for d0 in (-0.05, 0.0, 0.08) for b0 in (-0.3, 0.0, 0.3)), key=lambda r: r.fun)
    return float(best.x[0]), float(np.clip(best.x[1], *M.BETA_BOUNDS)), float(best.fun)


def shell_mass_for_angle(lens, dm, r_s, target):
    """Shell mass at r_s that makes the predicted Einstein angle equal the catalog value."""
    f = lambda lm: (angle_with_shell(lens, dm, lm * lens.Mstar, r_s) or -1e9) - target
    lo, hi = -50.0, 5000.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    return float(brentq(f, lo, hi, xtol=1e-12, rtol=1e-13)) * lens.Mstar


def companion_enclosed_mass(lens, dm, r):
    """Effective enclosed mass of the frozen R10 companion, r^2 g_chi / G."""
    _, A, rc, rt, q = base_terms(lens, dm)
    r = np.asarray(r, float)
    return r ** 2 * M.J.companion_force(r, A, rc, rt, q) / G


def run_lens(name):
    lens = M.make_lens(name)
    chol = cho_factor(np.asarray(lens.cov), lower=True)
    dm0, beta0, cost0 = kinematic_map(lens, chol)
    v0 = M.vrms(lens, dm0, beta0)
    d0 = v0 - lens.y
    chi2_0 = float(d0 @ cho_solve(chol, d0))
    theta0 = M.einstein_angle(lens, dm0)

    lo, hi = (x * lens.Re for x in SHELL_RANGE_IN_RE)
    n = int(SHELL_RADII_PER_DECADE * np.log10(hi / lo)) + 1
    radii = np.geomspace(lo, hi, n)

    # ---- linear-response kernels, per unit of added mass -----------------------
    probe = 1e-3 * lens.Mstar
    kernels = []
    for r_s in radii:
        a = angle_with_shell(lens, dm0, probe, r_s)
        v = vrms_with_shell(lens, dm0, beta0, probe, r_s)
        if a is None or v is None:
            kernels.append(None)
            continue
        dv = (v - v0) / probe
        kernels.append(dict(
            r_kpc=float(r_s), r_over_Re=float(r_s / lens.Re),
            dtheta_dM_arcsec_per_Msun=float((a - theta0) / probe),
            dlnTheta_dlnMstar=float((a - theta0) / theta0 * lens.Mstar / probe),
            dvrms_dM_kms_per_Msun=dv.tolist(),
            max_abs_dlnVrms_dlnMstar=float(np.max(np.abs(dv * lens.Mstar / lens.y))),
            # leverage: fractional lensing gain per unit of worst fractional V_rms disturbance
            leverage=float(abs((a - theta0) / theta0) / max(np.max(np.abs(dv / lens.y)), 1e-300))))

    # ---- cost of actually closing the angle gap at each radius ----------------
    costs = []
    for r_s in radii:
        m = shell_mass_for_angle(lens, dm0, r_s, lens.theta)
        if m is None:
            costs.append(None)
            continue

        def penalty(x):
            dm, beta = float(x[0]), float(np.clip(x[1], *M.BETA_BOUNDS))
            v = vrms_with_shell(lens, dm, beta, m, r_s)
            if v is None:
                return 1e9
            d = v - lens.y
            return float(d @ cho_solve(chol, d)) + (dm / lens.mass_log_error) ** 2 \
                + (beta / M.BETA_PRIOR_SIGMA) ** 2
        # The shell is held fixed; only ordinary nuisance re-adjusts in its defence.
        best = min((minimize(penalty, [dm0 + s, beta0], method='Nelder-Mead',
                             options=dict(xatol=1e-6, fatol=1e-8, maxiter=600))
                    for s in (-0.05, 0.0, 0.05)), key=lambda r: r.fun)
        dm_b, beta_b = float(best.x[0]), float(np.clip(best.x[1], *M.BETA_BOUNDS))
        v = vrms_with_shell(lens, dm_b, beta_b, m, r_s)
        d = v - lens.y
        costs.append(dict(
            r_kpc=float(r_s), r_over_Re=float(r_s / lens.Re),
            shell_mass_Msun=float(m), shell_over_stellar_mass=float(m / lens.Mstar),
            companion_enclosed_mass_Msun=float(companion_enclosed_mass(lens, dm0, r_s)),
            shell_over_companion_enclosed=float(m / companion_enclosed_mass(lens, dm0, r_s)),
            within_redistribution_budget=bool(
                m >= 0 or abs(m) <= companion_enclosed_mass(lens, dm0, r_s)),
            dm_after_dex=dm_b, beta_after=beta_b,
            kinematic_chi2=float(d @ cho_solve(chol, d)),
            delta_chi2_vs_kinematic_best=float(d @ cho_solve(chol, d) - chi2_0),
            total_penalty=float(best.fun), delta_penalty=float(best.fun - cost0),
            vrms_fractional_RMS=float(np.sqrt(np.mean((d / lens.y) ** 2)))))

    live = [c for c in costs if c]
    cheapest = min(live, key=lambda c: c['delta_penalty']) if live else None
    return dict(
        name=name, theta_observed_arcsec=float(lens.theta), theta_at_kinematic_best=float(theta0),
        angle_fractional_error_at_kinematic_best=float(theta0 / lens.theta - 1),
        Re_kpc=float(lens.Re), Mstar_Msun=float(lens.Mstar), n_bins=int(len(lens.y)),
        einstein_radius_kpc=float(lens.b), outer_kinematic_radius_kpc=float(
            lens.outer_arcsec[-1] * lens.Dl / M.J.ARCSEC),
        kinematic_best=dict(dm_dex=dm0, beta=beta0, chi2=chi2_0, penalty=cost0),
        kernels=[k for k in kernels if k], shell_costs=live,
        cheapest_repair=cheapest)


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory; old results are never overwritten')
    out.mkdir(parents=True)
    t0 = time.monotonic()
    rows = [run_lens(n) for n in M.LENSES]
    for r in rows:
        c = r['cheapest_repair']
        print(f"{r['name']}: gap {r['angle_fractional_error_at_kinematic_best']:+.2%}"
              f"  cheapest shell at r/Re={c['r_over_Re']:.3f}"
              f"  m/M*={c['shell_over_stellar_mass']:+.3f}"
              f"  dchi2={c['delta_chi2_vs_kinematic_best']:+.2f}"
              f"  dpenalty={c['delta_penalty']:+.2f}", flush=True)
    payload = dict(
        experiment='JR-9D', protocol='../PROTOCOL.md',
        scope='Linear response kernels of the two observables to added spherical mass, '
              'and the kinematic chi-square price of closing each Einstein-angle gap.',
        shell_radius_range_in_Re=list(SHELL_RANGE_IN_RE),
        shell_radii_per_decade=SHELL_RADII_PER_DECADE,
        note='The shell basis spans every spherically symmetric redistribution because both '
             'the deflection integral and the Jeans pressure integral are linear in the force.',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        lenses=rows, seconds=time.monotonic() - t0)
    (out / 'response-kernels.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    print('WROTE', out / 'response-kernels.json', flush=True)


if __name__ == '__main__':
    main()
