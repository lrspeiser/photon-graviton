#!/usr/bin/env python3
"""JR-10: the companion's stress state, not its amount, sets lensing versus kinematics.

The energy ledger forbids making the companion's gravity out of stored starlight by
four to eight orders of magnitude (RC-1, and the per-lens starlight check in JR-9).
So the companion cannot be static dust whose mass we must pay for. This module tests
the alternative that costs no extra energy at all: the same stress-energy, arranged
anisotropically.

In linearised GR with ds^2 = -(1+2 Phi) dt^2 + (1-2 Psi) dx^2,

    laplacian(Phi) = 4 pi G (rho + T^k_k)      slow stars feel Phi
    laplacian(Psi) = 4 pi G  rho               light feels Phi + Psi

so for a companion with radial and transverse stress ratios w_r = p_r/rho and
w_t = p_t/rho,

    gamma_chi  ==  Psi_chi / Phi_chi  =  1 / (1 + w_r + 2 w_t)

    static dust            w_r=0,  w_t=0   -> gamma = 1      (what R10 assumed)
    radial null flux       w_r=1,  w_t=0   -> gamma = 1/2
    radial tension         w_r=-1/2,w_t=0  -> gamma = 2
    static field gradient  w_r=1,  w_t=-1  -> gamma = infinite

The last case is exact and is the one the data has been asking for: T^k_k = -rho
cancels rho in the Tolman active mass, so a static radial field gradient exerts NO
force on slow stars while its energy density still deflects light. That is lensing
without kinematics, paid for out of stress rather than out of new energy.

R10's g_chi was fitted to stellar kinematics, so it already *is* the Phi gradient.
Nothing about the kinematic prediction changes. Only the lensing weight moves, from
the dust value 1 to (1 + gamma_chi)/2 -- which is exactly the R11 form, now with a
derived meaning and a derived allowed band instead of a free coefficient.

HONESTY ABOUT WHAT THIS CAN AND CANNOT TEST
-------------------------------------------
A per-lens gamma_chi is six numbers fitted to six angles and is NOT a test. It is
reported as a measurement of what the data demands. The test is whether those six
numbers (a) fall inside the energy-condition band and (b) follow one universal law
of a source property. Both are declared before the law search runs.

    python stress.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.linalg import cho_factor, cho_solve

import model as M

# Energy-condition band, declared before any number is computed.
# Dominant energy condition |p| <= rho gives each w in [-1, 1], so the Tolman
# factor 1 + w_r + 2 w_t lies in [-2, 3]; requiring the companion to attract stars
# forces it positive, hence gamma_chi >= 1/3. gamma_chi = 1 is static dust.
GAMMA_MIN_DEC = 1.0 / 3.0
GAMMA_DUST = 1.0
OBS = {x['Name']: x for x in json.loads(
    (M.ROOT / 'research_work/results/lensing-data-readiness/'
     'lens-observations-and-image-models.json').read_text())}
LIGHT = {x['Name']: x for x in json.loads(
    (M.ROOT / 'research_work/results/slacs-light-profile-audit/results.json').read_text())['rows']}


def companion_terms(lens, dm):
    """The frozen R10 companion force pieces plus the baryon scale at this offset."""
    local = lens.local_parameters(M.params(dm, 0.), M.SPEC)
    A, rc, rt = M.J.source_parameters(lens.Mstar, 0., lens.Re, local, M.SPEC['saturation'], 1.)
    return 10 ** local['logu'], float(A), float(rc), float(rt), float(local['q_sph'])


def angle_with_gamma(lens, dm, gamma):
    """Einstein angle when the companion's lensing weight is (1+gamma)/2, not 1.

    Baryons keep the standard weight. Only the companion's stress state changes.
    """
    scale, A, rc, rt, q = companion_terms(lens, dm)
    weight = 0.5 * (1.0 + gamma)

    def f(b):
        r = b / np.cos(lens.lens_t)
        g = scale * lens.baryon_force(r) + weight * M.J.companion_force(r, A, rc, rt, q)
        return lens.ratio * 4 / M.J.CLIGHT ** 2 * np.dot(lens.lens_w, g * r) - b / lens.Dl
    grid = lens.Re * np.geomspace(1e-6, 1e4, 200)
    vals = np.array([f(b) for b in grid])
    ii = np.flatnonzero((vals[:-1] > 0) & (vals[1:] <= 0))
    if not len(ii):
        return None
    return float(brentq(f, grid[ii[-1]], grid[ii[-1] + 1], xtol=1e-10) / lens.Dl * M.J.ARCSEC)


def kinematic_only_nuisance(lens):
    """dm and beta from the resolved V_rms alone. The lensing weight never enters."""
    chol = cho_factor(np.asarray(lens.cov), lower=True)
    lim = 5.0 * lens.mass_log_error

    def cost(x):
        dm = float(np.clip(x[0], -lim, lim))
        beta = float(np.clip(x[1], *M.BETA_BOUNDS))
        v = M.vrms(lens, dm, beta)
        if v is None:
            return 1e9
        d = v - lens.y
        return float(d @ cho_solve(chol, d)) + (dm / lens.mass_log_error) ** 2 \
            + (beta / M.BETA_PRIOR_SIGMA) ** 2
    best = min((minimize(cost, [d0, b0], method='Nelder-Mead',
                         options=dict(xatol=1e-7, fatol=1e-9, maxiter=900))
                for d0 in (-0.05, 0.0, 0.08) for b0 in (-0.3, 0.0, 0.3)), key=lambda r: r.fun)
    dm = float(np.clip(best.x[0], -lim, lim))
    beta = float(np.clip(best.x[1], *M.BETA_BOUNDS))
    v = M.vrms(lens, dm, beta)
    d = v - lens.y
    return dm, beta, float(d @ cho_solve(chol, d)), v, chol


def required_gamma(lens, dm):
    """The single stress ratio that reproduces the catalog Einstein angle exactly."""
    f = lambda g: (angle_with_gamma(lens, dm, g) or -1e9) - lens.theta
    lo, hi = -0.99, 60.0
    if f(lo) * f(hi) > 0:
        return None
    return float(brentq(f, lo, hi, xtol=1e-12, rtol=1e-14))


def source_properties(lens, dm, v_pred):
    """Candidate predictors, every one independent of the lensing measurement."""
    o = OBS[lens.name]
    scale, A, rc, rt, q = companion_terms(lens, dm)
    Mb = lens.Mstar * scale
    comps = LIGHT[lens.name]['components']
    # projected stellar surface density inside the half-light radius, Msun/kpc^2
    sigma_e = Mb / (2 * np.pi * lens.Re ** 2)
    return dict(
        zFG=o['zFG'], zBG=o['zBG'],
        Mstar_Msun=float(lens.Mstar), Mbaryon_Msun=float(Mb), Re_kpc=float(lens.Re),
        compactness=float((Mb / 1e10) / lens.Re ** 2),
        surface_density_Msun_kpc2=float(sigma_e),
        aperture_sigma_kms=float(o['sigma']),
        mean_vrms_kms=float(np.mean(v_pred)),
        central_vrms_kms=float(v_pred[0]),
        vrms_gradient=float((v_pred[-1] - v_pred[0]) / v_pred[0]),
        axis_ratio=float(o['b/a']), sersic_n_inner=float(comps[0]['n']),
        inner_light_fraction=float(comps[0]['total_light_fraction']),
        companion_core_kpc=float(rc), companion_taper_kpc=float(rt),
        einstein_radius_kpc=float(lens.b), theta_E_over_Re=float(lens.b / lens.Re),
        # dimensionless depth of the baryonic potential well
        potential_depth=float(M.J.G * Mb / (lens.Re * M.J.CLIGHT ** 2)))


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    rows = []
    for name in M.LENSES:
        lens = M.make_lens(name)
        dm, beta, chi2, v, _ = kinematic_only_nuisance(lens)
        g_dust = angle_with_gamma(lens, dm, GAMMA_DUST)
        g_req = required_gamma(lens, dm)
        tolman = None if g_req is None else 1.0 / g_req       # 1 + w_r + 2 w_t
        rows.append(dict(
            name=name, dm_from_kinematics_dex=dm, beta_from_kinematics=beta,
            kinematic_chi2=chi2, n_bins=int(len(lens.y)),
            vrms_fractional_RMS=float(np.sqrt(np.mean(((v - lens.y) / lens.y) ** 2))),
            theta_observed_arcsec=float(lens.theta),
            theta_at_dust_gamma1_arcsec=g_dust,
            angle_error_at_dust=float(g_dust / lens.theta - 1),
            required_gamma_chi=g_req,
            implied_tolman_factor=tolman,
            implied_wr_if_transverse_zero=None if tolman is None else float(tolman - 1),
            inside_energy_condition_band=None if g_req is None else bool(g_req >= GAMMA_MIN_DEC),
            **source_properties(lens, dm, v)))
        r = rows[-1]
        print(f"{name}: dust error {r['angle_error_at_dust']:+.2%}  "
              f"required gamma_chi = {r['required_gamma_chi']:.4f}  "
              f"Tolman factor = {r['implied_tolman_factor']:+.4f}  "
              f"{'IN band' if r['inside_energy_condition_band'] else 'OUT of band'}", flush=True)

    payload = dict(
        experiment='JR-10 stage 1',
        scope='What stress state each lens demands, measured under the JR-9 cross-prediction '
              'protocol. Six numbers from six angles: a measurement of the requirement, not '
              'a test of a model.',
        derivation=dict(
            relation='gamma_chi = Psi/Phi = 1/(1 + w_r + 2 w_t)',
            lensing_weight='(1 + gamma_chi)/2 on the companion force; baryons keep weight 1',
            declared_band=dict(gamma_min_from_energy_conditions=GAMMA_MIN_DEC,
                               gamma_static_dust=GAMMA_DUST,
                               note='Dominant energy condition plus attraction to stars.'),
            reference_states={'static dust': 1.0, 'radial null flux': 0.5,
                              'radial tension w_r=-1/2': 2.0, 'static field gradient': 'infinite'}),
        kinematics_untouched='R10 g_chi was fitted to stellar kinematics and IS the Phi gradient. '
                             'No V_rms prediction changes for any gamma_chi.',
        rows=rows,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'required-stress.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    print('WROTE', out / 'required-stress.json')


if __name__ == '__main__':
    main()
